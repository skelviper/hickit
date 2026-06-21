#!/usr/bin/env python3
"""Audit P9016 contacts.seg sources with duplicate-aware read statistics."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    rest = chrom[3:] if chrom.startswith("chr") else chrom
    if rest.isdigit():
        return (int(rest), "")
    if rest == "X":
        return (1000, "")
    if rest == "Y":
        return (1001, "")
    return (2000, rest)


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.9g}"
    return str(value)


def ratio(num: int | float, den: int | float) -> float:
    return float(num) / float(den) if den else float("nan")


def bucket(n: int) -> str:
    if n <= 10:
        return str(n)
    if n <= 20:
        return "11-20"
    if n <= 50:
        return "21-50"
    if n <= 100:
        return "51-100"
    return ">100"


@dataclass(frozen=True)
class Segment:
    chrom: str
    start: int
    end: int
    strand: str
    phase: str
    mapq: int | None
    count: int | None
    mate: str

    @property
    def midpoint(self) -> int:
        return (self.start + self.end) // 2

    def exact_key(self) -> tuple[str, int, int, str, str, int | None, int | None, str]:
        return (self.chrom, self.start, self.end, self.strand, self.phase, self.mapq, self.count, self.mate)

    def interval_key(self) -> tuple[str, int, int, str]:
        return (self.chrom, self.start, self.end, self.strand)

    def bin_key(self, bin_size: int) -> tuple[str, int]:
        return (self.chrom, self.midpoint // bin_size)


def parse_segment(field: str) -> Segment | None:
    parts = field.split("!")
    if len(parts) < 7:
        return None
    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError:
        return None
    try:
        mapq = int(parts[5])
    except ValueError:
        mapq = None
    try:
        count = int(parts[6])
    except ValueError:
        count = None
    return Segment(
        chrom=parts[0],
        start=start,
        end=end,
        strand=parts[3] if len(parts) > 3 else ".",
        phase=parts[4] if len(parts) > 4 else ".",
        mapq=mapq,
        count=count,
        mate=parts[7] if len(parts) > 7 else ".",
    )


def canonical_segment_pair(a: Segment, b: Segment) -> tuple[Segment, Segment]:
    ak = (chrom_sort_key(a.chrom), a.midpoint, a.start, a.end, a.strand, a.mate)
    bk = (chrom_sort_key(b.chrom), b.midpoint, b.start, b.end, b.strand, b.mate)
    return (a, b) if ak <= bk else (b, a)


def exact_pair_key(a: Segment, b: Segment) -> tuple[tuple[str, int, int, str, str, int | None, int | None, str], tuple[str, int, int, str, str, int | None, int | None, str]]:
    ca, cb = canonical_segment_pair(a, b)
    return ca.exact_key(), cb.exact_key()


def interval_pair_key(a: Segment, b: Segment) -> tuple[tuple[str, int, int, str], tuple[str, int, int, str]]:
    ca, cb = canonical_segment_pair(a, b)
    return ca.interval_key(), cb.interval_key()


def bpair_key(a: Segment, b: Segment, bin_size: int) -> tuple[str, int, str, int]:
    ka = a.bin_key(bin_size)
    kb = b.bin_key(bin_size)
    if (chrom_sort_key(ka[0]), ka[1]) <= (chrom_sort_key(kb[0]), kb[1]):
        return ka[0], ka[1], kb[0], kb[1]
    return kb[0], kb[1], ka[0], ka[1]


def write_kv(path: Path, values: list[tuple[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        for key, value in values:
            writer.writerow([key, fmt(value)])


def write_rows(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in fields})


class SourceAudit:
    def __init__(self, label: str, path: Path, bin_size: int, min_mapq: int, max_examples: int) -> None:
        self.label = label
        self.path = path
        self.bin_size = bin_size
        self.min_mapq = min_mapq
        self.max_examples = max_examples
        self.summary: dict[str, object] = {}
        self.segment_hist: Counter[str] = Counter()
        self.usable_segment_hist: Counter[str] = Counter()
        self.unique_segment_hist: Counter[str] = Counter()
        self.read_chrom_hist: Counter[str] = Counter()
        self.duplicate_read_hist: Counter[str] = Counter()
        self.bpair_count_hist: dict[str, Counter[str]] = {"all": Counter(), "cis": Counter(), "trans": Counter()}
        self.examples: list[dict[str, object]] = []

    def run(self) -> None:
        n_reads = 0
        total_segments = 0
        usable_segments = 0
        unique_exact_segments_total = 0
        unique_interval_segments_total = 0
        unique_bin_segments_total = 0
        duplicate_exact_segments_total = 0
        duplicate_interval_segments_total = 0
        duplicate_bin_segments_total = 0
        n_reads_with_exact_duplicate_segments = 0
        n_reads_with_interval_duplicate_segments = 0
        n_reads_with_bin_duplicate_segments = 0
        n_usable_reads_ge2_segments = 0
        n_pairwise_reads = 0
        n_multiseg_reads = 0
        n_multichrom_reads = 0
        n_multiseg_multichrom_reads = 0
        n_reads_with_duplicate_exact_pairs = 0
        n_reads_with_duplicate_interval_pairs = 0
        n_reads_with_duplicate_bpairs = 0
        n_reads_multiway_after_exact_dedup = 0
        n_reads_multiway_after_interval_dedup = 0
        n_reads_multiway_after_bin_dedup = 0
        phase_labeled_usable_segments = 0
        total_pair_candidates = 0
        total_pair_candidates_cis = 0
        total_pair_candidates_trans = 0
        unique_pair_candidates_within_read = 0
        duplicate_pair_candidates_within_read = 0
        unique_interval_pairs_within_read = 0
        duplicate_interval_pairs_within_read = 0
        unique_bpairs_within_read = 0
        duplicate_bpairs_within_read = 0
        bpair_counts: Counter[tuple[str, int, str, int]] = Counter()
        read_id_counts: Counter[str] = Counter()

        with open_text(self.path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    continue
                read_id = fields[0]
                raw_segments = [parse_segment(field) for field in fields[1:]]
                segments = [seg for seg in raw_segments if seg is not None]
                usable = [seg for seg in segments if seg.mapq is None or seg.mapq >= self.min_mapq]
                n_reads += 1
                read_id_counts[read_id] += 1
                total_segments += len(segments)
                usable_segments += len(usable)
                self.segment_hist[bucket(len(segments))] += 1
                self.usable_segment_hist[bucket(len(usable))] += 1
                exact_keys = [seg.exact_key() for seg in usable]
                interval_keys = [seg.interval_key() for seg in usable]
                bin_keys = [seg.bin_key(self.bin_size) for seg in usable]
                unique_exact = len(set(exact_keys))
                unique_interval = len(set(interval_keys))
                unique_bin = len(set(bin_keys))
                unique_exact_segments_total += unique_exact
                unique_interval_segments_total += unique_interval
                unique_bin_segments_total += unique_bin
                duplicate_exact_segments_total += len(usable) - unique_exact
                duplicate_interval_segments_total += len(usable) - unique_interval
                duplicate_bin_segments_total += len(usable) - unique_bin
                if unique_exact < len(usable):
                    n_reads_with_exact_duplicate_segments += 1
                if unique_interval < len(usable):
                    n_reads_with_interval_duplicate_segments += 1
                if unique_bin < len(usable):
                    n_reads_with_bin_duplicate_segments += 1
                self.unique_segment_hist[bucket(unique_exact)] += 1
                self.duplicate_read_hist[bucket(len(usable) - unique_exact)] += 1
                if unique_exact >= 3:
                    n_reads_multiway_after_exact_dedup += 1
                if unique_interval >= 3:
                    n_reads_multiway_after_interval_dedup += 1
                if unique_bin >= 3:
                    n_reads_multiway_after_bin_dedup += 1
                phase_labeled_usable_segments += sum(1 for seg in usable if seg.phase in {"0", "1"})
                chroms = {seg.chrom for seg in usable}
                self.read_chrom_hist[bucket(len(chroms))] += 1
                if len(usable) < 2:
                    continue
                n_usable_reads_ge2_segments += 1
                if len(usable) == 2:
                    n_pairwise_reads += 1
                else:
                    n_multiseg_reads += 1
                if len(chroms) > 1:
                    n_multichrom_reads += 1
                    if len(usable) > 2:
                        n_multiseg_multichrom_reads += 1
                read_exact_pairs: list[object] = []
                read_interval_pairs: list[object] = []
                read_bpairs: list[tuple[str, int, str, int]] = []
                for a, b in itertools.combinations(usable, 2):
                    total_pair_candidates += 1
                    if a.chrom == b.chrom:
                        total_pair_candidates_cis += 1
                    else:
                        total_pair_candidates_trans += 1
                    ep = exact_pair_key(a, b)
                    ip = interval_pair_key(a, b)
                    bp = bpair_key(a, b, self.bin_size)
                    read_exact_pairs.append(ep)
                    read_interval_pairs.append(ip)
                    read_bpairs.append(bp)
                    bpair_counts[bp] += 1
                    if len(self.examples) < self.max_examples and a.chrom != b.chrom:
                        ca, cb = canonical_segment_pair(a, b)
                        self.examples.append(
                            {
                                "source_label": self.label,
                                "readID": read_id,
                                "read_n_segments": len(usable),
                                "read_n_unique_exact_segments": unique_exact,
                                "read_n_unique_bin_segments": unique_bin,
                                "chr1": ca.chrom,
                                "pos1": ca.midpoint,
                                "chr2": cb.chrom,
                                "pos2": cb.midpoint,
                                "mate1": ca.mate,
                                "mate2": cb.mate,
                            }
                        )
                n_exact_pairs = len(read_exact_pairs)
                n_unique_exact_pairs = len(set(read_exact_pairs))
                n_unique_interval_pairs = len(set(read_interval_pairs))
                n_unique_bpairs = len(set(read_bpairs))
                unique_pair_candidates_within_read += n_unique_exact_pairs
                duplicate_pair_candidates_within_read += n_exact_pairs - n_unique_exact_pairs
                unique_interval_pairs_within_read += n_unique_interval_pairs
                duplicate_interval_pairs_within_read += n_exact_pairs - n_unique_interval_pairs
                unique_bpairs_within_read += n_unique_bpairs
                duplicate_bpairs_within_read += n_exact_pairs - n_unique_bpairs
                if n_unique_exact_pairs < n_exact_pairs:
                    n_reads_with_duplicate_exact_pairs += 1
                if n_unique_interval_pairs < n_exact_pairs:
                    n_reads_with_duplicate_interval_pairs += 1
                if n_unique_bpairs < n_exact_pairs:
                    n_reads_with_duplicate_bpairs += 1

        duplicate_read_ids = sum(1 for count in read_id_counts.values() if count > 1)
        repeated_read_id_rows = sum(count - 1 for count in read_id_counts.values() if count > 1)
        for key, count in bpair_counts.items():
            scope = "cis" if key[0] == key[2] else "trans"
            self.bpair_count_hist["all"][bucket(count)] += 1
            self.bpair_count_hist[scope][bucket(count)] += 1

        n_unique_bpairs_all = len(bpair_counts)
        n_unique_bpairs_cis = sum(1 for key in bpair_counts if key[0] == key[2])
        n_unique_bpairs_trans = sum(1 for key in bpair_counts if key[0] != key[2])
        self.summary = {
            "source_label": self.label,
            "seg_path": str(self.path),
            "file_size_bytes": self.path.stat().st_size if self.path.exists() else "NA",
            "bin_size_bp": self.bin_size,
            "min_mapq": self.min_mapq,
            "n_reads": n_reads,
            "n_duplicate_read_ids": duplicate_read_ids,
            "n_repeated_read_id_rows": repeated_read_id_rows,
            "total_segments": total_segments,
            "usable_segments": usable_segments,
            "n_usable_reads_ge2_segments": n_usable_reads_ge2_segments,
            "n_pairwise_reads": n_pairwise_reads,
            "n_multiseg_reads": n_multiseg_reads,
            "n_multichrom_reads": n_multichrom_reads,
            "n_multiseg_multichrom_reads": n_multiseg_multichrom_reads,
            "frac_multiseg_of_usable_ge2": ratio(n_multiseg_reads, n_usable_reads_ge2_segments),
            "frac_multichrom_of_usable_ge2": ratio(n_multichrom_reads, n_usable_reads_ge2_segments),
            "frac_multiseg_multichrom_of_usable_ge2": ratio(n_multiseg_multichrom_reads, n_usable_reads_ge2_segments),
            "unique_exact_segments_total": unique_exact_segments_total,
            "unique_interval_segments_total": unique_interval_segments_total,
            "unique_1mb_segment_bins_total": unique_bin_segments_total,
            "duplicate_exact_segments_within_read": duplicate_exact_segments_total,
            "duplicate_interval_segments_within_read": duplicate_interval_segments_total,
            "duplicate_1mb_segment_bins_within_read": duplicate_bin_segments_total,
            "frac_duplicate_exact_segments_within_read": ratio(duplicate_exact_segments_total, usable_segments),
            "frac_duplicate_interval_segments_within_read": ratio(duplicate_interval_segments_total, usable_segments),
            "frac_duplicate_1mb_segment_bins_within_read": ratio(duplicate_bin_segments_total, usable_segments),
            "n_reads_with_exact_duplicate_segments": n_reads_with_exact_duplicate_segments,
            "n_reads_with_interval_duplicate_segments": n_reads_with_interval_duplicate_segments,
            "n_reads_with_1mb_duplicate_segments": n_reads_with_bin_duplicate_segments,
            "n_reads_multiway_after_exact_dedup": n_reads_multiway_after_exact_dedup,
            "n_reads_multiway_after_interval_dedup": n_reads_multiway_after_interval_dedup,
            "n_reads_multiway_after_1mb_bin_dedup": n_reads_multiway_after_bin_dedup,
            "phase_labeled_usable_segments": phase_labeled_usable_segments,
            "phase_labels_present_in_source_for_eval_only": int(phase_labeled_usable_segments > 0),
            "training_export_should_drop_phase_labels": 1,
            "pair_candidates_all": total_pair_candidates,
            "pair_candidates_cis": total_pair_candidates_cis,
            "pair_candidates_trans": total_pair_candidates_trans,
            "unique_exact_pair_candidates_within_read": unique_pair_candidates_within_read,
            "duplicate_exact_pair_candidates_within_read": duplicate_pair_candidates_within_read,
            "unique_interval_pair_candidates_within_read": unique_interval_pairs_within_read,
            "duplicate_interval_pair_candidates_within_read": duplicate_interval_pairs_within_read,
            "unique_1mb_bpair_candidates_within_read": unique_bpairs_within_read,
            "duplicate_1mb_bpair_candidates_within_read": duplicate_bpairs_within_read,
            "frac_duplicate_exact_pair_candidates_within_read": ratio(duplicate_pair_candidates_within_read, total_pair_candidates),
            "frac_duplicate_interval_pair_candidates_within_read": ratio(duplicate_interval_pairs_within_read, total_pair_candidates),
            "frac_duplicate_1mb_bpair_candidates_within_read": ratio(duplicate_bpairs_within_read, total_pair_candidates),
            "n_reads_with_duplicate_exact_pairs": n_reads_with_duplicate_exact_pairs,
            "n_reads_with_duplicate_interval_pairs": n_reads_with_duplicate_interval_pairs,
            "n_reads_with_duplicate_1mb_bpairs": n_reads_with_duplicate_bpairs,
            "unique_1mb_bpairs_all": n_unique_bpairs_all,
            "unique_1mb_bpairs_cis": n_unique_bpairs_cis,
            "unique_1mb_bpairs_trans": n_unique_bpairs_trans,
            "bpair_raw_to_unique_ratio_all": ratio(total_pair_candidates, n_unique_bpairs_all),
            "bpair_raw_to_unique_ratio_cis": ratio(total_pair_candidates_cis, n_unique_bpairs_cis),
            "bpair_raw_to_unique_ratio_trans": ratio(total_pair_candidates_trans, n_unique_bpairs_trans),
            "read_level_anchor_status": "POTENTIAL_READ_LEVEL_MULTISEG_ANCHOR"
            if n_multiseg_multichrom_reads > 0
            else "NO_MULTISEG_MULTICHROM_ANCHOR",
        }

    def write_outputs(self, outdir: Path) -> None:
        source_dir = outdir / self.label
        write_kv(source_dir / "source_summary.tsv", list(self.summary.items()))
        write_rows(
            source_dir / "segment_count_hist.tsv",
            [{"bucket": k, "n_reads": v} for k, v in sorted(self.segment_hist.items())],
            ["bucket", "n_reads"],
        )
        write_rows(
            source_dir / "usable_segment_count_hist.tsv",
            [{"bucket": k, "n_reads": v} for k, v in sorted(self.usable_segment_hist.items())],
            ["bucket", "n_reads"],
        )
        write_rows(
            source_dir / "unique_exact_segment_count_hist.tsv",
            [{"bucket": k, "n_reads": v} for k, v in sorted(self.unique_segment_hist.items())],
            ["bucket", "n_reads"],
        )
        write_rows(
            source_dir / "read_chrom_count_hist.tsv",
            [{"bucket": k, "n_reads": v} for k, v in sorted(self.read_chrom_hist.items())],
            ["bucket", "n_reads"],
        )
        write_rows(
            source_dir / "duplicate_segment_count_hist.tsv",
            [{"bucket": k, "n_reads": v} for k, v in sorted(self.duplicate_read_hist.items())],
            ["bucket", "n_reads"],
        )
        bpair_rows: list[dict[str, object]] = []
        for scope, hist in self.bpair_count_hist.items():
            for key, value in sorted(hist.items()):
                bpair_rows.append({"scope": scope, "bpair_contact_count_bucket": key, "n_1mb_bpairs": value})
        write_rows(source_dir / "bpair_contact_count_hist.tsv", bpair_rows, ["scope", "bpair_contact_count_bucket", "n_1mb_bpairs"])
        write_rows(
            source_dir / "trans_pair_examples.tsv",
            self.examples,
            [
                "source_label",
                "readID",
                "read_n_segments",
                "read_n_unique_exact_segments",
                "read_n_unique_bin_segments",
                "chr1",
                "pos1",
                "chr2",
                "pos2",
                "mate1",
                "mate2",
            ],
        )


def plot_summary(outdir: Path, summaries: list[dict[str, object]]) -> None:
    plots = outdir / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        (plots / "PLOT_SKIPPED.txt").write_text(f"matplotlib unavailable: {exc}\n")
        return
    labels = [str(s["source_label"]) for s in summaries]
    metrics = [
        ("frac_multiseg_of_usable_ge2", "multi-seg"),
        ("frac_multichrom_of_usable_ge2", "multi-chrom"),
        ("frac_multiseg_multichrom_of_usable_ge2", "multi-seg+multi-chrom"),
        ("frac_duplicate_exact_segments_within_read", "dup exact seg"),
        ("frac_duplicate_1mb_bpair_candidates_within_read", "dup 1Mb bpair"),
    ]
    fig, ax = plt.subplots(figsize=(7.5, 3.5), dpi=300)
    width = 0.8 / len(metrics)
    x = list(range(len(labels)))
    for i, (key, name) in enumerate(metrics):
        vals = [float(s.get(key, float("nan"))) for s in summaries]
        offsets = [xx - 0.4 + width * (i + 0.5) for xx in x]
        ax.bar(offsets, vals, width=width, label=name)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_ylabel("fraction", fontsize=7)
    ax.set_title("P9016 contacts.seg source audit", fontsize=7)
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(plots / "seg_source_fraction_summary.png")
    plt.close(fig)


def write_readme(outdir: Path, summaries: list[dict[str, object]]) -> None:
    lines = [
        "# 053 P9016 contacts.seg Source Audit 1Mb",
        "",
        "This is a source audit only. No reconstruction was trained here.",
        "",
        "## Purpose",
        "",
        "The goal is to decide whether upstream `contacts.seg` contains read-level multiway information that is absent from the stripped P9016 `.pairs.gz` file, while explicitly accounting for the fact that the segment file is not deduplicated.",
        "",
        "## Why Duplicate Handling Matters",
        "",
        "If one read contains repeated or overlapping segment records, expanding all segment pairs can create artificial multiway contacts. This audit therefore reports exact segment duplicates, interval-level duplicates, 1Mb bin duplicates, exact pair duplicates, and collapsed 1Mb bpair duplicates before any training experiment is considered.",
        "",
        "## Source Summary",
        "",
        "| source | reads | usable >=2 seg | multi-seg frac | multi-chrom frac | multi-seg+multi-chrom frac | exact seg dup frac | 1Mb bpair dup frac | phase labels present |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for summary in summaries:
        lines.append(
            "| {source} | {reads} | {usable} | {multi:.4g} | {mc:.4g} | {both:.4g} | {dupseg:.4g} | {dupbp:.4g} | {phase} |".format(
                source=summary["source_label"],
                reads=summary["n_reads"],
                usable=summary["n_usable_reads_ge2_segments"],
                multi=float(summary["frac_multiseg_of_usable_ge2"]),
                mc=float(summary["frac_multichrom_of_usable_ge2"]),
                both=float(summary["frac_multiseg_multichrom_of_usable_ge2"]),
                dupseg=float(summary["frac_duplicate_exact_segments_within_read"]),
                dupbp=float(summary["frac_duplicate_1mb_bpair_candidates_within_read"]),
                phase=summary["phase_labels_present_in_source_for_eval_only"],
            )
        )
    lines.extend(
        [
            "",
            "## Training Boundary",
            "",
            "- This audit reads upstream segment sources but does not train a model.",
            "- Phase-like segment fields, when present, must be dropped for any blind training export. In this audit the sharec source has phase labels, while the archive source did not show phase labels.",
            "- The next training experiment should treat this as an explicit source-boundary change, not as the strict P9016 `.pairs.gz` baseline.",
            "",
            "## Proposed Next Step",
            "",
            "If these source-level duplicate statistics are acceptable, the next experiment should be a read-level posterior diagnostic: use current model posteriors plus read grouping to test whether read-consistent multi-segment decoding improves trans identity before changing EM/FDG training.",
        ]
    )
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", action="append", required=True, help="label=/path/to/contacts.seg.gz")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--bin-size", type=int, default=1_000_000)
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--max-examples", type=int, default=30)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    audits: list[SourceAudit] = []
    for source in args.source:
        if "=" not in source:
            raise SystemExit(f"--source must be label=/path form: {source}")
        label, raw_path = source.split("=", 1)
        audit = SourceAudit(label=label, path=Path(raw_path), bin_size=args.bin_size, min_mapq=args.min_mapq, max_examples=args.max_examples)
        audit.run()
        audit.write_outputs(args.outdir)
        audits.append(audit)

    fields = sorted({key for audit in audits for key in audit.summary})
    if "source_label" in fields:
        fields.remove("source_label")
    fields = ["source_label"] + fields
    write_rows(args.outdir / "seg_source_comparison.tsv", [audit.summary for audit in audits], fields)
    write_kv(
        args.outdir / "summary.tsv",
        [
            ("n_sources", len(audits)),
            ("bin_size_bp", args.bin_size),
            ("min_mapq", args.min_mapq),
            ("has_source_with_multiseg_multichrom", int(any(int(a.summary["n_multiseg_multichrom_reads"]) > 0 for a in audits))),
            ("training_run", 0),
            ("recommendation", "run_read_level_posterior_diagnostic_before_training"),
        ],
    )
    plot_summary(args.outdir, [audit.summary for audit in audits])
    write_readme(args.outdir, [audit.summary for audit in audits])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
