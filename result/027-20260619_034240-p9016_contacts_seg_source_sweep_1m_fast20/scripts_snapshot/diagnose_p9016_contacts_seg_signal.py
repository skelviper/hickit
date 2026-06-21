#!/usr/bin/env python3
"""Audit read-level/multi-segment signal in a CHARM contacts.seg.gz file."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO


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


def segment_midpoint(field: str) -> tuple[str, int, str, str | None, int | None, int | None] | None:
    parts = field.split("!")
    if len(parts) < 7:
        return None
    chrom = parts[0]
    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError:
        return None
    strand = parts[3] if len(parts) > 3 else "."
    phase = parts[4] if len(parts) > 4 else "."
    try:
        mapq = int(parts[5])
    except ValueError:
        mapq = None
    try:
        count = int(parts[6])
    except ValueError:
        count = None
    return chrom, (start + end) // 2, strand, phase, mapq, count


def pair_key(seg_a: tuple[str, int, str, str | None, int | None, int | None], seg_b: tuple[str, int, str, str | None, int | None, int | None]) -> tuple[str, int, str, int]:
    a = (chrom_sort_key(seg_a[0]), seg_a[1], seg_a[0])
    b = (chrom_sort_key(seg_b[0]), seg_b[1], seg_b[0])
    if a <= b:
        return seg_a[0], seg_a[1], seg_b[0], seg_b[1]
    return seg_b[0], seg_b[1], seg_a[0], seg_a[1]


def bpair_key(seg_a: tuple[str, int, str, str | None, int | None, int | None], seg_b: tuple[str, int, str, str | None, int | None, int | None], bin_size: int) -> tuple[str, int, str, int]:
    a = (seg_a[0], seg_a[1] // bin_size)
    b = (seg_b[0], seg_b[1] // bin_size)
    if (chrom_sort_key(a[0]), a[1]) <= (chrom_sort_key(b[0]), b[1]):
        return a[0], a[1], b[0], b[1]
    return b[0], b[1], a[0], a[1]


def write_kv(path: Path, values: list[tuple[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        for key, value in values:
            writer.writerow([key, fmt(value)])


def write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in fieldnames})


def bucket(n: int) -> str:
    if n <= 10:
        return str(n)
    if n <= 20:
        return "11-20"
    if n <= 50:
        return "21-50"
    return ">50"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seg", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--bin-size", type=int, default=1_000_000)
    ap.add_argument("--sample", default="P9016")
    ap.add_argument("--min-mapq", type=int, default=30)
    args = ap.parse_args()

    n_reads = 0
    n_usable_reads = 0
    n_pairwise_reads = 0
    n_multiseg_reads = 0
    n_multichrom_reads = 0
    n_multiseg_multichrom_reads = 0
    total_segments = 0
    usable_segments = 0
    total_export_pairs = 0
    trans_export_pairs = 0
    cis_export_pairs = 0
    phase_labeled_segments = 0
    segment_hist: Counter[str] = Counter()
    usable_segment_hist: Counter[str] = Counter()
    read_chrom_hist: Counter[str] = Counter()
    unique_bpair_counts: Counter[tuple[str, int, str, int]] = Counter()
    unique_pair_examples: list[dict[str, object]] = []

    with open_text(args.seg) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            read_id = fields[0]
            raw_segments = [segment_midpoint(field) for field in fields[1:]]
            segments = [seg for seg in raw_segments if seg is not None]
            n_reads += 1
            total_segments += len(segments)
            segment_hist[bucket(len(segments))] += 1
            usable = [seg for seg in segments if seg[4] is None or seg[4] >= args.min_mapq]
            usable_segments += len(usable)
            usable_segment_hist[bucket(len(usable))] += 1
            phase_labeled_segments += sum(1 for seg in usable if seg[3] in {"0", "1"})
            chroms = {seg[0] for seg in usable}
            read_chrom_hist[bucket(len(chroms))] += 1
            if len(usable) < 2:
                continue
            n_usable_reads += 1
            if len(usable) == 2:
                n_pairwise_reads += 1
            else:
                n_multiseg_reads += 1
            if len(chroms) > 1:
                n_multichrom_reads += 1
                if len(usable) > 2:
                    n_multiseg_multichrom_reads += 1
            for seg_a, seg_b in itertools.combinations(usable, 2):
                total_export_pairs += 1
                if seg_a[0] == seg_b[0]:
                    cis_export_pairs += 1
                else:
                    trans_export_pairs += 1
                key = bpair_key(seg_a, seg_b, args.bin_size)
                unique_bpair_counts[key] += 1
                if len(unique_pair_examples) < 20 and seg_a[0] != seg_b[0]:
                    chrom1, pos1, chrom2, pos2 = pair_key(seg_a, seg_b)
                    unique_pair_examples.append(
                        {
                            "readID": read_id,
                            "chr1": chrom1,
                            "pos1": pos1,
                            "chr2": chrom2,
                            "pos2": pos2,
                            "strand1": "+",
                            "strand2": "+",
                        }
                    )

    bpair_contact_hist: Counter[str] = Counter()
    trans_bpair_contact_hist: Counter[str] = Counter()
    cis_bpair_contact_hist: Counter[str] = Counter()
    for (chrom1, _, chrom2, _), count in unique_bpair_counts.items():
        bpair_contact_hist[bucket(count)] += 1
        if chrom1 == chrom2:
            cis_bpair_contact_hist[bucket(count)] += 1
        else:
            trans_bpair_contact_hist[bucket(count)] += 1

    summary = [
        ("sample", args.sample),
        ("seg_path", str(args.seg)),
        ("bin_size_bp", args.bin_size),
        ("min_mapq", args.min_mapq),
        ("n_reads", n_reads),
        ("n_usable_reads_ge2_segments", n_usable_reads),
        ("n_pairwise_reads", n_pairwise_reads),
        ("n_multiseg_reads", n_multiseg_reads),
        ("n_multichrom_reads", n_multichrom_reads),
        ("n_multiseg_multichrom_reads", n_multiseg_multichrom_reads),
        ("total_segments", total_segments),
        ("usable_segments", usable_segments),
        ("phase_labeled_usable_segments", phase_labeled_segments),
        ("phase_labels_present_in_source_for_eval_only", 1 if phase_labeled_segments else 0),
        ("training_export_should_drop_phase_labels", 1),
        ("export_pair_candidates_all", total_export_pairs),
        ("export_pair_candidates_cis", cis_export_pairs),
        ("export_pair_candidates_trans", trans_export_pairs),
        ("unique_1mb_bpairs_all", len(unique_bpair_counts)),
        ("unique_1mb_bpairs_cis", sum(1 for key in unique_bpair_counts if key[0] == key[2])),
        ("unique_1mb_bpairs_trans", sum(1 for key in unique_bpair_counts if key[0] != key[2])),
        ("read_level_anchor_status", "POTENTIAL_READ_LEVEL_MULTISEG_ANCHOR" if n_multiseg_multichrom_reads > 0 else "NO_MULTISEG_MULTICHROM_ANCHOR"),
    ]
    write_kv(args.outdir / "contacts_seg_signal_summary.tsv", summary)
    write_rows(
        args.outdir / "segment_count_hist.tsv",
        [{"n_segments_bucket": k, "n_reads": v} for k, v in sorted(segment_hist.items())],
        ["n_segments_bucket", "n_reads"],
    )
    write_rows(
        args.outdir / "usable_segment_count_hist.tsv",
        [{"n_usable_segments_bucket": k, "n_reads": v} for k, v in sorted(usable_segment_hist.items())],
        ["n_usable_segments_bucket", "n_reads"],
    )
    write_rows(
        args.outdir / "read_chrom_count_hist.tsv",
        [{"n_chroms_bucket": k, "n_reads": v} for k, v in sorted(read_chrom_hist.items())],
        ["n_chroms_bucket", "n_reads"],
    )
    rows = []
    for scope, hist in (("all", bpair_contact_hist), ("cis", cis_bpair_contact_hist), ("trans", trans_bpair_contact_hist)):
        for key, value in sorted(hist.items()):
            rows.append({"scope": scope, "bpair_contact_count_bucket": key, "n_1mb_bpairs": value})
    write_rows(args.outdir / "bpair_contact_count_hist.tsv", rows, ["scope", "bpair_contact_count_bucket", "n_1mb_bpairs"])
    write_rows(args.outdir / "trans_pair_examples.tsv", unique_pair_examples, ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"])


if __name__ == "__main__":
    main()
