#!/usr/bin/env python3
"""Audit whether a P9016-style pairs file contains blind source signal.

This diagnostic separates two questions:

1. Does the pairs file retain read/group identifiers that could support a
   blind multi-contact or molecule-level constraint?
2. How stable is the measured SNP-labeled dominant state after 1 Mb binning?

Phase labels are used only for the second, eval-only stability audit. They are
not used to construct any training prior here.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, TextIO


MISSING_VALUES = {"", ".", "NA", "N/A", "nan", "NaN", "None", "none"}
DEFAULT_COLUMNS = ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2", "phase0", "phase1"]


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    if chrom.startswith("chr"):
        rest = chrom[3:]
    else:
        rest = chrom
    if rest.isdigit():
        return (int(rest), "")
    if rest == "X":
        return (1000, "")
    if rest == "Y":
        return (1001, "")
    if rest in {"M", "MT"}:
        return (1002, "")
    return (2000, rest)


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def parse_phase(value: str) -> int | None:
    if value in MISSING_VALUES:
        return None
    try:
        phase = int(value)
    except ValueError:
        return None
    if phase not in (0, 1):
        return None
    return phase


def state_from_phase(phase0: int, phase1: int) -> int:
    return phase0 * 2 + phase1


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.9g}"
    return str(value)


def safe_div(num: float, den: float) -> float:
    if den == 0:
        return float("nan")
    return num / den


def bucket_count(n: int) -> str:
    if n <= 5:
        return str(n)
    if n <= 10:
        return "6-10"
    if n <= 20:
        return "11-20"
    if n <= 50:
        return "21-50"
    if n <= 100:
        return "51-100"
    return ">100"


def stable_half(parts: Iterable[object]) -> int:
    h = hashlib.blake2b(digest_size=8)
    for part in parts:
        h.update(str(part).encode())
        h.update(b"\0")
    return h.digest()[0] & 1


def iter_rows(path: Path):
    columns: list[str] | None = None
    with open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#columns:"):
                text = line[len("#columns:") :].strip()
                columns = [x for x in text.split("\t") if x]
                continue
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if columns is None:
                columns = DEFAULT_COLUMNS[: len(parts)]
            if len(parts) < len(columns):
                parts = parts + [""] * (len(columns) - len(parts))
            yield dict(zip(columns, parts))


def bpair_key(row: dict[str, str], bin_size: int) -> tuple[str, int, str, int]:
    chrom1 = row["chr1"]
    chrom2 = row["chr2"]
    bin1 = int(row["pos1"]) // bin_size
    bin2 = int(row["pos2"]) // bin_size
    a = (chrom_sort_key(chrom1), bin1, chrom1)
    b = (chrom_sort_key(chrom2), bin2, chrom2)
    if a <= b:
        return (chrom1, bin1, chrom2, bin2)
    return (chrom2, bin2, chrom1, bin1)


def write_kv(path: Path, values: list[tuple[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        for key, value in values:
            writer.writerow([key, fmt(value)])


def write_dict_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in fieldnames})


def readid_status(total_rows: int, nonmissing_rows: int, repeated_rows: int, repeated_groups: int) -> str:
    if total_rows == 0:
        return "NO_ROWS"
    if nonmissing_rows == 0:
        return "UNUSABLE_ALL_READID_MISSING"
    if repeated_rows == 0 or repeated_groups == 0:
        return "UNUSABLE_READIDS_UNIQUE_ONLY"
    frac = repeated_rows / total_rows
    if frac < 0.01:
        return "WEAK_REPEATED_READID_SIGNAL_LT_1PCT_ROWS"
    return "POTENTIALLY_USABLE_REPEATED_READID_GROUPS"


def aggregate_readid_groups(pairs: Path, repeated_readids: set[str]) -> tuple[list[dict[str, object]], dict[str, int]]:
    if not repeated_readids:
        return [], {}
    group_counts: Counter[str] = Counter()
    group_cis: Counter[str] = Counter()
    group_trans: Counter[str] = Counter()
    group_chroms: dict[str, set[str]] = defaultdict(set)
    group_bin_pairs: dict[str, set[tuple[str, int, str, int]]] = defaultdict(set)
    for row in iter_rows(pairs):
        rid = row.get("readID", "")
        if rid not in repeated_readids:
            continue
        group_counts[rid] += 1
        chrom1 = row["chr1"]
        chrom2 = row["chr2"]
        group_chroms[rid].add(chrom1)
        group_chroms[rid].add(chrom2)
        if chrom1 == chrom2:
            group_cis[rid] += 1
        else:
            group_trans[rid] += 1
        try:
            group_bin_pairs[rid].add(bpair_key(row, 1_000_000))
        except Exception:
            pass

    rows: list[dict[str, object]] = []
    summary = {
        "repeated_readid_groups_with_trans": 0,
        "repeated_readid_groups_with_cis_and_trans": 0,
        "repeated_readid_groups_with_ge3_chroms": 0,
        "max_repeated_readid_chroms": 0,
        "max_repeated_readid_unique_1mb_bpairs": 0,
    }
    hist: Counter[tuple[str, str, str]] = Counter()
    for rid, n in group_counts.items():
        n_chroms = len(group_chroms[rid])
        n_unique_bpairs = len(group_bin_pairs[rid])
        has_cis = group_cis[rid] > 0
        has_trans = group_trans[rid] > 0
        if has_trans:
            summary["repeated_readid_groups_with_trans"] += 1
        if has_cis and has_trans:
            summary["repeated_readid_groups_with_cis_and_trans"] += 1
        if n_chroms >= 3:
            summary["repeated_readid_groups_with_ge3_chroms"] += 1
        summary["max_repeated_readid_chroms"] = max(summary["max_repeated_readid_chroms"], n_chroms)
        summary["max_repeated_readid_unique_1mb_bpairs"] = max(summary["max_repeated_readid_unique_1mb_bpairs"], n_unique_bpairs)
        hist[(bucket_count(n), bucket_count(n_chroms), bucket_count(n_unique_bpairs))] += 1

    for (group_size_bucket, chrom_bucket, bpair_bucket), n_groups in sorted(hist.items()):
        rows.append(
            {
                "group_size_bucket": group_size_bucket,
                "n_chroms_bucket": chrom_bucket,
                "n_unique_1mb_bpairs_bucket": bpair_bucket,
                "n_readids": n_groups,
            }
        )
    return rows, summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", required=True, type=Path)
    ap.add_argument("--bin-size", type=int, default=1_000_000)
    ap.add_argument("--sample", default="P9016")
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    n_rows = 0
    n_cis = 0
    n_trans = 0
    n_phase_both = 0
    n_phase_both_cis = 0
    n_phase_both_trans = 0
    n_readid_missing = 0
    readid_counts: Counter[str] = Counter()
    bpair_counts: dict[tuple[str, int, str, int], list[int]] = defaultdict(lambda: [0, 0, 0, 0])
    bpair_halves: dict[tuple[str, int, str, int], list[int]] = defaultdict(lambda: [0, 0, 0, 0, 0, 0, 0, 0])

    for row_index, row in enumerate(iter_rows(args.pairs)):
        n_rows += 1
        chrom1 = row["chr1"]
        chrom2 = row["chr2"]
        is_cis = chrom1 == chrom2
        if is_cis:
            n_cis += 1
        else:
            n_trans += 1

        rid = row.get("readID", "")
        if rid in MISSING_VALUES:
            n_readid_missing += 1
        else:
            readid_counts[rid] += 1

        phase0 = parse_phase(row.get("phase0", "."))
        phase1 = parse_phase(row.get("phase1", "."))
        if phase0 is None or phase1 is None:
            continue
        n_phase_both += 1
        if is_cis:
            n_phase_both_cis += 1
        else:
            n_phase_both_trans += 1
        key = bpair_key(row, args.bin_size)
        state = state_from_phase(phase0, phase1)
        bpair_counts[key][state] += 1
        half = stable_half((row_index, row["chr1"], row["pos1"], row["chr2"], row["pos2"], phase0, phase1))
        bpair_halves[key][half * 4 + state] += 1

    n_readid_nonmissing = sum(readid_counts.values())
    repeated_readids = {rid for rid, n in readid_counts.items() if n >= 2}
    rows_in_repeated_readids = sum(readid_counts[rid] for rid in repeated_readids)
    max_readid_group = max(readid_counts.values(), default=0)
    status = readid_status(n_rows, n_readid_nonmissing, rows_in_repeated_readids, len(repeated_readids))

    group_size_hist: Counter[str] = Counter()
    group_size_row_hist: Counter[str] = Counter()
    for count in readid_counts.values():
        bucket = bucket_count(count)
        group_size_hist[bucket] += 1
        group_size_row_hist[bucket] += count
    readid_hist_rows = [
        {
            "group_size_bucket": bucket,
            "n_readids": group_size_hist[bucket],
            "n_rows": group_size_row_hist[bucket],
        }
        for bucket in sorted(group_size_hist, key=lambda x: (0 if x.isdigit() else 1, int(x) if x.isdigit() else x))
    ]

    repeated_group_rows, repeated_group_summary = aggregate_readid_groups(args.pairs, repeated_readids)

    truth_hist: dict[tuple[str, str], dict[str, float]] = defaultdict(
        lambda: {"n_bpairs": 0, "contact_weight": 0, "dominant_contact_weight": 0, "dominant_fraction_sum": 0.0}
    )
    truth_scope_totals: dict[str, dict[str, float]] = defaultdict(
        lambda: {"n_bpairs": 0, "contact_weight": 0, "singleton_contact_weight": 0, "n_singleton_bpairs": 0}
    )
    split_stats: dict[tuple[str, str], dict[str, float]] = defaultdict(
        lambda: {
            "n_bpairs_with_both_halves": 0,
            "contact_weight": 0,
            "state_agree_weight": 0,
            "same_cross_agree_weight": 0,
        }
    )

    for key, counts in bpair_counts.items():
        chrom1, _, chrom2, _ = key
        scope = "cis" if chrom1 == chrom2 else "trans"
        scopes = ("all", scope)
        n = sum(counts)
        dominant = max(counts)
        bucket = bucket_count(n)
        for s in scopes:
            entry = truth_hist[(s, bucket)]
            entry["n_bpairs"] += 1
            entry["contact_weight"] += n
            entry["dominant_contact_weight"] += dominant
            entry["dominant_fraction_sum"] += dominant / n if n else 0.0
            total = truth_scope_totals[s]
            total["n_bpairs"] += 1
            total["contact_weight"] += n
            if n == 1:
                total["n_singleton_bpairs"] += 1
                total["singleton_contact_weight"] += n

        halves = bpair_halves[key]
        a = halves[:4]
        b = halves[4:]
        na = sum(a)
        nb = sum(b)
        if na <= 0 or nb <= 0:
            continue
        top_a = max(range(4), key=lambda idx: a[idx])
        top_b = max(range(4), key=lambda idx: b[idx])
        state_agree = 1 if top_a == top_b else 0
        same_cross_agree = 1 if (top_a in (0, 3)) == (top_b in (0, 3)) else 0
        weight = n
        for s in scopes:
            entry = split_stats[(s, bucket)]
            entry["n_bpairs_with_both_halves"] += 1
            entry["contact_weight"] += weight
            entry["state_agree_weight"] += weight * state_agree
            entry["same_cross_agree_weight"] += weight * same_cross_agree

    truth_hist_rows: list[dict[str, object]] = []
    for (scope, bucket), entry in sorted(truth_hist.items()):
        n_bpairs = entry["n_bpairs"]
        contact_weight = entry["contact_weight"]
        truth_hist_rows.append(
            {
                "scope": scope,
                "truth_n_bucket": bucket,
                "n_bpairs": int(n_bpairs),
                "contact_weight": int(contact_weight),
                "dominant_contact_weight": int(entry["dominant_contact_weight"]),
                "dominant_contact_fraction_weighted": safe_div(entry["dominant_contact_weight"], contact_weight),
                "dominant_fraction_mean_per_bpair": safe_div(entry["dominant_fraction_sum"], n_bpairs),
            }
        )

    split_rows: list[dict[str, object]] = []
    for (scope, bucket), entry in sorted(split_stats.items()):
        weight = entry["contact_weight"]
        split_rows.append(
            {
                "scope": scope,
                "truth_n_bucket": bucket,
                "n_bpairs_with_both_halves": int(entry["n_bpairs_with_both_halves"]),
                "contact_weight": int(weight),
                "dominant_state_agreement_contact_weighted": safe_div(entry["state_agree_weight"], weight),
                "same_cross_agreement_contact_weighted": safe_div(entry["same_cross_agree_weight"], weight),
                "eval_only_uses_phase_labels": 1,
            }
        )

    summary: list[tuple[str, object]] = [
        ("sample", args.sample),
        ("pairs_path", str(args.pairs)),
        ("bin_size_bp", args.bin_size),
        ("n_rows", n_rows),
        ("n_cis_rows", n_cis),
        ("n_trans_rows", n_trans),
        ("n_phase_both_rows", n_phase_both),
        ("n_phase_both_cis_rows", n_phase_both_cis),
        ("n_phase_both_trans_rows", n_phase_both_trans),
        ("frac_phase_both_rows", safe_div(n_phase_both, n_rows)),
        ("n_readid_missing_rows", n_readid_missing),
        ("n_readid_nonmissing_rows", n_readid_nonmissing),
        ("frac_readid_missing_rows", safe_div(n_readid_missing, n_rows)),
        ("unique_nonmissing_readids", len(readid_counts)),
        ("repeated_nonmissing_readid_groups", len(repeated_readids)),
        ("rows_in_repeated_nonmissing_readids", rows_in_repeated_readids),
        ("frac_rows_in_repeated_nonmissing_readids", safe_div(rows_in_repeated_readids, n_rows)),
        ("max_nonmissing_readid_group_size", max_readid_group),
        ("readid_signal_status", status),
        ("n_1mb_truth_bpairs", len(bpair_counts)),
        ("eval_only_phase_labels_used_for_truth_stability", 1),
        ("training_uses_phase_labels", 0),
        ("training_uses_charm_or_reference", 0),
    ]
    for key, value in repeated_group_summary.items():
        summary.append((key, value))
    for scope in ("all", "cis", "trans"):
        total = truth_scope_totals[scope]
        summary.extend(
            [
                (f"{scope}_truth_bpairs", int(total["n_bpairs"])),
                (f"{scope}_truth_contact_weight", int(total["contact_weight"])),
                (f"{scope}_truth_singleton_bpairs", int(total["n_singleton_bpairs"])),
                (f"{scope}_truth_singleton_contact_weight", int(total["singleton_contact_weight"])),
                (f"{scope}_truth_singleton_contact_fraction", safe_div(total["singleton_contact_weight"], total["contact_weight"])),
            ]
        )

    write_kv(args.outdir / "pairs_source_signal_summary.tsv", summary)
    write_dict_rows(args.outdir / "readid_group_size_hist.tsv", readid_hist_rows, ["group_size_bucket", "n_readids", "n_rows"])
    write_dict_rows(
        args.outdir / "readid_repeated_group_diag.tsv",
        repeated_group_rows,
        ["group_size_bucket", "n_chroms_bucket", "n_unique_1mb_bpairs_bucket", "n_readids"],
    )
    write_dict_rows(
        args.outdir / "bpair_truth_count_hist.tsv",
        truth_hist_rows,
        [
            "scope",
            "truth_n_bucket",
            "n_bpairs",
            "contact_weight",
            "dominant_contact_weight",
            "dominant_contact_fraction_weighted",
            "dominant_fraction_mean_per_bpair",
        ],
    )
    write_dict_rows(
        args.outdir / "split_half_truth_reproducibility.tsv",
        split_rows,
        [
            "scope",
            "truth_n_bucket",
            "n_bpairs_with_both_halves",
            "contact_weight",
            "dominant_state_agreement_contact_weighted",
            "same_cross_agreement_contact_weighted",
            "eval_only_uses_phase_labels",
        ],
    )


if __name__ == "__main__":
    main()
