#!/usr/bin/env python3
"""Audit raw P9016 pairs-only identifiability for trans copy-state calls.

The scoring features in this script are split into:
- pairs-visible metadata: columns physically present in the P9016 pairs file;
- model-visible diagnostics: posterior/confidence values from a finished blind
  reconstruction trained on raw pairs;
- eval-only labels: SNP phase and CHARM/3DG reference, used only after scoring.

No output from this script is a training input.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np


STATE_COLUMNS = ("p00", "p01", "p10", "p11")


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def load_eval_module() -> Any:
    repo_root = Path(__file__).resolve().parents[1]
    eval_path = repo_root / "eval" / "evaluate_p9016_baseline.py"
    spec = importlib.util.spec_from_file_location("evaluate_p9016_baseline", eval_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load eval module from {eval_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.9g}"
    return str(value)


def write_rows(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field)) for field in fields})


def read_kv(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header != ["key", "value"]:
            raise ValueError(f"{path} is not a key/value TSV")
        for row in reader:
            if len(row) >= 2:
                out[row[0]] = row[1]
    return out


def entropy(p4: np.ndarray) -> float:
    vals = p4[p4 > 0]
    if vals.size == 0:
        return float("nan")
    return float(-(vals * np.log(vals)).sum())


def margin(p4: np.ndarray) -> float:
    vals = np.sort(p4)
    return float(vals[-1] - vals[-2])


def distance_gap(distances: np.ndarray | None) -> float:
    if distances is None or len(distances) < 2:
        return float("nan")
    vals = np.sort(np.asarray(distances, dtype=float))
    return float(vals[1] - vals[0])


def canonical_bpair(eval_mod: Any, chrom1: str, start1: int, chrom2: str, start2: int) -> tuple[str, int, str, int]:
    return eval_mod.canonical_bpair_key(chrom1, start1, chrom2, start2)


def inspect_pairs(path: Path, bin_size: int) -> tuple[list[dict[str, object]], dict[tuple[str, int, str, int], dict[str, object]]]:
    n_total = n_cis = n_trans = 0
    read_ids: Counter[str] = Counter()
    strand_pairs: Counter[tuple[str, str]] = Counter()
    phase_pairs: Counter[tuple[str, str]] = Counter()
    n_with_phase0 = n_with_phase1 = n_with_both_phase = n_with_any_phase = 0
    bpair_meta: dict[tuple[str, int, str, int], dict[str, object]] = {}
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            read_id, chrom1, pos1_s, chrom2, pos2_s, strand1, strand2, phase0, phase1 = fields[:9]
            pos1 = int(pos1_s)
            pos2 = int(pos2_s)
            start1 = (pos1 // bin_size) * bin_size
            start2 = (pos2 // bin_size) * bin_size
            key = (chrom1, start1, chrom2, start2)
            if (chrom2, start2, chrom1, start1) < key:
                key = (chrom2, start2, chrom1, start1)
            n_total += 1
            is_trans = chrom1 != chrom2
            n_trans += int(is_trans)
            n_cis += int(not is_trans)
            read_ids[read_id] += 1
            strand_pairs[(strand1, strand2)] += 1
            phase_pairs[(phase0, phase1)] += 1
            has0 = phase0 in ("0", "1")
            has1 = phase1 in ("0", "1")
            n_with_phase0 += int(has0)
            n_with_phase1 += int(has1)
            n_with_both_phase += int(has0 and has1)
            n_with_any_phase += int(has0 or has1)
            meta = bpair_meta.setdefault(
                key,
                {
                    "n_raw_pairs": 0,
                    "n_phase_both": 0,
                    "n_phase_any": 0,
                    "n_readid_non_dot": 0,
                    "read_ids": Counter(),
                    "strand_pairs": Counter(),
                },
            )
            meta["n_raw_pairs"] = int(meta["n_raw_pairs"]) + 1
            meta["n_phase_both"] = int(meta["n_phase_both"]) + int(has0 and has1)
            meta["n_phase_any"] = int(meta["n_phase_any"]) + int(has0 or has1)
            meta["n_readid_non_dot"] = int(meta["n_readid_non_dot"]) + int(read_id != ".")
            meta["read_ids"][read_id] += 1  # type: ignore[index]
            meta["strand_pairs"][(strand1, strand2)] += 1  # type: ignore[index]
    rows = [
        {"field": "n_total_pairs", "value": n_total},
        {"field": "n_cis_pairs", "value": n_cis},
        {"field": "n_trans_pairs", "value": n_trans},
        {"field": "n_unique_read_ids", "value": len(read_ids)},
        {"field": "most_common_read_id", "value": read_ids.most_common(1)[0][0] if read_ids else "NA"},
        {"field": "most_common_read_id_count", "value": read_ids.most_common(1)[0][1] if read_ids else "NA"},
        {"field": "n_non_dot_read_id_pairs", "value": n_total - read_ids.get(".", 0)},
        {"field": "n_unique_strand_pairs", "value": len(strand_pairs)},
        {"field": "strand_pair_distribution", "value": ";".join(f"{a}/{b}:{n}" for (a, b), n in strand_pairs.most_common())},
        {"field": "n_phase0_labeled", "value": n_with_phase0},
        {"field": "n_phase1_labeled", "value": n_with_phase1},
        {"field": "n_both_phase_labeled", "value": n_with_both_phase},
        {"field": "n_any_phase_labeled", "value": n_with_any_phase},
        {"field": "phase_pair_distribution", "value": ";".join(f"{a}/{b}:{n}" for (a, b), n in phase_pairs.most_common())},
        {"field": "pairs_visible_read_linkage_available", "value": int(len(read_ids) > 1 or read_ids.get(".", 0) != n_total)},
        {"field": "pairs_visible_strand_information_available", "value": int(len(strand_pairs) > 1)},
        {"field": "phase_columns_present_eval_only", "value": 1},
    ]
    return rows, bpair_meta


def load_swaps(path: Path) -> dict[str, int]:
    swaps: dict[str, int] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if int(row.get("selected_for_eval", "0")) == 1:
                swaps[row["chrom"]] = int(row["copy_swap"])
    return swaps


def build_records(
    eval_mod: Any,
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    contact_counts: dict[tuple[str, int, str, int], np.ndarray],
    reconstruction: dict[tuple[str, int, int], np.ndarray],
    reference: dict[tuple[str, int, int], np.ndarray],
    swaps: dict[str, int],
    bpair_meta: dict[tuple[str, int, str, int], dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for key, counts_raw in contact_counts.items():
        item = posterior.get(key)
        if item is None:
            continue
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        if chrom1 == chrom2:
            continue
        if eval_mod.p4_from_coords_for_bpair(reference, item) is None:
            continue
        n_contacts = int(np.asarray(counts_raw, dtype=np.int64).sum())
        if n_contacts <= 0:
            continue
        p4_raw = np.asarray(item["p4"], dtype=float)
        p4 = eval_mod.align_p4_to_truth_gauge(p4_raw, int(swaps.get(chrom1, 0)), int(swaps.get(chrom2, 0)))
        pred = int(np.argmax(p4))
        top1_correct = int(counts_raw[pred])
        pred_same = pred in (0, 3)
        same_cross_correct = int(counts_raw[0] + counts_raw[3]) if pred_same else int(counts_raw[1] + counts_raw[2])
        distances = eval_mod.all_copy_pair_distances(reconstruction, chrom1, int(item["start1"]), chrom2, int(item["start2"]))
        dgap = distance_gap(distances)
        meta_key = canonical_bpair(eval_mod, chrom1, int(item["start1"]), chrom2, int(item["start2"]))
        meta = bpair_meta.get(meta_key, {})
        rows.append(
            {
                "chrom_pair": "__".join(sorted((chrom1, chrom2), key=eval_mod.chrom_sort_key)),
                "chrom1": chrom1,
                "start1": int(item["start1"]),
                "chrom2": chrom2,
                "start2": int(item["start2"]),
                "n_eval_contacts": n_contacts,
                "n_raw_pairs": int(item.get("n_raw", meta.get("n_raw_pairs", n_contacts))),
                "n_phase_any_pairs": int(meta.get("n_phase_any", 0)),
                "n_phase_both_pairs": int(meta.get("n_phase_both", 0)),
                "n_readid_non_dot_pairs": int(meta.get("n_readid_non_dot", 0)),
                "pmax": float(np.max(p4_raw)),
                "margin": margin(p4_raw),
                "neg_entropy": -entropy(p4_raw),
                "distance_gap": dgap,
                "top1_correct_contacts": top1_correct,
                "same_cross_correct_contacts": same_cross_correct,
                "top1_accuracy_contact_weighted": top1_correct / n_contacts,
                "same_cross_accuracy_contact_weighted": same_cross_correct / n_contacts,
            }
        )
    return rows


def curve(records: list[dict[str, object]], score: str, fractions: list[float]) -> list[dict[str, object]]:
    total = sum(int(r["n_eval_contacts"]) for r in records)
    total_top1 = sum(int(r["top1_correct_contacts"]) for r in records)
    total_same = sum(int(r["same_cross_correct_contacts"]) for r in records)
    ordered = sorted(records, key=lambda r: float(r[score]), reverse=True)
    out: list[dict[str, object]] = []
    n = top1 = same = idx = 0
    for frac in fractions:
        target = math.ceil(total * frac)
        while idx < len(ordered) and n < target:
            rec = ordered[idx]
            n += int(rec["n_eval_contacts"])
            top1 += int(rec["top1_correct_contacts"])
            same += int(rec["same_cross_correct_contacts"])
            idx += 1
        threshold = float(ordered[idx - 1][score]) if idx else float("nan")
        out.append(
            {
                "score": score,
                "target_fraction": frac,
                "called_contacts": n,
                "called_fraction": n / total if total else float("nan"),
                "threshold": threshold,
                "top1_accuracy": top1 / n if n else float("nan"),
                "top1_recall": top1 / total if total else float("nan"),
                "same_cross_accuracy": same / n if n else float("nan"),
                "same_cross_recall": same / total if total else float("nan"),
                "baseline_top1_accuracy": total_top1 / total if total else float("nan"),
                "baseline_same_cross_accuracy": total_same / total if total else float("nan"),
                "delta_top1_vs_baseline": top1 / n - total_top1 / total if n and total else float("nan"),
                "score_uses_phase_labels": int(score.startswith("n_phase")),
                "score_uses_charm_or_reference": int(score == "distance_gap"),
                "eval_uses_phase_labels": 1,
            }
        )
    return out


def raw_count_bins(records: list[dict[str, object]]) -> list[dict[str, object]]:
    bins = [(1, 1), (2, 2), (3, 5), (6, 10), (11, 20), (21, 50), (51, 10**12)]
    out: list[dict[str, object]] = []
    for lo, hi in bins:
        subset = [r for r in records if lo <= int(r["n_raw_pairs"]) <= hi]
        n = sum(int(r["n_eval_contacts"]) for r in subset)
        if not subset:
            out.append({"raw_count_bin": f"{lo}-{hi if hi < 10**12 else 'inf'}", "n_bpairs": 0, "n_eval_contacts": 0})
            continue
        top1 = sum(int(r["top1_correct_contacts"]) for r in subset)
        same = sum(int(r["same_cross_correct_contacts"]) for r in subset)
        out.append(
            {
                "raw_count_bin": f"{lo}-{hi if hi < 10**12 else 'inf'}",
                "n_bpairs": len(subset),
                "n_eval_contacts": n,
                "top1_accuracy": top1 / n if n else float("nan"),
                "same_cross_accuracy": same / n if n else float("nan"),
                "mean_pmax": sum(float(r["pmax"]) * int(r["n_eval_contacts"]) for r in subset) / n if n else float("nan"),
            }
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", required=True, type=Path)
    ap.add_argument("--reference-3dg", required=True, type=Path)
    ap.add_argument("--config-output-dir", required=True, type=Path)
    ap.add_argument("--eval-output-dir", required=True, type=Path)
    ap.add_argument("--bin-size", type=int, default=1_000_000)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()

    eval_mod = load_eval_module()
    args.outdir.mkdir(parents=True, exist_ok=True)

    pair_rows, bpair_meta = inspect_pairs(args.pairs, args.bin_size)
    write_rows(args.outdir / "pairs_visible_fields.tsv", pair_rows, ["field", "value"])

    posterior = eval_mod.read_posterior(args.config_output_dir / "p9016_full.bpair_posterior.tsv")
    reconstruction = eval_mod.read_reconstruction(args.config_output_dir / "p9016_full.coords.tsv")
    reference = eval_mod.read_reference_3dg(args.reference_3dg, args.bin_size)
    contact_counts = eval_mod.read_contact_truth_counts(args.pairs, args.bin_size, posterior)
    swaps = load_swaps(args.eval_output_dir / "whole_chrom_snp_oracle_swaps.tsv")
    records = build_records(eval_mod, posterior, contact_counts, reconstruction, reference, swaps, bpair_meta)

    record_fields = [
        "chrom_pair", "chrom1", "start1", "chrom2", "start2", "n_eval_contacts", "n_raw_pairs",
        "n_phase_any_pairs", "n_phase_both_pairs", "n_readid_non_dot_pairs", "pmax", "margin",
        "neg_entropy", "distance_gap", "top1_correct_contacts", "same_cross_correct_contacts",
        "top1_accuracy_contact_weighted", "same_cross_accuracy_contact_weighted",
    ]
    write_rows(args.outdir / "trans_bpair_records.tsv", records, record_fields)

    fractions = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]
    scores = ["n_raw_pairs", "pmax", "margin", "neg_entropy", "distance_gap", "n_phase_any_pairs", "n_phase_both_pairs"]
    curve_rows: list[dict[str, object]] = []
    for score in scores:
        curve_rows.extend(curve(records, score, fractions))
    curve_fields = [
        "score", "target_fraction", "called_contacts", "called_fraction", "threshold",
        "top1_accuracy", "top1_recall", "same_cross_accuracy", "same_cross_recall",
        "baseline_top1_accuracy", "baseline_same_cross_accuracy", "delta_top1_vs_baseline",
        "score_uses_phase_labels", "score_uses_charm_or_reference", "eval_uses_phase_labels",
    ]
    write_rows(args.outdir / "coverage_accuracy_curve.tsv", curve_rows, curve_fields)
    write_rows(
        args.outdir / "raw_count_bin_accuracy.tsv",
        raw_count_bins(records),
        ["raw_count_bin", "n_bpairs", "n_eval_contacts", "top1_accuracy", "same_cross_accuracy", "mean_pmax"],
    )

    pairs_info = {str(r["field"]): str(r["value"]) for r in pair_rows}
    baseline_top1 = float(curve_rows[0]["baseline_top1_accuracy"]) if curve_rows else float("nan")
    best_blind = max(
        (r for r in curve_rows if int(r["score_uses_phase_labels"]) == 0 and int(r["score_uses_charm_or_reference"]) == 0 and float(r["called_fraction"]) <= 0.5),
        key=lambda r: float(r["top1_accuracy"]),
        default=None,
    )
    rows = [
        {"key": "config_output_dir", "value": args.config_output_dir},
        {"key": "eval_output_dir", "value": args.eval_output_dir},
        {"key": "n_trans_bpairs_with_eval_contacts", "value": len(records)},
        {"key": "n_trans_eval_contacts", "value": sum(int(r["n_eval_contacts"]) for r in records)},
        {"key": "baseline_trans_top1_accuracy", "value": baseline_top1},
        {"key": "pairs_visible_read_linkage_available", "value": pairs_info.get("pairs_visible_read_linkage_available", "NA")},
        {"key": "pairs_visible_strand_information_available", "value": pairs_info.get("pairs_visible_strand_information_available", "NA")},
        {"key": "phase_columns_present_eval_only", "value": pairs_info.get("phase_columns_present_eval_only", "NA")},
        {"key": "best_blind_score_le_50pct", "value": best_blind["score"] if best_blind else "NA"},
        {"key": "best_blind_called_fraction_le_50pct", "value": best_blind["called_fraction"] if best_blind else "NA"},
        {"key": "best_blind_top1_accuracy_le_50pct", "value": best_blind["top1_accuracy"] if best_blind else "NA"},
        {"key": "best_blind_delta_top1_vs_baseline", "value": best_blind["delta_top1_vs_baseline"] if best_blind else "NA"},
        {"key": "best_blind_recall_le_50pct", "value": best_blind["top1_recall"] if best_blind else "NA"},
        {"key": "full_denominator_plus_0p1_supported_by_pairs_visible_fields", "value": 0},
        {"key": "training_uses_phase_labels", "value": 0},
        {"key": "training_uses_charm_or_reference", "value": 0},
    ]
    with (args.outdir / "summary.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=["key", "value"], lineterminator="\n")
        writer.writeheader()
        writer.writerows({"key": r["key"], "value": fmt(r["value"])} for r in rows)


if __name__ == "__main__":
    main()
