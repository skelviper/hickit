#!/usr/bin/env python3
"""Evaluate blind callable-trans confidence scores for a finished P9016 run.

All call scores are computed without SNP truth or CHARM/3DG. SNP labels are
used only after scoring to compute contact accuracy under the standard
whole-chromosome cis SNP gauge. CHARM/3DG is used only to keep the same shared
denominator as the standard evaluator.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from pathlib import Path
from typing import Any

import numpy as np


def load_eval_module() -> Any:
    repo_root = Path(__file__).resolve().parents[1]
    eval_path = repo_root / "eval" / "evaluate_p9016_baseline.py"
    spec = importlib.util.spec_from_file_location("evaluate_p9016_baseline", eval_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load eval module from {eval_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_kv(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header != ["key", "value"]:
            raise ValueError(f"{path} is not a key/value TSV")
        for row in reader:
            if len(row) >= 2:
                values[row[0]] = row[1]
    return values


def read_summary_value(path: Path, key: str) -> str:
    if not path.exists():
        return "NA"
    return read_kv(path).get(key, "NA")


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.9g}"
    return str(value)


def write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in fieldnames})


def entropy_score(p4: np.ndarray) -> float:
    vals = p4[p4 > 0]
    if vals.size == 0:
        return float("-inf")
    return -float(-(vals * np.log(vals)).sum())


def margin_score(p4: np.ndarray) -> float:
    ordered = np.sort(p4)
    return float(ordered[-1] - ordered[-2])


def distance_gap_score(distances: np.ndarray) -> float:
    ordered = np.sort(distances)
    if ordered.size < 2:
        return float("-inf")
    return float(ordered[1] - ordered[0])


def canonical_pair(eval_mod: Any, chrom1: str, chrom2: str) -> tuple[str, str]:
    return tuple(sorted((chrom1, chrom2), key=eval_mod.chrom_sort_key))  # type: ignore[return-value]


def collect_trans_records(
    eval_mod: Any,
    contact_counts: dict[tuple[str, int, str, int], np.ndarray],
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    reconstruction: dict[tuple[str, int, int], np.ndarray],
    reference: dict[tuple[str, int, int], np.ndarray],
    swaps: dict[str, int],
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for key, counts in contact_counts.items():
        item = posterior.get(key)
        if item is None:
            continue
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        if chrom1 == chrom2:
            continue
        n = int(np.asarray(counts, dtype=np.int64).sum())
        if n <= 0:
            continue
        # Match standard contact_accuracy.tsv denominator.
        if eval_mod.p4_from_coords_for_bpair(reference, item) is None:
            continue
        distances = eval_mod.all_copy_pair_distances(reconstruction, chrom1, int(item["start1"]), chrom2, int(item["start2"]))
        if distances is None:
            continue
        p4 = np.asarray(item["p4"], dtype=float)
        aligned = eval_mod.align_p4_to_truth_gauge(p4, int(swaps.get(chrom1, 0)), int(swaps.get(chrom2, 0)))
        pred = int(np.argmax(aligned))
        correct = int(counts[pred])
        pred_same = pred in (0, 3)
        same_counts = int(counts[0] + counts[3])
        cross_counts = int(counts[1] + counts[2])
        same_cross_correct = same_counts if pred_same else cross_counts
        nearest_state = int(np.argmin(np.asarray(distances, dtype=float)))
        posterior_top = int(np.argmax(p4))
        scores = {
            "pmax": float(np.max(p4)),
            "margin": margin_score(p4),
            "neg_entropy": entropy_score(p4),
            "distance_gap": distance_gap_score(np.asarray(distances, dtype=float)),
            "posterior_nearest_agreement_gap": distance_gap_score(np.asarray(distances, dtype=float))
            if posterior_top == nearest_state
            else -distance_gap_score(np.asarray(distances, dtype=float)),
            "n_raw": float(item.get("n_raw", n)),
        }
        records.append(
            {
                "chrom_pair": "__".join(canonical_pair(eval_mod, chrom1, chrom2)),
                "chrom1": chrom1,
                "start1": int(item["start1"]),
                "chrom2": chrom2,
                "start2": int(item["start2"]),
                "n_contacts": n,
                "top1_correct": correct,
                "same_cross_correct": same_cross_correct,
                "scores": scores,
            }
        )
    return records


def curve_for_score(
    records: list[dict[str, object]],
    score_name: str,
    fractions: list[float],
) -> list[dict[str, object]]:
    total_contacts = sum(int(r["n_contacts"]) for r in records)
    total_correct = sum(int(r["top1_correct"]) for r in records)
    total_same_cross_correct = sum(int(r["same_cross_correct"]) for r in records)
    ordered = sorted(records, key=lambda r: float(r["scores"][score_name]), reverse=True)
    rows: list[dict[str, object]] = []
    cumulative_n = 0
    cumulative_correct = 0
    cumulative_same_cross = 0
    idx = 0
    for target_fraction in fractions:
        target = math.ceil(total_contacts * target_fraction)
        while idx < len(ordered) and cumulative_n < target:
            rec = ordered[idx]
            cumulative_n += int(rec["n_contacts"])
            cumulative_correct += int(rec["top1_correct"])
            cumulative_same_cross += int(rec["same_cross_correct"])
            idx += 1
        threshold = float(ordered[idx - 1]["scores"][score_name]) if idx > 0 else float("nan")
        rows.append(
            {
                "score_name": score_name,
                "target_called_fraction": target_fraction,
                "actual_called_contacts": cumulative_n,
                "actual_called_fraction": cumulative_n / total_contacts if total_contacts else float("nan"),
                "score_threshold": threshold,
                "top1_accuracy": cumulative_correct / cumulative_n if cumulative_n else float("nan"),
                "top1_recall": cumulative_correct / total_contacts if total_contacts else float("nan"),
                "same_cross_accuracy": cumulative_same_cross / cumulative_n if cumulative_n else float("nan"),
                "same_cross_recall": cumulative_same_cross / total_contacts if total_contacts else float("nan"),
                "baseline_top1_accuracy": total_correct / total_contacts if total_contacts else float("nan"),
                "baseline_same_cross_accuracy": total_same_cross_correct / total_contacts if total_contacts else float("nan"),
                "delta_top1_vs_all_trans": (cumulative_correct / cumulative_n - total_correct / total_contacts)
                if cumulative_n and total_contacts
                else float("nan"),
                "eval_only_uses_phase_labels_for_scoring": 1,
                "call_score_uses_phase_labels": 0,
                "call_score_uses_charm_or_reference": 0,
            }
        )
    return rows


def chrom_pair_curve(records: list[dict[str, object]], score_name: str) -> list[dict[str, object]]:
    by_pair: dict[str, list[dict[str, object]]] = {}
    for rec in records:
        by_pair.setdefault(str(rec["chrom_pair"]), []).append(rec)
    rows: list[dict[str, object]] = []
    for pair, pair_records in sorted(by_pair.items()):
        n = sum(int(r["n_contacts"]) for r in pair_records)
        correct = sum(int(r["top1_correct"]) for r in pair_records)
        same = sum(int(r["same_cross_correct"]) for r in pair_records)
        mean_score = sum(float(r["scores"][score_name]) * int(r["n_contacts"]) for r in pair_records) / n if n else float("nan")
        rows.append(
            {
                "chrom_pair": pair,
                "score_name": score_name,
                "n_contacts": n,
                "top1_accuracy": correct / n if n else float("nan"),
                "same_cross_accuracy": same / n if n else float("nan"),
                "contact_weighted_mean_score": mean_score,
                "eval_only_uses_phase_labels_for_scoring": 1,
                "call_score_uses_phase_labels": 0,
            }
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config-output-dir", required=True, type=Path)
    ap.add_argument("--eval-output-dir", required=True, type=Path)
    ap.add_argument("--pairs", required=True, type=Path)
    ap.add_argument("--reference-3dg", required=True, type=Path)
    ap.add_argument("--bin-size", type=int, default=1_000_000)
    ap.add_argument("--label", default=None)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()

    eval_mod = load_eval_module()
    out_dir = args.outdir
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = args.config_output_dir / "p9016_full.manifest.tsv"
    coords_path = args.config_output_dir / "p9016_full.coords.tsv"
    posterior_path = args.config_output_dir / "p9016_full.bpair_posterior.tsv"
    swaps_path = args.eval_output_dir / "whole_chrom_snp_oracle_swaps.tsv"
    if not swaps_path.exists():
        # Backward-compatible name from older eval outputs.
        swaps_path = args.eval_output_dir / "model_whole_chrom_snp_oracle_swaps.tsv"
    if not swaps_path.exists():
        raise FileNotFoundError(f"missing whole-chrom SNP swaps in {args.eval_output_dir}")

    manifest = read_kv(manifest_path)
    posterior = eval_mod.read_posterior(posterior_path)
    reconstruction = eval_mod.read_reconstruction(coords_path)
    reference = eval_mod.read_reference_3dg(args.reference_3dg, args.bin_size)
    contact_counts = eval_mod.read_contact_truth_counts(args.pairs, args.bin_size, posterior)
    swaps = {}
    with swaps_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if int(row.get("selected_for_eval", "0")) == 1:
                swaps[row["chrom"]] = int(row["copy_swap"])

    records = collect_trans_records(eval_mod, contact_counts, posterior, reconstruction, reference, swaps)
    fractions = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]
    score_names = ["pmax", "margin", "neg_entropy", "distance_gap", "posterior_nearest_agreement_gap", "n_raw"]
    curve_rows: list[dict[str, object]] = []
    for score_name in score_names:
        curve_rows.extend(curve_for_score(records, score_name, fractions))
    curve_fields = [
        "config_name",
        "score_name",
        "target_called_fraction",
        "actual_called_contacts",
        "actual_called_fraction",
        "score_threshold",
        "top1_accuracy",
        "top1_recall",
        "same_cross_accuracy",
        "same_cross_recall",
        "baseline_top1_accuracy",
        "baseline_same_cross_accuracy",
        "delta_top1_vs_all_trans",
        "eval_only_uses_phase_labels_for_scoring",
        "call_score_uses_phase_labels",
        "call_score_uses_charm_or_reference",
    ]
    config_name = args.label or manifest.get("config_name", args.config_output_dir.name)
    for row in curve_rows:
        row["config_name"] = config_name
    write_rows(out_dir / "callable_trans_curve.tsv", curve_rows, curve_fields)

    record_rows: list[dict[str, object]] = []
    for rec in records:
        row = {
            "config_name": config_name,
            "chrom_pair": rec["chrom_pair"],
            "chrom1": rec["chrom1"],
            "start1": rec["start1"],
            "chrom2": rec["chrom2"],
            "start2": rec["start2"],
            "n_contacts": rec["n_contacts"],
            "top1_correct": rec["top1_correct"],
            "same_cross_correct": rec["same_cross_correct"],
            "top1_is_correct": int(int(rec["top1_correct"]) > 0),
        }
        for score_name, score in rec["scores"].items():  # type: ignore[union-attr]
            row[score_name] = score
        record_rows.append(row)
    write_rows(
        out_dir / "callable_trans_records.tsv",
        record_rows,
        [
            "config_name",
            "chrom_pair",
            "chrom1",
            "start1",
            "chrom2",
            "start2",
            "n_contacts",
            "top1_correct",
            "same_cross_correct",
            "top1_is_correct",
            "pmax",
            "margin",
            "neg_entropy",
            "distance_gap",
            "posterior_nearest_agreement_gap",
            "n_raw",
        ],
    )

    best_rows: list[dict[str, object]] = []
    for score_name in score_names:
        candidates = [r for r in curve_rows if r["score_name"] == score_name and float(r["actual_called_fraction"]) <= 0.5]
        if not candidates:
            continue
        best = max(candidates, key=lambda r: (float(r["top1_accuracy"]), float(r["actual_called_fraction"])))
        best_rows.append(dict(best))
    write_rows(out_dir / "callable_trans_best.tsv", best_rows, curve_fields)

    pair_rows = chrom_pair_curve(records, "posterior_nearest_agreement_gap")
    for row in pair_rows:
        row["config_name"] = config_name
    write_rows(
        out_dir / "callable_trans_chrom_pair.tsv",
        pair_rows,
        [
            "config_name",
            "chrom_pair",
            "score_name",
            "n_contacts",
            "top1_accuracy",
            "same_cross_accuracy",
            "contact_weighted_mean_score",
            "eval_only_uses_phase_labels_for_scoring",
            "call_score_uses_phase_labels",
        ],
    )

    summary_rows = [
        ("config_name", config_name),
        ("n_trans_records", len(records)),
        ("n_trans_contacts", sum(int(r["n_contacts"]) for r in records)),
        ("baseline_top1_accuracy", curve_rows[0]["baseline_top1_accuracy"] if curve_rows else "NA"),
        ("baseline_same_cross_accuracy", curve_rows[0]["baseline_same_cross_accuracy"] if curve_rows else "NA"),
        ("best_score_le_50pct_coverage", best_rows[0]["score_name"] if best_rows else "NA"),
        ("best_top1_accuracy_le_50pct_coverage", max((float(r["top1_accuracy"]) for r in best_rows), default=float("nan"))),
        ("eval_only_uses_phase_labels_for_scoring", 1),
        ("call_score_uses_phase_labels", 0),
        ("call_score_uses_charm_or_reference", 0),
        ("selection_uses_truth", 0),
        ("selection_uses_charm", 0),
        ("scoring_uses_truth", 1),
        ("denominator_uses_charm_filter", 1),
        ("train_manifest_input_contact_source", manifest.get("input_contact_source", "NA")),
        ("train_manifest_uses_phase_labels", manifest.get("uses_phase_labels", "NA")),
        ("train_manifest_uses_charm_or_reference", manifest.get("uses_charm_or_reference", "NA")),
        ("train_manifest_uses_charm_for_training", manifest.get("uses_charm_for_training", "NA")),
        ("train_manifest_reference_derived_positive_control", manifest.get("reference_derived_positive_control", "NA")),
        ("train_manifest_approved_p9016_raw_pairs_realpath", manifest.get("approved_p9016_raw_pairs_realpath", "NA")),
    ]
    with (out_dir / "summary.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        for key, value in summary_rows:
            writer.writerow([key, fmt(value)])


if __name__ == "__main__":
    main()
