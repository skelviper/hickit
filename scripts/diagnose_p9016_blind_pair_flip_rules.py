#!/usr/bin/env python3
"""Compare blind chromosome-pair flip rules against eval-only trans oracle.

The flip rules labeled ``blind_*`` choose pair-specific copy flips without SNP
truth or CHARM/3DG. SNP labels are used only afterward to score the chosen rules.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from collections import defaultdict
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


def read_selected_swaps(path: Path) -> dict[str, int]:
    swaps: dict[str, int] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"chrom", "copy_swap", "selected_for_eval"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} missing required columns {sorted(required)}")
        for row in reader:
            if int(row["selected_for_eval"]) == 1:
                swaps[row["chrom"]] = int(row["copy_swap"])
    return swaps


def format_value(value: object) -> str:
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
            writer.writerow({key: format_value(row.get(key)) for key in fieldnames})


def canonical_pair(eval_mod: Any, chrom1: str, chrom2: str) -> tuple[str, str]:
    return tuple(sorted((chrom1, chrom2), key=eval_mod.chrom_sort_key))  # type: ignore[return-value]


def rel_flip_score_for_state(state: int, rel_flip: int) -> int:
    # A relative flip on endpoint 2 maps model state to the alternate pair gauge.
    return state ^ rel_flip


def load_trans_records(
    eval_mod: Any,
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    contact_counts: dict[tuple[str, int, str, int], np.ndarray],
    coords: dict[tuple[str, int, int], np.ndarray],
    reference: dict[tuple[str, int, int], np.ndarray],
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for key, counts in contact_counts.items():
        item = posterior.get(key)
        if item is None:
            continue
        if int(counts.sum()) <= 0:
            continue
        # Match the standard contact_accuracy.tsv denominator.
        if eval_mod.p4_from_coords_for_bpair(reference, item) is None:
            continue
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        if chrom1 == chrom2:
            continue
        distances = eval_mod.all_copy_pair_distances(coords, chrom1, int(item["start1"]), chrom2, int(item["start2"]))
        if distances is None:
            continue
        p4 = np.asarray(item["p4"], dtype=float)
        records.append(
            {
                "pair": canonical_pair(eval_mod, chrom1, chrom2),
                "chrom1": chrom1,
                "chrom2": chrom2,
                "counts": np.asarray(counts, dtype=np.int64),
                "p4": p4,
                "pmax": float(np.max(p4)),
                "margin": float(item.get("margin", 0.0)),
                "distances": np.asarray(distances, dtype=float),
                "posterior_top": int(np.argmax(p4)),
                "nearest_state": int(np.argmin(np.asarray(distances, dtype=float))),
                "n_contacts": int(counts.sum()),
            }
        )
    return records


def choose_pair_oracle(records: list[dict[str, object]], cis_swaps: dict[str, int], eval_mod: Any) -> dict[tuple[str, str], int]:
    scores: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0])
    for rec in records:
        chrom1 = str(rec["chrom1"])
        chrom2 = str(rec["chrom2"])
        counts = rec["counts"]
        p4 = rec["p4"]
        pair = rec["pair"]
        swap1 = cis_swaps.get(chrom1, 0)
        swap2 = cis_swaps.get(chrom2, 0)
        for rel in (0, 1):
            aligned = eval_mod.align_p4_to_truth_gauge(p4, swap1, swap2 ^ rel)
            pred = int(np.argmax(aligned))
            scores[pair][rel] += int(counts[pred])
    return {pair: 1 if vals[1] > vals[0] else 0 for pair, vals in scores.items()}


def choose_blind_rules(records: list[dict[str, object]]) -> dict[str, dict[tuple[str, str], int]]:
    pair_scores: dict[str, dict[tuple[str, str], list[float]]] = {
        "blind_geometry_nearest_pair_flip": defaultdict(lambda: [0.0, 0.0]),
        "blind_posterior_geometry_agreement_pair_flip": defaultdict(lambda: [0.0, 0.0]),
        "blind_high_pmax_geometry_pair_flip": defaultdict(lambda: [0.0, 0.0]),
        "blind_margin_geometry_pair_flip": defaultdict(lambda: [0.0, 0.0]),
    }
    for rec in records:
        pair = rec["pair"]
        n = float(rec["n_contacts"])
        p4 = np.asarray(rec["p4"], dtype=float)
        distances = np.asarray(rec["distances"], dtype=float)
        posterior_top = int(rec["posterior_top"])
        nearest = int(rec["nearest_state"])
        pmax = float(rec["pmax"])
        margin = max(float(rec["margin"]), 0.0)
        # Rule family: choose the relative flip that makes posterior top1 align
        # with the current reconstruction's nearest copy-pair state.
        rel = (posterior_top ^ nearest) & 1
        gap = float(np.partition(distances, 1)[1] - np.partition(distances, 0)[0])
        pair_scores["blind_geometry_nearest_pair_flip"][pair][rel] += n
        pair_scores["blind_posterior_geometry_agreement_pair_flip"][pair][rel] += n * max(gap, 0.0)
        if pmax >= 0.9:
            pair_scores["blind_high_pmax_geometry_pair_flip"][pair][rel] += n * max(gap, 0.0)
        pair_scores["blind_margin_geometry_pair_flip"][pair][rel] += n * margin * max(gap, 0.0)
    out: dict[str, dict[tuple[str, str], int]] = {}
    for policy, scores in pair_scores.items():
        out[policy] = {pair: 1 if vals[1] > vals[0] else 0 for pair, vals in scores.items()}
    return out


def _dist(a: np.ndarray, b: np.ndarray) -> float:
    d = a - b
    return float(np.sqrt(float(np.dot(d, d))))


def choose_center_blind_rules(
    eval_mod: Any,
    records: list[dict[str, object]],
    coords: dict[tuple[str, int, int], np.ndarray],
) -> dict[str, dict[tuple[str, str], int]]:
    by_chrom_copy: dict[tuple[str, int], list[np.ndarray]] = defaultdict(list)
    for (chrom, _start, copy), xyz in coords.items():
        by_chrom_copy[(str(chrom), int(copy))].append(np.asarray(xyz, dtype=float))
    centers: dict[tuple[str, int], np.ndarray] = {}
    for key, vals in by_chrom_copy.items():
        if vals:
            centers[key] = np.vstack(vals).mean(axis=0)

    pair_set = sorted({rec["pair"] for rec in records})
    same_center: dict[tuple[str, str], int] = {}
    vector_dot: dict[tuple[str, str], int] = {}
    for pair in pair_set:
        chrom1, chrom2 = pair
        c10 = centers.get((chrom1, 0))
        c11 = centers.get((chrom1, 1))
        c20 = centers.get((chrom2, 0))
        c21 = centers.get((chrom2, 1))
        if c10 is None or c11 is None or c20 is None or c21 is None:
            same_center[pair] = 0
            vector_dot[pair] = 0
            continue
        same = _dist(c10, c20) + _dist(c11, c21)
        cross = _dist(c10, c21) + _dist(c11, c20)
        same_center[pair] = 1 if cross < same else 0
        v1 = c11 - c10
        v2 = c21 - c20
        vector_dot[pair] = 1 if float(np.dot(v1, v2)) < 0.0 else 0

    distance_scores: dict[str, dict[tuple[str, str], list[float]]] = {
        "blind_contact_distance_same_closer_pair_flip": defaultdict(lambda: [0.0, 0.0]),
        "blind_contact_posterior_same_mass_pair_flip": defaultdict(lambda: [0.0, 0.0]),
    }
    for rec in records:
        pair = rec["pair"]
        n = float(rec["n_contacts"])
        distances = np.asarray(rec["distances"], dtype=float)
        p4 = np.asarray(rec["p4"], dtype=float)
        same_dist = float(distances[0] + distances[3])
        cross_dist = float(distances[1] + distances[2])
        distance_scores["blind_contact_distance_same_closer_pair_flip"][pair][0 if same_dist <= cross_dist else 1] += n * abs(cross_dist - same_dist)
        same_mass = float(p4[0] + p4[3])
        cross_mass = float(p4[1] + p4[2])
        distance_scores["blind_contact_posterior_same_mass_pair_flip"][pair][0 if same_mass >= cross_mass else 1] += n * abs(same_mass - cross_mass)

    out = {
        "blind_chrom_center_same_closer_pair_flip": same_center,
        "blind_chrom_homolog_vector_dot_pair_flip": vector_dot,
    }
    for policy, scores in distance_scores.items():
        out[policy] = {pair: 1 if vals[1] > vals[0] else 0 for pair, vals in scores.items()}
    return out


def evaluate_policy(
    eval_mod: Any,
    records: list[dict[str, object]],
    cis_swaps: dict[str, int],
    pair_flips: dict[tuple[str, str], int],
    policy: str,
    flip_source: str,
) -> dict[str, object]:
    stats = eval_mod.empty_accuracy_stats()
    for rec in records:
        chrom1 = str(rec["chrom1"])
        chrom2 = str(rec["chrom2"])
        p4 = rec["p4"]
        counts = rec["counts"]
        rel = int(pair_flips.get(rec["pair"], 0))
        swap1 = cis_swaps.get(chrom1, 0)
        swap2 = cis_swaps.get(chrom2, 0) ^ rel
        aligned = eval_mod.align_p4_to_truth_gauge(p4, swap1, swap2)
        eval_mod.update_accuracy_stats(stats, counts, aligned)
    row = eval_mod.finalize_accuracy_stats(stats)
    row.update(
        {
            "policy": policy,
            "scope": "genome_trans",
            "flip_source": flip_source,
            "eval_only_uses_snp_for_flip_selection": 1 if flip_source == "snp_truth_oracle" else 0,
            "eval_only_uses_snp_for_scoring": 1,
            "eval_only_uses_charm_shared_denominator": 1,
        }
    )
    return row


def pair_rows_for_policies(
    eval_mod: Any,
    records: list[dict[str, object]],
    cis_swaps: dict[str, int],
    policies: dict[str, tuple[dict[tuple[str, str], int], str]],
) -> list[dict[str, object]]:
    by_pair: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for rec in records:
        by_pair[rec["pair"]].append(rec)
    rows: list[dict[str, object]] = []
    for pair, pair_records in sorted(by_pair.items(), key=lambda item: (eval_mod.chrom_sort_key(item[0][0]), eval_mod.chrom_sort_key(item[0][1]))):
        for policy, (flips, source) in policies.items():
            stats = eval_mod.empty_accuracy_stats()
            rel = int(flips.get(pair, 0))
            for rec in pair_records:
                chrom1 = str(rec["chrom1"])
                chrom2 = str(rec["chrom2"])
                swap1 = cis_swaps.get(chrom1, 0)
                swap2 = cis_swaps.get(chrom2, 0) ^ rel
                aligned = eval_mod.align_p4_to_truth_gauge(rec["p4"], swap1, swap2)
                eval_mod.update_accuracy_stats(stats, rec["counts"], aligned)
            final = eval_mod.finalize_accuracy_stats(stats)
            rows.append(
                {
                    "chr1": pair[0],
                    "chr2": pair[1],
                    "policy": policy,
                    "relative_flip": rel,
                    "flip_source": source,
                    "n_eval_contacts": final["n_eval_contacts"],
                    "top1_accuracy": final["top1_accuracy"],
                    "same_cross_accuracy": final["same_cross_accuracy"],
                    "pmax_threshold_accuracy": final["pmax_threshold_accuracy"],
                    "pmax_threshold_recall": final["pmax_threshold_recall"],
                }
            )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config-output-dir", type=Path, required=True)
    parser.add_argument("--eval-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, default=None)
    parser.add_argument("--reference-3dg", type=Path, default=None)
    parser.add_argument("--bin-size", type=int, default=None)
    args = parser.parse_args()

    eval_mod = load_eval_module()
    summary = read_kv(args.eval_dir / "summary.tsv")
    pairs = args.pairs or Path(summary["pairs_path"])
    reference_path = args.reference_3dg or Path(summary["reference_3dg_path"])
    bin_size = args.bin_size or int(summary["bin_size_bp"])
    posterior = eval_mod.read_posterior(args.config_output_dir / "p9016_full.bpair_posterior.tsv")
    coords = eval_mod.read_reconstruction(args.config_output_dir / "p9016_full.coords.tsv")
    reference = eval_mod.read_reference_3dg(reference_path, bin_size)
    contact_counts = eval_mod.read_contact_truth_counts(pairs, bin_size, posterior)
    cis_swaps = read_selected_swaps(args.eval_dir / "whole_chrom_snp_oracle_swaps.tsv")
    records = load_trans_records(eval_mod, posterior, contact_counts, coords, reference)
    oracle_flips = choose_pair_oracle(records, cis_swaps, eval_mod)
    blind_rules = choose_blind_rules(records)
    policies: dict[str, tuple[dict[tuple[str, str], int], str]] = {
        "cis_selected_whole_chrom": ({}, "cis_snp_truth_whole_chrom"),
        "pair_independent_trans_oracle": (oracle_flips, "snp_truth_oracle"),
    }
    for policy, flips in blind_rules.items():
        policies[policy] = (flips, "blind_posterior_reconstruction_geometry")
    for policy, flips in choose_center_blind_rules(eval_mod, records, coords).items():
        policies[policy] = (flips, "blind_reconstruction_geometry")
    summary_rows = [
        evaluate_policy(eval_mod, records, cis_swaps, flips, policy, source)
        for policy, (flips, source) in policies.items()
    ]
    pair_rows = pair_rows_for_policies(eval_mod, records, cis_swaps, policies)
    summary_fields = [
        "policy",
        "scope",
        "flip_source",
        "n_eval_contacts",
        "top1_correct_contacts",
        "top1_accuracy",
        "same_cross_correct_contacts",
        "same_cross_accuracy",
        "pmax_threshold",
        "n_called_contacts",
        "called_contact_fraction",
        "pmax_threshold_accuracy",
        "pmax_threshold_recall",
        "eval_only_uses_snp_for_flip_selection",
        "eval_only_uses_snp_for_scoring",
        "eval_only_uses_charm_shared_denominator",
    ]
    pair_fields = [
        "chr1",
        "chr2",
        "policy",
        "relative_flip",
        "flip_source",
        "n_eval_contacts",
        "top1_accuracy",
        "same_cross_accuracy",
        "pmax_threshold_accuracy",
        "pmax_threshold_recall",
    ]
    write_rows(args.outdir / "blind_pair_flip_rule_summary.tsv", summary_rows, summary_fields)
    write_rows(args.outdir / "blind_pair_flip_rule_chr_pair.tsv", pair_rows, pair_fields)
    print(f"wrote\t{args.outdir / 'blind_pair_flip_rule_summary.tsv'}")
    print(f"wrote\t{args.outdir / 'blind_pair_flip_rule_chr_pair.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
