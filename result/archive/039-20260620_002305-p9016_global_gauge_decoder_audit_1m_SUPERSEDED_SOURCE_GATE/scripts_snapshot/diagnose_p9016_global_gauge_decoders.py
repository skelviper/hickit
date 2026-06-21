#!/usr/bin/env python3
"""Audit raw-only global chromosome gauge decoders for trans contacts.

This diagnostic does not train a model. It chooses one copy-swap bit per
chromosome from raw posterior/reconstruction-only objectives, then scores the
result with SNP labels as eval-only truth.
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
                swaps[str(row["chrom"])] = int(row["copy_swap"])
    return swaps


def read_posterior_with_extra(eval_mod: Any, path: Path) -> dict[tuple[str, int, str, int], dict[str, object]]:
    posterior = eval_mod.read_posterior(path)
    with eval_mod.open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a header")
        required_extra = {"margin"}
        missing_extra = required_extra.difference(reader.fieldnames)
        if missing_extra:
            raise ValueError(f"{path} is missing required gauge diagnostic columns: {sorted(missing_extra)}")
        numeric_fields = ("pU", "psame_raw", "pcross_raw", "entropy", "margin", "rho_output")
        for row in reader:
            key = eval_mod.canonical_bpair_key(row["chr1"], int(row["start1"]), row["chr2"], int(row["start2"]))
            item = posterior.get(key)
            if item is None:
                continue
            for field in numeric_fields:
                if field in row and row[field] not in {"", "NA"}:
                    item[field] = float(row[field])
    return posterior


def format_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.9g}"
    return str(value)


def write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(row.get(key)) for key in fieldnames})


def canonical_pair(eval_mod: Any, chrom1: str, chrom2: str) -> tuple[str, str]:
    return tuple(sorted((chrom1, chrom2), key=eval_mod.chrom_sort_key))  # type: ignore[return-value]


def swap_code_to_canonical(pair: tuple[str, str], chrom1: str, chrom2: str, swap: int) -> int:
    if (chrom1, chrom2) == pair:
        return swap
    return ((swap & 1) << 1) | ((swap >> 1) & 1)


def aligned_p4(eval_mod: Any, p4: np.ndarray, swap_code: int) -> np.ndarray:
    return eval_mod.align_p4_to_truth_gauge(p4, (swap_code >> 1) & 1, swap_code & 1)


def all_copy_distances(
    eval_mod: Any,
    coords: dict[tuple[str, int, int], np.ndarray],
    item: dict[str, object],
) -> np.ndarray | None:
    distances = eval_mod.all_copy_pair_distances(
        coords,
        str(item["chrom1"]),
        int(item["start1"]),
        str(item["chrom2"]),
        int(item["start2"]),
    )
    if distances is None:
        return None
    return np.asarray(distances, dtype=float)


def load_eval_records(
    eval_mod: Any,
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    contact_counts: dict[tuple[str, int, str, int], np.ndarray],
    coords: dict[tuple[str, int, int], np.ndarray],
    reference: dict[tuple[str, int, int], np.ndarray],
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for key, counts_raw in contact_counts.items():
        item = posterior.get(key)
        if item is None:
            continue
        counts = np.asarray(counts_raw, dtype=np.int64)
        n_contacts = int(counts.sum())
        if n_contacts <= 0:
            continue
        if eval_mod.p4_from_coords_for_bpair(reference, item) is None:
            continue
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        if chrom1 == chrom2:
            continue
        p4 = np.asarray(item["p4"], dtype=float)
        records.append(
            {
                "pair": canonical_pair(eval_mod, chrom1, chrom2),
                "chrom1": chrom1,
                "chrom2": chrom2,
                "counts": counts,
                "p4": p4,
                "n_contacts": n_contacts,
            }
        )
    return records


def load_raw_selection_records(
    eval_mod: Any,
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    coords: dict[tuple[str, int, int], np.ndarray],
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for item in posterior.values():
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        if chrom1 == chrom2 or str(item.get("contact_class", "")) != "trans":
            continue
        n_raw = int(item.get("n_raw", 0))
        if n_raw <= 0:
            continue
        distances = all_copy_distances(eval_mod, coords, item)
        if distances is None:
            continue
        p4 = np.asarray(item["p4"], dtype=float)
        records.append(
            {
                "pair": canonical_pair(eval_mod, chrom1, chrom2),
                "chrom1": chrom1,
                "chrom2": chrom2,
                "p4": p4,
                "distances": distances,
                "n_raw": n_raw,
                "pmax": float(item.get("pmax", np.max(p4))),
                "margin": float(item.get("margin", 0.0)),
            }
        )
    return records


def chroms_from_coords(eval_mod: Any, coords: dict[tuple[str, int, int], np.ndarray]) -> list[str]:
    chroms = sorted({chrom for chrom, _start, _copy in coords}, key=eval_mod.chrom_sort_key)
    return chroms


def add_parity_score(scores: np.ndarray, parity: int, value: float) -> None:
    if parity == 0:
        scores[0] += value
        scores[3] += value
    else:
        scores[1] += value
        scores[2] += value


def raw_edge_scores(records: list[dict[str, object]]) -> dict[str, dict[tuple[str, str], np.ndarray]]:
    policies: dict[str, dict[tuple[str, str], np.ndarray]] = {
        "raw_global_geom_nearest_to_posterior_top": defaultdict(lambda: np.zeros(4, dtype=float)),
        "raw_global_geom_nearest_gap_weighted": defaultdict(lambda: np.zeros(4, dtype=float)),
        "raw_global_margin_geom_gap_weighted": defaultdict(lambda: np.zeros(4, dtype=float)),
        "raw_global_pmax_geom_gap_weighted": defaultdict(lambda: np.zeros(4, dtype=float)),
        "raw_global_same_cross_distance_sum": defaultdict(lambda: np.zeros(4, dtype=float)),
        "raw_global_posterior_same_cross_mass": defaultdict(lambda: np.zeros(4, dtype=float)),
    }
    for rec in records:
        pair = rec["pair"]  # type: ignore[assignment]
        chrom1 = str(rec["chrom1"])
        chrom2 = str(rec["chrom2"])
        n = float(rec["n_raw"])
        p4 = np.asarray(rec["p4"], dtype=float)
        distances = np.asarray(rec["distances"], dtype=float)
        posterior_top = int(np.argmax(p4))
        nearest = int(np.argmin(distances))
        gap = float(np.partition(distances, 1)[1] - np.partition(distances, 0)[0])
        weight = n * max(gap, 0.0)
        target = swap_code_to_canonical(pair, chrom1, chrom2, posterior_top ^ nearest)
        policies["raw_global_geom_nearest_to_posterior_top"][pair][target] += n
        policies["raw_global_geom_nearest_gap_weighted"][pair][target] += weight
        policies["raw_global_margin_geom_gap_weighted"][pair][target] += weight * max(float(rec["margin"]), 0.0)
        policies["raw_global_pmax_geom_gap_weighted"][pair][target] += weight * float(rec["pmax"])

        same_distance = float(distances[0] + distances[3])
        cross_distance = float(distances[1] + distances[2])
        distance_gap = abs(cross_distance - same_distance)
        distance_parity = 0 if same_distance <= cross_distance else 1
        add_parity_score(policies["raw_global_same_cross_distance_sum"][pair], distance_parity, n * distance_gap)

        same_mass = float(p4[0] + p4[3])
        cross_mass = float(p4[1] + p4[2])
        mass_gap = abs(same_mass - cross_mass)
        mass_parity = 0 if same_mass >= cross_mass else 1
        add_parity_score(policies["raw_global_posterior_same_cross_mass"][pair], mass_parity, n * mass_gap)
    return {policy: dict(by_pair) for policy, by_pair in policies.items()}


def center_edge_scores(
    coords: dict[tuple[str, int, int], np.ndarray],
    chroms: list[str],
) -> dict[str, dict[tuple[str, str], np.ndarray]]:
    by_copy: dict[tuple[str, int], list[np.ndarray]] = defaultdict(list)
    for (chrom, _start, copy), xyz in coords.items():
        by_copy[(str(chrom), int(copy))].append(np.asarray(xyz, dtype=float))
    centers: dict[tuple[str, int], np.ndarray] = {}
    for key, values in by_copy.items():
        if values:
            centers[key] = np.vstack(values).mean(axis=0)

    same_center: dict[tuple[str, str], np.ndarray] = {}
    vector_dot: dict[tuple[str, str], np.ndarray] = {}
    for i, chrom1 in enumerate(chroms):
        for chrom2 in chroms[i + 1 :]:
            c10 = centers.get((chrom1, 0))
            c11 = centers.get((chrom1, 1))
            c20 = centers.get((chrom2, 0))
            c21 = centers.get((chrom2, 1))
            if c10 is None or c11 is None or c20 is None or c21 is None:
                continue
            pair = (chrom1, chrom2)
            same = float(np.linalg.norm(c10 - c20) + np.linalg.norm(c11 - c21))
            cross = float(np.linalg.norm(c10 - c21) + np.linalg.norm(c11 - c20))
            vals = np.zeros(4, dtype=float)
            add_parity_score(vals, 0 if same <= cross else 1, abs(cross - same))
            same_center[pair] = vals

            v1 = c11 - c10
            v2 = c21 - c20
            dot = float(np.dot(v1, v2))
            vals2 = np.zeros(4, dtype=float)
            add_parity_score(vals2, 0, dot)
            add_parity_score(vals2, 1, -dot)
            vector_dot[pair] = vals2
    return {
        "raw_global_chrom_center_same_closer": same_center,
        "raw_global_homolog_vector_dot": vector_dot,
    }


def eval_oracle_edge_scores(
    eval_mod: Any,
    records: list[dict[str, object]],
    cis_swaps: dict[str, int],
) -> dict[tuple[str, str], np.ndarray]:
    scores: dict[tuple[str, str], np.ndarray] = defaultdict(lambda: np.zeros(4, dtype=float))
    for rec in records:
        pair = rec["pair"]  # type: ignore[assignment]
        chrom1 = str(rec["chrom1"])
        chrom2 = str(rec["chrom2"])
        p4 = np.asarray(rec["p4"], dtype=float)
        counts = np.asarray(rec["counts"], dtype=np.int64)
        base_swap = ((int(cis_swaps.get(chrom1, 0)) & 1) << 1) | (int(cis_swaps.get(chrom2, 0)) & 1)
        for extra in range(4):
            p4_aligned = aligned_p4(eval_mod, p4, base_swap ^ extra)
            extra_canonical = swap_code_to_canonical(pair, chrom1, chrom2, extra)
            scores[pair][extra_canonical] += float(counts[int(np.argmax(p4_aligned))])
    return dict(scores)


def solve_global_swaps(
    chroms: list[str],
    edge_scores: dict[tuple[str, str], np.ndarray],
) -> tuple[dict[str, int], float]:
    if not chroms:
        return {}, 0.0
    if len(chroms) == 1 or not edge_scores:
        return {chrom: 0 for chrom in chroms}, 0.0
    if len(chroms) > 30:
        raise ValueError(f"too many chromosomes for exact global gauge search: {len(chroms)}")

    index = {chrom: i for i, chrom in enumerate(chroms)}
    n_free = len(chroms) - 1
    assignments = np.arange(1 << n_free, dtype=np.uint32)
    total = np.zeros(assignments.shape[0], dtype=np.float64)
    bit_cache: dict[int, np.ndarray] = {0: np.zeros(assignments.shape[0], dtype=np.uint8)}

    def bits_for(chrom: str) -> np.ndarray:
        idx = index[chrom]
        if idx not in bit_cache:
            bit_cache[idx] = ((assignments >> (idx - 1)) & 1).astype(np.uint8)
        return bit_cache[idx]

    for pair, vals_raw in edge_scores.items():
        chrom1, chrom2 = pair
        if chrom1 not in index or chrom2 not in index:
            continue
        vals = np.asarray(vals_raw, dtype=np.float64)
        b1 = bits_for(chrom1)
        b2 = bits_for(chrom2)
        states = (b1 << 1) | b2
        total += vals[states]

    best_idx = int(np.argmax(total))
    best_score = float(total[best_idx])
    swaps = {chroms[0]: 0}
    for idx, chrom in enumerate(chroms[1:], start=1):
        swaps[chrom] = int((best_idx >> (idx - 1)) & 1)
    return swaps, best_score


def evaluate_global_swaps(
    eval_mod: Any,
    config_name: str,
    records: list[dict[str, object]],
    cis_swaps: dict[str, int],
    policies: dict[str, dict[str, int]],
    source_by_policy: dict[str, str],
    score_by_policy: dict[str, float],
    meta_by_policy: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for policy, chrom_swaps in policies.items():
        stats = eval_mod.empty_accuracy_stats()
        for rec in records:
            chrom1 = str(rec["chrom1"])
            chrom2 = str(rec["chrom2"])
            selected = ((int(chrom_swaps.get(chrom1, 0)) & 1) << 1) | (int(chrom_swaps.get(chrom2, 0)) & 1)
            base_swap = ((int(cis_swaps.get(chrom1, 0)) & 1) << 1) | (int(cis_swaps.get(chrom2, 0)) & 1)
            p4 = aligned_p4(eval_mod, np.asarray(rec["p4"], dtype=float), base_swap ^ selected)
            eval_mod.update_accuracy_stats(stats, np.asarray(rec["counts"], dtype=np.int64), p4)
        final = eval_mod.finalize_accuracy_stats(stats)
        meta = meta_by_policy.get(policy, {})
        rows.append(
            {
                "config_name": config_name,
                "policy": policy,
                "scope": "genome_trans",
                "flip_source": source_by_policy.get(policy, "blind"),
                "objective_score": score_by_policy.get(policy, "NA"),
                "n_eval_contacts": final["n_eval_contacts"],
                "top1_correct_contacts": final["top1_correct_contacts"],
                "top1_accuracy": final["top1_accuracy"],
                "same_cross_correct_contacts": final["same_cross_correct_contacts"],
                "same_cross_accuracy": final["same_cross_accuracy"],
                "pmax_threshold": final["pmax_threshold"],
                "n_called_contacts": final["n_called_contacts"],
                "called_contact_fraction": final["called_contact_fraction"],
                "pmax_threshold_accuracy": final["pmax_threshold_accuracy"],
                "pmax_threshold_recall": final["pmax_threshold_recall"],
                "uses_snp_for_base_gauge": 1,
                "uses_snp_for_global_gauge_selection": meta.get("uses_snp_for_global_gauge_selection", 0),
                "uses_snp_labeled_denominator_for_selection_weight": meta.get("uses_snp_labeled_denominator_for_selection_weight", 0),
                "uses_charm_for_denominator_filter": 1,
                "decoder_selection_denominator": meta.get("decoder_selection_denominator", "NA"),
                "decoder_selection_n_edges": meta.get("decoder_selection_n_edges", "NA"),
                "decoder_selection_total_weight": meta.get("decoder_selection_total_weight", "NA"),
                "eval_only_uses_snp_for_selection": meta.get("uses_snp_for_global_gauge_selection", 0),
                "eval_only_uses_snp_for_scoring": 1,
                "eval_only_uses_charm_shared_denominator": 1,
            }
        )
    baseline = next((float(r["top1_accuracy"]) for r in rows if r["policy"] == "cis_selected_whole_chrom"), None)
    for row in rows:
        if baseline is None:
            row["delta_top1_vs_cis_selected"] = "NA"
            row["target_plus_0p1_met"] = 0
        else:
            delta = float(row["top1_accuracy"]) - baseline
            row["delta_top1_vs_cis_selected"] = delta
            row["target_plus_0p1_met"] = int(delta >= 0.1)
    return rows


def chrom_swap_rows(
    config_name: str,
    chroms: list[str],
    policies: dict[str, dict[str, int]],
    source_by_policy: dict[str, str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for policy, swaps in policies.items():
        for chrom in chroms:
            rows.append(
                {
                    "config_name": config_name,
                    "policy": policy,
                    "flip_source": source_by_policy.get(policy, "blind"),
                    "chrom": chrom,
                    "copy_swap": int(swaps.get(chrom, 0)),
                }
            )
    return rows


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config-output-dir", type=Path, required=True)
    ap.add_argument("--eval-dir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--pairs", type=Path, default=None)
    ap.add_argument("--reference-3dg", type=Path, default=None)
    ap.add_argument("--bin-size", type=int, default=None)
    ap.add_argument("--config-name", default=None)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    eval_mod = load_eval_module()
    summary = read_kv(args.eval_dir / "summary.tsv")
    pairs = args.pairs or Path(summary["pairs_path"])
    reference_path = args.reference_3dg or Path(summary["reference_3dg_path"])
    bin_size = args.bin_size or int(summary["bin_size_bp"])
    config_name = args.config_name or args.config_output_dir.name

    posterior = read_posterior_with_extra(eval_mod, args.config_output_dir / "p9016_full.bpair_posterior.tsv")
    coords = eval_mod.read_reconstruction(args.config_output_dir / "p9016_full.coords.tsv")
    reference = eval_mod.read_reference_3dg(reference_path, bin_size)
    contact_counts = eval_mod.read_contact_truth_counts(pairs, bin_size, posterior)
    cis_swaps = read_selected_swaps(args.eval_dir / "whole_chrom_snp_oracle_swaps.tsv")
    eval_records = load_eval_records(eval_mod, posterior, contact_counts, coords, reference)
    raw_records = load_raw_selection_records(eval_mod, posterior, coords)
    chroms = chroms_from_coords(eval_mod, coords)

    policies: dict[str, dict[str, int]] = {
        "cis_selected_whole_chrom": {chrom: 0 for chrom in chroms},
    }
    source_by_policy = {"cis_selected_whole_chrom": "cis_snp_truth_whole_chrom"}
    score_by_policy: dict[str, float] = {"cis_selected_whole_chrom": 0.0}
    meta_by_policy: dict[str, dict[str, object]] = {
        "cis_selected_whole_chrom": {
            "uses_snp_for_global_gauge_selection": 0,
            "uses_snp_labeled_denominator_for_selection_weight": 0,
            "decoder_selection_denominator": "none_cis_snp_gauge_only",
            "decoder_selection_n_edges": 0,
            "decoder_selection_total_weight": 0,
        }
    }

    oracle_scores = eval_oracle_edge_scores(eval_mod, eval_records, cis_swaps)
    swaps, score = solve_global_swaps(chroms, oracle_scores)
    policies["global_chrom_snp_oracle"] = swaps
    source_by_policy["global_chrom_snp_oracle"] = "snp_truth_global_chrom_oracle"
    score_by_policy["global_chrom_snp_oracle"] = score
    meta_by_policy["global_chrom_snp_oracle"] = {
        "uses_snp_for_global_gauge_selection": 1,
        "uses_snp_labeled_denominator_for_selection_weight": 1,
        "decoder_selection_denominator": "snp_truth_eval_contacts",
        "decoder_selection_n_edges": len(oracle_scores),
        "decoder_selection_total_weight": sum(int(np.asarray(rec["counts"], dtype=np.int64).sum()) for rec in eval_records),
    }

    raw_policy_scores = raw_edge_scores(raw_records)
    raw_policy_scores.update(center_edge_scores(coords, chroms))
    raw_total_weight = sum(int(rec["n_raw"]) for rec in raw_records)
    for policy, edge_scores in raw_policy_scores.items():
        swaps, score = solve_global_swaps(chroms, edge_scores)
        policies[policy] = swaps
        source_by_policy[policy] = "blind_raw_posterior_reconstruction_global_gauge"
        score_by_policy[policy] = score
        meta_by_policy[policy] = {
            "uses_snp_for_global_gauge_selection": 0,
            "uses_snp_labeled_denominator_for_selection_weight": 0,
            "decoder_selection_denominator": "raw_posterior_trans_bpair" if policy.startswith("raw_global_") and "chrom_center" not in policy and "homolog_vector" not in policy else "raw_reconstruction_chrom_centers",
            "decoder_selection_n_edges": len(edge_scores),
            "decoder_selection_total_weight": raw_total_weight if "chrom_center" not in policy and "homolog_vector" not in policy else "NA",
        }

    summary_rows = evaluate_global_swaps(
        eval_mod,
        config_name,
        eval_records,
        cis_swaps,
        policies,
        source_by_policy,
        score_by_policy,
        meta_by_policy,
    )
    swap_rows = chrom_swap_rows(config_name, chroms, policies, source_by_policy)
    summary_fields = [
        "config_name",
        "policy",
        "scope",
        "flip_source",
        "objective_score",
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
        "delta_top1_vs_cis_selected",
        "target_plus_0p1_met",
        "uses_snp_for_base_gauge",
        "uses_snp_for_global_gauge_selection",
        "uses_snp_labeled_denominator_for_selection_weight",
        "uses_charm_for_denominator_filter",
        "decoder_selection_denominator",
        "decoder_selection_n_edges",
        "decoder_selection_total_weight",
        "eval_only_uses_snp_for_selection",
        "eval_only_uses_snp_for_scoring",
        "eval_only_uses_charm_shared_denominator",
    ]
    swap_fields = ["config_name", "policy", "flip_source", "chrom", "copy_swap"]
    write_rows(args.outdir / "global_gauge_decoder_summary.tsv", summary_rows, summary_fields)
    write_rows(args.outdir / "global_gauge_decoder_chrom_swaps.tsv", swap_rows, swap_fields)
    print(f"wrote\t{args.outdir / 'global_gauge_decoder_summary.tsv'}")
    print(f"wrote\t{args.outdir / 'global_gauge_decoder_chrom_swaps.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
