#!/usr/bin/env python3
"""Audit blind post-training trans copy-state decoders.

This diagnostic does not train a model.  It reads a finished reconstruction and
posterior, chooses chromosome-pair decoder rules using reconstruction/posterior
quantities only, and then scores those rules with SNP labels as eval-only truth.
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
    """Read posterior and keep diagnostic fields used by blind decoder rules."""
    posterior = eval_mod.read_posterior(path)
    with eval_mod.open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a header")
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


def oriented_swap_for_item(pair: tuple[str, str], chrom1: str, chrom2: str, swap: int) -> int:
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
        distances = all_copy_distances(eval_mod, coords, item)
        if distances is None:
            continue
        p4 = np.asarray(item["p4"], dtype=float)
        pair = canonical_pair(eval_mod, chrom1, chrom2)
        records.append(
            {
                "pair": pair,
                "chrom1": chrom1,
                "chrom2": chrom2,
                "counts": counts,
                "p4": p4,
                "distances": distances,
                "n_contacts": n_contacts,
                "pmax": float(np.max(p4)),
                "margin": float(item.get("margin", 0.0)),
            }
        )
    return records


def load_raw_selection_records(
    eval_mod: Any,
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    coords: dict[tuple[str, int, int], np.ndarray],
) -> list[dict[str, object]]:
    """Build raw-only records for decoder selection.

    These records use only the reconstruction and posterior table emitted by
    blind training.  SNP-labeled contact counts and CHARM/3DG reference
    coordinates are intentionally excluded here; they enter only in scoring.
    """
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
        pair = canonical_pair(eval_mod, chrom1, chrom2)
        records.append(
            {
                "pair": pair,
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


def swap_scores_for_record(eval_mod: Any, rec: dict[str, object], cis_swaps: dict[str, int]) -> np.ndarray:
    p4 = np.asarray(rec["p4"], dtype=float)
    counts = np.asarray(rec["counts"], dtype=np.int64)
    chrom1 = str(rec["chrom1"])
    chrom2 = str(rec["chrom2"])
    scores = np.zeros(4, dtype=np.int64)
    base_swap = ((int(cis_swaps.get(chrom1, 0)) & 1) << 1) | (int(cis_swaps.get(chrom2, 0)) & 1)
    for extra in range(4):
        p4_aligned = aligned_p4(eval_mod, p4, base_swap ^ extra)
        scores[extra] = int(counts[int(np.argmax(p4_aligned))])
    return scores


def evaluate_decoders(
    eval_mod: Any,
    config_name: str,
    records: list[dict[str, object]],
    cis_swaps: dict[str, int],
    decoders: dict[str, dict[tuple[str, str], int]],
    source_by_policy: dict[str, str],
    policy_meta: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for policy, pair_swaps in decoders.items():
        stats = eval_mod.empty_accuracy_stats()
        for rec in records:
            pair = rec["pair"]  # type: ignore[assignment]
            chrom1 = str(rec["chrom1"])
            chrom2 = str(rec["chrom2"])
            selected = int(pair_swaps.get(pair, 0))
            selected = oriented_swap_for_item(pair, chrom1, chrom2, selected)
            base_swap = ((int(cis_swaps.get(chrom1, 0)) & 1) << 1) | (int(cis_swaps.get(chrom2, 0)) & 1)
            p4 = aligned_p4(eval_mod, np.asarray(rec["p4"], dtype=float), base_swap ^ selected)
            eval_mod.update_accuracy_stats(stats, np.asarray(rec["counts"], dtype=np.int64), p4)
        final = eval_mod.finalize_accuracy_stats(stats)
        rows.append(
            {
                "config_name": config_name,
                "policy": policy,
                "scope": "genome_trans",
                "flip_source": source_by_policy.get(policy, "blind"),
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
                "uses_snp_for_base_gauge": policy_meta.get(policy, {}).get("uses_snp_for_base_gauge", 1),
                "uses_snp_for_pair_decoder_selection": policy_meta.get(policy, {}).get("uses_snp_for_pair_decoder_selection", 0),
                "uses_snp_labeled_denominator_for_selection_weight": policy_meta.get(policy, {}).get("uses_snp_labeled_denominator_for_selection_weight", 0),
                "uses_charm_for_denominator_filter": 1,
                "decoder_selection_denominator": policy_meta.get(policy, {}).get("decoder_selection_denominator", "NA"),
                "decoder_selection_n_records": policy_meta.get(policy, {}).get("decoder_selection_n_records", "NA"),
                "decoder_selection_total_weight": policy_meta.get(policy, {}).get("decoder_selection_total_weight", "NA"),
                "eval_only_uses_snp_for_selection": policy_meta.get(policy, {}).get("uses_snp_for_pair_decoder_selection", 0),
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


def pair_detail_rows(
    eval_mod: Any,
    records: list[dict[str, object]],
    cis_swaps: dict[str, int],
    decoders: dict[str, dict[tuple[str, str], int]],
    source_by_policy: dict[str, str],
) -> list[dict[str, object]]:
    by_pair: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for rec in records:
        by_pair[rec["pair"]].append(rec)  # type: ignore[index]
    rows: list[dict[str, object]] = []
    for pair, pair_records in sorted(by_pair.items(), key=lambda item: (eval_mod.chrom_sort_key(item[0][0]), eval_mod.chrom_sort_key(item[0][1]))):
        for policy, pair_swaps in decoders.items():
            stats = eval_mod.empty_accuracy_stats()
            selected = int(pair_swaps.get(pair, 0))
            for rec in pair_records:
                chrom1 = str(rec["chrom1"])
                chrom2 = str(rec["chrom2"])
                oriented = oriented_swap_for_item(pair, chrom1, chrom2, selected)
                base_swap = ((int(cis_swaps.get(chrom1, 0)) & 1) << 1) | (int(cis_swaps.get(chrom2, 0)) & 1)
                p4 = aligned_p4(eval_mod, np.asarray(rec["p4"], dtype=float), base_swap ^ oriented)
                eval_mod.update_accuracy_stats(stats, np.asarray(rec["counts"], dtype=np.int64), p4)
            final = eval_mod.finalize_accuracy_stats(stats)
            rows.append(
                {
                    "chrom1": pair[0],
                    "chrom2": pair[1],
                    "policy": policy,
                    "selected_swap_code": selected,
                    "selected_swap1": (selected >> 1) & 1,
                    "selected_swap2": selected & 1,
                    "flip_source": source_by_policy.get(policy, "blind"),
                    "n_eval_contacts": final["n_eval_contacts"],
                    "top1_accuracy": final["top1_accuracy"],
                    "same_cross_accuracy": final["same_cross_accuracy"],
                    "pmax_threshold_accuracy": final["pmax_threshold_accuracy"],
                    "pmax_threshold_recall": final["pmax_threshold_recall"],
                }
            )
    return rows


def choose_oracle_decoders(
    eval_mod: Any,
    records: list[dict[str, object]],
    cis_swaps: dict[str, int],
) -> dict[str, dict[tuple[str, str], int]]:
    scores: dict[tuple[str, str], np.ndarray] = defaultdict(lambda: np.zeros(4, dtype=np.int64))
    for rec in records:
        scores[rec["pair"]] += swap_scores_for_record(eval_mod, rec, cis_swaps)  # type: ignore[index]
    return {"pair_independent_4state_oracle": {pair: int(np.argmax(vals)) for pair, vals in scores.items()}}


def swap_code_to_canonical(pair: tuple[str, str], chrom1: str, chrom2: str, swap: int) -> int:
    if (chrom1, chrom2) == pair:
        return swap
    return ((swap & 1) << 1) | ((swap >> 1) & 1)


def choose_blind_decoders(records: list[dict[str, object]]) -> dict[str, dict[tuple[str, str], int]]:
    scores: dict[str, dict[tuple[str, str], np.ndarray]] = {
        "blind_geom_nearest_to_posterior_top": defaultdict(lambda: np.zeros(4, dtype=float)),
        "blind_geom_nearest_gap_weighted": defaultdict(lambda: np.zeros(4, dtype=float)),
        "blind_margin_geom_gap_weighted": defaultdict(lambda: np.zeros(4, dtype=float)),
        "blind_pmax_geom_gap_weighted": defaultdict(lambda: np.zeros(4, dtype=float)),
        "blind_same_cross_distance_sum": defaultdict(lambda: np.zeros(4, dtype=float)),
        "blind_posterior_mass_same_cross": defaultdict(lambda: np.zeros(4, dtype=float)),
    }
    for rec in records:
        pair = rec["pair"]  # type: ignore[assignment]
        n = float(rec["n_raw"])
        chrom1 = str(rec["chrom1"])
        chrom2 = str(rec["chrom2"])
        p4 = np.asarray(rec["p4"], dtype=float)
        distances = np.asarray(rec["distances"], dtype=float)
        posterior_top = int(np.argmax(p4))
        nearest = int(np.argmin(distances))
        gap = float(np.partition(distances, 1)[1] - np.partition(distances, 0)[0])
        margin = max(float(rec["margin"]), 0.0)
        pmax = float(rec["pmax"])
        target_swap = swap_code_to_canonical(pair, chrom1, chrom2, posterior_top ^ nearest)
        weight = n * max(gap, 0.0)
        scores["blind_geom_nearest_to_posterior_top"][pair][target_swap] += n
        scores["blind_geom_nearest_gap_weighted"][pair][target_swap] += weight
        scores["blind_margin_geom_gap_weighted"][pair][target_swap] += weight * margin
        scores["blind_pmax_geom_gap_weighted"][pair][target_swap] += weight * pmax
        same_distance = distances[0] + distances[3]
        cross_distance = distances[1] + distances[2]
        same_mass = p4[0] + p4[3]
        cross_mass = p4[1] + p4[2]
        rel = 0 if same_distance <= cross_distance else 1
        scores["blind_same_cross_distance_sum"][pair][rel] += n * abs(cross_distance - same_distance)
        rel2 = 0 if same_mass >= cross_mass else 1
        scores["blind_posterior_mass_same_cross"][pair][rel2] += n * abs(same_mass - cross_mass)
    out: dict[str, dict[tuple[str, str], int]] = {}
    for policy, by_pair in scores.items():
        out[policy] = {pair: int(np.argmax(vals)) for pair, vals in by_pair.items()}
    return out


def choose_center_decoders(coords: dict[tuple[str, int, int], np.ndarray], records: list[dict[str, object]]) -> dict[str, dict[tuple[str, str], int]]:
    by_copy: dict[tuple[str, int], list[np.ndarray]] = defaultdict(list)
    for (chrom, _start, copy), xyz in coords.items():
        by_copy[(str(chrom), int(copy))].append(np.asarray(xyz, dtype=float))
    centers: dict[tuple[str, int], np.ndarray] = {}
    for key, values in by_copy.items():
        if values:
            centers[key] = np.vstack(values).mean(axis=0)
    pairs = sorted({rec["pair"] for rec in records})
    same_center: dict[tuple[str, str], int] = {}
    vector_dot: dict[tuple[str, str], int] = {}
    for pair in pairs:
        chrom1, chrom2 = pair
        c10 = centers.get((chrom1, 0))
        c11 = centers.get((chrom1, 1))
        c20 = centers.get((chrom2, 0))
        c21 = centers.get((chrom2, 1))
        if c10 is None or c11 is None or c20 is None or c21 is None:
            same_center[pair] = 0
            vector_dot[pair] = 0
            continue
        same = float(np.linalg.norm(c10 - c20) + np.linalg.norm(c11 - c21))
        cross = float(np.linalg.norm(c10 - c21) + np.linalg.norm(c11 - c20))
        same_center[pair] = 1 if cross < same else 0
        v1 = c11 - c10
        v2 = c21 - c20
        vector_dot[pair] = 1 if float(np.dot(v1, v2)) < 0.0 else 0
    return {
        "blind_chrom_center_same_closer": same_center,
        "blind_chrom_homolog_vector_dot": vector_dot,
    }


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
    raw_selection_records = load_raw_selection_records(eval_mod, posterior, coords)

    decoders: dict[str, dict[tuple[str, str], int]] = {
        "cis_selected_whole_chrom": {},
    }
    source_by_policy = {
        "cis_selected_whole_chrom": "cis_snp_truth_whole_chrom",
    }
    policy_meta: dict[str, dict[str, object]] = {
        "cis_selected_whole_chrom": {
            "uses_snp_for_base_gauge": 1,
            "uses_snp_for_pair_decoder_selection": 0,
            "uses_snp_labeled_denominator_for_selection_weight": 0,
            "decoder_selection_denominator": "none_cis_snp_gauge_only",
            "decoder_selection_n_records": 0,
            "decoder_selection_total_weight": 0,
        }
    }
    for policy, decoder in choose_oracle_decoders(eval_mod, eval_records, cis_swaps).items():
        decoders[policy] = decoder
        source_by_policy[policy] = "snp_truth_oracle"
        policy_meta[policy] = {
            "uses_snp_for_base_gauge": 1,
            "uses_snp_for_pair_decoder_selection": 1,
            "uses_snp_labeled_denominator_for_selection_weight": 1,
            "decoder_selection_denominator": "snp_truth_eval_contacts",
            "decoder_selection_n_records": len(eval_records),
            "decoder_selection_total_weight": sum(int(np.asarray(rec["counts"], dtype=np.int64).sum()) for rec in eval_records),
        }
    for policy, decoder in choose_blind_decoders(raw_selection_records).items():
        decoders[policy] = decoder
        source_by_policy[policy] = "blind_posterior_reconstruction_geometry"
        policy_meta[policy] = {
            "uses_snp_for_base_gauge": 1,
            "uses_snp_for_pair_decoder_selection": 0,
            "uses_snp_labeled_denominator_for_selection_weight": 0,
            "decoder_selection_denominator": "raw_posterior_trans_bpair",
            "decoder_selection_n_records": len(raw_selection_records),
            "decoder_selection_total_weight": sum(int(rec["n_raw"]) for rec in raw_selection_records),
        }
    for policy, decoder in choose_center_decoders(coords, raw_selection_records).items():
        decoders[policy] = decoder
        source_by_policy[policy] = "blind_reconstruction_geometry"
        policy_meta[policy] = {
            "uses_snp_for_base_gauge": 1,
            "uses_snp_for_pair_decoder_selection": 0,
            "uses_snp_labeled_denominator_for_selection_weight": 0,
            "decoder_selection_denominator": "raw_reconstruction_all_trans_bpair",
            "decoder_selection_n_records": len(raw_selection_records),
            "decoder_selection_total_weight": sum(int(rec["n_raw"]) for rec in raw_selection_records),
        }

    summary_rows = evaluate_decoders(eval_mod, config_name, eval_records, cis_swaps, decoders, source_by_policy, policy_meta)
    pair_rows = pair_detail_rows(eval_mod, eval_records, cis_swaps, decoders, source_by_policy)
    summary_fields = [
        "config_name",
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
        "delta_top1_vs_cis_selected",
        "target_plus_0p1_met",
        "uses_snp_for_base_gauge",
        "uses_snp_for_pair_decoder_selection",
        "uses_snp_labeled_denominator_for_selection_weight",
        "uses_charm_for_denominator_filter",
        "decoder_selection_denominator",
        "decoder_selection_n_records",
        "decoder_selection_total_weight",
        "eval_only_uses_snp_for_selection",
        "eval_only_uses_snp_for_scoring",
        "eval_only_uses_charm_shared_denominator",
    ]
    pair_fields = [
        "chrom1",
        "chrom2",
        "policy",
        "selected_swap_code",
        "selected_swap1",
        "selected_swap2",
        "flip_source",
        "n_eval_contacts",
        "top1_accuracy",
        "same_cross_accuracy",
        "pmax_threshold_accuracy",
        "pmax_threshold_recall",
    ]
    write_rows(args.outdir / "trans_decoder_rule_summary.tsv", summary_rows, summary_fields)
    write_rows(args.outdir / "trans_decoder_rule_chr_pair.tsv", pair_rows, pair_fields)
    print(f"wrote\t{args.outdir / 'trans_decoder_rule_summary.tsv'}")
    print(f"wrote\t{args.outdir / 'trans_decoder_rule_chr_pair.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
