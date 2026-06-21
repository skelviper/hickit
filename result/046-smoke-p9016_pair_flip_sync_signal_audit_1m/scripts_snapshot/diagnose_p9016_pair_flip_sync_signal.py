#!/usr/bin/env python3
"""Audit chromosome-pair copy-flip stability from blind post-training signals.

This diagnostic does not train a model.  It reads a completed blind P9016
reconstruction plus posterior table, chooses chromosome-pair copy flips from raw
posterior/reconstruction geometry only, and then scores those fixed choices with
SNP labels as eval-only truth.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np


STATE_COLUMNS = ("p00", "p01", "p10", "p11")


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
        numeric_fields = ("pU", "psame_raw", "pcross_raw", "entropy", "margin", "rho_output")
        for row in reader:
            key = eval_mod.canonical_bpair_key(row["chr1"], int(row["start1"]), row["chr2"], int(row["start2"]))
            item = posterior.get(key)
            if item is None:
                continue
            for field in numeric_fields:
                if field in row and row[field] not in {"", "NA"}:
                    item[field] = float(row[field])
            item["end1"] = int(row.get("end1", int(row["start1"]) + 1))
            item["end2"] = int(row.get("end2", int(row["start2"]) + 1))
    return posterior


def canonical_pair(eval_mod: Any, chrom1: str, chrom2: str) -> tuple[str, str]:
    return tuple(sorted((chrom1, chrom2), key=eval_mod.chrom_sort_key))  # type: ignore[return-value]


def swap_code_to_canonical(pair: tuple[str, str], chrom1: str, chrom2: str, swap: int) -> int:
    if (chrom1, chrom2) == pair:
        return swap
    return ((swap & 1) << 1) | ((swap >> 1) & 1)


def oriented_swap_for_item(pair: tuple[str, str], chrom1: str, chrom2: str, swap: int) -> int:
    if (chrom1, chrom2) == pair:
        return swap
    return ((swap & 1) << 1) | ((swap >> 1) & 1)


def split_id_for_record(chrom1: str, start1: int, chrom2: str, start2: int, seed: int) -> int:
    raw = f"{seed}|{chrom1}|{start1}|{chrom2}|{start2}".encode()
    digest = hashlib.blake2b(raw, digest_size=8).digest()
    return int.from_bytes(digest, "little") & 1


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
        # Match the standard eval denominator; CHARM/3DG is eval-only here.
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
        records.append(
            {
                "pair": canonical_pair(eval_mod, chrom1, chrom2),
                "chrom1": chrom1,
                "chrom2": chrom2,
                "start1": int(item["start1"]),
                "start2": int(item["start2"]),
                "counts": counts,
                "p4": p4,
                "distances": distances,
                "n_contacts": n_contacts,
                "pmax": float(item.get("pmax", np.max(p4))),
                "margin": float(item.get("margin", 0.0)),
            }
        )
    return records


def load_raw_records(
    eval_mod: Any,
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    coords: dict[tuple[str, int, int], np.ndarray],
    split_seed: int,
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
        start1 = int(item["start1"])
        start2 = int(item["start2"])
        records.append(
            {
                "pair": canonical_pair(eval_mod, chrom1, chrom2),
                "chrom1": chrom1,
                "chrom2": chrom2,
                "start1": start1,
                "start2": start2,
                "p4": p4,
                "distances": distances,
                "n_raw": n_raw,
                "pmax": float(item.get("pmax", np.max(p4))),
                "margin": float(item.get("margin", 0.0)),
                "posterior_top": int(np.argmax(p4)),
                "nearest_state": int(np.argmin(distances)),
                "split_id": split_id_for_record(chrom1, start1, chrom2, start2, split_seed),
            }
        )
    return records


def add_parity_score(scores: np.ndarray, parity: int, value: float) -> None:
    if parity == 0:
        scores[0] += value
        scores[3] += value
    else:
        scores[1] += value
        scores[2] += value


def compute_raw_policy_scores(records: list[dict[str, object]]) -> dict[str, dict[tuple[str, str], np.ndarray]]:
    policies: dict[str, dict[tuple[str, str], np.ndarray]] = {
        "geom_top_to_nearest": defaultdict(lambda: np.zeros(4, dtype=float)),
        "geom_nearest_gap": defaultdict(lambda: np.zeros(4, dtype=float)),
        "margin_geom_gap": defaultdict(lambda: np.zeros(4, dtype=float)),
        "pmax_geom_gap": defaultdict(lambda: np.zeros(4, dtype=float)),
        "same_cross_distance_sum": defaultdict(lambda: np.zeros(4, dtype=float)),
        "posterior_same_cross_mass": defaultdict(lambda: np.zeros(4, dtype=float)),
    }
    for rec in records:
        pair = rec["pair"]  # type: ignore[assignment]
        chrom1 = str(rec["chrom1"])
        chrom2 = str(rec["chrom2"])
        n = float(rec["n_raw"])
        p4 = np.asarray(rec["p4"], dtype=float)
        distances = np.asarray(rec["distances"], dtype=float)
        posterior_top = int(rec["posterior_top"])
        nearest = int(rec["nearest_state"])
        target = swap_code_to_canonical(pair, chrom1, chrom2, posterior_top ^ nearest)
        sorted_d = np.sort(distances)
        gap = float(sorted_d[1] - sorted_d[0]) if len(sorted_d) >= 2 else 0.0
        weight = n * max(gap, 0.0)
        policies["geom_top_to_nearest"][pair][target] += n
        policies["geom_nearest_gap"][pair][target] += weight
        policies["margin_geom_gap"][pair][target] += weight * max(float(rec["margin"]), 0.0)
        policies["pmax_geom_gap"][pair][target] += weight * float(rec["pmax"])

        same_distance = float(distances[0] + distances[3])
        cross_distance = float(distances[1] + distances[2])
        add_parity_score(
            policies["same_cross_distance_sum"][pair],
            0 if same_distance <= cross_distance else 1,
            n * abs(cross_distance - same_distance),
        )

        same_mass = float(p4[0] + p4[3])
        cross_mass = float(p4[1] + p4[2])
        add_parity_score(
            policies["posterior_same_cross_mass"][pair],
            0 if same_mass >= cross_mass else 1,
            n * abs(same_mass - cross_mass),
        )
    return {policy: dict(by_pair) for policy, by_pair in policies.items()}


def choose_from_scores(scores: dict[tuple[str, str], np.ndarray]) -> dict[tuple[str, str], int]:
    return {pair: int(np.argmax(np.asarray(vals, dtype=float))) for pair, vals in scores.items()}


def choose_oracle_scores(
    eval_mod: Any,
    records: list[dict[str, object]],
    cis_swaps: dict[str, int],
) -> dict[tuple[str, str], np.ndarray]:
    scores: dict[tuple[str, str], np.ndarray] = defaultdict(lambda: np.zeros(4, dtype=float))
    for rec in records:
        pair = rec["pair"]  # type: ignore[assignment]
        chrom1 = str(rec["chrom1"])
        chrom2 = str(rec["chrom2"])
        counts = np.asarray(rec["counts"], dtype=np.int64)
        p4 = np.asarray(rec["p4"], dtype=float)
        base_swap = ((int(cis_swaps.get(chrom1, 0)) & 1) << 1) | (int(cis_swaps.get(chrom2, 0)) & 1)
        for extra in range(4):
            oriented = oriented_swap_for_item(pair, chrom1, chrom2, extra)
            combined = base_swap ^ oriented
            aligned = eval_mod.align_p4_to_truth_gauge(p4, (combined >> 1) & 1, combined & 1)
            scores[pair][extra] += float(counts[int(np.argmax(aligned))])
    return dict(scores)


def evaluate_policy(
    eval_mod: Any,
    config_name: str,
    records: list[dict[str, object]],
    cis_swaps: dict[str, int],
    pair_swaps: dict[tuple[str, str], int],
    policy: str,
    policy_family: str,
    selection_source: str,
    meta: dict[str, object],
) -> dict[str, object]:
    stats = eval_mod.empty_accuracy_stats()
    for rec in records:
        pair = rec["pair"]  # type: ignore[assignment]
        chrom1 = str(rec["chrom1"])
        chrom2 = str(rec["chrom2"])
        selected = int(pair_swaps.get(pair, 0))
        oriented = oriented_swap_for_item(pair, chrom1, chrom2, selected)
        base_swap = ((int(cis_swaps.get(chrom1, 0)) & 1) << 1) | (int(cis_swaps.get(chrom2, 0)) & 1)
        combined = base_swap ^ oriented
        aligned = eval_mod.align_p4_to_truth_gauge(np.asarray(rec["p4"], dtype=float), (combined >> 1) & 1, combined & 1)
        eval_mod.update_accuracy_stats(stats, np.asarray(rec["counts"], dtype=np.int64), aligned)
    row = eval_mod.finalize_accuracy_stats(stats)
    row.update(
        {
            "config_name": config_name,
            "policy": policy,
            "policy_family": policy_family,
            "selection_source": selection_source,
            "scope": "genome_trans",
            "uses_snp_for_base_gauge": 1,
            "uses_snp_for_pair_selection": int(selection_source == "snp_truth_oracle"),
            "uses_snp_labeled_denominator_for_selection_weight": int(selection_source == "snp_truth_oracle"),
            "uses_charm_for_denominator_filter": 1,
            "eval_only_uses_snp_for_scoring": 1,
            "eval_only_uses_charm_shared_denominator": 1,
        }
    )
    row.update(meta)
    return row


def pair_eval_contact_counts(records: list[dict[str, object]]) -> dict[tuple[str, str], int]:
    counts: dict[tuple[str, str], int] = defaultdict(int)
    for rec in records:
        counts[rec["pair"]] += int(rec["n_contacts"])  # type: ignore[index]
    return dict(counts)


def raw_pair_weights(records: list[dict[str, object]]) -> dict[tuple[str, str], int]:
    weights: dict[tuple[str, str], int] = defaultdict(int)
    for rec in records:
        weights[rec["pair"]] += int(rec["n_raw"])  # type: ignore[index]
    return dict(weights)


def split_agreement_stats(
    all_pairs: set[tuple[str, str]],
    split_a: dict[tuple[str, str], int],
    split_b: dict[tuple[str, str], int],
    eval_weights: dict[tuple[str, str], int],
    raw_weights: dict[tuple[str, str], int],
) -> dict[str, object]:
    n_pairs = len(all_pairs)
    agree_pairs = {pair for pair in all_pairs if pair in split_a and pair in split_b and split_a[pair] == split_b[pair]}
    raw_total = sum(raw_weights.get(pair, 0) for pair in all_pairs)
    eval_total = sum(eval_weights.get(pair, 0) for pair in all_pairs)
    raw_agree = sum(raw_weights.get(pair, 0) for pair in agree_pairs)
    eval_agree = sum(eval_weights.get(pair, 0) for pair in agree_pairs)
    return {
        "selected_pair_count": n_pairs,
        "split_agree_pair_count": len(agree_pairs),
        "split_agree_pair_fraction": len(agree_pairs) / n_pairs if n_pairs else float("nan"),
        "split_agree_raw_weight": raw_agree,
        "split_agree_raw_weight_fraction": raw_agree / raw_total if raw_total else float("nan"),
        "split_agree_eval_contacts": eval_agree,
        "split_agree_eval_contact_fraction": eval_agree / eval_total if eval_total else float("nan"),
    }


def oracle_match_stats(
    all_pairs: set[tuple[str, str]],
    selected: dict[tuple[str, str], int],
    oracle: dict[tuple[str, str], int],
    eval_weights: dict[tuple[str, str], int],
    raw_weights: dict[tuple[str, str], int],
) -> dict[str, object]:
    comparable = {pair for pair in all_pairs if pair in oracle}
    matches = {pair for pair in comparable if int(selected.get(pair, 0)) == int(oracle[pair])}
    raw_total = sum(raw_weights.get(pair, 0) for pair in comparable)
    eval_total = sum(eval_weights.get(pair, 0) for pair in comparable)
    raw_match = sum(raw_weights.get(pair, 0) for pair in matches)
    eval_match = sum(eval_weights.get(pair, 0) for pair in matches)
    return {
        "oracle_comparable_pair_count": len(comparable),
        "oracle_match_pair_count": len(matches),
        "oracle_match_pair_fraction": len(matches) / len(comparable) if comparable else float("nan"),
        "oracle_match_raw_weight_fraction": raw_match / raw_total if raw_total else float("nan"),
        "oracle_match_eval_contact_fraction": eval_match / eval_total if eval_total else float("nan"),
    }


def make_split_agree_swaps(
    all_swaps: dict[tuple[str, str], int],
    split_a: dict[tuple[str, str], int],
    split_b: dict[tuple[str, str], int],
) -> dict[tuple[str, str], int]:
    out: dict[tuple[str, str], int] = {}
    for pair, selected in all_swaps.items():
        if pair in split_a and pair in split_b and split_a[pair] == split_b[pair]:
            out[pair] = int(selected)
        else:
            out[pair] = 0
    return out


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


def pair_rows_for_policies(
    eval_mod: Any,
    config_name: str,
    records: list[dict[str, object]],
    cis_swaps: dict[str, int],
    policies: dict[str, dict[tuple[str, str], int]],
    policy_family: dict[str, str],
    selection_source: dict[str, str],
    split_a_by_policy: dict[str, dict[tuple[str, str], int]],
    split_b_by_policy: dict[str, dict[tuple[str, str], int]],
    oracle_swaps: dict[tuple[str, str], int],
) -> list[dict[str, object]]:
    by_pair: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for rec in records:
        by_pair[rec["pair"]].append(rec)  # type: ignore[index]
    rows: list[dict[str, object]] = []
    for pair, pair_records in sorted(by_pair.items(), key=lambda x: (eval_mod.chrom_sort_key(x[0][0]), eval_mod.chrom_sort_key(x[0][1]))):
        for policy, swaps in policies.items():
            selected = int(swaps.get(pair, 0))
            stats = eval_mod.empty_accuracy_stats()
            for rec in pair_records:
                chrom1 = str(rec["chrom1"])
                chrom2 = str(rec["chrom2"])
                oriented = oriented_swap_for_item(pair, chrom1, chrom2, selected)
                base_swap = ((int(cis_swaps.get(chrom1, 0)) & 1) << 1) | (int(cis_swaps.get(chrom2, 0)) & 1)
                combined = base_swap ^ oriented
                aligned = eval_mod.align_p4_to_truth_gauge(np.asarray(rec["p4"], dtype=float), (combined >> 1) & 1, combined & 1)
                eval_mod.update_accuracy_stats(stats, np.asarray(rec["counts"], dtype=np.int64), aligned)
            final = eval_mod.finalize_accuracy_stats(stats)
            split_a = split_a_by_policy.get(policy, {}).get(pair)
            split_b = split_b_by_policy.get(policy, {}).get(pair)
            oracle = oracle_swaps.get(pair)
            rows.append(
                {
                    "config_name": config_name,
                    "chrom1": pair[0],
                    "chrom2": pair[1],
                    "policy": policy,
                    "policy_family": policy_family.get(policy, "NA"),
                    "selection_source": selection_source.get(policy, "NA"),
                    "selected_swap_code": selected,
                    "selected_swap1": (selected >> 1) & 1,
                    "selected_swap2": selected & 1,
                    "split_a_swap_code": split_a if split_a is not None else "NA",
                    "split_b_swap_code": split_b if split_b is not None else "NA",
                    "split_agree": int(split_a is not None and split_b is not None and split_a == split_b),
                    "oracle_swap_code": oracle if oracle is not None else "NA",
                    "oracle_matches_selected": int(oracle is not None and int(oracle) == selected),
                    "n_eval_contacts": final["n_eval_contacts"],
                    "top1_accuracy": final["top1_accuracy"],
                    "same_cross_accuracy": final["same_cross_accuracy"],
                    "pmax_threshold_accuracy": final["pmax_threshold_accuracy"],
                    "pmax_threshold_recall": final["pmax_threshold_recall"],
                }
            )
    return rows


def score_rows_for_policy(
    config_name: str,
    policy: str,
    scores_all: dict[tuple[str, str], np.ndarray],
    scores_a: dict[tuple[str, str], np.ndarray],
    scores_b: dict[tuple[str, str], np.ndarray],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    pairs = set(scores_all) | set(scores_a) | set(scores_b)
    for pair in sorted(pairs):
        all_vals = np.asarray(scores_all.get(pair, np.zeros(4)), dtype=float)
        a_vals = np.asarray(scores_a.get(pair, np.zeros(4)), dtype=float)
        b_vals = np.asarray(scores_b.get(pair, np.zeros(4)), dtype=float)
        rows.append(
            {
                "config_name": config_name,
                "policy": policy,
                "chrom1": pair[0],
                "chrom2": pair[1],
                "all_score0": all_vals[0],
                "all_score1": all_vals[1],
                "all_score2": all_vals[2],
                "all_score3": all_vals[3],
                "split_a_score0": a_vals[0],
                "split_a_score1": a_vals[1],
                "split_a_score2": a_vals[2],
                "split_a_score3": a_vals[3],
                "split_b_score0": b_vals[0],
                "split_b_score1": b_vals[1],
                "split_b_score2": b_vals[2],
                "split_b_score3": b_vals[3],
                "all_selected_swap_code": int(np.argmax(all_vals)),
                "split_a_selected_swap_code": int(np.argmax(a_vals)),
                "split_b_selected_swap_code": int(np.argmax(b_vals)),
                "split_agree": int(int(np.argmax(a_vals)) == int(np.argmax(b_vals))),
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
    ap.add_argument("--split-seed", type=int, default=17)
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
    raw_records = load_raw_records(eval_mod, posterior, coords, args.split_seed)

    raw_all_scores = compute_raw_policy_scores(raw_records)
    raw_a_scores = compute_raw_policy_scores([rec for rec in raw_records if int(rec["split_id"]) == 0])
    raw_b_scores = compute_raw_policy_scores([rec for rec in raw_records if int(rec["split_id"]) == 1])
    oracle_scores = choose_oracle_scores(eval_mod, eval_records, cis_swaps)
    oracle_swaps = choose_from_scores(oracle_scores)

    eval_weights = pair_eval_contact_counts(eval_records)
    raw_weights = raw_pair_weights(raw_records)
    all_pairs = set(raw_weights) | set(eval_weights)
    raw_total_records = len(raw_records)
    raw_total_weight = sum(raw_weights.values())

    policies: dict[str, dict[tuple[str, str], int]] = {
        "cis_selected_whole_chrom": {},
        "pair_independent_4state_snp_oracle": oracle_swaps,
    }
    policy_family = {
        "cis_selected_whole_chrom": "baseline_eval_gauge",
        "pair_independent_4state_snp_oracle": "eval_only_oracle",
    }
    selection_source = {
        "cis_selected_whole_chrom": "cis_snp_truth_whole_chrom",
        "pair_independent_4state_snp_oracle": "snp_truth_oracle",
    }
    split_a_by_policy: dict[str, dict[tuple[str, str], int]] = {}
    split_b_by_policy: dict[str, dict[tuple[str, str], int]] = {}
    audit_rows: list[dict[str, object]] = []
    score_rows: list[dict[str, object]] = []

    for base_policy in sorted(raw_all_scores):
        all_swaps = choose_from_scores(raw_all_scores[base_policy])
        split_a_swaps = choose_from_scores(raw_a_scores.get(base_policy, {}))
        split_b_swaps = choose_from_scores(raw_b_scores.get(base_policy, {}))
        agree_swaps = make_split_agree_swaps(all_swaps, split_a_swaps, split_b_swaps)

        all_name = f"raw_all_{base_policy}"
        a_name = f"raw_split_a_{base_policy}"
        b_name = f"raw_split_b_{base_policy}"
        agree_name = f"raw_split_agree_{base_policy}"
        for name, swaps, family in (
            (all_name, all_swaps, "blind_raw_all"),
            (a_name, split_a_swaps, "blind_raw_split_a"),
            (b_name, split_b_swaps, "blind_raw_split_b"),
            (agree_name, agree_swaps, "blind_raw_split_agree_else_zero"),
        ):
            policies[name] = swaps
            policy_family[name] = family
            selection_source[name] = "blind_raw_posterior_reconstruction_geometry"
            split_a_by_policy[name] = split_a_swaps
            split_b_by_policy[name] = split_b_swaps

        split_stats = split_agreement_stats(all_pairs, split_a_swaps, split_b_swaps, eval_weights, raw_weights)
        for name, swaps, family in (
            (all_name, all_swaps, "blind_raw_all"),
            (agree_name, agree_swaps, "blind_raw_split_agree_else_zero"),
        ):
            row = {
                "config_name": config_name,
                "base_policy": base_policy,
                "policy": name,
                "policy_family": family,
                "selection_source": "blind_raw_posterior_reconstruction_geometry",
                "raw_selection_n_records": raw_total_records,
                "raw_selection_total_weight": raw_total_weight,
            }
            row.update(split_stats)
            row.update(oracle_match_stats(all_pairs, swaps, oracle_swaps, eval_weights, raw_weights))
            audit_rows.append(row)
        score_rows.extend(
            score_rows_for_policy(
                config_name,
                base_policy,
                raw_all_scores[base_policy],
                raw_a_scores.get(base_policy, {}),
                raw_b_scores.get(base_policy, {}),
            )
        )

    summary_rows = []
    for policy, swaps in policies.items():
        meta = {
            "raw_selection_n_records": raw_total_records if policy.startswith("raw_") else 0,
            "raw_selection_total_weight": raw_total_weight if policy.startswith("raw_") else 0,
        }
        if policy.startswith("raw_"):
            base_policy = policy
            for prefix in ("raw_all_", "raw_split_a_", "raw_split_b_", "raw_split_agree_"):
                if base_policy.startswith(prefix):
                    base_policy = base_policy[len(prefix) :]
                    break
            meta.update(
                split_agreement_stats(
                    all_pairs,
                    split_a_by_policy.get(policy, {}),
                    split_b_by_policy.get(policy, {}),
                    eval_weights,
                    raw_weights,
                )
            )
            meta.update(oracle_match_stats(all_pairs, swaps, oracle_swaps, eval_weights, raw_weights))
            meta["base_policy"] = base_policy
        else:
            meta.update(
                {
                    "base_policy": "NA",
                    "selected_pair_count": len(all_pairs),
                    "split_agree_pair_count": "NA",
                    "split_agree_pair_fraction": "NA",
                    "split_agree_raw_weight": "NA",
                    "split_agree_raw_weight_fraction": "NA",
                    "split_agree_eval_contacts": "NA",
                    "split_agree_eval_contact_fraction": "NA",
                    "oracle_comparable_pair_count": len(oracle_swaps),
                    "oracle_match_pair_count": "NA",
                    "oracle_match_pair_fraction": "NA",
                    "oracle_match_raw_weight_fraction": "NA",
                    "oracle_match_eval_contact_fraction": "NA",
                }
            )
        summary_rows.append(
            evaluate_policy(
                eval_mod,
                config_name,
                eval_records,
                cis_swaps,
                swaps,
                policy,
                policy_family.get(policy, "NA"),
                selection_source.get(policy, "NA"),
                meta,
            )
        )

    baseline = next((float(row["top1_accuracy"]) for row in summary_rows if row["policy"] == "cis_selected_whole_chrom"), None)
    for row in summary_rows:
        if baseline is None:
            row["delta_top1_vs_cis_selected"] = "NA"
            row["target_plus_0p1_met"] = 0
        else:
            delta = float(row["top1_accuracy"]) - baseline
            row["delta_top1_vs_cis_selected"] = delta
            row["target_plus_0p1_met"] = int(delta >= 0.1)

    pair_rows = pair_rows_for_policies(
        eval_mod,
        config_name,
        eval_records,
        cis_swaps,
        policies,
        policy_family,
        selection_source,
        split_a_by_policy,
        split_b_by_policy,
        oracle_swaps,
    )

    summary_fields = [
        "config_name",
        "policy",
        "base_policy",
        "policy_family",
        "selection_source",
        "scope",
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
        "raw_selection_n_records",
        "raw_selection_total_weight",
        "selected_pair_count",
        "split_agree_pair_count",
        "split_agree_pair_fraction",
        "split_agree_raw_weight",
        "split_agree_raw_weight_fraction",
        "split_agree_eval_contacts",
        "split_agree_eval_contact_fraction",
        "oracle_comparable_pair_count",
        "oracle_match_pair_count",
        "oracle_match_pair_fraction",
        "oracle_match_raw_weight_fraction",
        "oracle_match_eval_contact_fraction",
        "uses_snp_for_base_gauge",
        "uses_snp_for_pair_selection",
        "uses_snp_labeled_denominator_for_selection_weight",
        "uses_charm_for_denominator_filter",
        "eval_only_uses_snp_for_scoring",
        "eval_only_uses_charm_shared_denominator",
    ]
    pair_fields = [
        "config_name",
        "chrom1",
        "chrom2",
        "policy",
        "policy_family",
        "selection_source",
        "selected_swap_code",
        "selected_swap1",
        "selected_swap2",
        "split_a_swap_code",
        "split_b_swap_code",
        "split_agree",
        "oracle_swap_code",
        "oracle_matches_selected",
        "n_eval_contacts",
        "top1_accuracy",
        "same_cross_accuracy",
        "pmax_threshold_accuracy",
        "pmax_threshold_recall",
    ]
    audit_fields = [
        "config_name",
        "base_policy",
        "policy",
        "policy_family",
        "selection_source",
        "raw_selection_n_records",
        "raw_selection_total_weight",
        "selected_pair_count",
        "split_agree_pair_count",
        "split_agree_pair_fraction",
        "split_agree_raw_weight",
        "split_agree_raw_weight_fraction",
        "split_agree_eval_contacts",
        "split_agree_eval_contact_fraction",
        "oracle_comparable_pair_count",
        "oracle_match_pair_count",
        "oracle_match_pair_fraction",
        "oracle_match_raw_weight_fraction",
        "oracle_match_eval_contact_fraction",
    ]
    score_fields = [
        "config_name",
        "policy",
        "chrom1",
        "chrom2",
        "all_score0",
        "all_score1",
        "all_score2",
        "all_score3",
        "split_a_score0",
        "split_a_score1",
        "split_a_score2",
        "split_a_score3",
        "split_b_score0",
        "split_b_score1",
        "split_b_score2",
        "split_b_score3",
        "all_selected_swap_code",
        "split_a_selected_swap_code",
        "split_b_selected_swap_code",
        "split_agree",
    ]
    write_rows(args.outdir / "pair_flip_sync_summary.tsv", summary_rows, summary_fields)
    write_rows(args.outdir / "pair_flip_sync_chr_pair.tsv", pair_rows, pair_fields)
    write_rows(args.outdir / "pair_flip_sync_audit.tsv", audit_rows, audit_fields)
    write_rows(args.outdir / "pair_flip_signal_scores.tsv", score_rows, score_fields)
    print(f"wrote\t{args.outdir / 'pair_flip_sync_summary.tsv'}")
    print(f"wrote\t{args.outdir / 'pair_flip_sync_chr_pair.tsv'}")
    print(f"wrote\t{args.outdir / 'pair_flip_sync_audit.tsv'}")
    print(f"wrote\t{args.outdir / 'pair_flip_signal_scores.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
