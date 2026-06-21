#!/usr/bin/env python3
"""Evaluate trans copy-gauge oracle ceilings for a finished P9016 run.

This script is eval-only. It reads SNP phase labels and CHARM/3DG coordinates
only to diagnose the metric ceiling after blind training has finished.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
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


def read_selected_swaps(path: Path) -> dict[str, int]:
    swaps: dict[str, int] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"chrom", "copy_swap", "selected_for_eval"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"{path} is missing required columns: {sorted(required)}")
        for row in reader:
            if int(row["selected_for_eval"]) == 1:
                swaps[row["chrom"]] = int(row["copy_swap"])
    return swaps


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fmt_float(value: float) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.9g}"


def valid_shared_denominator(eval_mod: Any, reference: dict[tuple[str, int, int], np.ndarray], item: dict[str, object]) -> bool:
    # Match contact_accuracy.tsv: require the CHARM/3DG FDG p4 to be computable.
    return eval_mod.p4_from_coords_for_bpair(reference, item) is not None


def iter_eval_contacts(
    eval_mod: Any,
    contact_counts: dict[tuple[str, int, str, int], np.ndarray],
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    reference: dict[tuple[str, int, int], np.ndarray],
):
    for key, counts in contact_counts.items():
        item = posterior.get(key)
        if item is None:
            continue
        if int(np.asarray(counts, dtype=np.int64).sum()) <= 0:
            continue
        if not valid_shared_denominator(eval_mod, reference, item):
            continue
        yield key, np.asarray(counts, dtype=np.int64), item


def evaluate_fixed_swaps(
    eval_mod: Any,
    config_name: str,
    policy: str,
    swaps: dict[str, int],
    contact_counts: dict[tuple[str, int, str, int], np.ndarray],
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    reference: dict[tuple[str, int, int], np.ndarray],
    optimized_on: str,
    swap_json_limit: int = 4096,
) -> list[dict[str, object]]:
    stats_by_scope: dict[str, dict[str, float]] = {
        "genome_all": eval_mod.empty_accuracy_stats(),
        "genome_cis": eval_mod.empty_accuracy_stats(),
        "genome_trans": eval_mod.empty_accuracy_stats(),
    }
    for _, counts, item in iter_eval_contacts(eval_mod, contact_counts, posterior, reference):
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        p4 = eval_mod.align_p4_to_truth_gauge(
            np.asarray(item["p4"], dtype=float),
            int(swaps.get(chrom1, 0)),
            int(swaps.get(chrom2, 0)),
        )
        eval_mod.update_accuracy_stats(stats_by_scope["genome_all"], counts, p4)
        scope = "genome_cis" if chrom1 == chrom2 else "genome_trans"
        eval_mod.update_accuracy_stats(stats_by_scope[scope], counts, p4)

    swaps_json = json.dumps(swaps, sort_keys=True)
    if len(swaps_json) > swap_json_limit:
        swaps_json = swaps_json[:swap_json_limit] + "...truncated"
    rows: list[dict[str, object]] = []
    for scope in ("genome_all", "genome_cis", "genome_trans"):
        finalized = eval_mod.finalize_accuracy_stats(stats_by_scope[scope])
        rows.append(
            {
                "config_name": config_name,
                "policy": policy,
                "scope": scope,
                "optimized_on": optimized_on,
                "n_eval_contacts": finalized["n_eval_contacts"],
                "top1_correct_contacts": finalized["top1_correct_contacts"],
                "top1_accuracy": fmt_float(float(finalized["top1_accuracy"])),
                "same_cross_correct_contacts": finalized["same_cross_correct_contacts"],
                "same_cross_accuracy": fmt_float(float(finalized["same_cross_accuracy"])),
                "pmax_threshold": finalized["pmax_threshold"],
                "n_called_contacts": finalized["n_called_contacts"],
                "called_contact_fraction": fmt_float(float(finalized["called_contact_fraction"])),
                "pmax_threshold_accuracy": fmt_float(float(finalized["pmax_threshold_accuracy"])),
                "pmax_threshold_recall": fmt_float(float(finalized["pmax_threshold_recall"])),
                "copy_swaps_json": swaps_json,
                "eval_only_uses_snp_truth": 1,
                "eval_only_uses_charm_shared_denominator": 1,
            }
        )
    return rows


def build_pair_scores(
    eval_mod: Any,
    contact_counts: dict[tuple[str, int, str, int], np.ndarray],
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    reference: dict[tuple[str, int, int], np.ndarray],
    chroms: list[str],
) -> tuple[dict[tuple[int, int], np.ndarray], dict[tuple[int, int], np.ndarray], dict[tuple[int, int], int]]:
    chrom_index = {chrom: idx for idx, chrom in enumerate(chroms)}
    top1_scores: dict[tuple[int, int], np.ndarray] = defaultdict(lambda: np.zeros((2, 2), dtype=np.int64))
    same_cross_scores: dict[tuple[int, int], np.ndarray] = defaultdict(lambda: np.zeros((2, 2), dtype=np.int64))
    contacts: dict[tuple[int, int], int] = defaultdict(int)
    for _, counts, item in iter_eval_contacts(eval_mod, contact_counts, posterior, reference):
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        if chrom1 == chrom2:
            continue
        if chrom1 not in chrom_index or chrom2 not in chrom_index:
            continue
        i = chrom_index[chrom1]
        j = chrom_index[chrom2]
        if i > j:
            i, j = j, i
            # Posterior keys are canonical, so this should not happen for sorted
            # chromosome order, but keep the pair key stable if it does.
        pair = (i, j)
        contacts[pair] += int(counts.sum())
        base_p4 = np.asarray(item["p4"], dtype=float)
        for swap1 in (0, 1):
            for swap2 in (0, 1):
                p4 = eval_mod.align_p4_to_truth_gauge(base_p4, swap1, swap2)
                pred = int(np.argmax(p4))
                top1_scores[pair][swap1, swap2] += int(counts[pred])
                pred_same = pred in (0, 3)
                same_counts = int(counts[0] + counts[3])
                cross_counts = int(counts[1] + counts[2])
                same_cross_scores[pair][swap1, swap2] += same_counts if pred_same else cross_counts
    return dict(top1_scores), dict(same_cross_scores), dict(contacts)


def optimize_global_trans_swaps(pair_scores: dict[tuple[int, int], np.ndarray], n_chroms: int) -> tuple[dict[int, int], int]:
    if n_chroms <= 0:
        return {}, 0
    if n_chroms > 24:
        raise ValueError(f"too many chromosomes for exact global oracle: {n_chroms}")
    masks = np.arange(1 << n_chroms, dtype=np.uint32)
    total = np.zeros(len(masks), dtype=np.int64)
    for (i, j), score in pair_scores.items():
        flat = np.asarray(score, dtype=np.int64).reshape(4)
        code = (((masks >> np.uint32(i)) & 1) << 1) | ((masks >> np.uint32(j)) & 1)
        total += flat[code]
    best_idx = int(np.argmax(total))
    best_mask = int(masks[best_idx])
    swaps = {idx: (best_mask >> idx) & 1 for idx in range(n_chroms)}
    return swaps, int(total[best_idx])


def evaluate_pair_independent_policy(
    eval_mod: Any,
    config_name: str,
    policy: str,
    pair_swaps: dict[tuple[str, str], tuple[int, int]],
    contact_counts: dict[tuple[str, int, str, int], np.ndarray],
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    reference: dict[tuple[str, int, int], np.ndarray],
) -> dict[str, object]:
    stats = eval_mod.empty_accuracy_stats()
    for _, counts, item in iter_eval_contacts(eval_mod, contact_counts, posterior, reference):
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        if chrom1 == chrom2:
            continue
        pair_key = tuple(sorted((chrom1, chrom2), key=eval_mod.chrom_sort_key))
        swap1, swap2 = pair_swaps[pair_key]
        if (chrom1, chrom2) != pair_key:
            swap1, swap2 = swap2, swap1
        p4 = eval_mod.align_p4_to_truth_gauge(np.asarray(item["p4"], dtype=float), swap1, swap2)
        eval_mod.update_accuracy_stats(stats, counts, p4)
    finalized = eval_mod.finalize_accuracy_stats(stats)
    return {
        "config_name": config_name,
        "policy": policy,
        "scope": "genome_trans",
        "optimized_on": "trans_top1_pair_independent",
        "n_eval_contacts": finalized["n_eval_contacts"],
        "top1_correct_contacts": finalized["top1_correct_contacts"],
        "top1_accuracy": fmt_float(float(finalized["top1_accuracy"])),
        "same_cross_correct_contacts": finalized["same_cross_correct_contacts"],
        "same_cross_accuracy": fmt_float(float(finalized["same_cross_accuracy"])),
        "pmax_threshold": finalized["pmax_threshold"],
        "n_called_contacts": finalized["n_called_contacts"],
        "called_contact_fraction": fmt_float(float(finalized["called_contact_fraction"])),
        "pmax_threshold_accuracy": fmt_float(float(finalized["pmax_threshold_accuracy"])),
        "pmax_threshold_recall": fmt_float(float(finalized["pmax_threshold_recall"])),
        "copy_swaps_json": "pair_independent_in_pair_tsv",
        "eval_only_uses_snp_truth": 1,
        "eval_only_uses_charm_shared_denominator": 1,
    }


def pair_policy_rows(
    chroms: list[str],
    pair_scores: dict[tuple[int, int], np.ndarray],
    pair_same_cross_scores: dict[tuple[int, int], np.ndarray],
    pair_contacts: dict[tuple[int, int], int],
    cis_swaps: dict[str, int],
    global_swaps_by_index: dict[int, int],
) -> tuple[list[dict[str, object]], dict[tuple[str, str], tuple[int, int]]]:
    rows: list[dict[str, object]] = []
    selected_pair_swaps: dict[tuple[str, str], tuple[int, int]] = {}
    for pair in sorted(pair_contacts, key=lambda p: (chroms[p[0]], chroms[p[1]])):
        i, j = pair
        chrom1 = chroms[i]
        chrom2 = chroms[j]
        n = int(pair_contacts[pair])
        score = pair_scores[pair]
        same_score = pair_same_cross_scores[pair]
        best_flat = int(np.argmax(score.reshape(4)))
        best_swap1 = best_flat >> 1
        best_swap2 = best_flat & 1
        selected_pair_swaps[(chrom1, chrom2)] = (best_swap1, best_swap2)
        global_swap1 = int(global_swaps_by_index.get(i, 0))
        global_swap2 = int(global_swaps_by_index.get(j, 0))
        cis_swap1 = int(cis_swaps.get(chrom1, 0))
        cis_swap2 = int(cis_swaps.get(chrom2, 0))
        pair_correct = int(score[best_swap1, best_swap2])
        global_correct = int(score[global_swap1, global_swap2])
        cis_correct = int(score[cis_swap1, cis_swap2])
        no_correct = int(score[0, 0])
        rows.append(
            {
                "chrom1": chrom1,
                "chrom2": chrom2,
                "n_eval_contacts": n,
                "pair_best_swap1": best_swap1,
                "pair_best_swap2": best_swap2,
                "global_swap1": global_swap1,
                "global_swap2": global_swap2,
                "cis_selected_swap1": cis_swap1,
                "cis_selected_swap2": cis_swap2,
                "top1_correct_no_swap": no_correct,
                "top1_accuracy_no_swap": fmt_float(no_correct / n if n else float("nan")),
                "top1_correct_cis_selected": cis_correct,
                "top1_accuracy_cis_selected": fmt_float(cis_correct / n if n else float("nan")),
                "top1_correct_global_trans_oracle": global_correct,
                "top1_accuracy_global_trans_oracle": fmt_float(global_correct / n if n else float("nan")),
                "top1_correct_pair_independent_oracle": pair_correct,
                "top1_accuracy_pair_independent_oracle": fmt_float(pair_correct / n if n else float("nan")),
                "pair_minus_global_top1_accuracy": fmt_float((pair_correct - global_correct) / n if n else float("nan")),
                "pair_minus_cis_selected_top1_accuracy": fmt_float((pair_correct - cis_correct) / n if n else float("nan")),
                "same_cross_accuracy_pair_independent_oracle": fmt_float(
                    int(same_score[best_swap1, best_swap2]) / n if n else float("nan")
                ),
                "score00": int(score[0, 0]),
                "score01": int(score[0, 1]),
                "score10": int(score[1, 0]),
                "score11": int(score[1, 1]),
            }
        )
    return rows, selected_pair_swaps


def add_deltas(rows: list[dict[str, object]]) -> None:
    baseline_by_config_scope: dict[tuple[str, str], float] = {}
    for row in rows:
        if row["policy"] == "cis_selected_whole_chrom_swaps":
            baseline_by_config_scope[(str(row["config_name"]), str(row["scope"]))] = float(row["top1_accuracy"])
    for row in rows:
        baseline = baseline_by_config_scope.get((str(row["config_name"]), str(row["scope"])))
        value = float(row["top1_accuracy"])
        row["delta_top1_accuracy_vs_cis_selected"] = fmt_float(value - baseline) if baseline is not None else "NA"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True, help="Pairs file with eval-only SNP phase columns.")
    parser.add_argument("--posterior", type=Path, required=True, help="Finished p9016_full.bpair_posterior.tsv.")
    parser.add_argument("--reference-3dg", type=Path, required=True, help="CHARM/3DG reference coordinates for denominator matching.")
    parser.add_argument("--whole-chrom-swaps", type=Path, required=True, help="Existing whole_chrom_snp_oracle_swaps.tsv.")
    parser.add_argument("--bin-size", type=int, default=1_000_000)
    parser.add_argument("--config-name", required=True)
    parser.add_argument("--out", type=Path, required=True, help="Main oracle summary TSV.")
    parser.add_argument("--pair-out", type=Path, default=None, help="Chromosome-pair detail TSV.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    eval_mod = load_eval_module()
    posterior = eval_mod.read_posterior(args.posterior)
    reference = eval_mod.read_reference_3dg(args.reference_3dg, args.bin_size)
    contact_counts = eval_mod.read_contact_truth_counts(args.pairs, args.bin_size, posterior)
    chroms = sorted(
        {
            str(item["chrom1"])
            for item in posterior.values()
        }
        | {str(item["chrom2"]) for item in posterior.values()},
        key=eval_mod.chrom_sort_key,
    )
    cis_swaps = read_selected_swaps(args.whole_chrom_swaps)

    rows: list[dict[str, object]] = []
    rows.extend(
        evaluate_fixed_swaps(
            eval_mod,
            args.config_name,
            "no_swaps",
            {},
            contact_counts,
            posterior,
            reference,
            "none",
        )
    )
    rows.extend(
        evaluate_fixed_swaps(
            eval_mod,
            args.config_name,
            "cis_selected_whole_chrom_swaps",
            cis_swaps,
            contact_counts,
            posterior,
            reference,
            "per_chrom_cis_top1",
        )
    )

    pair_scores, pair_same_cross_scores, pair_contacts = build_pair_scores(
        eval_mod, contact_counts, posterior, reference, chroms
    )
    global_swaps_by_index, _ = optimize_global_trans_swaps(pair_scores, len(chroms))
    global_swaps = {chroms[idx]: int(swap) for idx, swap in global_swaps_by_index.items()}
    rows.extend(
        evaluate_fixed_swaps(
            eval_mod,
            args.config_name,
            "global_trans_oracle_whole_chrom_swaps",
            global_swaps,
            contact_counts,
            posterior,
            reference,
            "genome_trans_top1",
        )
    )
    pair_rows, pair_swaps = pair_policy_rows(chroms, pair_scores, pair_same_cross_scores, pair_contacts, cis_swaps, global_swaps_by_index)
    rows.append(
        evaluate_pair_independent_policy(
            eval_mod,
            args.config_name,
            "pair_independent_trans_oracle_swaps",
            pair_swaps,
            contact_counts,
            posterior,
            reference,
        )
    )
    add_deltas(rows)
    write_table(args.out, rows)
    pair_out = args.pair_out if args.pair_out is not None else args.out.with_name(args.out.stem + ".pairs.tsv")
    write_table(pair_out, pair_rows)
    print(f"wrote\t{args.out}")
    print(f"wrote\t{pair_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
