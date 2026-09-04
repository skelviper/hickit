#!/usr/bin/env python3
"""Logical Experiment 058: held-out leakage, stability, and blind selection audit."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import math
import os
import sys
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from phase_diagnostics.metrics import (
    contact_metrics,
    jensen_shannon_divergence,
    model_selection_metrics,
    structure_metrics,
)
from phase_diagnostics.provenance import collect_provenance, file_sha256, write_provenance
from phase_diagnostics.splits import Observation, SPLIT_MODES, blocked_split, leakage_audit
from phase_diagnostics.state import parity_probabilities, relabel_probabilities
from phase_diagnostics.synthetic import (
    InferenceConfig,
    SyntheticCondition,
    evaluate_synthetic,
    generate_synthetic,
    geometry_selected_chromosome_flips,
    infer_synthetic,
    relabel_contact_probabilities_to_geometry_gauge,
    score_coordinates,
    subset_synthetic,
)


DEFAULT_PAIRS = Path("/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz")
HISTORICAL_078 = Path(
    "/work/phase3/hickit/test_res/078-20260712_211825-p9016_corrected_combined_candidate_1m"
)


def require_analysis() -> None:
    if os.environ.get("CONDA_DEFAULT_ENV") != "analysis":
        raise RuntimeError("Experiment 058 must run inside conda environment 'analysis'")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    names = fields or list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def pairs_columns(path: Path) -> dict[str, int]:
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#columns:"):
                return {name: index for index, name in enumerate(line.split(":", 1)[1].strip().split())}
            if not line.startswith("#"):
                break
    raise ValueError(f"{path}: missing #columns header")


def read_real_observations(
    path: Path, max_observations: int | None, seed: int,
) -> tuple[list[Observation], dict[str, object]]:
    """Read raw-safe fields only; phase columns are deliberately ignored."""
    columns = pairs_columns(path)
    required = {"chr1", "pos1", "chr2", "pos2"}
    if not required.issubset(columns):
        raise ValueError(f"{path}: missing required endpoint fields")
    read_column = columns.get("readID", columns.get("read_id"))
    observations: list[Observation] = []
    rng = np.random.default_rng(seed)
    n_nonempty_read_id = 0
    n_source = 0
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            read_id = fields[read_column] if read_column is not None else "."
            molecule = None if read_id in {"", ".", "NA", "-"} else read_id
            n_nonempty_read_id += molecule is not None
            observation = Observation(
                observation_id=f"P9016:{n_source}",
                chrom1=fields[columns["chr1"]],
                pos1=int(fields[columns["pos1"]]),
                chrom2=fields[columns["chr2"]],
                pos2=int(fields[columns["pos2"]]),
                molecule_id=molecule,
            )
            n_source += 1
            if max_observations is None or len(observations) < max_observations:
                observations.append(observation)
            else:
                replace_index = int(rng.integers(0, n_source))
                if replace_index < max_observations:
                    observations[replace_index] = observation
    return observations, {
        "n_source_observations_scanned": n_source,
        "n_observations_retained": len(observations),
        "n_nonempty_read_or_molecule_ids": n_nonempty_read_id,
        "source_scan_is_complete": 1,
        "retained_sample_is_complete": int(max_observations is None),
        "sampling_policy": "all" if max_observations is None else "deterministic_reservoir",
        "phase_columns_read": 0,
    }


def synthetic_observations(dataset) -> list[Observation]:
    rows: list[Observation] = []
    for index, (left, right, molecule) in enumerate(
        zip(dataset.pair_i, dataset.pair_j, dataset.molecule_id, strict=True)
    ):
        rows.append(Observation(
            observation_id=f"synthetic:{index}",
            chrom1=f"chr{int(dataset.chrom_index[left]) + 1}",
            pos1=int(dataset.bin_offset[left]) * dataset.condition.resolution_bp,
            chrom2=f"chr{int(dataset.chrom_index[right]) + 1}",
            pos2=int(dataset.bin_offset[right]) * dataset.condition.resolution_bp,
            molecule_id=(f"mol:{int(molecule)}" if molecule >= 0 else None),
        ))
    return rows


def split_rows(
    observations: list[Observation], *, source: str, seed: int, block_sizes: list[int],
) -> tuple[list[dict[str, object]], dict[str, dict[str, str]]]:
    rows: list[dict[str, object]] = []
    assignments: dict[str, dict[str, str]] = {}
    specifications = [(mode, 5_000_000) for mode in SPLIT_MODES if mode != "genomic_block"]
    specifications.extend(("genomic_block", block_size) for block_size in block_sizes)
    for mode, block_size in specifications:
        label = mode if mode != "genomic_block" else f"genomic_block_{block_size // 1_000_000}Mb"
        assignment = blocked_split(
            observations, mode=mode, heldout_fraction=0.2, seed=seed,
            bin_size=1_000_000, block_size=block_size,
        )
        assignments[label] = assignment
        audit = leakage_audit(observations, assignment, bin_size=1_000_000, neighborhood_bp=block_size)
        train = sum(value == "train" for value in assignment.values())
        heldout = len(assignment) - train
        rows.append({
            "source": source,
            "split_mode": label,
            "block_size_bp": block_size if mode == "genomic_block" else "NA",
            "seed": seed,
            "n_observations": len(observations),
            "n_train": train,
            "n_heldout": heldout,
            "heldout_fraction_observed": heldout / max(1, len(observations)),
            **audit,
            "molecule_atomicity_enforced": 1,
        })
    return rows, assignments


def optimal_gauge_alignment(dataset, left: np.ndarray, right: np.ndarray) -> tuple[np.ndarray, dict[int, int], float]:
    chromosomes = range(dataset.condition.chromosome_count)
    best_agreement = -1.0
    best = {}
    best_q = right
    for bits in itertools.product((0, 1), repeat=max(0, dataset.condition.chromosome_count - 1)):
        flips = {0: 0, **{chrom: bit for chrom, bit in zip(range(1, dataset.condition.chromosome_count), bits, strict=True)}}
        aligned = np.vstack([
            relabel_probabilities(
                row,
                flip_left=flips[int(dataset.chrom_index[dataset.pair_i[index]])],
                flip_right=flips[int(dataset.chrom_index[dataset.pair_j[index]])],
            )
            for index, row in enumerate(right)
        ])
        agreement = float(np.mean(np.argmax(left, axis=1) == np.argmax(aligned, axis=1)))
        if agreement > best_agreement:
            best_agreement, best, best_q = agreement, flips, aligned
    return best_q, best, best_agreement


def run_synthetic_heldout(
    output: Path, *, quick: bool, seed: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    condition = SyntheticCondition(
        generator_family="model_mismatched", homolog_separation=1.25,
        contacts_per_cell=1200 if quick else 3000, cis_fraction=0.65,
        background_fraction=0.10, chromosome_count=3, chromosome_length_bins=10,
        readchain_length=3, molecule_count=40 if quick else 100, seed=seed,
    )
    dataset = generate_synthetic(condition)
    observations = synthetic_observations(dataset)
    manifest, assignments = split_rows(
        observations, source="synthetic_model_mismatched", seed=seed,
        block_sizes=[1_000_000, 5_000_000, 10_000_000],
    )
    offsets = [0, 12, 26] if quick else [0, 12, 26, 42, 54, 72, 84]
    seeds = [seed + offset for offset in offsets]
    heldout_rows: list[dict[str, object]] = []
    results_by_seed: dict[int, object] = {}
    posterior_by_seed: dict[int, np.ndarray] = {}
    selection_groups: dict[str, list[dict[str, object]]] = {}
    for split_label, assignment in assignments.items():
        train_indices = np.asarray([
            index for index, observation in enumerate(observations)
            if assignment[observation.observation_id] == "train"
        ], dtype=int)
        heldout_indices = np.asarray([
            index for index, observation in enumerate(observations)
            if assignment[observation.observation_id] == "heldout"
        ], dtype=int)
        if train_indices.size < 20 or heldout_indices.size < 10:
            continue
        train = subset_synthetic(dataset, train_indices)
        heldout = subset_synthetic(dataset, heldout_indices)
        group_rows: list[dict[str, object]] = []
        for run_seed in seeds:
            config = InferenceConfig(
                seed=run_seed, iterations=55 if quick else 140, learning_rate=0.025,
                likelihood_mode="monotone_log", staged_cis_then_trans=False,
            )
            result = infer_synthetic(train, config)
            q_heldout, log_marginal = score_coordinates(
                heldout, result.coordinates, likelihood_mode=config.likelihood_mode,
                alpha=config.alpha, d0=config.d0,
                observation_model="contrastive",
                signal_fraction=result.signal_fraction,
            )
            signal = ~heldout.is_background
            raw_metrics, _, _ = contact_metrics(heldout.truth_state[signal], q_heldout[signal])
            evaluated = evaluate_synthetic(train, result)
            # Truth geometry is opened only after blind held-out scores close.
            geometry_flips = geometry_selected_chromosome_flips(train, result.coordinates)
            aligned_q = relabel_contact_probabilities_to_geometry_gauge(
                heldout, q_heldout, geometry_flips
            )
            metrics, _, selective = contact_metrics(
                heldout.truth_state[signal], aligned_q[signal]
            )
            row: dict[str, object] = {
                "source": "synthetic_model_mismatched",
                "configuration": "blind_midpoint_delta_pairwise",
                "split_mode": split_label,
                "seed": run_seed,
                "n_train": int(train_indices.size),
                "n_heldout": int(heldout_indices.size),
                "n_heldout_signal": int(signal.sum()),
                "blind_heldout_pseudo_nll": float(-log_marginal.mean()),
                "blind_heldout_score_kernel_energy": float(-log_marginal.mean()),
                "blind_heldout_pmax_surprisal": float(-np.max(np.log(np.maximum(q_heldout, 1e-15)), axis=1).mean()),
                "raw_gauge_exact_top1_accuracy_eval_only": raw_metrics["exact_top1_accuracy"],
                "raw_gauge_same_cross_accuracy_eval_only": raw_metrics["same_cross_accuracy"],
                "exact_top1_accuracy_eval_only": metrics["exact_top1_accuracy"],
                "same_cross_accuracy_eval_only": metrics["same_cross_accuracy"],
                "exact_log_loss_eval_only": metrics["exact_log_loss"],
                "same_cross_log_loss_eval_only": metrics["same_cross_log_loss"],
                "exact_ece_eval_only": metrics["exact_ece"],
                "same_cross_ece_eval_only": metrics["same_cross_ece"],
                "unordered_structure_distance_spearman_eval_only": evaluated["unordered_structure_distance_spearman"],
                "runtime_seconds": result.runtime_seconds,
                "converged": int(result.converged),
                "truth_used_in_training_or_selection": "none",
                "contact_metric_gauge_policy": "ORACLE_GEOMETRY_ALIGNMENT_EVAL_ONLY",
                "heldout_fit_status": "RUN_SYNTHETIC_ONLY",
            }
            for selective_row in selective:
                metric = str(selective_row["metric"])
                target = float(selective_row["target_precision"])
                row[f"{metric}_coverage_at_precision_{target:.2f}_eval_only"] = selective_row["callable_fraction"]
            heldout_rows.append(row)
            group_rows.append(row)
            if split_label == "exact_haploid_bin_pair":
                q_all, _ = score_coordinates(
                    dataset, result.coordinates, observation_model="contrastive",
                    signal_fraction=result.signal_fraction,
                )
                results_by_seed[run_seed] = result
                posterior_by_seed[run_seed] = q_all
        selection_groups[split_label] = group_rows

    model_selection: list[dict[str, object]] = []
    blind_metrics = (
        ("blind_heldout_pseudo_nll", False),
        ("blind_heldout_pmax_surprisal", False),
    )
    truth_metrics = (
        "exact_top1_accuracy_eval_only",
        "same_cross_accuracy_eval_only",
        "unordered_structure_distance_spearman_eval_only",
    )
    for split_label, group in selection_groups.items():
        if len(group) < 2:
            continue
        for blind_name, higher_blind in blind_metrics:
            for truth_name in truth_metrics:
                values = model_selection_metrics(
                    [float(row[blind_name]) for row in group],
                    [float(row[truth_name]) for row in group],
                    higher_blind_is_better=higher_blind,
                )
                model_selection.append({
                    "source": "synthetic_model_mismatched",
                    "split_mode": split_label,
                    "blind_selection_metric": blind_name,
                    "truth_metric_eval_only": truth_name,
                    "n_candidates": len(group),
                    **values,
                    "truth_used_in_selection": "none",
                })

    stability: list[dict[str, object]] = []
    for left_seed, right_seed in itertools.combinations(sorted(results_by_seed), 2):
        left_result, right_result = results_by_seed[left_seed], results_by_seed[right_seed]
        left_q, right_q = posterior_by_seed[left_seed], posterior_by_seed[right_seed]
        aligned_q, flips, exact_agreement = optimal_gauge_alignment(dataset, left_q, right_q)
        aligned_coords = right_result.coordinates.copy()
        for bin_index, chrom in enumerate(dataset.chrom_index):
            if flips[int(chrom)]:
                aligned_coords[bin_index] = aligned_coords[bin_index, ::-1]
        parity_left = parity_probabilities(left_q)
        parity_right_raw = parity_probabilities(right_q)
        parity_right_aligned = parity_probabilities(aligned_q)
        parity_prediction_left = np.argmax(parity_left, axis=1)
        parity_prediction_right_raw = np.argmax(parity_right_raw, axis=1)
        parity_prediction_right_aligned = np.argmax(parity_right_aligned, axis=1)
        parity_raw_agreement = float(np.mean(parity_prediction_left == parity_prediction_right_raw))
        parity_aligned_agreement = float(np.mean(parity_prediction_left == parity_prediction_right_aligned))
        contact_cis = dataset.chrom_index[dataset.pair_i] == dataset.chrom_index[dataset.pair_j]
        raw_agreement = float(np.mean(np.argmax(left_q, axis=1) == np.argmax(right_q, axis=1)))
        per_chrom_structure = []
        for chrom in range(dataset.condition.chromosome_count):
            selected = dataset.chrom_index == chrom
            per_chrom_structure.append(float(structure_metrics(
                left_result.coordinates[selected], aligned_coords[selected]
            )["distance_spearman"]))
        callable_left = left_q.max(axis=1) >= 0.70
        callable_right = aligned_q.max(axis=1) >= 0.70
        stability.append({
            "configuration": "blind_midpoint_delta_pairwise",
            "split_mode": "exact_haploid_bin_pair",
            "seed_left": left_seed,
            "seed_right": right_seed,
            "same_cross_prediction_agreement_raw": parity_raw_agreement,
            "same_cross_prediction_agreement_after_optimal_chromosome_flips": parity_aligned_agreement,
            "same_cross_prediction_agreement_after_optimal_chromosome_flips_cis": float(np.mean(
                parity_prediction_left[contact_cis] == parity_prediction_right_aligned[contact_cis]
            )),
            "same_cross_prediction_agreement_after_optimal_chromosome_flips_trans": float(np.mean(
                parity_prediction_left[~contact_cis] == parity_prediction_right_aligned[~contact_cis]
            )),
            "exact_state_raw_agreement": raw_agreement,
            "exact_state_agreement_after_optimal_chromosome_flips": exact_agreement,
            "gauge_equivalent_agreement_gain": exact_agreement - raw_agreement,
            "distance_matrix_spearman": structure_metrics(
                left_result.coordinates, aligned_coords
            )["distance_spearman"],
            "mean_chromosome_structure_spearman": float(np.nanmean(per_chrom_structure)),
            "posterior_js_divergence_after_gauge_alignment": jensen_shannon_divergence(left_q, aligned_q),
            "fraction_contacts_callable_both_pmax_ge_0p70": float(np.mean(callable_left & callable_right)),
            "fraction_contacts_consistently_callable_same_state": float(np.mean(
                callable_left & callable_right & (np.argmax(left_q, axis=1) == np.argmax(aligned_q, axis=1))
            )),
            "optimal_flip_pattern": ",".join(f"chr{chrom + 1}:{flip}" for chrom, flip in sorted(flips.items())),
            "truth_used": "none",
        })
    return manifest, heldout_rows, stability, model_selection


def import_historical_078() -> list[dict[str, object]]:
    path = HISTORICAL_078 / "phase_metric_summary.tsv"
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as handle:
        return [{
            "source": "historical_078_registered_final_test",
            "configuration": row["arm"],
            "split_mode": "registered_exact_bpair_plus_source_unit",
            "seed": 17,
            "n_train": 1_379_753,
            "n_heldout": int(row["n_contacts"]),
            "n_heldout_signal": int(row["n_contacts"]),
            "blind_heldout_pseudo_nll": "NA_not_reported_as_normalized_likelihood",
            "blind_heldout_score_kernel_energy": "NA",
            "blind_heldout_pmax_surprisal": "NA",
            "exact_top1_accuracy_eval_only": float(row["top1_accuracy"]),
            "same_cross_accuracy_eval_only": float(row["same_cross_accuracy"]),
            "exact_log_loss_eval_only": float(row["four_state_nll"]),
            "same_cross_log_loss_eval_only": "NA",
            "exact_ece_eval_only": float(row["ece"]),
            "same_cross_ece_eval_only": "NA",
            "unordered_structure_distance_spearman_eval_only": "see_summary_geometry_selected",
            "runtime_seconds": "NA",
            "converged": 1,
            "truth_used_in_training_or_selection": "none; phase opened postclosure for evaluation",
            "heldout_fit_status": "HISTORICAL_EVAL_IMPORT; NORMALIZED_BLIND_HELDOUT_SCORE_NOT_AVAILABLE",
            "historical_scope": row["scope"],
            "pmax_ge_0p9_accuracy_eval_only": float(row["callable_accuracy_pmax_ge_0p9"]),
            "pmax_ge_0p9_recall_eval_only": float(row["callable_recall_pmax_ge_0p9"]),
        } for row in csv.DictReader(handle, delimiter="\t")]


def summarize_heldout(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Replicate mean/SD for the synthetic held-out ensemble only."""
    synthetic = [row for row in rows if row["source"] == "synthetic_model_mismatched"]
    metrics = (
        "blind_heldout_pseudo_nll",
        "blind_heldout_score_kernel_energy",
        "blind_heldout_pmax_surprisal",
        "raw_gauge_exact_top1_accuracy_eval_only",
        "raw_gauge_same_cross_accuracy_eval_only",
        "exact_top1_accuracy_eval_only",
        "same_cross_accuracy_eval_only",
        "exact_ece_eval_only",
        "same_cross_ece_eval_only",
        "unordered_structure_distance_spearman_eval_only",
    )
    summary: list[dict[str, object]] = []
    for split_mode in sorted({str(row["split_mode"]) for row in synthetic}):
        selected = [row for row in synthetic if row["split_mode"] == split_mode]
        output: dict[str, object] = {
            "source": "synthetic_model_mismatched",
            "split_mode": split_mode,
            "n_seeds": len(selected),
            "n_cells": 1,
            "n_chromosomes": 3,
            "truth_used_in_training_or_selection": "none",
            "uncertainty_definition": "sample standard deviation across fixed seeds",
        }
        for metric in metrics:
            values = np.asarray([float(row[metric]) for row in selected], dtype=float)
            output[f"{metric}_mean"] = float(np.nanmean(values))
            output[f"{metric}_sd"] = float(np.nanstd(values, ddof=1)) if values.size > 1 else math.nan
        summary.append(output)
    return summary


def plot_results(output: Path, heldout: list[dict[str, object]], stability: list[dict[str, object]]) -> None:
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    synthetic = [row for row in heldout if row["source"] == "synthetic_model_mismatched"]
    labels = sorted(set(str(row["split_mode"]) for row in synthetic))
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.0), constrained_layout=True)
    for axis, metric, title in (
        (axes[0], "exact_top1_accuracy_eval_only", "Exact four-state (geometry gauge aligned)"),
        (axes[1], "same_cross_accuracy_eval_only", "Same/cross (geometry gauge aligned)"),
    ):
        values = [[float(row[metric]) for row in synthetic if row["split_mode"] == label] for label in labels]
        axis.boxplot(values, tick_labels=labels, showfliers=False)
        axis.set_title(title)
        axis.set_ylabel("Held-out accuracy")
        axis.tick_params(axis="x", rotation=55)
        axis.axhline(0.25 if "exact" in metric else 0.5, color="0.5", linestyle="--", linewidth=0.8)
    fig.savefig(figures / "blocked_split_accuracy.pdf")
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(3.4, 3.0), constrained_layout=True)
    x = [float(row["exact_state_agreement_after_optimal_chromosome_flips"]) for row in stability]
    y = [float(row["same_cross_prediction_agreement_after_optimal_chromosome_flips"]) for row in stability]
    axis.scatter(x, y, s=24, color="#0072B2")
    axis.set_xlabel("Exact agreement after chromosome flips")
    axis.set_ylabel("Same/cross agreement")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    fig.savefig(figures / "seed_stability.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--output-dir", type=Path, default=Path("results/phase_diagnostics/058"))
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    require_analysis()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    real, real_info = read_real_observations(args.pairs, 40_000 if args.quick else None, args.seed)
    real_manifest, _ = split_rows(
        real, source="P9016_raw_safe_fields_only", seed=args.seed,
        block_sizes=[1_000_000, 5_000_000, 10_000_000],
    )
    for row in real_manifest:
        row.update(real_info)
        row["input_sha256"] = file_sha256(args.pairs)
        row["readID_observation"] = "P9016 readID is '.' in inspected rows; no recoverable molecule grouping"
        row["molecule_atomicity_enforced"] = "NA_NOT_AUDITABLE_NO_MOLECULE_IDS"
        row["molecule_leakage_groups"] = "NA_NOT_AUDITABLE_NO_MOLECULE_IDS"
        row["heldout_fit_status"] = "NOT_RUN_FOR_P9016; SPLIT_AND_LEAKAGE_AUDIT_ONLY"

    synthetic_manifest, heldout, stability, selection = run_synthetic_heldout(
        args.output_dir, quick=args.quick, seed=args.seed
    )
    heldout.extend(import_historical_078())
    write_tsv(args.output_dir / "split_manifest.tsv", real_manifest + synthetic_manifest)
    write_tsv(args.output_dir / "heldout_metrics.tsv", heldout)
    write_tsv(args.output_dir / "heldout_summary.tsv", summarize_heldout(heldout))
    write_tsv(args.output_dir / "seed_stability.tsv", stability)
    write_tsv(args.output_dir / "model_selection.tsv", selection)
    config_payload = {
        "logical_experiment_id": "058",
        "resolution_bp": 1_000_000,
        "seed": args.seed,
        "quick": args.quick,
        "P9016_input": str(args.pairs),
        "P9016_fit_status": "NOT_RUN; split/leakage audit only",
        "synthetic_seed_offsets": [0, 12, 26] if args.quick else [0, 12, 26, 42, 54, 72, 84],
        "split_modes": list(SPLIT_MODES),
        "genomic_block_sizes_bp": [1_000_000, 5_000_000, 10_000_000],
        "synthetic_observation_score": "contrastive local four-state pseudo-likelihood",
        "truth_boundary": "truth metrics opened only after each blind score exists",
    }
    (args.output_dir / "config.json").write_text(
        json.dumps(config_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    plot_results(args.output_dir, heldout, stability)
    provenance = collect_provenance(
        Path(__file__).resolve().parents[1], sys.argv, args.seed,
        {
            **config_payload, "experiment": "058", "mode": "quick" if args.quick else "full", **real_info,
            "P9016_heldout_fit_status": "NOT_RUN; synthetic diagnostic only",
            "historical_078_sha256": file_sha256(HISTORICAL_078 / "phase_metric_summary.tsv"),
        },
        (args.pairs, HISTORICAL_078 / "phase_metric_summary.tsv"),
    )
    write_provenance(args.output_dir / "provenance.json", provenance)


if __name__ == "__main__":
    main()
