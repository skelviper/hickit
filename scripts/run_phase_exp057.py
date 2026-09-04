#!/usr/bin/env python3
"""Logical Experiment 057: ORACLE decomposition of E-step and M-step failure."""

from __future__ import annotations

import argparse
import csv
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

from phase_diagnostics.likelihood import SCORE_MODES, ScoreConfig, log_score
from phase_diagnostics.metrics import contact_metrics, selective_threshold_curve
from phase_diagnostics.oracle_mstep import run_native_oracle_mstep
from phase_diagnostics.p9016 import PhaseContactData, load_phase_contacts
from phase_diagnostics.provenance import collect_provenance, file_sha256, write_provenance
from phase_diagnostics.synthetic import (
    InferenceConfig,
    SyntheticCondition,
    evaluate_synthetic,
    generate_synthetic,
    infer_synthetic,
)


DEFAULT_PAIRS = Path("/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz")
DEFAULT_REFERENCE = Path("/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz")
HISTORICAL_ORACLE = Path(
    "/work/phase3/hickit/test_res/083-20260713_111134-p9016_phased_hickit_oracle_gate_1m"
)
PHASE_KNOWN_TRAIN = HISTORICAL_ORACLE / "inputs/train.true_phase.onehot.pairs.gz"


def require_analysis() -> None:
    if os.environ.get("CONDA_DEFAULT_ENV") != "analysis":
        raise RuntimeError("Experiment 057 must run inside conda environment 'analysis'")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    names = fields or list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def posterior_matrix(data: PhaseContactData, config: ScoreConfig) -> np.ndarray:
    # Current Hickit stores base_d_scale=n_raw^(-1/3). Normalize distances by it.
    normalized_distance = data.distances * np.cbrt(data.multiplicity[:, None])
    scores = log_score(normalized_distance.ravel(), config).reshape(-1, 4)
    scores -= scores.max(axis=1, keepdims=True)
    values = np.exp(np.clip(scores, -745, 0))
    return values / values.sum(axis=1, keepdims=True)


def append_metrics(
    long_rows: list[dict[str, object]],
    calibration_rows: list[dict[str, object]],
    *,
    component: str,
    configuration: str,
    scope: str,
    truth_used: str,
    truth: np.ndarray,
    probabilities: np.ndarray,
) -> dict[str, float | int]:
    metrics, calibration, selective = contact_metrics(truth, probabilities)
    for metric, value in metrics.items():
        long_rows.append({
            "experiment": "057",
            "component": component,
            "configuration": configuration,
            "scope": scope,
            "metric": metric,
            "value": value,
            "n": truth.size,
            "truth_used_in_inference": truth_used,
        })
    for row in calibration:
        calibration_rows.append({
            "experiment": "057",
            "component": component,
            "configuration": configuration,
            "scope": scope,
            **asdict(row),
        })
    for row in selective:
        long_rows.append({
            "experiment": "057",
            "component": component,
            "configuration": configuration,
            "scope": scope,
            "metric": f"{row['metric']}_coverage_at_precision_{float(row['target_precision']):.2f}",
            "value": row["callable_fraction"],
            "n": truth.size,
            "truth_used_in_inference": truth_used,
        })
    return metrics


def genomic_strata(data: PhaseContactData) -> np.ndarray:
    result = np.full(data.n_contacts, "trans", dtype=object)
    distance = data.genomic_separation
    result[data.is_cis & (distance < 5_000_000)] = "cis_lt5Mb"
    result[data.is_cis & (distance >= 5_000_000) & (distance < 20_000_000)] = "cis_5_20Mb"
    result[data.is_cis & (distance >= 20_000_000) & (distance < 50_000_000)] = "cis_20_50Mb"
    result[data.is_cis & (distance >= 50_000_000)] = "cis_ge50Mb"
    return result


def margin_strata(data: PhaseContactData) -> np.ndarray:
    truth_distance = data.distances[np.arange(data.n_contacts), data.truth_state]
    other = data.distances.copy()
    other[np.arange(data.n_contacts), data.truth_state] = np.inf
    margin = other.min(axis=1) - truth_distance
    positive = margin[margin > 0]
    q1, q2 = np.quantile(positive, (1 / 3, 2 / 3)) if positive.size else (0.0, 0.0)
    result = np.full(data.n_contacts, "truth_not_nearest", dtype=object)
    result[(margin > 0) & (margin <= q1)] = "truth_nearest_low_margin"
    result[(margin > q1) & (margin <= q2)] = "truth_nearest_mid_margin"
    result[margin > q2] = "truth_nearest_high_margin"
    return result


def distance_diagnostics(data: PhaseContactData, probabilities: np.ndarray) -> dict[str, float]:
    order = np.argsort(data.distances, axis=1)
    truth_rank = np.empty(data.n_contacts, dtype=int)
    for rank in range(4):
        selected = order[:, rank]
        truth_rank[selected == data.truth_state] = rank + 1
    nearest = order[:, 0]
    top = np.argmax(probabilities, axis=1)
    return {
        "mean_truth_state_distance_rank": float(truth_rank.mean()),
        "truth_state_geometrically_nearest_fraction": float(np.mean(nearest == data.truth_state)),
        "posterior_top1_equals_geometrically_nearest_fraction": float(np.mean(top == nearest)),
    }


def prior_baselines(data: PhaseContactData, rng: np.random.Generator) -> dict[str, np.ndarray]:
    uniform = np.full((data.n_contacts, 4), 0.25)
    cis_parity = uniform.copy()
    cis_parity[data.is_cis] = np.asarray([0.45, 0.05, 0.05, 0.45])
    random_prediction = rng.dirichlet(np.ones(4), size=data.n_contacts)
    return {
        "prior_uniform_four_state": uniform,
        "prior_cis_same_cross_configured_0p9": cis_parity,
        "prior_chromosome_pair_unlabeled_gauge_symmetric": uniform.copy(),
        "prior_state_frequency_unlabeled_max_entropy": uniform.copy(),
        "prior_random_seeded": random_prediction,
    }


def shuffled_within_strata(data: PhaseContactData, rng: np.random.Generator) -> np.ndarray:
    result = data.distances.copy()
    labels = genomic_strata(data)
    trans_pair = np.asarray([
        f"{min(str(a), str(b))}|{max(str(a), str(b))}" if not cis else stratum
        for a, b, cis, stratum in zip(data.chrom1, data.chrom2, data.is_cis, labels, strict=True)
    ], dtype=object)
    for label in sorted(set(trans_pair.tolist())):
        selected = np.flatnonzero(trans_pair == label)
        if selected.size > 1:
            result[selected] = result[rng.permutation(selected)]
    return result


def run_estep(
    data: PhaseContactData,
    output: Path,
    seed: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], dict[str, np.ndarray]]:
    long_rows: list[dict[str, object]] = []
    calibration_rows: list[dict[str, object]] = []
    per_chromosome: list[dict[str, object]] = []
    posterior_by_config: dict[str, np.ndarray] = {}
    rng = np.random.default_rng(seed)
    strata_sets = {
        "all": np.ones(data.n_contacts, dtype=bool),
        "cis": data.is_cis,
        "trans": ~data.is_cis,
    }
    for mode in SCORE_MODES:
        config = ScoreConfig(mode=mode, temperature=1.0, k=1.0, alpha=2.0, d0=0.05, r0=1.0, scale=0.25)
        posterior = posterior_matrix(data, config)
        posterior_by_config[f"oracle_coords_{mode}"] = posterior
        for scope, selected in strata_sets.items():
            if not selected.any():
                continue
            append_metrics(
                long_rows, calibration_rows, component="057A_oracle_coordinates_score_kernel_isolation",
                configuration=f"oracle_coords_{mode}", scope=scope, truth_used="ORACLE_COORDINATES",
                truth=data.truth_state[selected], probabilities=posterior[selected],
            )
            for metric, value in distance_diagnostics(
                PhaseContactData(
                    chrom1=data.chrom1[selected], start1=data.start1[selected], chrom2=data.chrom2[selected],
                    start2=data.start2[selected], truth_state=data.truth_state[selected], distances=data.distances[selected],
                    multiplicity=data.multiplicity[selected], genomic_separation=data.genomic_separation[selected],
                    is_cis=data.is_cis[selected],
                ), posterior[selected]
            ).items():
                long_rows.append({
                    "experiment": "057", "component": "057A_oracle_coordinates_score_kernel_isolation",
                    "configuration": f"oracle_coords_{mode}", "scope": scope, "metric": metric,
                    "value": value, "n": int(selected.sum()), "truth_used_in_inference": "ORACLE_COORDINATES",
                })
        for chrom in sorted(set(data.chrom1.tolist()) | set(data.chrom2.tolist())):
            selected = (data.chrom1 == chrom) | (data.chrom2 == chrom)
            metrics, _, _ = contact_metrics(data.truth_state[selected], posterior[selected])
            per_chromosome.append({
                "component": "057A", "configuration": f"oracle_coords_{mode}", "chromosome": chrom,
                "n_contacts": int(selected.sum()), "exact_top1_accuracy": metrics["exact_top1_accuracy"],
                "same_cross_accuracy": metrics["same_cross_accuracy"], "exact_log_loss": metrics["exact_log_loss"],
                "same_cross_log_loss": metrics["same_cross_log_loss"], "structure_metric": "NA",
                "truth_used_in_inference": "ORACLE_COORDINATES",
            })

    primary = posterior_by_config["oracle_coords_monotone_log"]
    for labels, family in ((genomic_strata(data), "genomic_separation"), (margin_strata(data), "truth_distance_margin")):
        for label in sorted(set(labels.tolist())):
            selected = labels == label
            append_metrics(
                long_rows, calibration_rows, component=f"057A_stratified_{family}",
                configuration="oracle_coords_monotone_log", scope=str(label), truth_used="ORACLE_COORDINATES",
                truth=data.truth_state[selected], probabilities=primary[selected],
            )

    for name, posterior in prior_baselines(data, rng).items():
        posterior_by_config[name] = posterior
        for scope, selected in strata_sets.items():
            append_metrics(
                long_rows, calibration_rows, component="057D_prior_only", configuration=name, scope=scope,
                truth_used="none", truth=data.truth_state[selected], probabilities=posterior[selected],
            )

    normalized = data.distances * np.cbrt(data.multiplicity[:, None])
    control_distances = {
        "negative_shuffled_endpoints_stratum_preserving": shuffled_within_strata(data, rng),
        "negative_randomized_geometry": rng.lognormal(mean=0.0, sigma=1.0, size=data.distances.shape),
        "negative_permuted_chromosome_labels_trans": data.distances.copy(),
    }
    trans = np.flatnonzero(~data.is_cis)
    control_distances["negative_permuted_chromosome_labels_trans"][trans] = data.distances[rng.permutation(trans)]
    for name, distances in control_distances.items():
        score = log_score((distances * np.cbrt(data.multiplicity[:, None])).ravel(), ScoreConfig(mode="monotone_log")).reshape(-1, 4)
        score -= score.max(axis=1, keepdims=True)
        posterior = np.exp(score)
        posterior /= posterior.sum(axis=1, keepdims=True)
        posterior_by_config[name] = posterior
        for scope, selected in strata_sets.items():
            append_metrics(
                long_rows, calibration_rows, component="057E_randomized_negative_control",
                configuration=name, scope=scope, truth_used="ORACLE_COORDINATES_RANDOMIZED_CONTROL",
                truth=data.truth_state[selected], probabilities=posterior[selected],
            )
    random_p4 = rng.dirichlet(np.ones(4), size=data.n_contacts)
    posterior_by_config["negative_randomized_p4_initialization"] = random_p4
    for scope, selected in strata_sets.items():
        append_metrics(
            long_rows, calibration_rows, component="057E_randomized_negative_control",
            configuration="negative_randomized_p4_initialization", scope=scope, truth_used="none",
            truth=data.truth_state[selected], probabilities=random_p4[selected],
        )
    return long_rows, calibration_rows, per_chromosome, posterior_by_config


def selective_calling_rows(
    data: PhaseContactData, posterior_by_config: dict[str, np.ndarray],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    configurations = (
        "oracle_coords_fdg_flat",
        "oracle_coords_monotone_log",
        "prior_uniform_four_state",
        "negative_randomized_geometry",
    )
    scopes = {
        "all": np.ones(data.n_contacts, dtype=bool),
        "cis": data.is_cis,
        "trans": ~data.is_cis,
    }
    for configuration in configurations:
        probability = posterior_by_config[configuration]
        for scope, selected in scopes.items():
            truth = data.truth_state[selected]
            conditional = probability[selected]
            conditional /= conditional.sum(axis=1, keepdims=True)
            exact_prediction = np.argmax(conditional, axis=1)
            exact_confidence = np.max(conditional, axis=1)
            parity_probability = np.column_stack((
                conditional[:, 0] + conditional[:, 3],
                conditional[:, 1] + conditional[:, 2],
            ))
            parity_truth = np.isin(truth, (1, 2)).astype(int)
            for target, correct, confidence in (
                ("exact", exact_prediction == truth, exact_confidence),
                (
                    "same_cross",
                    np.argmax(parity_probability, axis=1) == parity_truth,
                    np.max(parity_probability, axis=1),
                ),
            ):
                for row in selective_threshold_curve(correct, confidence):
                    rows.append({
                        "experiment": "057",
                        "configuration": configuration,
                        "scope": scope,
                        "target": target,
                        **row,
                        "truth_used_in_inference": (
                            "ORACLE_COORDINATES"
                            if configuration.startswith("oracle_coords_") else "none"
                        ),
                    })
    return rows


def import_oracle_mstep(
    long_rows: list[dict[str, object]], per_chromosome: list[dict[str, object]],
) -> None:
    summary_path = HISTORICAL_ORACLE / "strict_eval/summary.tsv"
    per_chrom_path = HISTORICAL_ORACLE / "strict_eval/per_chromosome.tsv"
    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row["arm"].startswith("true_seed"):
                continue
            for metric in (
                "strict_unordered_pair_spearman_mean", "strict_worst_copy_spearman_mean",
                "assignment_gap_mean", "d0_minus_d1_abs_spearman_mean",
                "pooled_two_copy_spearman_mean_diagnostic",
            ):
                long_rows.append({
                    "experiment": "057", "component": "057B_historical_joint_oracle_phase_hickit_partial",
                    "configuration": row["arm"], "scope": "all_chromosomes", "metric": metric,
                    "value": row[metric], "n": 20, "truth_used_in_inference": "ORACLE_CONTACT_PHASE",
                })
    with per_chrom_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row["arm"].startswith("true_seed"):
                continue
            per_chromosome.append({
                "component": "057B_historical_joint_partial", "configuration": row["arm"], "chromosome": row["chrom"],
                "n_contacts": "NA", "exact_top1_accuracy": "NA", "same_cross_accuracy": "NA",
                "exact_log_loss": "NA", "same_cross_log_loss": "NA",
                "structure_metric": row["strict_unordered_pair_spearman"],
                "truth_used_in_inference": "ORACLE_CONTACT_PHASE",
            })
def run_near_truth(
    long_rows: list[dict[str, object]], output: Path, *, quick: bool, base_seed: int,
) -> list[dict[str, object]]:
    levels = (0.0, 0.10, 0.50, 1.00, 2.00) if quick else (0.0, 0.05, 0.10, 0.25, 0.50, 1.00, 2.00)
    seeds = (base_seed, base_seed + 12) if quick else (base_seed, base_seed + 12, base_seed + 26)
    rows: list[dict[str, object]] = []
    for seed in seeds:
        dataset = generate_synthetic(SyntheticCondition(
            generator_family="well_specified", homolog_separation=1.5,
            contacts_per_cell=1000 if quick else 2500, chromosome_count=2,
            chromosome_length_bins=8, background_fraction=0.05, seed=seed,
        ))
        adjacent = []
        for chrom in range(dataset.condition.chromosome_count):
            selected = np.flatnonzero(dataset.chrom_index == chrom)
            adjacent.extend(np.linalg.norm(
                dataset.truth_coordinates[selected[1:]] - dataset.truth_coordinates[selected[:-1]], axis=-1
            ).ravel())
        scale = float(np.median(adjacent))
        for level in levels:
            rng = np.random.default_rng(seed * 1009 + int(level * 1000))
            initial = dataset.truth_coordinates + rng.normal(scale=level * scale, size=dataset.truth_coordinates.shape)
            result = infer_synthetic(
                dataset,
                InferenceConfig(seed=seed, iterations=70 if quick else 140, learning_rate=0.02),
                initial_coordinates=initial,
            )
            metrics = evaluate_synthetic(dataset, result)
            row = {"noise_level_adjacent_median": level, "seed": seed, **metrics}
            rows.append(row)
            for metric in (
                "exact_top1_accuracy", "same_cross_accuracy",
                "geometry_gauge_aligned_exact_top1_accuracy",
                "geometry_gauge_aligned_same_cross_accuracy",
                "unordered_structure_rmsd",
                "unordered_structure_distance_spearman", "mean_exact_confidence", "mean_entropy",
            ):
                long_rows.append({
                    "experiment": "057", "component": "057C_prototype_near_truth_initialization",
                    "configuration": f"noise_{level:g}_seed_{seed}", "scope": "synthetic_all",
                    "metric": metric, "value": metrics[metric], "n": metrics["n_contacts"],
                    "truth_used_in_inference": "ORACLE_NEAR_TRUTH_INITIALIZATION",
                })
    write_tsv(output / "near_truth_initialization.tsv", rows)
    return rows


def plot_results(
    output: Path, long_rows: list[dict[str, object]], near_truth: list[dict[str, object]],
    selective_rows: list[dict[str, object]],
) -> None:
    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7})
    configs = [f"oracle_coords_{mode}" for mode in SCORE_MODES]
    figure, axes = plt.subplots(1, 2, figsize=(6.4, 3.0), constrained_layout=True)
    for axis, metric, title in zip(
        axes, ("exact_top1_accuracy", "same_cross_accuracy"), ("Exact four-state", "Same/cross"), strict=True
    ):
        cis = [float(next(row["value"] for row in long_rows if row["configuration"] == config and row["scope"] == "cis" and row["metric"] == metric)) for config in configs]
        trans = [float(next(row["value"] for row in long_rows if row["configuration"] == config and row["scope"] == "trans" and row["metric"] == metric)) for config in configs]
        x = np.arange(len(configs))
        axis.bar(x - 0.18, cis, width=0.36, color="#0072B2", label="cis")
        axis.bar(x + 0.18, trans, width=0.36, color="#D55E00", label="trans")
        axis.axhline(0.25 if metric.startswith("exact") else 0.5, color="black", linewidth=0.6, linestyle=":")
        axis.set_xticks(x, [item.replace("oracle_coords_", "") for item in configs], rotation=35, ha="right")
        axis.set_ylim(0, 1)
        axis.set_ylabel("accuracy")
        axis.set_title(title)
        axis.legend(frameon=False)
    figure.savefig(figure_dir / "oracle_coordinate_estep_accuracy.pdf", dpi=300)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(6.4, 3.0), constrained_layout=True)
    seed_values = sorted({int(row["seed"]) for row in near_truth})
    colors = ("#0072B2", "#D55E00", "#009E73", "#CC79A7")
    for seed, color in zip(seed_values, colors, strict=False):
        selected = sorted((row for row in near_truth if row["seed"] == seed), key=lambda row: float(row["noise_level_adjacent_median"]))
        if not selected:
            continue
        x = [float(row["noise_level_adjacent_median"]) for row in selected]
        axes[0].plot(
            x,
            [float(row["geometry_gauge_aligned_same_cross_accuracy"]) for row in selected],
            marker="o", ms=2.5, color=color, label=f"seed {seed}",
        )
        axes[1].plot(x, [float(row["unordered_structure_distance_spearman"]) for row in selected], marker="o", ms=2.5, color=color, label=f"seed {seed}")
    axes[0].axhline(0.5, color="black", linewidth=0.6, linestyle=":")
    axes[0].set_ylabel("same/cross accuracy (geometry gauge aligned)")
    axes[1].set_ylabel("unordered distance Spearman")
    for axis in axes:
        axis.set_xlabel("initial noise / median adjacent distance")
        axis.legend(frameon=False)
    figure.savefig(figure_dir / "near_truth_initialization.pdf", dpi=300)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(6.4, 3.0), constrained_layout=True)
    styles = {
        "oracle_coords_fdg_flat": ("#D55E00", "FDG-flat"),
        "oracle_coords_monotone_log": ("#0072B2", "monotone-log"),
        "prior_uniform_four_state": ("#666666", "uniform prior"),
    }
    for axis, target in zip(axes, ("exact", "same_cross"), strict=True):
        for configuration, (color, label) in styles.items():
            selected = [
                row for row in selective_rows
                if row["configuration"] == configuration and row["scope"] == "all"
                and row["target"] == target and int(row["callable_count"]) > 0
            ]
            selected.sort(key=lambda row: float(row["coverage"]))
            axis.plot(
                [float(row["coverage"]) for row in selected],
                [float(row["precision"]) for row in selected],
                marker="o", markersize=2, linewidth=0.8, color=color, label=label,
            )
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.set_xlabel("callable coverage")
        axis.set_ylabel("accuracy / precision")
        axis.set_title("Exact state" if target == "exact" else "Same/cross")
        axis.legend(frameon=False)
    figure.savefig(figure_dir / "selective_calling.pdf", dpi=300)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--reference-3dg", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output-dir", type=Path, default=Path("results/phase_diagnostics/057"))
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    require_analysis()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "configs").mkdir(exist_ok=True)
    config_payload = {
        "logical_experiment_id": "057", "formal_sequence_alias": "assigned_at_packaging",
        "resolution_bp": 1_000_000, "pairs": str(args.pairs), "reference_3dg": str(args.reference_3dg),
        "seed": args.seed, "quick": args.quick, "score_modes": list(SCORE_MODES),
        "oracle_boundary": "phase and CHARM/3DG enter only this explicitly ORACLE experiment",
        "057A_scope": "ORACLE-coordinate score-kernel isolation, not the full production C E-step",
        "057A_omitted_terms": "production neighbor-derived base_k/unit, dynamic priors, and rate weights are unavailable on frozen reference coordinates",
        "057B_scope": (
            "quick mode imports PARTIAL historical joint evidence; full mode runs current CPU Hickit "
            "cis-only and cis+trans plus an experimental fixed-cis rigid trans-placement arm"
        ),
        "057C_scope": "Python midpoint/Delta prototype, not the production C full inference loop",
    }
    (args.output_dir / "configs/experiment_057.json").write_text(
        json.dumps(config_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    data = load_phase_contacts(
        args.pairs, args.reference_3dg, resolution=1_000_000,
        max_contacts=50_000 if args.quick else None,
    )
    long_rows, calibration_rows, per_chromosome, posterior_by_config = run_estep(
        data, args.output_dir, args.seed
    )
    import_oracle_mstep(long_rows, per_chromosome)
    native_manifest: dict[str, object] | None = None
    if not args.quick:
        native_rows, native_chromosome_rows, native_manifest = run_native_oracle_mstep(
            Path(__file__).resolve().parents[1], args.output_dir, data,
            args.reference_3dg, PHASE_KNOWN_TRAIN,
            (args.seed, args.seed + 12, args.seed + 26),
        )
        write_tsv(args.output_dir / "oracle_mstep_runs.tsv", native_rows)
        per_chromosome.extend(native_chromosome_rows)
        metadata = {"arm", "seed", "truth_used_in_inference"}
        for row in native_rows:
            for metric, value in row.items():
                if metric in metadata or isinstance(value, str):
                    continue
                long_rows.append({
                    "experiment": "057",
                    "component": "057B_explicit_oracle_phase_mstep",
                    "configuration": f"{row['arm']}_seed{row['seed']}",
                    "scope": "all_chromosomes",
                    "metric": metric,
                    "value": value,
                    "n": row["n_chromosomes"],
                    "truth_used_in_inference": "ORACLE_CONTACT_PHASE",
                })
    near_truth = run_near_truth(long_rows, args.output_dir, quick=args.quick, base_seed=args.seed)
    selective_rows = selective_calling_rows(data, posterior_by_config)
    write_tsv(args.output_dir / "metrics_long.tsv", long_rows)
    write_tsv(args.output_dir / "per_chromosome.tsv", per_chromosome)
    write_tsv(args.output_dir / "calibration.tsv", calibration_rows)
    write_tsv(args.output_dir / "selective_calling.tsv", selective_rows)
    manifest_rows = [
        {"field": "logical_experiment_id", "value": "057", "notes": "requested scientific label"},
        {"field": "formal_test_res_sequence", "value": "assigned_at_packaging", "notes": "packager creates a new non-overwriting sequence"},
        {"field": "resolution_bp", "value": 1_000_000, "notes": "P9016 and synthetic final resolution"},
        {"field": "evaluated_phase_contacts", "value": data.n_contacts, "notes": "double-end phase plus coordinate overlap"},
        {"field": "pairs_sha256", "value": file_sha256(args.pairs), "notes": "blind source; phase read only here"},
        {"field": "reference_sha256", "value": file_sha256(args.reference_3dg), "notes": "ORACLE coordinates"},
        {"field": "oracle_mstep_source", "value": str(HISTORICAL_ORACLE), "notes": "hash-closed native Hickit phase-known experiment 083"},
        {"field": "oracle_mstep_source_summary_sha256", "value": file_sha256(HISTORICAL_ORACLE / "strict_eval/summary.tsv"), "notes": "imported evidence"},
        {"field": "057B_status", "value": native_manifest["status"] if native_manifest else "PARTIAL_QUICK_MODE", "notes": "full mode has three explicit arms; quick mode imports historical joint evidence only"},
        {"field": "057B_cis_only_arm", "value": "RUN_CURRENT_HICKIT_CPU" if native_manifest else "NOT_RUN_QUICK_MODE", "notes": "ORACLE phase; independent input containing cis contacts only"},
        {"field": "057B_cis_plus_trans_arm", "value": "RUN_CURRENT_HICKIT_CPU" if native_manifest else "HISTORICAL_JOINT_ONLY", "notes": "ORACLE phase"},
        {"field": "057B_fixed_cis_trans_placement_arm", "value": "RUN_EXPERIMENTAL_RIGID_BODY" if native_manifest else "NOT_RUN_QUICK_MODE", "notes": "cis coordinates immutable; separate from native Hickit"},
        {"field": "truth_used_in_blind_inference", "value": 0, "notes": "all truth-dependent rows explicitly ORACLE"},
        {"field": "quick", "value": int(args.quick), "notes": "bounded diagnostic mode"},
    ]
    if native_manifest is not None:
        manifest_rows.extend(
            {"field": f"057B_{key}", "value": value, "notes": "explicit ORACLE-phase M-step provenance"}
            for key, value in native_manifest.items()
            if key != "status"
        )
    write_tsv(args.output_dir / "run_manifest.tsv", manifest_rows, ["field", "value", "notes"])
    plot_results(args.output_dir, long_rows, near_truth, selective_rows)
    provenance = collect_provenance(
        Path(__file__).resolve().parents[1], sys.argv, args.seed, config_payload,
        (args.pairs, args.reference_3dg, PHASE_KNOWN_TRAIN),
    )
    write_provenance(args.output_dir / "provenance.json", provenance)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
