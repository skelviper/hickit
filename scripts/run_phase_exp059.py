#!/usr/bin/env python3
"""Logical Experiment 059: synthetic identifiability phase diagrams."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from phase_diagnostics.metrics import contact_metrics
from phase_diagnostics.provenance import collect_provenance, write_provenance
from phase_diagnostics.synthetic import (
    InferenceConfig,
    SyntheticCondition,
    evaluate_synthetic,
    generate_synthetic,
    infer_synthetic,
    oracle_coordinate_posteriors,
    oracle_generator_conditional_posteriors,
)


def require_analysis() -> None:
    if os.environ.get("CONDA_DEFAULT_ENV") != "analysis":
        raise RuntimeError("Experiment 059 must run inside conda environment 'analysis'")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    names = fields or list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def condition_grid(quick: bool, base_seed: int) -> list[dict[str, object]]:
    families = ("well_specified", "model_mismatched")
    rows: list[dict[str, object]] = []
    seen: set[tuple[object, ...]] = set()

    def add(family: str, diagram: str, *, grid_stage: str = "coarse", **updates: object) -> None:
        base = SyntheticCondition(
            generator_family=family, homolog_separation=1.0,
            contacts_per_cell=1200, cis_fraction=0.65, background_fraction=0.05,
            resolution_bp=1_000_000, chromosome_count=3, chromosome_length_bins=10,
            readchain_length=2, molecule_count=0,
            count_overdispersion=(0.0 if family == "well_specified" else 0.2),
            anchor_fraction=0.0, seed=base_seed,
        )
        condition = replace(base, **updates)
        key = tuple(asdict(condition).items())
        identity = (diagram, key)
        if identity in seen:
            return
        seen.add(identity)
        rows.append({"diagram": diagram, "condition": condition, "grid_stage": grid_stage})

    depths = (400, 1600) if quick else (300, 900, 3000, 12000)
    separations = (0.25, 1.5) if quick else (0.1, 0.75, 2.0)
    trans_fractions = (0.2, 0.8) if quick else (0.1, 0.5, 0.9)
    backgrounds = (0.0, 0.3) if quick else (0.0, 0.2, 0.5)
    anchor_depths = (1000,) if quick else (500, 2500, 10000)
    anchors = (0.0, 0.001, 0.005, 0.01, 0.05)
    chain_lengths = (2, 4) if quick else (2, 4, 6)
    molecule_counts = (0, 40) if quick else (0, 20, 80)
    for family in families:
        for depth in depths:
            for separation in separations:
                add(family, "contacts_x_homolog_separation", contacts_per_cell=depth, homolog_separation=separation)
        if not quick:
            # A preliminary coarse run localized the largest parity change to
            # this neighborhood. Refine it without opening a full factorial grid.
            for depth in (1500, 6000):
                for separation in (1.25, 1.6):
                    add(
                        family, "contacts_x_homolog_separation",
                        grid_stage="transition_refinement",
                        contacts_per_cell=depth, homolog_separation=separation,
                    )
        for trans_fraction in trans_fractions:
            for background in backgrounds:
                add(family, "trans_fraction_x_background", cis_fraction=1.0 - trans_fraction, background_fraction=background)
        for anchor in anchors:
            for depth in anchor_depths:
                add(family, "anchor_fraction_x_contacts", anchor_fraction=anchor, contacts_per_cell=depth)
        for chain_length in chain_lengths:
            for molecule_count in molecule_counts:
                budget = max(1200, molecule_count * chain_length * (chain_length - 1) // 2)
                add(
                    family, "readchain_length_x_molecule_count", readchain_length=chain_length,
                    molecule_count=molecule_count, contacts_per_cell=budget,
                )
    return rows


def condition_identifier(diagram: str, condition: SyntheticCondition) -> str:
    return (
        f"{diagram}__{condition.generator_family}__d{condition.contacts_per_cell}"
        f"__sep{condition.homolog_separation:g}__cis{condition.cis_fraction:g}"
        f"__bg{condition.background_fraction:g}__a{condition.anchor_fraction:g}"
        f"__r{condition.readchain_length}__m{condition.molecule_count}"
    )


def condition_replicate_seed(
    base_seed: int, replicate: int, diagram: str, condition: SyntheticCondition,
) -> int:
    """Derive stable seeds while sharing data across nested anchor fractions."""
    fields = asdict(condition)
    fields.pop("seed")
    if diagram == "anchor_fraction_x_contacts":
        fields.pop("anchor_fraction")
    payload = json.dumps(
        {"diagram": diagram, "condition": fields}, sort_keys=True, separators=(",", ":")
    )
    offset = int.from_bytes(hashlib.sha256(payload.encode("utf-8")).digest()[:4], "big")
    return int(base_seed + 1009 * replicate + offset)


def long_metrics(
    condition_id: str, diagram: str, configuration: str, replicate: int,
    condition: SyntheticCondition, metrics: dict[str, object], truth_use: str,
    likelihood_match_status: str | None = None,
) -> list[dict[str, object]]:
    rows = []
    for metric, value in metrics.items():
        if isinstance(value, (int, float, np.integer, np.floating)):
            rows.append({
                "condition_id": condition_id,
                "diagram": diagram,
                "generator_family": condition.generator_family,
                "configuration": configuration,
                "replicate": replicate,
                "seed": condition.seed,
                "metric": metric,
                "value": value,
                "truth_used_in_inference": truth_use,
                "likelihood_match_status": likelihood_match_status or (
                    "MATCHED_PAIRWISE_ENDPOINT_LIKELIHOOD"
                    if condition.generator_family == "well_specified"
                    and condition.molecule_count == 0
                    and condition.count_overdispersion == 0
                    else "MISMATCHED_OR_MOLECULE_COMPOSITE"
                ),
            })
    return rows


def oracle_contact_metrics(dataset, probabilities: np.ndarray) -> dict[str, object]:
    """Evaluate an explicitly ORACLE state distribution on non-anchor contacts."""
    evaluation_signal = ~dataset.is_background & ~dataset.anchor_mask
    metrics, _, selective = contact_metrics(
        dataset.truth_state[evaluation_signal], probabilities[evaluation_signal]
    )
    cis = dataset.chrom_index[dataset.pair_i] == dataset.chrom_index[dataset.pair_j]
    scoped: dict[str, object] = {}
    for scope, selected in (
        ("cis", evaluation_signal & cis),
        ("trans", evaluation_signal & ~cis),
    ):
        if selected.any():
            values, _, _ = contact_metrics(dataset.truth_state[selected], probabilities[selected])
            for name in (
                "n_contacts", "exact_top1_accuracy", "exact_log_loss", "exact_brier",
                "same_cross_accuracy", "same_cross_log_loss", "same_cross_brier",
                "exact_ece", "same_cross_ece",
            ):
                scoped[f"{scope}_{name}"] = values[name]
    output: dict[str, object] = {
        **metrics,
        **{
            f"geometry_gauge_aligned_{name}": value
            for name, value in metrics.items()
        },
        **scoped,
        **{
            f"geometry_gauge_aligned_{name}": value
            for name, value in scoped.items()
        },
        "n_anchor_contacts": int(dataset.anchor_mask.sum()),
        "actual_anchor_fraction": float(dataset.anchor_mask.mean()),
        "runtime_seconds": 0.0,
        "converged": 1,
        "failure": 0,
        "unordered_structure_rmsd": 0.0,
        "unordered_structure_distance_spearman": 1.0,
        "globally_aligned_structure_rmsd": 0.0,
        "gauge_recovery": 1.0,
        "relative_gauge_recovery": 1.0,
        "mean_homolog_separation_spearman": 1.0,
    }
    for row in selective:
        metric_name = (
            f"{row['metric']}_coverage_at_precision_{float(row['target_precision']):.2f}"
        )
        output[metric_name] = float(row["callable_fraction"])
        output[f"geometry_gauge_aligned_{metric_name}"] = float(row["callable_fraction"])
    return output


def run_grid(quick: bool, base_seed: int) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    grid = condition_grid(quick, base_seed)
    replicates = 2 if quick else 3
    condition_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    wide_rows: list[dict[str, object]] = []
    for item in grid:
        diagram = str(item["diagram"])
        grid_stage = str(item["grid_stage"])
        base_condition = item["condition"]
        assert isinstance(base_condition, SyntheticCondition)
        condition_id = condition_identifier(diagram, base_condition)
        condition_fields = asdict(base_condition)
        condition_fields.pop("seed")
        replicate_seeds = [
            condition_replicate_seed(base_seed, replicate, diagram, base_condition)
            for replicate in range(replicates)
        ]
        condition_rows.append({
            "condition_id": condition_id,
            "diagram": diagram,
            **condition_fields,
            "base_seed": base_seed,
            "replicate_seeds": ",".join(map(str, replicate_seeds)),
            "replicates": replicates,
            "grid_policy": "scientific_slices_not_cartesian",
            "grid_stage": grid_stage,
            "truth_used_in_inference": (
                "ORACLE_EXTERNAL_PHASE_ANCHORS"
                if base_condition.anchor_fraction > 0 else "none"
            ),
            "likelihood_match_status": (
                "MATCHED_PAIRWISE_ENDPOINT_LIKELIHOOD"
                if base_condition.generator_family == "well_specified"
                and base_condition.molecule_count == 0
                and base_condition.count_overdispersion == 0
                else "MISMATCHED_OR_MOLECULE_COMPOSITE"
            ),
        })
        for replicate in range(replicates):
            seed = condition_replicate_seed(
                base_seed, replicate, diagram, base_condition
            )
            condition = replace(base_condition, seed=seed)
            dataset = generate_synthetic(condition)
            config = InferenceConfig(
                seed=seed, iterations=55 if quick else 130, learning_rate=0.025,
                likelihood_mode="monotone_log", staged_cis_then_trans=False,
                molecule_likelihood=False,
                observation_model=(
                    "normalized_endpoint"
                    if condition.molecule_count == 0
                    else "contrastive"
                ),
            )
            result = infer_synthetic(dataset, config)
            evaluated = evaluate_synthetic(dataset, result)
            uses_oracle_anchors = bool(dataset.anchor_mask.any())
            if uses_oracle_anchors:
                configuration = "ORACLE_external_phase_anchors_monotone_log_midpoint_delta_pairwise"
            elif config.observation_model == "normalized_endpoint":
                configuration = "blind_normalized_endpoint_monotone_log_midpoint_delta_pairwise"
            else:
                configuration = "blind_contrastive_monotone_log_midpoint_delta_pairwise"
            row = {
                "condition_id": condition_id, "diagram": diagram, "grid_stage": grid_stage,
                "generator_family": condition.generator_family,
                "configuration": configuration,
                "replicate": replicate, "seed": seed,
                **asdict(condition), **evaluated,
                "truth_used_in_inference": (
                    "ORACLE_EXTERNAL_PHASE_ANCHORS" if uses_oracle_anchors else "none"
                ),
                "likelihood_match_status": (
                    "MATCHED_PAIRWISE_ENDPOINT_LIKELIHOOD"
                    if condition.generator_family == "well_specified"
                    and condition.molecule_count == 0
                    and condition.count_overdispersion == 0
                    and config.observation_model == "normalized_endpoint"
                    else "MISMATCHED_OR_MOLECULE_COMPOSITE"
                ),
            }
            wide_rows.append(row)
            metric_rows.extend(long_metrics(
                condition_id, diagram, str(row["configuration"]), replicate,
                condition, evaluated, str(row["truth_used_in_inference"]),
            ))

            oracle_q = oracle_coordinate_posteriors(dataset, "monotone_log")
            oracle_metrics = oracle_contact_metrics(dataset, oracle_q)
            oracle_configuration = "ORACLE_coordinates_inference_monotone_log"
            metric_rows.extend(long_metrics(
                condition_id, diagram, oracle_configuration, replicate,
                condition, oracle_metrics, "ORACLE_COORDINATES",
                (
                    "MATCHED_PAIRWISE_ENDPOINT_LIKELIHOOD"
                    if condition.generator_family == "well_specified"
                    and condition.molecule_count == 0
                    and condition.count_overdispersion == 0
                    else "MISMATCHED_INFERENCE_SCORE_OR_MOLECULE_COMPOSITE"
                ),
            ))
            oracle_row = {
                "condition_id": condition_id, "diagram": diagram, "grid_stage": grid_stage,
                "generator_family": condition.generator_family,
                "configuration": oracle_configuration,
                "replicate": replicate, "seed": seed,
                **asdict(condition), **oracle_metrics,
                "truth_used_in_inference": "ORACLE_COORDINATES",
                "likelihood_match_status": (
                    "MATCHED_PAIRWISE_ENDPOINT_LIKELIHOOD"
                    if condition.generator_family == "well_specified"
                    and condition.molecule_count == 0
                    and condition.count_overdispersion == 0
                    else "MISMATCHED_INFERENCE_SCORE_OR_MOLECULE_COMPOSITE"
                ),
            }
            wide_rows.append(oracle_row)

            if condition.molecule_count == 0:
                generator_q = oracle_generator_conditional_posteriors(dataset)
                generator_metrics = oracle_contact_metrics(dataset, generator_q)
                generator_configuration = "ORACLE_coordinates_generator_conditional"
                metric_rows.extend(long_metrics(
                    condition_id, diagram, generator_configuration, replicate,
                    condition, generator_metrics, "ORACLE_COORDINATES_AND_GENERATOR",
                    "ORACLE_GENERATOR_CONDITIONAL",
                ))
                wide_rows.append({
                    "condition_id": condition_id, "diagram": diagram,
                    "grid_stage": grid_stage,
                    "generator_family": condition.generator_family,
                    "configuration": generator_configuration,
                    "replicate": replicate, "seed": seed,
                    **asdict(condition), **generator_metrics,
                    "truth_used_in_inference": "ORACLE_COORDINATES_AND_GENERATOR",
                    "likelihood_match_status": "ORACLE_GENERATOR_CONDITIONAL",
                })

            if diagram == "readchain_length_x_molecule_count" and condition.molecule_count > 0:
                molecule_result = infer_synthetic(dataset, replace(config, molecule_likelihood=True))
                molecule_evaluated = evaluate_synthetic(dataset, molecule_result)
                molecule_row = {
                    **row,
                    "configuration": "blind_monotone_log_molecule_normalized_composite",
                    **molecule_evaluated,
                }
                wide_rows.append(molecule_row)
                metric_rows.extend(long_metrics(
                    condition_id, diagram, str(molecule_row["configuration"]), replicate,
                    condition, molecule_evaluated, str(row["truth_used_in_inference"]),
                ))
    return condition_rows, metric_rows, wide_rows


def summarize(wide: list[dict[str, object]]) -> list[dict[str, object]]:
    metrics = (
        "exact_top1_accuracy", "same_cross_accuracy", "exact_log_loss", "same_cross_log_loss",
        "exact_ece", "same_cross_ece", "unordered_structure_rmsd",
        "geometry_gauge_aligned_exact_top1_accuracy",
        "geometry_gauge_aligned_same_cross_accuracy",
        "geometry_gauge_aligned_exact_log_loss",
        "geometry_gauge_aligned_same_cross_log_loss",
        "geometry_gauge_aligned_exact_ece", "geometry_gauge_aligned_same_cross_ece",
        "unordered_structure_distance_spearman", "globally_aligned_structure_rmsd",
        "gauge_recovery", "relative_gauge_recovery", "mean_homolog_separation_spearman",
        "cis_exact_top1_accuracy", "cis_same_cross_accuracy",
        "trans_exact_top1_accuracy", "trans_same_cross_accuracy",
        "geometry_gauge_aligned_cis_exact_top1_accuracy",
        "geometry_gauge_aligned_cis_same_cross_accuracy",
        "geometry_gauge_aligned_trans_exact_top1_accuracy",
        "geometry_gauge_aligned_trans_same_cross_accuracy",
        "runtime_seconds", "converged", "failure",
        "actual_anchor_fraction", "n_anchor_contacts",
        "background_accuracy", "background_log_loss", "background_brier", "background_ece",
        "exact_coverage_at_precision_0.55", "exact_coverage_at_precision_0.70",
        "exact_coverage_at_precision_0.90", "same_cross_coverage_at_precision_0.70",
        "same_cross_coverage_at_precision_0.90",
        "geometry_gauge_aligned_exact_coverage_at_precision_0.55",
        "geometry_gauge_aligned_exact_coverage_at_precision_0.70",
        "geometry_gauge_aligned_exact_coverage_at_precision_0.90",
        "geometry_gauge_aligned_same_cross_coverage_at_precision_0.70",
        "geometry_gauge_aligned_same_cross_coverage_at_precision_0.90",
        "exact_signal_aware_coverage_at_precision_0.55",
        "exact_signal_aware_coverage_at_precision_0.70",
        "exact_signal_aware_coverage_at_precision_0.90",
        "same_cross_signal_aware_coverage_at_precision_0.70",
        "same_cross_signal_aware_coverage_at_precision_0.90",
        "fitted_signal_fraction",
    )
    groups: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in wide:
        groups.setdefault((str(row["condition_id"]), str(row["configuration"])), []).append(row)
    output = []
    for (condition_id, configuration), rows in sorted(groups.items()):
        first = rows[0]
        result: dict[str, object] = {
            "condition_id": condition_id, "diagram": first["diagram"],
            "grid_stage": first.get("grid_stage", "coarse"),
            "generator_family": first["generator_family"], "configuration": configuration,
            "n_replicates": len(rows), "contacts_per_cell": first["contacts_per_cell"],
            "homolog_separation": first["homolog_separation"],
            "trans_fraction": 1.0 - float(first["cis_fraction"]),
            "background_fraction": first["background_fraction"], "anchor_fraction": first["anchor_fraction"],
            "readchain_length": first["readchain_length"], "molecule_count": first["molecule_count"],
            "truth_used_in_inference": first["truth_used_in_inference"],
            "likelihood_match_status": first["likelihood_match_status"],
        }
        for metric in metrics:
            values = np.asarray([float(row.get(metric, math.nan)) for row in rows], dtype=float)
            finite = values[np.isfinite(values)]
            result[f"{metric}_mean"] = float(finite.mean()) if finite.size else math.nan
            result[f"{metric}_sd"] = float(finite.std(ddof=1)) if finite.size > 1 else (0.0 if finite.size else math.nan)
            result[f"{metric}_n_finite_replicates"] = int(finite.size)
        output.append(result)
    return output


def phase_diagrams(output: Path, summary: list[dict[str, object]]) -> None:
    specs = {
        "contacts_x_homolog_separation": ("contacts_per_cell", "homolog_separation", "Contacts per cell", "Homolog separation"),
        "trans_fraction_x_background": ("trans_fraction", "background_fraction", "Trans fraction", "Background fraction"),
        "anchor_fraction_x_contacts": ("contacts_per_cell", "anchor_fraction", "Contacts per cell", "Anchor fraction"),
        "readchain_length_x_molecule_count": ("molecule_count", "readchain_length", "Molecule count", "Readchain length"),
    }
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for diagram, (x_name, y_name, x_label, y_label) in specs.items():
        figure, axes = plt.subplots(4, 2, figsize=(6.4, 11.6), constrained_layout=True)
        for column, family in enumerate(("well_specified", "model_mismatched")):
            rows = [
                row for row in summary
                if row["diagram"] == diagram and row["generator_family"] == family
                and str(row["configuration"]).endswith("monotone_log_midpoint_delta_pairwise")
            ]
            xs = sorted(set(float(row[x_name]) for row in rows))
            ys = sorted(set(float(row[y_name]) for row in rows))
            plot_metrics = (
                ("exact_top1_accuracy_mean", "raw-gauge exact state"),
                ("same_cross_accuracy_mean", "raw-gauge same/cross"),
                ("geometry_gauge_aligned_exact_top1_accuracy_mean", "geometry-gauge-aligned exact state"),
                ("geometry_gauge_aligned_same_cross_accuracy_mean", "geometry-gauge-aligned same/cross"),
            )
            for row_index, (metric, metric_label) in enumerate(plot_metrics):
                matrix = np.full((len(ys), len(xs)), np.nan)
                for row in rows:
                    matrix[ys.index(float(row[y_name])), xs.index(float(row[x_name]))] = float(row[metric])
                image = axes[row_index, column].imshow(
                    matrix, origin="lower", aspect="auto", vmin=0.0, vmax=1.0, cmap="viridis"
                )
                axes[row_index, column].set_xticks(range(len(xs)), [f"{value:g}" for value in xs])
                axes[row_index, column].set_yticks(range(len(ys)), [f"{value:g}" for value in ys])
                axes[row_index, column].set_xlabel(x_label)
                axes[row_index, column].set_ylabel(y_label)
                axes[row_index, column].set_title(
                    f"{family.replace('_', ' ')}: {metric_label}"
                )
                figure.colorbar(image, ax=axes[row_index, column], fraction=0.046, pad=0.04, label="Accuracy")
        if diagram == "anchor_fraction_x_contacts":
            figure.suptitle("Anchor fraction > 0: ORACLE_EXTERNAL_PHASE_ANCHORS")
        figure.savefig(figures / f"{diagram}.pdf", dpi=300)
        plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results/phase_diagnostics/059"))
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    require_analysis()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    condition_rows, metrics, wide = run_grid(args.quick, args.seed)
    replicate_summary = summarize(wide)
    write_tsv(args.output_dir / "conditions.tsv", condition_rows)
    write_tsv(args.output_dir / "metrics_long.tsv", metrics)
    write_tsv(args.output_dir / "replicate_summary.tsv", replicate_summary)
    phase_diagrams(args.output_dir, replicate_summary)
    config = {
        "logical_experiment_id": "059", "formal_sequence_alias": "assigned_at_packaging",
        "mode": "quick" if args.quick else "full", "seed": args.seed,
        "replicates": 2 if args.quick else 3,
        "grid_policy": "four low-dimensional scientific slices; no broad Cartesian sweep",
        "transition_refinement": "contacts x separation only: depths 1500/6000 and separations 1.25/1.6 after preliminary coarse transition",
        "anchor_fractions": [0, 0.001, 0.005, 0.01, 0.05],
        "anchor_truth_boundary": "anchor fraction > 0 is ORACLE_EXTERNAL_PHASE_ANCHORS",
        "anchor_mass_policy": "one joint endpoint-and-state likelihood term per anchored observation; anchors excluded from contact evaluation",
        "anchor_pairing_policy": "within generator/depth/replicate, anchor masks are nested prefixes on one identical geometry and contact dataset",
        "minimum_nonzero_anchor_count": 1,
        "regularization_depth_policy": "fixed prior calibrated at effective observation mass 1200; relative prior scale is 1200/effective_mass",
        "resolution_bp": 1_000_000,
    }
    (args.output_dir / "config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_provenance(
        args.output_dir / "provenance.json",
        collect_provenance(Path(__file__).resolve().parents[1], sys.argv, args.seed, config, ()),
    )


if __name__ == "__main__":
    main()
