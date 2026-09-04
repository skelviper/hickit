#!/usr/bin/env python3
"""Build a scoped comparison matrix for current, corrected, prior, and ORACLE models."""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from phase_diagnostics.metrics import contact_metrics
from phase_diagnostics.provenance import collect_provenance, write_provenance
from phase_diagnostics.splits import Observation, blocked_split
from phase_diagnostics.synthetic import (
    InferenceConfig,
    SyntheticCondition,
    evaluate_synthetic,
    generate_synthetic,
    geometry_selected_chromosome_flips,
    infer_synthetic,
    oracle_coordinate_posteriors,
    relabel_contact_probabilities_to_geometry_gauge,
    score_coordinates,
    subset_synthetic,
    synchronize_result_z2,
)
from scripts.package_phase_formal_results import assess_057b


HISTORICAL_056 = Path(
    "/work/phase3/test_res/056-20260620_230016-p9016_readchain_dedup_training_1m"
)
HISTORICAL_078 = Path(
    "/work/phase3/hickit/test_res/078-20260712_211825-p9016_corrected_combined_candidate_1m"
)
HISTORICAL_083 = Path(
    "/work/phase3/hickit/test_res/083-20260713_111134-p9016_phased_hickit_oracle_gate_1m"
)


def require_analysis() -> None:
    if os.environ.get("CONDA_DEFAULT_ENV") != "analysis":
        raise RuntimeError("Model comparison must run inside conda environment 'analysis'")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty comparison: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def synthetic_observations(dataset) -> list[Observation]:
    return [
        Observation(
            observation_id=str(index),
            chrom1=f"chr{int(dataset.chrom_index[left]) + 1}",
            pos1=int(dataset.bin_offset[left]) * dataset.condition.resolution_bp,
            chrom2=f"chr{int(dataset.chrom_index[right]) + 1}",
            pos2=int(dataset.bin_offset[right]) * dataset.condition.resolution_bp,
            molecule_id=(f"mol:{int(dataset.molecule_id[index])}" if dataset.molecule_id[index] >= 0 else None),
        )
        for index, (left, right) in enumerate(zip(dataset.pair_i, dataset.pair_j, strict=True))
    ]


def metric_row(
    configuration: str, dataset, heldout, result, q4: np.ndarray, log_marginal: np.ndarray,
    *, notes: str, z2: dict[str, float | int | str] | None = None,
) -> dict[str, object]:
    signal = ~heldout.is_background
    cis = heldout.chrom_index[heldout.pair_i] == heldout.chrom_index[heldout.pair_j]
    raw_all_metrics, _, _ = contact_metrics(heldout.truth_state[signal], q4[signal])
    geometry_flips = geometry_selected_chromosome_flips(dataset, result.coordinates)
    aligned_q4 = relabel_contact_probabilities_to_geometry_gauge(
        heldout, q4, geometry_flips
    )
    all_metrics, _, selective = contact_metrics(
        heldout.truth_state[signal], aligned_q4[signal]
    )
    scoped = {}
    raw_scoped = {}
    for scope, selected in (("cis", signal & cis), ("trans", signal & ~cis)):
        raw_values, _, _ = contact_metrics(heldout.truth_state[selected], q4[selected])
        values, _, _ = contact_metrics(heldout.truth_state[selected], aligned_q4[selected])
        scoped[f"{scope}_exact"] = values["exact_top1_accuracy"]
        scoped[f"{scope}_same_cross"] = values["same_cross_accuracy"]
        raw_scoped[f"raw_gauge_{scope}_exact"] = raw_values["exact_top1_accuracy"]
        raw_scoped[f"raw_gauge_{scope}_same_cross"] = raw_values["same_cross_accuracy"]
    structure = evaluate_synthetic(dataset, result)
    coverage = {
        f"{entry['metric']}_coverage_at_precision_{float(entry['target_precision']):.2f}": entry["callable_fraction"]
        for entry in selective
    }
    return {
        "configuration": configuration,
        "evidence_scope": "synthetic_model_mismatched_blocked_heldout",
        "truth_used_in_inference": "none",
        **scoped,
        **raw_scoped,
        "exact_all": all_metrics["exact_top1_accuracy"],
        "same_cross_all": all_metrics["same_cross_accuracy"],
        "raw_gauge_exact_all": raw_all_metrics["exact_top1_accuracy"],
        "raw_gauge_same_cross_all": raw_all_metrics["same_cross_accuracy"],
        "exact_log_loss": all_metrics["exact_log_loss"],
        "same_cross_log_loss": all_metrics["same_cross_log_loss"],
        "calibration_ece": all_metrics["exact_ece"],
        "callable_fraction": coverage.get("exact_signal_aware_coverage_at_precision_0.70", math.nan),
        "callable_definition": "eval-only maximum signal-aware p(state) threshold coverage achieving exact precision >=0.70; confidence ties intact",
        "callable_fraction_precision_0.55": coverage.get("exact_signal_aware_coverage_at_precision_0.55", math.nan),
        "callable_fraction_precision_0.80": coverage.get("exact_signal_aware_coverage_at_precision_0.80", math.nan),
        "callable_fraction_precision_0.90": coverage.get("exact_signal_aware_coverage_at_precision_0.90", math.nan),
        "structure_metric": structure["unordered_structure_distance_spearman"],
        "structure_metric_name": "per_chromosome_unordered_distance_spearman",
        "heldout_metric": float(-log_marginal.mean()),
        "heldout_metric_name": "blind_normalized_endpoint_nll_or_labeled_contrastive_pseudo_nll_lower_better",
        "runtime_seconds": result.runtime_seconds,
        "n_evaluated_contacts": int(signal.sum()),
        "n_cells": 1,
        "n_chromosomes": dataset.condition.chromosome_count,
        "copy_swap_policy": "per-chromosome geometry-selected ORACLE alignment for evaluation metrics only",
        "notes": f"{notes}; raw arbitrary-gauge contact metrics retained in raw_gauge_* columns",
        **(z2 or {}),
    }


def run_synthetic(quick: bool, seed: int) -> list[dict[str, object]]:
    condition = SyntheticCondition(
        generator_family="model_mismatched", homolog_separation=1.25,
        contacts_per_cell=1400 if quick else 3500, cis_fraction=0.65,
        background_fraction=0.15, chromosome_count=3, chromosome_length_bins=10,
        readchain_length=4, molecule_count=45 if quick else 100,
        count_overdispersion=0.2, anchor_fraction=0.0, seed=seed,
    )
    dataset = generate_synthetic(condition)
    observations = synthetic_observations(dataset)
    assignment = blocked_split(
        observations, mode="exact_haploid_bin_pair", heldout_fraction=0.2,
        seed=seed, bin_size=condition.resolution_bp, block_size=5_000_000,
    )
    train_index = [index for index, row in enumerate(observations) if assignment[row.observation_id] == "train"]
    heldout_index = [index for index, row in enumerate(observations) if assignment[row.observation_id] == "heldout"]
    train = subset_synthetic(dataset, train_index)
    heldout = subset_synthetic(dataset, heldout_index)
    offsets = [0, 12] if quick else [0, 12, 26]
    seeds = [seed + offset for offset in offsets]
    specifications = [
        ("current_pipeline_FDG_flat_synthetic_proxy", InferenceConfig(likelihood_mode="fdg_flat", background_component=False, observation_model="normalized_endpoint"), "FDG shell score used as endpoint likelihood proxy; not a real-data rerun"),
        ("current_pipeline_DIST2_synthetic_proxy", InferenceConfig(likelihood_mode="dist2", background_component=False, observation_model="normalized_endpoint"), "DIST2 endpoint likelihood; same synthetic optimizer"),
        ("current_pipeline_LOGDIST2_synthetic_proxy", InferenceConfig(likelihood_mode="logdist2", background_component=False, observation_model="normalized_endpoint"), "LOGDIST2 endpoint likelihood; same synthetic optimizer"),
        ("monotone_contact_likelihood_posterior_once", InferenceConfig(likelihood_mode="monotone_log", background_component=False, observation_model="normalized_endpoint"), "normalized endpoint marginal applies latent posterior once; no duplicate alias row"),
        (
            "monotone_plus_background_legacy_regularizers_on",
            InferenceConfig(
                likelihood_mode="monotone_log", background_component=True,
                observation_model="normalized_endpoint",
                global_copy_track_weight=0.2, minimum_homolog_separation_weight=0.2,
                minimum_homolog_separation=1.5,
            ),
            "background contrastive component plus gauge-forcing copy-track and hard-separation prototypes",
        ),
        ("monotone_plus_background_conservative_regularizers_off", InferenceConfig(likelihood_mode="monotone_log", background_component=True, observation_model="normalized_endpoint"), "same endpoint/background mixture with global copy-track and hard minimum separation disabled; Delta may shrink"),
        ("staged_cis_then_trans", InferenceConfig(likelihood_mode="monotone_log", background_component=True, observation_model="normalized_endpoint", staged_cis_then_trans=True), "cis first, shape-preserving trans placement, weak joint refinement"),
        ("staged_plus_Z2_synchronization", InferenceConfig(likelihood_mode="monotone_log", background_component=True, observation_model="normalized_endpoint", staged_cis_then_trans=True), "no calibrated blind relative-gauge edges are available; all chromosomes remain unresolved and synchronization is a no-op"),
        ("molecule_level_readchain_normalized_composite", InferenceConfig(likelihood_mode="monotone_log", background_component=True, observation_model="contrastive", molecule_likelihood=True), "exact 2^m shared-copy enumeration under a unit-mass tempered composite; contrastive normalization, not an exact joint likelihood"),
    ]
    replicate_rows: list[dict[str, object]] = []
    for run_seed in seeds:
        cache: dict[InferenceConfig, object] = {}
        for name, template, note in specifications:
            config = replace(template, seed=run_seed, iterations=65 if quick else 150, learning_rate=0.025)
            result = cache.get(config)
            if result is None:
                result = infer_synthetic(train, config)
                cache[config] = result
            z2 = None
            if name == "staged_plus_Z2_synchronization":
                result, z2 = synchronize_result_z2(train, result)
            q4, log_marginal = score_coordinates(
                heldout, result.coordinates, likelihood_mode=config.likelihood_mode,
                alpha=config.alpha, d0=config.d0,
                observation_model=config.observation_model,
                signal_fraction=result.signal_fraction,
            )
            row = metric_row(name, train, heldout, result, q4, log_marginal, notes=note, z2=z2)
            row["seed"] = run_seed
            replicate_rows.append(row)
    aggregate: list[dict[str, object]] = []
    for name, _, note in specifications:
        rows = [row for row in replicate_rows if row["configuration"] == name]
        first = rows[0]
        result = {key: value for key, value in first.items() if key not in {"seed"}}
        result["n_replicates"] = len(rows)
        for key in (
            "cis_exact", "cis_same_cross", "trans_exact", "trans_same_cross", "exact_all", "same_cross_all",
            "raw_gauge_cis_exact", "raw_gauge_cis_same_cross",
            "raw_gauge_trans_exact", "raw_gauge_trans_same_cross",
            "raw_gauge_exact_all", "raw_gauge_same_cross_all",
            "exact_log_loss", "same_cross_log_loss", "calibration_ece", "callable_fraction",
            "callable_fraction_precision_0.55", "callable_fraction_precision_0.80",
            "callable_fraction_precision_0.90", "structure_metric", "heldout_metric", "runtime_seconds",
        ):
            values = np.asarray([float(row[key]) for row in rows], dtype=float)
            result[key] = float(np.mean(values))
            result[f"{key}_sd"] = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
        result["notes"] = (
            f"{note}; mean of {len(rows)} seeds; contact metrics use per-chromosome "
            "geometry-selected ORACLE alignment for evaluation only; raw arbitrary-gauge "
            "values are retained in raw_gauge_* columns"
        )
        aggregate.append(result)
    aggregate.extend(prior_and_oracle_synthetic(dataset, heldout_index))
    return aggregate


def prior_and_oracle_synthetic(dataset, heldout_index: list[int]) -> list[dict[str, object]]:
    heldout = subset_synthetic(dataset, heldout_index)
    signal = ~heldout.is_background
    cis = heldout.chrom_index[heldout.pair_i] == heldout.chrom_index[heldout.pair_j]
    result = []
    for name, q4, truth_use in (
        ("prior_only_uniform", np.full((heldout.pair_i.size, 4), 0.25), "none"),
        ("ORACLE_coordinates_monotone_log", oracle_coordinate_posteriors(heldout), "ORACLE_COORDINATES"),
    ):
        metrics, _, selective = contact_metrics(heldout.truth_state[signal], q4[signal])
        scoped = {}
        for scope, selected in (("cis", signal & cis), ("trans", signal & ~cis)):
            values, _, _ = contact_metrics(heldout.truth_state[selected], q4[selected])
            scoped[f"{scope}_exact"] = values["exact_top1_accuracy"]
            scoped[f"{scope}_same_cross"] = values["same_cross_accuracy"]
        coverage = next(
            float(row["callable_fraction"]) for row in selective
            if row["metric"] == "exact" and float(row["target_precision"]) == 0.70
        )
        result.append({
            "configuration": name, "evidence_scope": "synthetic_model_mismatched_blocked_heldout",
            "truth_used_in_inference": truth_use, **scoped,
            "exact_all": metrics["exact_top1_accuracy"], "same_cross_all": metrics["same_cross_accuracy"],
            "exact_log_loss": metrics["exact_log_loss"], "same_cross_log_loss": metrics["same_cross_log_loss"],
            "calibration_ece": metrics["exact_ece"], "callable_fraction": coverage,
            "callable_definition": "eval-only maximum threshold coverage achieving exact precision >=0.70; confidence ties intact",
            "structure_metric": "NA", "structure_metric_name": "NA",
            "heldout_metric": "NA", "heldout_metric_name": "NA", "runtime_seconds": 0.0,
            "n_evaluated_contacts": int(signal.sum()), "n_cells": 1,
            "n_chromosomes": dataset.condition.chromosome_count,
            "copy_swap_policy": "truth coordinate labels only for ORACLE row",
            "notes": "geometry-free prior" if truth_use == "none" else "ORACLE positive control",
            "n_replicates": 1,
        })
    return result


def current_oracle_mstep_rows(results_root: Path) -> list[dict[str, object]]:
    """Summarize the current 057B ORACLE M-step arms without scoring phase calls."""
    status = assess_057b(results_root / "057")
    if not status.complete:
        details = "; ".join(status.problems) or "057B closure is incomplete"
        raise ValueError(f"cannot summarize current 057B results: {details}")
    table = list(status.rows)
    expected = {
        "native_hickit_cis_only": (
            "ORACLE_contact_phase_native_hickit_cis_only",
            "mean_copy_swap_invariant_distance_spearman",
            "mean_per_chromosome_copy_swap_invariant_distance_spearman",
            "native Hickit cis-only reconstruction; chromosome placement is not identified",
        ),
        "native_hickit_cis_plus_trans": (
            "ORACLE_contact_phase_native_hickit_cis_plus_trans",
            "mean_copy_swap_invariant_distance_spearman",
            "mean_per_chromosome_copy_swap_invariant_distance_spearman",
            "native Hickit joint cis+trans reconstruction",
        ),
        "experimental_fixed_cis_rigid_trans_placement": (
            "ORACLE_contact_phase_experimental_fixed_cis_rigid_trans",
            "chromosome_centroid_distance_spearman",
            "chromosome_centroid_distance_spearman",
            "experimental rigid-body trans placement after native cis-only shapes were frozen",
        ),
    }
    observed_arms = {row.get("arm", "") for row in table}
    if observed_arms != set(expected):
        raise ValueError(
            "057B comparison requires exactly the three declared arms; "
            f"observed {sorted(observed_arms)}"
        )
    expected_seeds = set(status.seeds)
    if len(expected_seeds) < 2:
        raise ValueError(
            "057B comparison requires at least two declared seeds to quantify variability"
        )
    observed_chromosome_counts = {int(row["n_chromosomes"]) for row in table}
    if len(observed_chromosome_counts) != 1:
        raise ValueError(
            "057B comparison requires one consistent chromosome count; "
            f"observed {sorted(observed_chromosome_counts)}"
        )
    cis_runtime_by_seed = {
        int(row["seed"]): float(row["runtime_seconds"])
        for row in table if row["arm"] == "native_hickit_cis_only"
    }
    rows: list[dict[str, object]] = []
    for arm, (configuration, metric_column, metric_name, note) in expected.items():
        selected = [row for row in table if row["arm"] == arm]
        seeds = [int(row["seed"]) for row in selected]
        if set(seeds) != expected_seeds or len(seeds) != len(expected_seeds):
            raise ValueError(
                f"057B arm {arm!r} seed set {sorted(seeds)} does not match "
                f"declared seeds {sorted(expected_seeds)}"
            )
        metric_values = np.asarray([float(row[metric_column]) for row in selected])
        runtime_values = np.asarray([float(row["runtime_seconds"]) for row in selected])
        if not np.all(np.isfinite(metric_values)) or not np.all(np.isfinite(runtime_values)):
            raise ValueError(f"057B arm {arm!r} contains non-finite comparison values")
        extra = ""
        if arm == "experimental_fixed_cis_rigid_trans_placement":
            shape_values = np.asarray([
                float(row["mean_copy_swap_invariant_distance_spearman"])
                for row in selected
            ])
            placement_runtime = runtime_values.copy()
            runtime_values = np.asarray([
                cis_runtime_by_seed[seed] + placement
                for seed, placement in zip(seeds, placement_runtime, strict=True)
            ])
            extra = (
                f"; frozen-shape mean distance Spearman={shape_values.mean():.6f}"
                f"; mean placement-only runtime={placement_runtime.mean():.6f}s"
                "; reported runtime includes its required cis-only reconstruction"
            )
        rows.append({
            "configuration": configuration,
            "evidence_scope": "P9016_1Mb_ORACLE_phase_current_057B",
            "truth_used_in_inference": "ORACLE_CONTACT_PHASE",
            "cis_exact": "NA_NOT_AN_INFERENCE_TARGET",
            "cis_same_cross": "NA_NOT_AN_INFERENCE_TARGET",
            "trans_exact": "NA_NOT_AN_INFERENCE_TARGET",
            "trans_same_cross": "NA_NOT_AN_INFERENCE_TARGET",
            "exact_all": "NA_NOT_AN_INFERENCE_TARGET",
            "same_cross_all": "NA_NOT_AN_INFERENCE_TARGET",
            "calibration_ece": "NA_NOT_AN_INFERENCE_TARGET",
            "callable_fraction": "NA_NOT_AN_INFERENCE_TARGET",
            "callable_definition": "NA_ORACLE_PHASE_FIXED_BEFORE_MSTEP",
            "structure_metric": float(metric_values.mean()),
            "structure_metric_sd": float(metric_values.std(ddof=1)),
            "structure_metric_name": metric_name,
            "heldout_metric": "NA_NOT_HELD_OUT",
            "heldout_metric_name": "NA_NOT_HELD_OUT",
            "runtime_seconds": float(runtime_values.mean()),
            "runtime_seconds_sd": float(runtime_values.std(ddof=1)),
            "n_evaluated_contacts": "NA_STRUCTURE_ONLY",
            "n_cells": 1,
            "n_chromosomes": int(selected[0]["n_chromosomes"]),
            "n_replicates": len(selected),
            "copy_swap_policy": "per-chromosome unordered homolog assignment; geometry-selected for evaluation only",
            "notes": f"{note}; seeds={','.join(map(str, sorted(seeds)))}{extra}",
        })
    return rows


def historical_rows(results_root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    contact_path = HISTORICAL_056 / "eval/approved_baseline_pairs_only/contact_accuracy.tsv"
    if contact_path.exists():
        with contact_path.open(encoding="utf-8") as handle:
            table = list(csv.DictReader(handle, delimiter="\t"))
        selected = {
            row["scope"]: row for row in table
            if row["source"] == "reconstruction_posterior"
            and row["copy_swap_policy"] == "whole_chrom_snp_cis_top1_oracle_eval_only"
        }
        rows.append({
            "configuration": "historical_056_approved_pairs_only",
            "evidence_scope": "P9016_1Mb_historical_not_current_correctness_valid",
            "truth_used_in_inference": "none; SNP cis used only for evaluation copy alignment",
            "cis_exact": selected["genome_cis"]["top1_accuracy"],
            "cis_same_cross": "NA_NOT_COMPARABLE_LEGACY_ARGMAX_STATE_PARITY",
            "trans_exact": selected["genome_trans"]["top1_accuracy"],
            "trans_same_cross": "NA_NOT_COMPARABLE_LEGACY_ARGMAX_STATE_PARITY",
            "exact_all": selected["genome_all"]["top1_accuracy"],
            "same_cross_all": "NA_NOT_COMPARABLE_LEGACY_ARGMAX_STATE_PARITY",
            "calibration_ece": "NA", "callable_fraction": selected["genome_all"]["called_contact_fraction"],
            "callable_definition": "historical arm-specific called-contact fraction; NOT_COMPARABLE to target-precision coverage",
            "structure_metric": 0.8009042076870312, "structure_metric_name": "mean_per_chrom_cis_distance_spearman",
            "heldout_metric": "NA", "runtime_seconds": "NA", "n_evaluated_contacts": selected["genome_all"]["n_eval_contacts"],
            "n_cells": 1, "n_chromosomes": 20,
            "copy_swap_policy": "ORACLE SNP cis top1 for exact-state evaluation",
            "notes": "legacy force implementation later invalidated; legacy same/cross used argmax-state parity",
        })

    phase_path = HISTORICAL_078 / "phase_metric_summary.tsv"
    summary_path = HISTORICAL_078 / "summary.tsv"
    if phase_path.exists() and summary_path.exists():
        with phase_path.open(encoding="utf-8") as handle:
            phase = [row for row in csv.DictReader(handle, delimiter="\t") if row["arm"] == "K1_corrected_combined_expected"]
        scopes = {row["scope"]: row for row in phase}
        with summary_path.open(encoding="utf-8") as handle:
            summary = next(row for row in csv.DictReader(handle, delimiter="\t") if row["arm"] == "K1_corrected_combined_expected")
        callable_accuracy = float(scopes["all"]["callable_accuracy_pmax_ge_0p9"])
        callable_recall = float(scopes["all"]["callable_recall_pmax_ge_0p9"])
        callable_coverage = callable_recall / callable_accuracy if callable_accuracy > 0 else 0.0
        rows.append({
            "configuration": "current_corrected_078_K1",
            "evidence_scope": "P9016_1Mb_registered_train_disjoint_reused_final_test",
            "truth_used_in_inference": "none; phase opened only after blind scorer closure",
            "cis_exact": scopes["cis"]["top1_accuracy"], "cis_same_cross": scopes["cis"]["same_cross_accuracy"],
            "trans_exact": scopes["trans"]["top1_accuracy"], "trans_same_cross": scopes["trans"]["same_cross_accuracy"],
            "exact_all": scopes["all"]["top1_accuracy"], "same_cross_all": scopes["all"]["same_cross_accuracy"],
            "exact_log_loss": scopes["all"]["four_state_nll"], "calibration_ece": scopes["all"]["ece"],
            "callable_fraction": callable_coverage,
            "callable_definition": "fixed posterior pmax>=0.9 coverage (recall/accuracy)",
            "structure_metric": summary["final_selected_spearman_equal_chrom_mean"],
            "structure_metric_name": "mean_per_chrom_cis_distance_spearman",
            "heldout_metric": "NA_unnormalized_score_only", "runtime_seconds": "NA",
            "n_evaluated_contacts": scopes["all"]["n_contacts"], "n_cells": 1, "n_chromosomes": 20,
            "copy_swap_policy": "phase evaluation uses postclosure SNP cis swap; structure uses geometry swap",
            "notes": "best current correctness-audited real-data evidence; callable fraction is pmax>=0.9 coverage; final-test is reused, not model selection",
        })

    oracle_path = HISTORICAL_083 / "strict_eval/summary.tsv"
    if oracle_path.exists():
        with oracle_path.open(encoding="utf-8") as handle:
            oracle = [row for row in csv.DictReader(handle, delimiter="\t") if row["arm"].startswith("true_seed")]
        rows.append({
            "configuration": "ORACLE_contact_phase_historical_Hickit_joint_partial",
            "evidence_scope": "P9016_1Mb_ORACLE_phase_historical_083",
            "truth_used_in_inference": "ORACLE_CONTACT_PHASE",
            "cis_exact": "NA", "cis_same_cross": "NA", "trans_exact": "NA", "trans_same_cross": "NA",
            "exact_all": "NA", "same_cross_all": "NA", "calibration_ece": "NA", "callable_fraction": "NA",
            "callable_definition": "NA",
            "structure_metric": float(np.mean([float(row["strict_unordered_pair_spearman_mean"]) for row in oracle])),
            "structure_metric_name": "strict_unordered_pair_spearman_mean",
            "heldout_metric": "NA", "runtime_seconds": "NA", "n_evaluated_contacts": "NA",
            "n_cells": 1, "n_chromosomes": 20,
            "copy_swap_policy": "per-chromosome unordered homolog assignment",
            "notes": "ORACLE historical full cis+trans phase-known run retained for lineage; current explicit 057B arms are separate rows",
        })

    exp057 = results_root / "057/metrics_long.tsv"
    if exp057.exists():
        with exp057.open(encoding="utf-8") as handle:
            table = list(csv.DictReader(handle, delimiter="\t"))
        configurations = sorted({
            row["configuration"] for row in table
            if row["component"] in {"057A_oracle_coordinates_score_kernel_isolation", "057D_prior_only"}
        })
        for configuration in configurations:
            subset = [row for row in table if row["configuration"] == configuration]
            lookup = {(row["scope"], row["metric"]): row["value"] for row in subset}
            evidence_mode = "full" if int(subset[0]["n"]) > 50_000 else "quick"
            rows.append({
                "configuration": configuration,
                "evidence_scope": (
                    f"P9016_1Mb_ORACLE_coordinate_E_step_{evidence_mode}"
                    if "oracle_coords" in configuration
                    else f"P9016_1Mb_prior_only_{evidence_mode}"
                ),
                "truth_used_in_inference": subset[0]["truth_used_in_inference"],
                "cis_exact": lookup.get(("cis", "exact_top1_accuracy"), "NA"),
                "cis_same_cross": lookup.get(("cis", "same_cross_accuracy"), "NA"),
                "trans_exact": lookup.get(("trans", "exact_top1_accuracy"), "NA"),
                "trans_same_cross": lookup.get(("trans", "same_cross_accuracy"), "NA"),
                "exact_all": lookup.get(("all", "exact_top1_accuracy"), "NA"),
                "same_cross_all": lookup.get(("all", "same_cross_accuracy"), "NA"),
                "exact_log_loss": lookup.get(("all", "exact_log_loss"), "NA"),
                "same_cross_log_loss": lookup.get(("all", "same_cross_log_loss"), "NA"),
                "calibration_ece": lookup.get(("all", "exact_ece"), "NA"),
                "callable_fraction": lookup.get(("all", "exact_coverage_at_precision_0.70"), "NA"),
                "callable_definition": "eval-only maximum threshold coverage achieving exact precision >=0.70; confidence ties intact",
                "structure_metric": "NA", "heldout_metric": "NA", "runtime_seconds": 0,
                "n_evaluated_contacts": subset[0]["n"], "n_cells": 1, "n_chromosomes": 20,
                "copy_swap_policy": "CHARM/3DG reference copy labels for ORACLE coordinate rows",
                "notes": "fixed-coordinate E-step or geometry-free prior; not a full reconstruction",
            })
    rows.extend(current_oracle_mstep_rows(results_root))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, default=Path("results/phase_diagnostics"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/phase_diagnostics/model_comparison"))
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    require_analysis()
    rows = historical_rows(args.results_root) + run_synthetic(args.quick, args.seed)
    write_tsv(args.output_dir / "comparison_matrix.tsv", rows)
    write_provenance(
        args.output_dir / "provenance.json",
        collect_provenance(
            Path(__file__).resolve().parents[1], sys.argv, args.seed,
            {"mode": "quick" if args.quick else "full", "comparison_is_scoped": True}, (),
        ),
    )


if __name__ == "__main__":
    main()
