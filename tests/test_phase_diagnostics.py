from __future__ import annotations

import itertools
import csv
import hashlib
import json
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from phase_diagnostics import provenance as provenance_module
from phase_diagnostics.likelihood import ScoreConfig, is_nonincreasing, log_score, posterior_from_distances
from phase_diagnostics.mass import assert_mass_conserved, mstep_contribution, normalized_pair_expansion
from phase_diagnostics.p9016 import read_3dg
from phase_diagnostics.metrics import (
    contact_metrics,
    selective_call_curve,
    selective_coverage,
    selective_threshold_curve,
    structure_metrics,
)
from phase_diagnostics.splits import Observation, blocked_split, leakage_audit
from phase_diagnostics.state import (
    CanonicalAccumulator,
    STATE_LABELS,
    apply_chromosome_gauge,
    read_p4,
    relabel_probabilities,
    swap_endpoint_state,
    write_p4,
)
from phase_diagnostics.synthetic import (
    InferenceConfig,
    InferenceResult,
    SyntheticCondition,
    evaluate_synthetic,
    generate_synthetic,
    geometry_selected_chromosome_flips,
    infer_synthetic,
    observation_masses,
    oracle_coordinate_posteriors,
    oracle_generator_conditional_posteriors,
    regularization_scale,
    relabel_contact_probabilities_to_geometry_gauge,
    score_coordinates,
    synchronize_result_z2,
)
from phase_diagnostics.z2 import RelativeFlipEdge, synchronize_z2
from scripts import package_phase_formal_results as packaging
from scripts import run_phase_model_comparison as model_comparison


def _write_test_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_complete_057b_fixture(
    root: Path,
    *,
    truth_boundary: str = "ORACLE_CONTACT_PHASE",
    declared_seeds: tuple[int, ...] = (17, 29),
    arm_seeds: dict[str, tuple[int, ...]] | None = None,
    chromosome_overrides: dict[tuple[str, int], int] | None = None,
) -> None:
    arm_seeds = arm_seeds or {
        arm: declared_seeds for arm in packaging.EXPECTED_057B_ARMS
    }
    chromosome_overrides = chromosome_overrides or {}
    _write_test_tsv(root / "057/run_manifest.tsv", [
        {"field": "057B_status", "value": "COMPLETE_THREE_EXPLICIT_ARMS"},
        {"field": "057B_seeds", "value": ",".join(map(str, declared_seeds))},
    ])
    run_rows: list[dict[str, object]] = []
    long_rows: list[dict[str, object]] = []
    chromosome_rows: list[dict[str, object]] = []
    for arm in packaging.EXPECTED_057B_ARMS:
        for seed in arm_seeds[arm]:
            offset = 0.0 if seed == 17 else 0.1
            n_chromosomes = chromosome_overrides.get((arm, seed), 1)
            numeric = {
                "n_chromosomes": n_chromosomes,
                "mean_copy_swap_invariant_procrustes_rmsd": 0.2 + offset,
                "mean_copy_swap_invariant_rigid_rmsd": 0.3 + offset,
                "mean_copy_swap_invariant_distance_pearson": 0.75 + offset,
                "mean_copy_swap_invariant_distance_spearman": 0.8 + offset,
                "mean_chromosome_shape_preservation": 0.77 + offset,
                "mean_homolog_separation_spearman": 0.7 + offset,
                "chromosome_centroid_distance_spearman": 0.4 + offset,
                "runtime_seconds": 2.0 + offset,
            }
            configuration = f"{arm}_seed{seed}"
            run_rows.append({
                "arm": arm,
                "seed": seed,
                **numeric,
                "truth_used_in_inference": truth_boundary,
            })
            for metric, value in numeric.items():
                long_rows.append({
                    "component": "057B_explicit_oracle_phase_mstep",
                    "configuration": configuration,
                    "metric": metric,
                    "value": value,
                    "truth_used_in_inference": truth_boundary,
                })
            for chromosome in range(n_chromosomes):
                chromosome_rows.append({
                    "component": "057B_oracle_phase_mstep",
                    "configuration": configuration,
                    "chromosome": f"chr{chromosome + 1}",
                    "procrustes_rmsd": numeric["mean_copy_swap_invariant_procrustes_rmsd"],
                    "rigid_rmsd": numeric["mean_copy_swap_invariant_rigid_rmsd"],
                    "distance_pearson": numeric["mean_copy_swap_invariant_distance_pearson"],
                    "distance_spearman": numeric["mean_copy_swap_invariant_distance_spearman"],
                    "chromosome_shape_preservation": numeric["mean_chromosome_shape_preservation"],
                    "homolog_separation_spearman": numeric["mean_homolog_separation_spearman"],
                    "truth_used_in_inference": truth_boundary,
                })
    _write_test_tsv(root / "057/oracle_mstep_runs.tsv", run_rows)
    _write_test_tsv(root / "057/metrics_long.tsv", long_rows)
    _write_test_tsv(root / "057/per_chromosome.tsv", chromosome_rows)


def test_packaging_git_output_preserves_porcelain_status_columns(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        packaging.subprocess,
        "check_output",
        lambda *args, **kwargs: " M results/phase_diagnostics/057/metrics_long.tsv\n",
    )
    assert packaging.git_output("status", "--porcelain").startswith(" M ")


def test_run_provenance_preserves_porcelain_status_columns(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        provenance_module.subprocess,
        "check_output",
        lambda *args, **kwargs: " M results/phase_diagnostics/057/metrics_long.tsv\n",
    )
    assert provenance_module._git(tmp_path, "status", "--porcelain").startswith(" M ")


def test_packaging_dirty_filter_distinguishes_outputs_from_source() -> None:
    output_lines = [
        " M results/phase_diagnostics/057/metrics_long.tsv",
        "?? test_res/119-example/",
        "R  results/old.tsv -> results/new.tsv",
    ]
    assert packaging.non_output_git_status(output_lines) == []
    source_line = " M scripts/package_phase_formal_results.py"
    assert packaging.non_output_git_status([source_line]) == [source_line]


def test_formal_057b_closure_rejects_undeclared_seed(tmp_path: Path) -> None:
    _write_test_tsv(tmp_path / "run_manifest.tsv", [
        {"field": "057B_status", "value": "COMPLETE_THREE_EXPLICIT_ARMS"},
        {"field": "057B_seeds", "value": "17"},
    ])
    run_rows: list[dict[str, object]] = []
    long_rows: list[dict[str, object]] = []
    chromosome_rows: list[dict[str, object]] = []
    numeric = {
        "n_chromosomes": 1,
        "mean_copy_swap_invariant_procrustes_rmsd": 0.2,
        "mean_copy_swap_invariant_rigid_rmsd": 0.3,
        "mean_copy_swap_invariant_distance_pearson": 0.75,
        "mean_copy_swap_invariant_distance_spearman": 0.8,
        "mean_chromosome_shape_preservation": 0.77,
        "mean_homolog_separation_spearman": 0.7,
        "chromosome_centroid_distance_spearman": 0.6,
        "runtime_seconds": 1.0,
    }
    for arm in packaging.EXPECTED_057B_ARMS:
        configuration = f"{arm}_seed17"
        run_rows.append({
            "arm": arm, "seed": 17, **numeric,
            "truth_used_in_inference": "ORACLE_CONTACT_PHASE",
        })
        for metric, value in numeric.items():
            long_rows.append({
                "component": "057B_explicit_oracle_phase_mstep",
                "configuration": configuration,
                "metric": metric,
                "value": value,
                "truth_used_in_inference": "ORACLE_CONTACT_PHASE",
            })
        chromosome_rows.append({
            "component": "057B_oracle_phase_mstep",
            "configuration": configuration,
            "chromosome": "chr1",
            "procrustes_rmsd": 0.2,
            "rigid_rmsd": 0.3,
            "distance_pearson": 0.75,
            "distance_spearman": 0.8,
            "chromosome_shape_preservation": 0.77,
            "homolog_separation_spearman": 0.7,
            "truth_used_in_inference": "ORACLE_CONTACT_PHASE",
        })
    _write_test_tsv(tmp_path / "oracle_mstep_runs.tsv", run_rows)
    _write_test_tsv(tmp_path / "metrics_long.tsv", long_rows)
    _write_test_tsv(tmp_path / "per_chromosome.tsv", chromosome_rows)
    assert packaging.assess_057b(tmp_path).complete

    run_rows.append({
        "arm": packaging.EXPECTED_057B_ARMS[0], "seed": 99, **numeric,
        "truth_used_in_inference": "ORACLE_CONTACT_PHASE",
    })
    _write_test_tsv(tmp_path / "oracle_mstep_runs.tsv", run_rows)
    status = packaging.assess_057b(tmp_path)
    assert not status.complete
    assert any("undeclared arm/seed" in problem for problem in status.problems)

    _write_test_tsv(tmp_path / "oracle_mstep_runs.tsv", run_rows[:-1])
    long_rows[1]["truth_used_in_inference"] = "none"
    _write_test_tsv(tmp_path / "metrics_long.tsv", long_rows)
    status = packaging.assess_057b(tmp_path)
    assert not status.complete
    assert any("non-ORACLE explicit 057B" in problem for problem in status.problems)

    long_rows[1]["truth_used_in_inference"] = "ORACLE_CONTACT_PHASE"
    chromosome_rows[0]["homolog_separation_spearman"] = 0.1
    _write_test_tsv(tmp_path / "metrics_long.tsv", long_rows)
    _write_test_tsv(tmp_path / "per_chromosome.tsv", chromosome_rows)
    status = packaging.assess_057b(tmp_path)
    assert not status.complete
    assert any("mean_homolog_separation_spearman" in problem for problem in status.problems)


def test_model_comparison_includes_current_oracle_mstep_arms(tmp_path: Path) -> None:
    _write_complete_057b_fixture(tmp_path)
    summaries = model_comparison.current_oracle_mstep_rows(tmp_path)
    assert len(summaries) == 3
    assert {row["truth_used_in_inference"] for row in summaries} == {
        "ORACLE_CONTACT_PHASE"
    }
    assert {row["n_replicates"] for row in summaries} == {2}
    fixed = next(
        row for row in summaries
        if row["configuration"].endswith("experimental_fixed_cis_rigid_trans")
    )
    assert fixed["structure_metric_name"] == "chromosome_centroid_distance_spearman"
    assert math.isclose(float(fixed["structure_metric"]), 0.45)
    assert math.isclose(float(fixed["runtime_seconds"]), 4.1)
    assert "placement-only runtime" in str(fixed["notes"])


def test_model_comparison_oracle_mstep_fails_closed(tmp_path: Path) -> None:
    with pytest.raises(OSError):
        model_comparison.current_oracle_mstep_rows(tmp_path)

    wrong_truth = tmp_path / "wrong_truth"
    _write_complete_057b_fixture(wrong_truth, truth_boundary="none_WRONG")
    with pytest.raises(ValueError, match="ORACLE"):
        model_comparison.current_oracle_mstep_rows(wrong_truth)

    wrong_seeds = tmp_path / "wrong_seeds"
    _write_complete_057b_fixture(wrong_seeds, arm_seeds={
        packaging.EXPECTED_057B_ARMS[0]: (17, 29),
        packaging.EXPECTED_057B_ARMS[1]: (17, 43),
        packaging.EXPECTED_057B_ARMS[2]: (17, 29),
    })
    with pytest.raises(ValueError, match="seed"):
        model_comparison.current_oracle_mstep_rows(wrong_seeds)

    wrong_chromosomes = tmp_path / "wrong_chromosomes"
    _write_complete_057b_fixture(
        wrong_chromosomes,
        chromosome_overrides={(packaging.EXPECTED_057B_ARMS[0], 29): 2},
    )
    with pytest.raises(ValueError, match="consistent chromosome count"):
        model_comparison.current_oracle_mstep_rows(wrong_chromosomes)

    single_seed = tmp_path / "single_seed"
    _write_complete_057b_fixture(single_seed, declared_seeds=(17,))
    with pytest.raises(ValueError, match="at least two declared seeds"):
        model_comparison.current_oracle_mstep_rows(single_seed)


def test_formal_059_standard_cross_file_identity(tmp_path: Path) -> None:
    plot_dir = tmp_path / "plots"
    plot_dir.mkdir()
    for name in packaging.STANDARD_PLOT_NAMES:
        (plot_dir / name).write_bytes(b"plot fixture")
    condition = {
        "resolution_bp": 1_000_000,
        "seed": 17,
        "anchor_fraction": 0.0,
        "chromosome_count": 2,
        "generator_family": "well_specified",
    }
    inference = {"seed": 17, "likelihood_mode": "monotone_log"}
    contract = {
        "condition": condition,
        "inference": inference,
        "truth_used_in_inference": "none",
    }
    serialized = json.dumps(contract, sort_keys=True, separators=(",", ":"))
    contract_hash = hashlib.sha256(serialized.encode("utf-8")).hexdigest()
    run_id = f"phase059_standard_{contract_hash[:16]}"
    (tmp_path / "standard_plot_config.json").write_text(json.dumps({
        **contract,
        "run_contract_sha256": contract_hash,
        "run_id": run_id,
    }), encoding="utf-8")
    common = {
        "run_id": run_id,
        "run_contract_sha256": contract_hash,
        "seed": 17,
        "generator_family": "well_specified",
    }
    summary = {
        **common,
        "truth_used_in_inference": "none",
        "copy_swap_policy": "ORACLE_GEOMETRY_ALIGNMENT_EVAL_ONLY",
        "resolution_bp": 1_000_000,
        "n_contacts": 100,
        "converged": 1,
    }
    _write_test_tsv(tmp_path / "summary.tsv", [summary])
    _write_test_tsv(tmp_path / "per_chromosome.tsv", [
        {"chromosome": "chr1", "geometry_selected_copy_flip": 0, **common},
        {"chromosome": "chr2", "geometry_selected_copy_flip": 1, **common},
    ])
    assert packaging.validate_059_standard(tmp_path, 17) == []

    summary["run_id"] = "mixed_run"
    _write_test_tsv(tmp_path / "summary.tsv", [summary])
    assert any("run_id" in problem for problem in packaging.validate_059_standard(tmp_path, 17))


def test_hickit_ab_coordinate_suffix_and_rigid_rmsd(tmp_path: Path) -> None:
    path = tmp_path / "hickit.3dg"
    path.write_text(
        "chr1a\t0\t0\t0\t0\nchr1b\t0\t0\t1\t0\n"
        "chr1a\t1000000\t1\t0\t0\nchr1b\t1000000\t1\t1\t0\n",
        encoding="utf-8",
    )
    coordinates = read_3dg(path, 1_000_000)
    assert set(coordinates) == {
        ("chr1", 0, 0), ("chr1", 0, 1),
        ("chr1", 1_000_000, 0), ("chr1", 1_000_000, 1),
    }
    reference = np.asarray(list(coordinates.values()))
    estimate = reference * 2.0
    metrics = structure_metrics(reference, estimate)
    assert metrics["procrustes_rmsd"] < 1e-12
    assert metrics["rigid_rmsd"] > 0.1


@pytest.mark.parametrize("state,expected", [(0, 0), (1, 2), (2, 1), (3, 3)])
def test_endpoint_swap_state_order(state: int, expected: int) -> None:
    assert swap_endpoint_state(state) == expected
    pure = np.eye(4)[state]
    assert np.argmax(relabel_probabilities(pure, swap_endpoints=True)) == expected


def test_full_canonicalization_serialization_roundtrip(tmp_path: Path) -> None:
    accumulator = CanonicalAccumulator()
    for state in range(4):
        pure = np.eye(4)[state]
        accumulator.add("chr1", 10, "chr2", 20, pure)
        accumulator.add("chr2", 20, "chr1", 10, relabel_probabilities(pure, swap_endpoints=True))
    combined = accumulator.normalized()
    assert list(combined) == [("chr1", 10, "chr2", 20)]
    np.testing.assert_allclose(combined[("chr1", 10, "chr2", 20)], np.full(4, 0.25))
    path = tmp_path / "p4.tsv"
    write_p4(path, combined)
    reloaded = read_p4(path)
    np.testing.assert_array_equal(reloaded[("chr1", 10, "chr2", 20)], combined[("chr1", 10, "chr2", 20)])


def test_p4_normalization_and_background_mixture() -> None:
    q4, q_background, _ = posterior_from_distances(
        [0.1, 0.5, 1.0, 2.0],
        ScoreConfig(mode="monotone_log", d0=0.1),
        background_prior=0.2,
        background_log_score=-1.0,
    )
    assert np.isclose(q4.sum() + q_background, 1.0)
    assert 0 < q_background < 1


def test_contact_metrics_condition_on_signal_mass() -> None:
    truth = np.asarray([0, 1, 2, 3])
    conditional = np.asarray([
        [0.7, 0.1, 0.1, 0.1],
        [0.1, 0.7, 0.1, 0.1],
        [0.1, 0.1, 0.7, 0.1],
        [0.1, 0.1, 0.1, 0.7],
    ])
    full, _, _ = contact_metrics(truth, conditional)
    subnormalized, _, _ = contact_metrics(truth, conditional * 0.2)
    assert subnormalized["exact_log_loss"] == pytest.approx(full["exact_log_loss"])
    assert subnormalized["same_cross_log_loss"] == pytest.approx(full["same_cross_log_loss"])
    assert subnormalized["exact_ece"] == pytest.approx(full["exact_ece"])
    assert subnormalized["mean_signal_probability"] == pytest.approx(0.2)


def test_uniform_prior_reports_random_tie_expectation() -> None:
    truth = np.asarray([0, 0, 0, 1, 2, 3])
    metrics, _, _ = contact_metrics(truth, np.full((truth.size, 4), 0.25))
    assert metrics["exact_top1_accuracy"] == pytest.approx(0.5)
    assert metrics["exact_random_tie_expected_accuracy"] == pytest.approx(0.25)
    assert metrics["same_cross_random_tie_expected_accuracy"] == pytest.approx(0.5)


def test_selective_coverage_never_splits_confidence_ties() -> None:
    correct = np.asarray([True, True, False, False])
    confidence = np.full(4, 0.25)
    curve = selective_call_curve(correct, confidence)
    assert curve == [{
        "callable_count": 4,
        "correct_count": 2,
        "coverage": 1.0,
        "precision": 0.5,
        "correct_recall": 1.0,
        "minimum_confidence": 0.25,
    }]
    rows = {row["target_precision"]: row for row in selective_coverage(correct, confidence)}
    assert rows[0.55]["callable_count"] == 0
    reversed_rows = {
        row["target_precision"]: row
        for row in selective_coverage(correct[::-1], confidence[::-1])
    }
    for target in rows:
        assert reversed_rows[target]["callable_count"] == rows[target]["callable_count"]
        assert reversed_rows[target]["callable_fraction"] == rows[target]["callable_fraction"]
    threshold_rows = selective_threshold_curve(correct, confidence, thresholds=(0.0, 0.25, 0.3))
    assert [row["callable_count"] for row in threshold_rows] == [4, 4, 0]
    assert threshold_rows[0]["precision"] == pytest.approx(0.5)


@pytest.mark.parametrize("posterior", [
    [0.25, 0.25, 0.25, 0.25],
    [0.5, 1 / 6, 1 / 6, 1 / 6],
    [0.75, 1 / 12, 1 / 12, 1 / 12],
    [0.95, 1 / 60, 1 / 60, 1 / 60],
    [0.999, 1 / 3000, 1 / 3000, 1 / 3000],
])
def test_source_observation_mass_conservation(posterior: list[float]) -> None:
    result = mstep_contribution(1.0, posterior, target_mode="posterior_once")
    assert_mass_conserved(result.state_mass, 1.0)


@pytest.mark.parametrize("segments", [2, 3, 4, 8])
def test_molecule_mass_conservation(segments: int) -> None:
    weights = normalized_pair_expansion(segments, source_mass=1.0)
    assert weights.size == segments * (segments - 1) // 2
    assert_mass_conserved(weights, 1.0)


def test_synthetic_molecule_observation_masses_sum_to_one() -> None:
    dataset = generate_synthetic(SyntheticCondition(
        contacts_per_cell=100, chromosome_count=2, chromosome_length_bins=6,
        readchain_length=4, molecule_count=5, seed=19,
    ))
    masses = observation_masses(dataset)
    for molecule_id in np.unique(dataset.molecule_id[dataset.molecule_id >= 0]):
        assert masses[dataset.molecule_id == molecule_id].sum() == pytest.approx(1.0)
    assert np.all(masses[dataset.molecule_id < 0] == 1.0)


def test_synthetic_anchor_mass_and_evaluation_exclusion() -> None:
    dataset = generate_synthetic(SyntheticCondition(
        contacts_per_cell=500, chromosome_count=2, chromosome_length_bins=6,
        background_fraction=0.0, anchor_fraction=0.001, seed=43,
    ))
    assert dataset.anchor_mask.sum() == 1
    masses = observation_masses(dataset)
    assert masses[dataset.anchor_mask].sum() == pytest.approx(1.0)
    perfect = np.eye(4)[dataset.truth_state]
    result = InferenceResult(
        coordinates=dataset.truth_coordinates.copy(),
        state_probabilities=perfect,
        background_probability=np.zeros(dataset.pair_i.size),
        objective_history=[0.0],
        runtime_seconds=0.0,
        converged=True,
        anchor_flips={},
    )
    metrics = evaluate_synthetic(dataset, result)
    assert metrics["n_anchor_contacts"] == 1
    assert metrics["n_contacts"] == 499
    assert metrics["exact_top1_accuracy"] == pytest.approx(1.0)


def test_synthetic_anchor_masks_are_nested_on_identical_data() -> None:
    common = dict(
        contacts_per_cell=1000, chromosome_count=2, chromosome_length_bins=6,
        background_fraction=0.05, seed=59,
    )
    low = generate_synthetic(SyntheticCondition(anchor_fraction=0.001, **common))
    high = generate_synthetic(SyntheticCondition(anchor_fraction=0.05, **common))
    np.testing.assert_array_equal(low.truth_coordinates, high.truth_coordinates)
    np.testing.assert_array_equal(low.pair_i, high.pair_i)
    np.testing.assert_array_equal(low.pair_j, high.pair_j)
    np.testing.assert_array_equal(low.truth_state, high.truth_state)
    assert np.all(~low.anchor_mask | high.anchor_mask)


def test_regularization_prior_weakens_with_observation_mass() -> None:
    assert regularization_scale(1200.0, 300.0) == pytest.approx(4.0)
    assert regularization_scale(1200.0, 1200.0) == pytest.approx(1.0)
    assert regularization_scale(1200.0, 12000.0) == pytest.approx(0.1)


def test_duplicate_empirical_distribution_doubles_mass_and_halves_prior_scale() -> None:
    base_mass = 400.0
    duplicated_mass = 2.0 * base_mass
    assert regularization_scale(1200.0, duplicated_mass) == pytest.approx(
        0.5 * regularization_scale(1200.0, base_mass)
    )


def test_geometry_gauge_aligned_contact_metrics_are_separate() -> None:
    dataset = generate_synthetic(SyntheticCondition(
        contacts_per_cell=600, chromosome_count=2, chromosome_length_bins=6,
        background_fraction=0.0, seed=47,
    ))
    flips = {0: 1, 1: 0}
    coordinates = dataset.truth_coordinates.copy()
    coordinates[dataset.chrom_index == 0] = coordinates[dataset.chrom_index == 0, ::-1]
    truth_probability = np.eye(4)[dataset.truth_state]
    arbitrary_gauge_probability = relabel_contact_probabilities_to_geometry_gauge(
        dataset, truth_probability, flips
    )
    result = InferenceResult(
        coordinates=coordinates,
        state_probabilities=arbitrary_gauge_probability,
        background_probability=np.zeros(dataset.pair_i.size),
        objective_history=[0.0],
        runtime_seconds=0.0,
        converged=True,
        anchor_flips={},
    )
    selected_flips = geometry_selected_chromosome_flips(dataset, coordinates)
    assert selected_flips == flips
    metrics = evaluate_synthetic(dataset, result)
    assert float(metrics["exact_top1_accuracy"]) < 1.0
    assert metrics["geometry_gauge_aligned_exact_top1_accuracy"] == pytest.approx(1.0)
    assert metrics["geometry_gauge_aligned_same_cross_accuracy"] == pytest.approx(1.0)


def test_mismatched_generator_oracle_uses_logistic_conditional() -> None:
    dataset = generate_synthetic(SyntheticCondition(
        generator_family="model_mismatched", contacts_per_cell=300,
        chromosome_count=2, chromosome_length_bins=6, molecule_count=0, seed=53,
    ))
    posterior = oracle_generator_conditional_posteriors(dataset)
    delta = (
        dataset.truth_coordinates[dataset.pair_i, :, None, :]
        - dataset.truth_coordinates[dataset.pair_j, None, :, :]
    )
    distance = np.linalg.norm(delta, axis=-1).reshape(-1, 4)
    expected = 1.0 / (1.0 + np.exp(np.clip((distance - 1.25) / 0.22, -60, 60)))
    expected += 0.015
    expected /= expected.sum(axis=1, keepdims=True)
    np.testing.assert_allclose(posterior, expected)
    assert not np.allclose(posterior, oracle_coordinate_posteriors(dataset))


def test_molecule_output_is_joint_assignment_pair_marginal() -> None:
    dataset = generate_synthetic(SyntheticCondition(
        contacts_per_cell=6, chromosome_count=2, chromosome_length_bins=6,
        readchain_length=4, molecule_count=1, background_fraction=0.0, seed=61,
    ))
    config = InferenceConfig(
        seed=61, iterations=2, observation_model="contrastive",
        molecule_likelihood=True, background_component=False,
    )
    result = infer_synthetic(dataset, config)
    molecule_id = 0
    segments = dataset.molecule_segments[molecule_id]
    assignments = np.asarray(list(itertools.product((0, 1), repeat=segments.size)))
    pair_positions = list(itertools.combinations(range(segments.size), 2))
    pair_scores = []
    for left_position, right_position in pair_positions:
        left = int(segments[left_position])
        right = int(segments[right_position])
        distances = np.linalg.norm(
            result.coordinates[left, :, None, :]
            - result.coordinates[right, None, :, :],
            axis=-1,
        )
        pair_scores.append(log_score(
            distances.ravel(),
            ScoreConfig(mode="monotone_log", alpha=config.alpha, d0=config.d0),
        ))
    pair_scores_array = np.asarray(pair_scores)
    assignment_scores = np.zeros(assignments.shape[0])
    for pair_index, (left_position, right_position) in enumerate(pair_positions):
        states = 2 * assignments[:, left_position] + assignments[:, right_position]
        assignment_scores += pair_scores_array[pair_index, states]
    assignment_scores /= len(pair_positions)
    assignment_weights = np.exp(assignment_scores - assignment_scores.max())
    assignment_weights /= assignment_weights.sum()
    segment_lookup = {int(bin_index): index for index, bin_index in enumerate(segments)}
    for observation_index in np.flatnonzero(dataset.molecule_id == molecule_id):
        left_position = segment_lookup[int(dataset.pair_i[observation_index])]
        right_position = segment_lookup[int(dataset.pair_j[observation_index])]
        states = 2 * assignments[:, left_position] + assignments[:, right_position]
        expected = np.asarray([
            assignment_weights[states == state].sum() for state in range(4)
        ])
        np.testing.assert_allclose(result.state_probabilities[observation_index], expected)


def _objective(
    coordinates: dict[tuple[str, int, int], np.ndarray],
    contacts: list[tuple[str, int, str, int, np.ndarray]],
) -> tuple[float, float]:
    estep = mstep = 0.0
    score_config = ScoreConfig(mode="monotone_log", d0=0.1)
    for chrom1, bin1, chrom2, bin2, q4 in contacts:
        distances = np.asarray([
            np.linalg.norm(coordinates[(chrom1, bin1, left)] - coordinates[(chrom2, bin2, right)])
            for left, right in itertools.product((0, 1), repeat=2)
        ])
        score = log_score(distances, score_config)
        estep += float(np.log(np.exp(score).sum() / 4.0))
        mstep += float(np.dot(q4, -score))
    return estep, mstep


@pytest.mark.parametrize("flipped", [
    {"chr1": 1},
    {"chr1": 1, "chr2": 1},
    {"chr1": 1, "chr2": 1, "chr3": 1},
])
def test_chromosome_copy_flip_symmetry(flipped: dict[str, int]) -> None:
    rng = np.random.default_rng(7)
    coordinates = {
        (chrom, bin_index, copy): rng.normal(size=3)
        for chrom in ("chr1", "chr2", "chr3")
        for bin_index in range(2)
        for copy in (0, 1)
    }
    contacts = [
        ("chr1", 0, "chr1", 1, np.asarray([0.1, 0.2, 0.3, 0.4])),
        ("chr1", 0, "chr2", 1, np.asarray([0.4, 0.1, 0.2, 0.3])),
        ("chr2", 0, "chr3", 1, np.asarray([0.2, 0.3, 0.1, 0.4])),
    ]
    before = _objective(coordinates, contacts)
    transformed = dict(coordinates)
    for chrom, should_flip in flipped.items():
        if should_flip:
            for bin_index in range(2):
                transformed[(chrom, bin_index, 0)], transformed[(chrom, bin_index, 1)] = (
                    coordinates[(chrom, bin_index, 1)], coordinates[(chrom, bin_index, 0)]
                )
    transformed_contacts = [
        (chrom1, bin1, chrom2, bin2, apply_chromosome_gauge(q4, chrom1, chrom2, flipped))
        for chrom1, bin1, chrom2, bin2, q4 in contacts
    ]
    after = _objective(transformed, transformed_contacts)
    np.testing.assert_allclose(before, after, atol=1e-12)
    before_parity = [q[0] + q[3] for *_, q in contacts]
    after_parity = [q[0] + q[3] for *_, q in transformed_contacts]
    # Same/cross itself is invariant only when both endpoints receive the same flip;
    # exact labels always relabel. The total objective remains invariant in every case.
    if all(flipped.get(c1, 0) == flipped.get(c2, 0) for c1, _, c2, _, _ in contacts):
        np.testing.assert_allclose(before_parity, after_parity)


def test_monotonicity_and_near_zero_stability() -> None:
    distance = np.concatenate(([0.0, 1e-15, 1e-12, 1e-9], np.linspace(1e-6, 5.0, 1000)))
    for mode in ("dist2", "logdist2", "monotonic_powerlaw", "monotone_log", "monotone_logistic"):
        config = ScoreConfig(mode=mode)
        values = log_score(distance, config)
        assert np.isfinite(values).all()
        assert is_nonincreasing(config, distance)
    assert not is_nonincreasing(ScoreConfig(mode="fdg_flat"), distance)


def test_blocked_split_prevents_molecule_and_pair_leakage() -> None:
    observations = [
        Observation("a", "chr1", 10, "chr2", 20, "m1"),
        Observation("b", "chr1", 12, "chr2", 22, "m1"),
        Observation("c", "chr1", 10, "chr2", 20, "m2"),
        Observation("d", "chr1", 2_000_010, "chr2", 3_000_020, "m3"),
        Observation("e", "chr2", 8_000_000, "chr3", 9_000_000, "m4"),
    ]
    for mode in ("raw_observation", "molecule_or_readchain", "exact_haploid_bin_pair", "genomic_block", "chromosome_pair"):
        assignment = blocked_split(observations, mode=mode, heldout_fraction=0.4, seed=11)
        assert assignment["a"] == assignment["b"]
        if mode in {"exact_haploid_bin_pair", "genomic_block", "chromosome_pair"}:
            assert assignment["a"] == assignment["c"]
    exact = blocked_split(observations, mode="exact_haploid_bin_pair", heldout_fraction=0.4, seed=11)
    audit = leakage_audit(observations, exact)
    assert audit["molecule_leakage_groups"] == 0
    assert audit["exact_haploid_bin_pair_leakage_groups"] == 0


def test_z2_connected_components_and_frustration() -> None:
    result = synchronize_z2([
        RelativeFlipEdge("chr1", "chr2", 0.9, 3.0),
        RelativeFlipEdge("chr2", "chr3", -0.8, 2.0),
        RelativeFlipEdge("chr1", "chr3", -0.7, 2.0),
    ], chromosomes=("chr1", "chr2", "chr3", "chr4"))
    assert result.spins["chr1"] == result.spins["chr2"]
    assert result.spins["chr2"] != result.spins["chr3"]
    assert "chr4" in result.unresolved
    assert result.frustrated_edge_fraction == 0


def test_z2_detects_frustrated_four_cycle_without_triangles() -> None:
    result = synchronize_z2([
        RelativeFlipEdge("chr1", "chr2", 0.9),
        RelativeFlipEdge("chr2", "chr3", 0.9),
        RelativeFlipEdge("chr3", "chr4", 0.9),
        RelativeFlipEdge("chr1", "chr4", -0.9),
    ])
    assert result.cycle_inconsistency == pytest.approx(1.0)


def test_blind_z2_does_not_read_background_truth() -> None:
    dataset = generate_synthetic(SyntheticCondition(
        contacts_per_cell=120, chromosome_count=3, chromosome_length_bins=5,
        background_fraction=0.4, seed=13,
    ))
    rng = np.random.default_rng(13)
    q4 = rng.dirichlet(np.ones(4), size=dataset.pair_i.size) * 0.8
    result = InferenceResult(
        coordinates=dataset.truth_coordinates.copy(), state_probabilities=q4,
        background_probability=np.full(dataset.pair_i.size, 0.2),
        objective_history=[1.0], runtime_seconds=0.0, converged=True,
        anchor_flips={},
    )
    left, left_diagnostics = synchronize_result_z2(dataset, result)
    changed_truth = replace(dataset, is_background=~dataset.is_background)
    right, right_diagnostics = synchronize_result_z2(changed_truth, result)
    np.testing.assert_allclose(left.coordinates, right.coordinates)
    np.testing.assert_allclose(left.state_probabilities, right.state_probabilities)
    assert left_diagnostics == right_diagnostics


def test_deterministic_synthetic_generation() -> None:
    condition = SyntheticCondition(contacts_per_cell=200, chromosome_count=2, chromosome_length_bins=6, seed=23)
    left = generate_synthetic(condition)
    right = generate_synthetic(condition)
    np.testing.assert_array_equal(left.truth_coordinates, right.truth_coordinates)
    np.testing.assert_array_equal(left.pair_i, right.pair_i)
    np.testing.assert_array_equal(left.truth_state, right.truth_state)


def test_high_depth_exposure_sampling_remains_defined() -> None:
    dataset = generate_synthetic(SyntheticCondition(
        generator_family="well_specified", contacts_per_cell=3000,
        chromosome_count=2, chromosome_length_bins=6, seed=31,
    ))
    result = infer_synthetic(dataset, InferenceConfig(seed=31, iterations=2))
    assert result.state_probabilities.shape == (dataset.pair_i.size, 4)
    assert np.isfinite(result.state_probabilities).all()


def test_well_specified_endpoint_likelihood_prefers_truth_geometry() -> None:
    dataset = generate_synthetic(SyntheticCondition(
        generator_family="well_specified", contacts_per_cell=3000,
        background_fraction=0.0, count_overdispersion=0.0,
        chromosome_count=2, chromosome_length_bins=6, seed=5,
    ))
    _, truth_log_probability = score_coordinates(
        dataset, dataset.truth_coordinates, observation_model="normalized_endpoint"
    )
    randomized = np.random.default_rng(5).normal(size=dataset.truth_coordinates.shape)
    _, random_log_probability = score_coordinates(
        dataset, randomized, observation_model="normalized_endpoint"
    )
    assert truth_log_probability.mean() > random_log_probability.mean()


def test_oracle_geometry_tiny_recovery_and_null_control() -> None:
    dataset = generate_synthetic(SyntheticCondition(
        generator_family="well_specified", homolog_separation=2.0, contacts_per_cell=2500,
        chromosome_count=2, chromosome_length_bins=8, background_fraction=0.0, seed=29,
    ))
    posterior = oracle_coordinate_posteriors(dataset)
    metrics, _, _ = contact_metrics(dataset.truth_state, posterior)
    assert metrics["same_cross_accuracy"] > 0.58
    rng = np.random.default_rng(29)
    null = rng.dirichlet(np.ones(4), size=dataset.truth_state.size)
    null_metrics, _, _ = contact_metrics(dataset.truth_state, null)
    assert null_metrics["exact_top1_accuracy"] < 0.30


@pytest.mark.slow
def test_blind_tiny_problem_is_finite_and_above_random_parity() -> None:
    dataset = generate_synthetic(SyntheticCondition(
        generator_family="well_specified", homolog_separation=3.0, contacts_per_cell=5000,
        cis_fraction=1.0, chromosome_count=1, chromosome_length_bins=6,
        background_fraction=0.0, anchor_fraction=0.0, seed=31,
    ))
    result = infer_synthetic(dataset, InferenceConfig(
        seed=31, iterations=180, learning_rate=0.025, background_component=False,
        observation_model="normalized_endpoint",
    ))
    metrics = evaluate_synthetic(dataset, result)
    assert dataset.anchor_mask.sum() == 0
    assert math.isfinite(float(metrics["unordered_structure_rmsd"]))
    assert float(metrics["same_cross_accuracy"]) > 0.70
