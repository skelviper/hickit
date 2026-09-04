#!/usr/bin/env python3
"""Generate the standard plot set for one representative synthetic condition."""

from __future__ import annotations

import argparse
import csv
import hashlib
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
from scipy.spatial import ConvexHull, QhullError
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from phase_diagnostics.metrics import chromosome_swap_invariant_metrics, contact_metrics
from phase_diagnostics.state import relabel_probabilities
from phase_diagnostics.synthetic import (
    InferenceConfig,
    SyntheticCondition,
    generate_synthetic,
    infer_synthetic,
    oracle_coordinate_posteriors,
)


COPY_COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00")


def require_analysis() -> None:
    if os.environ.get("CONDA_DEFAULT_ENV") != "analysis":
        raise RuntimeError("standard plots must run inside conda environment 'analysis'")


def align_to_reference(reference: np.ndarray, estimate: np.ndarray) -> np.ndarray:
    ref = reference.reshape(-1, 3)
    est = estimate.reshape(-1, 3)
    ref_center = ref.mean(axis=0, keepdims=True)
    est_center = est.mean(axis=0, keepdims=True)
    ref_zero = ref - ref_center
    est_zero = est - est_center
    left, singular, right_t = np.linalg.svd(est_zero.T @ ref_zero)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    denominator = float(np.square(est_zero).sum())
    scale = float(singular.sum() / denominator) if denominator > 1e-15 else 1.0
    return (est_zero @ rotation * scale + ref_center).reshape(reference.shape)


def volume(points: np.ndarray) -> float:
    try:
        return float(ConvexHull(points).volume)
    except QhullError:
        return math.nan


def relabel_to_reference(dataset, probabilities: np.ndarray, flips: dict[int, int]) -> np.ndarray:
    return np.vstack([
        relabel_probabilities(
            row,
            flip_left=flips[int(dataset.chrom_index[left])],
            flip_right=flips[int(dataset.chrom_index[right])],
        )
        for row, left, right in zip(
            probabilities, dataset.pair_i, dataset.pair_j, strict=True
        )
    ])


def common_limits(reference: np.ndarray, estimate: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    joined = np.concatenate((reference.reshape(-1, 3), estimate.reshape(-1, 3)))
    lower = joined.min(axis=0)
    upper = joined.max(axis=0)
    center = 0.5 * (lower + upper)
    radius = max(1e-6, 0.55 * float(np.max(upper - lower)))
    return center - radius, center + radius


def style_3d(axis, lower: np.ndarray, upper: np.ndarray, title: str) -> None:
    axis.set_xlim(lower[0], upper[0])
    axis.set_ylim(lower[1], upper[1])
    axis.set_zlim(lower[2], upper[2])
    axis.set_title(title)
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_zlabel("z")
    axis.set_box_aspect((1, 1, 1))


def plot_all_chromosomes(path: Path, dataset, reference: np.ndarray, estimate: np.ndarray) -> None:
    lower, upper = common_limits(reference, estimate)
    figure = plt.figure(figsize=(7.2, 3.0), constrained_layout=True)
    axes = [figure.add_subplot(1, 2, index + 1, projection="3d") for index in range(2)]
    for axis, coordinates, title in zip(axes, (reference, estimate), ("Synthetic truth", "Reconstruction"), strict=True):
        for chrom in range(dataset.condition.chromosome_count):
            selected = dataset.chrom_index == chrom
            for copy in (0, 1):
                color = COPY_COLORS[2 * chrom + copy]
                axis.scatter(
                    *coordinates[selected, copy].T, s=8, alpha=0.85, color=color,
                    label=f"chr{chrom + 1} copy{copy}", rasterized=True,
                )
        style_3d(axis, lower, upper, title)
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    figure.savefig(path, dpi=300)
    plt.close(figure)


def plot_chr1(path: Path, dataset, reference: np.ndarray, estimate: np.ndarray) -> None:
    lower, upper = common_limits(reference, estimate)
    selected = dataset.chrom_index == 0
    background = ~selected
    figure = plt.figure(figsize=(6.0, 3.0), constrained_layout=True)
    axes = [figure.add_subplot(1, 2, index + 1, projection="3d") for index in range(2)]
    for axis, coordinates, title in zip(axes, (reference, estimate), ("Synthetic truth", "Reconstruction"), strict=True):
        for copy in (0, 1):
            axis.scatter(
                *coordinates[background, copy].T, s=5, alpha=0.08, color="#777777",
                rasterized=True,
            )
            axis.scatter(
                *coordinates[selected, copy].T, s=13, alpha=0.95,
                color=COPY_COLORS[copy], label=f"chr1 copy{copy}", rasterized=True,
            )
        style_3d(axis, lower, upper, title)
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    figure.savefig(path, dpi=300)
    plt.close(figure)


def plot_chr1_distance_maps(path: Path, dataset, reference: np.ndarray, estimate: np.ndarray) -> None:
    selected = dataset.chrom_index == 0
    matrices = [
        squareform(pdist(coordinates[selected, copy]))
        for coordinates in (reference, estimate) for copy in (0, 1)
    ]
    vmax = max(float(matrix.max()) for matrix in matrices)
    figure, axes = plt.subplots(2, 2, figsize=(6.5, 6.0), constrained_layout=True)
    image = None
    for axis, matrix, title in zip(
        axes.ravel(), matrices,
        ("Truth copy0", "Truth copy1", "Reconstruction copy0", "Reconstruction copy1"),
        strict=True,
    ):
        image = axis.imshow(matrix, cmap="coolwarm_r", vmin=0, vmax=vmax, origin="lower")
        axis.set_title(title)
        axis.set_xlabel("chr1 bin")
        axis.set_ylabel("chr1 bin")
    assert image is not None
    figure.colorbar(image, ax=axes, fraction=0.025, pad=0.02, label="Euclidean distance (larger is blue)")
    figure.savefig(path, dpi=300)
    plt.close(figure)


def write_tables(
    output: Path, dataset, result, estimate: np.ndarray, flips: dict[int, int],
    *, run_id: str, run_contract_sha256: str,
) -> None:
    aligned_probability = relabel_to_reference(dataset, result.state_probabilities, flips)
    signal = ~dataset.is_background
    conditional = aligned_probability[signal]
    conditional /= conditional.sum(axis=1, keepdims=True)
    metrics, _, _ = contact_metrics(dataset.truth_state[signal], conditional)
    prediction = np.argmax(conditional, axis=1)
    confidence = np.max(conditional, axis=1)
    calls = confidence >= 0.9
    correct = prediction == dataset.truth_state[signal]
    oracle_probability = oracle_coordinate_posteriors(dataset)[signal]
    flat_left = conditional.ravel()
    flat_right = oracle_probability.ravel()
    summary = {
        **metrics,
        "pmax_ge_0p9_call_count": int(calls.sum()),
        "pmax_ge_0p9_coverage": float(calls.mean()),
        "pmax_ge_0p9_accuracy": float(correct[calls].mean()) if calls.any() else math.nan,
        "pmax_ge_0p9_recall": float((correct & calls).mean()),
        "reference_coordinate_p4_pearson": float(pearsonr(flat_left, flat_right).statistic),
        "reference_coordinate_p4_spearman": float(spearmanr(flat_left, flat_right).statistic),
        "copy_swap_policy": "ORACLE_GEOMETRY_ALIGNMENT_EVAL_ONLY",
        "truth_used_in_inference": "none",
        "resolution_bp": dataset.condition.resolution_bp,
        "runtime_seconds": result.runtime_seconds,
        "converged": int(result.converged),
        "run_id": run_id,
        "run_contract_sha256": run_contract_sha256,
        "seed": dataset.condition.seed,
        "generator_family": dataset.condition.generator_family,
    }
    with (output / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(summary)

    rows = []
    for chrom in range(dataset.condition.chromosome_count):
        selected = dataset.chrom_index == chrom
        row = {
            "chromosome": f"chr{chrom + 1}",
            "geometry_selected_copy_flip": flips[chrom],
            "run_id": run_id,
            "run_contract_sha256": run_contract_sha256,
            "seed": dataset.condition.seed,
            "generator_family": dataset.condition.generator_family,
        }
        correlations = []
        for copy in (0, 1):
            truth_distance = pdist(dataset.truth_coordinates[selected, copy])
            reconstruction_distance = pdist(estimate[selected, copy])
            correlation = float(spearmanr(truth_distance, reconstruction_distance).statistic)
            correlations.append(correlation)
            row[f"copy{copy}_cis_distance_spearman"] = correlation
            row[f"truth_copy{copy}_volume"] = volume(dataset.truth_coordinates[selected, copy])
            row[f"reconstruction_copy{copy}_volume"] = volume(estimate[selected, copy])
        row["mean_cis_distance_spearman"] = float(np.mean(correlations))
        row["truth_mean_homolog_separation"] = float(np.linalg.norm(
            dataset.truth_coordinates[selected, 0] - dataset.truth_coordinates[selected, 1], axis=1
        ).mean())
        row["reconstruction_mean_homolog_separation"] = float(np.linalg.norm(
            estimate[selected, 0] - estimate[selected, 1], axis=1
        ).mean())
        rows.append(row)
    with (output / "per_chromosome.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    require_analysis()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    })
    condition = SyntheticCondition(
        generator_family="well_specified", homolog_separation=2.0,
        contacts_per_cell=3000, cis_fraction=0.65, background_fraction=0.05,
        resolution_bp=1_000_000, chromosome_count=3, chromosome_length_bins=10,
        seed=args.seed,
    )
    dataset = generate_synthetic(condition)
    inference = InferenceConfig(
        likelihood_mode="monotone_log", seed=args.seed, iterations=130,
        learning_rate=0.025,
    )
    run_contract = {
        "condition": asdict(condition),
        "inference": asdict(inference),
        "truth_used_in_inference": "none",
    }
    serialized_contract = json.dumps(run_contract, sort_keys=True, separators=(",", ":"))
    run_contract_sha256 = hashlib.sha256(serialized_contract.encode("utf-8")).hexdigest()
    run_id = f"phase059_standard_{run_contract_sha256[:16]}"
    result = infer_synthetic(dataset, inference)
    estimate = result.coordinates.copy()
    flips: dict[int, int] = {}
    for chrom in range(condition.chromosome_count):
        selected = dataset.chrom_index == chrom
        _, flip = chromosome_swap_invariant_metrics(
            dataset.truth_coordinates[selected], estimate[selected]
        )
        flips[chrom] = flip
        if flip:
            estimate[selected] = estimate[selected, ::-1]
    estimate = align_to_reference(dataset.truth_coordinates, estimate)

    plot_all_chromosomes(
        args.output_dir / "all_chrom_3d_scatter.png", dataset,
        dataset.truth_coordinates, estimate,
    )
    plot_chr1(
        args.output_dir / "chr1_copy_3d_scatter.png", dataset,
        dataset.truth_coordinates, estimate,
    )
    plot_chr1_distance_maps(
        args.output_dir / "chr1_distance_maps.png", dataset,
        dataset.truth_coordinates, estimate,
    )
    write_tables(
        args.output_dir.parent, dataset, result, estimate, flips,
        run_id=run_id, run_contract_sha256=run_contract_sha256,
    )
    (args.output_dir.parent / "standard_plot_config.json").write_text(
        json.dumps({
            **run_contract,
            "run_id": run_id,
            "run_contract_sha256": run_contract_sha256,
            "alignment": "per-chromosome ORACLE geometry swap, then global Procrustes; evaluation only",
        }, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
