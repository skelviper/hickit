"""Central metric definitions for contact, structure, gauge, and selection audits."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np
from scipy.spatial.distance import pdist
from scipy.stats import kendalltau, pearsonr, spearmanr

from .state import parity_probabilities, relabel_probabilities


EPS = 1e-15


def _wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total < 1:
        return math.nan, math.nan
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    half = z * math.sqrt(
        proportion * (1.0 - proportion) / total + z * z / (4.0 * total * total)
    ) / denominator
    return center - half, center + half


@dataclass(frozen=True)
class ReliabilityBin:
    metric: str
    bin_index: int
    lower: float
    upper: float
    count: int
    coverage: float
    mean_confidence: float
    accuracy: float


def _validate_contact_inputs(
    truth: Sequence[int], probabilities: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(truth, dtype=int)
    q = np.asarray(probabilities, dtype=float)
    if y.ndim != 1 or q.shape != (y.size, 4):
        raise ValueError("truth must be N and probabilities must be Nx4")
    if y.size == 0:
        raise ValueError("at least one evaluated contact is required")
    if (y < 0).any() or (y >= 4).any():
        raise ValueError("truth states must be in [0,3]")
    if not np.isfinite(q).all() or (q < 0).any() or (q.sum(axis=1) > 1 + 1e-8).any():
        raise ValueError("state probabilities must be finite, non-negative, and have row sum <= 1")
    return y, q


def reliability_curve(
    correct: Sequence[bool], confidence: Sequence[float], *, bins: int = 10,
    metric: str = "exact",
) -> list[ReliabilityBin]:
    observed = np.asarray(correct, dtype=bool)
    conf = np.asarray(confidence, dtype=float)
    if observed.shape != conf.shape or observed.ndim != 1:
        raise ValueError("correct and confidence must be matching vectors")
    if not np.isfinite(conf).all() or (conf < 0).any() or (conf > 1 + 1e-10).any():
        raise ValueError("confidence must be finite and in [0,1]")
    edges = np.linspace(0.0, 1.0, bins + 1)
    result: list[ReliabilityBin] = []
    for index in range(bins):
        lower, upper = float(edges[index]), float(edges[index + 1])
        selected = (conf >= lower) & ((conf < upper) if index + 1 < bins else (conf <= upper))
        count = int(selected.sum())
        result.append(
            ReliabilityBin(
                metric=metric,
                bin_index=index,
                lower=lower,
                upper=upper,
                count=count,
                coverage=count / max(1, observed.size),
                mean_confidence=float(conf[selected].mean()) if count else math.nan,
                accuracy=float(observed[selected].mean()) if count else math.nan,
            )
        )
    return result


def expected_calibration_error(curve: Iterable[ReliabilityBin]) -> float:
    return float(
        sum(row.coverage * abs(row.accuracy - row.mean_confidence) for row in curve if row.count)
    )


def selective_coverage(
    correct: Sequence[bool], confidence: Sequence[float],
    precision_targets: Sequence[float] = (0.55, 0.60, 0.70, 0.80, 0.90),
) -> list[dict[str, float | int]]:
    curve = selective_call_curve(correct, confidence)
    rows: list[dict[str, float | int]] = []
    for target in precision_targets:
        eligible = [row for row in curve if float(row["precision"]) >= target]
        if eligible:
            selected = max(eligible, key=lambda row: int(row["callable_count"]))
            rows.append(
                {
                    "target_precision": float(target),
                    "callable_count": int(selected["callable_count"]),
                    "callable_fraction": float(selected["coverage"]),
                    "achieved_precision": float(selected["precision"]),
                    "minimum_confidence": float(selected["minimum_confidence"]),
                }
            )
        else:
            rows.append(
                {
                    "target_precision": float(target),
                    "callable_count": 0,
                    "callable_fraction": 0.0,
                    "achieved_precision": math.nan,
                    "minimum_confidence": math.nan,
                }
            )
    return rows


def selective_call_curve(
    correct: Sequence[bool], confidence: Sequence[float],
) -> list[dict[str, float | int]]:
    """Return precision/coverage only at realizable confidence thresholds.

    Every row includes the complete tie group at its threshold. This prevents
    input order from creating apparently high precision by splitting calls
    that have exactly the same confidence.
    """
    observed = np.asarray(correct, dtype=bool)
    conf = np.asarray(confidence, dtype=float)
    if observed.ndim != 1 or observed.shape != conf.shape or observed.size == 0:
        raise ValueError("correct and confidence must be matching non-empty vectors")
    if not np.isfinite(conf).all() or (conf < 0).any() or (conf > 1 + 1e-10).any():
        raise ValueError("confidence must be finite and in [0,1]")
    order = np.argsort(-conf, kind="stable")
    ranked_correct = observed[order].astype(float)
    ranked_conf = conf[order]
    cumulative_correct = np.cumsum(ranked_correct)
    # A threshold cannot select only part of an exact confidence tie. Evaluate
    # the curve at the final index of each tie group.
    group_ends = np.flatnonzero(np.r_[ranked_conf[:-1] != ranked_conf[1:], True])
    total_correct = max(1.0, float(observed.sum()))
    rows = []
    for end in group_ends:
        stop = int(end) + 1
        correct_count = int(cumulative_correct[end])
        rows.append(
            {
                "callable_count": stop,
                "correct_count": correct_count,
                "coverage": stop / observed.size,
                "precision": correct_count / stop,
                "correct_recall": correct_count / total_correct,
                "minimum_confidence": float(ranked_conf[end]),
            }
        )
    return rows


def selective_threshold_curve(
    correct: Sequence[bool], confidence: Sequence[float],
    thresholds: Sequence[float] = tuple(np.linspace(0.0, 1.0, 21)),
) -> list[dict[str, float | int]]:
    """Evaluate selective-call accuracy and recall at fixed confidence thresholds."""
    observed = np.asarray(correct, dtype=bool)
    conf = np.asarray(confidence, dtype=float)
    if observed.ndim != 1 or observed.shape != conf.shape or observed.size == 0:
        raise ValueError("correct and confidence must be matching non-empty vectors")
    if not np.isfinite(conf).all() or (conf < 0).any() or (conf > 1 + 1e-10).any():
        raise ValueError("confidence must be finite and in [0,1]")
    total_correct = int(observed.sum())
    rows: list[dict[str, float | int]] = []
    for threshold in thresholds:
        selected = conf >= float(threshold)
        count = int(selected.sum())
        correct_count = int(observed[selected].sum())
        rows.append({
            "minimum_confidence": float(threshold),
            "callable_count": count,
            "correct_count": correct_count,
            "coverage": count / observed.size,
            "precision": correct_count / count if count else math.nan,
            "correct_recall": correct_count / total_correct if total_correct else math.nan,
        })
    return rows


def contact_metrics(
    truth: Sequence[int], probabilities: np.ndarray, *, calibration_bins: int = 10,
) -> tuple[dict[str, float | int], list[ReliabilityBin], list[dict[str, float | int]]]:
    y, q = _validate_contact_inputs(truth, probabilities)
    row_sums = q.sum(axis=1)
    if (row_sums <= EPS).any():
        raise ValueError("four-state signal probability must be positive")
    # Contact-state metrics are conditional on the signal component. A model
    # with background mass must calibrate that fifth state separately.
    conditional = q / row_sums[:, None]
    predicted = np.argmax(conditional, axis=1)
    confidence = np.max(conditional, axis=1)
    correct = predicted == y
    exact_ties = np.isclose(conditional, confidence[:, None], rtol=0.0, atol=1e-15)
    exact_fractional_correct = exact_ties[np.arange(y.size), y] / exact_ties.sum(axis=1)
    truth_prob = conditional[np.arange(y.size), y]
    one_hot = np.eye(4, dtype=float)[y]

    parity_truth = ((y == 1) | (y == 2)).astype(int)
    parity_q = parity_probabilities(conditional)
    parity_pred = np.argmax(parity_q, axis=1)
    parity_conf = np.max(parity_q, axis=1)
    parity_correct = parity_pred == parity_truth
    parity_ties = np.isclose(parity_q, parity_conf[:, None], rtol=0.0, atol=1e-15)
    parity_fractional_correct = (
        parity_ties[np.arange(y.size), parity_truth] / parity_ties.sum(axis=1)
    )
    parity_truth_prob = parity_q[np.arange(y.size), parity_truth]
    parity_one_hot = np.eye(2, dtype=float)[parity_truth]

    exact_curve = reliability_curve(correct, confidence, bins=calibration_bins, metric="exact")
    parity_curve = reliability_curve(
        parity_correct, parity_conf, bins=calibration_bins, metric="same_cross"
    )
    entropy = -np.sum(
        np.where(conditional > 0, conditional * np.log(np.maximum(conditional, EPS)), 0.0),
        axis=1,
    )
    exact_ci = _wilson_interval(int(correct.sum()), y.size)
    parity_ci = _wilson_interval(int(parity_correct.sum()), y.size)
    result: dict[str, float | int] = {
        "n_contacts": int(y.size),
        "callable_fraction": 1.0,
        "exact_top1_accuracy": float(correct.mean()),
        "exact_random_tie_expected_accuracy": float(exact_fractional_correct.mean()),
        "exact_accuracy_ci95_lower": exact_ci[0],
        "exact_accuracy_ci95_upper": exact_ci[1],
        "exact_log_loss": float(-np.log(np.maximum(truth_prob, EPS)).mean()),
        "exact_brier": float(np.square(conditional - one_hot).sum(axis=1).mean()),
        "same_cross_accuracy": float(parity_correct.mean()),
        "same_cross_random_tie_expected_accuracy": float(parity_fractional_correct.mean()),
        "same_cross_accuracy_ci95_lower": parity_ci[0],
        "same_cross_accuracy_ci95_upper": parity_ci[1],
        "same_cross_log_loss": float(-np.log(np.maximum(parity_truth_prob, EPS)).mean()),
        "same_cross_brier": float(np.square(parity_q - parity_one_hot).sum(axis=1).mean()),
        "mean_entropy": float(entropy.mean()),
        "exact_ece": expected_calibration_error(exact_curve),
        "same_cross_ece": expected_calibration_error(parity_curve),
        "mean_exact_confidence": float(confidence.mean()),
        "mean_same_cross_confidence": float(parity_conf.mean()),
        "mean_signal_probability": float(row_sums.mean()),
    }
    curves = [*exact_curve, *parity_curve]
    selective: list[dict[str, float | int]] = []
    for name, corr, conf in (
        ("exact", correct, confidence),
        ("same_cross", parity_correct, parity_conf),
        ("exact_signal_aware", correct, np.max(q, axis=1)),
        ("same_cross_signal_aware", parity_correct, np.max(parity_probabilities(q), axis=1)),
    ):
        for row in selective_coverage(corr, conf):
            selective.append({"metric": name, **row})
    return result, curves, selective


def stratified_contact_metrics(
    truth: Sequence[int], probabilities: np.ndarray, strata: Sequence[str],
) -> list[dict[str, float | int | str]]:
    y, q = _validate_contact_inputs(truth, probabilities)
    labels = np.asarray(strata, dtype=object)
    if labels.shape != y.shape:
        raise ValueError("strata must have one label per contact")
    rows: list[dict[str, float | int | str]] = []
    for label in sorted(set(str(item) for item in labels)):
        selected = labels == label
        metrics, _, _ = contact_metrics(y[selected], q[selected])
        rows.append({"stratum": label, **metrics})
    return rows


def _kabsch_align(reference: np.ndarray, estimate: np.ndarray, scale: bool = True) -> np.ndarray:
    if reference.shape != estimate.shape or reference.ndim != 2 or reference.shape[1] != 3:
        raise ValueError("reference and estimate must be matching Nx3 arrays")
    ref = reference - reference.mean(axis=0, keepdims=True)
    est = estimate - estimate.mean(axis=0, keepdims=True)
    covariance = est.T @ ref
    left, singular, right_t = np.linalg.svd(covariance)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    aligned = est @ rotation
    if scale:
        denominator = float(np.square(est).sum())
        factor = float(singular.sum() / denominator) if denominator > EPS else 1.0
        aligned *= factor
    return aligned + reference.mean(axis=0, keepdims=True)


def _safe_corr(left: np.ndarray, right: np.ndarray, method: str) -> float:
    if left.size < 2 or right.size != left.size or np.ptp(left) <= EPS or np.ptp(right) <= EPS:
        return math.nan
    value = pearsonr(left, right).statistic if method == "pearson" else spearmanr(left, right).statistic
    return float(value)


def structure_metrics(reference: np.ndarray, estimate: np.ndarray) -> dict[str, float | int]:
    ref = np.asarray(reference, dtype=float)
    est = np.asarray(estimate, dtype=float)
    aligned = _kabsch_align(ref.reshape(-1, 3), est.reshape(-1, 3)).reshape(ref.shape)
    rigid_aligned = _kabsch_align(
        ref.reshape(-1, 3), est.reshape(-1, 3), scale=False
    ).reshape(ref.shape)
    ref_dist = pdist(ref.reshape(-1, 3))
    est_dist = pdist(est.reshape(-1, 3))
    residual = aligned - ref
    return {
        "n_points": int(ref.reshape(-1, 3).shape[0]),
        "procrustes_rmsd": float(np.sqrt(np.square(residual).sum(axis=-1).mean())),
        "rigid_rmsd": float(np.sqrt(np.square(rigid_aligned - ref).sum(axis=-1).mean())),
        "distance_pearson": _safe_corr(ref_dist, est_dist, "pearson"),
        "distance_spearman": _safe_corr(ref_dist, est_dist, "spearman"),
    }


def chromosome_swap_invariant_metrics(
    reference: np.ndarray, estimate: np.ndarray,
) -> tuple[dict[str, float | int], int]:
    """Evaluate one chromosome with arrays shaped bins x copies x xyz."""
    ref = np.asarray(reference, dtype=float)
    est = np.asarray(estimate, dtype=float)
    if ref.ndim != 3 or ref.shape[1:] != (2, 3) or est.shape != ref.shape:
        raise ValueError("chromosome coordinates must be bins x 2 x 3")
    direct = structure_metrics(ref, est)
    swapped = structure_metrics(ref, est[:, ::-1, :])
    if float(swapped["procrustes_rmsd"]) < float(direct["procrustes_rmsd"]):
        return swapped, 1
    return direct, 0


def homolog_separation_correlation(reference: np.ndarray, estimate: np.ndarray) -> float:
    ref_sep = np.linalg.norm(reference[:, 0] - reference[:, 1], axis=-1)
    est_sep = np.linalg.norm(estimate[:, 0] - estimate[:, 1], axis=-1)
    return _safe_corr(ref_sep, est_sep, "spearman")


def gauge_metrics(
    truth_flips: Mapping[str, int], inferred_flips: Mapping[str, int],
    edges: Iterable[tuple[str, str]] | None = None,
) -> dict[str, float | int]:
    chromosomes = sorted(set(truth_flips) & set(inferred_flips))
    if not chromosomes:
        return {"n_chromosomes": 0, "per_chromosome_flip_accuracy": math.nan, "relative_flip_accuracy": math.nan}
    # A global flip of every chromosome is itself gauge-equivalent.
    direct = np.mean([truth_flips[c] == inferred_flips[c] for c in chromosomes])
    inverted = np.mean([truth_flips[c] != inferred_flips[c] for c in chromosomes])
    pairs = list(edges) if edges is not None else [
        (left, right) for index, left in enumerate(chromosomes) for right in chromosomes[index + 1 :]
    ]
    relative = [
        (truth_flips[left] ^ truth_flips[right]) == (inferred_flips[left] ^ inferred_flips[right])
        for left, right in pairs if left in truth_flips and right in truth_flips and left in inferred_flips and right in inferred_flips
    ]
    return {
        "n_chromosomes": len(chromosomes),
        "per_chromosome_flip_accuracy": float(max(direct, inverted)),
        "relative_flip_accuracy": float(np.mean(relative)) if relative else math.nan,
    }


def model_selection_metrics(
    blind_score: Sequence[float], truth_score: Sequence[float], *, higher_blind_is_better: bool = True,
    higher_truth_is_better: bool = True, top_k: int = 1,
) -> dict[str, float]:
    blind = np.asarray(blind_score, dtype=float)
    truth = np.asarray(truth_score, dtype=float)
    if blind.shape != truth.shape or blind.ndim != 1 or blind.size < 2:
        raise ValueError("blind and truth scores must be matching vectors with at least two seeds")
    blind_rank_value = blind if higher_blind_is_better else -blind
    truth_rank_value = truth if higher_truth_is_better else -truth
    selected = np.argsort(-blind_rank_value, kind="stable")[:top_k]
    best_selected_truth = float(np.max(truth_rank_value[selected]))
    optimal_truth = float(np.max(truth_rank_value))
    quartile = float(np.quantile(truth_rank_value, 0.75))
    return {
        "spearman": float(spearmanr(blind_rank_value, truth_rank_value).statistic),
        "kendall": float(kendalltau(blind_rank_value, truth_rank_value).statistic),
        "top_k_selection_regret": optimal_truth - best_selected_truth,
        "selected_in_truth_top_quartile": float(best_selected_truth >= quartile),
    }


def optimal_chromosome_flip_agreement(
    left: np.ndarray, right: np.ndarray, chrom1: Sequence[str], chrom2: Sequence[str],
    flips: Mapping[str, int],
) -> float:
    q_left = np.asarray(left, dtype=float)
    q_right = np.asarray(right, dtype=float)
    if q_left.shape != q_right.shape or q_left.shape[1] != 4:
        raise ValueError("posterior matrices must be matching Nx4 arrays")
    aligned = np.vstack([
        relabel_probabilities(q_right[index], flip_left=flips.get(str(c1), 0), flip_right=flips.get(str(c2), 0))
        for index, (c1, c2) in enumerate(zip(chrom1, chrom2, strict=True))
    ])
    return float(np.mean(np.argmax(q_left, axis=1) == np.argmax(aligned, axis=1)))


def jensen_shannon_divergence(left: np.ndarray, right: np.ndarray) -> float:
    p = np.asarray(left, dtype=float)
    q = np.asarray(right, dtype=float)
    if p.shape != q.shape or p.shape[-1] != 4:
        raise ValueError("posterior arrays must be matching and end in four states")
    p = p / np.maximum(p.sum(axis=-1, keepdims=True), EPS)
    q = q / np.maximum(q.sum(axis=-1, keepdims=True), EPS)
    midpoint = 0.5 * (p + q)
    kl_p = np.sum(np.where(p > 0, p * np.log(np.maximum(p / midpoint, EPS)), 0.0), axis=-1)
    kl_q = np.sum(np.where(q > 0, q * np.log(np.maximum(q / midpoint, EPS)), 0.0), axis=-1)
    return float(np.mean(0.5 * (kl_p + kl_q)))
