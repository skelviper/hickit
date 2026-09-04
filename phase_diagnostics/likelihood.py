"""Contact observation scores kept separate from the FDG layout potential."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np


SCORE_MODES = (
    "fdg_flat",
    "dist2",
    "logdist2",
    "monotonic_powerlaw",
    "monotone_log",
    "monotone_logistic",
)


@dataclass(frozen=True)
class ScoreConfig:
    mode: str = "monotone_log"
    temperature: float = 1.0
    k: float = 1.0
    alpha: float = 2.0
    d0: float = 0.05
    r0: float = 1.0
    scale: float = 0.25
    unit: float = 1.0
    d_scale: float = 1.0
    fdg_d_c1: float = 0.5
    fdg_d_c2: float = 1.5
    fdg_d_c3: float = 2.0

    def validate(self) -> None:
        if self.mode not in SCORE_MODES:
            raise ValueError(f"unsupported score mode: {self.mode}")
        for name in ("temperature", "k", "alpha", "d0", "scale", "unit", "d_scale"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if not math.isfinite(self.r0):
            raise ValueError("r0 must be finite")


def fdg_flat_energy(distance: np.ndarray | float, config: ScoreConfig) -> np.ndarray:
    """Exact shell-shaped contact term used by the current Hickit FDG code."""
    r = np.asarray(distance, dtype=float) / (config.unit * config.d_scale)
    d1, d2, d3 = config.fdg_d_c1, config.fdg_d_c2, config.fdg_d_c3
    c1 = 3.0 * (d3 - d2)
    c2 = (d3 - d2) ** 3
    energy = np.empty_like(r)
    close = r < d1
    flat = (r >= d1) & (r <= d2)
    shoulder = (r > d2) & (r <= d3)
    far = r > d3
    energy[close] = config.k * np.square(d1 - r[close])
    energy[flat] = 0.0
    energy[shoulder] = config.k * np.square(r[shoulder] - d2)
    energy[far] = config.k * (c1 * (r[far] - d3) + c2 / (r[far] - d2))
    return energy


def log_score(distance: np.ndarray | float, config: ScoreConfig) -> np.ndarray:
    """Return an unnormalized log contact score; higher means more likely."""
    config.validate()
    d = np.asarray(distance, dtype=float)
    if not np.isfinite(d).all() or (d < 0).any():
        raise ValueError("distances must be finite and non-negative")
    r = d / (config.unit * config.d_scale)
    if config.mode == "fdg_flat":
        value = -fdg_flat_energy(d, config)
    elif config.mode == "dist2":
        value = -config.k * np.square(r)
    elif config.mode == "logdist2":
        value = -config.k * np.log1p(np.square(r) / 1e-6)
    elif config.mode == "monotonic_powerlaw":
        value = -config.k * np.log1p(np.square(r))
    elif config.mode == "monotone_log":
        value = -config.k * config.alpha * np.log(r + config.d0)
    elif config.mode == "monotone_logistic":
        # log(sigmoid((r0-r)/scale)) = -log(1 + exp((r-r0)/scale)).
        value = -config.k * np.logaddexp(0.0, (r - config.r0) / config.scale)
    else:
        raise AssertionError(config.mode)
    return value / config.temperature


def _logsumexp(values: np.ndarray) -> float:
    maximum = float(np.max(values))
    if not math.isfinite(maximum):
        raise ValueError("all candidate scores are non-finite")
    return maximum + math.log(float(np.exp(values - maximum).sum()))


def posterior_from_distances(
    distances: Sequence[float],
    config: ScoreConfig,
    prior: Sequence[float] | None = None,
    background_prior: float = 0.0,
    background_log_score: float = 0.0,
) -> tuple[np.ndarray, float, float]:
    """Return four signal probabilities, background probability, and log normalizer.

    The four returned signal probabilities include their mixture mass and therefore
    sum to ``1 - q_background``. This prevents the background component from being
    silently renormalized away before the M-step.
    """
    d = np.asarray(distances, dtype=float)
    if d.shape != (4,):
        raise ValueError("exactly four state distances are required")
    if prior is None:
        p = np.full(4, 0.25, dtype=float)
    else:
        p = np.asarray(prior, dtype=float)
        if p.shape != (4,) or not np.isfinite(p).all() or (p < 0).any() or p.sum() <= 0:
            raise ValueError("prior must contain four finite non-negative values")
        p = p / p.sum()
    if not 0.0 <= background_prior < 1.0:
        raise ValueError("background_prior must be in [0, 1)")
    signal_prior_mass = 1.0 - background_prior
    scores = log_score(d, config) + np.log(np.maximum(p * signal_prior_mass, 1e-300)) / config.temperature
    if background_prior > 0.0:
        bg_score = math.log(background_prior) + float(background_log_score) / config.temperature
        all_scores = np.concatenate((scores, np.asarray([bg_score], dtype=float)))
    else:
        all_scores = scores
    normalizer = _logsumexp(all_scores)
    posterior = np.exp(scores - normalizer)
    q_background = math.exp(all_scores[4] - normalizer) if all_scores.size == 5 else 0.0
    if not np.isclose(float(posterior.sum()) + q_background, 1.0, atol=1e-12):
        raise AssertionError("posterior mixture is not normalized")
    return posterior, q_background, normalizer


def score_curve(
    config: ScoreConfig, distances: Sequence[float], reference_distance: float | None = None,
) -> list[dict[str, float | str]]:
    values = np.asarray(distances, dtype=float)
    scores = log_score(values, config)
    raw_energy = -scores * config.temperature
    if reference_distance is None:
        reference_index = 0
    else:
        reference_index = int(np.argmin(np.abs(values - reference_distance)))
    relative = np.exp(np.clip(scores - scores[reference_index], -700, 700))
    return [
        {
            "score_mode": config.mode,
            "temperature": config.temperature,
            "distance": float(distance),
            "raw_energy": float(energy),
            "log_likelihood_like_score": float(score),
            "normalized_relative_weight": float(weight),
        }
        for distance, energy, score, weight in zip(values, raw_energy, scores, relative, strict=True)
    ]


def is_nonincreasing(config: ScoreConfig, distances: Sequence[float], tolerance: float = 1e-12) -> bool:
    values = log_score(np.asarray(distances, dtype=float), config)
    return bool(np.all(np.diff(values) <= tolerance))
