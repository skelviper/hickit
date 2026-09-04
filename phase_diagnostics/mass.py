"""Observation and molecule mass accounting for the diploid M-step."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class MassResult:
    state_mass: np.ndarray
    target_distance: np.ndarray
    force_coefficient: np.ndarray

    @property
    def total_state_mass(self) -> float:
        return float(self.state_mass.sum())

    @property
    def total_force_coefficient(self) -> float:
        return float(self.force_coefficient.sum())


def mstep_contribution(
    raw_contact_mass: float,
    posterior: Sequence[float],
    *,
    contact_probability: float = 1.0,
    confidence_mass: float = 1.0,
    target_mode: str = "posterior_once",
    posterior_gamma: float = 1.0,
    epsilon_count: float = 1e-6,
) -> MassResult:
    """Calculate the contribution of one raw observation to four state edges.

    ``posterior_once`` matches the conservative raw-count target-distance mode.
    ``legacy_posterior_target`` reproduces the duplicated posterior dependence in
    which posterior controls both edge weight and target distance.
    """
    if not np.isfinite(raw_contact_mass) or raw_contact_mass < 0:
        raise ValueError("raw_contact_mass must be finite and non-negative")
    if not 0 <= contact_probability <= 1 or not 0 <= confidence_mass <= 1:
        raise ValueError("contact_probability and confidence_mass must be in [0,1]")
    p = np.asarray(posterior, dtype=float)
    if p.shape != (4,) or not np.isfinite(p).all() or (p < 0).any() or p.sum() <= 0:
        raise ValueError("posterior must contain four finite non-negative values")
    p = p / p.sum()
    state_mass = raw_contact_mass * contact_probability * confidence_mass * p
    if target_mode == "posterior_once":
        effective_count = np.full(4, max(raw_contact_mass, epsilon_count), dtype=float)
    elif target_mode == "legacy_posterior_target":
        effective_count = np.maximum(raw_contact_mass * np.power(p, posterior_gamma), epsilon_count)
    else:
        raise ValueError(f"unsupported target mode: {target_mode}")
    target = np.power(effective_count, -1.0 / 3.0)
    force = state_mass / target
    return MassResult(state_mass=state_mass, target_distance=target, force_coefficient=force)


def normalized_pair_expansion(number_of_segments: int, source_mass: float = 1.0) -> np.ndarray:
    if number_of_segments < 2:
        raise ValueError("a molecule must contain at least two segments")
    if not np.isfinite(source_mass) or source_mass < 0:
        raise ValueError("source_mass must be finite and non-negative")
    number_of_pairs = number_of_segments * (number_of_segments - 1) // 2
    return np.full(number_of_pairs, source_mass / number_of_pairs, dtype=float)


def assert_mass_conserved(values: Sequence[float], expected: float, tolerance: float = 1e-9) -> None:
    observed = float(np.asarray(values, dtype=float).sum())
    if not np.isclose(observed, expected, atol=tolerance, rtol=tolerance):
        raise AssertionError(f"mass is not conserved: observed={observed:.17g}, expected={expected:.17g}")
