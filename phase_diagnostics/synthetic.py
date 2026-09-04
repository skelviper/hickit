"""Compact synthetic phase-diagram generator and blind midpoint/delta solver."""

from __future__ import annotations

import itertools
import math
import time
from dataclasses import asdict, dataclass, replace
from typing import Sequence

import numpy as np

from .likelihood import ScoreConfig, log_score, posterior_from_distances
from .metrics import (
    chromosome_swap_invariant_metrics,
    contact_metrics,
    expected_calibration_error,
    homolog_separation_correlation,
    reliability_curve,
)
from .state import flip_state, relabel_probabilities
from .z2 import RelativeFlipEdge, synchronize_z2


@dataclass(frozen=True)
class SyntheticCondition:
    generator_family: str = "well_specified"
    homolog_separation: float = 1.0
    contacts_per_cell: int = 4000
    cis_fraction: float = 0.7
    background_fraction: float = 0.05
    resolution_bp: int = 1_000_000
    chromosome_count: int = 3
    chromosome_length_bins: int = 12
    readchain_length: int = 2
    molecule_count: int = 0
    count_overdispersion: float = 0.0
    anchor_fraction: float = 0.0
    seed: int = 17

    def validate(self) -> None:
        if self.generator_family not in {"well_specified", "model_mismatched"}:
            raise ValueError("generator_family must be well_specified or model_mismatched")
        if self.homolog_separation < 0 or self.contacts_per_cell < 1:
            raise ValueError("homolog separation and contact depth are invalid")
        for name in ("cis_fraction", "background_fraction", "anchor_fraction"):
            value = getattr(self, name)
            if not 0 <= value <= 1:
                raise ValueError(f"{name} must be in [0,1]")
        if self.resolution_bp < 1 or self.chromosome_count < 1 or self.chromosome_length_bins < 4:
            raise ValueError("resolution/chromosome dimensions are invalid")
        if self.readchain_length < 2 or self.molecule_count < 0 or self.count_overdispersion < 0:
            raise ValueError("readchain or overdispersion setting is invalid")


@dataclass
class SyntheticDataset:
    condition: SyntheticCondition
    truth_coordinates: np.ndarray
    chrom_index: np.ndarray
    bin_offset: np.ndarray
    pair_i: np.ndarray
    pair_j: np.ndarray
    truth_state: np.ndarray
    is_background: np.ndarray
    anchor_mask: np.ndarray
    molecule_id: np.ndarray
    molecule_segments: tuple[np.ndarray, ...]
    molecule_truth_copies: tuple[np.ndarray, ...]

    @property
    def n_bins(self) -> int:
        return int(self.truth_coordinates.shape[0])


@dataclass(frozen=True)
class InferenceConfig:
    seed: int = 17
    iterations: int = 250
    learning_rate: float = 0.03
    likelihood_mode: str = "monotone_log"
    alpha: float = 2.0
    d0: float = 0.1
    negative_ratio: float = 1.0
    regularization_reference_mass: float = 1200.0
    backbone_weight: float = 3.0
    bending_weight: float = 0.2
    delta_smooth_weight: float = 0.3
    delta_shrink_weight: float = 0.01
    scale_weight: float = 0.5
    observation_model: str = "auto"
    background_component: bool = True
    background_fraction_min: float = 0.001
    background_fraction_max: float = 0.5
    staged_cis_then_trans: bool = False
    stage_shape_weight: float = 2.0
    joint_refine_shape_weight: float = 0.25
    molecule_likelihood: bool = False
    global_copy_track_weight: float = 0.0
    minimum_homolog_separation_weight: float = 0.0
    minimum_homolog_separation: float = 1.0
    device: str = "cpu"


@dataclass
class InferenceResult:
    coordinates: np.ndarray
    state_probabilities: np.ndarray
    background_probability: np.ndarray
    objective_history: list[float]
    runtime_seconds: float
    converged: bool
    anchor_flips: dict[int, int]
    signal_fraction: float = 1.0


def synchronize_result_z2(
    dataset: SyntheticDataset,
    result: InferenceResult,
    relative_edges: Sequence[RelativeFlipEdge] | None = None,
) -> tuple[InferenceResult, dict[str, float | int | str]]:
    """Apply only calibrated relative-gauge edges supplied by a blind procedure.

    Four-state posteriors alone cannot calibrate a chromosome-pair flip: their
    same/cross balance changes under the very gauge being estimated. Therefore
    an ordinary SNP-free call supplies no edges and remains unresolved rather
    than imposing a hidden same-homolog trans-contact prior.
    """
    edges = list(relative_edges or ())
    solved = synchronize_z2(
        edges, chromosomes=[f"chr{chrom + 1}" for chrom in range(dataset.condition.chromosome_count)]
    )
    flips = {
        chrom: int(
            f"chr{chrom + 1}" not in solved.unresolved
            and solved.spins.get(f"chr{chrom + 1}", 1) < 0
        )
        for chrom in range(dataset.condition.chromosome_count)
    }
    q4 = np.vstack([
        relabel_probabilities(
            row,
            flip_left=flips[int(dataset.chrom_index[dataset.pair_i[index]])],
            flip_right=flips[int(dataset.chrom_index[dataset.pair_j[index]])],
        )
        for index, row in enumerate(result.state_probabilities)
    ])
    coordinates = result.coordinates.copy()
    for bin_index, chrom in enumerate(dataset.chrom_index):
        if flips[int(chrom)]:
            coordinates[bin_index] = coordinates[bin_index, ::-1]
    synchronized = InferenceResult(
        coordinates=coordinates,
        state_probabilities=q4,
        background_probability=result.background_probability.copy(),
        objective_history=list(result.objective_history),
        runtime_seconds=result.runtime_seconds,
        converged=result.converged,
        anchor_flips={**result.anchor_flips, **flips},
        signal_fraction=result.signal_fraction,
    )
    diagnostics: dict[str, float | int | str] = {
        "z2_evidence_status": (
            "CALIBRATED_BLIND_RELATIVE_EDGES"
            if edges else "UNRESOLVED_NO_CALIBRATED_BLIND_RELATIVE_GAUGE_EVIDENCE"
        ),
        "z2_n_edges": len(edges),
        "z2_n_components": len(set(solved.component.values())),
        "z2_largest_component": max(
            (sum(component == value for component in solved.component.values()) for value in set(solved.component.values())),
            default=0,
        ),
        "z2_frustrated_edge_fraction": solved.frustrated_edge_fraction,
        "z2_weighted_frustrated_edge_fraction": solved.weighted_frustrated_edge_fraction,
        "z2_cycle_inconsistency": solved.cycle_inconsistency,
        "z2_unresolved_chromosome_fraction": len(solved.unresolved) / max(1, dataset.condition.chromosome_count),
        "z2_min_spin_margin": min(solved.margins.values(), default=0.0),
    }
    return synchronized, diagnostics


def _unit_vector(rng: np.random.Generator) -> np.ndarray:
    while True:
        value = rng.normal(size=3)
        norm = float(np.linalg.norm(value))
        if norm > 1e-12:
            return value / norm


def _polymer(rng: np.random.Generator, n: int, origin: np.ndarray, persistence: float = 0.75) -> np.ndarray:
    values = np.empty((n, 3), dtype=float)
    values[0] = origin
    direction = _unit_vector(rng)
    for index in range(1, n):
        direction = persistence * direction + math.sqrt(1 - persistence**2) * _unit_vector(rng)
        direction /= np.linalg.norm(direction)
        values[index] = values[index - 1] + direction
    return values


def _geometry(condition: SyntheticCondition, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    coords: list[np.ndarray] = []
    chrom_index: list[int] = []
    offsets: list[int] = []
    for chrom in range(condition.chromosome_count):
        center = rng.normal(scale=4.0, size=3)
        midpoint = _polymer(rng, condition.chromosome_length_bins, center)
        direction = _unit_vector(rng)
        delta = np.empty_like(midpoint)
        current = direction
        for index in range(condition.chromosome_length_bins):
            current = 0.9 * current + 0.1 * _unit_vector(rng)
            current /= np.linalg.norm(current)
            delta[index] = 0.5 * condition.homolog_separation * current
        # Mild copy-specific shape variation prevents a pure rigid-copy toy.
        variation = np.cumsum(rng.normal(scale=0.03, size=midpoint.shape), axis=0)
        copies = np.stack((midpoint - delta + variation, midpoint + delta - variation), axis=1)
        coords.append(copies)
        chrom_index.extend([chrom] * condition.chromosome_length_bins)
        offsets.extend(range(condition.chromosome_length_bins))
    return np.concatenate(coords), np.asarray(chrom_index, dtype=int), np.asarray(offsets, dtype=int)


def _candidate_pairs(chrom_index: np.ndarray, offsets: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    left, right, cis = [], [], []
    for i in range(chrom_index.size - 1):
        for j in range(i + 1, chrom_index.size):
            left.append(i)
            right.append(j)
            cis.append(chrom_index[i] == chrom_index[j])
    return np.asarray(left), np.asarray(right), np.asarray(cis, dtype=bool)


def _state_weights(
    coordinates: np.ndarray, pair_i: np.ndarray, pair_j: np.ndarray, family: str,
) -> np.ndarray:
    delta = coordinates[pair_i, :, None, :] - coordinates[pair_j, None, :, :]
    distance = np.linalg.norm(delta, axis=-1).reshape(-1, 4)
    if family == "well_specified":
        weights = np.power(distance + 0.1, -2.0)
    else:
        # Mismatched capture adds a threshold, saturation, and random-ligation floor.
        weights = 1.0 / (1.0 + np.exp(np.clip((distance - 1.25) / 0.22, -60, 60)))
        weights += 0.015
    return weights


def _sample_pair_contacts(
    condition: SyntheticCondition,
    coordinates: np.ndarray,
    chrom_index: np.ndarray,
    offsets: np.ndarray,
    count: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    candidate_i, candidate_j, candidate_cis = _candidate_pairs(chrom_index, offsets)
    weights4 = _state_weights(coordinates, candidate_i, candidate_j, condition.generator_family)
    pair_weight = weights4.sum(axis=1)
    if condition.count_overdispersion > 0:
        shape = 1.0 / condition.count_overdispersion
        pair_weight *= rng.gamma(shape=shape, scale=1.0 / shape, size=pair_weight.size)
    output_i = np.empty(count, dtype=int)
    output_j = np.empty(count, dtype=int)
    state = np.empty(count, dtype=int)
    background = rng.random(count) < condition.background_fraction
    cis_draw = rng.random(count) < condition.cis_fraction
    for index in range(count):
        eligible = candidate_cis == cis_draw[index]
        candidate_indices = np.flatnonzero(eligible)
        if background[index]:
            selected = int(rng.choice(candidate_indices))
            state[index] = int(rng.integers(0, 4))
        else:
            probabilities = pair_weight[candidate_indices]
            probabilities = probabilities / probabilities.sum()
            selected = int(rng.choice(candidate_indices, p=probabilities))
            state_prob = weights4[selected] / weights4[selected].sum()
            state[index] = int(rng.choice(4, p=state_prob))
        output_i[index], output_j[index] = candidate_i[selected], candidate_j[selected]
    return output_i, output_j, state, background


def _sample_molecule(
    coordinates: np.ndarray, length: int, rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    n_bins = coordinates.shape[0]
    segments = [int(rng.integers(n_bins))]
    copies = [int(rng.integers(2))]
    while len(segments) < length:
        origin = coordinates[segments[-1], copies[-1]]
        distance = np.linalg.norm(coordinates - origin, axis=-1)
        weights = np.exp(-distance / 0.8)
        weights[np.asarray(segments), :] = 0
        flat = weights.ravel()
        if flat.sum() <= 0:
            break
        selected = int(rng.choice(flat.size, p=flat / flat.sum()))
        segments.append(selected // 2)
        copies.append(selected % 2)
    return np.asarray(segments, dtype=int), np.asarray(copies, dtype=int)


def generate_synthetic(condition: SyntheticCondition) -> SyntheticDataset:
    condition.validate()
    sequence = np.random.SeedSequence(condition.seed)
    geometry_rng, contact_rng, molecule_rng, anchor_rng = [np.random.default_rng(item) for item in sequence.spawn(4)]
    coordinates, chrom_index, offsets = _geometry(condition, geometry_rng)
    per_molecule_pairs = condition.readchain_length * (condition.readchain_length - 1) // 2
    molecule_contact_budget = min(condition.contacts_per_cell, condition.molecule_count * per_molecule_pairs)
    ordinary_count = condition.contacts_per_cell - molecule_contact_budget
    pair_i, pair_j, state, background = _sample_pair_contacts(
        condition, coordinates, chrom_index, offsets, ordinary_count, contact_rng
    )
    molecule_ids = np.full(ordinary_count, -1, dtype=int)
    molecule_segments: list[np.ndarray] = []
    molecule_truth: list[np.ndarray] = []
    expanded_i: list[int] = []
    expanded_j: list[int] = []
    expanded_state: list[int] = []
    expanded_molecule: list[int] = []
    for _ in range(condition.molecule_count):
        segments, copies = _sample_molecule(coordinates, condition.readchain_length, molecule_rng)
        if segments.size < 2:
            continue
        molecule_id = len(molecule_segments)
        molecule_segments.append(segments)
        molecule_truth.append(copies)
        for left, right in itertools.combinations(range(segments.size), 2):
            i, j = int(segments[left]), int(segments[right])
            a, b = int(copies[left]), int(copies[right])
            if j < i:
                i, j, a, b = j, i, b, a
            expanded_i.append(i)
            expanded_j.append(j)
            expanded_state.append(2 * a + b)
            expanded_molecule.append(molecule_id)
            if ordinary_count + len(expanded_i) >= condition.contacts_per_cell:
                break
        if ordinary_count + len(expanded_i) >= condition.contacts_per_cell:
            break
    if expanded_i:
        pair_i = np.concatenate((pair_i, np.asarray(expanded_i, dtype=int)))
        pair_j = np.concatenate((pair_j, np.asarray(expanded_j, dtype=int)))
        state = np.concatenate((state, np.asarray(expanded_state, dtype=int)))
        background = np.concatenate((background, np.zeros(len(expanded_i), dtype=bool)))
        molecule_ids = np.concatenate((molecule_ids, np.asarray(expanded_molecule, dtype=int)))
    anchor_mask = np.zeros(pair_i.size, dtype=bool)
    eligible = np.flatnonzero(~background)
    requested_anchor_count = int(round(condition.anchor_fraction * pair_i.size))
    if condition.anchor_fraction > 0 and eligible.size:
        requested_anchor_count = max(1, requested_anchor_count)
    n_anchor = min(eligible.size, requested_anchor_count)
    if n_anchor:
        # Prefixes of one seeded permutation make anchor fractions strictly nested.
        anchor_order = anchor_rng.permutation(eligible)
        anchor_mask[anchor_order[:n_anchor]] = True
    return SyntheticDataset(
        condition=condition,
        truth_coordinates=coordinates,
        chrom_index=chrom_index,
        bin_offset=offsets,
        pair_i=pair_i,
        pair_j=pair_j,
        truth_state=state,
        is_background=background,
        anchor_mask=anchor_mask,
        molecule_id=molecule_ids,
        molecule_segments=tuple(molecule_segments),
        molecule_truth_copies=tuple(molecule_truth),
    )


def subset_synthetic(dataset: SyntheticDataset, indices: Sequence[int]) -> SyntheticDataset:
    """Select observations while retaining only complete, densely relabeled molecules."""
    selected = np.asarray(indices, dtype=int)
    if selected.ndim != 1 or (selected < 0).any() or (selected >= dataset.pair_i.size).any():
        raise ValueError("indices must be a valid one-dimensional observation index")
    old_molecule = dataset.molecule_id[selected]
    retained = sorted(set(int(item) for item in old_molecule if item >= 0))
    mapping = {old: new for new, old in enumerate(retained)}
    relabeled = np.asarray([mapping.get(int(item), -1) for item in old_molecule], dtype=int)
    for old in retained:
        all_members = np.flatnonzero(dataset.molecule_id == old)
        if not np.all(np.isin(all_members, selected)):
            raise ValueError("a source molecule was divided across the requested subset")
    condition = replace(
        dataset.condition,
        contacts_per_cell=int(selected.size),
        molecule_count=len(retained),
    )
    return SyntheticDataset(
        condition=condition,
        truth_coordinates=dataset.truth_coordinates.copy(),
        chrom_index=dataset.chrom_index.copy(),
        bin_offset=dataset.bin_offset.copy(),
        pair_i=dataset.pair_i[selected].copy(),
        pair_j=dataset.pair_j[selected].copy(),
        truth_state=dataset.truth_state[selected].copy(),
        is_background=dataset.is_background[selected].copy(),
        anchor_mask=dataset.anchor_mask[selected].copy(),
        molecule_id=relabeled,
        molecule_segments=tuple(dataset.molecule_segments[old].copy() for old in retained),
        molecule_truth_copies=tuple(dataset.molecule_truth_copies[old].copy() for old in retained),
    )


def observation_masses(dataset: SyntheticDataset) -> np.ndarray:
    """Assign one mass unit per pair or per source molecule."""
    masses = np.ones(dataset.pair_i.size, dtype=float)
    for molecule_id in np.unique(dataset.molecule_id[dataset.molecule_id >= 0]):
        selected = dataset.molecule_id == molecule_id
        masses[selected] = 1.0 / selected.sum()
    return masses


def regularization_scale(reference_mass: float, effective_mass: float) -> float:
    """Scale a fixed structural prior relative to accumulated observation mass."""
    if not np.isfinite(reference_mass) or reference_mass <= 0:
        raise ValueError("reference_mass must be finite and positive")
    if not np.isfinite(effective_mass) or effective_mass <= 0:
        raise ValueError("effective_mass must be finite and positive")
    return float(reference_mass / effective_mass)


def _sample_negatives(dataset: SyntheticDataset, ratio: float, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    n = max(1, int(round(dataset.pair_i.size * ratio)))
    cis_fraction = float(np.mean(dataset.chrom_index[dataset.pair_i] == dataset.chrom_index[dataset.pair_j]))
    candidate_i, candidate_j, candidate_cis = _candidate_pairs(
        dataset.chrom_index, dataset.bin_offset
    )
    cis_candidates = np.flatnonzero(candidate_cis)
    trans_candidates = np.flatnonzero(~candidate_cis)
    target_cis = rng.random(n) < cis_fraction
    selected = np.empty(n, dtype=int)
    if target_cis.any():
        selected[target_cis] = rng.choice(cis_candidates, size=int(target_cis.sum()), replace=True)
    if (~target_cis).any():
        selected[~target_cis] = rng.choice(trans_candidates, size=int((~target_cis).sum()), replace=True)
    # These are exposure/noise samples, not a claim that a bin pair was never
    # observed. At high contact depth every finite candidate pair may occur.
    return candidate_i[selected], candidate_j[selected]


def infer_synthetic(
    dataset: SyntheticDataset,
    config: InferenceConfig,
    *,
    initial_coordinates: np.ndarray | None = None,
) -> InferenceResult:
    """Fit only raw endpoint contacts plus explicitly configured phased anchors."""
    import torch
    import torch.nn.functional as functional

    torch.manual_seed(config.seed)
    rng = np.random.default_rng(config.seed)
    device = torch.device(config.device)
    dtype = torch.float64
    n_bins = dataset.n_bins
    observation_model = config.observation_model
    if observation_model == "auto":
        observation_model = (
            "normalized_endpoint"
            if dataset.condition.generator_family == "well_specified"
            and dataset.condition.count_overdispersion == 0
            and not dataset.molecule_segments
            else "contrastive"
        )
    if observation_model not in {"normalized_endpoint", "contrastive"}:
        raise ValueError("observation_model must be auto, normalized_endpoint, or contrastive")
    if observation_model == "normalized_endpoint" and config.molecule_likelihood:
        raise ValueError(
            "normalized_endpoint is a pairwise likelihood; exact shared-copy molecule mode uses contrastive normalization"
        )
    if not 0 <= config.background_fraction_min < config.background_fraction_max < 1:
        raise ValueError("background fraction bounds must satisfy 0 <= min < max < 1")
    if config.regularization_reference_mass <= 0:
        raise ValueError("regularization_reference_mass must be positive")
    if config.molecule_likelihood and np.any(dataset.anchor_mask & (dataset.molecule_id >= 0)):
        raise ValueError(
            "phased molecule anchors require a joint constrained-assignment likelihood"
        )
    if initial_coordinates is None:
        midpoint_init = np.zeros((n_bins, 3), dtype=float)
        for chrom in range(dataset.condition.chromosome_count):
            selected = np.flatnonzero(dataset.chrom_index == chrom)
            midpoint_init[selected] = _polymer(rng, selected.size, rng.normal(scale=2.0, size=3))
        delta_init = rng.normal(scale=0.05, size=(n_bins, 3))
    else:
        initial = np.asarray(initial_coordinates, dtype=float)
        if initial.shape != (n_bins, 2, 3) or not np.isfinite(initial).all():
            raise ValueError("initial_coordinates must be a finite bins x 2 x 3 array")
        midpoint_init = initial.mean(axis=1)
        delta_init = 0.5 * (initial[:, 1] - initial[:, 0])
    midpoint = torch.nn.Parameter(torch.tensor(midpoint_init, dtype=dtype, device=device))
    delta = torch.nn.Parameter(torch.tensor(delta_init, dtype=dtype, device=device))
    intercept = torch.nn.Parameter(torch.tensor(0.0, dtype=dtype, device=device))
    optimizer = torch.optim.Adam((midpoint, delta, intercept), lr=config.learning_rate)
    pair_i = torch.tensor(dataset.pair_i, dtype=torch.long, device=device)
    pair_j = torch.tensor(dataset.pair_j, dtype=torch.long, device=device)
    truth_state = torch.tensor(dataset.truth_state, dtype=torch.long, device=device)
    anchor_mask = torch.tensor(dataset.anchor_mask, dtype=torch.bool, device=device)
    chrom = torch.tensor(dataset.chrom_index, dtype=torch.long, device=device)
    positive_weight = torch.tensor(observation_masses(dataset), dtype=dtype, device=device)
    pair_is_cis = chrom[pair_i] == chrom[pair_j]
    candidate_i_np, candidate_j_np, candidate_cis_np = _candidate_pairs(
        dataset.chrom_index, dataset.bin_offset
    )
    candidate_lookup = {
        (int(left), int(right)): index
        for index, (left, right) in enumerate(zip(candidate_i_np, candidate_j_np, strict=True))
    }
    observed_candidate_index = torch.tensor(
        [candidate_lookup[(int(left), int(right))] for left, right in zip(dataset.pair_i, dataset.pair_j, strict=True)],
        dtype=torch.long,
        device=device,
    )
    candidate_i = torch.tensor(candidate_i_np, dtype=torch.long, device=device)
    candidate_j = torch.tensor(candidate_j_np, dtype=torch.long, device=device)
    candidate_is_cis = torch.tensor(candidate_cis_np, dtype=torch.bool, device=device)
    if observation_model == "contrastive":
        neg_i_np, neg_j_np = _sample_negatives(dataset, config.negative_ratio, rng)
        neg_i = torch.tensor(neg_i_np, dtype=torch.long, device=device)
        neg_j = torch.tensor(neg_j_np, dtype=torch.long, device=device)
    else:
        neg_i = neg_j = None
    molecule_cache: list[tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, bool]] = []
    for segments_np in dataset.molecule_segments:
        segments = torch.tensor(segments_np, dtype=torch.long, device=device)
        n_segments = int(segments.numel())
        assignments = torch.tensor(
            list(itertools.product((0, 1), repeat=n_segments)), dtype=torch.long, device=device
        )
        pair_indices = list(itertools.combinations(range(n_segments), 2))
        pair_left = torch.tensor([left for left, _ in pair_indices], dtype=torch.long, device=device)
        pair_right = torch.tensor([right for _, right in pair_indices], dtype=torch.long, device=device)
        molecule_is_cis = bool(np.unique(dataset.chrom_index[segments_np]).size == 1)
        molecule_cache.append((segments, assignments, pair_left, pair_right, molecule_is_cis))

    def coordinates() -> torch.Tensor:
        return torch.stack((midpoint - delta, midpoint + delta), dim=1)

    def log_kernel(distances: torch.Tensor) -> torch.Tensor:
        if config.likelihood_mode == "monotone_log":
            return -config.alpha * torch.log(distances + config.d0)
        if config.likelihood_mode == "monotonic_powerlaw":
            return -torch.log1p(distances.square())
        if config.likelihood_mode == "monotone_logistic":
            return -functional.softplus((distances - 1.25) / 0.25)
        if config.likelihood_mode == "dist2":
            return -distances.square()
        if config.likelihood_mode == "logdist2":
            return -torch.log1p(distances.square() / 1e-6)
        if config.likelihood_mode == "fdg_flat":
            close = torch.square(torch.clamp(0.5 - distances, min=0.0))
            shoulder = torch.square(torch.clamp(distances - 1.5, min=0.0, max=0.5))
            far_r = torch.clamp(distances, min=2.0)
            far = 1.5 * (far_r - 2.0) + 0.125 / (far_r - 1.5)
            energy = torch.where(
                distances < 0.5, close,
                torch.where(distances <= 1.5, torch.zeros_like(distances),
                            torch.where(distances <= 2.0, shoulder, far)),
            )
            return -energy
        raise ValueError(f"unsupported synthetic inference likelihood: {config.likelihood_mode}")

    def log_components(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
        coord = coordinates()
        distances = torch.sqrt(
            (coord[left, :, None, :] - coord[right, None, :, :]).square().sum(dim=-1) + 1.0e-18
        )
        return log_kernel(distances)

    def normalized_endpoint_log_probabilities(
        candidate_components: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Return log P(endpoint pair | cis/trans stratum) and signal posterior."""
        candidate_log_weight = torch.logsumexp(candidate_components, dim=1) - math.log(4.0)
        log_probability = candidate_log_weight.new_empty(pair_i.numel())
        signal_probability = candidate_log_weight.new_ones(pair_i.numel())
        for stratum in (True, False):
            observed = pair_is_cis == stratum
            candidates = candidate_is_cis == stratum
            if not bool(observed.any()):
                continue
            if not bool(candidates.any()):
                raise ValueError("an observed cis/trans stratum has no candidate endpoint pairs")
            log_partition = torch.logsumexp(candidate_log_weight[candidates], dim=0)
            signal_log_probability = (
                candidate_log_weight[observed_candidate_index[observed]] - log_partition
            )
            if config.background_component:
                background_log_probability = -math.log(int(candidates.sum().item()))
                background_fraction = config.background_fraction_min + (
                    config.background_fraction_max - config.background_fraction_min
                ) * torch.sigmoid(-intercept)
                log_signal_fraction = torch.log1p(-background_fraction)
                log_background_fraction = torch.log(background_fraction)
                mixed = torch.logaddexp(
                    log_signal_fraction + signal_log_probability,
                    log_background_fraction + background_log_probability,
                )
                log_probability[observed] = mixed
                signal_probability[observed] = torch.exp(
                    log_signal_fraction + signal_log_probability - mixed
                )
            else:
                log_probability[observed] = signal_log_probability
        return log_probability, signal_probability

    history: list[float] = []
    staged_reference: torch.Tensor | None = None
    started = time.monotonic()
    for iteration in range(config.iterations):
        optimizer.zero_grad(set_to_none=True)
        components = log_components(pair_i, pair_j).reshape(-1, 4)
        if config.staged_cis_then_trans and iteration < config.iterations // 2:
            selected = pair_is_cis
        else:
            selected = torch.ones(pair_i.numel(), dtype=torch.bool, device=device)
        if observation_model == "normalized_endpoint":
            candidate_components = log_components(candidate_i, candidate_j).reshape(-1, 4)
            endpoint_log_probability, signal_probability = normalized_endpoint_log_probabilities(
                candidate_components
            )
            observation_log_probability = endpoint_log_probability.clone()
            selected_anchors = selected & anchor_mask
            if bool(selected_anchors.any()):
                state_log_probability = functional.log_softmax(components, dim=1)
                anchor_index = torch.nonzero(selected_anchors, as_tuple=False).flatten()
                observation_log_probability[anchor_index] = (
                    endpoint_log_probability[anchor_index]
                    + torch.log(signal_probability[anchor_index].clamp_min(1.0e-15))
                    + state_log_probability[anchor_index, truth_state[anchor_index]]
                )
            selected_mass = positive_weight[selected].sum().clamp_min(1.0)
            loss = -(
                observation_log_probability[selected] * positive_weight[selected]
            ).sum() / selected_mass
        else:
            assert neg_i is not None and neg_j is not None
            signal_score = torch.logsumexp(components, dim=1) - math.log(4.0) + intercept
            negative_components = log_components(neg_i, neg_j).reshape(-1, 4)
            negative_score = torch.logsumexp(negative_components, dim=1) - math.log(4.0) + intercept
            negative_selected = (
                chrom[neg_i] == chrom[neg_j]
                if config.staged_cis_then_trans and iteration < config.iterations // 2
                else torch.ones(neg_i.numel(), dtype=torch.bool, device=device)
            )
            if config.molecule_likelihood:
                ordinary = dataset.molecule_id < 0
                ordinary_selected = selected & torch.tensor(ordinary, dtype=torch.bool, device=device)
            else:
                ordinary_selected = selected
            if config.background_component:
                positive_terms = functional.softplus(-signal_score[ordinary_selected])
            else:
                positive_terms = -(
                    torch.logsumexp(components[ordinary_selected], dim=1) - math.log(4.0)
                )
            ordinary_anchor = anchor_mask[ordinary_selected]
            if bool(ordinary_anchor.any()):
                ordinary_components = components[ordinary_selected]
                ordinary_truth = truth_state[ordinary_selected]
                phase_nll = -functional.log_softmax(ordinary_components, dim=1)[
                    torch.arange(ordinary_components.shape[0], device=device), ordinary_truth
                ]
                positive_terms = positive_terms + torch.where(
                    ordinary_anchor, phase_nll, torch.zeros_like(phase_nll)
                )
            positive_loss_sum = (positive_terms * positive_weight[ordinary_selected]).sum()
            positive_mass = positive_weight[ordinary_selected].sum()
            if config.molecule_likelihood and molecule_cache:
                molecule_losses: list[torch.Tensor] = []
                coord = coordinates()
                for segments, assignments, pair_left, pair_right, molecule_is_cis in molecule_cache:
                    if config.staged_cis_then_trans and iteration < config.iterations // 2 and not molecule_is_cis:
                        continue
                    n_segments = int(segments.numel())
                    molecule_coordinates = coord[segments]
                    pair_distance = torch.sqrt(
                        (
                            molecule_coordinates[pair_left, :, None, :]
                            - molecule_coordinates[pair_right, None, :, :]
                        ).square().sum(dim=-1)
                        + 1.0e-18
                    )
                    pair_log_kernel = log_kernel(pair_distance)
                    assignment_pair_score = pair_log_kernel[
                        torch.arange(pair_left.numel(), device=device)[None, :],
                        assignments[:, pair_left],
                        assignments[:, pair_right],
                    ]
                    assignment_scores = assignment_pair_score.mean(dim=1)
                    molecule_score = torch.logsumexp(assignment_scores, dim=0)
                    molecule_score = molecule_score - n_segments * math.log(2.0)
                    if config.background_component:
                        molecule_score = molecule_score + intercept
                        molecule_losses.append(functional.softplus(-molecule_score))
                    else:
                        molecule_losses.append(-molecule_score)
                if molecule_losses:
                    positive_loss_sum = positive_loss_sum + torch.stack(molecule_losses).sum()
                    positive_mass = positive_mass + len(molecule_losses)
            positive_loss = positive_loss_sum / positive_mass.clamp_min(1.0)
            negative_loss = (
                functional.softplus(negative_score[negative_selected]).mean()
                if config.background_component else signal_score.new_zeros(())
            )
            loss = positive_loss + negative_loss
            selected_mass = positive_mass.clamp_min(1.0)
        backbone_terms = []
        bending_terms = []
        delta_smooth_terms = []
        for chrom_id in torch.unique(chrom):
            indices = torch.nonzero(chrom == chrom_id, as_tuple=False).flatten()
            step = midpoint[indices[1:]] - midpoint[indices[:-1]]
            backbone_terms.append((step.norm(dim=1) - 1.0).square().mean())
            if indices.numel() >= 3:
                bending_terms.append(
                    (midpoint[indices[2:]] - 2 * midpoint[indices[1:-1]] + midpoint[indices[:-2]]).square().sum(dim=1).mean()
                )
            delta_smooth_terms.append((delta[indices[1:]] - delta[indices[:-1]]).square().sum(dim=1).mean())
        centered = midpoint - midpoint.mean(dim=0, keepdim=True)
        rms = centered.square().sum(dim=1).mean().sqrt()
        regularizer = config.backbone_weight * torch.stack(backbone_terms).mean()
        regularizer = regularizer + config.bending_weight * torch.stack(bending_terms).mean()
        regularizer = regularizer + config.delta_smooth_weight * torch.stack(delta_smooth_terms).mean()
        regularizer = regularizer + config.delta_shrink_weight * delta.square().sum(dim=1).mean()
        if config.minimum_homolog_separation_weight > 0:
            homolog_distance = 2.0 * delta.norm(dim=1)
            regularizer = regularizer + config.minimum_homolog_separation_weight * torch.relu(
                config.minimum_homolog_separation - homolog_distance
            ).square().mean()
        if config.global_copy_track_weight > 0 and int(torch.unique(chrom).numel()) > 1:
            directions = torch.stack([
                delta[chrom == chrom_id].mean(dim=0) for chrom_id in torch.unique(chrom)
            ])
            directions = directions / directions.norm(dim=1, keepdim=True).clamp_min(1.0e-12)
            agreement = directions @ directions.T
            upper = torch.triu(torch.ones_like(agreement, dtype=torch.bool), diagonal=1)
            regularizer = regularizer + config.global_copy_track_weight * (
                1.0 - agreement[upper]
            ).square().mean()
        regularizer = regularizer + config.scale_weight * (rms - 3.0).square()
        if config.staged_cis_then_trans and iteration == config.iterations // 2:
            staged_reference = coordinates().detach().clone()
        if staged_reference is not None:
            current = coordinates()
            shape_terms = []
            for chrom_id in torch.unique(chrom):
                indices = torch.nonzero(chrom == chrom_id, as_tuple=False).flatten()
                for copy in (0, 1):
                    ref_distance = torch.pdist(staged_reference[indices, copy])
                    current_distance = torch.pdist(current[indices, copy])
                    shape_terms.append((current_distance - ref_distance).square().mean())
            refinement_start = 3 * config.iterations // 4
            shape_weight = (
                config.stage_shape_weight if iteration < refinement_start
                else config.joint_refine_shape_weight
            )
            regularizer = regularizer + shape_weight * torch.stack(shape_terms).mean()
        prior_scale = regularization_scale(
            config.regularization_reference_mass,
            float(selected_mass.detach().cpu()),
        )
        loss = loss + prior_scale * regularizer
        if not torch.isfinite(loss):
            break
        loss.backward()
        torch.nn.utils.clip_grad_norm_((midpoint, delta, intercept), max_norm=10.0)
        optimizer.step()
        history.append(float(loss.detach().cpu()))
    runtime = time.monotonic() - started
    with torch.no_grad():
        final_coordinates = coordinates().cpu().numpy()
        components = log_components(pair_i, pair_j).reshape(-1, 4)
        q4_conditional = torch.softmax(components, dim=1)
        if config.molecule_likelihood and molecule_cache:
            coord = coordinates()
            q4_conditional = q4_conditional.clone()
            for molecule_id, (segments, assignments, pair_left, pair_right, _) in enumerate(
                molecule_cache
            ):
                molecule_coordinates = coord[segments]
                pair_distance = torch.sqrt(
                    (
                        molecule_coordinates[pair_left, :, None, :]
                        - molecule_coordinates[pair_right, None, :, :]
                    ).square().sum(dim=-1)
                    + 1.0e-18
                )
                pair_log_kernel = log_kernel(pair_distance)
                assignment_pair_score = pair_log_kernel[
                    torch.arange(pair_left.numel(), device=device)[None, :],
                    assignments[:, pair_left],
                    assignments[:, pair_right],
                ]
                assignment_posterior = torch.softmax(
                    assignment_pair_score.mean(dim=1), dim=0
                )
                segment_lookup = {
                    int(bin_index): position
                    for position, bin_index in enumerate(segments.detach().cpu().numpy())
                }
                for observation_index in np.flatnonzero(dataset.molecule_id == molecule_id):
                    left_position = segment_lookup[int(dataset.pair_i[observation_index])]
                    right_position = segment_lookup[int(dataset.pair_j[observation_index])]
                    states = (
                        2 * assignments[:, left_position] + assignments[:, right_position]
                    )
                    marginal = torch.stack([
                        assignment_posterior[states == state].sum()
                        for state in range(4)
                    ])
                    q4_conditional[observation_index] = marginal
        if observation_model == "normalized_endpoint":
            candidate_components = log_components(candidate_i, candidate_j).reshape(-1, 4)
            _, signal_probability = normalized_endpoint_log_probabilities(candidate_components)
        else:
            signal_score = torch.logsumexp(components, dim=1) - math.log(4.0) + intercept
            signal_probability = torch.sigmoid(signal_score)
        if config.background_component:
            q4 = (q4_conditional * signal_probability[:, None]).cpu().numpy()
            qbg = (1.0 - signal_probability).cpu().numpy()
            if config.molecule_likelihood:
                molecule_observations = dataset.molecule_id >= 0
                q4[molecule_observations] = q4_conditional.cpu().numpy()[molecule_observations]
                qbg[molecule_observations] = 0.0
            if observation_model == "normalized_endpoint":
                fitted_background_fraction = config.background_fraction_min + (
                    config.background_fraction_max - config.background_fraction_min
                ) * torch.sigmoid(-intercept)
                fitted_signal_fraction = float((1.0 - fitted_background_fraction).cpu())
            else:
                fitted_signal_fraction = float(torch.sigmoid(intercept).cpu())
        else:
            q4 = q4_conditional.cpu().numpy()
            qbg = np.zeros(dataset.pair_i.size, dtype=float)
            fitted_signal_fraction = 1.0
    anchor_flips = infer_anchor_gauge(dataset, q4)
    if anchor_flips:
        q4 = np.vstack([
            relabel_probabilities(
                row,
                flip_left=anchor_flips.get(int(dataset.chrom_index[dataset.pair_i[index]]), 0),
                flip_right=anchor_flips.get(int(dataset.chrom_index[dataset.pair_j[index]]), 0),
            )
            for index, row in enumerate(q4)
        ])
        for bin_index, chrom_id in enumerate(dataset.chrom_index):
            if anchor_flips.get(int(chrom_id), 0):
                final_coordinates[bin_index] = final_coordinates[bin_index, ::-1]
    converged = bool(history) and np.isfinite(final_coordinates).all() and (
        len(history) < 20 or abs(history[-1] - history[-20]) / max(1.0, abs(history[-20])) < 0.02
    )
    return InferenceResult(
        coordinates=final_coordinates,
        state_probabilities=q4,
        background_probability=qbg,
        objective_history=history,
        runtime_seconds=runtime,
        converged=converged,
        anchor_flips=anchor_flips,
        signal_fraction=fitted_signal_fraction,
    )


def infer_anchor_gauge(dataset: SyntheticDataset, probabilities: np.ndarray) -> dict[int, int]:
    """Select chromosome flips using only explicitly phased anchor observations."""
    anchors = np.flatnonzero(dataset.anchor_mask)
    if anchors.size == 0:
        return {}
    chromosomes = sorted(set(dataset.chrom_index[dataset.pair_i[anchors]]) | set(dataset.chrom_index[dataset.pair_j[anchors]]))
    best_score = -math.inf
    best: dict[int, int] = {}
    for bits in itertools.product((0, 1), repeat=len(chromosomes)):
        flips = dict(zip(chromosomes, bits, strict=True))
        score = 0.0
        for index in anchors:
            c1 = int(dataset.chrom_index[dataset.pair_i[index]])
            c2 = int(dataset.chrom_index[dataset.pair_j[index]])
            state = flip_state(int(dataset.truth_state[index]), flips.get(c1, 0), flips.get(c2, 0))
            score += math.log(max(float(probabilities[index, state]), 1e-15))
        if score > best_score:
            best_score, best = score, flips
    return best


def geometry_selected_chromosome_flips(
    dataset: SyntheticDataset, coordinates: np.ndarray,
) -> dict[int, int]:
    """Select per-chromosome copy flips from reference geometry for evaluation only."""
    values = np.asarray(coordinates, dtype=float)
    if values.shape != dataset.truth_coordinates.shape or not np.isfinite(values).all():
        raise ValueError("coordinates must match finite synthetic truth coordinates")
    flips: dict[int, int] = {}
    for chrom in range(dataset.condition.chromosome_count):
        selected = dataset.chrom_index == chrom
        if not selected.any():
            raise ValueError(f"synthetic chromosome {chrom} has no bins")
        _, flips[chrom] = chromosome_swap_invariant_metrics(
            dataset.truth_coordinates[selected], values[selected]
        )
    return flips


def relabel_contact_probabilities_to_geometry_gauge(
    dataset: SyntheticDataset,
    probabilities: np.ndarray,
    flips: dict[int, int],
) -> np.ndarray:
    """Relabel contact probabilities using evaluation-only geometry-selected flips."""
    values = np.asarray(probabilities, dtype=float)
    if values.shape != (dataset.pair_i.size, 4):
        raise ValueError("probabilities must contain one four-state row per observation")
    expected = set(range(dataset.condition.chromosome_count))
    if set(flips) != expected or any(value not in (0, 1) for value in flips.values()):
        raise ValueError("flips must define one binary value for every chromosome")
    return np.vstack([
        relabel_probabilities(
            row,
            flip_left=flips[int(dataset.chrom_index[dataset.pair_i[index]])],
            flip_right=flips[int(dataset.chrom_index[dataset.pair_j[index]])],
        )
        for index, row in enumerate(values)
    ])


def evaluate_synthetic(dataset: SyntheticDataset, result: InferenceResult) -> dict[str, float | int | str]:
    evaluation = ~dataset.anchor_mask
    signal = ~dataset.is_background & evaluation
    if not signal.any():
        raise ValueError("synthetic evaluation has no non-anchor signal observations")
    raw_metrics, _, raw_selective = contact_metrics(
        dataset.truth_state[signal], result.state_probabilities[signal]
    )
    cis = dataset.chrom_index[dataset.pair_i] == dataset.chrom_index[dataset.pair_j]
    raw_scoped_metrics: dict[str, float | int] = {}
    for scope, selected in (("cis", signal & cis), ("trans", signal & ~cis)):
        if selected.any():
            scoped, _, _ = contact_metrics(dataset.truth_state[selected], result.state_probabilities[selected])
            for name in (
                "n_contacts", "exact_top1_accuracy", "exact_log_loss", "exact_brier",
                "same_cross_accuracy", "same_cross_log_loss", "same_cross_brier",
                "exact_ece", "same_cross_ece",
            ):
                raw_scoped_metrics[f"{scope}_{name}"] = scoped[name]
    per_chrom = []
    swaps = []
    for chrom in range(dataset.condition.chromosome_count):
        selected = dataset.chrom_index == chrom
        row, swap = chromosome_swap_invariant_metrics(
            dataset.truth_coordinates[selected], result.coordinates[selected]
        )
        row["homolog_separation_spearman"] = homolog_separation_correlation(
            dataset.truth_coordinates[selected], result.coordinates[selected]
        )
        per_chrom.append(row)
        swaps.append(swap)
    geometry_flips = {chrom: int(swap) for chrom, swap in enumerate(swaps)}
    aligned_probabilities = relabel_contact_probabilities_to_geometry_gauge(
        dataset, result.state_probabilities, geometry_flips
    )
    aligned_metrics, _, aligned_selective = contact_metrics(
        dataset.truth_state[signal], aligned_probabilities[signal]
    )
    aligned_scoped_metrics: dict[str, float | int] = {}
    for scope, selected in (("cis", signal & cis), ("trans", signal & ~cis)):
        if selected.any():
            scoped, _, _ = contact_metrics(dataset.truth_state[selected], aligned_probabilities[selected])
            for name in (
                "n_contacts", "exact_top1_accuracy", "exact_log_loss", "exact_brier",
                "same_cross_accuracy", "same_cross_log_loss", "same_cross_brier",
                "exact_ece", "same_cross_ece",
            ):
                aligned_scoped_metrics[f"geometry_gauge_aligned_{scope}_{name}"] = scoped[name]
    raw_exact = float(raw_metrics["exact_top1_accuracy"])
    raw_same_cross = float(raw_metrics["same_cross_accuracy"])
    output: dict[str, float | int | str] = {
        **raw_metrics,
        **{
            f"geometry_gauge_aligned_{name}": value
            for name, value in aligned_metrics.items()
        },
        "generator_family": dataset.condition.generator_family,
        "n_chromosomes": dataset.condition.chromosome_count,
        "n_cells": 1,
        "n_background_contacts": int(dataset.is_background.sum()),
        "n_anchor_contacts": int(dataset.anchor_mask.sum()),
        "n_evaluation_contacts": int(evaluation.sum()),
        "actual_anchor_fraction": float(dataset.anchor_mask.mean()),
        "effective_observation_mass": float(observation_masses(dataset).sum()),
        "effective_anchor_mass": float(observation_masses(dataset)[dataset.anchor_mask].sum()),
        "unordered_structure_rmsd": float(np.mean([float(row["procrustes_rmsd"]) for row in per_chrom])),
        "unordered_structure_distance_spearman": float(np.nanmean([float(row["distance_spearman"]) for row in per_chrom])),
        "globally_aligned_structure_rmsd": float(
            chromosome_swap_invariant_metrics(dataset.truth_coordinates, result.coordinates)[0]["procrustes_rmsd"]
        ),
        "gauge_recovery": float(max(np.mean(np.asarray(swaps) == 0), np.mean(np.asarray(swaps) == 1))),
        "relative_gauge_recovery": float(np.mean(
            (np.asarray(swaps[1:]) ^ swaps[0]) == 0
        )) if len(swaps) > 1 else 1.0,
        "mean_homolog_separation_spearman": float(np.nanmean([
            float(row["homolog_separation_spearman"]) for row in per_chrom
        ])),
        "runtime_seconds": result.runtime_seconds,
        "fitted_signal_fraction": result.signal_fraction,
        "converged": int(result.converged),
        "failure": int(not result.converged or not np.isfinite(result.coordinates).all()),
        "raw_exact_accuracy": raw_exact,
        "raw_same_cross_accuracy": raw_same_cross,
        "raw_contact_gauge_policy": (
            "EXTERNAL_ANCHOR_ALIGNED_INFERENCE_GAUGE"
            if dataset.anchor_mask.any() else "ARBITRARY_INFERENCE_GAUGE"
        ),
        "geometry_gauge_alignment_policy": "ORACLE_GEOMETRY_ALIGNMENT_EVAL_ONLY",
        "geometry_gauge_flip_pattern": ",".join(
            f"chr{chrom + 1}:{flip}" for chrom, flip in sorted(geometry_flips.items())
        ),
        **raw_scoped_metrics,
        **aligned_scoped_metrics,
    }
    background_truth = dataset.is_background[evaluation].astype(int)
    background_probability = np.asarray(result.background_probability, dtype=float)[evaluation]
    if background_probability.shape == background_truth.shape:
        background_probability = np.clip(background_probability, 1.0e-15, 1.0 - 1.0e-15)
        background_prediction = background_probability >= 0.5
        background_correct = background_prediction == background_truth
        background_confidence = np.where(
            background_prediction, background_probability, 1.0 - background_probability
        )
        background_truth_probability = np.where(
            background_truth == 1, background_probability, 1.0 - background_probability
        )
        background_curve = reliability_curve(
            background_correct, background_confidence, metric="background"
        )
        output.update({
            "background_accuracy": float(background_correct.mean()),
            "background_log_loss": float(-np.log(background_truth_probability).mean()),
            "background_brier": float(np.square(background_probability - background_truth).mean()),
            "background_ece": expected_calibration_error(background_curve),
        })
    for row in raw_selective:
        output[
            f"{row['metric']}_coverage_at_precision_{float(row['target_precision']):.2f}"
        ] = float(row["callable_fraction"])
    for row in aligned_selective:
        output[
            "geometry_gauge_aligned_"
            f"{row['metric']}_coverage_at_precision_{float(row['target_precision']):.2f}"
        ] = float(row["callable_fraction"])
    return output


def oracle_coordinate_posteriors(dataset: SyntheticDataset, mode: str = "monotone_log") -> np.ndarray:
    config = ScoreConfig(mode=mode, alpha=2.0, d0=0.1)
    rows = []
    for left, right in zip(dataset.pair_i, dataset.pair_j, strict=True):
        distance = np.linalg.norm(
            dataset.truth_coordinates[left, :, None, :] - dataset.truth_coordinates[right, None, :, :], axis=-1
        ).reshape(4)
        q4, _, _ = posterior_from_distances(distance, config)
        rows.append(q4)
    return np.asarray(rows)


def oracle_generator_conditional_posteriors(dataset: SyntheticDataset) -> np.ndarray:
    """Return the true pairwise state conditional for non-molecule synthetic data.

    This explicitly ORACLE diagnostic uses both reference coordinates and the
    known generator family. It is not available to blind inference.
    """
    if np.any(dataset.molecule_id >= 0):
        raise ValueError("generator-conditional pair posterior is undefined for molecule-derived contacts")
    weights = _state_weights(
        dataset.truth_coordinates, dataset.pair_i, dataset.pair_j,
        dataset.condition.generator_family,
    )
    totals = weights.sum(axis=1, keepdims=True)
    if not np.isfinite(weights).all() or (weights < 0).any() or (totals <= 0).any():
        raise AssertionError("synthetic generator produced invalid state weights")
    return weights / totals


def score_coordinates(
    dataset: SyntheticDataset,
    coordinates: np.ndarray,
    *,
    likelihood_mode: str = "monotone_log",
    alpha: float = 2.0,
    d0: float = 0.1,
    observation_model: str = "contrastive",
    signal_fraction: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Score observations without consulting phase labels or reference coordinates.

    In ``normalized_endpoint`` mode, returns the exact cis/trans-stratified
    endpoint-mixture log likelihood, including an optional uniform background.
    ``contrastive`` mode returns the local four-state kernel score and is labeled
    pseudo-likelihood by callers.
    """
    values = np.asarray(coordinates, dtype=float)
    if values.shape != (dataset.n_bins, 2, 3) or not np.isfinite(values).all():
        raise ValueError("coordinates must be a finite bins x 2 x 3 array")
    config = ScoreConfig(mode=likelihood_mode, alpha=alpha, d0=d0)
    delta = values[dataset.pair_i, :, None, :] - values[dataset.pair_j, None, :, :]
    distances = np.linalg.norm(delta, axis=-1).reshape(-1, 4)
    scores = log_score(distances.ravel(), config).reshape(-1, 4)
    maximum = scores.max(axis=1, keepdims=True)
    weights = np.exp(np.clip(scores - maximum, -745, 0))
    posterior = weights / weights.sum(axis=1, keepdims=True)
    log_marginal = maximum[:, 0] + np.log(weights.sum(axis=1)) - math.log(4.0)
    if observation_model == "normalized_endpoint":
        if not 0 < signal_fraction <= 1:
            raise ValueError("signal_fraction must be in (0,1]")
        candidate_i, candidate_j, candidate_cis = _candidate_pairs(
            dataset.chrom_index, dataset.bin_offset
        )
        candidate_delta = (
            values[candidate_i, :, None, :] - values[candidate_j, None, :, :]
        )
        candidate_distance = np.linalg.norm(candidate_delta, axis=-1).reshape(-1, 4)
        candidate_score = log_score(candidate_distance.ravel(), config).reshape(-1, 4)
        candidate_maximum = candidate_score.max(axis=1, keepdims=True)
        candidate_log_weight = (
            candidate_maximum[:, 0]
            + np.log(np.exp(np.clip(candidate_score - candidate_maximum, -745, 0)).sum(axis=1))
            - math.log(4.0)
        )
        lookup = {
            (int(left), int(right)): index
            for index, (left, right) in enumerate(zip(candidate_i, candidate_j, strict=True))
        }
        observed_index = np.asarray([
            lookup[(int(left), int(right))]
            for left, right in zip(dataset.pair_i, dataset.pair_j, strict=True)
        ])
        observed_cis = dataset.chrom_index[dataset.pair_i] == dataset.chrom_index[dataset.pair_j]
        for stratum in (True, False):
            observed = observed_cis == stratum
            candidates = candidate_cis == stratum
            if not observed.any():
                continue
            candidate_values = candidate_log_weight[candidates]
            partition_max = float(candidate_values.max())
            log_partition = partition_max + math.log(
                float(np.exp(candidate_values - partition_max).sum())
            )
            signal_log_probability = candidate_log_weight[observed_index[observed]] - log_partition
            if signal_fraction < 1.0:
                background_log_probability = -math.log(int(candidates.sum()))
                signal_term = math.log(signal_fraction) + signal_log_probability
                background_term = math.log1p(-signal_fraction) + background_log_probability
                mixed = np.logaddexp(signal_term, background_term)
                posterior[observed] *= np.exp(signal_term - mixed)[:, None]
                log_marginal[observed] = mixed
            else:
                log_marginal[observed] = signal_log_probability
    elif observation_model != "contrastive":
        raise ValueError("observation_model must be normalized_endpoint or contrastive")
    return posterior, log_marginal


def condition_dict(condition: SyntheticCondition) -> dict[str, object]:
    return asdict(condition)
