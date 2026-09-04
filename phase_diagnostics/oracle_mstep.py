"""Explicit ORACLE-phase M-step arms for logical experiment 057B."""

from __future__ import annotations

import csv
import gzip
import io
import math
import subprocess
import time
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr

from .likelihood import ScoreConfig, log_score
from .metrics import (
    chromosome_swap_invariant_metrics,
    homolog_separation_correlation,
    structure_metrics,
)
from .p9016 import PhaseContactData, read_3dg
from .provenance import file_sha256


STATE_COLUMNS = ("phase_prob00", "phase_prob01", "phase_prob10", "phase_prob11")


def _chrom_key(chrom: str) -> tuple[int, str]:
    token = chrom[3:] if chrom.startswith("chr") else chrom
    return (int(token), "") if token.isdigit() else (10_000, token)


def _deterministic_gzip(path: Path):
    raw = path.open("wb")
    zipped = gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0)
    text = io.TextIOWrapper(zipped, encoding="utf-8", newline="")
    return raw, zipped, text


def filter_onehot_cis(source: Path, destination: Path) -> int:
    """Preserve headers and measured one-hot phase while retaining cis rows only."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] | None = None
    count = 0
    raw, zipped, output = _deterministic_gzip(destination)
    try:
        with gzip.open(source, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    output.write(line)
                    if line.startswith("#columns:"):
                        columns = line.split(":", 1)[1].strip().split()
                    continue
                if not line.strip():
                    continue
                if columns is None:
                    raise ValueError(f"{source}: missing #columns header")
                fields = line.rstrip("\n").split("\t")
                by_name = dict(zip(columns, fields, strict=True))
                if by_name["chr1"] == by_name["chr2"]:
                    output.write(line)
                    count += 1
    finally:
        output.close()
        if not zipped.closed:
            zipped.close()
        if not raw.closed:
            raw.close()
    if count == 0:
        raise ValueError("cis-only ORACLE input is empty")
    return count


def _run_hickit(
    binary: Path, source: Path, destination: Path, log_path: Path, seed: int,
) -> tuple[float, str]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(binary), "--min-leg-dist=0", "--dup-dist=0", "-P", "2",
        "-s", str(seed), "-n", "1000", "-M", "--fdg-backend=cpu",
        "-i", str(source), "-p", "1.0", "-S",
        "-r", "1m", "-c", "1", "-r", "10m", "-c", "5",
        "-b", "4m", "-b", "1m", "-O", str(destination),
    ]
    start = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    runtime = time.perf_counter() - start
    if completed.returncode != 0 or not destination.is_file():
        raise RuntimeError(f"Hickit ORACLE M-step failed; see {log_path}")
    return runtime, " ".join(command)


def _coordinate_tensor(
    reference: dict[tuple[str, int, int], np.ndarray],
    estimate: dict[tuple[str, int, int], np.ndarray],
    chrom: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    starts = sorted({
        start for candidate_chrom, start, copy in reference
        if candidate_chrom == chrom and copy == 0
        and all((chrom, start, candidate_copy) in reference and (chrom, start, candidate_copy) in estimate
                for candidate_copy in (0, 1))
    })
    truth = np.asarray([[reference[(chrom, start, copy)] for copy in (0, 1)] for start in starts])
    inferred = np.asarray([[estimate[(chrom, start, copy)] for copy in (0, 1)] for start in starts])
    return np.asarray(starts, dtype=np.int64), truth, inferred


def _mean_adjacent_distance(coordinates: dict[tuple[str, int, int], np.ndarray]) -> float:
    values: list[float] = []
    for chrom in sorted({key[0] for key in coordinates}, key=_chrom_key):
        for copy in (0, 1):
            starts = sorted(start for candidate, start, candidate_copy in coordinates if candidate == chrom and candidate_copy == copy)
            values.extend(
                float(np.linalg.norm(coordinates[(chrom, right, copy)] - coordinates[(chrom, left, copy)]))
                for left, right in zip(starts, starts[1:], strict=False)
            )
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite) & (finite > 0)]
    return float(np.median(finite)) if finite.size else 1.0


def _contact_energy(
    data: PhaseContactData,
    estimate: dict[tuple[str, int, int], np.ndarray],
    selected: np.ndarray,
) -> tuple[float, int]:
    distance: list[float] = []
    multiplicity: list[int] = []
    for index in np.flatnonzero(selected):
        state = int(data.truth_state[index])
        left = (str(data.chrom1[index]), int(data.start1[index]), state >> 1)
        right = (str(data.chrom2[index]), int(data.start2[index]), state & 1)
        if left not in estimate or right not in estimate:
            continue
        distance.append(float(np.linalg.norm(estimate[left] - estimate[right])))
        multiplicity.append(int(data.multiplicity[index]))
    if not distance:
        return math.nan, 0
    unit = _mean_adjacent_distance(estimate)
    normalized = np.asarray(distance) / max(unit, 1e-12) * np.cbrt(np.asarray(multiplicity))
    energy = -log_score(normalized, ScoreConfig(mode="fdg_flat"))
    return float(np.mean(energy)), len(distance)


def evaluate_structure(
    reference: dict[tuple[str, int, int], np.ndarray],
    estimate: dict[tuple[str, int, int], np.ndarray],
    contacts: PhaseContactData,
    *,
    arm: str,
    seed: int,
    runtime_seconds: float,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    per_chromosome: list[dict[str, object]] = []
    centroids_reference: list[np.ndarray] = []
    centroids_estimate: list[np.ndarray] = []
    for chrom in sorted({key[0] for key in reference} & {key[0] for key in estimate}, key=_chrom_key):
        starts, truth, inferred = _coordinate_tensor(reference, estimate, chrom)
        if starts.size < 4:
            continue
        metrics, swap = chromosome_swap_invariant_metrics(truth, inferred)
        aligned_gauge = inferred[:, ::-1] if swap else inferred
        per_copy = [
            structure_metrics(truth[:, copy], aligned_gauge[:, copy])["distance_spearman"]
            for copy in (0, 1)
        ]
        per_chromosome.append({
            "component": "057B_oracle_phase_mstep",
            "configuration": f"{arm}_seed{seed}",
            "chromosome": chrom,
            "n_contacts": int(np.sum((contacts.chrom1 == chrom) | (contacts.chrom2 == chrom))),
            "n_bins": int(starts.size),
            "geometry_selected_copy_flip_eval_only": swap,
            "procrustes_rmsd": metrics["procrustes_rmsd"],
            "rigid_rmsd": metrics["rigid_rmsd"],
            "distance_pearson": metrics["distance_pearson"],
            "distance_spearman": metrics["distance_spearman"],
            "copy0_distance_spearman": per_copy[0],
            "copy1_distance_spearman": per_copy[1],
            "chromosome_shape_preservation": float(np.nanmean(per_copy)),
            "homolog_separation_spearman": homolog_separation_correlation(truth, aligned_gauge),
            "truth_mean_homolog_separation": float(np.linalg.norm(truth[:, 0] - truth[:, 1], axis=1).mean()),
            "reconstruction_mean_homolog_separation": float(np.linalg.norm(
                aligned_gauge[:, 0] - aligned_gauge[:, 1], axis=1
            ).mean()),
            "truth_used_in_inference": "ORACLE_CONTACT_PHASE",
        })
        centroids_reference.append(truth.reshape(-1, 3).mean(axis=0))
        centroids_estimate.append(aligned_gauge.reshape(-1, 3).mean(axis=0))
    if not per_chromosome:
        raise ValueError(f"no chromosomes could be evaluated for {arm}")
    cis_energy, n_cis = _contact_energy(contacts, estimate, contacts.is_cis)
    trans_energy, n_trans = _contact_energy(contacts, estimate, ~contacts.is_cis)
    ref_centroid_dist = pdist(np.asarray(centroids_reference))
    est_centroid_dist = pdist(np.asarray(centroids_estimate))
    summary: dict[str, object] = {
        "arm": arm,
        "seed": seed,
        "n_chromosomes": len(per_chromosome),
        "mean_copy_swap_invariant_procrustes_rmsd": float(np.mean([float(row["procrustes_rmsd"]) for row in per_chromosome])),
        "mean_copy_swap_invariant_rigid_rmsd": float(np.mean([float(row["rigid_rmsd"]) for row in per_chromosome])),
        "mean_copy_swap_invariant_distance_pearson": float(np.nanmean([float(row["distance_pearson"]) for row in per_chromosome])),
        "mean_copy_swap_invariant_distance_spearman": float(np.nanmean([float(row["distance_spearman"]) for row in per_chromosome])),
        "mean_chromosome_shape_preservation": float(np.nanmean([float(row["chromosome_shape_preservation"]) for row in per_chromosome])),
        "mean_homolog_separation_spearman": float(np.nanmean([float(row["homolog_separation_spearman"]) for row in per_chromosome])),
        "cis_truth_state_fdg_energy_proxy": cis_energy,
        "trans_truth_state_fdg_energy_proxy": trans_energy,
        "n_cis_contacts_scored": n_cis,
        "n_trans_contacts_scored": n_trans,
        "chromosome_centroid_distance_spearman": float(spearmanr(ref_centroid_dist, est_centroid_dist).statistic),
        "runtime_seconds": runtime_seconds,
        "truth_used_in_inference": "ORACLE_CONTACT_PHASE",
    }
    return summary, per_chromosome


def _read_trans_onehot(source: Path, resolution: int) -> Counter[tuple[str, int, str, int, int]]:
    counts: Counter[tuple[str, int, str, int, int]] = Counter()
    columns: list[str] | None = None
    with gzip.open(source, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#columns:"):
                columns = line.split(":", 1)[1].strip().split()
                continue
            if line.startswith("#") or not line.strip():
                continue
            if columns is None:
                raise ValueError(f"{source}: missing #columns header")
            values = dict(zip(columns, line.rstrip("\n").split("\t"), strict=True))
            if values["chr1"] == values["chr2"]:
                continue
            probabilities = [float(values[name]) for name in STATE_COLUMNS]
            state = int(np.argmax(probabilities))
            if probabilities[state] != 1.0 or not np.isclose(sum(probabilities), 1.0):
                raise ValueError("ORACLE placement requires pure one-hot contact states")
            counts[(
                values["chr1"], int(values["pos1"]) // resolution * resolution,
                values["chr2"], int(values["pos2"]) // resolution * resolution, state,
            )] += 1
    return counts


def rigid_trans_placement(
    source: Path,
    cis_coordinates: dict[tuple[str, int, int], np.ndarray],
    *,
    seed: int,
    iterations: int = 350,
) -> tuple[dict[tuple[str, int, int], np.ndarray], dict[str, float | int | str]]:
    """Optimize chromosome rotations/translations while keeping every cis shape fixed."""
    import torch

    torch.manual_seed(seed)
    torch.set_num_threads(1)
    chromosomes = sorted({key[0] for key in cis_coordinates}, key=_chrom_key)
    chrom_to_index = {chrom: index for index, chrom in enumerate(chromosomes)}
    centroids = {
        chrom: np.mean([value for (candidate, _, _), value in cis_coordinates.items() if candidate == chrom], axis=0)
        for chrom in chromosomes
    }
    local = {key: value - centroids[key[0]] for key, value in cis_coordinates.items()}
    contacts = _read_trans_onehot(source, 1_000_000)
    left_values: list[np.ndarray] = []
    right_values: list[np.ndarray] = []
    left_chrom: list[int] = []
    right_chrom: list[int] = []
    weights: list[float] = []
    dscale: list[float] = []
    for (chrom1, start1, chrom2, start2, state), count in contacts.items():
        left = (chrom1, start1, state >> 1)
        right = (chrom2, start2, state & 1)
        if left not in local or right not in local:
            continue
        left_values.append(local[left])
        right_values.append(local[right])
        left_chrom.append(chrom_to_index[chrom1])
        right_chrom.append(chrom_to_index[chrom2])
        weights.append(float(count))
        dscale.append(float(count) ** (-1.0 / 3.0))
    if not weights:
        raise ValueError("no trans ORACLE contacts overlap cis-only coordinates")
    dtype = torch.float64
    left_base = torch.tensor(np.asarray(left_values), dtype=dtype)
    right_base = torch.tensor(np.asarray(right_values), dtype=dtype)
    left_index = torch.tensor(left_chrom, dtype=torch.long)
    right_index = torch.tensor(right_chrom, dtype=torch.long)
    contact_weight = torch.tensor(weights, dtype=dtype)
    contact_dscale = torch.tensor(dscale, dtype=dtype)
    unit = _mean_adjacent_distance(cis_coordinates)
    initial_translation = np.asarray([centroids[chrom] for chrom in chromosomes])
    rotation_free = torch.nn.Parameter(torch.zeros((len(chromosomes) - 1, 3), dtype=dtype))
    translation_free = torch.nn.Parameter(torch.tensor(initial_translation[1:], dtype=dtype))
    anchor_translation = torch.tensor(initial_translation[:1], dtype=dtype)
    optimizer = torch.optim.Adam((rotation_free, translation_free), lr=0.025)
    identity = torch.eye(3, dtype=dtype)

    def rotation_matrices(vector: torch.Tensor) -> torch.Tensor:
        angle = torch.linalg.norm(vector, dim=1, keepdim=True).clamp_min(1e-12)
        axis = vector / angle
        zero = torch.zeros_like(axis[:, 0])
        x, y, z = axis.unbind(dim=1)
        skew = torch.stack((zero, -z, y, z, zero, -x, -y, x, zero), dim=1).reshape(-1, 3, 3)
        sine = torch.sin(angle)[:, :, None]
        cosine = torch.cos(angle)[:, :, None]
        return identity[None] + sine * skew + (1.0 - cosine) * (skew @ skew)

    def energy() -> torch.Tensor:
        rotation = torch.cat((identity[None], rotation_matrices(rotation_free)), dim=0)
        translation = torch.cat((anchor_translation, translation_free), dim=0)
        left = torch.einsum("nij,nj->ni", rotation[left_index], left_base) + translation[left_index]
        right = torch.einsum("nij,nj->ni", rotation[right_index], right_base) + translation[right_index]
        ratio = torch.linalg.norm(left - right, dim=1).clamp_min(1e-12) / (unit * contact_dscale)
        shell = torch.where(
            ratio < 0.5,
            torch.square(0.5 - ratio),
            torch.where(
                ratio <= 1.5,
                torch.zeros_like(ratio),
                torch.where(
                    ratio <= 2.0,
                    torch.square(ratio - 1.5),
                    1.5 * (ratio - 2.0) + 0.125 / (ratio - 1.5).clamp_min(1e-6),
                ),
            ),
        )
        contact_term = torch.sum(contact_weight * shell) / torch.sum(contact_weight)
        center_distance = torch.pdist(translation)
        center_repulsion = torch.square(torch.relu(2.0 * unit - center_distance) / max(unit, 1e-12)).mean()
        confinement = torch.square(translation - translation.mean(dim=0)).mean() / max(unit * unit, 1e-12)
        return contact_term + 0.02 * center_repulsion + 1e-5 * confinement

    with torch.no_grad():
        initial_objective = float(energy())
    history: list[float] = []
    start_time = time.perf_counter()
    for _ in range(iterations):
        optimizer.zero_grad()
        loss = energy()
        if not bool(torch.isfinite(loss)):
            raise FloatingPointError("non-finite rigid trans-placement objective")
        loss.backward()
        optimizer.step()
        history.append(float(loss.detach()))
    runtime = time.perf_counter() - start_time
    with torch.no_grad():
        rotation = torch.cat((identity[None], rotation_matrices(rotation_free)), dim=0).numpy()
        translation = torch.cat((anchor_translation, translation_free), dim=0).numpy()
    result = {
        key: rotation[chrom_to_index[key[0]]] @ value + translation[chrom_to_index[key[0]]]
        for key, value in local.items()
    }
    final_objective = history[-1]
    return result, {
        "optimizer": "fixed_shape_rigid_body_adam_fdg_layout",
        "iterations": iterations,
        "n_unique_trans_state_contacts": len(weights),
        "n_raw_trans_contact_mass": int(sum(weights)),
        "initial_objective": initial_objective,
        "final_objective": final_objective,
        "runtime_seconds": runtime,
        "converged": int(final_objective < initial_objective and np.isfinite(final_objective)),
        "seed": seed,
    }


def run_native_oracle_mstep(
    repo: Path,
    output: Path,
    phase_contacts: PhaseContactData,
    reference_path: Path,
    phase_known_input: Path,
    seeds: Iterable[int],
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    work = output / "work" / "057B"
    work.mkdir(parents=True, exist_ok=True)
    build_log = work / "build_hickit.log"
    with build_log.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            ("make", "-j4", "hickit"), cwd=repo, stdout=log, stderr=subprocess.STDOUT, text=True
        )
    binary = repo / "hickit"
    if completed.returncode != 0 or not binary.is_file():
        raise RuntimeError(f"failed to build current Hickit; see {build_log}")
    cis_input = work / "train.true_phase.cis_only.onehot.pairs.gz"
    n_cis_train = filter_onehot_cis(phase_known_input, cis_input)
    reference = read_3dg(reference_path, 1_000_000)
    summaries: list[dict[str, object]] = []
    chromosome_rows: list[dict[str, object]] = []
    commands: list[str] = []
    seed_values = list(seeds)
    for seed in seed_values:
        cis_path = work / f"cis_only_seed{seed}.1m.3dg"
        joint_path = work / f"cis_plus_trans_seed{seed}.1m.3dg"
        cis_runtime, cis_command = _run_hickit(
            binary, cis_input, cis_path, work / f"cis_only_seed{seed}.log", seed
        )
        joint_runtime, joint_command = _run_hickit(
            binary, phase_known_input, joint_path, work / f"cis_plus_trans_seed{seed}.log", seed
        )
        commands.extend((cis_command, joint_command))
        cis_map = read_3dg(cis_path, 1_000_000)
        joint_map = read_3dg(joint_path, 1_000_000)
        for arm, coordinate_map, runtime in (
            ("native_hickit_cis_only", cis_map, cis_runtime),
            ("native_hickit_cis_plus_trans", joint_map, joint_runtime),
        ):
            summary, per_chrom = evaluate_structure(
                reference, coordinate_map, phase_contacts, arm=arm, seed=seed,
                runtime_seconds=runtime,
            )
            summaries.append(summary)
            chromosome_rows.extend(per_chrom)
        placed, placement = rigid_trans_placement(phase_known_input, cis_map, seed=seed)
        summary, per_chrom = evaluate_structure(
            reference, placed, phase_contacts,
            arm="experimental_fixed_cis_rigid_trans_placement", seed=seed,
            runtime_seconds=float(placement["runtime_seconds"]),
        )
        summary.update({f"placement_{key}": value for key, value in placement.items()})
        summaries.append(summary)
        chromosome_rows.extend(per_chrom)
    manifest: dict[str, object] = {
        "status": "COMPLETE_THREE_EXPLICIT_ARMS",
        "backend": "cpu",
        "resolution_bp": 1_000_000,
        "seeds": ",".join(str(seed) for seed in seed_values),
        "phase_known_input": str(phase_known_input),
        "phase_known_input_sha256": file_sha256(phase_known_input),
        "cis_only_input_sha256": file_sha256(cis_input),
        "n_cis_train_contacts": n_cis_train,
        "hickit_binary_sha256": file_sha256(binary),
        "commands": " || ".join(commands),
        "fixed_shape_arm": "experimental rigid-body FDG layout; cis coordinates are immutable",
    }
    return summaries, chromosome_rows, manifest
