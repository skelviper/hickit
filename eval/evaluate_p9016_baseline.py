#!/usr/bin/env python3
"""Evaluate the P9016 softall baseline after blind training has finished."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.spatial import ConvexHull, QhullError
from scipy.stats import rankdata


COPY_RE = re.compile(r"^(?P<chrom>.+)\((?P<copy>mat|pat|copy0|copy1|0|1)\)$")
STATE_LABELS = ("00", "01", "10", "11")
STATE_COLUMNS = ("p00", "p01", "p10", "p11")
PANEL_SIZE_IN = 3.0
PLOT_FONT_SIZE = 7
SCATTER_POINT_SIZE = 3
MAX_AXIS_TICKS = 5
RETINA_DPI = 300
SCATTER_PANEL_GAP_IN = 0.80
SCATTER_LEGEND_WIDTH_IN = 2.15
SCATTER_CHR1_LEGEND_WIDTH_IN = 1.15
BACKGROUND_POINT_SIZE = 1.2
BACKGROUND_POINT_ALPHA = 0.14
PMAX_THRESHOLD = 0.90
FDG_D_C1 = 0.5
FDG_D_C2 = 1.5
FDG_D_C3 = 2.0
FDG_C_C1 = 3.0 * (FDG_D_C3 - FDG_D_C2)
FDG_C_C2 = (FDG_D_C3 - FDG_D_C2) ** 3
COPY0_COLORS = [
    "#4E79A7", "#F28E2B", "#59A14F", "#B6992D", "#499894",
    "#E15759", "#79706E", "#D37295", "#B07AA1", "#9D7660",
    "#5B8DB8", "#EFA65A", "#79A95B", "#C6A54A", "#65A8A5",
    "#D66B6B", "#8C857D", "#C785A8", "#8F7BB8", "#B08A73",
]
COPY1_COLORS = [
    "#2A9D8F", "#C77DFF", "#E9C46A", "#4361EE", "#F28482",
    "#6A994E", "#FF9F1C", "#4CC9F0", "#BC4749", "#90BE6D",
    "#9B5DE5", "#F15BB5", "#00B4D8", "#F4A261", "#577590",
    "#B56576", "#80B918", "#3A86FF", "#C1121F", "#8AC926",
]


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    if chrom.startswith("chr"):
        tail = chrom[3:]
    else:
        tail = chrom
    if tail.isdigit():
        return (int(tail), "")
    if tail == "X":
        return (10_000, "")
    if tail == "Y":
        return (10_001, "")
    if tail in {"M", "MT"}:
        return (10_002, "")
    return (20_000, tail)


def copy_to_int(value: str) -> int:
    if value in {"mat", "copy0", "0"}:
        return 0
    if value in {"pat", "copy1", "1"}:
        return 1
    raise ValueError(f"unsupported copy label: {value}")


def split_chrom_copy(raw_chrom: str) -> tuple[str, int]:
    match = COPY_RE.match(raw_chrom)
    if match:
        return match.group("chrom"), copy_to_int(match.group("copy"))
    return raw_chrom, 0


def finite_float(value: str) -> float:
    x = float(value)
    if not math.isfinite(x):
        raise ValueError(f"non-finite coordinate: {value}")
    return x


def read_reconstruction(path: Path) -> dict[tuple[str, int, int], np.ndarray]:
    coords: dict[tuple[str, int, int], np.ndarray] = {}
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"chr", "start", "copy", "x", "y", "z"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"{path} is missing required columns: {sorted(required)}")
        for row in reader:
            chrom = row["chr"]
            start = int(row["start"])
            copy = int(row["copy"])
            if copy not in (0, 1):
                continue
            coords[(chrom, start, copy)] = np.array(
                [finite_float(row["x"]), finite_float(row["y"]), finite_float(row["z"])],
                dtype=float,
            )
    return coords


def read_reference_3dg(path: Path, bin_size: int) -> dict[tuple[str, int, int], np.ndarray]:
    grouped: dict[tuple[str, int, int], list[np.ndarray]] = defaultdict(list)
    with open_text(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split()
            if len(fields) < 5:
                raise ValueError(f"{path}:{line_no} expected at least 5 columns")
            chrom, copy = split_chrom_copy(fields[0])
            start = (int(fields[1]) // bin_size) * bin_size
            xyz = np.array([finite_float(fields[2]), finite_float(fields[3]), finite_float(fields[4])], dtype=float)
            grouped[(chrom, start, copy)].append(xyz)
    coords: dict[tuple[str, int, int], np.ndarray] = {}
    for key, values in grouped.items():
        coords[key] = np.vstack(values).mean(axis=0)
    return coords


def read_train_manifest(path: Path | None) -> dict[str, str]:
    if path is None:
        raise ValueError("--train-manifest is required for a discipline-compliant eval")
    if not path.exists():
        raise FileNotFoundError(path)
    values: dict[str, str] = {}
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header != ["key", "value"]:
            raise ValueError(f"{path} is not a key/value manifest")
        for row in reader:
            if len(row) >= 2:
                values[row[0]] = row[1]
    return values


def read_posterior(path: Path) -> dict[tuple[str, int, str, int], dict[str, object]]:
    posterior: dict[tuple[str, int, str, int], dict[str, object]] = {}
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"chr1", "start1", "chr2", "start2", "base_d_scale", "base_k", "p00", "p01", "p10", "p11", "pmax", "contact_class"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"{path} is missing required columns: {sorted(required)}")
        for row in reader:
            chrom1 = row["chr1"]
            chrom2 = row["chr2"]
            start1 = int(row["start1"])
            start2 = int(row["start2"])
            p4 = np.array([float(row[col]) for col in STATE_COLUMNS], dtype=float)
            key = canonical_bpair_key(chrom1, start1, chrom2, start2)
            posterior[key] = {
                "chrom1": chrom1,
                "start1": start1,
                "chrom2": chrom2,
                "start2": start2,
                "base_d_scale": float(row["base_d_scale"]),
                "base_k": float(row["base_k"]),
                "p4": p4,
                "pmax": float(row["pmax"]),
                "contact_class": row["contact_class"],
            }
    return posterior


def manifest_int(manifest: dict[str, str], key: str) -> int:
    if key not in manifest:
        raise ValueError(f"training manifest missing {key}")
    try:
        return int(manifest[key])
    except ValueError as exc:
        raise ValueError(f"training manifest {key} is not an integer: {manifest[key]}") from exc


def validate_train_manifest(manifest: dict[str, str], bin_size: int) -> None:
    required = {
        "sample": "P9016",
        "runner_family": "p9016_minimal",
        "baseline": "softall",
        "mstep_graph_mode": "raw_expected_soft_all",
        "input_contact_source": "raw_pairs",
        "status": "OK",
    }
    for key, expected in required.items():
        got = manifest.get(key)
        if got != expected:
            raise ValueError(f"training manifest {key}={got!r}; expected {expected!r}")
    if manifest_int(manifest, "resolution") != bin_size:
        raise ValueError("training manifest resolution disagrees with --bin-size")
    if manifest_int(manifest, "bin_size_bp") != bin_size:
        raise ValueError("training manifest bin_size_bp disagrees with --bin-size")
    for key in ("uses_phase_labels", "uses_charm_or_reference", "uses_charm_for_training"):
        if manifest_int(manifest, key) != 0:
            raise ValueError(f"training manifest {key} must be 0")
    if manifest_int(manifest, "copy_labels_are_gauge_only") != 1:
        raise ValueError("training manifest copy_labels_are_gauge_only must be 1")
    for key in ("n_raw", "n_bpair", "n_beads"):
        if manifest_int(manifest, key) <= 0:
            raise ValueError(f"training manifest {key} must be positive")


def summarize_pairs(path: Path) -> tuple[dict[str, object], dict[str, int]]:
    chrom_sizes: dict[str, int] = {}
    columns: list[str] | None = None
    stats = {
        "pairs_total": 0,
        "pairs_cis": 0,
        "pairs_trans": 0,
        "pairs_with_any_phase": 0,
        "pairs_with_both_phases": 0,
        "phase0_labeled": 0,
        "phase1_labeled": 0,
    }
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#chromosome:"):
                parts = line.split()
                if len(parts) >= 3:
                    chrom_sizes[parts[1]] = int(parts[2])
                continue
            if line.startswith("#columns:"):
                columns = line.split(":", 1)[1].strip().split()
                continue
            if line.startswith("#"):
                continue
            fields = line.split()
            if columns is None:
                columns = ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]
            row = dict(zip(columns, fields))
            stats["pairs_total"] += 1
            if row.get("chr1") == row.get("chr2"):
                stats["pairs_cis"] += 1
            else:
                stats["pairs_trans"] += 1
            p0 = row.get("phase0", ".")
            p1 = row.get("phase1", ".")
            has0 = p0 in {"0", "1"}
            has1 = p1 in {"0", "1"}
            if has0:
                stats["phase0_labeled"] += 1
            if has1:
                stats["phase1_labeled"] += 1
            if has0 or has1:
                stats["pairs_with_any_phase"] += 1
            if has0 and has1:
                stats["pairs_with_both_phases"] += 1
    stats["pairs_phase_columns_present"] = int(columns is not None and "phase0" in columns and "phase1" in columns)
    return stats, chrom_sizes


def canonical_bpair_key(chrom1: str, start1: int, chrom2: str, start2: int) -> tuple[str, int, str, int]:
    a = (chrom_sort_key(chrom1), start1, chrom1)
    b = (chrom_sort_key(chrom2), start2, chrom2)
    if a <= b:
        return chrom1, start1, chrom2, start2
    return chrom2, start2, chrom1, start1


def canonicalize_contact(
    chrom1: str,
    start1: int,
    phase1: int,
    chrom2: str,
    start2: int,
    phase2: int,
) -> tuple[tuple[str, int, str, int], int]:
    key = canonical_bpair_key(chrom1, start1, chrom2, start2)
    if key == (chrom1, start1, chrom2, start2):
        return key, 2 * phase1 + phase2
    return key, 2 * phase2 + phase1


def read_contact_truth_counts(
    path: Path,
    bin_size: int,
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
) -> dict[tuple[str, int, str, int], np.ndarray]:
    columns: list[str] | None = None
    counts: dict[tuple[str, int, str, int], np.ndarray] = defaultdict(lambda: np.zeros(4, dtype=np.int64))
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#columns:"):
                columns = line.split(":", 1)[1].strip().split()
                continue
            if line.startswith("#"):
                continue
            fields = line.split()
            if columns is None:
                columns = ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]
            row = dict(zip(columns, fields))
            p0 = row.get("phase0", ".")
            p1 = row.get("phase1", ".")
            if p0 not in {"0", "1"} or p1 not in {"0", "1"}:
                continue
            chrom1 = row["chr1"]
            chrom2 = row["chr2"]
            start1 = (int(row["pos1"]) // bin_size) * bin_size
            start2 = (int(row["pos2"]) // bin_size) * bin_size
            if chrom1 == chrom2 and start1 == start2:
                continue
            key, state = canonicalize_contact(chrom1, start1, int(p0), chrom2, start2, int(p1))
            if key in posterior:
                counts[key][state] += 1
    return dict(counts)


def ordered_keys(
    reference: dict[tuple[str, int, int], np.ndarray],
    reconstruction: dict[tuple[str, int, int], np.ndarray],
    swaps_by_chrom: dict[str, int],
    chrom: str | None = None,
) -> list[tuple[str, int, int]]:
    keys: list[tuple[str, int, int]] = []
    for key in reference:
        k_chrom, start, ref_copy = key
        if chrom is not None and k_chrom != chrom:
            continue
        swap = swaps_by_chrom.get(k_chrom, 0)
        recon_key = (k_chrom, start, ref_copy ^ swap)
        if recon_key in reconstruction:
            keys.append(key)
    keys.sort(key=lambda k: (chrom_sort_key(k[0]), k[1], k[2]))
    return keys


def radius_of_gyration(x: np.ndarray) -> float:
    if len(x) == 0:
        return float("nan")
    centered = x - x.mean(axis=0, keepdims=True)
    return float(np.sqrt(np.mean(np.sum(centered * centered, axis=1))))


def condensed_distances(x: np.ndarray) -> np.ndarray:
    n = x.shape[0]
    if n < 2:
        return np.array([], dtype=float)
    diff = x[:, None, :] - x[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    iu = np.triu_indices(n, k=1)
    return dist[iu]


def pearson_corr(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def spearman_corr(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    return pearson_corr(rankdata(a), rankdata(b))


def fit_procrustes(ref: np.ndarray, query: np.ndarray) -> tuple[float, float, np.ndarray, np.ndarray, np.ndarray]:
    """Fit rigid Procrustes alignment and report similarity scale as a diagnostic.

    The standard alignment intentionally applies no scale change. The returned
    scale is the optimal similarity-Procrustes scale that would have been used by
    a scaled fit, and is kept only to diagnose global scale mismatch.
    """
    if len(ref) == 0:
        identity = np.eye(3)
        origin = np.zeros((1, 3), dtype=float)
        return float("nan"), float("nan"), identity, origin, origin
    ref_center = ref.mean(axis=0, keepdims=True)
    query_center = query.mean(axis=0, keepdims=True)
    ref0 = ref - ref_center
    query0 = query - query_center
    covariance = query0.T @ ref0
    u, singular, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = u @ vt
    denom = float(np.sum(query0 * query0))
    similarity_scale = float(np.sum(singular) / denom) if denom > 0 else float("nan")
    aligned = apply_procrustes(query, 1.0, rotation, query_center, ref_center)
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - ref) ** 2, axis=1))))
    ref_rg = radius_of_gyration(ref)
    if ref_rg > 0:
        rmsd = rmsd / ref_rg
    return rmsd, similarity_scale, rotation, query_center, ref_center


def apply_procrustes(
    query: np.ndarray,
    scale: float,
    rotation: np.ndarray,
    query_center: np.ndarray,
    ref_center: np.ndarray,
) -> np.ndarray:
    if not math.isfinite(scale):
        return query.copy()
    return scale * (query - query_center) @ rotation + ref_center


def procrustes(ref: np.ndarray, query: np.ndarray) -> tuple[float, float, np.ndarray]:
    rmsd, similarity_scale, rotation, query_center, ref_center = fit_procrustes(ref, query)
    aligned = apply_procrustes(query, 1.0, rotation, query_center, ref_center)
    return rmsd, similarity_scale, aligned


def points_for_keys(
    keys: list[tuple[str, int, int]],
    coords: dict[tuple[str, int, int], np.ndarray],
    query_coords: dict[tuple[str, int, int], np.ndarray] | None = None,
    swaps_by_chrom: dict[str, int] | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    left = np.vstack([coords[k] for k in keys]) if keys else np.zeros((0, 3), dtype=float)
    if query_coords is None or swaps_by_chrom is None:
        return left, None
    right = np.vstack([query_coords[(k[0], k[1], k[2] ^ swaps_by_chrom.get(k[0], 0))] for k in keys]) if keys else np.zeros((0, 3), dtype=float)
    return left, right


def cis_distance_correlation_row(
    chrom: str,
    swap: int,
    reference: dict[tuple[str, int, int], np.ndarray],
    reconstruction: dict[tuple[str, int, int], np.ndarray],
) -> dict[str, object]:
    swaps = {chrom: swap}
    keys = ordered_keys(reference, reconstruction, swaps, chrom=chrom)
    ref, rec = points_for_keys(keys, reference, reconstruction, swaps)
    assert rec is not None
    ref_dist = condensed_distances(ref)
    rec_dist = condensed_distances(rec)
    ref_rg = radius_of_gyration(ref)
    rec_rg = radius_of_gyration(rec)
    ref_norm = ref_dist / ref_rg if ref_rg > 0 else ref_dist * float("nan")
    rec_norm = rec_dist / rec_rg if rec_rg > 0 else rec_dist * float("nan")
    return {
        "chrom": chrom,
        "copy_swap": swap,
        "n_points": len(keys),
        "n_pairwise_distances": len(ref_dist),
        "cis_distance_pearson": pearson_corr(ref_norm, rec_norm),
        "cis_distance_spearman": spearman_corr(ref_norm, rec_norm),
    }


def choose_distance_swaps(
    chroms: list[str],
    reference: dict[tuple[str, int, int], np.ndarray],
    reconstruction: dict[tuple[str, int, int], np.ndarray],
) -> tuple[dict[str, int], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    swaps: dict[str, int] = {}
    for chrom in chroms:
        candidates = [cis_distance_correlation_row(chrom, swap, reference, reconstruction) for swap in (0, 1)]
        best = max(candidates, key=lambda row: (nan_to_neg_inf(row["cis_distance_spearman"]), int(row["n_points"])))
        swaps[chrom] = int(best["copy_swap"])
        for row in candidates:
            out = dict(row)
            out["selected_for_eval"] = int(int(row["copy_swap"]) == swaps[chrom])
            rows.append(out)
    return swaps, rows


def nan_to_neg_inf(value: object) -> float:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return float("-inf")
    return x if math.isfinite(x) else float("-inf")


def swap_state(state: int, swap1: int, swap2: int) -> int:
    return ((state >> 1) ^ swap1) * 2 + ((state & 1) ^ swap2)


def align_p4_to_truth_gauge(p4: np.ndarray, swap1: int, swap2: int) -> np.ndarray:
    aligned = np.zeros(4, dtype=float)
    for model_state in range(4):
        aligned[swap_state(model_state, swap1, swap2)] = p4[model_state]
    return aligned


def fdg_contact_energy_r(r: float, k: float) -> float:
    if r < FDG_D_C1:
        t = FDG_D_C1 - r
        return k * t * t
    if r <= FDG_D_C2:
        return 0.0
    if r <= FDG_D_C3:
        t = r - FDG_D_C2
        return k * t * t
    t = r - FDG_D_C2
    return k * (FDG_C_C1 * (r - FDG_D_C3) + FDG_C_C2 / t)


def p4_from_coords_for_bpair(
    coords: dict[tuple[str, int, int], np.ndarray],
    item: dict[str, object],
) -> np.ndarray | None:
    chrom1 = str(item["chrom1"])
    chrom2 = str(item["chrom2"])
    start1 = int(item["start1"])
    start2 = int(item["start2"])
    d_scale = float(item["base_d_scale"])
    base_k = float(item["base_k"])
    if d_scale <= 0:
        return None
    distances: list[float] = []
    for state in range(4):
        copy1 = state >> 1
        copy2 = state & 1
        key1 = (chrom1, start1, copy1)
        key2 = (chrom2, start2, copy2)
        if key1 not in coords or key2 not in coords:
            return None
        distances.append(float(np.linalg.norm(coords[key1] - coords[key2])))
    energy = np.array([fdg_contact_energy_r(distance / d_scale, base_k) for distance in distances], dtype=float)
    score = -energy
    score -= np.max(score)
    exp_score = np.exp(score)
    return exp_score / exp_score.sum()


def update_accuracy_stats(stats: dict[str, float], counts: np.ndarray, p4_truth_gauge: np.ndarray) -> None:
    n = int(counts.sum())
    if n <= 0:
        return
    pred = int(np.argmax(p4_truth_gauge))
    pmax = float(np.max(p4_truth_gauge))
    correct = int(counts[pred])
    stats["contacts"] += n
    stats["correct_top1"] += correct
    if pmax >= PMAX_THRESHOLD:
        stats["called_contacts"] += n
        stats["called_correct"] += correct


def empty_accuracy_stats() -> dict[str, float]:
    return {"contacts": 0.0, "correct_top1": 0.0, "called_contacts": 0.0, "called_correct": 0.0}


def finalize_accuracy_stats(stats: dict[str, float]) -> dict[str, object]:
    contacts = int(stats["contacts"])
    called = int(stats["called_contacts"])
    correct = int(stats["correct_top1"])
    called_correct = int(stats["called_correct"])
    return {
        "n_eval_contacts": contacts,
        "top1_correct_contacts": correct,
        "top1_accuracy": correct / contacts if contacts else float("nan"),
        "pmax_threshold": PMAX_THRESHOLD,
        "n_called_contacts": called,
        "called_contact_fraction": called / contacts if contacts else float("nan"),
        "pmax_threshold_correct_contacts": called_correct,
        "pmax_threshold_accuracy": called_correct / called if called else float("nan"),
        "pmax_threshold_recall": called_correct / contacts if contacts else float("nan"),
    }


def contact_accuracy_for_swaps(
    contact_counts: dict[tuple[str, int, str, int], np.ndarray],
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    reference: dict[tuple[str, int, int], np.ndarray],
    swaps_by_chrom: dict[str, int],
) -> list[dict[str, object]]:
    stats_by_source_scope: dict[tuple[str, str], dict[str, float]] = {}
    for key, counts in contact_counts.items():
        item = posterior.get(key)
        if item is None:
            continue
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        is_cis = chrom1 == chrom2
        model_p4 = align_p4_to_truth_gauge(
            np.asarray(item["p4"], dtype=float),
            swaps_by_chrom.get(chrom1, 0),
            swaps_by_chrom.get(chrom2, 0),
        )
        charm_p4 = p4_from_coords_for_bpair(reference, item)
        if charm_p4 is None:
            continue
        for source, p4 in [("reconstruction_posterior", model_p4), ("charm3dg_uniform_prior_fdg", charm_p4)]:
            scopes = ["genome_all", "genome_cis" if is_cis else "genome_trans"]
            if is_cis:
                scopes.append(f"{chrom1}_cis")
            for scope in scopes:
                stats = stats_by_source_scope.setdefault((source, scope), empty_accuracy_stats())
                update_accuracy_stats(stats, counts, p4)
    rows: list[dict[str, object]] = []
    def scope_sort_key(scope: str) -> tuple[int, tuple[int, str]]:
        if scope == "genome_all":
            return (0, (0, ""))
        if scope == "genome_cis":
            return (1, (0, ""))
        if scope == "genome_trans":
            return (2, (0, ""))
        return (3, chrom_sort_key(scope.removesuffix("_cis")))

    copy_swap_policy_by_source = {
        "reconstruction_posterior": "per_chrom_cis_distance_spearman_fixed",
        "charm3dg_uniform_prior_fdg": "charm3dg_reference_copy_labels",
    }
    for (source, scope), stats in sorted(stats_by_source_scope.items(), key=lambda x: (x[0][0], scope_sort_key(x[0][1]))):
        row = {
            "source": source,
            "scope": scope,
            "copy_swap_policy": copy_swap_policy_by_source.get(source, "unknown"),
        }
        row.update(finalize_accuracy_stats(stats))
        rows.append(row)
    return rows


def copy_separation_rows(
    coords_by_source: dict[str, dict[tuple[str, int, int], np.ndarray]],
    swaps_by_chrom: dict[str, int],
    chroms: list[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for source, coords in coords_by_source.items():
        genome_distances: list[float] = []
        for chrom in chroms:
            starts = sorted({start for c, start, copy in coords if c == chrom})
            distances: list[float] = []
            for start in starts:
                key0 = (chrom, start, swaps_by_chrom.get(chrom, 0) if source == "reconstruction" else 0)
                key1 = (chrom, start, (1 ^ swaps_by_chrom.get(chrom, 0)) if source == "reconstruction" else 1)
                if key0 in coords and key1 in coords:
                    distances.append(float(np.linalg.norm(coords[key0] - coords[key1])))
            genome_distances.extend(distances)
            rows.append({
                "source": source,
                "chrom": chrom,
                "n_bins": len(distances),
                "mean_copy01_separation": float(np.mean(distances)) if distances else float("nan"),
                "median_copy01_separation": float(np.median(distances)) if distances else float("nan"),
            })
        rows.append({
            "source": source,
            "chrom": "genome",
            "n_bins": len(genome_distances),
            "mean_copy01_separation": float(np.mean(genome_distances)) if genome_distances else float("nan"),
            "median_copy01_separation": float(np.median(genome_distances)) if genome_distances else float("nan"),
        })
    return rows


def point_cloud_volume(points: np.ndarray) -> float:
    if len(points) < 4:
        return float("nan")
    try:
        return float(ConvexHull(points).volume)
    except QhullError:
        return 0.0


def volume_rows(
    coords_by_source: dict[str, dict[tuple[str, int, int], np.ndarray]],
    chroms: list[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for source, coords in coords_by_source.items():
        all_points = np.vstack(list(coords.values())) if coords else np.zeros((0, 3), dtype=float)
        genome_volume = point_cloud_volume(all_points)
        for chrom in chroms:
            points = np.vstack([xyz for (c, _, _), xyz in coords.items() if c == chrom]) if any(c == chrom for c, _, _ in coords) else np.zeros((0, 3), dtype=float)
            volume = point_cloud_volume(points)
            rows.append({
                "source": source,
                "chrom": chrom,
                "n_points": len(points),
                "convex_hull_volume": volume,
                "volume_fraction_of_source_genome": volume / genome_volume if genome_volume > 0 and math.isfinite(volume) else float("nan"),
            })
        rows.append({
            "source": source,
            "chrom": "genome",
            "n_points": len(all_points),
            "convex_hull_volume": genome_volume,
            "volume_fraction_of_source_genome": 1.0 if genome_volume > 0 else float("nan"),
        })
    return rows


def format_value(value: object) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.6g}"
    return str(value)


def write_summary(path: Path, values: dict[str, object]) -> None:
    with path.open("w") as fh:
        fh.write("key\tvalue\n")
        for key in sorted(values):
            fh.write(f"{key}\t{format_value(values[key])}\n")


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def first_row(rows: list[dict[str, object]], **matches: object) -> dict[str, object] | None:
    for row in rows:
        if all(row.get(key) == value for key, value in matches.items()):
            return row
    return None


def distance_matrix_for_plot(points: np.ndarray) -> np.ndarray:
    rg = radius_of_gyration(points)
    diff = points[:, None, :] - points[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    return dist / rg if rg > 0 else dist


def color_for_chrom_copy(chrom: str, copy: int, chrom_index: dict[str, int]) -> tuple[float, float, float, float]:
    palette = COPY0_COLORS if copy == 0 else COPY1_COLORS
    return matplotlib.colors.to_rgba(palette[chrom_index[chrom] % len(palette)])


def style_axis_text(ax) -> None:
    ax.title.set_fontsize(PLOT_FONT_SIZE)
    ax.xaxis.label.set_size(PLOT_FONT_SIZE)
    ax.yaxis.label.set_size(PLOT_FONT_SIZE)
    ax.tick_params(axis="both", labelsize=PLOT_FONT_SIZE)
    if hasattr(ax, "zaxis"):
        ax.zaxis.label.set_size(PLOT_FONT_SIZE)
        ax.tick_params(axis="z", labelsize=PLOT_FONT_SIZE)


def evenly_spaced_indices(n: int, max_ticks: int = MAX_AXIS_TICKS) -> list[int]:
    if n <= 0:
        return []
    if n <= max_ticks:
        return list(range(n))
    return sorted(set(int(round(x)) for x in np.linspace(0, n - 1, max_ticks)))


def limit_3d_ticks(ax, max_ticks: int = MAX_AXIS_TICKS) -> None:
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        lo, hi = axis.get_view_interval()
        axis.set_ticks(np.linspace(lo, hi, max_ticks))


def shared_3d_limits(*point_sets: np.ndarray) -> list[tuple[float, float]]:
    points = np.vstack([points for points in point_sets if len(points)])
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    center = (mins + maxs) / 2.0
    half_width = float(np.max(maxs - mins) / 2.0)
    if half_width <= 0:
        half_width = 1.0
    half_width *= 1.04
    return [(float(center[i] - half_width), float(center[i] + half_width)) for i in range(3)]


def apply_3d_limits(ax, limits: list[tuple[float, float]]) -> None:
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_zlim(*limits[2])
    ax.set_box_aspect((1, 1, 1))


def plot_chr1_distance_maps(
    path: Path,
    reference: dict[tuple[str, int, int], np.ndarray],
    reconstruction: dict[tuple[str, int, int], np.ndarray],
    swaps_by_chrom: dict[str, int],
) -> int:
    chr1_swap = swaps_by_chrom.get("chr1", 0)
    copy_panels: list[tuple[int, list[tuple[str, int, int]], np.ndarray, np.ndarray]] = []
    all_distances: list[np.ndarray] = []
    for copy in (0, 1):
        keys = ordered_keys(reference, reconstruction, swaps_by_chrom, chrom="chr1")
        keys = [k for k in keys if k[2] == copy]
        if not keys:
            raise ValueError(f"no shared chr1 copy{copy} bins for chr1 distance-map plot")
        ref = np.vstack([reference[k] for k in keys])
        rec = np.vstack([reconstruction[(k[0], k[1], k[2] ^ chr1_swap)] for k in keys])
        ref_d = distance_matrix_for_plot(ref)
        rec_d = distance_matrix_for_plot(rec)
        copy_panels.append((copy, keys, ref_d, rec_d))
        all_distances.extend([ref_d.ravel(), rec_d.ravel()])
    vmax = float(np.nanpercentile(np.concatenate(all_distances), 99))
    cmap = plt.get_cmap("coolwarm").reversed()

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(2 * PANEL_SIZE_IN + 0.8, 2 * PANEL_SIZE_IN),
        constrained_layout=True,
    )
    im = None
    for col, (copy, keys, ref_d, rec_d) in enumerate(copy_panels):
        labels = [f"{k[1] // 1_000_000}M" for k in keys]
        tick_idx = evenly_spaced_indices(len(keys))
        model_copy = copy ^ chr1_swap
        for row, mat, title in [
            (0, ref_d, f"CHARM/3DG chr1 copy{copy}"),
            (1, rec_d, f"Reconstruction chr1 copy{copy} (model copy{model_copy})"),
        ]:
            ax = axes[row, col]
            im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=vmax, interpolation="nearest")
            ax.set_title(title, fontsize=PLOT_FONT_SIZE)
            ax.set_xticks(tick_idx, [labels[i] for i in tick_idx], rotation=45, ha="right", fontsize=PLOT_FONT_SIZE)
            ax.set_yticks(tick_idx, [labels[i] for i in tick_idx], fontsize=PLOT_FONT_SIZE)
            ax.set_xlabel(f"chr1 copy{copy} bins")
            ax.set_ylabel(f"chr1 copy{copy} bins")
            style_axis_text(ax)
    assert im is not None
    cbar = fig.colorbar(im, ax=axes, shrink=0.82)
    cbar.set_label("RG-normalized distance; larger distance is blue", fontsize=PLOT_FONT_SIZE)
    cbar.ax.tick_params(labelsize=PLOT_FONT_SIZE)
    fig.savefig(path, dpi=RETINA_DPI)
    plt.close(fig)
    return sum(len(keys) for _, keys, _, _ in copy_panels)


def add_background_scatter(
    ax,
    points: np.ndarray,
    keys: list[tuple[str, int, int]],
    highlight_chrom: str | None,
) -> None:
    if highlight_chrom is None:
        return
    idx = [i for i, key in enumerate(keys) if key[0] != highlight_chrom]
    if not idx:
        return
    sub = points[idx]
    ax.scatter(
        sub[:, 0],
        sub[:, 1],
        sub[:, 2],
        s=BACKGROUND_POINT_SIZE,
        color="#8A8A8A",
        alpha=BACKGROUND_POINT_ALPHA,
        depthshade=False,
        linewidths=0,
    )


def add_chrom_copy_scatter(
    ax,
    points: np.ndarray,
    keys: list[tuple[str, int, int]],
    chroms: list[str],
    *,
    highlight_chrom: str | None = None,
) -> None:
    chrom_index = {chrom: i for i, chrom in enumerate(chroms)}
    for chrom in chroms:
        if highlight_chrom is not None and chrom != highlight_chrom:
            continue
        for copy in (0, 1):
            idx = [i for i, key in enumerate(keys) if key[0] == chrom and key[2] == copy]
            if idx:
                sub = points[idx]
                ax.scatter(
                    sub[:, 0],
                    sub[:, 1],
                    sub[:, 2],
                    s=SCATTER_POINT_SIZE,
                    color=color_for_chrom_copy(chrom, copy, chrom_index),
                    alpha=0.84,
                )


def chrom_copy_legend(chroms: list[str]) -> list[Line2D]:
    chrom_index = {chrom: i for i, chrom in enumerate(chroms)}
    handles: list[Line2D] = []
    for copy in (0, 1):
        for chrom in chroms:
            handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    markersize=5,
                    markerfacecolor=color_for_chrom_copy(chrom, copy, chrom_index),
                    markeredgecolor=color_for_chrom_copy(chrom, copy, chrom_index),
                    label=f"{chrom} copy{copy}",
                )
            )
    return handles


def plot_3d_scatter(
    path: Path,
    reference: dict[tuple[str, int, int], np.ndarray],
    reconstruction: dict[tuple[str, int, int], np.ndarray],
    swaps_by_chrom: dict[str, int],
    *,
    chrom: str | None = None,
) -> int:
    all_keys = ordered_keys(reference, reconstruction, swaps_by_chrom)
    keys = [key for key in all_keys if chrom is None or key[0] == chrom]
    if not keys:
        raise ValueError("no shared bins for 3D scatter plot")
    all_ref = np.vstack([reference[k] for k in all_keys])
    all_rec = np.vstack([reconstruction[(k[0], k[1], k[2] ^ swaps_by_chrom.get(k[0], 0))] for k in all_keys])
    key_index = {key: i for i, key in enumerate(all_keys)}
    highlight_idx = [key_index[key] for key in keys]
    ref = all_ref[highlight_idx]

    _, _, rotation, query_center, ref_center = fit_procrustes(all_ref, all_rec)
    all_rec_aligned = apply_procrustes(all_rec, 1.0, rotation, query_center, ref_center)
    rec_aligned = all_rec_aligned[highlight_idx]

    chroms = sorted({k[0] for k in all_keys}, key=chrom_sort_key)
    legend_chroms = [chrom] if chrom else chroms
    limits = shared_3d_limits(all_ref, all_rec_aligned)
    legend_width = SCATTER_CHR1_LEGEND_WIDTH_IN if chrom else SCATTER_LEGEND_WIDTH_IN

    fig = plt.figure(
        figsize=(2 * PANEL_SIZE_IN + SCATTER_PANEL_GAP_IN + legend_width, PANEL_SIZE_IN),
        constrained_layout=False,
    )
    grid = fig.add_gridspec(
        1,
        4,
        width_ratios=[PANEL_SIZE_IN, SCATTER_PANEL_GAP_IN, PANEL_SIZE_IN, legend_width],
        wspace=0.0,
    )
    title_suffix = f" {chrom}" if chrom else ""
    panels = [(all_ref, f"CHARM/3DG{title_suffix}"), (all_rec_aligned, f"Reconstruction{title_suffix}")]
    for panel_idx, (points, title) in enumerate(panels, 1):
        grid_col = 0 if panel_idx == 1 else 2
        ax = fig.add_subplot(grid[0, grid_col], projection="3d")
        add_background_scatter(ax, points, all_keys, chrom)
        add_chrom_copy_scatter(ax, points, all_keys, chroms, highlight_chrom=chrom)
        ax.set_title(title, fontsize=PLOT_FONT_SIZE)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        apply_3d_limits(ax, limits)
        limit_3d_ticks(ax)
        style_axis_text(ax)
        ax.view_init(elev=22, azim=35)
    spacer_ax = fig.add_subplot(grid[0, 1])
    spacer_ax.axis("off")
    legend_ax = fig.add_subplot(grid[0, 3])
    legend_ax.axis("off")
    legend_ax.legend(
        handles=chrom_copy_legend(legend_chroms),
        loc="center left",
        fontsize=PLOT_FONT_SIZE,
        ncol=2,
        columnspacing=0.7,
        handlelength=0.8,
        handletextpad=0.3,
        labelspacing=0.18,
        borderpad=0.0,
        frameon=False,
    )
    fig.savefig(path, dpi=RETINA_DPI)
    plt.close(fig)
    return len(keys)


def write_readme(
    path: Path,
    args: argparse.Namespace,
    summary: dict[str, object],
    cis_rows: list[dict[str, object]],
    accuracy_rows: list[dict[str, object]],
    separation_rows: list[dict[str, object]],
    volume_rows_data: list[dict[str, object]],
) -> None:
    model_all = first_row(accuracy_rows, source="reconstruction_posterior", scope="genome_all")
    model_cis = first_row(accuracy_rows, source="reconstruction_posterior", scope="genome_cis")
    model_trans = first_row(accuracy_rows, source="reconstruction_posterior", scope="genome_trans")
    charm_all = first_row(accuracy_rows, source="charm3dg_uniform_prior_fdg", scope="genome_all")
    charm_cis = first_row(accuracy_rows, source="charm3dg_uniform_prior_fdg", scope="genome_cis")
    charm_trans = first_row(accuracy_rows, source="charm3dg_uniform_prior_fdg", scope="genome_trans")
    chr1_cis = next((row for row in cis_rows if row["chrom"] == "chr1" and int(row["selected_for_eval"]) == 1), None)
    rec_sep = first_row(separation_rows, source="reconstruction", chrom="genome")
    ref_sep = first_row(separation_rows, source="charm3dg", chrom="genome")
    rec_vol = first_row(volume_rows_data, source="reconstruction", chrom="genome")
    ref_vol = first_row(volume_rows_data, source="charm3dg", chrom="genome")
    rows_for_readme = [
        ("bin_size_bp", args.bin_size),
        ("pairs_total", summary["pairs_total"]),
        ("pairs_cis", summary["pairs_cis"]),
        ("pairs_trans", summary["pairs_trans"]),
        ("pairs_with_any_phase_eval_only", summary["pairs_with_any_phase"]),
        ("pairs_with_both_phases_eval_only", summary["pairs_with_both_phases"]),
        ("reference_points_aggregated", summary["reference_points_aggregated"]),
        ("reconstruction_points", summary["reconstruction_points"]),
        ("shared_points_per_chrom_best", summary["shared_points_per_chrom_best"]),
        ("chr1_shared_points_per_chrom_best", summary["chr1_shared_points_per_chrom_best"]),
        ("mean_per_chrom_cis_distance_spearman", summary["mean_per_chrom_cis_distance_spearman"]),
        ("chr1_cis_distance_spearman", chr1_cis["cis_distance_spearman"] if chr1_cis else float("nan")),
        ("model_top1_accuracy_genome_all", model_all["top1_accuracy"] if model_all else float("nan")),
        ("model_pmax90_accuracy_genome_all", model_all["pmax_threshold_accuracy"] if model_all else float("nan")),
        ("model_pmax90_recall_genome_all", model_all["pmax_threshold_recall"] if model_all else float("nan")),
        ("model_top1_accuracy_genome_cis", model_cis["top1_accuracy"] if model_cis else float("nan")),
        ("model_pmax90_accuracy_genome_cis", model_cis["pmax_threshold_accuracy"] if model_cis else float("nan")),
        ("model_pmax90_recall_genome_cis", model_cis["pmax_threshold_recall"] if model_cis else float("nan")),
        ("model_top1_accuracy_genome_trans", model_trans["top1_accuracy"] if model_trans else float("nan")),
        ("model_pmax90_accuracy_genome_trans", model_trans["pmax_threshold_accuracy"] if model_trans else float("nan")),
        ("model_pmax90_recall_genome_trans", model_trans["pmax_threshold_recall"] if model_trans else float("nan")),
        ("charm3dg_top1_accuracy_genome_all", charm_all["top1_accuracy"] if charm_all else float("nan")),
        ("charm3dg_pmax90_accuracy_genome_all", charm_all["pmax_threshold_accuracy"] if charm_all else float("nan")),
        ("charm3dg_pmax90_recall_genome_all", charm_all["pmax_threshold_recall"] if charm_all else float("nan")),
        ("charm3dg_top1_accuracy_genome_cis", charm_cis["top1_accuracy"] if charm_cis else float("nan")),
        ("charm3dg_pmax90_accuracy_genome_cis", charm_cis["pmax_threshold_accuracy"] if charm_cis else float("nan")),
        ("charm3dg_pmax90_recall_genome_cis", charm_cis["pmax_threshold_recall"] if charm_cis else float("nan")),
        ("charm3dg_top1_accuracy_genome_trans", charm_trans["top1_accuracy"] if charm_trans else float("nan")),
        ("charm3dg_pmax90_accuracy_genome_trans", charm_trans["pmax_threshold_accuracy"] if charm_trans else float("nan")),
        ("charm3dg_pmax90_recall_genome_trans", charm_trans["pmax_threshold_recall"] if charm_trans else float("nan")),
        ("reconstruction_mean_copy01_separation", rec_sep["mean_copy01_separation"] if rec_sep else float("nan")),
        ("charm3dg_mean_copy01_separation", ref_sep["mean_copy01_separation"] if ref_sep else float("nan")),
        ("reconstruction_genome_volume", rec_vol["convex_hull_volume"] if rec_vol else float("nan")),
        ("charm3dg_genome_volume", ref_vol["convex_hull_volume"] if ref_vol else float("nan")),
    ]
    with path.open("w") as fh:
        fh.write(f"# {args.outdir.name}\n\n")
        fh.write(f"- label: `{args.label}`\n")
        fh.write(f"- created_at: `{summary['created_at']}`\n")
        fh.write(f"- training input: `{args.pairs}`\n")
        fh.write(f"- reconstruction: `{args.reconstruction}`\n")
        fh.write(f"- CHARM/3DG eval reference: `{args.reference_3dg}`\n")
        fh.write(f"- train manifest: `{args.train_manifest}`\n")
        fh.write("- boundary: training used raw P9016 contact information only; phase labels and CHARM/3DG were read only by this post-training evaluator.\n")
        fh.write("- copy gauge: evaluation first selects one copy swap per chromosome by per-chromosome cis distance-matrix Spearman correlation. Contact top1 and pmax metrics then use this fixed geometry-selected gauge when comparing four-state probabilities to SNP phase truth.\n")
        fh.write("- contact denominator: contact accuracy uses eval-only raw contacts with both `phase0` and `phase1`, excluding same-bin contacts, and requiring a matching posterior bpair.\n")
        fh.write("- CHARM/3DG probability baseline: `charm3dg_uniform_prior_fdg` recomputes four-state probabilities from CHARM/3DG distances using hickit FDG contact energy, posterior `base_d_scale/base_k`, CHARM/3DG reference copy labels, and uniform four-state prior because posterior log-priors are not exported.\n\n")
        fh.write("- alignment: 3D scatter plots and Procrustes RMSD use rigid alignment only: translation and rotation are fitted, reconstruction scale is not fitted to CHARM/3DG. The reported similarity scale is diagnostic only and is not applied.\n\n")
        fh.write("## Quantitative Results\n\n")
        fh.write("| metric | value |\n")
        fh.write("| --- | ---: |\n")
        for key, value in rows_for_readme:
            fh.write(f"| {key} | {format_value(value)} |\n")
        fh.write("\n## Output Tables\n\n")
        fh.write("- `cis_distance_correlations.tsv`: per-chromosome cis distance-matrix Pearson/Spearman for both copy swaps, with the selected per-chrom swap marked.\n")
        fh.write("- `contact_accuracy.tsv`: four-state top1 accuracy, pmax >= 0.9 accuracy, called fraction, and recall for all/cis/trans contacts plus per-chromosome cis contacts. The `copy_swap_policy` column records whether rows use the reconstruction's fixed cis-distance-selected gauge or CHARM/3DG reference copy labels.\n")
        fh.write("- `copy_separation.tsv`: per-chromosome and genome mean/median distance between copy0 and copy1 of the same bin.\n")
        fh.write("- `per_chrom_volume.tsv`: per-chromosome and genome convex-hull volumes for CHARM/3DG and reconstruction.\n")
        fh.write("\n## Plots\n\n")
        fh.write("- `plots/chr1_distance_maps.png`: rows are CHARM/3DG and reconstruction; columns are copy0 and copy1; larger distances are blue.\n")
        fh.write("- `plots/all_chrom_3d_scatter.png`: each point is one bin; chromosome+copy states are colored separately; reconstruction is rigid-Procrustes-aligned into the CHARM/3DG coordinate frame with no scale fitting, and both panels share one 3D coordinate range.\n")
        fh.write("- `plots/chr1_copy_3d_scatter.png`: chr1 copy0/copy1 bins highlighted on top of low-alpha non-chr1 background points, using the same global rigid transform as the all-chromosome scatter plot.\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True, help="P9016 pairs file; phase columns are eval-only.")
    parser.add_argument("--reference-3dg", type=Path, required=True, help="CHARM/3DG reference coordinates.")
    parser.add_argument("--reconstruction", type=Path, required=True, help="Baseline reconstruction coords TSV.")
    parser.add_argument("--posterior", type=Path, default=None, help="Baseline bpair posterior TSV. Defaults to p9016_full.bpair_posterior.tsv next to reconstruction.")
    parser.add_argument("--train-manifest", type=Path, default=None, help="Training manifest TSV.")
    parser.add_argument("--outdir", type=Path, required=True, help="Experiment output directory.")
    parser.add_argument("--bin-size", type=int, default=4_000_000, help="Evaluation bin size in bp.")
    parser.add_argument("--label", default="P9016 softall baseline", help="Human-readable run label.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    plots_dir = args.outdir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    pair_stats, chrom_sizes = summarize_pairs(args.pairs)
    reference = read_reference_3dg(args.reference_3dg, args.bin_size)
    reconstruction = read_reconstruction(args.reconstruction)
    posterior_path = args.posterior if args.posterior else args.reconstruction.with_name("p9016_full.bpair_posterior.tsv")
    posterior = read_posterior(posterior_path)
    contact_counts = read_contact_truth_counts(args.pairs, args.bin_size, posterior)
    train_manifest = read_train_manifest(args.train_manifest)
    validate_train_manifest(train_manifest, args.bin_size)
    chroms = sorted({key[0] for key in reference} | {key[0] for key in reconstruction}, key=chrom_sort_key)

    distance_swaps, cis_rows = choose_distance_swaps(chroms, reference, reconstruction)
    accuracy_rows = contact_accuracy_for_swaps(contact_counts, posterior, reference, distance_swaps)
    separation_rows = copy_separation_rows({"charm3dg": reference, "reconstruction": reconstruction}, distance_swaps, chroms)
    volume_rows_data = volume_rows({"charm3dg": reference, "reconstruction": reconstruction}, chroms)
    model_all = first_row(accuracy_rows, source="reconstruction_posterior", scope="genome_all")
    model_cis = first_row(accuracy_rows, source="reconstruction_posterior", scope="genome_cis")
    model_trans = first_row(accuracy_rows, source="reconstruction_posterior", scope="genome_trans")
    charm_all = first_row(accuracy_rows, source="charm3dg_uniform_prior_fdg", scope="genome_all")
    charm_cis = first_row(accuracy_rows, source="charm3dg_uniform_prior_fdg", scope="genome_cis")
    charm_trans = first_row(accuracy_rows, source="charm3dg_uniform_prior_fdg", scope="genome_trans")

    best_keys = ordered_keys(reference, reconstruction, distance_swaps)
    best_chr1_keys = ordered_keys(reference, reconstruction, distance_swaps, chrom="chr1")
    if len(best_keys) < 2:
        raise ValueError("fewer than 2 shared genome bins after per-chrom copy swap")
    if len(best_chr1_keys) < 2:
        raise ValueError("fewer than 2 shared chr1 bins after per-chrom copy swap")
    selected_cis_rows = [row for row in cis_rows if int(row["selected_for_eval"]) == 1 and math.isfinite(float(row["cis_distance_spearman"]))]
    mean_cis_spearman = float(np.mean([float(row["cis_distance_spearman"]) for row in selected_cis_rows])) if selected_cis_rows else float("nan")

    summary: dict[str, object] = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "label": args.label,
        "bin_size_bp": args.bin_size,
        "pairs_path": str(args.pairs),
        "reference_3dg_path": str(args.reference_3dg),
        "reconstruction_path": str(args.reconstruction),
        "posterior_path": str(posterior_path),
        "train_manifest_path": str(args.train_manifest) if args.train_manifest else "",
        "chromosomes_in_pairs_header": len(chrom_sizes),
        "reference_points_aggregated": len(reference),
        "reconstruction_points": len(reconstruction),
        "posterior_bpair_rows": len(posterior),
        "truth_bpair_rows_with_eval_contacts": len(contact_counts),
        "truth_eval_contacts": int(sum(int(counts.sum()) for counts in contact_counts.values())),
        "shared_points_per_chrom_best": len(best_keys),
        "chr1_shared_points_per_chrom_best": len(best_chr1_keys),
        "distance_per_chrom_copy_swaps_json": json.dumps(distance_swaps, sort_keys=True),
        "contact_eval_per_chrom_copy_swaps_json": json.dumps(distance_swaps, sort_keys=True),
        "contact_copy_swap_source": "distance_per_chrom_copy_swaps_json",
        "copy_swap_policy": "per_chrom_best only; one geometry gauge is selected by cis distance Spearman and reused for contact top1/pmax metrics",
        "mean_per_chrom_cis_distance_spearman": mean_cis_spearman,
        "charm3dg_probability_baseline": "uniform_prior_fdg_energy_from_charm3dg_distances",
        "pmax_threshold": PMAX_THRESHOLD,
        "model_top1_accuracy_genome_all": model_all["top1_accuracy"] if model_all else float("nan"),
        "model_pmax90_accuracy_genome_all": model_all["pmax_threshold_accuracy"] if model_all else float("nan"),
        "model_pmax90_recall_genome_all": model_all["pmax_threshold_recall"] if model_all else float("nan"),
        "model_top1_accuracy_genome_cis": model_cis["top1_accuracy"] if model_cis else float("nan"),
        "model_pmax90_accuracy_genome_cis": model_cis["pmax_threshold_accuracy"] if model_cis else float("nan"),
        "model_pmax90_recall_genome_cis": model_cis["pmax_threshold_recall"] if model_cis else float("nan"),
        "model_top1_accuracy_genome_trans": model_trans["top1_accuracy"] if model_trans else float("nan"),
        "model_pmax90_accuracy_genome_trans": model_trans["pmax_threshold_accuracy"] if model_trans else float("nan"),
        "model_pmax90_recall_genome_trans": model_trans["pmax_threshold_recall"] if model_trans else float("nan"),
        "charm3dg_top1_accuracy_genome_all": charm_all["top1_accuracy"] if charm_all else float("nan"),
        "charm3dg_pmax90_accuracy_genome_all": charm_all["pmax_threshold_accuracy"] if charm_all else float("nan"),
        "charm3dg_pmax90_recall_genome_all": charm_all["pmax_threshold_recall"] if charm_all else float("nan"),
        "charm3dg_top1_accuracy_genome_cis": charm_cis["top1_accuracy"] if charm_cis else float("nan"),
        "charm3dg_pmax90_accuracy_genome_cis": charm_cis["pmax_threshold_accuracy"] if charm_cis else float("nan"),
        "charm3dg_pmax90_recall_genome_cis": charm_cis["pmax_threshold_recall"] if charm_cis else float("nan"),
        "charm3dg_top1_accuracy_genome_trans": charm_trans["top1_accuracy"] if charm_trans else float("nan"),
        "charm3dg_pmax90_accuracy_genome_trans": charm_trans["pmax_threshold_accuracy"] if charm_trans else float("nan"),
        "charm3dg_pmax90_recall_genome_trans": charm_trans["pmax_threshold_recall"] if charm_trans else float("nan"),
        "phase_used_eval_only": 1,
        "charm_3dg_used_eval_only": 1,
    }
    summary.update(pair_stats)
    for key in ["n_raw", "n_bpair", "n_beads", "resolution", "bin_size_bp", "mstep_graph_mode", "input_contact_source", "uses_charm_or_reference", "uses_phase_labels", "copy_labels_are_gauge_only"]:
        if key in train_manifest:
            summary[f"train_manifest_{key}"] = train_manifest[key]

    write_table(args.outdir / "cis_distance_correlations.tsv", cis_rows)
    write_table(args.outdir / "contact_accuracy.tsv", accuracy_rows)
    write_table(args.outdir / "copy_separation.tsv", separation_rows)
    write_table(args.outdir / "per_chrom_volume.tsv", volume_rows_data)
    metrics_path = args.outdir / "metrics.tsv"
    if metrics_path.exists():
        metrics_path.unlink()
    write_summary(args.outdir / "summary.tsv", summary)
    plot_chr1_distance_maps(plots_dir / "chr1_distance_maps.png", reference, reconstruction, distance_swaps)
    plot_3d_scatter(plots_dir / "all_chrom_3d_scatter.png", reference, reconstruction, distance_swaps)
    plot_3d_scatter(plots_dir / "chr1_copy_3d_scatter.png", reference, reconstruction, distance_swaps, chrom="chr1")
    write_readme(args.outdir / "README.md", args, summary, cis_rows, accuracy_rows, separation_rows, volume_rows_data)
    with (args.outdir / "eval_manifest.json").open("w") as fh:
        json.dump(
            {
                "summary": summary,
                "cis_distance_correlations": cis_rows,
                "contact_accuracy": accuracy_rows,
                "copy_separation": separation_rows,
                "per_chrom_volume": volume_rows_data,
            },
            fh,
            indent=2,
        )

    print(f"wrote\t{args.outdir / 'README.md'}")
    print(f"wrote\t{args.outdir / 'cis_distance_correlations.tsv'}")
    print(f"wrote\t{args.outdir / 'contact_accuracy.tsv'}")
    print(f"wrote\t{args.outdir / 'copy_separation.tsv'}")
    print(f"wrote\t{args.outdir / 'per_chrom_volume.tsv'}")
    print(f"wrote\t{plots_dir / 'chr1_distance_maps.png'}")
    print(f"wrote\t{plots_dir / 'all_chrom_3d_scatter.png'}")
    print(f"wrote\t{plots_dir / 'chr1_copy_3d_scatter.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
