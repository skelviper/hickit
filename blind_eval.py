#!/usr/bin/env python3
"""Eval-only helpers for blind diploid P9016 coordinate diagnostics.

This module intentionally has no training entry points.  CHARM/3DG haplotype
labels are reference-only evaluation labels and must not be imported by C
training runners or used for model selection inside a reconstruction run.
"""

from __future__ import annotations

import gzip
import hashlib
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


MODEL_COORD_COLUMNS = ["chr", "start", "end", "bid", "copy", "diploid_bid", "x", "y", "z"]
REF_COLUMNS = ["chr", "hap", "start", "x", "y", "z"]
HAPS = ("mat", "pat")
COPIES = (0, 1)

Permutation = Tuple[Tuple[int, str], Tuple[int, str]]
PERMUTATIONS: Tuple[Permutation, Permutation] = (
    ((0, "mat"), (1, "pat")),
    ((0, "pat"), (1, "mat")),
)
DISTANCE_STRATA: Tuple[Tuple[str, int, Optional[int]], ...] = (
    ("1_2Mb", 1_000_000, 2_000_000),
    ("2_5Mb", 2_000_000, 5_000_000),
    ("5_10Mb", 5_000_000, 10_000_000),
    ("10_20Mb", 10_000_000, 20_000_000),
    ("20_50Mb", 20_000_000, 50_000_000),
    ("50_100Mb", 50_000_000, 100_000_000),
    ("gt100Mb", 100_000_000, None),
)


def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _as_path(path: str | Path) -> Path:
    return path if isinstance(path, Path) else Path(path)


def load_model_coords(path: str | Path) -> pd.DataFrame:
    """Load model coordinates from Hickit blind TSV/TSV.GZ output."""
    path = _as_path(path)
    frame = pd.read_csv(path, sep="\t")
    missing = [col for col in MODEL_COORD_COLUMNS if col not in frame.columns]
    if missing:
        raise ValueError(f"{path} missing model coordinate columns: {missing}")
    frame = frame[MODEL_COORD_COLUMNS].copy()
    frame["chr"] = frame["chr"].astype(str)
    for col in ["start", "end", "bid", "copy", "diploid_bid"]:
        frame[col] = frame[col].astype(int)
    for col in ["x", "y", "z"]:
        frame[col] = frame[col].astype(float)
    bad_copy = sorted(set(frame["copy"]) - set(COPIES))
    if bad_copy:
        raise ValueError(f"{path} has non-diploid copy labels: {bad_copy}")
    if frame[["x", "y", "z"]].isna().any().any() or not np.isfinite(frame[["x", "y", "z"]].to_numpy()).all():
        raise ValueError(f"{path} has non-finite coordinates")
    dup = frame.duplicated(["chr", "copy", "start"], keep=False)
    if dup.any():
        examples = frame.loc[dup, ["chr", "copy", "start"]].drop_duplicates().head().to_dict("records")
        raise ValueError(f"{path} has duplicate model bins: {examples}")
    return frame.sort_values(["chr", "start", "copy"], kind="mergesort").reset_index(drop=True)


def load_charm_3dg(path: str | Path) -> pd.DataFrame:
    """Load CHARM 3DG coordinates with first column formatted as chr(hap)."""
    path = _as_path(path)
    records: List[Tuple[str, str, int, float, float, float]] = []
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                raise ValueError(f"{path} line {line_number}: expected at least 5 columns, got {len(fields)}")
            chrom_hap = fields[0]
            if not chrom_hap.endswith(")") or "(" not in chrom_hap:
                raise ValueError(f"{path} line {line_number}: expected chr(hap), got {chrom_hap!r}")
            chrom, hap = chrom_hap[:-1].split("(", 1)
            records.append((chrom, hap, int(fields[1]), float(fields[2]), float(fields[3]), float(fields[4])))
    frame = pd.DataFrame(records, columns=REF_COLUMNS)
    if frame.empty:
        raise ValueError(f"{path} contains no 3DG records")
    unexpected = sorted(set(frame["hap"]) - set(HAPS))
    if unexpected:
        raise ValueError(f"{path} has unexpected reference hap labels: {unexpected}")
    if frame[["x", "y", "z"]].isna().any().any() or not np.isfinite(frame[["x", "y", "z"]].to_numpy()).all():
        raise ValueError(f"{path} has non-finite reference coordinates")
    dup = frame.duplicated(["chr", "hap", "start"], keep=False)
    if dup.any():
        examples = frame.loc[dup, ["chr", "hap", "start"]].drop_duplicates().head().to_dict("records")
        raise ValueError(f"{path} has duplicate reference bins: {examples}")
    return frame.sort_values(["chr", "start", "hap"], kind="mergesort").reset_index(drop=True)


def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    if int(finite.sum()) < 3:
        return math.nan
    x = x[finite] - np.mean(x[finite])
    y = y[finite] - np.mean(y[finite])
    denom = math.sqrt(float(np.dot(x, x) * np.dot(y, y)))
    if denom == 0.0:
        return math.nan
    return float(np.dot(x, y) / denom)


def rankdata_average(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.shape[0], dtype=np.float64)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + j - 1) + 1.0
        i = j
    return ranks


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    if int(finite.sum()) < 3:
        return math.nan
    return pearson_corr(rankdata_average(x[finite]), rankdata_average(y[finite]))


def distance_vector(coords: np.ndarray) -> np.ndarray:
    coords = np.asarray(coords, dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"expected Nx3 coordinates, got {coords.shape}")
    if coords.shape[0] < 2:
        return np.empty(0, dtype=np.float64)
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2, dtype=np.float64))
    return dist[np.triu_indices(coords.shape[0], k=1)]


def distance_vector_with_delta(starts: Sequence[int], coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    starts_arr = np.asarray(starts, dtype=np.int64)
    coords = np.asarray(coords, dtype=np.float64)
    if coords.shape[0] < 2:
        return np.empty(0), np.empty(0, dtype=np.int64)
    upper = np.triu_indices(coords.shape[0], k=1)
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2, dtype=np.float64))[upper]
    delta = np.abs(starts_arr[upper[1]] - starts_arr[upper[0]])
    return dist, delta


def distance_regression(model_dist: np.ndarray, ref_dist: np.ndarray) -> Tuple[float, float]:
    model_dist = np.asarray(model_dist, dtype=np.float64)
    ref_dist = np.asarray(ref_dist, dtype=np.float64)
    finite = np.isfinite(model_dist) & np.isfinite(ref_dist)
    if int(finite.sum()) < 3:
        return math.nan, math.nan
    x = ref_dist[finite]
    y = model_dist[finite]
    x_center = x - x.mean()
    denom = float(np.dot(x_center, x_center))
    if denom == 0.0:
        return math.nan, math.nan
    slope = float(np.dot(x_center, y - y.mean()) / denom)
    intercept = float(y.mean() - slope * x.mean())
    return slope, intercept


def distance_ratio_quantiles(model_dist: np.ndarray, ref_dist: np.ndarray) -> Dict[str, float]:
    model_dist = np.asarray(model_dist, dtype=np.float64)
    ref_dist = np.asarray(ref_dist, dtype=np.float64)
    finite = np.isfinite(model_dist) & np.isfinite(ref_dist) & (ref_dist > 0)
    out = {f"ratio_q{q}": math.nan for q in [10, 25, 50, 75, 90]}
    out["median_distance_ratio"] = math.nan
    if int(finite.sum()) == 0:
        return out
    ratio = model_dist[finite] / ref_dist[finite]
    for q in [10, 25, 50, 75, 90]:
        out[f"ratio_q{q}"] = float(np.quantile(ratio, q / 100.0))
    out["median_distance_ratio"] = out["ratio_q50"]
    return out


def radius_of_gyration(coords: np.ndarray) -> float:
    coords = np.asarray(coords, dtype=np.float64)
    if coords.ndim != 2 or coords.shape[0] == 0:
        return math.nan
    centered = coords - coords.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum(centered * centered, axis=1))))


def centroid(coords: np.ndarray) -> np.ndarray:
    coords = np.asarray(coords, dtype=np.float64)
    if coords.ndim != 2 or coords.shape[0] == 0:
        return np.array([math.nan, math.nan, math.nan], dtype=np.float64)
    return coords.mean(axis=0)


def procrustes_rmsd(model_coords: np.ndarray, ref_coords: np.ndarray) -> float:
    """Similarity Procrustes RMSD after translation, rotation and uniform scale."""
    x = np.asarray(model_coords, dtype=np.float64)
    y = np.asarray(ref_coords, dtype=np.float64)
    finite = np.isfinite(x).all(axis=1) & np.isfinite(y).all(axis=1)
    x = x[finite]
    y = y[finite]
    if x.shape[0] < 3:
        return math.nan
    x0 = x - x.mean(axis=0)
    y0 = y - y.mean(axis=0)
    denom = float(np.sum(x0 * x0))
    if denom == 0.0:
        return math.nan
    u, s, vt = np.linalg.svd(x0.T @ y0, full_matrices=False)
    sign = np.ones_like(s)
    r = u @ vt
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1.0
        sign[-1] = -1.0
        r = u @ vt
    scale = float(np.sum(s * sign) / denom)
    aligned = scale * x0 @ r
    diff = aligned - y0
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def permutation_name(permutation: Permutation) -> str:
    return "_".join(f"copy{copy}_{hap}" for copy, hap in permutation)


def inverse_permutation(permutation: Permutation) -> Dict[int, str]:
    return {copy: hap for copy, hap in permutation}


def chromosome_sort_key(chrom: str) -> Tuple[int, str]:
    if chrom.startswith("chr"):
        tail = chrom[3:]
    else:
        tail = chrom
    if tail.isdigit():
        return int(tail), ""
    special = {"X": 1000, "Y": 1001, "M": 1002, "MT": 1002}
    return special.get(tail, 2000), tail


def stable_seed(*parts: object) -> int:
    text = "|".join(str(part) for part in parts)
    digest = hashlib.sha256(text.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "little", signed=False) & 0x7FFFFFFF


@dataclass
class CopyDistance:
    copy: int
    hap: str
    chrom: str
    starts: List[int]
    model_coords: np.ndarray
    ref_coords: np.ndarray
    model_dist: np.ndarray
    ref_dist: np.ndarray
    genomic_delta: np.ndarray
    pearson: float
    spearman: float
    slope: float
    intercept: float


def _coords_by_key(frame: pd.DataFrame, key_cols: Sequence[str]) -> Dict[Tuple[object, ...], np.ndarray]:
    return {
        tuple(getattr(row, col) for col in key_cols): np.array([row.x, row.y, row.z], dtype=np.float64)
        for row in frame.itertuples(index=False)
    }


def matched_copy_distance(model_df: pd.DataFrame, ref_df: pd.DataFrame, chrom: str, copy: int, hap: str) -> Optional[CopyDistance]:
    m = model_df.loc[(model_df["chr"] == chrom) & (model_df["copy"] == copy), ["start", "x", "y", "z"]]
    r = ref_df.loc[(ref_df["chr"] == chrom) & (ref_df["hap"] == hap), ["start", "x", "y", "z"]]
    m_by = {int(row.start): np.array([row.x, row.y, row.z], dtype=np.float64) for row in m.itertuples(index=False)}
    r_by = {int(row.start): np.array([row.x, row.y, row.z], dtype=np.float64) for row in r.itertuples(index=False)}
    starts = sorted(set(m_by).intersection(r_by))
    if len(starts) < 3:
        return None
    model_coords = np.vstack([m_by[start] for start in starts])
    ref_coords = np.vstack([r_by[start] for start in starts])
    model_dist, delta = distance_vector_with_delta(starts, model_coords)
    ref_dist = distance_vector(ref_coords)
    slope, intercept = distance_regression(model_dist, ref_dist)
    return CopyDistance(
        copy=copy,
        hap=hap,
        chrom=chrom,
        starts=starts,
        model_coords=model_coords,
        ref_coords=ref_coords,
        model_dist=model_dist,
        ref_dist=ref_dist,
        genomic_delta=delta,
        pearson=pearson_corr(model_dist, ref_dist),
        spearman=spearman_corr(model_dist, ref_dist),
        slope=slope,
        intercept=intercept,
    )


def evaluate_chromosome_permutation(model_df: pd.DataFrame, ref_df: pd.DataFrame, chrom: str, permutation: Permutation) -> Dict[str, object]:
    copy_results: List[CopyDistance] = []
    for copy, hap in permutation:
        result = matched_copy_distance(model_df, ref_df, chrom, copy, hap)
        if result is None:
            return {
                "chrom": chrom,
                "mapping": permutation_name(permutation),
                "valid": False,
                "n_common_bins_min": 0,
                "n_distance_pairs_total": 0,
                "min_copy_pearson": math.nan,
                "mean_copy_pearson": math.nan,
                "pooled_pearson": math.nan,
                "min_copy_spearman": math.nan,
                "mean_copy_spearman": math.nan,
                "pooled_spearman": math.nan,
            }
        copy_results.append(result)
    pearsons = np.array([r.pearson for r in copy_results], dtype=np.float64)
    spearmans = np.array([r.spearman for r in copy_results], dtype=np.float64)
    pooled_model = np.concatenate([r.model_dist for r in copy_results])
    pooled_ref = np.concatenate([r.ref_dist for r in copy_results])
    slopes = np.array([r.slope for r in copy_results], dtype=np.float64)
    ratios = distance_ratio_quantiles(pooled_model, pooled_ref)
    return {
        "chrom": chrom,
        "mapping": permutation_name(permutation),
        "valid": bool(np.isfinite(pearsons).all() and np.isfinite(spearmans).all()),
        "n_common_bins_min": min(len(r.starts) for r in copy_results),
        "n_distance_pairs_total": int(sum(r.model_dist.size for r in copy_results)),
        "min_copy_pearson": float(np.min(pearsons)) if np.isfinite(pearsons).all() else math.nan,
        "mean_copy_pearson": float(np.mean(pearsons)) if np.isfinite(pearsons).all() else math.nan,
        "pooled_pearson": pearson_corr(pooled_model, pooled_ref),
        "min_copy_spearman": float(np.min(spearmans)) if np.isfinite(spearmans).all() else math.nan,
        "mean_copy_spearman": float(np.mean(spearmans)) if np.isfinite(spearmans).all() else math.nan,
        "pooled_spearman": spearman_corr(pooled_model, pooled_ref),
        "slope_median": float(np.nanmedian(slopes)),
        "intercept_median": float(np.nanmedian([r.intercept for r in copy_results])),
        **ratios,
        "copy_results": copy_results,
    }


def best_chromosome_mapping(model_df: pd.DataFrame, ref_df: pd.DataFrame, chrom: str, metric: str = "min_copy_pearson") -> Dict[str, object]:
    evaluated = [evaluate_chromosome_permutation(model_df, ref_df, chrom, p) for p in PERMUTATIONS]
    valid = [row for row in evaluated if row.get("valid")]
    if not valid:
        return evaluated[0]
    best = max(valid, key=lambda row: (float(row.get(metric, math.nan)), float(row.get("min_copy_spearman", math.nan))))
    # Copy identity margins use the four one-copy correlations.
    r00 = matched_copy_distance(model_df, ref_df, chrom, 0, "mat")
    r01 = matched_copy_distance(model_df, ref_df, chrom, 0, "pat")
    r10 = matched_copy_distance(model_df, ref_df, chrom, 1, "mat")
    r11 = matched_copy_distance(model_df, ref_df, chrom, 1, "pat")
    if all(r is not None for r in [r00, r01, r10, r11]):
        fit_a = min(r00.pearson, r11.pearson)
        fit_b = min(r01.pearson, r10.pearson)
        margin_a = min(r00.pearson - r01.pearson, r11.pearson - r10.pearson)
        margin_b = min(r01.pearson - r00.pearson, r10.pearson - r11.pearson)
        sfit_a = min(r00.spearman, r11.spearman)
        sfit_b = min(r01.spearman, r10.spearman)
        smargin_a = min(r00.spearman - r01.spearman, r11.spearman - r10.spearman)
        smargin_b = min(r01.spearman - r00.spearman, r10.spearman - r11.spearman)
        use_a = best["mapping"] == "copy0_mat_copy1_pat"
        best.update({
            "copy_fit_min_pearson": fit_a if use_a else fit_b,
            "copy_identity_margin_min_pearson": margin_a if use_a else margin_b,
            "copy_contrast_score_pearson": (fit_a if use_a else fit_b) * max(margin_a if use_a else margin_b, 0.0),
            "copy_fit_min_spearman": sfit_a if use_a else sfit_b,
            "copy_identity_margin_min_spearman": smargin_a if use_a else smargin_b,
            "copy_contrast_score_spearman": (sfit_a if use_a else sfit_b) * max(smargin_a if use_a else smargin_b, 0.0),
            "copy0_mat_pearson": r00.pearson,
            "copy0_pat_pearson": r01.pearson,
            "copy1_mat_pearson": r10.pearson,
            "copy1_pat_pearson": r11.pearson,
        })
    ref_mat = matched_copy_distance(
        ref_df.rename(columns={"hap": "copy"}).assign(copy=lambda x: x["copy"].map({"mat": 0, "pat": 1})),
        ref_df,
        chrom,
        0,
        "pat",
    )
    if ref_mat is not None:
        best["ref_mat_pat_corr"] = ref_mat.pearson
        best["ref_separability"] = 1.0 - ref_mat.pearson if math.isfinite(ref_mat.pearson) else math.nan
    return best


def chromosome_list(model_df: pd.DataFrame, ref_df: pd.DataFrame) -> List[str]:
    return sorted(set(model_df["chr"]).intersection(set(ref_df["chr"])), key=chromosome_sort_key)


def mapping_from_name(name: str) -> Dict[int, str]:
    if name == "copy0_mat_copy1_pat":
        return {0: "mat", 1: "pat"}
    if name == "copy0_pat_copy1_mat":
        return {0: "pat", 1: "mat"}
    raise ValueError(f"unknown mapping {name}")


def mapped_point_table(model_df: pd.DataFrame, ref_df: pd.DataFrame, per_chr_mapping: Mapping[str, str]) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for chrom, mapping_name in per_chr_mapping.items():
        mapping = mapping_from_name(mapping_name)
        for copy, hap in mapping.items():
            m = model_df.loc[
                (model_df["chr"] == chrom) & (model_df["copy"] == copy),
                ["start", "x", "y", "z"],
            ].rename(columns={"x": "model_x", "y": "model_y", "z": "model_z"})
            r = ref_df.loc[
                (ref_df["chr"] == chrom) & (ref_df["hap"] == hap),
                ["start", "x", "y", "z"],
            ].rename(columns={"x": "ref_x", "y": "ref_y", "z": "ref_z"})
            merged = m.merge(r, on="start", how="inner", sort=True)
            if merged.empty:
                continue
            for item in merged.itertuples(index=False):
                rows.append({
                    "chr": chrom,
                    "start": int(item.start),
                    "copy": copy,
                    "hap": hap,
                    "model_x": float(item.model_x),
                    "model_y": float(item.model_y),
                    "model_z": float(item.model_z),
                    "ref_x": float(item.ref_x),
                    "ref_y": float(item.ref_y),
                    "ref_z": float(item.ref_z),
                })
    return pd.DataFrame(rows)


def centroid_rows(points: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row_key, group in points.groupby(["chr", "copy", "hap"], sort=False):
        chrom, copy, hap = row_key
        rows.append({
            "chr": chrom,
            "copy": int(copy),
            "hap": hap,
            "model_x": float(group["model_x"].mean()),
            "model_y": float(group["model_y"].mean()),
            "model_z": float(group["model_z"].mean()),
            "ref_x": float(group["ref_x"].mean()),
            "ref_y": float(group["ref_y"].mean()),
            "ref_z": float(group["ref_z"].mean()),
            "n_bins": int(len(group)),
        })
    return pd.DataFrame(rows)


def distance_metrics_from_points(points: pd.DataFrame) -> Dict[str, float]:
    if len(points) < 3:
        return {"pearson": math.nan, "spearman": math.nan, "slope": math.nan, "median_distance_ratio": math.nan}
    model = points[["model_x", "model_y", "model_z"]].to_numpy(float)
    ref = points[["ref_x", "ref_y", "ref_z"]].to_numpy(float)
    md = distance_vector(model)
    rd = distance_vector(ref)
    slope, _ = distance_regression(md, rd)
    ratios = distance_ratio_quantiles(md, rd)
    return {
        "pearson": pearson_corr(md, rd),
        "spearman": spearman_corr(md, rd),
        "slope": slope,
        "median_distance_ratio": ratios["median_distance_ratio"],
    }


def deterministic_pair_indices(n_a: int, n_b: int, max_pairs: int, seed: int) -> Tuple[np.ndarray, np.ndarray]:
    total = n_a * n_b
    if total == 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    if total <= max_pairs:
        flat = np.arange(total, dtype=np.int64)
    else:
        rng = np.random.default_rng(seed)
        flat = rng.choice(total, size=max_pairs, replace=False)
    return flat // n_b, flat % n_b


def sampled_cross_distances(
    left_model: np.ndarray,
    right_model: np.ndarray,
    left_ref: np.ndarray,
    right_ref: np.ndarray,
    max_pairs: int,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray]:
    ia, ib = deterministic_pair_indices(left_model.shape[0], right_model.shape[0], max_pairs, seed)
    if ia.size == 0:
        return np.empty(0), np.empty(0)
    model_diff = left_model[ia] - right_model[ib]
    ref_diff = left_ref[ia] - right_ref[ib]
    return (
        np.sqrt(np.sum(model_diff * model_diff, axis=1)),
        np.sqrt(np.sum(ref_diff * ref_diff, axis=1)),
    )
