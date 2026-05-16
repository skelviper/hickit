#!/usr/bin/env python3
"""Eval-only whole-genome diagnostics for blind P9016 diploid coordinates.

This runner reads finished model coordinate files and an external CHARM 3DG
reference.  It never invokes training and never writes back into model outputs.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from blind_eval import (
        COPIES,
        DISTANCE_STRATA,
        HAPS,
        PERMUTATIONS,
        chromosome_list,
        chromosome_sort_key,
        centroid_rows,
        deterministic_pair_indices,
        distance_metrics_from_points,
        distance_ratio_quantiles,
        distance_regression,
        distance_vector,
        distance_vector_with_delta,
        evaluate_chromosome_permutation,
        load_charm_3dg,
        load_model_coords,
        mapped_point_table,
        matched_copy_distance,
        pearson_corr,
        permutation_name,
        procrustes_rmsd,
        radius_of_gyration,
        sampled_cross_distances,
        spearman_corr,
        stable_seed,
    )
except ImportError:  # pragma: no cover - supports package-style imports.
    from hickit.blind_eval import (  # type: ignore
        COPIES,
        DISTANCE_STRATA,
        HAPS,
        PERMUTATIONS,
        chromosome_list,
        chromosome_sort_key,
        centroid_rows,
        deterministic_pair_indices,
        distance_metrics_from_points,
        distance_ratio_quantiles,
        distance_regression,
        distance_vector,
        distance_vector_with_delta,
        evaluate_chromosome_permutation,
        load_charm_3dg,
        load_model_coords,
        mapped_point_table,
        matched_copy_distance,
        pearson_corr,
        permutation_name,
        procrustes_rmsd,
        radius_of_gyration,
        sampled_cross_distances,
        spearman_corr,
        stable_seed,
    )


DEFAULT_ROOT = Path("/tmp/hk_blind_p9016_rep1_fine_n_rs_grid_135846_1777436403_0")
DEFAULT_SUMMARY = DEFAULT_ROOT / "scan_summary.tsv"
DEFAULT_MATRIX_SUMMARY = DEFAULT_ROOT / "matrix_summary.tsv"
DEFAULT_TDG = Path("/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz")
LOG4 = math.log(4.0)
ABS_RATIO_MIN = 0.7
ABS_RATIO_MAX = 1.3
TRANS_NULL_MARGIN = 0.05
TRANS_CIS_MIN_PEARSON = 0.30
NULL_PERMUTATIONS = 20
NULL_BOOTSTRAPS = 30
NULL_MAX_POINTS = 20_000


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def finite_mean(values: Sequence[float]) -> float:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return float(arr.mean()) if arr.size else math.nan


def finite_median(values: Sequence[float]) -> float:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return float(np.median(arr)) if arr.size else math.nan


def fraction_lt(values: Sequence[float], threshold: float) -> float:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return float(np.mean(arr < threshold)) if arr.size else math.nan


def finite_values(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    return arr[np.isfinite(arr)]


def metric_pass(value: object, threshold: object) -> int:
    val = safe_float(value)
    thr = safe_float(threshold)
    return int(math.isfinite(val) and math.isfinite(thr) and val > thr)


def in_range(value: object, low: float, high: float) -> int:
    val = safe_float(value)
    return int(math.isfinite(val) and low <= val <= high)


def ratio_not_compressed(value: object, low: float = ABS_RATIO_MIN) -> int:
    val = safe_float(value)
    return int(math.isfinite(val) and val >= low)


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(list(rows))
    frame.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)


def existing_path(value: object) -> Path | None:
    if value is None:
        return None
    text = str(value)
    if not text or text == "nan":
        return None
    path = Path(text)
    return path if path.exists() else None


def coords_path_from_summary(summary_row: pd.Series) -> Path:
    for column in ["coords_gz", "output_coords", "coords", "coords_path"]:
        if column in summary_row.index:
            path = existing_path(summary_row[column])
            if path is not None:
                return path
    if "output_dir" in summary_row.index:
        out_dir = existing_path(summary_row["output_dir"])
        if out_dir is not None:
            for name in ["p9016_final.coords.tsv.gz", "p9016_full.coords.tsv", "p9016_full.coords.tsv.gz"]:
                path = out_dir / name
                if path.exists():
                    return path
    raise ValueError("summary row has no readable coords path column")


def posterior_path_from_summary(summary_row: pd.Series) -> Path | None:
    for column in ["output_bpair_posterior", "bpair_posterior", "posterior_path"]:
        if column in summary_row.index:
            path = existing_path(summary_row[column])
            if path is not None:
                return path
    if "output_dir" in summary_row.index:
        out_dir = existing_path(summary_row["output_dir"])
        if out_dir is not None:
            path = out_dir / "p9016_full.bpair_posterior.tsv"
            if path.exists():
                return path
    return None


def force_diag_path_from_summary(summary_row: pd.Series) -> Path | None:
    for column in ["output_force_class_diag", "force_class_diag", "force_diag_path"]:
        if column in summary_row.index:
            path = existing_path(summary_row[column])
            if path is not None:
                return path
    if "output_dir" in summary_row.index:
        out_dir = existing_path(summary_row["output_dir"])
        if out_dir is not None:
            for name in ["p9016_full.force_class_diag.tsv", "eval.force_diag.tsv"]:
                path = out_dir / name
                if path.exists():
                    return path
    return None


def read_force_class_diag(path: Path | None) -> Dict[str, object]:
    if path is None:
        return {"force_class_diag_available": 0}
    frame = pd.read_csv(path, sep="\t")
    if "class" not in frame.columns:
        return {"force_class_diag_available": 0}
    by_class = {str(row["class"]): row for _, row in frame.iterrows()}
    cis = by_class.get("cis", {})
    trans = by_class.get("trans", {})
    total = by_class.get("total", {})
    cis_force = safe_float(cis.get("contact_force_l1", math.nan))
    trans_force = safe_float(trans.get("contact_force_l1", math.nan))
    repulsion_force = safe_float(total.get("repulsion_force_l1", math.nan))
    return {
        "force_class_diag_available": 1,
        "final_contact_energy_cis": safe_float(cis.get("contact_energy", math.nan)),
        "final_contact_energy_trans": safe_float(trans.get("contact_energy", math.nan)),
        "final_contact_force_l1_cis": cis_force,
        "final_contact_force_l1_trans": trans_force,
        "final_backbone_force_l1": safe_float(total.get("backbone_force_l1", math.nan)),
        "final_repulsion_force_l1": repulsion_force,
        "final_homolog_sep_force_l1": safe_float(total.get("homolog_sep_force_l1", math.nan)),
        "final_contact_wedges_cis": safe_float(cis.get("n_wedges", math.nan)),
        "final_contact_wedges_trans": safe_float(trans.get("n_wedges", math.nan)),
        "force_l1_cis_trans_ratio": cis_force / trans_force if math.isfinite(cis_force) and math.isfinite(trans_force) and trans_force > 0.0 else math.nan,
        "repulsion_to_trans_force_ratio": repulsion_force / trans_force if math.isfinite(repulsion_force) and math.isfinite(trans_force) and trans_force > 0.0 else math.nan,
    }


def best_mapping_for_chrom(model_df: pd.DataFrame, ref_df: pd.DataFrame, chrom: str) -> Dict[str, object]:
    rows = [evaluate_chromosome_permutation(model_df, ref_df, chrom, perm) for perm in PERMUTATIONS]
    valid = [row for row in rows if bool(row.get("valid"))]
    if not valid:
        row = rows[0]
        row["best_by"] = "invalid"
        return row
    row = max(valid, key=lambda x: (safe_float(x.get("min_copy_pearson")), safe_float(x.get("min_copy_spearman"))))
    row = dict(row)
    row["best_by"] = "min_copy_pearson"

    r00 = matched_copy_distance(model_df, ref_df, chrom, 0, "mat")
    r01 = matched_copy_distance(model_df, ref_df, chrom, 0, "pat")
    r10 = matched_copy_distance(model_df, ref_df, chrom, 1, "mat")
    r11 = matched_copy_distance(model_df, ref_df, chrom, 1, "pat")
    if all(x is not None for x in [r00, r01, r10, r11]):
        assert r00 is not None and r01 is not None and r10 is not None and r11 is not None
        use_a = row["mapping"] == "copy0_mat_copy1_pat"
        fit_a = min(r00.pearson, r11.pearson)
        fit_b = min(r01.pearson, r10.pearson)
        margin_a = min(r00.pearson - r01.pearson, r11.pearson - r10.pearson)
        margin_b = min(r01.pearson - r00.pearson, r10.pearson - r11.pearson)
        sfit_a = min(r00.spearman, r11.spearman)
        sfit_b = min(r01.spearman, r10.spearman)
        smargin_a = min(r00.spearman - r01.spearman, r11.spearman - r10.spearman)
        smargin_b = min(r01.spearman - r00.spearman, r10.spearman - r11.spearman)
        row.update({
            "copy0_mat_pearson": r00.pearson,
            "copy0_pat_pearson": r01.pearson,
            "copy1_mat_pearson": r10.pearson,
            "copy1_pat_pearson": r11.pearson,
            "copy0_mat_spearman": r00.spearman,
            "copy0_pat_spearman": r01.spearman,
            "copy1_mat_spearman": r10.spearman,
            "copy1_pat_spearman": r11.spearman,
            "copy_fit_min_pearson": fit_a if use_a else fit_b,
            "copy_identity_margin_min_pearson": margin_a if use_a else margin_b,
            "copy_contrast_score_pearson": (fit_a if use_a else fit_b) * max(margin_a if use_a else margin_b, 0.0),
            "copy_fit_min_spearman": sfit_a if use_a else sfit_b,
            "copy_identity_margin_min_spearman": smargin_a if use_a else smargin_b,
            "copy_contrast_score_spearman": (sfit_a if use_a else sfit_b) * max(smargin_a if use_a else smargin_b, 0.0),
        })
    ref_mat = matched_copy_distance(
        ref_df.assign(copy=ref_df["hap"].map({"mat": 0, "pat": 1})),
        ref_df,
        chrom,
        0,
        "pat",
    )
    if ref_mat is not None:
        row["ref_mat_pat_corr"] = ref_mat.pearson
        row["ref_separability"] = 1.0 - ref_mat.pearson if math.isfinite(ref_mat.pearson) else math.nan
    return row


def mapping_copy_results(row: Mapping[str, object]) -> List[object]:
    return [x for x in row.get("copy_results", []) if x is not None]


def cis_strata_rows(config_name: str, best_rows: Sequence[Mapping[str, object]]) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    for row in best_rows:
        chrom = str(row["chrom"])
        mapping = str(row["mapping"])
        for copy_result in mapping_copy_results(row):
            for stratum, low, high in DISTANCE_STRATA:
                delta = copy_result.genomic_delta
                mask = delta >= low
                if high is not None:
                    mask &= delta < high
                model = copy_result.model_dist[mask]
                ref = copy_result.ref_dist[mask]
                slope, intercept = distance_regression(model, ref)
                ratios = distance_ratio_quantiles(model, ref)
                out.append({
                    "config_name": config_name,
                    "chr": chrom,
                    "mapping": mapping,
                    "copy": copy_result.copy,
                    "hap": copy_result.hap,
                    "stratum": stratum,
                    "min_delta_bp": low,
                    "max_delta_bp": high if high is not None else "",
                    "n_pairs": int(model.size),
                    "pearson": pearson_corr(model, ref),
                    "spearman": spearman_corr(model, ref),
                    "slope": slope,
                    "intercept": intercept,
                    **ratios,
                })
    return out


def rg_and_scale_rows(config_name: str, best_rows: Sequence[Mapping[str, object]]) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], List[Dict[str, object]]]:
    rg_rows: List[Dict[str, object]] = []
    quantile_rows: List[Dict[str, object]] = []
    scale_rows: List[Dict[str, object]] = []
    for row in best_rows:
        chrom = str(row["chrom"])
        mapping = str(row["mapping"])
        for copy_result in mapping_copy_results(row):
            model_rg = radius_of_gyration(copy_result.model_coords)
            ref_rg = radius_of_gyration(copy_result.ref_coords)
            ratio = model_rg / ref_rg if ref_rg > 0 else math.nan
            rg_rows.append({
                "config_name": config_name,
                "chr": chrom,
                "mapping": mapping,
                "copy": copy_result.copy,
                "hap": copy_result.hap,
                "n_bins": len(copy_result.starts),
                "rg_model": model_rg,
                "rg_ref": ref_rg,
                "rg_ratio": ratio,
            })
            qs = {}
            for q in [5, 10, 25, 50, 75, 90, 95]:
                if copy_result.model_dist.size and copy_result.ref_dist.size:
                    qs[f"model_q{q}"] = float(np.quantile(copy_result.model_dist, q / 100.0))
                    qs[f"ref_q{q}"] = float(np.quantile(copy_result.ref_dist, q / 100.0))
                    qs[f"ratio_q{q}"] = qs[f"model_q{q}"] / qs[f"ref_q{q}"] if qs[f"ref_q{q}"] > 0 else math.nan
                else:
                    qs[f"model_q{q}"] = math.nan
                    qs[f"ref_q{q}"] = math.nan
                    qs[f"ratio_q{q}"] = math.nan
            quantile_rows.append({
                "config_name": config_name,
                "chr": chrom,
                "mapping": mapping,
                "copy": copy_result.copy,
                "hap": copy_result.hap,
                "n_pairs": int(copy_result.model_dist.size),
                **qs,
            })
            slope, intercept = distance_regression(copy_result.model_dist, copy_result.ref_dist)
            scale_rows.append({
                "config_name": config_name,
                "chr": chrom,
                "mapping": mapping,
                "copy": copy_result.copy,
                "hap": copy_result.hap,
                "n_pairs": int(copy_result.model_dist.size),
                "pearson": copy_result.pearson,
                "spearman": copy_result.spearman,
                "slope": slope,
                "intercept": intercept,
                **distance_ratio_quantiles(copy_result.model_dist, copy_result.ref_dist),
            })
    return rg_rows, quantile_rows, scale_rows


def points_for_mapping(model_df: pd.DataFrame, ref_df: pd.DataFrame, mapping_by_chr: Mapping[str, str]) -> pd.DataFrame:
    return mapped_point_table(model_df, ref_df, mapping_by_chr)


def global_centroid_mapping(
    model_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    chroms: Sequence[str],
    initial_mapping: Mapping[str, str] | None = None,
) -> Tuple[Dict[str, str], Dict[str, object]]:
    # Brute force 2^N is exact for small fixtures.  For P9016 (~20 chromosomes)
    # 2^N per config is too costly across a full grid, so use deterministic
    # coordinate descent while recording the mode in eval.global_gauge.tsv.
    n = len(chroms)
    if n == 0:
        return {}, {"global_mapping_mode": "none", "centroid_distance_spearman": math.nan, "centroid_distance_pearson": math.nan}

    centroid_cache: Dict[Tuple[str, str], pd.DataFrame] = {}
    for chrom in chroms:
        for mapping_name in ["copy0_mat_copy1_pat", "copy0_pat_copy1_mat"]:
            centroid_cache[(chrom, mapping_name)] = centroid_rows(points_for_mapping(model_df, ref_df, {chrom: mapping_name}))

    best_mapping: Dict[str, str] = {}
    best_metrics: Dict[str, object] = {}
    best_key = (-math.inf, -math.inf)

    def score_mapping(mapping: Mapping[str, str]) -> Tuple[Tuple[float, float], Dict[str, object]]:
        frames = [centroid_cache[(chrom, mapping[chrom])] for chrom in chroms if (chrom, mapping[chrom]) in centroid_cache]
        centroids = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        metrics = distance_metrics_from_points(centroids)
        return (safe_float(metrics["spearman"]), safe_float(metrics["pearson"])), metrics

    if n <= 12:
        for mask in range(1 << n):
            mapping = {
                chrom: ("copy0_pat_copy1_mat" if (mask >> i) & 1 else "copy0_mat_copy1_pat")
                for i, chrom in enumerate(chroms)
            }
            key, metrics = score_mapping(mapping)
            if key > best_key:
                best_key = key
                best_mapping = dict(mapping)
                best_metrics = metrics
        mode = "bruteforce_centroid"
    else:
        # Start from independent per-chromosome cis-best mappings, then
        # greedily flip one chromosome at a time until centroid score stops
        # improving.  This preserves deterministic behavior without reference
        # metrics entering any training path.
        mapping = {chrom: str((initial_mapping or {}).get(chrom, "copy0_mat_copy1_pat")) for chrom in chroms}
        best_key, best_metrics = score_mapping(mapping)
        changed = True
        while changed:
            changed = False
            for chrom in chroms:
                trial = dict(mapping)
                trial[chrom] = "copy0_pat_copy1_mat" if trial[chrom] == "copy0_mat_copy1_pat" else "copy0_mat_copy1_pat"
                key, metrics = score_mapping(trial)
                if key > best_key:
                    mapping = trial
                    best_key = key
                    best_metrics = metrics
                    changed = True
        best_mapping = dict(mapping)
        mode = "coordinate_descent_centroid"
    return best_mapping, {
        "global_mapping_mode": mode,
        "centroid_distance_pearson": best_metrics.get("pearson", math.nan),
        "centroid_distance_spearman": best_metrics.get("spearman", math.nan),
        "centroid_distance_slope": best_metrics.get("slope", math.nan),
        "centroid_median_distance_ratio": best_metrics.get("median_distance_ratio", math.nan),
        "centroid_relative_metric": best_metrics.get("spearman", math.nan),
        "centroid_absolute_scale_metric": best_metrics.get("median_distance_ratio", math.nan),
    }


def trans_metrics(
    config_name: str,
    points: pd.DataFrame,
    max_pairs: int,
    sample_seed: int,
) -> Tuple[List[Dict[str, object]], Dict[str, object], Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    pooled_model: List[np.ndarray] = []
    pooled_ref: List[np.ndarray] = []
    pooled_model_centroid: List[np.ndarray] = []
    pooled_ref_centroid: List[np.ndarray] = []
    chroms = sorted(points["chr"].unique(), key=chromosome_sort_key)
    for i, chr_a in enumerate(chroms):
        for chr_b in chroms[i + 1:]:
            for copy_a in COPIES:
                for copy_b in COPIES:
                    left = points.loc[(points["chr"] == chr_a) & (points["copy"] == copy_a)]
                    right = points.loc[(points["chr"] == chr_b) & (points["copy"] == copy_b)]
                    if left.empty or right.empty:
                        continue
                    md, rd = sampled_cross_distances(
                        left[["model_x", "model_y", "model_z"]].to_numpy(float),
                        right[["model_x", "model_y", "model_z"]].to_numpy(float),
                        left[["ref_x", "ref_y", "ref_z"]].to_numpy(float),
                        right[["ref_x", "ref_y", "ref_z"]].to_numpy(float),
                        max_pairs=max_pairs,
                        seed=stable_seed(config_name, chr_a, chr_b, copy_a, copy_b, sample_seed),
                    )
                    if md.size == 0:
                        continue
                    model_centroid_distance = float(np.linalg.norm(
                        left[["model_x", "model_y", "model_z"]].to_numpy(float).mean(axis=0) -
                        right[["model_x", "model_y", "model_z"]].to_numpy(float).mean(axis=0)
                    ))
                    ref_centroid_distance = float(np.linalg.norm(
                        left[["ref_x", "ref_y", "ref_z"]].to_numpy(float).mean(axis=0) -
                        right[["ref_x", "ref_y", "ref_z"]].to_numpy(float).mean(axis=0)
                    ))
                    pooled_model.append(md)
                    pooled_ref.append(rd)
                    pooled_model_centroid.append(np.full(md.shape, model_centroid_distance, dtype=np.float64))
                    pooled_ref_centroid.append(np.full(rd.shape, ref_centroid_distance, dtype=np.float64))
                    metrics = distance_metric_summary(md, rd, "")
                    rows.append({
                        "config_name": config_name,
                        "chr_a": chr_a,
                        "chr_b": chr_b,
                        "copy_a": copy_a,
                        "copy_b": copy_b,
                        "n_pairs": int(md.size),
                        "model_centroid_distance": model_centroid_distance,
                        "ref_centroid_distance": ref_centroid_distance,
                        **metrics,
                    })
    if pooled_model:
        model = np.concatenate(pooled_model)
        ref = np.concatenate(pooled_ref)
        model_centroid = np.concatenate(pooled_model_centroid)
        ref_centroid = np.concatenate(pooled_ref_centroid)
    else:
        model = np.empty(0)
        ref = np.empty(0)
        model_centroid = np.empty(0)
        ref_centroid = np.empty(0)
    metrics = distance_metric_summary(model, ref, "trans_")
    model_resid = residualize_against_covariate(model, model_centroid)
    ref_resid = residualize_against_covariate(ref, ref_centroid)
    summary = {
        "config_name": config_name,
        "trans_n_pairs": int(model.size),
        "trans_distance_pearson": metrics["trans_pearson"],
        "trans_distance_spearman": metrics["trans_spearman"],
        "trans_distance_slope": metrics["trans_slope"],
        "trans_distance_intercept": metrics["trans_intercept"],
        "trans_median_distance_ratio": metrics["trans_median_distance_ratio"],
        "trans_relative_metric": metrics["trans_spearman"],
        "trans_absolute_scale_metric": metrics["trans_median_distance_ratio"],
        "trans_centroid_residual_pearson": pearson_corr(model_resid, ref_resid),
        "trans_centroid_residual_spearman": spearman_corr(model_resid, ref_resid),
        **metrics,
    }
    return rows, summary, {"trans_model_dist": model, "trans_ref_dist": ref}


def trans_chrpair_mean_rows(config_name: str, trans_rows: Sequence[Mapping[str, object]]) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    if not trans_rows:
        return [], {
            "config_name": config_name,
            "trans_chrpair_mean_spearman": math.nan,
            "trans_chrpair_mean_pearson": math.nan,
            "trans_chrpair_mean_slope": math.nan,
            "trans_chrpair_mean_median_distance_ratio": math.nan,
        }
    frame = pd.DataFrame(list(trans_rows))
    rows: List[Dict[str, object]] = []
    for key, group in frame.groupby(["chr_a", "chr_b"], sort=False):
        chr_a, chr_b = key
        row = {
            "config_name": config_name,
            "chr_a": chr_a,
            "chr_b": chr_b,
            "n_copy_pairs": int(len(group)),
            "n_pairs_total": int(group["n_pairs"].sum()),
            "pearson_mean": finite_mean(group["pearson"].to_numpy(float)),
            "spearman_mean": finite_mean(group["spearman"].to_numpy(float)),
            "slope_mean": finite_mean(group["slope"].to_numpy(float)),
            "median_distance_ratio_mean": finite_mean(group["median_distance_ratio"].to_numpy(float)),
            "abs_log_ratio_median_mean": finite_mean(group["abs_log_ratio_median"].to_numpy(float)),
            "residual_mad_mean": finite_mean(group["residual_mad"].to_numpy(float)),
        }
        rows.append(row)
    return rows, {
        "config_name": config_name,
        "trans_chrpair_mean_spearman": finite_mean([safe_float(r["spearman_mean"]) for r in rows]),
        "trans_chrpair_mean_pearson": finite_mean([safe_float(r["pearson_mean"]) for r in rows]),
        "trans_chrpair_mean_slope": finite_mean([safe_float(r["slope_mean"]) for r in rows]),
        "trans_chrpair_mean_median_distance_ratio": finite_mean([safe_float(r["median_distance_ratio_mean"]) for r in rows]),
        "trans_chrpair_mean_abs_log_ratio_median": finite_mean([safe_float(r["abs_log_ratio_median_mean"]) for r in rows]),
        "trans_chrpair_mean_residual_mad": finite_mean([safe_float(r["residual_mad_mean"]) for r in rows]),
    }


def replace_model_coords(points: pd.DataFrame, coords: np.ndarray) -> pd.DataFrame:
    out = points.copy()
    out[["model_x", "model_y", "model_z"]] = coords
    return out


def transformed_baseline_points(points: pd.DataFrame, baseline_name: str, config_name: str, sample_seed: int) -> pd.DataFrame:
    points = points.reset_index(drop=True).copy()
    coords = points[["model_x", "model_y", "model_z"]].to_numpy(float)
    center = coords.mean(axis=0)
    scale = float(np.std(coords))
    if not math.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    rng = np.random.default_rng(stable_seed(config_name, baseline_name, sample_seed))

    if baseline_name == "random_diploid":
        return replace_model_coords(points, center + rng.normal(size=coords.shape) * scale)

    if baseline_name == "shuffled_chromosome_labels":
        perm = rng.permutation(len(coords))
        return replace_model_coords(points, coords[perm])

    if baseline_name == "centroid_only":
        out = points.copy()
        for _, idx in out.groupby(["chr", "copy"], sort=False).groups.items():
            centroid = coords[list(idx)].mean(axis=0)
            out.loc[list(idx), ["model_x", "model_y", "model_z"]] = centroid
        return out

    if baseline_name == "compact_ball_per_chromosome":
        out = points.copy()
        for key, idx in out.groupby(["chr", "copy"], sort=False).groups.items():
            group_idx = list(idx)
            group = coords[group_idx]
            centroid = group.mean(axis=0)
            rg = float(np.sqrt(np.mean(np.sum((group - centroid) ** 2, axis=1)))) if len(group_idx) else 0.0
            radius = max(rg * 0.05, scale * 0.01)
            local_rng = np.random.default_rng(stable_seed(config_name, baseline_name, key[0], key[1], sample_seed))
            out.loc[group_idx, ["model_x", "model_y", "model_z"]] = centroid + local_rng.normal(size=(len(group_idx), 3)) * radius
        return out

    if baseline_name == "per_chromosome_translated_scaffold":
        out = points.copy()
        for key, idx in out.groupby(["chr", "copy"], sort=False).groups.items():
            group_idx = list(idx)
            local_rng = np.random.default_rng(stable_seed(config_name, baseline_name, key[0], key[1], sample_seed))
            shift = local_rng.normal(size=3) * scale * 2.0
            out.loc[group_idx, ["model_x", "model_y", "model_z"]] = coords[group_idx] + shift
        return out

    raise ValueError(f"unknown baseline {baseline_name}")


def baseline_diagnostics(
    config_name: str,
    points: pd.DataFrame,
    max_pairs: int,
    sample_seed: int,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    baseline_names = [
        "random_diploid",
        "shuffled_chromosome_labels",
        "centroid_only",
        "compact_ball_per_chromosome",
        "per_chromosome_translated_scaffold",
    ]
    rows: List[Dict[str, object]] = []
    summary: Dict[str, object] = {"config_name": config_name}
    for name in baseline_names:
        bpoints = transformed_baseline_points(points, name, config_name, sample_seed)
        _, trans_summary, _ = trans_metrics(f"{config_name}:{name}", bpoints, max_pairs, sample_seed)
        centroid_metrics = distance_metrics_from_points(centroid_rows(bpoints))
        row = {
            "config_name": config_name,
            "baseline_name": name,
            "is_baseline": 1,
            "trans_relative_metric": trans_summary.get("trans_relative_metric", math.nan),
            "trans_distance_spearman": trans_summary.get("trans_distance_spearman", math.nan),
            "trans_distance_slope": trans_summary.get("trans_distance_slope", math.nan),
            "trans_median_distance_ratio": trans_summary.get("trans_median_distance_ratio", math.nan),
            "trans_abs_log_ratio_median": trans_summary.get("trans_abs_log_ratio_median", math.nan),
            "trans_centroid_residual_spearman": trans_summary.get("trans_centroid_residual_spearman", math.nan),
            "centroid_relative_metric": centroid_metrics.get("spearman", math.nan),
            "centroid_distance_spearman": centroid_metrics.get("spearman", math.nan),
            "centroid_distance_slope": centroid_metrics.get("slope", math.nan),
            "centroid_median_distance_ratio": centroid_metrics.get("median_distance_ratio", math.nan),
        }
        rows.append(row)
        prefix = f"baseline_{name}_"
        for key, value in row.items():
            if key not in {"config_name", "baseline_name", "is_baseline"}:
                summary[prefix + key] = value
    return rows, summary


def sampled_trans_null_diagnostics(
    config_name: str,
    trans_payload: Mapping[str, object],
    trans_chrpair_rows_in: Sequence[Mapping[str, object]],
    sample_seed: int,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    model = np.asarray(trans_payload.get("trans_model_dist", np.empty(0)), dtype=np.float64)
    ref = np.asarray(trans_payload.get("trans_ref_dist", np.empty(0)), dtype=np.float64)
    mask = np.isfinite(model) & np.isfinite(ref)
    model = model[mask]
    ref = ref[mask]
    rng = np.random.default_rng(stable_seed(config_name, "sampled_trans_null", sample_seed))
    if model.size > NULL_MAX_POINTS:
        idx = rng.choice(model.size, size=NULL_MAX_POINTS, replace=False)
        model = model[idx]
        ref = ref[idx]

    rows: List[Dict[str, object]] = []
    perm_values: List[float] = []
    if model.size >= 3:
        for i in range(NULL_PERMUTATIONS):
            shuffled = ref[rng.permutation(ref.size)]
            value = spearman_corr(model, shuffled)
            perm_values.append(value)
            rows.append({
                "config_name": config_name,
                "null_type": "sampled_trans_ref_distance_permutation",
                "iteration": i,
                "metric": "spearman",
                "value": value,
            })

    boot_values: List[float] = []
    if model.size >= 3:
        for i in range(NULL_BOOTSTRAPS):
            idx = rng.integers(0, model.size, size=model.size)
            value = spearman_corr(model[idx], ref[idx])
            boot_values.append(value)
            rows.append({
                "config_name": config_name,
                "null_type": "sampled_trans_pair_bootstrap",
                "iteration": i,
                "metric": "spearman",
                "value": value,
            })

    chrpair_values: List[float] = [safe_float(r.get("spearman_mean")) for r in trans_chrpair_rows_in]
    chrpair_values = [v for v in chrpair_values if math.isfinite(v)]
    chrpair_boot: List[float] = []
    if chrpair_values:
        arr = np.asarray(chrpair_values, dtype=np.float64)
        for i in range(NULL_BOOTSTRAPS):
            value = float(np.mean(arr[rng.integers(0, arr.size, size=arr.size)]))
            chrpair_boot.append(value)
            rows.append({
                "config_name": config_name,
                "null_type": "chromosome_pair_mean_bootstrap",
                "iteration": i,
                "metric": "spearman_mean",
                "value": value,
            })

    def quantile(values: Sequence[float], q: float) -> float:
        arr = np.asarray([v for v in values if math.isfinite(v)], dtype=np.float64)
        return float(np.quantile(arr, q)) if arr.size else math.nan

    summary = {
        "config_name": config_name,
        "trans_null_permutation_spearman_p50": quantile(perm_values, 0.50),
        "trans_null_permutation_spearman_p95": quantile(perm_values, 0.95),
        "trans_bootstrap_spearman_q025": quantile(boot_values, 0.025),
        "trans_bootstrap_spearman_q975": quantile(boot_values, 0.975),
        "trans_chrpair_bootstrap_spearman_mean_q025": quantile(chrpair_boot, 0.025),
        "trans_chrpair_bootstrap_spearman_mean_q975": quantile(chrpair_boot, 0.975),
        "trans_null_n_points": int(model.size),
        "trans_null_n_permutations": NULL_PERMUTATIONS,
        "trans_null_n_bootstraps": NULL_BOOTSTRAPS,
    }
    return rows, summary


def rho_train_for_bpair(mode: str, schedule_value: float, rho_output: float, contact_class: str, floor: float) -> float:
    rho_conf = max(0.0, min(1.0, rho_output))
    if not math.isfinite(schedule_value):
        schedule_value = 1.0
    if not math.isfinite(floor) or floor <= 0.0:
        floor = 0.5
    if mode == "constant":
        return schedule_value
    if mode == "entropy":
        return schedule_value * rho_conf
    if mode == "entropy_with_floor":
        return schedule_value * max(rho_conf, floor)
    if mode == "entropy_cis_constant_trans":
        return schedule_value if contact_class == "trans" else schedule_value * rho_conf
    if mode == "entropy_cis_floor_trans":
        return schedule_value * max(rho_conf, floor) if contact_class == "trans" else schedule_value * rho_conf
    return math.nan


def training_state_weight_sum(
    mode: str,
    contact_class: str,
    p: np.ndarray,
    pmax: float,
    margin: float,
    trans_margin_min: float,
    trans_pmax_min: float,
    gamma: float,
) -> float:
    if contact_class != "trans" or mode in {"", "posterior"}:
        return 1.0
    if mode == "top_if_confident":
        if not math.isfinite(trans_margin_min):
            trans_margin_min = 0.0
        if not math.isfinite(trans_pmax_min):
            trans_pmax_min = 0.0
        return 1.0 if margin >= trans_margin_min and pmax >= trans_pmax_min else 0.0
    if mode == "power_sharpen":
        if not math.isfinite(gamma) or gamma <= 0.0:
            gamma = 1.0
        total = float(np.sum(np.power(np.clip(p, 0.0, 1.0), gamma)))
        return 1.0 if total > 0.0 and math.isfinite(total) else 0.0
    return 1.0


def high_count_short_distance_precision(frame: pd.DataFrame, fraction: float, distance_column: str = "expected_distance") -> float:
    if frame.empty or "n_raw" not in frame.columns or distance_column not in frame.columns:
        return math.nan
    n = max(1, int(math.ceil(len(frame) * fraction)))
    high_count = set(frame.sort_values("n_raw", ascending=False).head(n).index)
    short_distance = set(frame.sort_values(distance_column, ascending=True).head(n).index)
    return len(high_count & short_distance) / float(n) if n > 0 else math.nan


def fraction_le(values: Sequence[float], threshold: float) -> float:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return float(np.mean(arr <= threshold)) if arr.size else math.nan


def distance_residual_summary(model: np.ndarray, ref: np.ndarray, prefix: str = "") -> Dict[str, object]:
    model = finite_values(model)
    ref = finite_values(ref)
    if model.size != ref.size:
        n = min(model.size, ref.size)
        model = model[:n]
        ref = ref[:n]
    if model.size == 0:
        return {
            f"{prefix}residual_rmse": math.nan,
            f"{prefix}residual_mae": math.nan,
            f"{prefix}residual_median_abs": math.nan,
            f"{prefix}residual_mad": math.nan,
            f"{prefix}abs_log_ratio_median": math.nan,
            f"{prefix}abs_log_ratio_q90": math.nan,
            f"{prefix}fraction_ratio_lt_0_7": math.nan,
            f"{prefix}fraction_ratio_gt_1_3": math.nan,
            f"{prefix}fraction_abs_log_ratio_gt_log_1_3": math.nan,
        }
    slope, intercept = distance_regression(model, ref)
    fitted = slope * ref + intercept if math.isfinite(slope) and math.isfinite(intercept) else np.full_like(model, math.nan)
    residual = model - fitted
    abs_residual = np.abs(residual[np.isfinite(residual)])
    positive = (model > 0.0) & (ref > 0.0)
    ratio = model[positive] / ref[positive]
    log_abs = np.abs(np.log(ratio)) if ratio.size else np.empty(0)
    med_abs = float(np.median(abs_residual)) if abs_residual.size else math.nan
    mad = float(np.median(np.abs(abs_residual - med_abs))) if abs_residual.size and math.isfinite(med_abs) else math.nan
    return {
        f"{prefix}residual_rmse": float(math.sqrt(np.mean(residual[np.isfinite(residual)] ** 2))) if np.isfinite(residual).any() else math.nan,
        f"{prefix}residual_mae": float(np.mean(abs_residual)) if abs_residual.size else math.nan,
        f"{prefix}residual_median_abs": med_abs,
        f"{prefix}residual_mad": mad,
        f"{prefix}abs_log_ratio_median": float(np.median(log_abs)) if log_abs.size else math.nan,
        f"{prefix}abs_log_ratio_q90": float(np.quantile(log_abs, 0.9)) if log_abs.size else math.nan,
        f"{prefix}fraction_ratio_lt_0_7": float(np.mean(ratio < ABS_RATIO_MIN)) if ratio.size else math.nan,
        f"{prefix}fraction_ratio_gt_1_3": float(np.mean(ratio > ABS_RATIO_MAX)) if ratio.size else math.nan,
        f"{prefix}fraction_abs_log_ratio_gt_log_1_3": float(np.mean(log_abs > math.log(ABS_RATIO_MAX))) if log_abs.size else math.nan,
    }


def distance_metric_summary(model: np.ndarray, ref: np.ndarray, prefix: str = "") -> Dict[str, object]:
    slope, intercept = distance_regression(model, ref)
    ratios = distance_ratio_quantiles(model, ref)
    return {
        f"{prefix}pearson": pearson_corr(model, ref),
        f"{prefix}spearman": spearman_corr(model, ref),
        f"{prefix}slope": slope,
        f"{prefix}intercept": intercept,
        f"{prefix}median_distance_ratio": ratios["median_distance_ratio"],
        **distance_residual_summary(model, ref, prefix),
    }


def quantile_edges(values: Sequence[float], quantiles: Sequence[float]) -> List[float]:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return []
    edges = [float(np.quantile(arr, q)) for q in quantiles]
    edges = sorted(set(edges))
    if len(edges) < 2:
        value = edges[0]
        return [value, np.nextafter(value, math.inf)]
    edges[-1] = np.nextafter(edges[-1], math.inf)
    return edges


def finite_fixed_edges(values: Sequence[float], candidate_edges: Sequence[float]) -> List[float]:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return []
    lo = float(np.min(arr))
    hi = float(np.max(arr))
    edges = [x for x in candidate_edges if math.isfinite(x) and lo <= x <= hi]
    edges = sorted(set([lo, *edges, np.nextafter(hi, math.inf)]))
    if len(edges) < 2:
        return [lo, np.nextafter(lo, math.inf)]
    return edges


def bin_label(left: float, right: float) -> str:
    right_text = f"{right:.6g}"
    if math.isfinite(right):
        previous = np.nextafter(right, -math.inf)
        right_text = f"{previous:.6g}" if previous > left else f"{right:.6g}"
    return f"[{left:.6g},{right_text}]"


def dscale_group_summary(config_name: str, frame: pd.DataFrame) -> Dict[str, object]:
    return {
        "config_name": config_name,
        "available": 1,
        "n_bpair": int(len(frame)),
        "n_raw_sum": float(frame["n_raw"].sum()) if "n_raw" in frame.columns else math.nan,
        "mean_n_raw": float(frame["n_raw"].mean()) if "n_raw" in frame.columns and len(frame) else math.nan,
        "mean_pU": float(frame["pU"].mean()) if "pU" in frame.columns and len(frame) else math.nan,
        "median_margin": float(frame["margin"].median()) if "margin" in frame.columns and len(frame) else math.nan,
        "mean_base_d_scale": float(frame["base_d_scale"].mean()) if "base_d_scale" in frame.columns and len(frame) else math.nan,
        "median_base_d_scale": float(frame["base_d_scale"].median()) if "base_d_scale" in frame.columns and len(frame) else math.nan,
        "median_expected_distance": float(frame["expected_distance"].median()) if "expected_distance" in frame.columns and len(frame) else math.nan,
        "median_top_distance": float(frame["top_distance"].median()) if "top_distance" in frame.columns and len(frame) else math.nan,
        "median_expected_residual": float(frame["expected_residual"].median()) if "expected_residual" in frame.columns and len(frame) else math.nan,
        "median_top_residual": float(frame["top_residual"].median()) if "top_residual" in frame.columns and len(frame) else math.nan,
        "contact_count_vs_inv_expected_distance_spearman": (
            spearman_corr(np.log1p(frame["n_raw"].to_numpy(float)), -frame["expected_distance"].to_numpy(float))
            if {"n_raw", "expected_distance"}.issubset(frame.columns) else math.nan
        ),
    }


def dscale_stratified_diagnostics(
    config_name: str,
    diag: pd.DataFrame,
    config_dir: Path,
    dscale_source_column: str | None,
) -> Dict[str, object]:
    if dscale_source_column is None:
        row = {"config_name": config_name, "available": 0, "reason": "missing_base_d_scale_or_d_scale"}
        write_tsv(config_dir / "eval.dscale_stratified.tsv", [row])
        return {"dscale_stratified_available": 0, "dscale_stratified_reason": row["reason"]}
    if diag.empty:
        row = {"config_name": config_name, "available": 0, "reason": "no_bpair_rows_with_coords"}
        write_tsv(config_dir / "eval.dscale_stratified.tsv", [row])
        return {"dscale_stratified_available": 0, "dscale_stratified_reason": row["reason"]}
    finite_d = finite_values(diag["base_d_scale"].to_numpy(float))
    if finite_d.size == 0:
        row = {"config_name": config_name, "available": 0, "reason": f"no_finite_{dscale_source_column}"}
        write_tsv(config_dir / "eval.dscale_stratified.tsv", [row])
        return {"dscale_stratified_available": 0, "dscale_stratified_reason": row["reason"]}

    rows: List[Dict[str, object]] = []

    def append_row(stratum_type: str, stratum_value: object, frame: pd.DataFrame, contact_class: str = "all", left: object = "", right: object = "") -> None:
        if frame.empty:
            return
        rows.append({
            "stratum_type": stratum_type,
            "stratum_value": stratum_value,
            "contact_class": contact_class,
            "bin_left": left,
            "bin_right": right,
            "d_scale_source_column": dscale_source_column,
            **dscale_group_summary(config_name, frame),
        })

    append_row("all", "all", diag)
    for contact_class, group in diag.groupby("contact_class", sort=False):
        append_row("contact_class", contact_class, group, str(contact_class))

    bin_specs = [
        ("n_raw_bin", "n_raw", quantile_edges(diag["n_raw"].to_numpy(float), [0.0, 0.25, 0.5, 0.75, 0.9, 1.0])),
        ("base_d_scale_bin", "base_d_scale", quantile_edges(diag["base_d_scale"].to_numpy(float), [0.0, 0.25, 0.5, 0.75, 0.9, 1.0])),
        ("pU_bin", "pU", finite_fixed_edges(diag["pU"].to_numpy(float), [0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0])),
        ("margin_bin", "margin", finite_fixed_edges(diag["margin"].to_numpy(float), [0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0])),
    ]
    for stratum_type, column, edges in bin_specs:
        if len(edges) < 2 or column not in diag.columns:
            continue
        values = diag[column].to_numpy(float)
        finite_mask = np.isfinite(values)
        for left, right in zip(edges[:-1], edges[1:]):
            mask = finite_mask & (values >= left) & (values < right)
            group = diag.loc[mask]
            label = bin_label(float(left), float(right))
            append_row(stratum_type, label, group, "all", float(left), float(right))
            for contact_class, class_group in group.groupby("contact_class", sort=False):
                append_row(stratum_type, label, class_group, str(contact_class), float(left), float(right))

    if not rows:
        row = {"config_name": config_name, "available": 0, "reason": "no_nonempty_strata"}
        write_tsv(config_dir / "eval.dscale_stratified.tsv", [row])
        return {"dscale_stratified_available": 0, "dscale_stratified_reason": row["reason"]}
    write_tsv(config_dir / "eval.dscale_stratified.tsv", rows)
    return {
        "dscale_stratified_available": 1,
        "dscale_stratified_reason": "",
        "dscale_stratified_n_rows": len(rows),
        "dscale_source_column": dscale_source_column,
    }


def residualize_against_covariate(values: np.ndarray, covariate: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    covariate = np.asarray(covariate, dtype=np.float64)
    if values.size != covariate.size or values.size == 0:
        return np.empty(0)
    mask = np.isfinite(values) & np.isfinite(covariate)
    out = np.full(values.shape, math.nan, dtype=np.float64)
    if mask.sum() < 3:
        return out
    slope, intercept = distance_regression(values[mask], covariate[mask])
    if not math.isfinite(slope) or not math.isfinite(intercept):
        return out
    out[mask] = values[mask] - (slope * covariate[mask] + intercept)
    return out


def posterior_contact_diagnostics(
    config_name: str,
    model_df: pd.DataFrame,
    posterior_path: Path | None,
    summary_row: pd.Series,
    config_dir: Path,
) -> Dict[str, object]:
    if posterior_path is None:
        missing = {"config_name": config_name, "available": 0, "reason": "bpair_posterior_not_found"}
        write_tsv(config_dir / "eval.posterior_cis_trans_summary.tsv", [missing])
        write_tsv(config_dir / "eval.posterior_by_chrpair.tsv", [missing])
        write_tsv(config_dir / "eval.posterior_pu_hist.tsv", [missing])
        write_tsv(config_dir / "eval.contact_distance_summary.tsv", [missing])
        write_tsv(config_dir / "eval.trans_contact_by_chrpair.tsv", [missing])
        dscale_diag = dscale_stratified_diagnostics(config_name, pd.DataFrame(), config_dir, None)
        return {"posterior_diag_available": 0, **dscale_diag}

    posterior = pd.read_csv(posterior_path, sep="\t")
    dscale_source_column = next((col for col in ["base_d_scale", "d_scale"] if col in posterior.columns), None)
    required = ["chr1", "start1", "chr2", "start2", "n_raw", "base_k", "p00", "p01", "p10", "p11", "pU"]
    missing_cols = [col for col in required if col not in posterior.columns]
    if missing_cols:
        missing = {"config_name": config_name, "available": 0, "reason": f"missing_columns:{','.join(missing_cols)}"}
        write_tsv(config_dir / "eval.posterior_cis_trans_summary.tsv", [missing])
        write_tsv(config_dir / "eval.contact_distance_summary.tsv", [missing])
        dscale_diag = dscale_stratified_diagnostics(config_name, pd.DataFrame(), config_dir, dscale_source_column)
        return {"posterior_diag_available": 0, **dscale_diag}

    coords = {
        (str(row.chr), int(row.start), int(row.copy)): np.array([row.x, row.y, row.z], dtype=np.float64)
        for row in model_df.itertuples(index=False)
    }
    rho_mode = str(summary_row.get("rho_train_mode", "constant"))
    rho_schedule = safe_float(summary_row.get("rho_train_end", summary_row.get("final_rho_train", 1.0)))
    rho_floor = safe_float(summary_row.get("rho_train_floor", summary_row.get("rho_floor", 0.5)))
    cis_multiplier = safe_float(summary_row.get("contact_k_multiplier_cis", 1.0))
    trans_multiplier = safe_float(summary_row.get("contact_k_multiplier_trans", 1.0))
    state_weight_mode = str(summary_row.get("state_weight_mode", "posterior"))
    trans_margin_min = safe_float(summary_row.get("trans_margin_min", 0.0))
    trans_pmax_min = safe_float(summary_row.get("trans_pmax_min", 0.0))
    trans_posterior_power_gamma = safe_float(summary_row.get("trans_posterior_power_gamma", 1.0))
    if not math.isfinite(cis_multiplier) or cis_multiplier == 0.0:
        cis_multiplier = 1.0
    if not math.isfinite(trans_multiplier) or trans_multiplier == 0.0:
        trans_multiplier = 1.0

    rows: List[Dict[str, object]] = []
    n_missing_coords = 0
    for item in posterior.itertuples(index=False):
        chr1 = str(getattr(item, "chr1"))
        chr2 = str(getattr(item, "chr2"))
        start1 = int(getattr(item, "start1"))
        start2 = int(getattr(item, "start2"))
        xyz10 = coords.get((chr1, start1, 0))
        xyz11 = coords.get((chr1, start1, 1))
        xyz20 = coords.get((chr2, start2, 0))
        xyz21 = coords.get((chr2, start2, 1))
        if xyz10 is None or xyz11 is None or xyz20 is None or xyz21 is None:
            n_missing_coords += 1
            continue
        d00 = float(np.linalg.norm(xyz10 - xyz20))
        d01 = float(np.linalg.norm(xyz10 - xyz21))
        d10 = float(np.linalg.norm(xyz11 - xyz20))
        d11 = float(np.linalg.norm(xyz11 - xyz21))
        p = np.array([float(getattr(item, "p00")), float(getattr(item, "p01")), float(getattr(item, "p10")), float(getattr(item, "p11"))], dtype=np.float64)
        d = np.array([d00, d01, d10, d11], dtype=np.float64)
        top = int(np.argmax(p))
        class_from_chr = "cis" if chr1 == chr2 else "trans"
        contact_class = str(getattr(item, "contact_class", class_from_chr))
        if contact_class not in {"cis", "trans"}:
            contact_class = class_from_chr
        rho_output = float(getattr(item, "rho_output", max(0.0, min(1.0, 1.0 - float(getattr(item, "entropy", LOG4)) / LOG4))))
        rho = rho_train_for_bpair(rho_mode, rho_schedule, rho_output, contact_class, rho_floor)
        multiplier = trans_multiplier if contact_class == "trans" else cis_multiplier
        base_k = float(getattr(item, "base_k"))
        base_d = float(getattr(item, dscale_source_column)) if dscale_source_column is not None else math.nan
        pmax = float(getattr(item, "pmax", float(np.max(p))))
        margin = float(getattr(item, "margin", math.nan))
        if not math.isfinite(margin):
            sorted_p = np.sort(p)
            margin = float(sorted_p[-1] - sorted_p[-2]) if sorted_p.size >= 2 else math.nan
        state_weight_sum = training_state_weight_sum(
            state_weight_mode, contact_class, p, pmax, margin,
            trans_margin_min, trans_pmax_min, trans_posterior_power_gamma,
        )
        eff_k = base_k * rho * multiplier * state_weight_sum if math.isfinite(rho) else math.nan
        e_dist = float(np.dot(p, d))
        rows.append({
            "config_name": config_name,
            "chr1": chr1,
            "chr2": chr2,
            "chrpair": f"{chr1}:{chr2}" if chr1 <= chr2 else f"{chr2}:{chr1}",
            "contact_class": contact_class,
            "n_raw": float(getattr(item, "n_raw")),
            "pU": float(getattr(item, "pU")),
            "entropy": float(getattr(item, "entropy", math.nan)),
            "pmax": pmax,
            "margin": margin,
            "psame": float(p[0] + p[3]),
            "pcross": float(p[1] + p[2]),
            "rho_output": rho_output,
            "rho_train_bpair": rho,
            "base_k": base_k,
            "base_d_scale": base_d,
            "state_weight_mode": state_weight_mode,
            "training_state_weight_sum": state_weight_sum,
            "total_effective_k": eff_k,
            "d00": d00,
            "d01": d01,
            "d10": d10,
            "d11": d11,
            "expected_distance": e_dist,
            "top_distance": float(d[top]),
            "min_distance": float(np.min(d)),
            "top_state": ["00", "01", "10", "11"][top],
            "expected_residual": e_dist / base_d if base_d > 0.0 else math.nan,
            "top_residual": float(d[top]) / base_d if base_d > 0.0 else math.nan,
            "min_residual": float(np.min(d)) / base_d if base_d > 0.0 else math.nan,
        })

    diag = pd.DataFrame(rows)
    if diag.empty:
        missing = {"config_name": config_name, "available": 0, "reason": "no_bpair_rows_with_coords", "n_missing_coords": n_missing_coords}
        write_tsv(config_dir / "eval.posterior_cis_trans_summary.tsv", [missing])
        write_tsv(config_dir / "eval.contact_distance_summary.tsv", [missing])
        dscale_diag = dscale_stratified_diagnostics(config_name, diag, config_dir, dscale_source_column)
        return {"posterior_diag_available": 0, "posterior_n_missing_coords": n_missing_coords, **dscale_diag}

    def class_summary(label: str, frame: pd.DataFrame) -> Dict[str, object]:
        count_corr_expected = spearman_corr(np.log1p(frame["n_raw"].to_numpy(float)), -frame["expected_distance"].to_numpy(float))
        count_corr_top = spearman_corr(np.log1p(frame["n_raw"].to_numpy(float)), -frame["top_distance"].to_numpy(float))
        count_corr_min = spearman_corr(np.log1p(frame["n_raw"].to_numpy(float)), -frame["min_distance"].to_numpy(float))
        return {
            "config_name": config_name,
            "class": label,
            "available": 1,
            "n_bpair": int(len(frame)),
            "n_raw_sum": float(frame["n_raw"].sum()),
            "mean_pU": float(frame["pU"].mean()),
            "median_pU": float(frame["pU"].median()),
            "mean_pmax": float(frame["pmax"].mean()),
            "median_pmax": float(frame["pmax"].median()),
            "mean_margin": float(frame["margin"].mean()),
            "median_margin": float(frame["margin"].median()),
            "mean_rho_train_bpair": float(frame["rho_train_bpair"].mean()),
            "median_rho_train_bpair": float(frame["rho_train_bpair"].median()),
            "fraction_contacts_with_rho_train_near_zero": fraction_le(frame["rho_train_bpair"].to_numpy(float), 1e-4),
            "sum_effective_k": float(frame["total_effective_k"].sum()),
            "mean_psame": float(frame["psame"].mean()),
            "median_psame": float(frame["psame"].median()),
            "mean_pcross": float(frame["pcross"].mean()),
            "median_pcross": float(frame["pcross"].median()),
            "contact_count_vs_inv_expected_distance_spearman": count_corr_expected,
            "contact_count_vs_inv_top_distance_spearman": count_corr_top,
            "contact_count_vs_inv_min_distance_spearman": count_corr_min,
            "median_expected_residual": float(frame["expected_residual"].median()),
            "median_top_residual": float(frame["top_residual"].median()),
            "median_min_residual": float(frame["min_residual"].median()),
            "high_count_short_expected_precision_at_1pct": high_count_short_distance_precision(frame, 0.01, "expected_distance"),
            "high_count_short_expected_precision_at_5pct": high_count_short_distance_precision(frame, 0.05, "expected_distance"),
            "high_count_short_expected_precision_at_10pct": high_count_short_distance_precision(frame, 0.10, "expected_distance"),
        }

    posterior_rows = [class_summary(label, group) for label, group in diag.groupby("contact_class", sort=False)]
    by_class = {row["class"]: row for row in posterior_rows}
    cis_k = safe_float(by_class.get("cis", {}).get("sum_effective_k", math.nan))
    trans_k = safe_float(by_class.get("trans", {}).get("sum_effective_k", math.nan))
    total_k = cis_k + trans_k if math.isfinite(cis_k) and math.isfinite(trans_k) else math.nan
    write_tsv(config_dir / "eval.posterior_cis_trans_summary.tsv", posterior_rows)

    chrpair_rows = []
    for key, group in diag.groupby(["contact_class", "chrpair"], sort=False):
        contact_class, chrpair = key
        row = class_summary(str(contact_class), group)
        row["chrpair"] = chrpair
        chrpair_rows.append(row)
    write_tsv(config_dir / "eval.posterior_by_chrpair.tsv", chrpair_rows)
    write_tsv(config_dir / "eval.trans_contact_by_chrpair.tsv", [row for row in chrpair_rows if row.get("class") == "trans"])

    hist_rows = []
    bins = [0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.000001]
    for contact_class, group in diag.groupby("contact_class", sort=False):
        pU = group["pU"].to_numpy(float)
        for left, right in zip(bins[:-1], bins[1:]):
            mask = (pU >= left) & (pU < right)
            hist_rows.append({
                "config_name": config_name,
                "class": contact_class,
                "pU_left": left,
                "pU_right": right,
                "n_bpair": int(mask.sum()),
            })
    write_tsv(config_dir / "eval.posterior_pu_hist.tsv", hist_rows)

    contact_rows = []
    for row in posterior_rows:
        contact_rows.append({
            "config_name": config_name,
            "class": row["class"],
            "available": 1,
            "contact_count_vs_inv_expected_distance_spearman": row["contact_count_vs_inv_expected_distance_spearman"],
            "contact_count_vs_inv_top_distance_spearman": row["contact_count_vs_inv_top_distance_spearman"],
            "contact_count_vs_inv_min_distance_spearman": row["contact_count_vs_inv_min_distance_spearman"],
            "median_expected_residual": row["median_expected_residual"],
            "median_top_residual": row["median_top_residual"],
            "median_min_residual": row["median_min_residual"],
            "high_count_short_expected_precision_at_1pct": row["high_count_short_expected_precision_at_1pct"],
            "high_count_short_expected_precision_at_5pct": row["high_count_short_expected_precision_at_5pct"],
            "high_count_short_expected_precision_at_10pct": row["high_count_short_expected_precision_at_10pct"],
        })
    write_tsv(config_dir / "eval.contact_distance_summary.tsv", contact_rows)
    dscale_diag = dscale_stratified_diagnostics(config_name, diag, config_dir, dscale_source_column)

    return {
        "posterior_diag_available": 1,
        **dscale_diag,
        "posterior_n_bpair_with_coords": int(len(diag)),
        "posterior_n_missing_coords": n_missing_coords,
        "mean_pU_cis": by_class.get("cis", {}).get("mean_pU", math.nan),
        "mean_pU_trans": by_class.get("trans", {}).get("mean_pU", math.nan),
        "median_pU_cis": by_class.get("cis", {}).get("median_pU", math.nan),
        "median_pU_trans": by_class.get("trans", {}).get("median_pU", math.nan),
        "mean_pmax_cis": by_class.get("cis", {}).get("mean_pmax", math.nan),
        "mean_pmax_trans": by_class.get("trans", {}).get("mean_pmax", math.nan),
        "median_pmax_cis": by_class.get("cis", {}).get("median_pmax", math.nan),
        "median_pmax_trans": by_class.get("trans", {}).get("median_pmax", math.nan),
        "mean_margin_cis": by_class.get("cis", {}).get("mean_margin", math.nan),
        "mean_margin_trans": by_class.get("trans", {}).get("mean_margin", math.nan),
        "median_margin_cis": by_class.get("cis", {}).get("median_margin", math.nan),
        "median_margin_trans": by_class.get("trans", {}).get("median_margin", math.nan),
        "mean_psame_cis": by_class.get("cis", {}).get("mean_psame", math.nan),
        "mean_psame_trans": by_class.get("trans", {}).get("mean_psame", math.nan),
        "mean_pcross_cis": by_class.get("cis", {}).get("mean_pcross", math.nan),
        "mean_pcross_trans": by_class.get("trans", {}).get("mean_pcross", math.nan),
        "mean_rho_train_cis": by_class.get("cis", {}).get("mean_rho_train_bpair", math.nan),
        "mean_rho_train_trans": by_class.get("trans", {}).get("mean_rho_train_bpair", math.nan),
        "fraction_trans_contacts_with_rho_train_near_zero": by_class.get("trans", {}).get("fraction_contacts_with_rho_train_near_zero", math.nan),
        "sum_effective_k_cis": cis_k,
        "sum_effective_k_trans": trans_k,
        "effective_k_cis_trans_ratio": cis_k / trans_k if math.isfinite(cis_k) and math.isfinite(trans_k) and trans_k > 0 else math.nan,
        "trans_fraction_effective_k": trans_k / total_k if math.isfinite(total_k) and total_k > 0 else math.nan,
        "cis_contact_count_vs_inv_distance_spearman": by_class.get("cis", {}).get("contact_count_vs_inv_expected_distance_spearman", math.nan),
        "trans_contact_count_vs_inv_distance_spearman": by_class.get("trans", {}).get("contact_count_vs_inv_expected_distance_spearman", math.nan),
        "cis_contact_residual_median": by_class.get("cis", {}).get("median_top_residual", math.nan),
        "trans_contact_residual_median": by_class.get("trans", {}).get("median_top_residual", math.nan),
        "cis_high_count_short_expected_precision_at_5pct": by_class.get("cis", {}).get("high_count_short_expected_precision_at_5pct", math.nan),
        "trans_high_count_short_expected_precision_at_5pct": by_class.get("trans", {}).get("high_count_short_expected_precision_at_5pct", math.nan),
    }


def collect_cis_distance_arrays(best_rows: Sequence[Mapping[str, object]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    pooled_model = []
    pooled_ref = []
    long_model = []
    long_ref = []
    for row in [r for r in best_rows if bool(r.get("valid"))]:
        for cr in mapping_copy_results(row):
            pooled_model.append(cr.model_dist)
            pooled_ref.append(cr.ref_dist)
            mask = cr.genomic_delta >= 20_000_000
            long_model.append(cr.model_dist[mask])
            long_ref.append(cr.ref_dist[mask])
    pm = np.concatenate(pooled_model) if pooled_model else np.empty(0)
    pr = np.concatenate(pooled_ref) if pooled_ref else np.empty(0)
    lm = np.concatenate(long_model) if long_model else np.empty(0)
    lr = np.concatenate(long_ref) if long_ref else np.empty(0)
    return pm, pr, lm, lr


def whole_genome_summary(best_rows: Sequence[Mapping[str, object]], trans_payload: Mapping[str, object]) -> Dict[str, object]:
    cis_model, cis_ref, _, _ = collect_cis_distance_arrays(best_rows)
    trans_model = np.asarray(trans_payload.get("trans_model_dist", np.empty(0)), dtype=np.float64)
    trans_ref = np.asarray(trans_payload.get("trans_ref_dist", np.empty(0)), dtype=np.float64)
    if cis_model.size and trans_model.size:
        all_model = np.concatenate([cis_model, trans_model])
        all_ref = np.concatenate([cis_ref, trans_ref])
    elif cis_model.size:
        all_model, all_ref = cis_model, cis_ref
    else:
        all_model, all_ref = trans_model, trans_ref
    slope, _ = distance_regression(all_model, all_ref)
    ratios = distance_ratio_quantiles(all_model, all_ref)
    return {
        "whole_genome_distance_pearson": pearson_corr(all_model, all_ref),
        "whole_genome_distance_spearman": spearman_corr(all_model, all_ref),
        "whole_genome_distance_slope": slope,
        "whole_genome_median_distance_ratio": ratios["median_distance_ratio"],
        "all_cis_pairs_used": int(cis_model.size),
        "all_trans_pairs_used": int(trans_model.size),
    }


def procrustes_rows(config_name: str, points: pd.DataFrame) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    model_all = points[["model_x", "model_y", "model_z"]].to_numpy(float)
    ref_all = points[["ref_x", "ref_y", "ref_z"]].to_numpy(float)
    global_rmsd = procrustes_rmsd(model_all, ref_all)
    per_chr = []
    per_chr_copy = []
    for chrom, group in points.groupby("chr", sort=False):
        rmsd = procrustes_rmsd(group[["model_x", "model_y", "model_z"]].to_numpy(float), group[["ref_x", "ref_y", "ref_z"]].to_numpy(float))
        per_chr.append(rmsd)
        rows.append({"config_name": config_name, "level": "chr", "chr": chrom, "copy": "", "hap": "", "n_bins": len(group), "rmsd": rmsd})
    for key, group in points.groupby(["chr", "copy", "hap"], sort=False):
        chrom, copy, hap = key
        rmsd = procrustes_rmsd(group[["model_x", "model_y", "model_z"]].to_numpy(float), group[["ref_x", "ref_y", "ref_z"]].to_numpy(float))
        per_chr_copy.append(rmsd)
        rows.append({"config_name": config_name, "level": "chr_copy", "chr": chrom, "copy": copy, "hap": hap, "n_bins": len(group), "rmsd": rmsd})
    summary = {
        "config_name": config_name,
        "procrustes_global_rmsd": global_rmsd,
        "procrustes_per_chr_rmsd_mean": finite_mean(per_chr),
        "procrustes_per_chr_rmsd_median": finite_median(per_chr),
        "procrustes_per_chr_copy_rmsd_mean": finite_mean(per_chr_copy),
        "procrustes_global_over_per_chr_ratio": global_rmsd / finite_mean(per_chr) if finite_mean(per_chr) > 0 else math.nan,
    }
    return rows, summary


def evaluate_config(
    summary_row: pd.Series,
    ref_df: pd.DataFrame,
    out_dir: Path,
    max_trans_pairs: int,
    sample_seed: int,
    baseline_trans_pairs: int,
    run_baselines: bool = True,
) -> Dict[str, object]:
    config_name = str(summary_row["config_name"])
    model_df = load_model_coords(coords_path_from_summary(summary_row))
    chroms = chromosome_list(model_df, ref_df)
    config_dir = out_dir / config_name
    config_dir.mkdir(parents=True, exist_ok=True)
    posterior_diag = posterior_contact_diagnostics(config_name, model_df, posterior_path_from_summary(summary_row), summary_row, config_dir)

    all_perm_rows: List[Dict[str, object]] = []
    best_rows: List[Dict[str, object]] = []
    per_chr_mapping: Dict[str, str] = {}
    for chrom in chroms:
        chrom_perm_rows = [evaluate_chromosome_permutation(model_df, ref_df, chrom, p) for p in PERMUTATIONS]
        for row in chrom_perm_rows:
            out_row = {k: v for k, v in row.items() if k != "copy_results"}
            out_row["config_name"] = config_name
            all_perm_rows.append(out_row)
        best = best_mapping_for_chrom(model_df, ref_df, chrom)
        best_clean = {k: v for k, v in best.items() if k != "copy_results"}
        best_clean["config_name"] = config_name
        best_rows.append(best)
        per_chr_mapping[chrom] = str(best["mapping"])
    write_tsv(config_dir / "eval.cis_per_chr.tsv", [{k: v for k, v in r.items() if k != "copy_results"} | {"config_name": config_name} for r in best_rows])
    write_tsv(config_dir / "eval.cis_per_chr_all_mappings.tsv", all_perm_rows)
    write_tsv(config_dir / "eval.cis_distance_strata.tsv", cis_strata_rows(config_name, best_rows))

    rg_rows, quantile_rows, scale_rows = rg_and_scale_rows(config_name, best_rows)
    write_tsv(config_dir / "eval.rg_by_chr.tsv", rg_rows)
    write_tsv(config_dir / "eval.distance_quantiles_by_chr.tsv", quantile_rows)
    write_tsv(config_dir / "eval.scale_by_chr.tsv", scale_rows)
    write_tsv(config_dir / "eval.compaction_summary.tsv", [compaction_summary(config_name, rg_rows, scale_rows, quantile_rows)])

    global_mapping, centroid_summary = global_centroid_mapping(model_df, ref_df, chroms, per_chr_mapping)
    points = points_for_mapping(model_df, ref_df, global_mapping)
    centroids = centroid_rows(points)
    write_tsv(config_dir / "eval.centroid_matrix.tsv", centroid_pair_rows(config_name, centroids))
    write_tsv(config_dir / "eval.centroid_summary.tsv", [{"config_name": config_name, **centroid_summary}])
    write_tsv(config_dir / "eval.global_gauge.tsv", [
        {"config_name": config_name, "chr": chrom, "global_mapping": global_mapping.get(chrom, ""), "per_chr_cis_mapping": per_chr_mapping.get(chrom, "")}
        for chrom in chroms
    ])

    trans_rows, trans_summary, trans_payload = trans_metrics(config_name, points, max_trans_pairs, sample_seed)
    trans_chrpair_rows_out, trans_chrpair_summary = trans_chrpair_mean_rows(config_name, trans_rows)
    trans_summary.update(trans_chrpair_summary)
    write_tsv(config_dir / "eval.trans_by_chrpair.tsv", trans_rows)
    write_tsv(config_dir / "eval.trans_chrpair_mean.tsv", trans_chrpair_rows_out)
    write_tsv(config_dir / "eval.trans_summary.tsv", [trans_summary])
    if run_baselines:
        baseline_rows, baseline_summary = baseline_diagnostics(config_name, points, baseline_trans_pairs, sample_seed)
    else:
        baseline_rows, baseline_summary = [], {"config_name": config_name, "baseline_diagnostics_available": 0}
    baseline_summary["baseline_diagnostics_available"] = 1 if run_baselines else 0
    write_tsv(config_dir / "eval.baseline_diagnostics.tsv", baseline_rows or [{"config_name": config_name, "available": 0, "reason": "disabled"}])
    null_rows, null_summary = sampled_trans_null_diagnostics(config_name, trans_payload, trans_chrpair_rows_out, sample_seed)
    write_tsv(config_dir / "eval.sampled_trans_null.tsv", null_rows or [{"config_name": config_name, "available": 0, "reason": "insufficient_trans_samples"}])
    wg_summary = whole_genome_summary(best_rows, trans_payload)
    write_tsv(config_dir / "eval.whole_genome_distance_summary.tsv", [{"config_name": config_name, **wg_summary}])
    pro_rows, pro_summary = procrustes_rows(config_name, points)
    write_tsv(config_dir / "eval.procrustes_by_chr.tsv", pro_rows)
    write_tsv(config_dir / "eval.procrustes_summary.tsv", [pro_summary])

    cis_summary = aggregate_cis_summary(config_name, best_rows, rg_rows, scale_rows, quantile_rows)
    write_tsv(config_dir / "eval.cis_summary.tsv", [cis_summary])

    out = dict(summary_row.to_dict())
    out.update(cis_summary)
    out.update(compaction_summary(config_name, rg_rows, scale_rows, quantile_rows))
    out.update({k: v for k, v in centroid_summary.items() if k != "global_mapping_mode"})
    out["global_mapping_mode"] = centroid_summary.get("global_mapping_mode", "")
    out.update(trans_summary)
    out.update(baseline_summary)
    out.update(null_summary)
    out.update(wg_summary)
    out.update(pro_summary)
    out.update(posterior_diag)
    out.update(read_force_class_diag(force_diag_path_from_summary(summary_row)))
    out["global_vs_per_chr_gauge_gain"] = safe_float(out.get("all_chr_cis_min_copy_spearman_mean")) - safe_float(out.get("centroid_distance_spearman"))
    out["eval_dir"] = str(config_dir)
    return out


def aggregate_cis_summary(
    config_name: str,
    best_rows: Sequence[Mapping[str, object]],
    rg_rows: Sequence[Mapping[str, object]],
    scale_rows: Sequence[Mapping[str, object]],
    quantile_rows: Sequence[Mapping[str, object]],
) -> Dict[str, object]:
    valid = [r for r in best_rows if bool(r.get("valid"))]
    chr1 = next((r for r in valid if r["chrom"] == "chr1"), {})
    pm, pr, lm, lr = collect_cis_distance_arrays(best_rows)
    long_slope, _ = distance_regression(lm, lr)
    return {
        "config_name": config_name,
        "n_chr_eval": len(valid),
        "chr1_cis_min_copy_pearson": chr1.get("min_copy_pearson", math.nan),
        "chr1_cis_min_copy_spearman": chr1.get("min_copy_spearman", math.nan),
        "chr1_copy_identity_margin_min_pearson": chr1.get("copy_identity_margin_min_pearson", math.nan),
        "chr1_copy_contrast_score_pearson": chr1.get("copy_contrast_score_pearson", math.nan),
        "all_chr_cis_min_copy_pearson_mean": finite_mean([safe_float(r.get("min_copy_pearson")) for r in valid]),
        "all_chr_cis_min_copy_pearson_median": finite_median([safe_float(r.get("min_copy_pearson")) for r in valid]),
        "all_chr_cis_min_copy_spearman_mean": finite_mean([safe_float(r.get("min_copy_spearman")) for r in valid]),
        "all_chr_cis_pooled_pearson": pearson_corr(pm, pr),
        "all_chr_cis_pooled_spearman": spearman_corr(pm, pr),
        "long_range_cis_pooled_pearson": pearson_corr(lm, lr),
        "long_range_cis_pooled_spearman": spearman_corr(lm, lr),
        "long_range_cis_slope_median": long_slope,
    }


def compaction_summary(
    config_name: str,
    rg_rows: Sequence[Mapping[str, object]],
    scale_rows: Sequence[Mapping[str, object]],
    quantile_rows: Sequence[Mapping[str, object]],
) -> Dict[str, object]:
    rg = [safe_float(r.get("rg_ratio")) for r in rg_rows]
    slopes = [safe_float(r.get("slope")) for r in scale_rows]
    q90 = [safe_float(r.get("ratio_q90")) for r in quantile_rows]
    return {
        "config_name": config_name,
        "median_rg_ratio": finite_median(rg),
        "min_rg_ratio": float(np.nanmin(rg)) if np.isfinite(rg).any() else math.nan,
        "median_q90_distance_ratio": finite_median(q90),
        "median_cis_distance_slope": finite_median(slopes),
        "fraction_chr_rg_ratio_lt_0_7": fraction_lt(rg, 0.7),
        "fraction_chr_slope_lt_0_7": fraction_lt(slopes, 0.7),
    }


def centroid_pair_rows(config_name: str, centroids: pd.DataFrame) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    if centroids.empty:
        return rows
    model = centroids[["model_x", "model_y", "model_z"]].to_numpy(float)
    ref = centroids[["ref_x", "ref_y", "ref_z"]].to_numpy(float)
    for i in range(len(centroids)):
        for j in range(i + 1, len(centroids)):
            md = float(np.linalg.norm(model[i] - model[j]))
            rd = float(np.linalg.norm(ref[i] - ref[j]))
            rows.append({
                "config_name": config_name,
                "left_chr": centroids.iloc[i]["chr"],
                "left_copy": centroids.iloc[i]["copy"],
                "left_hap": centroids.iloc[i]["hap"],
                "right_chr": centroids.iloc[j]["chr"],
                "right_copy": centroids.iloc[j]["copy"],
                "right_hap": centroids.iloc[j]["hap"],
                "model_centroid_distance": md,
                "ref_centroid_distance": rd,
                "distance_ratio": md / rd if rd > 0 else math.nan,
            })
    return rows


def annotate_acceptability(summary: pd.DataFrame) -> pd.DataFrame:
    frame = summary.copy()
    if "trans_relative_metric" not in frame.columns and "trans_distance_spearman" in frame.columns:
        frame["trans_relative_metric"] = frame["trans_distance_spearman"]
    if "centroid_relative_metric" not in frame.columns and "centroid_distance_spearman" in frame.columns:
        frame["centroid_relative_metric"] = frame["centroid_distance_spearman"]
    if "trans_absolute_scale_metric" not in frame.columns and "trans_median_distance_ratio" in frame.columns:
        frame["trans_absolute_scale_metric"] = frame["trans_median_distance_ratio"]
    if "centroid_absolute_scale_metric" not in frame.columns and "centroid_median_distance_ratio" in frame.columns:
        frame["centroid_absolute_scale_metric"] = frame["centroid_median_distance_ratio"]
    if "baseline_diagnostics_available" not in frame.columns:
        frame["baseline_diagnostics_available"] = 0
    for col in [
        "baseline_random_diploid_trans_relative_metric",
        "baseline_shuffled_chromosome_labels_trans_relative_metric",
        "baseline_centroid_only_trans_relative_metric",
        "baseline_compact_ball_per_chromosome_trans_relative_metric",
        "baseline_per_chromosome_translated_scaffold_trans_relative_metric",
        "trans_null_permutation_spearman_p50",
        "trans_null_permutation_spearman_p95",
        "trans_bootstrap_spearman_q025",
        "trans_bootstrap_spearman_q975",
        "trans_chrpair_bootstrap_spearman_mean_q025",
        "trans_chrpair_bootstrap_spearman_mean_q975",
    ]:
        if col not in frame.columns:
            frame[col] = math.nan

    if "trans_relative_metric" in frame.columns:
        frame["trans_metric_rank"] = frame["trans_relative_metric"].rank(method="min", ascending=False, na_option="bottom").astype(int)
        frame["trans_metric_rank_best"] = (frame["trans_metric_rank"] == 1).astype(int)
    else:
        frame["trans_metric_rank"] = math.nan
        frame["trans_metric_rank_best"] = 0
    if "centroid_relative_metric" in frame.columns:
        frame["centroid_metric_rank"] = frame["centroid_relative_metric"].rank(method="min", ascending=False, na_option="bottom").astype(int)
        frame["centroid_metric_rank_best"] = (frame["centroid_metric_rank"] == 1).astype(int)
    else:
        frame["centroid_metric_rank"] = math.nan
        frame["centroid_metric_rank_best"] = 0

    trans_null = frame.get("trans_null_permutation_spearman_p95", pd.Series([math.nan] * len(frame), index=frame.index))
    frame["trans_relative_exceeds_null_margin"] = [
        metric_pass(v, safe_float(n) + TRANS_NULL_MARGIN if math.isfinite(safe_float(n)) else math.nan)
        for v, n in zip(frame.get("trans_relative_metric", pd.Series([math.nan] * len(frame))), trans_null)
    ]
    frame["trans_slope_not_severely_compressed"] = [ratio_not_compressed(v) for v in frame.get("trans_distance_slope", pd.Series([math.nan] * len(frame)))]
    frame["trans_distance_ratio_not_severely_compressed"] = [ratio_not_compressed(v) for v in frame.get("trans_median_distance_ratio", pd.Series([math.nan] * len(frame)))]
    frame["trans_centroid_metric_agreement_pass"] = [
        int(
            math.isfinite(safe_float(t)) and math.isfinite(safe_float(c)) and
            abs(safe_float(t) - safe_float(c)) <= 0.15 and safe_float(c) >= 0.50
        )
        for t, c in zip(frame.get("trans_relative_metric", pd.Series([math.nan] * len(frame))),
                        frame.get("centroid_relative_metric", pd.Series([math.nan] * len(frame))))
    ]
    frame["trans_contact_consistency_available"] = frame.get("posterior_diag_available", pd.Series([0] * len(frame))).fillna(0).astype(int)
    frame["trans_contact_consistency_pass"] = [
        int(avail and math.isfinite(safe_float(v)) and safe_float(v) > 0.10)
        for avail, v in zip(frame["trans_contact_consistency_available"],
                            frame.get("trans_contact_count_vs_inv_distance_spearman", pd.Series([math.nan] * len(frame))))
    ]
    frame["trans_all_chr_cis_nontrivial_pass"] = [
        int(math.isfinite(safe_float(v)) and safe_float(v) >= TRANS_CIS_MIN_PEARSON)
        for v in frame.get("all_chr_cis_min_copy_pearson_mean", pd.Series([math.nan] * len(frame)))
    ]
    frame["trans_rg_scale_no_collapse_pass"] = [
        int(
            math.isfinite(safe_float(rg)) and math.isfinite(safe_float(slope)) and
            safe_float(rg) >= ABS_RATIO_MIN and safe_float(slope) >= ABS_RATIO_MIN
        )
        for rg, slope in zip(frame.get("median_rg_ratio", pd.Series([math.nan] * len(frame))),
                             frame.get("median_cis_distance_slope", pd.Series([math.nan] * len(frame))))
    ]
    trans_components = [
        "trans_relative_exceeds_null_margin",
        "trans_slope_not_severely_compressed",
        "trans_distance_ratio_not_severely_compressed",
        "trans_centroid_metric_agreement_pass",
        "trans_contact_consistency_pass",
        "trans_all_chr_cis_nontrivial_pass",
        "trans_rg_scale_no_collapse_pass",
    ]
    frame["trans_metric_absolute_pass"] = frame[trans_components].all(axis=1).astype(int)

    frame["centroid_slope_not_severely_compressed"] = [ratio_not_compressed(v) for v in frame.get("centroid_distance_slope", pd.Series([math.nan] * len(frame)))]
    frame["centroid_distance_ratio_not_severely_compressed"] = [ratio_not_compressed(v) for v in frame.get("centroid_median_distance_ratio", pd.Series([math.nan] * len(frame)))]
    frame["centroid_relative_nontrivial_pass"] = [
        int(math.isfinite(safe_float(v)) and safe_float(v) >= 0.50)
        for v in frame.get("centroid_relative_metric", pd.Series([math.nan] * len(frame)))
    ]
    centroid_components = [
        "centroid_relative_nontrivial_pass",
        "centroid_slope_not_severely_compressed",
        "centroid_distance_ratio_not_severely_compressed",
        "trans_all_chr_cis_nontrivial_pass",
        "trans_rg_scale_no_collapse_pass",
    ]
    frame["centroid_metric_absolute_pass"] = frame[centroid_components].all(axis=1).astype(int)

    def failed_components(row: pd.Series, components: Sequence[str]) -> str:
        return ",".join([name for name in components if int(row.get(name, 0)) == 0])

    frame["trans_metric_absolute_fail_reasons"] = [failed_components(row, trans_components) for _, row in frame.iterrows()]
    frame["centroid_metric_absolute_fail_reasons"] = [failed_components(row, centroid_components) for _, row in frame.iterrows()]
    frame["trans_metric_rank_note"] = "relative_rank_only_not_absolute_success"
    frame["centroid_metric_rank_note"] = "relative_rank_only_not_absolute_success"
    return frame


def warning_rows(summary: pd.DataFrame) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = [{
        "warning_type": "global_interpretation",
        "config_name": "",
        "message": (
            "The top trans/centroid-ranked configs are not absolute-pass trans candidates. "
            "They mostly indicate that current chr1/cis optimization and current trans/centroid "
            "ranking metrics pull in different directions. Absolute trans reconstruction remains poor."
        ),
    }]
    for _, row in summary.iterrows():
        name = str(row.get("config_name", ""))
        if int(row.get("trans_metric_rank_best", 0)) == 1 and int(row.get("trans_metric_absolute_pass", 0)) == 0:
            rows.append({
                "warning_type": "relative_top_trans_not_absolute_pass",
                "config_name": name,
                "message": "Relative top trans metric only; not biologically acceptable trans reconstruction.",
                "fail_reasons": row.get("trans_metric_absolute_fail_reasons", ""),
            })
        if int(row.get("centroid_metric_rank_best", 0)) == 1 and int(row.get("centroid_metric_absolute_pass", 0)) == 0:
            rows.append({
                "warning_type": "relative_top_centroid_not_absolute_pass",
                "config_name": name,
                "message": "Relative top centroid metric only; not an absolute-pass centroid candidate.",
                "fail_reasons": row.get("centroid_metric_absolute_fail_reasons", ""),
            })
        if safe_float(row.get("chr1_cis_min_copy_pearson")) >= 0.55 and int(row.get("trans_metric_absolute_pass", 0)) == 0:
            rows.append({
                "warning_type": "cis_high_trans_absolute_fail",
                "config_name": name,
                "message": "High chr1/cis metric does not imply acceptable global trans layout.",
                "fail_reasons": row.get("trans_metric_absolute_fail_reasons", ""),
            })
    return rows


def write_rankings(root: Path, summary: pd.DataFrame) -> None:
    for stale in ["rank_by_trans.tsv", "rank_by_centroid.tsv"]:
        stale_path = root / stale
        if stale_path.exists():
            stale_path.unlink()
    ranking_specs = [
        ("rank_by_chr1_cis.tsv", "chr1_cis_min_copy_pearson", False),
        ("rank_by_all_chr_cis.tsv", "all_chr_cis_min_copy_pearson_mean", False),
        ("rank_by_trans_relative_metric.tsv", "trans_relative_metric", False),
        ("rank_by_centroid_relative_metric.tsv", "centroid_relative_metric", False),
        ("rank_by_compaction.tsv", "median_rg_ratio", False),
        ("pareto_summary.tsv", "pareto_score", False),
    ]
    frame = summary.copy()
    frame["pareto_score"] = (
        frame["all_chr_cis_min_copy_spearman_mean"].fillna(-1)
        + frame["trans_relative_metric"].fillna(-1)
        + frame["centroid_relative_metric"].fillna(-1)
        - (frame["median_rg_ratio"].fillna(0) - 1.0).abs()
    )
    keep = [
        "config_name", "mode_group", "n_iter", "relax_steps", "multiplier",
        "chr1_cis_min_copy_pearson", "chr1_copy_identity_margin_min_pearson",
        "all_chr_cis_min_copy_pearson_mean", "all_chr_cis_pooled_spearman",
        "trans_relative_metric", "trans_metric_rank", "trans_metric_rank_best",
        "trans_metric_absolute_pass", "trans_metric_absolute_fail_reasons",
        "trans_distance_spearman", "trans_distance_slope", "trans_median_distance_ratio",
        "centroid_relative_metric", "centroid_metric_rank", "centroid_metric_rank_best",
        "centroid_metric_absolute_pass", "centroid_metric_absolute_fail_reasons",
        "centroid_distance_spearman", "centroid_distance_slope", "centroid_median_distance_ratio",
        "median_rg_ratio",
        "median_cis_distance_slope", "procrustes_global_over_per_chr_ratio",
        "pareto_score", "eval_dir",
    ]
    keep = [col for col in keep if col in frame.columns]
    for name, metric, ascending in ranking_specs:
        ranked = frame.sort_values(metric, ascending=ascending, na_position="last").copy()
        ranked["rank_metric_name"] = metric
        ranked["rank_metric_kind"] = "relative_spearman" if metric in {"trans_relative_metric", "centroid_relative_metric"} else "diagnostic_metric"
        ranked[["rank_metric_name", "rank_metric_kind", *keep]].to_csv(root / name, sep="\t", index=False)


def write_mechanism_tables(root: Path, summary: pd.DataFrame) -> None:
    scatter_cols = [
        "config_name", "mode_group", "n_iter", "relax_steps", "multiplier",
        "chr1_cis_min_copy_pearson", "chr1_copy_identity_margin_min_pearson",
        "all_chr_cis_min_copy_pearson_mean",
        "trans_relative_metric", "trans_metric_rank", "trans_metric_rank_best",
        "trans_metric_absolute_pass", "trans_distance_slope", "trans_median_distance_ratio",
        "trans_abs_log_ratio_median", "trans_centroid_residual_spearman",
        "centroid_relative_metric", "centroid_metric_rank", "centroid_metric_rank_best",
        "centroid_metric_absolute_pass", "centroid_distance_slope", "centroid_median_distance_ratio",
        "median_rg_ratio", "median_cis_distance_slope",
        "mean_pU", "final_mean_rho_train_bpair", "mean_pU_trans",
        "effective_k_cis_trans_ratio", "force_l1_cis_trans_ratio",
    ]
    scatter_cols = [col for col in scatter_cols if col in summary.columns]
    summary[scatter_cols].to_csv(root / "eval_mechanism_scatter.tsv", sep="\t", index=False)

    metrics = [
        "chr1_cis_min_copy_pearson",
        "all_chr_cis_min_copy_pearson_mean",
        "trans_relative_metric",
        "trans_distance_slope",
        "trans_median_distance_ratio",
        "trans_abs_log_ratio_median",
        "trans_centroid_residual_spearman",
        "centroid_relative_metric",
        "centroid_distance_slope",
        "centroid_median_distance_ratio",
        "median_rg_ratio",
        "median_cis_distance_slope",
        "mean_pU_trans",
        "effective_k_cis_trans_ratio",
        "force_l1_cis_trans_ratio",
    ]
    rows = []
    for metric in metrics:
        if metric not in summary.columns:
            continue
        for row in summary.itertuples(index=False):
            rows.append({
                "metric": metric,
                "mode_group": getattr(row, "mode_group", ""),
                "n_iter": getattr(row, "n_iter", ""),
                "relax_steps": getattr(row, "relax_steps", ""),
                "multiplier": getattr(row, "multiplier", ""),
                "value": getattr(row, metric),
                "config_name": getattr(row, "config_name"),
            })
    pd.DataFrame(rows).to_csv(root / "eval_heatmap_matrix.tsv", sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--tdg", type=Path, default=DEFAULT_TDG)
    parser.add_argument("--out-root", type=Path, default=None)
    parser.add_argument("--annotate-existing-summary", type=Path, default=None, help="Only add corrected flags/rank tables to an existing eval_summary.tsv")
    parser.add_argument("--trans-sample-per-chrpair", type=int, default=2048)
    parser.add_argument("--baseline-trans-sample-per-chrpair", type=int, default=256)
    parser.add_argument("--sample-seed", type=int, default=17)
    parser.add_argument("--max-configs", type=int, default=0)
    parser.add_argument("--config", action="append", default=None, help="Evaluate only named config(s)")
    parser.add_argument("--no-baselines", action="store_true", help="Skip eval-only baseline coordinate transforms")
    args = parser.parse_args()

    if args.out_root is None:
        args.out_root = args.annotate_existing_summary.parent if args.annotate_existing_summary is not None else args.summary.parent / "eval"
    args.out_root.mkdir(parents=True, exist_ok=True)

    if args.annotate_existing_summary is not None:
        eval_summary = pd.read_csv(args.annotate_existing_summary, sep="\t")
        eval_summary = annotate_acceptability(eval_summary)
        eval_summary.to_csv(args.annotate_existing_summary, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
        write_rankings(args.out_root, eval_summary)
        write_mechanism_tables(args.out_root, eval_summary)
        write_tsv(args.out_root / "eval_acceptability_warnings.tsv", warning_rows(eval_summary))
        print(f"annotated {len(eval_summary)} rows in {args.annotate_existing_summary}")
        return 0

    summary_path = args.summary
    if summary_path == DEFAULT_SUMMARY and not summary_path.exists() and DEFAULT_MATRIX_SUMMARY.exists():
        summary_path = DEFAULT_MATRIX_SUMMARY
    summary = pd.read_csv(summary_path, sep="\t")
    if "config_name" not in summary.columns and "config_id" in summary.columns:
        summary = summary.rename(columns={"config_id": "config_name"})
    if "status" in summary.columns:
        summary = summary.loc[summary["status"].astype(str) == "OK"].copy()
    if args.config:
        summary = summary.loc[summary["config_name"].isin(args.config)].copy()
    if args.max_configs and args.max_configs > 0:
        summary = summary.head(args.max_configs).copy()
    ref_df = load_charm_3dg(args.tdg)

    rows: List[Dict[str, object]] = []
    for idx, row in enumerate(summary.itertuples(index=False), start=1):
        series = pd.Series(row._asdict())
        config_name = str(series["config_name"])
        print(f"[{idx}/{len(summary)}] evaluating {config_name}", flush=True)
        rows.append(evaluate_config(
            series,
            ref_df,
            args.out_root,
            args.trans_sample_per_chrpair,
            args.sample_seed,
            args.baseline_trans_sample_per_chrpair,
            run_baselines=not args.no_baselines,
        ))

    eval_summary = annotate_acceptability(pd.DataFrame(rows))
    eval_summary_path = args.out_root / "eval_summary.tsv"
    eval_summary.to_csv(eval_summary_path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    write_rankings(args.out_root, eval_summary)
    write_mechanism_tables(args.out_root, eval_summary)
    write_tsv(args.out_root / "eval_acceptability_warnings.tsv", warning_rows(eval_summary))
    print(f"wrote {len(eval_summary)} rows to {eval_summary_path}")
    for metric in ["chr1_cis_min_copy_pearson", "all_chr_cis_min_copy_pearson_mean", "trans_relative_metric", "centroid_relative_metric", "median_rg_ratio"]:
        if metric in eval_summary.columns:
            top = eval_summary.sort_values(metric, ascending=False).head(5)
            label = f"relative top by {metric}" if metric in {"trans_relative_metric", "centroid_relative_metric"} else f"top by {metric}"
            print(f"\n{label}")
            cols = [c for c in ["config_name", "mode_group", "state_weight_mode", "n_iter", "relax_steps", metric] if c in top.columns]
            if metric in {"trans_relative_metric", "centroid_relative_metric"}:
                cols += [c for c in ["trans_metric_absolute_pass", "centroid_metric_absolute_pass"] if c in top.columns]
            print(top[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
