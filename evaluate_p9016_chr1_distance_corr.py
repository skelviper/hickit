#!/usr/bin/env python3
"""Evaluate chr1 distance-map correlation against CHARM 3DG.

The metric compares within-copy chr1 pairwise Euclidean distance vectors from a
blind reconstruction to the corresponding CHARM 3DG haplotype distance vectors.
Because model copy labels have gauge symmetry, both copy-to-haplotype mappings
are evaluated and reported.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd


DEFAULT_SUMMARY = Path(
    "/tmp/hk_blind_p9016_main_random_repel_grid_4135715_1777429703_0/scan_summary.tsv"
)
DEFAULT_TDG = Path("/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz")
DEFAULT_OUT_NAME = "chr1_distance_corr_summary.tsv"

CopyHap = Tuple[int, str]
Permutation = Tuple[CopyHap, CopyHap]

COPY_HAP_PERMUTATIONS: Tuple[Permutation, Permutation] = (
    ((0, "mat"), (1, "pat")),
    ((0, "pat"), (1, "mat")),
)


def read_tdg(path: Path) -> pd.DataFrame:
    records: List[Tuple[str, str, int, float, float, float]] = []
    with gzip.open(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                raise ValueError(f"Expected at least 5 columns in {path} line {line_number}")
            chrom_hap = fields[0]
            if not chrom_hap.endswith(")") or "(" not in chrom_hap:
                raise ValueError(f"Expected chr(hap) label in {path} line {line_number}: {chrom_hap!r}")
            chrom, hap = chrom_hap[:-1].split("(", 1)
            records.append((chrom, hap, int(fields[1]), float(fields[2]), float(fields[3]), float(fields[4])))
    frame = pd.DataFrame(records, columns=["chr", "hap", "start", "x", "y", "z"])
    if frame.empty:
        raise ValueError(f"No 3DG records read from {path}")
    duplicate_mask = frame.duplicated(["chr", "hap", "start"], keep=False)
    if duplicate_mask.any():
        examples = frame.loc[duplicate_mask, ["chr", "hap", "start"]].drop_duplicates().head().to_dict("records")
        raise ValueError(f"Duplicate 3DG bins found: {examples}")
    return frame


def load_model_chr(path: Path, target_chr: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", usecols=["chr", "start", "copy", "x", "y", "z"])
    frame = frame.loc[frame["chr"] == target_chr].copy()
    frame["copy"] = frame["copy"].astype(int)
    duplicate_mask = frame.duplicated(["chr", "copy", "start"], keep=False)
    if duplicate_mask.any():
        examples = frame.loc[duplicate_mask, ["chr", "copy", "start"]].drop_duplicates().head().to_dict("records")
        raise ValueError(f"Duplicate model bins found in {path}: {examples}")
    return frame


def distance_vector(coords: np.ndarray) -> np.ndarray:
    coords = np.asarray(coords, dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"Expected Nx3 coordinate array, got {coords.shape}")
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2, dtype=np.float64))
    upper = np.triu_indices(coords.shape[0], k=1)
    return dist[upper]


def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    if int(finite.sum()) < 3:
        return math.nan
    x = x[finite]
    y = y[finite]
    x = x - x.mean()
    y = y - y.mean()
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
        rank = 0.5 * (i + j - 1) + 1.0
        ranks[order[i:j]] = rank
        i = j
    return ranks


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    if int(finite.sum()) < 3:
        return math.nan
    return pearson_corr(rankdata_average(x[finite]), rankdata_average(y[finite]))


def permutation_name(permutation: Permutation) -> str:
    return "_".join(f"copy{copy}_{hap}" for copy, hap in permutation)


def build_reference_lookup(tdg_df: pd.DataFrame, target_chr: str) -> Dict[str, pd.DataFrame]:
    lookup: Dict[str, pd.DataFrame] = {}
    for hap in sorted(tdg_df["hap"].unique()):
        part = (
            tdg_df.loc[(tdg_df["chr"] == target_chr) & (tdg_df["hap"] == hap), ["start", "x", "y", "z"]]
            .sort_values("start", kind="mergesort")
            .reset_index(drop=True)
        )
        lookup[hap] = part
    return lookup


def coords_by_start(frame: pd.DataFrame) -> Mapping[int, np.ndarray]:
    return {
        int(row.start): np.array([row.x, row.y, row.z], dtype=np.float64)
        for row in frame.itertuples(index=False)
    }


def evaluate_copy_pair(
    model_df: pd.DataFrame,
    ref_lookup: Mapping[str, pd.DataFrame],
    copy: int,
    hap: str,
) -> Dict[str, object]:
    model_part = model_df.loc[model_df["copy"] == copy, ["start", "x", "y", "z"]]
    ref_part = ref_lookup[hap]
    model_by_start = coords_by_start(model_part)
    ref_by_start = coords_by_start(ref_part)
    common_starts = sorted(set(model_by_start).intersection(ref_by_start))
    if len(common_starts) < 3:
        return {
            "copy": copy,
            "hap": hap,
            "n_common_bins": len(common_starts),
            "n_distance_pairs": 0,
            "pearson": math.nan,
            "spearman": math.nan,
        }
    model_coords = np.vstack([model_by_start[start] for start in common_starts])
    ref_coords = np.vstack([ref_by_start[start] for start in common_starts])
    model_dist = distance_vector(model_coords)
    ref_dist = distance_vector(ref_coords)
    return {
        "copy": copy,
        "hap": hap,
        "n_common_bins": len(common_starts),
        "n_distance_pairs": int(model_dist.shape[0]),
        "pearson": pearson_corr(model_dist, ref_dist),
        "spearman": spearman_corr(model_dist, ref_dist),
        "model_dist": model_dist,
        "ref_dist": ref_dist,
    }


def evaluate_permutation(
    model_df: pd.DataFrame,
    ref_lookup: Mapping[str, pd.DataFrame],
    permutation: Permutation,
) -> Dict[str, object]:
    copy_results = [evaluate_copy_pair(model_df, ref_lookup, copy, hap) for copy, hap in permutation]
    pearsons = [float(result["pearson"]) for result in copy_results]
    spearmans = [float(result["spearman"]) for result in copy_results]
    valid = all("model_dist" in result for result in copy_results) and all(math.isfinite(x) for x in pearsons + spearmans)
    if valid:
        pooled_model = np.concatenate([result["model_dist"] for result in copy_results])
        pooled_ref = np.concatenate([result["ref_dist"] for result in copy_results])
        min_copy_pearson = float(np.min(pearsons))
        mean_copy_pearson = float(np.mean(pearsons))
        pooled_pearson = pearson_corr(pooled_model, pooled_ref)
        min_copy_spearman = float(np.min(spearmans))
        mean_copy_spearman = float(np.mean(spearmans))
        pooled_spearman = spearman_corr(pooled_model, pooled_ref)
    else:
        min_copy_pearson = math.nan
        mean_copy_pearson = math.nan
        pooled_pearson = math.nan
        min_copy_spearman = math.nan
        mean_copy_spearman = math.nan
        pooled_spearman = math.nan
    return {
        "name": permutation_name(permutation),
        "valid": valid,
        "copy_results": copy_results,
        "min_copy_pearson": min_copy_pearson,
        "mean_copy_pearson": mean_copy_pearson,
        "pooled_pearson": pooled_pearson,
        "min_copy_spearman": min_copy_spearman,
        "mean_copy_spearman": mean_copy_spearman,
        "pooled_spearman": pooled_spearman,
        "n_common_bins_min": min(int(result["n_common_bins"]) for result in copy_results),
        "n_distance_pairs_total": sum(int(result["n_distance_pairs"]) for result in copy_results),
    }


def fmt(value: object) -> object:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.9g}"
    return value


def evaluate_summary(summary_path: Path, tdg_path: Path, target_chr: str) -> pd.DataFrame:
    summary_df = pd.read_csv(summary_path, sep="\t")
    tdg_df = read_tdg(tdg_path)
    ref_lookup = build_reference_lookup(tdg_df, target_chr)
    rows: List[Dict[str, object]] = []
    for summary_row in summary_df.itertuples(index=False):
        coords_path = Path(summary_row.coords_gz)
        model_df = load_model_chr(coords_path, target_chr)
        permutation_results = [
            evaluate_permutation(model_df, ref_lookup, permutation)
            for permutation in COPY_HAP_PERMUTATIONS
        ]
        valid_results = [result for result in permutation_results if result["valid"]]
        if valid_results:
            best_by_min = max(valid_results, key=lambda result: float(result["min_copy_pearson"]))
            best_by_pooled = max(valid_results, key=lambda result: float(result["pooled_pearson"]))
            literal_min = min(float(result["min_copy_pearson"]) for result in valid_results)
        else:
            best_by_min = permutation_results[0]
            best_by_pooled = permutation_results[0]
            literal_min = math.nan

        out: Dict[str, object] = {
            "config_name": summary_row.config_name,
            "mode_group": getattr(summary_row, "mode_group", ""),
            "multiplier": summary_row.multiplier,
            "k_rel_rep": summary_row.k_rel_rep,
            "n_iter": summary_row.n_iter,
            "relax_steps": summary_row.relax_steps,
            "mean_pU": summary_row.mean_pU,
            "mean_pmax": summary_row.mean_pmax,
            "mean_margin": summary_row.mean_margin,
            "sep_mean": summary_row.sep_mean,
            "contact_energy_per_wedge_k": summary_row.contact_energy_per_wedge_k,
            "coords_gz": str(coords_path),
            "target_chr": target_chr,
            "label_invariant_min_copy_pearson": best_by_min["min_copy_pearson"],
            "label_invariant_best_mapping": best_by_min["name"],
            "label_invariant_mean_copy_pearson": best_by_min["mean_copy_pearson"],
            "label_invariant_pooled_pearson_for_best_min": best_by_min["pooled_pearson"],
            "label_invariant_min_copy_spearman": best_by_min["min_copy_spearman"],
            "label_invariant_pooled_spearman_for_best_min": best_by_min["pooled_spearman"],
            "best_pooled_pearson": best_by_pooled["pooled_pearson"],
            "best_pooled_mapping": best_by_pooled["name"],
            "literal_min_over_mappings_min_copy_pearson": literal_min,
            "n_common_bins_min": best_by_min["n_common_bins_min"],
            "n_distance_pairs_total": best_by_min["n_distance_pairs_total"],
        }
        for result in permutation_results:
            prefix = result["name"]
            out[f"{prefix}_min_copy_pearson"] = result["min_copy_pearson"]
            out[f"{prefix}_mean_copy_pearson"] = result["mean_copy_pearson"]
            out[f"{prefix}_pooled_pearson"] = result["pooled_pearson"]
            out[f"{prefix}_min_copy_spearman"] = result["min_copy_spearman"]
            for copy_result in result["copy_results"]:
                copy = copy_result["copy"]
                hap = copy_result["hap"]
                out[f"{prefix}_copy{copy}_{hap}_pearson"] = copy_result["pearson"]
                out[f"{prefix}_copy{copy}_{hap}_spearman"] = copy_result["spearman"]
        rows.append(out)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--tdg", type=Path, default=DEFAULT_TDG)
    parser.add_argument("--chr", default="chr1", dest="target_chr")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    if args.out is None:
        args.out = args.summary.parent / DEFAULT_OUT_NAME
    result = evaluate_summary(args.summary, args.tdg, args.target_chr)
    result.to_csv(args.out, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)

    print(f"wrote {len(result)} rows to {args.out}")
    for mode, group in result.groupby("mode_group", sort=True):
        top = group.sort_values("label_invariant_min_copy_pearson", ascending=False).head(5)
        print(f"\nmode={mode} top by label_invariant_min_copy_pearson")
        print(
            top[
                [
                    "config_name",
                    "label_invariant_min_copy_pearson",
                    "label_invariant_best_mapping",
                    "label_invariant_pooled_pearson_for_best_min",
                    "mean_pU",
                    "n_common_bins_min",
                ]
            ].to_string(index=False)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
