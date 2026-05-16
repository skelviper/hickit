#!/usr/bin/env python3
"""Post-process full-output P9016 candidate reruns.

This module is eval-only.  It normalizes the C rerun manifest/summary output,
invokes CHARM-based evaluation only after training has finished, and writes the
candidate-level diagnostic tables requested for the rerun audit.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence

import numpy as np
import pandas as pd

from blind_eval import load_model_coords, spearman_corr
from evaluate_blind_p9016_grid import (
    DEFAULT_TDG,
    LOG4,
    annotate_acceptability,
    evaluate_config,
    load_charm_3dg,
    rho_train_for_bpair,
    safe_float,
    warning_rows,
    write_mechanism_tables,
    write_rankings,
    write_tsv,
)


CANDIDATE_ORDER = [
    "main_rep1_n100_rs100",
    "main_rep1_n120_rs80",
    "main_rep1_n80_rs100",
    "main_rep1_n5_rs80",
    "main_rep1_n5_rs5",
    "random_rep1_n120_rs50",
]

PRIOR_RHO_2X2_ORDER = [
    "prior_rho2x2_n100_rs100_cis_entropy",
    "prior_rho2x2_n100_rs100_uniform_entropy",
    "prior_rho2x2_n100_rs100_cis_constant",
    "prior_rho2x2_n100_rs100_uniform_constant",
    "prior_rho2x2_n120_rs80_cis_entropy",
    "prior_rho2x2_n120_rs80_uniform_entropy",
    "prior_rho2x2_n120_rs80_cis_constant",
    "prior_rho2x2_n120_rs80_uniform_constant",
]

SCAFFOLD_STAGE_ORDER = [
    "main_rep1_n100_rs100",
    "main_rep1_n120_rs80",
    "prior_rho2x2_n100_rs100_uniform_entropy",
    "prior_rho2x2_n100_rs100_cis_constant",
]

DSCALE_ABLATION_ORDER = [
    "dscale_n100_rs100_raw",
    "dscale_n100_rs100_expected",
    "dscale_n120_rs80_raw",
    "dscale_n120_rs80_expected",
    "dscale_n80_rs100_raw",
    "dscale_n80_rs100_expected",
]

RUN_SET_CONFIG_ORDERS = {
    "candidate_rerun": CANDIDATE_ORDER,
    "prior_rho_2x2": PRIOR_RHO_2X2_ORDER,
    "scaffold_stage": SCAFFOLD_STAGE_ORDER,
    "dscale_ablation": DSCALE_ABLATION_ORDER,
}

POSTERIOR_REQUIRED_COLUMNS = {
    "chr1",
    "start1",
    "end1",
    "chr2",
    "start2",
    "end2",
    "bid1",
    "bid2",
    "n_raw",
    "base_d_scale",
    "base_k",
    "p00",
    "p01",
    "p10",
    "p11",
    "pU",
    "entropy",
    "margin",
    "pmax",
    "rho_output",
    "contact_class",
}

FORCE_REQUIRED_COLUMNS = {
    "class",
    "n_wedges",
    "sum_wedge_k",
    "contact_energy",
    "contact_force_l1",
    "backbone_force_l1",
    "repulsion_force_l1",
    "homolog_sep_force_l1",
    "n_nonfinite",
}

CONTACT_STATES = [
    ("00", 0, 0, "p00"),
    ("01", 0, 1, "p01"),
    ("10", 1, 0, "p10"),
    ("11", 1, 1, "p11"),
]

PU_BINS = [(0.0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 0.9), (0.9, 0.95), (0.95, 0.99), (0.99, 1.000001)]
PMAX_BINS = [(0.0, 0.5), (0.5, 0.75), (0.75, 0.9), (0.9, 0.99), (0.99, 1.000001)]
D_SCALE_BINS = [(0.0, 0.5), (0.5, 0.75), (0.75, 1.0), (1.0, 1.5), (1.5, math.inf)]
N_RAW_BINS = [(1.0, 2.0), (2.0, 5.0), (5.0, 10.0), (10.0, 25.0), (25.0, math.inf)]
STATE_P_BINS = [
    (0.0, 1e-6),
    (1e-6, 1e-4),
    (1e-4, 1e-3),
    (1e-3, 1e-2),
    (1e-2, 5e-2),
    (5e-2, 1e-1),
    (1e-1, 2.5e-1),
    (2.5e-1, 1.000001),
]
SCAFFOLD_STAGE_FILES = [
    ("init_haploid_scaffold", "init_haploid_scaffold.coords.tsv.gz", False),
    ("init_diploid_split", "init_diploid_split.coords.tsv.gz", False),
    ("after_iter1", "after_iter1.coords.tsv.gz", False),
    ("final", "final.coords.tsv.gz", True),
]

FDG_D_C1 = 0.5
FDG_D_C2 = 1.5
FDG_D_C3 = 2.0
FDG_C_C1 = 3.0 * (FDG_D_C3 - FDG_D_C2)
FDG_C_C2 = (FDG_D_C3 - FDG_D_C2) ** 3


def read_tsv_trimmed(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=object)
    frame.columns = [str(col).strip() for col in frame.columns]
    for col in frame.columns:
        if frame[col].dtype == object:
            frame[col] = frame[col].map(lambda x: x.strip() if isinstance(x, str) else x)
    return frame


def read_manifest(path: Path) -> Dict[str, str]:
    manifest: Dict[str, str] = {}
    with open(path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header != ["key", "value"]:
            raise ValueError(f"{path} is not a key/value manifest")
        for row in reader:
            if len(row) >= 2:
                manifest[row[0]] = row[1]
    return manifest


def first_existing_path(paths: Iterable[Path]) -> str:
    for path in paths:
        if path.exists():
            return str(path)
    return ""


def mode_group_for_config(config_name: str) -> str:
    return "random" if config_name.startswith("random_") else "main"


def normalize_candidate_summary(
    root: Path,
    summary_name: str = "matrix_summary.tsv",
    require_exact_candidates: bool = False,
    expected_order: Sequence[str] | None = None,
) -> pd.DataFrame:
    summary_path = root / summary_name
    if not summary_path.exists():
        raise FileNotFoundError(summary_path)
    frame = read_tsv_trimmed(summary_path)
    if "config_name" not in frame.columns:
        if "config_id" not in frame.columns:
            raise ValueError(f"{summary_path} has neither config_name nor config_id")
        frame = frame.rename(columns={"config_id": "config_name"})
    if require_exact_candidates:
        names = [str(name) for name in frame["config_name"]]
        if expected_order is None:
            expected_order = CANDIDATE_ORDER
        expected = set(expected_order)
        seen = set(names)
        missing = sorted(expected - seen)
        extra = sorted(seen - expected)
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if missing or extra or duplicates:
            raise ValueError(
                "candidate set mismatch: "
                f"missing={','.join(missing) or 'none'} "
                f"extra={','.join(extra) or 'none'} "
                f"duplicates={','.join(duplicates) or 'none'}"
            )

    rows: List[Dict[str, object]] = []
    for raw in frame.to_dict("records"):
        config_name = str(raw["config_name"])
        output_dir = Path(str(raw["output_dir"]))
        manifest_path = output_dir / "p9016_full.manifest.tsv"
        manifest = read_manifest(manifest_path)
        output_coords = manifest.get("output_coords", str(output_dir / "p9016_full.coords.tsv"))
        output_coords_gz = manifest.get("output_coords_gz", str(output_dir / "p9016_full.coords.tsv.gz"))
        coords_gz = output_coords_gz if Path(output_coords_gz).exists() else output_coords
        row: Dict[str, object] = dict(raw)
        row.update({
            "config_name": config_name,
            "mode_group": raw.get("mode_group", mode_group_for_config(config_name)),
            "output_dir": str(output_dir),
            "coords_gz": coords_gz,
            "output_coords": output_coords,
            "output_coords_gz": output_coords_gz if Path(output_coords_gz).exists() else "",
            "output_bpair_posterior": manifest.get("output_bpair_posterior", str(output_dir / "p9016_full.bpair_posterior.tsv")),
            "output_loop_diag": manifest.get("output_loop_diag", str(output_dir / "p9016_full.loop_diag.tsv")),
            "output_force_class_diag": manifest.get("output_force_class_diag", str(output_dir / "p9016_full.force_class_diag.tsv")),
            "output_manifest": str(manifest_path),
            "status": raw.get("status", raw.get("audit_status", manifest.get("status", ""))),
            "multiplier": manifest.get("repulsion_multiplier", raw.get("repulsion_multiplier", raw.get("multiplier", "1"))),
            "k_rel_rep": manifest.get("k_rel_rep_effective", raw.get("k_rel_rep_effective", raw.get("k_rel_rep", ""))),
            "relax_step": manifest.get("relax_step", raw.get("relax_step", "")),
            "temperature_start": manifest.get("temperature_start", raw.get("temperature_start", "")),
            "temperature_end": manifest.get("temperature_end", raw.get("temperature_end", "")),
            "rho_train_start": manifest.get("rho_train_start", raw.get("rho_train_start", "")),
            "rho_train_end": manifest.get("rho_train_end", raw.get("rho_train_end", "")),
            "unit": manifest.get("unit", raw.get("unit", "1")),
            "d_scale": manifest.get("d_scale", raw.get("d_scale", "1")),
            "legacy_base_k_unused": manifest.get("legacy_base_k_unused", raw.get("legacy_base_k_unused", "2")),
            "min_sep_unit": manifest.get("min_sep_unit", raw.get("min_sep_unit", "")),
            "lambda_sep": manifest.get("lambda_sep", raw.get("lambda_sep", "")),
            "init_eps": manifest.get("init_eps_effective", raw.get("init_eps", "")),
            "init_noise_scale": manifest.get("init_noise_scale_effective", raw.get("init_noise_scale", "")),
            "init_seed": manifest.get("init_seed", raw.get("seed", raw.get("init_seed", ""))),
            "init_seed_effective": manifest.get("init_seed", raw.get("seed", raw.get("init_seed_effective", ""))),
            "init_scale": manifest.get("init_scale", raw.get("init_scale", "")),
            "init_mode": manifest.get("init_mode", raw.get("init_mode", "")),
            "prior_mode": manifest.get("prior_mode", raw.get("prior_mode", "")),
            "rho_train_mode": manifest.get("rho_train_mode", raw.get("rho_train_mode", "")),
            "rho_train_floor": manifest.get("rho_train_floor", raw.get("rho_train_floor", "")),
            "d_scale_mode": manifest.get("d_scale_mode", raw.get("d_scale_mode", "")),
            "scaffold_source": manifest.get("scaffold_source", raw.get("scaffold_source", "")),
            "same_bin_filter_enabled": manifest.get("same_bin_filter_enabled", raw.get("same_bin_filter_enabled", "")),
            "n_raw_same_bin_excluded": manifest.get("n_raw_same_bin_excluded", raw.get("n_same_bin_excluded", "")),
            "n_bpair_same_bin_excluded": manifest.get("n_bpair_same_bin_excluded", raw.get("n_bpair_same_bin_excluded", "")),
            "posterior_refreshed_after_final_relax": manifest.get(
                "posterior_refreshed_after_final_relax",
                raw.get("posterior_refreshed_after_final_relax", ""),
            ),
        })
        rows.append(row)

    out = pd.DataFrame(rows)
    if expected_order is None:
        expected_order = CANDIDATE_ORDER
    order = {name: i for i, name in enumerate(expected_order)}
    out["_candidate_order"] = out["config_name"].map(lambda x: order.get(str(x), len(order)))
    out = out.sort_values(["_candidate_order", "config_name"], kind="mergesort").drop(columns=["_candidate_order"])
    return out


def expected_config_order(run_set: str) -> Sequence[str]:
    if run_set not in RUN_SET_CONFIG_ORDERS:
        raise ValueError(f"unknown run set: {run_set}")
    return RUN_SET_CONFIG_ORDERS[run_set]


def validate_candidate_file_schemas(summary: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for item in summary.itertuples(index=False):
        config_name = str(getattr(item, "config_name"))
        posterior_path = Path(str(getattr(item, "output_bpair_posterior")))
        force_path = Path(str(getattr(item, "output_force_class_diag")))
        posterior_cols = set(read_tsv_trimmed(posterior_path).columns) if posterior_path.exists() else set()
        force_cols = set(read_tsv_trimmed(force_path).columns) if force_path.exists() else set()
        posterior_missing = sorted(POSTERIOR_REQUIRED_COLUMNS - posterior_cols)
        force_missing = sorted(FORCE_REQUIRED_COLUMNS - force_cols)
        rows.append({
            "config_name": config_name,
            "posterior_path": str(posterior_path),
            "force_path": str(force_path),
            "posterior_schema_ok": int(not posterior_missing),
            "force_schema_ok": int(not force_missing),
            "posterior_missing_columns": ",".join(posterior_missing),
            "force_missing_columns": ",".join(force_missing),
        })
    schema = pd.DataFrame(rows)
    if not schema.empty and int(schema["posterior_schema_ok"].min()) == 0:
        bad = schema.loc[schema["posterior_schema_ok"] == 0, ["config_name", "posterior_missing_columns"]]
        raise ValueError(f"bpair posterior schema failure: {bad.to_dict('records')}")
    if not schema.empty and int(schema["force_schema_ok"].min()) == 0:
        bad = schema.loc[schema["force_schema_ok"] == 0, ["config_name", "force_missing_columns"]]
        raise ValueError(f"force diag schema failure: {bad.to_dict('records')}")
    return schema


def evaluate_candidates(
    summary: pd.DataFrame,
    tdg_path: Path,
    out_root: Path,
    trans_sample_per_chrpair: int,
    baseline_trans_sample_per_chrpair: int,
    sample_seed: int,
    run_baselines: bool,
) -> pd.DataFrame:
    out_root.mkdir(parents=True, exist_ok=True)
    ref_df = load_charm_3dg(tdg_path)
    rows = [
        evaluate_config(
            pd.Series(row),
            ref_df,
            out_root,
            trans_sample_per_chrpair,
            sample_seed,
            baseline_trans_sample_per_chrpair,
            run_baselines=run_baselines,
        )
        for row in summary.to_dict("records")
    ]
    eval_summary = annotate_acceptability(pd.DataFrame(rows))
    eval_summary.to_csv(out_root / "eval_summary.tsv", sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    write_rankings(out_root, eval_summary)
    write_mechanism_tables(out_root, eval_summary)
    write_tsv(out_root / "eval_acceptability_warnings.tsv", warning_rows(eval_summary))
    return eval_summary


def select_existing_columns(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    present = [col for col in columns if col in frame.columns]
    return frame[present].copy()


def output_path(root: Path, prefix: str, suffix: str) -> Path:
    return root / f"{prefix}_{suffix}"


def finite_ratio(numerator: object, denominator: object) -> float:
    num = safe_float(numerator)
    den = safe_float(denominator)
    return num / den if math.isfinite(num) and math.isfinite(den) and den > 0.0 else math.nan


def finite_mean(values: Sequence[object]) -> float:
    arr = np.asarray([safe_float(value) for value in values], dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return float(arr.mean()) if arr.size else math.nan


def finite_median(values: Sequence[object]) -> float:
    arr = np.asarray([safe_float(value) for value in values], dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return float(np.median(arr)) if arr.size else math.nan


def bin_label(name: str, left: float, right: float) -> str:
    right_text = "inf" if math.isinf(right) else f"{right:g}"
    return f"{name}[{left:g},{right_text})"


def contact_force_r(r: float, k: float) -> float:
    if r < FDG_D_C1:
        return 2.0 * k * (FDG_D_C1 - r)
    if r <= FDG_D_C2:
        return 0.0
    if r <= FDG_D_C3:
        return -2.0 * k * (r - FDG_D_C2)
    t = r - FDG_D_C2
    return -k * (FDG_C_C1 - FDG_C_C2 / (t * t))


def contact_energy_r(r: float, k: float) -> float:
    if r < FDG_D_C1:
        return k * (FDG_D_C1 - r) ** 2
    if r <= FDG_D_C2:
        return 0.0
    if r <= FDG_D_C3:
        return k * (r - FDG_D_C2) ** 2
    t = r - FDG_D_C2
    return k * (FDG_C_C1 * (r - FDG_D_C3) + FDG_C_C2 / t)


def class_force_row(force: pd.DataFrame, row_class: str) -> Mapping[str, object]:
    if force.empty or "class" not in force.columns:
        return {}
    rows = force.loc[force["class"].astype(str) == row_class]
    return {} if rows.empty else rows.iloc[0].to_dict()


def coords_by_key(coords_path: Path) -> Dict[tuple[str, int, int], np.ndarray]:
    coord_lookup, _, _ = coords_lookup_and_index(coords_path)
    return coord_lookup


def coords_lookup_and_index(coords_path: Path) -> tuple[Dict[tuple[str, int, int], np.ndarray], Dict[tuple[str, int, int], int], int]:
    coords = load_model_coords(coords_path)
    coord_lookup: Dict[tuple[str, int, int], np.ndarray] = {}
    index_lookup: Dict[tuple[str, int, int], int] = {}
    for idx, row in enumerate(coords.itertuples(index=False)):
        key = (str(row.chr), int(row.start), int(row.copy))
        coord_lookup[key] = np.array([row.x, row.y, row.z], dtype=np.float64)
        index_lookup[key] = idx
    return coord_lookup, index_lookup, len(coords)


def diagnostic_rows_from_posterior(item: object) -> tuple[pd.DataFrame, str]:
    config_name = str(getattr(item, "config_name"))
    posterior_path = Path(str(getattr(item, "output_bpair_posterior", "")))
    coords_path = Path(str(getattr(item, "coords_gz", getattr(item, "output_coords", ""))))
    if not posterior_path.exists():
        return pd.DataFrame(), f"posterior_not_found:{posterior_path}"
    if not coords_path.exists():
        return pd.DataFrame(), f"coords_not_found:{coords_path}"
    posterior = read_tsv_trimmed(posterior_path)
    missing = sorted(POSTERIOR_REQUIRED_COLUMNS - set(posterior.columns))
    if missing:
        return pd.DataFrame(), f"posterior_missing_columns:{','.join(missing)}"
    try:
        coord_lookup = coords_by_key(coords_path)
    except Exception as exc:  # pragma: no cover - defensive diagnostics path
        return pd.DataFrame(), f"coords_load_failed:{exc}"

    rho_mode = str(getattr(item, "rho_train_mode", "constant"))
    rho_schedule = safe_float(getattr(item, "rho_train_end", getattr(item, "final_rho_train", 1.0)))
    rho_floor = safe_float(getattr(item, "rho_train_floor", getattr(item, "rho_floor", 0.5)))
    unit = safe_float(getattr(item, "unit", 1.0))
    if not math.isfinite(unit) or unit <= 0.0:
        unit = 1.0
    d_scale_mode = str(getattr(item, "d_scale_mode", "raw_count"))
    d_scale_eps_count = safe_float(getattr(item, "d_scale_eps_count", 1e-6))
    if not math.isfinite(d_scale_eps_count) or d_scale_eps_count <= 0.0:
        d_scale_eps_count = 1e-6
    cis_multiplier = safe_float(getattr(item, "contact_k_multiplier_cis", 1.0))
    trans_multiplier = safe_float(getattr(item, "contact_k_multiplier_trans", 1.0))
    if not math.isfinite(cis_multiplier) or cis_multiplier == 0.0:
        cis_multiplier = 1.0
    if not math.isfinite(trans_multiplier) or trans_multiplier == 0.0:
        trans_multiplier = 1.0

    rows: List[Dict[str, object]] = []
    for raw in posterior.to_dict("records"):
        chr1 = str(raw["chr1"])
        chr2 = str(raw["chr2"])
        start1 = int(safe_float(raw["start1"]))
        start2 = int(safe_float(raw["start2"]))
        contact_class = str(raw.get("contact_class") or ("cis" if chr1 == chr2 else "trans"))
        if contact_class not in {"cis", "trans"}:
            contact_class = "cis" if chr1 == chr2 else "trans"
        rho_output = safe_float(raw.get("rho_output"))
        if not math.isfinite(rho_output):
            entropy = safe_float(raw.get("entropy", LOG4))
            rho_output = max(0.0, min(1.0, 1.0 - entropy / LOG4)) if math.isfinite(entropy) else math.nan
        rho_train = rho_train_for_bpair(rho_mode, rho_schedule, rho_output, contact_class, rho_floor)
        multiplier = trans_multiplier if contact_class == "trans" else cis_multiplier
        base_k = safe_float(raw["base_k"])
        base_d_scale = safe_float(raw["base_d_scale"])
        n_raw = safe_float(raw["n_raw"])
        pvals = {state: safe_float(raw[col]) for state, _, _, col in CONTACT_STATES}
        distances: Dict[str, float] = {}
        force_vectors: Dict[str, np.ndarray] = {}
        edge_force_l1 = 0.0
        edge_energy = 0.0
        sum_wedge_k = 0.0
        n_wedges = 0
        missing_coord = False
        for state, copy1, copy2, _ in CONTACT_STATES:
            xyz1 = coord_lookup.get((chr1, start1, copy1))
            xyz2 = coord_lookup.get((chr2, start2, copy2))
            if xyz1 is None or xyz2 is None:
                missing_coord = True
                break
            delta = xyz1 - xyz2
            distance = float(np.linalg.norm(delta))
            direction = delta / distance if distance > 0.0 else np.array([1.0, 0.0, 0.0], dtype=np.float64)
            p = pvals[state]
            d_scale = base_d_scale
            if d_scale_mode == "expected_count":
                n_eff = n_raw * p if math.isfinite(n_raw) and math.isfinite(p) else math.nan
                if math.isfinite(n_eff):
                    d_scale = max(n_eff, d_scale_eps_count) ** (-1.0 / 3.0)
            k = base_k * multiplier * rho_train * p
            r = (distance / unit) / d_scale if math.isfinite(d_scale) and d_scale > 0.0 else math.nan
            fmag = contact_force_r(r, k) if math.isfinite(r) and math.isfinite(k) and k >= 0.0 else math.nan
            force = direction * fmag if math.isfinite(fmag) else np.full(3, math.nan)
            distances[state] = distance
            force_vectors[state] = force
            if math.isfinite(k):
                sum_wedge_k += k
            if math.isfinite(fmag):
                edge_force_l1 += 2.0 * abs(fmag) * float(np.abs(direction).sum())
                edge_energy += contact_energy_r(r, k)
            n_wedges += 1
        if missing_coord:
            continue
        p_array = np.array([pvals[state] for state, _, _, _ in CONTACT_STATES], dtype=np.float64)
        d_array = np.array([distances[state] for state, _, _, _ in CONTACT_STATES], dtype=np.float64)
        top_index = int(np.nanargmax(p_array))
        expected_distance = float(np.dot(p_array, d_array))
        net_vectors = [
            force_vectors["00"] + force_vectors["01"],
            force_vectors["10"] + force_vectors["11"],
            -force_vectors["00"] - force_vectors["10"],
            -force_vectors["01"] - force_vectors["11"],
        ]
        net_force_norm = float(sum(np.linalg.norm(vec) for vec in net_vectors))
        rows.append({
            "config_name": config_name,
            "chr1": chr1,
            "start1": start1,
            "chr2": chr2,
            "start2": start2,
            "chrpair": f"{chr1}:{chr2}" if chr1 <= chr2 else f"{chr2}:{chr1}",
            "contact_class": contact_class,
            "genomic_sep": abs(start2 - start1) if chr1 == chr2 else math.nan,
            "n_raw": n_raw,
            "base_d_scale": base_d_scale,
            "base_k": base_k,
            "pU": safe_float(raw["pU"]),
            "pmax": safe_float(raw["pmax"]),
            "rho_output": rho_output,
            "rho_train_bpair": rho_train,
            "sum_effective_k": sum_wedge_k,
            "n_wedges": n_wedges,
            "expected_distance": expected_distance,
            "top_distance": float(d_array[top_index]),
            "min_distance": float(np.min(d_array)),
            "expected_residual": expected_distance / base_d_scale if base_d_scale > 0.0 else math.nan,
            "top_residual": float(d_array[top_index]) / base_d_scale if base_d_scale > 0.0 else math.nan,
            "min_residual": float(np.min(d_array)) / base_d_scale if base_d_scale > 0.0 else math.nan,
            "contact_energy": edge_energy,
            "edge_endpoint_contact_force_l1": edge_force_l1,
            "net_force_norm": net_force_norm,
            "net_to_endpoint_force_l1_ratio": finite_ratio(net_force_norm, edge_force_l1),
        })
    if not rows:
        return pd.DataFrame(), "no_rows_with_coords"
    return pd.DataFrame(rows), ""


def write_posterior_outputs(root: Path, eval_summary: pd.DataFrame, prefix: str = "candidate_rerun") -> None:
    columns = [
        "config_name",
        "posterior_diag_available",
        "mean_pU_cis",
        "mean_pU_trans",
        "median_pU_cis",
        "median_pU_trans",
        "mean_pmax_cis",
        "mean_pmax_trans",
        "median_pmax_cis",
        "median_pmax_trans",
        "mean_margin_cis",
        "mean_margin_trans",
        "median_margin_cis",
        "median_margin_trans",
        "mean_psame_cis",
        "mean_psame_trans",
        "mean_pcross_cis",
        "mean_pcross_trans",
        "mean_rho_train_cis",
        "mean_rho_train_trans",
        "sum_effective_k_cis",
        "sum_effective_k_trans",
        "effective_k_cis_trans_ratio",
        "trans_fraction_effective_k",
        "fraction_trans_contacts_with_rho_train_near_zero",
    ]
    posterior_summary = select_existing_columns(eval_summary, columns)
    posterior_summary.to_csv(
        output_path(root, prefix, "posterior_cis_trans.tsv"),
        sep="\t",
        index=False,
        quoting=csv.QUOTE_MINIMAL,
    )
    posterior_summary.to_csv(
        output_path(root, prefix, "posterior_summary.tsv"),
        sep="\t",
        index=False,
        quoting=csv.QUOTE_MINIMAL,
    )
    hist_frames = []
    for eval_dir in eval_summary.get("eval_dir", []):
        path = Path(str(eval_dir)) / "eval.posterior_pu_hist.tsv"
        if path.exists():
            hist_frames.append(read_tsv_trimmed(path))
    if hist_frames:
        pd.concat(hist_frames, ignore_index=True).to_csv(output_path(root, prefix, "pu_hist.tsv"), sep="\t", index=False)


def force_value(frame: pd.DataFrame, row_class: str, column: str) -> float:
    if frame.empty or "class" not in frame.columns or column not in frame.columns:
        return math.nan
    rows = frame.loc[frame["class"].astype(str) == row_class]
    if rows.empty:
        return math.nan
    return safe_float(rows.iloc[0].get(column))


def write_force_output(root: Path, summary: pd.DataFrame, prefix: str = "candidate_rerun") -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for item in summary.itertuples(index=False):
        config_name = str(getattr(item, "config_name"))
        force_path = Path(str(getattr(item, "output_force_class_diag")))
        force = read_tsv_trimmed(force_path) if force_path.exists() else pd.DataFrame()
        trans_force = force_value(force, "trans", "contact_force_l1")
        cis_force = force_value(force, "cis", "contact_force_l1")
        repulsion_force = force_value(force, "total", "repulsion_force_l1")
        if not math.isfinite(repulsion_force):
            repulsion_force = safe_float(getattr(item, "final_repulsion_force_l1", math.nan))
        backbone_force = force_value(force, "total", "backbone_force_l1")
        if not math.isfinite(backbone_force):
            backbone_force = safe_float(getattr(item, "final_backbone_force_l1", math.nan))
        sep_force = force_value(force, "total", "homolog_sep_force_l1")
        rows.append({
            "config_name": config_name,
            "force_diag_path": str(force_path),
            "contact_energy_cis": force_value(force, "cis", "contact_energy"),
            "contact_energy_trans": force_value(force, "trans", "contact_energy"),
            "contact_force_l1_cis": cis_force,
            "contact_force_l1_trans": trans_force,
            "backbone_force_l1": backbone_force,
            "repulsion_force_l1": repulsion_force,
            "homolog_sep_force_l1": sep_force,
            "force_l1_cis_trans_ratio": cis_force / trans_force if math.isfinite(cis_force) and math.isfinite(trans_force) and trans_force > 0 else math.nan,
            "repulsion_to_trans_force_ratio": repulsion_force / trans_force if math.isfinite(repulsion_force) and math.isfinite(trans_force) and trans_force > 0 else math.nan,
        })
    out = pd.DataFrame(rows)
    out.to_csv(output_path(root, prefix, "force_diag.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    return out


def posterior_class_stats(eval_summary: pd.DataFrame) -> Dict[str, Dict[str, Dict[str, object]]]:
    out: Dict[str, Dict[str, Dict[str, object]]] = {}
    for row in eval_summary.to_dict("records"):
        config_name = str(row.get("config_name"))
        eval_dir = Path(str(row.get("eval_dir", "")))
        class_path = eval_dir / "eval.posterior_cis_trans_summary.tsv"
        if class_path.exists():
            frame = read_tsv_trimmed(class_path)
            out[config_name] = {str(item.get("class")): item for item in frame.to_dict("records")}
    return out


def write_force_normalized_output(
    root: Path,
    summary: pd.DataFrame,
    eval_summary: pd.DataFrame,
    prefix: str = "candidate_rerun",
) -> pd.DataFrame:
    posterior_stats = posterior_class_stats(eval_summary)
    eval_by_config = {str(row["config_name"]): row for row in eval_summary.to_dict("records")}
    rows: List[Dict[str, object]] = []
    for item in summary.itertuples(index=False):
        config_name = str(getattr(item, "config_name"))
        force_path = Path(str(getattr(item, "output_force_class_diag")))
        force = read_tsv_trimmed(force_path) if force_path.exists() else pd.DataFrame()
        cis_force = class_force_row(force, "cis")
        trans_force = class_force_row(force, "trans")
        class_stats = posterior_stats.get(config_name, {})
        eval_row = eval_by_config.get(config_name, {})
        row: Dict[str, object] = {"config_name": config_name}
        same_bin_bpair_excluded = safe_float(getattr(item, "n_bpair_same_bin_excluded", 0.0))
        same_bin_raw_excluded = safe_float(getattr(item, "n_raw_same_bin_excluded", 0.0))
        if not math.isfinite(same_bin_bpair_excluded):
            same_bin_bpair_excluded = 0.0
        if not math.isfinite(same_bin_raw_excluded):
            same_bin_raw_excluded = 0.0
        for label, force_row in [("cis", cis_force), ("trans", trans_force)]:
            stats = class_stats.get(label, {})
            n_bpair_all = safe_float(stats.get("n_bpair", getattr(item, f"n_bpair_{label}", math.nan)))
            n_raw_all = safe_float(stats.get("n_raw_sum", getattr(item, f"n_raw_{label}", math.nan)))
            n_bpair = n_bpair_all
            n_raw = n_raw_all
            if label == "cis":
                if math.isfinite(n_bpair):
                    n_bpair = max(0.0, n_bpair - same_bin_bpair_excluded)
                if math.isfinite(n_raw):
                    n_raw = max(0.0, n_raw - same_bin_raw_excluded)
            n_wedge = safe_float(force_row.get("n_wedges", math.nan))
            sum_effective_k = safe_float(
                force_row.get("sum_wedge_k", stats.get("sum_effective_k", eval_row.get(f"sum_effective_k_{label}", math.nan)))
            )
            sum_effective_k_posterior = safe_float(stats.get("sum_effective_k", eval_row.get(f"sum_effective_k_{label}", math.nan)))
            force_l1 = safe_float(force_row.get("contact_force_l1", eval_row.get(f"final_contact_force_l1_{label}", math.nan)))
            row.update({
                f"n_bpair_{label}": n_bpair,
                f"n_raw_{label}": n_raw,
                f"n_bpair_posterior_{label}": n_bpair_all,
                f"n_raw_posterior_{label}": n_raw_all,
                f"n_wedge_{label}": n_wedge,
                f"sum_effective_k_{label}": sum_effective_k,
                f"sum_effective_k_posterior_{label}": sum_effective_k_posterior,
                f"mean_effective_k_per_bpair_{label}": finite_ratio(sum_effective_k, n_bpair),
                f"mean_effective_k_per_raw_{label}": finite_ratio(sum_effective_k, n_raw),
                f"mean_effective_k_per_wedge_{label}": finite_ratio(sum_effective_k, n_wedge),
                f"force_l1_{label}": force_l1,
                f"force_l1_per_bpair_{label}": finite_ratio(force_l1, n_bpair),
                f"force_l1_per_raw_{label}": finite_ratio(force_l1, n_raw),
                f"force_l1_per_wedge_{label}": finite_ratio(force_l1, n_wedge),
            })
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(output_path(root, prefix, "force_normalized.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    return out


def force_cancellation_from_posterior(item: object) -> Dict[str, object]:
    config_name = str(getattr(item, "config_name"))
    posterior_path = Path(str(getattr(item, "output_bpair_posterior", "")))
    coords_path = Path(str(getattr(item, "coords_gz", getattr(item, "output_coords", ""))))
    empty: Dict[str, object] = {
        "config_name": config_name,
        "available": 0,
        "reason": "",
    }
    if not posterior_path.exists():
        empty["reason"] = f"posterior_not_found:{posterior_path}"
        return empty
    if not coords_path.exists():
        empty["reason"] = f"coords_not_found:{coords_path}"
        return empty
    posterior = read_tsv_trimmed(posterior_path)
    missing = sorted(POSTERIOR_REQUIRED_COLUMNS - set(posterior.columns))
    if missing:
        empty["reason"] = f"posterior_missing_columns:{','.join(missing)}"
        return empty
    try:
        coord_lookup, index_lookup, n_coords = coords_lookup_and_index(coords_path)
    except Exception as exc:  # pragma: no cover - defensive diagnostics path
        empty["reason"] = f"coords_load_failed:{exc}"
        return empty

    rho_mode = str(getattr(item, "rho_train_mode", "constant"))
    rho_schedule = safe_float(getattr(item, "rho_train_end", getattr(item, "final_rho_train", 1.0)))
    rho_floor = safe_float(getattr(item, "rho_train_floor", getattr(item, "rho_floor", 0.5)))
    unit = safe_float(getattr(item, "unit", 1.0))
    if not math.isfinite(unit) or unit <= 0.0:
        unit = 1.0
    d_scale_mode = str(getattr(item, "d_scale_mode", "raw_count"))
    d_scale_eps_count = safe_float(getattr(item, "d_scale_eps_count", 1e-6))
    if not math.isfinite(d_scale_eps_count) or d_scale_eps_count <= 0.0:
        d_scale_eps_count = 1e-6
    cis_multiplier = safe_float(getattr(item, "contact_k_multiplier_cis", 1.0))
    trans_multiplier = safe_float(getattr(item, "contact_k_multiplier_trans", 1.0))
    if not math.isfinite(cis_multiplier) or cis_multiplier == 0.0:
        cis_multiplier = 1.0
    if not math.isfinite(trans_multiplier) or trans_multiplier == 0.0:
        trans_multiplier = 1.0

    forces = {
        "cis": np.zeros((n_coords, 3), dtype=np.float64),
        "trans": np.zeros((n_coords, 3), dtype=np.float64),
    }
    sums: Dict[str, Dict[str, float]] = {
        "cis": {"n_bpair": 0.0, "n_wedge": 0.0, "sum_effective_k": 0.0, "force_l1": 0.0, "contact_energy": 0.0, "n_nonfinite": 0.0},
        "trans": {"n_bpair": 0.0, "n_wedge": 0.0, "sum_effective_k": 0.0, "force_l1": 0.0, "contact_energy": 0.0, "n_nonfinite": 0.0},
    }
    skipped_same_bin = 0
    skipped_missing_coords = 0

    for raw in posterior.to_dict("records"):
        chr1 = str(raw["chr1"])
        chr2 = str(raw["chr2"])
        start1 = int(safe_float(raw["start1"]))
        start2 = int(safe_float(raw["start2"]))
        bid1 = int(safe_float(raw["bid1"]))
        bid2 = int(safe_float(raw["bid2"]))
        if bid1 == bid2:
            skipped_same_bin += 1
            continue
        contact_class = str(raw.get("contact_class") or ("cis" if chr1 == chr2 else "trans"))
        if contact_class not in {"cis", "trans"}:
            contact_class = "cis" if chr1 == chr2 else "trans"
        rho_output = safe_float(raw.get("rho_output"))
        if not math.isfinite(rho_output):
            entropy = safe_float(raw.get("entropy", LOG4))
            rho_output = max(0.0, min(1.0, 1.0 - entropy / LOG4)) if math.isfinite(entropy) else math.nan
        rho_train = rho_train_for_bpair(rho_mode, rho_schedule, rho_output, contact_class, rho_floor)
        multiplier = trans_multiplier if contact_class == "trans" else cis_multiplier
        base_k = safe_float(raw["base_k"])
        base_d_scale = safe_float(raw["base_d_scale"])
        n_raw = safe_float(raw["n_raw"])
        pvals = {state: safe_float(raw[col]) for state, _, _, col in CONTACT_STATES}

        bpair_has_force = False
        for state, copy1, copy2, _ in CONTACT_STATES:
            key1 = (chr1, start1, copy1)
            key2 = (chr2, start2, copy2)
            xyz1 = coord_lookup.get(key1)
            xyz2 = coord_lookup.get(key2)
            idx1 = index_lookup.get(key1)
            idx2 = index_lookup.get(key2)
            if xyz1 is None or xyz2 is None or idx1 is None or idx2 is None:
                skipped_missing_coords += 1
                continue
            if idx1 == idx2:
                continue
            p = pvals[state]
            d_scale = base_d_scale
            if d_scale_mode == "expected_count":
                n_eff = n_raw * p if math.isfinite(n_raw) and math.isfinite(p) else math.nan
                if math.isfinite(n_eff):
                    d_scale = max(n_eff, d_scale_eps_count) ** (-1.0 / 3.0)
            k = base_k * multiplier * rho_train * p
            delta = xyz1 - xyz2
            distance = float(np.linalg.norm(delta))
            direction = delta / distance if distance > 0.0 else np.array([1.0, 0.0, 0.0], dtype=np.float64)
            r = (distance / unit) / d_scale if math.isfinite(d_scale) and d_scale > 0.0 else math.nan
            fmag = contact_force_r(r, k) if math.isfinite(r) and math.isfinite(k) and k >= 0.0 else math.nan
            if not math.isfinite(fmag):
                sums[contact_class]["n_nonfinite"] += 1.0
                continue
            force = direction * fmag
            forces[contact_class][idx1] += force
            forces[contact_class][idx2] -= force
            sums[contact_class]["n_wedge"] += 1.0
            sums[contact_class]["sum_effective_k"] += k if math.isfinite(k) else 0.0
            sums[contact_class]["force_l1"] += 2.0 * abs(fmag) * float(np.abs(direction).sum())
            sums[contact_class]["contact_energy"] += contact_energy_r(r, k)
            bpair_has_force = True
        if bpair_has_force:
            sums[contact_class]["n_bpair"] += 1.0

    if not any(sums[label]["n_wedge"] > 0.0 for label in ("cis", "trans")):
        empty["reason"] = "no_force_rows_with_coords"
        empty["skipped_same_bin_bpairs"] = skipped_same_bin
        empty["skipped_missing_state_coords"] = skipped_missing_coords
        return empty

    all_force = forces["cis"] + forces["trans"]
    row: Dict[str, object] = {
        "config_name": config_name,
        "available": 1,
        "reason": "",
        "skipped_same_bin_bpairs": skipped_same_bin,
        "skipped_missing_state_coords": skipped_missing_coords,
    }
    all_force_l1 = 0.0
    all_net_norm = float(np.linalg.norm(all_force, axis=1).sum())
    for label in ("cis", "trans"):
        net_norm = float(np.linalg.norm(forces[label], axis=1).sum())
        force_l1 = sums[label]["force_l1"]
        all_force_l1 += force_l1
        row.update({
            f"n_force_bpair_{label}": int(sums[label]["n_bpair"]),
            f"n_force_wedge_{label}": int(sums[label]["n_wedge"]),
            f"sum_effective_k_{label}": sums[label]["sum_effective_k"],
            f"contact_energy_{label}": sums[label]["contact_energy"],
            f"force_l1_{label}": force_l1,
            f"net_force_norm_{label}": net_norm,
            f"force_cancellation_ratio_{label}": finite_ratio(net_norm, force_l1),
            f"n_nonfinite_state_force_{label}": int(sums[label]["n_nonfinite"]),
        })
    row.update({
        "force_l1_all": all_force_l1,
        "net_force_norm_all": all_net_norm,
        "force_cancellation_ratio_all": finite_ratio(all_net_norm, all_force_l1),
    })
    return row


def write_force_cancellation_output(root: Path, summary: pd.DataFrame, prefix: str = "candidate_rerun") -> pd.DataFrame:
    rows = [force_cancellation_from_posterior(item) for item in summary.itertuples(index=False)]
    out = pd.DataFrame(rows)
    out.to_csv(output_path(root, prefix, "force_cancellation.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    return out


def residual_stratum_row(config_name: str, stratification: str, stratum: str, frame: pd.DataFrame) -> Dict[str, object]:
    if frame.empty:
        return {
            "config_name": config_name,
            "stratification": stratification,
            "stratum": stratum,
            "available": 0,
            "n_bpair": 0,
        }
    return {
        "config_name": config_name,
        "stratification": stratification,
        "stratum": stratum,
        "available": 1,
        "n_bpair": int(len(frame)),
        "n_raw_sum": float(frame["n_raw"].sum()),
        "mean_n_raw": float(frame["n_raw"].mean()),
        "mean_pU": float(frame["pU"].mean()),
        "median_pU": float(frame["pU"].median()),
        "mean_pmax": float(frame["pmax"].mean()),
        "median_pmax": float(frame["pmax"].median()),
        "mean_d_scale": float(frame["base_d_scale"].mean()),
        "median_expected_residual": float(frame["expected_residual"].median()),
        "median_top_residual": float(frame["top_residual"].median()),
        "median_min_residual": float(frame["min_residual"].median()),
        "sum_effective_k": float(frame["sum_effective_k"].sum()),
        "contact_count_vs_inv_expected_distance_spearman": spearman_corr(
            np.log1p(frame["n_raw"].to_numpy(float)),
            -frame["expected_distance"].to_numpy(float),
        ),
    }


def append_bin_rows(
    rows: List[Dict[str, object]],
    config_name: str,
    diag: pd.DataFrame,
    stratification: str,
    column: str,
    bins: Sequence[tuple[float, float]],
) -> None:
    values = diag[column].to_numpy(float)
    for left, right in bins:
        mask = values >= left
        mask &= values < right if math.isfinite(right) else values >= left
        rows.append(residual_stratum_row(config_name, stratification, bin_label(column, left, right), diag.loc[mask]))


def state_d_scale(base_d_scale: float, n_raw: float, state_p: float, d_scale_mode: str, eps_count: float) -> tuple[float, int]:
    if d_scale_mode != "expected_count":
        return base_d_scale, 0
    n_eff = n_raw * state_p if math.isfinite(n_raw) and math.isfinite(state_p) else math.nan
    if not math.isfinite(n_eff):
        return math.nan, 0
    clipped = int(n_eff < eps_count)
    return max(n_eff, eps_count) ** (-1.0 / 3.0), clipped


def dscale_distribution_row(
    config_name: str,
    d_scale_mode: str,
    stratification: str,
    stratum: str,
    frame: pd.DataFrame,
) -> Dict[str, object]:
    if frame.empty:
        return {
            "config_name": config_name,
            "d_scale_mode": d_scale_mode,
            "stratification": stratification,
            "stratum": stratum,
            "available": 0,
            "n_state": 0,
        }
    weights = frame["state_p"].to_numpy(float)
    dscale = frame["state_d_scale"].to_numpy(float)
    finite_weight = np.isfinite(weights) & np.isfinite(dscale)
    weighted = float(np.average(dscale[finite_weight], weights=weights[finite_weight])) if finite_weight.any() and float(weights[finite_weight].sum()) > 0.0 else math.nan
    return {
        "config_name": config_name,
        "d_scale_mode": d_scale_mode,
        "stratification": stratification,
        "stratum": stratum,
        "available": 1,
        "n_state": int(len(frame)),
        "n_bpair_approx": float(len(frame)) / 4.0,
        "sum_state_p": float(frame["state_p"].sum()),
        "mean_state_p": float(frame["state_p"].mean()),
        "median_state_p": float(frame["state_p"].median()),
        "mean_n_raw": float(frame["n_raw"].mean()),
        "median_n_raw": float(frame["n_raw"].median()),
        "mean_base_d_scale": float(frame["base_d_scale"].mean()),
        "median_base_d_scale": float(frame["base_d_scale"].median()),
        "mean_state_d_scale": float(frame["state_d_scale"].mean()),
        "median_state_d_scale": float(frame["state_d_scale"].median()),
        "p10_state_d_scale": float(frame["state_d_scale"].quantile(0.10)),
        "p90_state_d_scale": float(frame["state_d_scale"].quantile(0.90)),
        "weighted_mean_state_d_scale_by_p": weighted,
        "fraction_clipped_eps_count": float(frame["d_scale_clipped"].mean()),
    }


def dscale_state_table(item: object) -> pd.DataFrame:
    posterior_path = Path(str(getattr(item, "output_bpair_posterior", "")))
    if not posterior_path.exists():
        return pd.DataFrame()
    posterior = read_tsv_trimmed(posterior_path)
    missing = sorted(POSTERIOR_REQUIRED_COLUMNS - set(posterior.columns))
    if missing:
        return pd.DataFrame()
    d_scale_mode = str(getattr(item, "d_scale_mode", "raw_count"))
    eps_count = safe_float(getattr(item, "d_scale_eps_count", 1e-6))
    if not math.isfinite(eps_count) or eps_count <= 0.0:
        eps_count = 1e-6

    rows: List[Dict[str, object]] = []
    for raw in posterior.to_dict("records"):
        chr1 = str(raw["chr1"])
        chr2 = str(raw["chr2"])
        contact_class = str(raw.get("contact_class") or ("cis" if chr1 == chr2 else "trans"))
        if contact_class not in {"cis", "trans"}:
            contact_class = "cis" if chr1 == chr2 else "trans"
        n_raw = safe_float(raw["n_raw"])
        base_d_scale = safe_float(raw["base_d_scale"])
        pU = safe_float(raw["pU"])
        pmax = safe_float(raw["pmax"])
        for state, _, _, column in CONTACT_STATES:
            state_p = safe_float(raw[column])
            dscale, clipped = state_d_scale(base_d_scale, n_raw, state_p, d_scale_mode, eps_count)
            if not math.isfinite(dscale):
                continue
            rows.append({
                "contact_class": contact_class,
                "chrpair": f"{chr1}:{chr2}" if chr1 <= chr2 else f"{chr2}:{chr1}",
                "state": state,
                "n_raw": n_raw,
                "base_d_scale": base_d_scale,
                "pU": pU,
                "pmax": pmax,
                "state_p": state_p,
                "state_d_scale": dscale,
                "d_scale_clipped": clipped,
            })
    return pd.DataFrame(rows)


def write_dscale_diag_output(root: Path, summary: pd.DataFrame, prefix: str = "dscale_ablation") -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for item in summary.itertuples(index=False):
        config_name = str(getattr(item, "config_name"))
        d_scale_mode = str(getattr(item, "d_scale_mode", "raw_count"))
        table = dscale_state_table(item)
        if table.empty:
            rows.append({
                "config_name": config_name,
                "d_scale_mode": d_scale_mode,
                "stratification": "unavailable",
                "stratum": "all",
                "available": 0,
                "n_state": 0,
            })
            continue
        rows.append(dscale_distribution_row(config_name, d_scale_mode, "all", "all", table))
        for label, group in table.groupby("contact_class", sort=False):
            rows.append(dscale_distribution_row(config_name, d_scale_mode, "contact_class", str(label), group))
        for left, right in N_RAW_BINS:
            mask = table["n_raw"].to_numpy(float) >= left
            mask &= table["n_raw"].to_numpy(float) < right if math.isfinite(right) else table["n_raw"].to_numpy(float) >= left
            rows.append(dscale_distribution_row(config_name, d_scale_mode, "n_raw_bin", bin_label("n_raw", left, right), table.loc[mask]))
        for left, right in PU_BINS:
            mask = table["pU"].to_numpy(float) >= left
            mask &= table["pU"].to_numpy(float) < right if math.isfinite(right) else table["pU"].to_numpy(float) >= left
            rows.append(dscale_distribution_row(config_name, d_scale_mode, "pU_bin", bin_label("pU", left, right), table.loc[mask]))
        for left, right in PMAX_BINS:
            mask = table["pmax"].to_numpy(float) >= left
            mask &= table["pmax"].to_numpy(float) < right if math.isfinite(right) else table["pmax"].to_numpy(float) >= left
            rows.append(dscale_distribution_row(config_name, d_scale_mode, "posterior_confidence_bin", bin_label("pmax", left, right), table.loc[mask]))
        for left, right in STATE_P_BINS:
            mask = table["state_p"].to_numpy(float) >= left
            mask &= table["state_p"].to_numpy(float) < right if math.isfinite(right) else table["state_p"].to_numpy(float) >= left
            rows.append(dscale_distribution_row(config_name, d_scale_mode, "state_p_bin", bin_label("state_p", left, right), table.loc[mask]))
        for left, right in D_SCALE_BINS:
            mask = table["state_d_scale"].to_numpy(float) >= left
            mask &= table["state_d_scale"].to_numpy(float) < right if math.isfinite(right) else table["state_d_scale"].to_numpy(float) >= left
            rows.append(dscale_distribution_row(config_name, d_scale_mode, "state_d_scale_bin", bin_label("state_d_scale", left, right), table.loc[mask]))
    out = pd.DataFrame(rows)
    out.to_csv(output_path(root, prefix, "dscale_diag.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    return out


def write_residual_stratified_output(root: Path, summary: pd.DataFrame, prefix: str = "candidate_rerun") -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for item in summary.itertuples(index=False):
        config_name = str(getattr(item, "config_name"))
        diag, reason = diagnostic_rows_from_posterior(item)
        if diag.empty:
            rows.append({
                "config_name": config_name,
                "stratification": "unavailable",
                "stratum": "all",
                "available": 0,
                "reason": reason,
                "n_bpair": 0,
            })
            continue
        rows.append(residual_stratum_row(config_name, "contact_class", "cis_short_range", diag.loc[(diag["contact_class"] == "cis") & (diag["genomic_sep"] <= 20_000_000)]))
        rows.append(residual_stratum_row(config_name, "contact_class", "cis_long_range", diag.loc[(diag["contact_class"] == "cis") & (diag["genomic_sep"] > 20_000_000)]))
        rows.append(residual_stratum_row(config_name, "contact_class", "trans", diag.loc[diag["contact_class"] == "trans"]))
        for chrpair, group in diag.loc[diag["contact_class"] == "trans"].groupby("chrpair", sort=False):
            rows.append(residual_stratum_row(config_name, "trans_chrpair", str(chrpair), group))
        append_bin_rows(rows, config_name, diag, "n_raw_bin", "n_raw", N_RAW_BINS)
        append_bin_rows(rows, config_name, diag, "pU_bin", "pU", PU_BINS)
        append_bin_rows(rows, config_name, diag, "pmax_bin", "pmax", PMAX_BINS)
        append_bin_rows(rows, config_name, diag, "d_scale_bin", "base_d_scale", D_SCALE_BINS)
    out = pd.DataFrame(rows)
    out.to_csv(output_path(root, prefix, "residual_stratified.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    return out


def write_contact_residual_output(root: Path, eval_summary: pd.DataFrame, prefix: str = "candidate_rerun") -> None:
    frames = []
    summary_columns = [
        "config_name",
        "cis_contact_count_vs_inv_distance_spearman",
        "trans_contact_count_vs_inv_distance_spearman",
        "cis_contact_residual_median",
        "trans_contact_residual_median",
        "cis_high_count_short_expected_precision_at_5pct",
        "trans_high_count_short_expected_precision_at_5pct",
    ]
    summary = select_existing_columns(eval_summary, summary_columns)
    if not summary.empty:
        summary.insert(0, "row_type", "summary")
        frames.append(summary)
    for eval_dir in eval_summary.get("eval_dir", []):
        eval_path = Path(str(eval_dir))
        class_path = eval_path / "eval.contact_distance_summary.tsv"
        chrpair_path = eval_path / "eval.trans_contact_by_chrpair.tsv"
        if class_path.exists():
            frame = read_tsv_trimmed(class_path)
            frame.insert(0, "row_type", "class")
            frames.append(frame)
        if chrpair_path.exists():
            frame = read_tsv_trimmed(chrpair_path)
            frame.insert(0, "row_type", "trans_chrpair")
            frames.append(frame)
    if frames:
        pd.concat(frames, ignore_index=True, sort=False).to_csv(
            output_path(root, prefix, "contact_residual.tsv"),
            sep="\t",
            index=False,
            quoting=csv.QUOTE_MINIMAL,
        )


def fmt(value: object, digits: int = 3) -> str:
    val = safe_float(value)
    return f"{val:.{digits}g}" if math.isfinite(val) else "NA"


def pass_label(value: object) -> str:
    val = safe_float(value)
    return "PASS" if math.isfinite(val) and int(val) == 1 else "FAIL"


def write_interpretation(root: Path, eval_summary: pd.DataFrame, force_summary: pd.DataFrame, prefix: str = "candidate_rerun") -> None:
    force_by_config = {str(row["config_name"]): row for row in force_summary.to_dict("records")}
    title = "P9016 prior/rho 2x2 interpretation" if prefix == "prior_rho_2x2" else "P9016 full-output candidate rerun interpretation"
    lines = [
        f"# {title}",
        "",
        "This is post-hoc eval-only analysis. CHARM/3DG is used only after rerun outputs exist, and no configuration is selected as a best model by CHARM.",
        "copy0/copy1 labels are gauge labels, not maternal/paternal haplotypes.",
        "",
    ]
    if prefix == "prior_rho_2x2":
        lines.extend([
            "## 2x2 Decision Frame",
            "",
            "The baseline is `cis_inter_ratio + entropy`. `uniform + entropy` isolates the prior effect. `cis_inter_ratio + constant` isolates entropy rho_train. `uniform + constant` tests whether a flat 0.25 prior with full scheduled contact force is sufficient.",
            "",
            "Interpretation rules: uniform prior improving scale while damaging cis means the cis/inter prior is informative but too strong; constant rho improving trans residual without damaging cis means entropy rho_train contributes to failure; uniformly poor `uniform + constant` means flat priors alone are not a solution.",
            "",
        ])
    lines.extend([
        "## Candidate summary",
        "",
    ])
    for row in eval_summary.to_dict("records"):
        name = str(row["config_name"])
        force = force_by_config.get(name, {})
        underweighted = (
            safe_float(row.get("mean_pU_trans")) > 0.8
            and safe_float(row.get("mean_rho_train_trans")) < 0.3
            and safe_float(row.get("trans_fraction_effective_k")) < 0.2
        )
        force_driven = (
            safe_float(force.get("force_l1_cis_trans_ratio")) > 5.0
            and safe_float(row.get("median_rg_ratio")) < 0.7
        )
        negative = name in {"main_rep1_n5_rs80", "main_rep1_n5_rs5", "random_rep1_n120_rs50"}
        if negative:
            role = "diagnostic negative/control, not an absolute-pass trans candidate"
        elif pass_label(row.get("trans_metric_absolute_pass")) == "PASS" and pass_label(row.get("centroid_metric_absolute_pass")) == "PASS":
            role = "valid Pareto follow-up candidate by absolute flags"
        else:
            role = "cis/compaction mechanism candidate, not an absolute trans candidate"
        lines.extend([
            f"### {name}",
            "",
            f"- cis quality: chr1 Pearson {fmt(row.get('chr1_cis_min_copy_pearson'))}, all-chr mean Pearson {fmt(row.get('all_chr_cis_min_copy_pearson_mean'))}, pooled Spearman {fmt(row.get('all_chr_cis_pooled_spearman'))}.",
            f"- compaction severity: median Rg ratio {fmt(row.get('median_rg_ratio'))}, median cis distance slope {fmt(row.get('median_cis_distance_slope'))}.",
            f"- trans/centroid absolute status: trans {pass_label(row.get('trans_metric_absolute_pass'))}; centroid {pass_label(row.get('centroid_metric_absolute_pass'))}.",
            f"- posterior/rho evidence: mean_pU_trans {fmt(row.get('mean_pU_trans'))}, mean_rho_train_trans {fmt(row.get('mean_rho_train_trans'))}, trans effective-k fraction {fmt(row.get('trans_fraction_effective_k'))}, near-zero trans rho fraction {fmt(row.get('fraction_trans_contacts_with_rho_train_near_zero'))}. Underweighting pattern: {'yes' if underweighted else 'not proven'}.",
            f"- force evidence: cis/trans contact-force ratio {fmt(force.get('force_l1_cis_trans_ratio'))}, repulsion/trans-force ratio {fmt(force.get('repulsion_to_trans_force_ratio'))}, backbone force {fmt(force.get('backbone_force_l1'))}. Force-driven compaction pattern: {'yes' if force_driven else 'not proven'}.",
            f"- contact residual evidence: cis residual median {fmt(row.get('cis_contact_residual_median'))}, trans residual median {fmt(row.get('trans_contact_residual_median'))}, trans count-vs-inv-distance Spearman {fmt(row.get('trans_contact_count_vs_inv_distance_spearman'))}.",
            f"- interpretation: {role}.",
            "",
        ])
    (root / f"{prefix}_interpretation.md").write_text("\n".join(lines), encoding="utf-8")


def write_mechanism_update(
    root: Path,
    eval_summary: pd.DataFrame,
    force_normalized: pd.DataFrame,
    force_cancellation: pd.DataFrame,
    residual_stratified: pd.DataFrame,
    prefix: str = "candidate_rerun",
) -> None:
    norm_by_config = {str(row["config_name"]): row for row in force_normalized.to_dict("records")}
    cancel_by_config = {str(row["config_name"]): row for row in force_cancellation.to_dict("records")} if not force_cancellation.empty else {}
    residual_available = int(residual_stratified.get("available", pd.Series(dtype=int)).fillna(0).astype(int).sum()) if not residual_stratified.empty else 0
    fn = force_normalized
    fc = force_cancellation
    high_cis_names = set(eval_summary.loc[
        (eval_summary.get("all_chr_cis_min_copy_pearson_mean", pd.Series(dtype=float)).map(safe_float) >= 0.35)
        & (eval_summary.get("median_rg_ratio", pd.Series(dtype=float)).map(safe_float) < 0.7),
        "config_name",
    ].astype(str)) if not eval_summary.empty else set()
    high_cis_norm = fn.loc[fn["config_name"].astype(str).isin(high_cis_names)] if not fn.empty and high_cis_names else pd.DataFrame()
    high_cis_cancel = fc.loc[fc["config_name"].astype(str).isin(high_cis_names)] if not fc.empty and high_cis_names else pd.DataFrame()
    high_cis_trans_force = finite_median(high_cis_norm["force_l1_per_wedge_trans"]) if not high_cis_norm.empty else math.nan
    high_cis_cis_force = finite_median(high_cis_norm["force_l1_per_wedge_cis"]) if not high_cis_norm.empty else math.nan
    high_cis_trans_cancel = finite_median(high_cis_cancel["force_cancellation_ratio_trans"]) if not high_cis_cancel.empty else math.nan
    high_cis_trans_residual = finite_median(eval_summary.loc[
        eval_summary["config_name"].astype(str).isin(high_cis_names),
        "trans_contact_residual_median",
    ]) if high_cis_names and "trans_contact_residual_median" in eval_summary.columns else math.nan
    lines = [
        "# Mechanism Update",
        "",
        "This update is postprocess-only: it reads existing full-output, posterior, coordinate, force-diagnostic, and eval files. It does not invoke training, config selection, or any CHARM-dependent rerun logic.",
        "",
        "## Normalized Force Interpretation",
        "",
        "Raw class force L1 is dominated by how many contacts and expanded state wedges each class contributes. The normalized table therefore reports force per bpair, raw count, and wedge alongside effective-k per bpair/raw/wedge. The bpair/raw denominators subtract same-bin cis rows because the C force wedge builder skips same-bin bpairs; posterior totals are retained in separate `*_posterior_*` columns.",
        "",
        "Conceptually, a class with 10x more wedges can have 10x larger force even if each wedge is ordinary. Per-wedge force separates the per-contact mechanical scale from the bookkeeping scale.",
        "",
        "## Cancellation Interpretation",
        "",
        "The cancellation table reconstructs posterior-weighted contact forces from posterior probabilities and final coordinates, skips same-bin bpairs to match the force wedge support, accumulates class-specific vectors on each diploid bead, and then measures bead-level L2 norm. `force_cancellation_ratio_* = net_force_norm_* / force_l1_*`; a small ratio means large endpoint force is cancelling across constraints or posterior states before it becomes coherent bead motion.",
        "",
        "## Residual Stratification",
        "",
        "The residual stratified table partitions reconstructed contact residuals by cis short range, cis long range, trans, individual trans chromosome pair, raw-count bins, posterior uncertainty bins, pmax bins, and d_scale bins. This checks whether residual patterns are broad or localized to a subset such as uncertain trans contacts.",
        "",
        f"Residual strata with available rows: {residual_available}.",
        "",
        "## Cross-Run Readout",
        "",
        f"- High-cis compressed configs: {', '.join(sorted(high_cis_names)) if high_cis_names else 'none under the current thresholds'}.",
        f"- Median per-wedge force in those configs, cis/trans: {fmt(high_cis_cis_force)} / {fmt(high_cis_trans_force)}.",
        f"- Median trans cancellation ratio in those configs: {fmt(high_cis_trans_cancel)}.",
        f"- Median trans top-state residual in those configs: {fmt(high_cis_trans_residual)}.",
        "- Mechanistic read: aggregate trans force is not absent. When trans force remains larger than cis after bpair/raw/wedge normalization but trans residual and absolute trans flags remain poor, the failure is better interpreted as incoherent/noisy or scale-incompatible trans constraints rather than simple trans underweighting.",
        "- Scale read: when cis short-range residual is better than trans residual but Rg ratio and cis distance slope remain low, the current objective is still compatible with compressed structures; this points to scale/anti-collapse/d_scale diagnostics rather than simply amplifying trans force.",
        "",
        "## Candidate-Level Notes",
        "",
    ]
    for row in eval_summary.to_dict("records"):
        name = str(row["config_name"])
        norm = norm_by_config.get(name, {})
        cancel = cancel_by_config.get(name, {})
        lines.extend([
            f"### {name}",
            "",
            f"- cis force per bpair/raw/wedge: {fmt(norm.get('force_l1_per_bpair_cis'))} / {fmt(norm.get('force_l1_per_raw_cis'))} / {fmt(norm.get('force_l1_per_wedge_cis'))}.",
            f"- trans force per bpair/raw/wedge: {fmt(norm.get('force_l1_per_bpair_trans'))} / {fmt(norm.get('force_l1_per_raw_trans'))} / {fmt(norm.get('force_l1_per_wedge_trans'))}.",
            f"- effective-k per bpair cis/trans: {fmt(norm.get('mean_effective_k_per_bpair_cis'))} / {fmt(norm.get('mean_effective_k_per_bpair_trans'))}.",
            f"- reconstructed cancellation ratio cis/trans/all: {fmt(cancel.get('force_cancellation_ratio_cis'))} / {fmt(cancel.get('force_cancellation_ratio_trans'))} / {fmt(cancel.get('force_cancellation_ratio_all'))}.",
            f"- posterior context: mean_pU_trans {fmt(row.get('mean_pU_trans'))}, trans effective-k fraction {fmt(row.get('trans_fraction_effective_k'))}.",
            "",
        ])
    (root / f"{prefix}_mechanism_update.md").write_text("\n".join(lines), encoding="utf-8")


def scaffold_stage_summary_rows(summary: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for raw in summary.to_dict("records"):
        source_name = str(raw["config_name"])
        source_dir = Path(str(raw["output_dir"]))
        for stage, filename, use_final_diagnostics in SCAFFOLD_STAGE_FILES:
            coords_path = source_dir / filename
            if not coords_path.exists():
                raise FileNotFoundError(coords_path)
            row = dict(raw)
            row["source_config_name"] = source_name
            row["stage"] = stage
            row["config_name"] = f"{source_name}__{stage}"
            row["coords_gz"] = str(coords_path)
            row["output_coords"] = str(coords_path)
            row["output_coords_gz"] = str(coords_path)
            if not use_final_diagnostics:
                # Initial and one-iteration snapshots intentionally do not
                # reuse the final posterior or force file; otherwise residual
                # diagnostics would mix stages.
                row["output_dir"] = ""
                row["output_bpair_posterior"] = ""
                row["output_force_class_diag"] = ""
            rows.append(row)
    return pd.DataFrame(rows)


def add_scaffold_stage_ratios(eval_summary: pd.DataFrame) -> pd.DataFrame:
    frame = eval_summary.copy()
    ratio_specs = [
        ("median_rg_ratio", "rg"),
        ("median_q90_distance_ratio", "q90"),
        ("long_range_cis_slope_median", "long_range_cis_slope"),
    ]
    for source_name, group in frame.groupby("source_config_name", sort=False):
        by_stage = {str(row["stage"]): row for row in group.to_dict("records")}
        init = by_stage.get("init_diploid_split")
        final = by_stage.get("final")
        after = by_stage.get("after_iter1")
        if init is None:
            continue
        idx = frame["source_config_name"].astype(str) == str(source_name)
        for metric, label in ratio_specs:
            init_value = safe_float(init.get(metric))
            final_value = safe_float(final.get(metric)) if final is not None else math.nan
            after_value = safe_float(after.get(metric)) if after is not None else math.nan
            frame.loc[idx, f"final_over_init_diploid_{label}_ratio"] = finite_ratio(final_value, init_value)
            frame.loc[idx, f"after_iter1_over_init_diploid_{label}_ratio"] = finite_ratio(after_value, init_value)
    return frame


def collect_scaffold_stage_table(eval_summary: pd.DataFrame, filename: str) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    meta = {
        str(row["config_name"]): {
            "source_config_name": row.get("source_config_name", ""),
            "stage": row.get("stage", ""),
        }
        for row in eval_summary.to_dict("records")
    }
    for row in eval_summary.to_dict("records"):
        path = Path(str(row.get("eval_dir", ""))) / filename
        if not path.exists():
            continue
        frame = read_tsv_trimmed(path)
        config_name = str(row["config_name"])
        frame.insert(0, "stage", meta.get(config_name, {}).get("stage", ""))
        frame.insert(0, "source_config_name", meta.get(config_name, {}).get("source_config_name", ""))
        frames.append(frame)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def write_scaffold_stage_interpretation(root: Path, eval_summary: pd.DataFrame, prefix: str = "scaffold_stage") -> None:
    lines = [
        "# P9016 Scaffold Stage Interpretation",
        "",
        "This is eval-only analysis. CHARM/3DG is used only after each diagnostic coordinate file exists. The diagnostic dump path is explicit and does not alter default training behavior.",
        "The `init_haploid_scaffold` file is evaluated as two duplicated copy rows per bead, so homolog separation at that stage is not biologically interpretable.",
        "copy0/copy1 labels remain gauge labels and are not maternal/paternal haplotypes.",
        "",
        "## Decision Rules",
        "",
        "- If `init_diploid_split` has reasonable Rg/slope but `final` is compressed, the EM/contact objective is the collapse source.",
        "- If `init_diploid_split` is already compressed, scaffold generation or global scale calibration is the first problem.",
        "- If `after_iter1` is already compressed while init is not, early E-step/contact targets are too aggressive.",
        "- If compression accumulates from init to after_iter1 to final, prioritize d_scale, repulsion balance, or anti-collapse regularization.",
        "",
        "## Per-Config Readout",
        "",
    ]
    for source_name, group in eval_summary.groupby("source_config_name", sort=False):
        by_stage = {str(row["stage"]): row for row in group.to_dict("records")}
        init = by_stage.get("init_diploid_split", {})
        after = by_stage.get("after_iter1", {})
        final = by_stage.get("final", {})
        init_rg = safe_float(init.get("median_rg_ratio"))
        init_slope = safe_float(init.get("median_cis_distance_slope"))
        after_rg = safe_float(after.get("median_rg_ratio"))
        after_slope = safe_float(after.get("median_cis_distance_slope"))
        final_rg = safe_float(final.get("median_rg_ratio"))
        final_slope = safe_float(final.get("median_cis_distance_slope"))
        init_reasonable = init_rg >= 0.7 and init_slope >= 0.7
        final_compressed = final_rg < 0.7 or final_slope < 0.7
        after_compressed = after_rg < 0.7 or after_slope < 0.7
        if not init_reasonable:
            origin = "init/scaffold-scale already compressed"
        elif after_compressed:
            origin = "collapse appears by after_iter1"
        elif final_compressed:
            origin = "collapse accumulates during EM/contact relaxation"
        else:
            origin = "no severe collapse by current thresholds"
        lines.extend([
            f"### {source_name}",
            "",
            f"- init_diploid_split: median Rg {fmt(init_rg)}, cis slope {fmt(init_slope)}, long-range cis Spearman {fmt(init.get('long_range_cis_pooled_spearman'))}.",
            f"- after_iter1: median Rg {fmt(after_rg)}, cis slope {fmt(after_slope)}, final/init Rg ratio proxy {fmt(after.get('after_iter1_over_init_diploid_rg_ratio'))}.",
            f"- final: median Rg {fmt(final_rg)}, cis slope {fmt(final_slope)}, final/init Rg ratio {fmt(final.get('final_over_init_diploid_rg_ratio'))}, final/init q90 ratio {fmt(final.get('final_over_init_diploid_q90_ratio'))}.",
            f"- trans/centroid absolute status at final: trans {pass_label(final.get('trans_metric_absolute_pass'))}; centroid {pass_label(final.get('centroid_metric_absolute_pass'))}.",
            f"- interpretation: {origin}.",
            "",
        ])
    (root / f"{prefix}_interpretation.md").write_text("\n".join(lines), encoding="utf-8")


def write_scaffold_stage_outputs(
    root: Path,
    summary: pd.DataFrame,
    tdg_path: Path,
    trans_sample_per_chrpair: int,
    baseline_trans_sample_per_chrpair: int,
    sample_seed: int,
    run_baselines: bool,
    prefix: str = "scaffold_stage",
) -> pd.DataFrame:
    stage_summary = scaffold_stage_summary_rows(summary)
    stage_summary.to_csv(output_path(root, prefix, "stage_input_summary.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    eval_summary = evaluate_candidates(
        stage_summary,
        tdg_path,
        root / "eval",
        trans_sample_per_chrpair,
        baseline_trans_sample_per_chrpair,
        sample_seed,
        run_baselines=run_baselines,
    )
    eval_summary = add_scaffold_stage_ratios(eval_summary)
    eval_summary.to_csv(output_path(root, prefix, "eval_summary.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    rg = collect_scaffold_stage_table(eval_summary, "eval.rg_by_chr.tsv")
    if not rg.empty:
        rg.to_csv(output_path(root, prefix, "rg_by_chr.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    quantiles = collect_scaffold_stage_table(eval_summary, "eval.distance_quantiles_by_chr.tsv")
    if not quantiles.empty:
        quantiles.to_csv(output_path(root, prefix, "distance_quantiles.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    write_scaffold_stage_interpretation(root, eval_summary, prefix)
    return eval_summary


def write_dscale_interpretation(
    root: Path,
    eval_summary: pd.DataFrame,
    force_normalized: pd.DataFrame,
    dscale_diag: pd.DataFrame,
    prefix: str = "dscale_ablation",
) -> None:
    force_by_config = {str(row["config_name"]): row for row in force_normalized.to_dict("records")}
    dscale_focus = dscale_diag.loc[
        (dscale_diag.get("stratification", pd.Series(dtype=str)).astype(str) == "state_p_bin")
        & (dscale_diag.get("available", pd.Series(dtype=int)).fillna(0).astype(int) == 1)
    ] if not dscale_diag.empty else pd.DataFrame()
    lines = [
        "# P9016 d_scale Ablation Interpretation",
        "",
        "This is post-hoc eval-only analysis. The training runner reads only the P9016 pairs file; CHARM/3DG is used after outputs exist.",
        "The ablation compares `raw_count` with `expected_count` under the same cis_inter_ratio prior, entropy rho_train, unphased scaffold initialization, and repulsion multiplier.",
        "",
        "## Decision Rules",
        "",
        "- If `expected_count` improves Rg/slope while preserving cis, state-specific d_scale geometry contributes to collapse.",
        "- If `expected_count` improves trans residual but not Rg, uncertain copy-state force geometry improved but collapse needs a separate anti-collapse term.",
        "- If `expected_count` does not help, prioritize scaffold-scale anti-collapse or aggregate trans constraints.",
        "- Reject any mode that improves trans Spearman only by destroying all-chr cis or causing homolog collapse.",
        "",
        "## Pairwise Raw vs Expected",
        "",
    ]
    rows = {str(row["config_name"]): row for row in eval_summary.to_dict("records")}
    for base in ["n100_rs100", "n120_rs80", "n80_rs100"]:
        raw = rows.get(f"dscale_{base}_raw", {})
        expected = rows.get(f"dscale_{base}_expected", {})
        if not raw or not expected:
            continue
        delta_rg = safe_float(expected.get("median_rg_ratio")) - safe_float(raw.get("median_rg_ratio"))
        delta_slope = safe_float(expected.get("median_cis_distance_slope")) - safe_float(raw.get("median_cis_distance_slope"))
        delta_cis = safe_float(expected.get("all_chr_cis_min_copy_pearson_mean")) - safe_float(raw.get("all_chr_cis_min_copy_pearson_mean"))
        delta_trans_resid = safe_float(expected.get("trans_contact_residual_median")) - safe_float(raw.get("trans_contact_residual_median"))
        if delta_rg > 0.05 and delta_slope > 0.05 and delta_cis > -0.05:
            readout = "expected_count improves scale while roughly preserving cis"
        elif delta_trans_resid < -0.05 and not (delta_rg > 0.05 or delta_slope > 0.05):
            readout = "expected_count improves contact residual more than global scale"
        elif delta_cis < -0.05:
            readout = "expected_count damages cis quality"
        else:
            readout = "no clear improvement from expected_count"
        lines.extend([
            f"### {base}",
            "",
            f"- raw: cis mean Pearson {fmt(raw.get('all_chr_cis_min_copy_pearson_mean'))}, Rg {fmt(raw.get('median_rg_ratio'))}, slope {fmt(raw.get('median_cis_distance_slope'))}, trans residual {fmt(raw.get('trans_contact_residual_median'))}, trans absolute {pass_label(raw.get('trans_metric_absolute_pass'))}.",
            f"- expected_count: cis mean Pearson {fmt(expected.get('all_chr_cis_min_copy_pearson_mean'))}, Rg {fmt(expected.get('median_rg_ratio'))}, slope {fmt(expected.get('median_cis_distance_slope'))}, trans residual {fmt(expected.get('trans_contact_residual_median'))}, trans absolute {pass_label(expected.get('trans_metric_absolute_pass'))}.",
            f"- deltas expected-raw: Rg {fmt(delta_rg)}, slope {fmt(delta_slope)}, cis Pearson {fmt(delta_cis)}, trans residual {fmt(delta_trans_resid)}.",
            f"- normalized force raw/expected trans per wedge: {fmt(force_by_config.get(f'dscale_{base}_raw', {}).get('force_l1_per_wedge_trans'))} / {fmt(force_by_config.get(f'dscale_{base}_expected', {}).get('force_l1_per_wedge_trans'))}.",
            f"- interpretation: {readout}.",
            "",
        ])
    if not dscale_focus.empty:
        lines.extend([
            "## Low-Posterior d_scale Check",
            "",
            "The d_scale diagnostic stratifies all four copy-state wedges per bpair. In expected_count mode, low posterior states should have larger state-specific d_scale and no NaN values.",
            "",
        ])
        for config_name, group in dscale_focus.groupby("config_name", sort=False):
            low = group.loc[group["stratum"].astype(str).isin(["state_p[0,1e-06)", "state_p[1e-06,0.0001)", "state_p[0.0001,0.001)", "state_p[0.001,0.01)"])]
            high = group.loc[group["stratum"].astype(str) == "state_p[0.25,1)"]
            lines.append(
                f"- {config_name}: low-state median d_scale {fmt(finite_median(low['median_state_d_scale']) if not low.empty else math.nan)}, high-state median d_scale {fmt(finite_median(high['median_state_d_scale']) if not high.empty else math.nan)}."
            )
    (root / f"{prefix}_interpretation.md").write_text("\n".join(lines), encoding="utf-8")


def run_postprocess(args: argparse.Namespace) -> None:
    root = args.root
    prefix = args.prefix or args.run_set
    order = expected_config_order(args.run_set)
    summary = normalize_candidate_summary(root, require_exact_candidates=True, expected_order=order)
    schema = validate_candidate_file_schemas(summary)
    summary.to_csv(output_path(root, prefix, "summary.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    schema.to_csv(output_path(root, prefix, "schema_check.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    if args.run_set == "scaffold_stage":
        write_scaffold_stage_outputs(
            root,
            summary,
            args.tdg,
            args.trans_sample_per_chrpair,
            args.baseline_trans_sample_per_chrpair,
            args.sample_seed,
            run_baselines=not args.no_baselines,
            prefix=prefix,
        )
        return
    eval_summary = evaluate_candidates(
        summary,
        args.tdg,
        root / "eval",
        args.trans_sample_per_chrpair,
        args.baseline_trans_sample_per_chrpair,
        args.sample_seed,
        run_baselines=not args.no_baselines,
    )
    eval_summary.to_csv(output_path(root, prefix, "eval_summary.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    write_posterior_outputs(root, eval_summary, prefix)
    force_summary = write_force_output(root, summary, prefix)
    force_normalized = write_force_normalized_output(root, summary, eval_summary, prefix)
    force_cancellation = write_force_cancellation_output(root, summary, prefix)
    residual_stratified = write_residual_stratified_output(root, summary, prefix)
    write_contact_residual_output(root, eval_summary, prefix)
    dscale_diag = pd.DataFrame()
    if args.run_set == "dscale_ablation":
        dscale_diag = write_dscale_diag_output(root, summary, prefix)
        write_dscale_interpretation(root, eval_summary, force_normalized, dscale_diag, prefix)
    else:
        write_interpretation(root, eval_summary, force_summary, prefix)
    write_mechanism_update(root, eval_summary, force_normalized, force_cancellation, residual_stratified, prefix)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--prefix", default=None)
    parser.add_argument("--run-set", choices=sorted(RUN_SET_CONFIG_ORDERS), default="candidate_rerun")
    parser.add_argument("--tdg", type=Path, default=DEFAULT_TDG)
    parser.add_argument("--trans-sample-per-chrpair", type=int, default=2048)
    parser.add_argument("--baseline-trans-sample-per-chrpair", type=int, default=256)
    parser.add_argument("--sample-seed", type=int, default=17)
    parser.add_argument("--no-baselines", action="store_true")
    args = parser.parse_args()
    run_postprocess(args)
    print(f"wrote {args.prefix or args.run_set} diagnostics under {args.root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
