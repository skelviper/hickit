#!/usr/bin/env python3
"""Summarize P9016 4Mb trans-compression diagnostics from finished eval TSVs.

This script is intentionally read-only with respect to the experiment root.  It
does not run training, reload coordinates, or rerun CHARM/3DG evaluation.  It
only gathers already-written TSV diagnostics into small tables that make the
absolute trans failure easier to inspect.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence


ABS_RATIO_MIN = 0.70
TRANS_CIS_MIN_PEARSON = 0.30
CENTROID_RELATIVE_MIN = 0.50


RUN_HEALTH_COLUMNS = [
    "config_name",
    "output_dir",
    "scan_config",
    "input_path",
    "sample",
    "runner_family",
    "default_profile",
    "resolution_label",
    "resolution",
    "bin_size_bp",
    "n_iter",
    "relax_steps",
    "base_k_mode",
    "init_mode",
    "prior_mode",
    "rho_train_mode",
    "state_weight_mode",
    "uses_phase_labels",
    "uses_charm_or_reference",
    "uses_charm_for_training",
    "audit_status",
    "n_bad_iter",
    "n_relax_nonfinite_iter",
    "n_coord_nonfinite",
    "n_repulsion_nonfinite_step",
    "posterior_refreshed_after_final_relax",
    "n_raw",
    "n_bpair",
    "n_beads",
    "n_raw_cis",
    "n_raw_trans",
    "n_bpair_cis",
    "n_bpair_trans",
    "final_mean_entropy",
    "final_mean_pU",
    "final_sum_wedge_k",
    "final_mean_rho_train_bpair",
    "final_contact_energy",
    "final_repulsion_energy",
    "final_backbone_energy",
]


ABSOLUTE_GATE_COLUMNS = [
    "config_name",
    "state_weight_mode",
    "trans_metric_rank",
    "trans_metric_rank_best",
    "trans_relative_metric",
    "trans_metric_absolute_pass",
    "trans_metric_absolute_fail_reasons",
    "centroid_metric_rank",
    "centroid_metric_rank_best",
    "centroid_relative_metric",
    "centroid_metric_absolute_pass",
    "centroid_metric_absolute_fail_reasons",
    "chr1_cis_min_copy_pearson",
    "all_chr_cis_min_copy_pearson_mean",
    "all_chr_cis_pooled_spearman",
    "median_rg_ratio",
    "median_cis_distance_slope",
]


TRANS_COMPRESSION_COLUMNS = [
    "config_name",
    "state_weight_mode",
    "trans_distance_spearman",
    "trans_distance_slope",
    "trans_slope_deficit_to_0_70",
    "trans_median_distance_ratio",
    "trans_ratio_deficit_to_0_70",
    "trans_abs_log_ratio_median",
    "trans_fraction_ratio_lt_0_7",
    "trans_fraction_ratio_gt_1_3",
    "trans_distance_intercept",
    "trans_residual_median_abs",
    "trans_residual_mad",
    "trans_centroid_residual_spearman",
    "centroid_distance_spearman",
    "centroid_distance_slope",
    "centroid_slope_deficit_to_0_70",
    "centroid_median_distance_ratio",
    "centroid_ratio_deficit_to_0_70",
    "median_rg_ratio",
    "rg_ratio_deficit_to_0_70",
    "median_cis_distance_slope",
    "cis_slope_deficit_to_0_70",
    "fraction_chr_rg_ratio_lt_0_7",
    "fraction_chr_slope_lt_0_7",
    "procrustes_global_over_per_chr_ratio",
    "whole_genome_distance_spearman",
    "whole_genome_median_distance_ratio",
]


POSTERIOR_FORCE_COLUMNS = [
    "config_name",
    "state_weight_mode",
    "mean_pU_cis",
    "mean_pU_trans",
    "median_pU_cis",
    "median_pU_trans",
    "mean_pmax_cis",
    "mean_pmax_trans",
    "mean_margin_cis",
    "mean_margin_trans",
    "mean_rho_train_cis",
    "mean_rho_train_trans",
    "fraction_trans_contacts_with_rho_train_near_zero",
    "sum_effective_k_cis",
    "sum_effective_k_trans",
    "effective_k_cis_trans_ratio",
    "trans_fraction_effective_k",
    "cis_contact_count_vs_inv_distance_spearman",
    "trans_contact_count_vs_inv_distance_spearman",
    "final_contact_energy_cis",
    "final_contact_energy_trans",
    "final_contact_force_l1_cis",
    "final_contact_force_l1_trans",
    "force_l1_cis_trans_ratio",
    "repulsion_to_trans_force_ratio",
    "final_backbone_force_l1",
    "final_repulsion_force_l1",
    "final_contact_wedges_cis",
    "final_contact_wedges_trans",
]


RANK_DIGEST_COLUMNS = [
    "rank_table",
    "rank",
    "config_name",
    "rank_metric_name",
    "rank_metric_kind",
    "rank_metric_value",
    "trans_metric_absolute_pass",
    "centroid_metric_absolute_pass",
    "trans_distance_spearman",
    "trans_distance_slope",
    "trans_median_distance_ratio",
    "centroid_distance_spearman",
    "centroid_distance_slope",
    "centroid_median_distance_ratio",
    "median_rg_ratio",
    "all_chr_cis_min_copy_pearson_mean",
    "note",
]


CHRPAIR_COLUMNS = [
    "config_name",
    "source_rank",
    "chr_a",
    "chr_b",
    "copy_a",
    "copy_b",
    "n_pairs",
    "model_centroid_distance",
    "ref_centroid_distance",
    "spearman",
    "slope",
    "median_distance_ratio",
    "ratio_deficit_to_0_70",
    "abs_log_ratio_median",
    "residual_median_abs",
    "residual_mad",
    "compression_reason",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize finished P9016 4Mb trans-compression eval TSVs."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("/tmp/hk_blind_p9016_mr4mb_inter_gate_1778044113"),
        help="Finished experiment root containing matrix_summary.tsv and eval outputs.",
    )
    parser.add_argument(
        "--eval-name",
        default="eval_rerun",
        help="Evaluation subdirectory name under --root.",
    )
    parser.add_argument(
        "--matrix-summary",
        type=Path,
        default=None,
        help=(
            "Optional matrix summary TSV. If omitted, use matrix_summary.tsv "
            "under --root, or auto-detect a single matrix_summary*.tsv file."
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("diagnostics/p9016_mr4mb_inter_gate_1778044113_trans_compression"),
        help="Directory where summary TSVs and README are written.",
    )
    parser.add_argument(
        "--top-configs",
        type=int,
        default=3,
        help="Number of trans-ranked configs to include in chrpair outlier summaries.",
    )
    parser.add_argument(
        "--top-chrpairs",
        type=int,
        default=20,
        help="Number of chrpair/copy outliers to keep per selected config.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def find_matrix_summary(root: Path, requested: Path | None) -> Path:
    if requested is not None:
        return requested.resolve()

    default = root / "matrix_summary.tsv"
    if default.exists():
        return default

    candidates = sorted(root.glob("matrix_summary*.tsv"))
    if len(candidates) == 1:
        return candidates[0]

    if not candidates:
        raise FileNotFoundError(default)

    listed = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(
        f"multiple matrix summaries found under {root}; pass --matrix-summary explicitly: {listed}"
    )


def read_manifest(path: Path) -> Dict[str, str]:
    if not path.exists():
        return {}
    values: Dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "key" not in reader.fieldnames or "value" not in reader.fieldnames:
            return values
        for row in reader:
            key = row.get("key", "")
            if key:
                values[key] = row.get("value", "")
    return values


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: clean_value(row.get(col, "")) for col in columns})


def clean_value(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return value


def safe_float(value: object) -> float:
    try:
        if value is None or value == "":
            return math.nan
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def safe_int(value: object, default: int = 0) -> int:
    try:
        if value is None or value == "":
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def finite_deficit(value: object, threshold: float = ABS_RATIO_MIN) -> float:
    numeric = safe_float(value)
    if not math.isfinite(numeric):
        return math.nan
    return max(0.0, threshold - numeric)


def select_columns(rows: Iterable[Mapping[str, object]], columns: Sequence[str]) -> List[Dict[str, object]]:
    return [{col: row.get(col, "") for col in columns} for row in rows]


def merge_manifest_fields(rows: Iterable[Mapping[str, object]]) -> List[Dict[str, object]]:
    manifest_keys = [
        "input_path",
        "sample",
        "runner_family",
        "default_profile",
        "resolution_label",
        "resolution",
        "bin_size_bp",
        "base_k_mode",
        "init_mode",
        "prior_mode",
        "rho_train_mode",
        "state_weight_mode",
        "uses_phase_labels",
        "uses_charm_or_reference",
        "uses_charm_for_training",
    ]
    out: List[Dict[str, object]] = []
    for row in rows:
        item = dict(row)
        output_dir = str(item.get("output_dir", ""))
        manifest = read_manifest(Path(output_dir) / "p9016_full.manifest.tsv") if output_dir else {}
        for key in manifest_keys:
            if not item.get(key, ""):
                item[key] = manifest.get(key, "")
        out.append(item)
    return out


def sort_by_float(rows: List[Dict[str, object]], key: str, reverse: bool = True) -> List[Dict[str, object]]:
    return sorted(
        rows,
        key=lambda row: safe_float(row.get(key, "")) if math.isfinite(safe_float(row.get(key, ""))) else -math.inf,
        reverse=reverse,
    )


def add_compression_fields(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    for row in rows:
        item = dict(row)
        item["trans_slope_deficit_to_0_70"] = finite_deficit(item.get("trans_distance_slope"))
        item["trans_ratio_deficit_to_0_70"] = finite_deficit(item.get("trans_median_distance_ratio"))
        item["centroid_slope_deficit_to_0_70"] = finite_deficit(item.get("centroid_distance_slope"))
        item["centroid_ratio_deficit_to_0_70"] = finite_deficit(item.get("centroid_median_distance_ratio"))
        item["rg_ratio_deficit_to_0_70"] = finite_deficit(item.get("median_rg_ratio"))
        item["cis_slope_deficit_to_0_70"] = finite_deficit(item.get("median_cis_distance_slope"))
        out.append(item)
    return out


def rank_metric_value(row: Mapping[str, object]) -> object:
    metric = str(row.get("rank_metric_name", ""))
    if metric == "trans_relative_metric":
        return row.get("trans_relative_metric", "")
    if metric == "centroid_relative_metric":
        return row.get("centroid_relative_metric", "")
    if metric == "median_rg_ratio":
        return row.get("median_rg_ratio", "")
    if metric == "all_chr_cis_min_copy_pearson_mean":
        return row.get("all_chr_cis_min_copy_pearson_mean", "")
    if metric == "pareto_score":
        return row.get("pareto_score", "")
    return row.get(metric, "")


def build_rank_digest(eval_dir: Path) -> List[Dict[str, object]]:
    rank_specs = [
        ("rank_by_trans_relative_metric.tsv", "trans_relative"),
        ("rank_by_centroid_relative_metric.tsv", "centroid_relative"),
        ("rank_by_compaction.tsv", "compaction"),
        ("rank_by_all_chr_cis.tsv", "all_chr_cis"),
        ("pareto_summary.tsv", "pareto"),
    ]
    digest: List[Dict[str, object]] = []
    for file_name, table_name in rank_specs:
        path = eval_dir / file_name
        if not path.exists():
            continue
        for rank, row in enumerate(read_tsv(path), start=1):
            item = dict(row)
            item["rank_table"] = table_name
            item["rank"] = rank
            item["rank_metric_value"] = rank_metric_value(item)
            trans_pass = safe_int(item.get("trans_metric_absolute_pass"))
            centroid_pass = safe_int(item.get("centroid_metric_absolute_pass"))
            item["note"] = "absolute_pass" if trans_pass and centroid_pass else "relative_only"
            digest.append(item)
    return digest


def compression_reason(row: Mapping[str, object]) -> str:
    reasons: List[str] = []
    slope = safe_float(row.get("slope"))
    ratio = safe_float(row.get("median_distance_ratio"))
    spearman = safe_float(row.get("spearman"))
    if math.isfinite(slope) and slope < ABS_RATIO_MIN:
        reasons.append("low_slope")
    if math.isfinite(ratio) and ratio < ABS_RATIO_MIN:
        reasons.append("low_median_ratio")
    if math.isfinite(spearman) and spearman < 0.10:
        reasons.append("weak_local_rank")
    if not reasons:
        reasons.append("high_abs_log_ratio")
    return ",".join(reasons)


def chrpair_outlier_score(row: Mapping[str, object]) -> float:
    ratio_deficit = finite_deficit(row.get("median_distance_ratio"))
    abs_log = safe_float(row.get("abs_log_ratio_median"))
    slope_deficit = finite_deficit(row.get("slope"))
    weak_rank = max(0.0, 0.10 - safe_float(row.get("spearman"))) if math.isfinite(safe_float(row.get("spearman"))) else 0.0
    score = 0.0
    for value in [ratio_deficit, abs_log, slope_deficit, weak_rank]:
        if math.isfinite(value):
            score += value
    return score


def build_chrpair_outliers(
    eval_dir: Path,
    ranked_rows: Sequence[Mapping[str, object]],
    top_configs: int,
    top_chrpairs: int,
) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    selected = ranked_rows[: max(0, top_configs)]
    for rank, row in enumerate(selected, start=1):
        config = str(row.get("config_name", ""))
        if not config:
            continue
        path = eval_dir / config / "eval.trans_by_chrpair.tsv"
        if not path.exists():
            continue
        chr_rows = read_tsv(path)
        scored = sorted(chr_rows, key=chrpair_outlier_score, reverse=True)
        for chr_row in scored[: max(0, top_chrpairs)]:
            item = dict(chr_row)
            item["source_rank"] = rank
            item["ratio_deficit_to_0_70"] = finite_deficit(item.get("median_distance_ratio"))
            item["compression_reason"] = compression_reason(item)
            out.append(item)
    return out


def write_warning_digest(eval_dir: Path, out_dir: Path) -> None:
    path = eval_dir / "eval_acceptability_warnings.tsv"
    if not path.exists():
        return
    rows = read_tsv(path)
    columns = ["warning_type", "config_name", "message", "fail_reasons"]
    write_tsv(out_dir / "acceptability_warnings.tsv", rows, columns)


def summarize_counts(summary: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    total = len(summary)
    trans_pass = sum(safe_int(row.get("trans_metric_absolute_pass")) for row in summary)
    centroid_pass = sum(safe_int(row.get("centroid_metric_absolute_pass")) for row in summary)
    top = sort_by_float([dict(row) for row in summary], "trans_relative_metric")[:1]
    top_row = top[0] if top else {}
    return {
        "n_configs": total,
        "n_trans_absolute_pass": trans_pass,
        "n_centroid_absolute_pass": centroid_pass,
        "top_trans_config": top_row.get("config_name", ""),
        "top_trans_spearman": top_row.get("trans_distance_spearman", ""),
        "top_trans_slope": top_row.get("trans_distance_slope", ""),
        "top_trans_median_distance_ratio": top_row.get("trans_median_distance_ratio", ""),
        "top_trans_fail_reasons": top_row.get("trans_metric_absolute_fail_reasons", ""),
    }


def write_readme(
    out_dir: Path,
    root: Path,
    eval_dir: Path,
    matrix_summary_path: Path,
    counts: Mapping[str, object],
) -> None:
    text = f"""# P9016 Trans Compression Diagnostic

Generated: {datetime.now(timezone.utc).isoformat()}

Experiment root: `{root}`

Eval directory: `{eval_dir}`

Matrix summary: `{matrix_summary_path}`

This report is read-only with respect to the experiment root.  It summarizes
finished TSV diagnostics and does not rerun training, coordinate loading, or
CHARM/3DG evaluation.

## Summary

- Number of configs: {counts.get("n_configs", "")}
- Trans absolute-pass configs: {counts.get("n_trans_absolute_pass", "")}
- Centroid absolute-pass configs: {counts.get("n_centroid_absolute_pass", "")}
- Top trans-relative config: `{counts.get("top_trans_config", "")}`
- Top trans Spearman: {counts.get("top_trans_spearman", "")}
- Top trans slope: {counts.get("top_trans_slope", "")}
- Top trans median distance ratio: {counts.get("top_trans_median_distance_ratio", "")}
- Top trans fail reasons: `{counts.get("top_trans_fail_reasons", "")}`

## Tables

- `run_health.tsv`: run/audit/numeric-health columns from the matrix summary TSV.
- `absolute_gate_summary.tsv`: relative rank versus absolute gate status.
- `trans_compression.tsv`: trans, centroid, Rg, and cis compression indicators.
- `posterior_force_balance.tsv`: posterior ambiguity, rho, effective k, and force balance.
- `rank_warning_digest.tsv`: rank-table digest marking relative-only winners.
- `chrpair_failure_top.tsv`: top chromosome-pair/copy outliers for the top trans-ranked configs.
- `acceptability_warnings.tsv`: copied acceptability warnings from the eval directory.

Interpretation rule: a high relative metric is not an absolute success unless
`trans_metric_absolute_pass` and the relevant centroid/cis/scale gates also pass.
"""
    (out_dir / "README.md").write_text(text)


def write_provenance(
    out_dir: Path,
    root: Path,
    eval_dir: Path,
    matrix_summary_path: Path,
    counts: Mapping[str, object],
) -> None:
    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(root),
        "eval_dir": str(eval_dir),
        "matrix_summary_path": str(matrix_summary_path),
        "thresholds": {
            "ABS_RATIO_MIN": ABS_RATIO_MIN,
            "TRANS_CIS_MIN_PEARSON": TRANS_CIS_MIN_PEARSON,
            "CENTROID_RELATIVE_MIN": CENTROID_RELATIVE_MIN,
        },
        "summary": dict(counts),
    }
    (out_dir / "provenance.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    eval_dir = root / args.eval_name
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    matrix_summary_path = find_matrix_summary(root, args.matrix_summary)
    matrix_summary = read_tsv(matrix_summary_path)
    eval_summary = read_tsv(eval_dir / "eval_summary.tsv")
    eval_ranked = sort_by_float(eval_summary, "trans_relative_metric")

    run_health = select_columns(merge_manifest_fields(matrix_summary), RUN_HEALTH_COLUMNS)
    absolute_gate = select_columns(eval_ranked, ABSOLUTE_GATE_COLUMNS)
    trans_compression = select_columns(add_compression_fields(eval_ranked), TRANS_COMPRESSION_COLUMNS)
    posterior_force = select_columns(eval_ranked, POSTERIOR_FORCE_COLUMNS)
    rank_digest = build_rank_digest(eval_dir)
    chrpair_outliers = build_chrpair_outliers(eval_dir, eval_ranked, args.top_configs, args.top_chrpairs)

    write_tsv(out_dir / "run_health.tsv", run_health, RUN_HEALTH_COLUMNS)
    write_tsv(out_dir / "absolute_gate_summary.tsv", absolute_gate, ABSOLUTE_GATE_COLUMNS)
    write_tsv(out_dir / "trans_compression.tsv", trans_compression, TRANS_COMPRESSION_COLUMNS)
    write_tsv(out_dir / "posterior_force_balance.tsv", posterior_force, POSTERIOR_FORCE_COLUMNS)
    write_tsv(out_dir / "rank_warning_digest.tsv", rank_digest, RANK_DIGEST_COLUMNS)
    write_tsv(out_dir / "chrpair_failure_top.tsv", chrpair_outliers, CHRPAIR_COLUMNS)
    write_warning_digest(eval_dir, out_dir)

    counts = summarize_counts(eval_summary)
    write_readme(out_dir, root, eval_dir, matrix_summary_path, counts)
    write_provenance(out_dir, root, eval_dir, matrix_summary_path, counts)

    print(f"wrote P9016 trans-compression diagnostics to {out_dir}")
    print(
        "summary: "
        f"configs={counts['n_configs']} "
        f"trans_absolute_pass={counts['n_trans_absolute_pass']} "
        f"centroid_absolute_pass={counts['n_centroid_absolute_pass']} "
        f"top_trans={counts['top_trans_config']} "
        f"slope={counts['top_trans_slope']} "
        f"ratio={counts['top_trans_median_distance_ratio']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
