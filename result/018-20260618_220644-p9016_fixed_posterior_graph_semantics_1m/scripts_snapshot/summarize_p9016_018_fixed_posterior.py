#!/usr/bin/env python3
"""Summarize experiment 018 fixed-posterior graph-semantics diagnostics."""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path


COLUMNS = [
    "config_name",
    "graph_mode",
    "dscale_mode",
    "dscale_gamma",
    "trans_dscale_multiplier",
    "repulsion_multiplier",
    "native_softall_replay",
    "min_sep_unit",
    "lambda_sep",
    "n_fixed_edges",
    "final_contact_energy",
    "final_repulsion_energy",
    "final_backbone_energy",
    "final_sep_energy",
    "final_min_sep",
    "sep_p05",
    "final_mean_sep",
    "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans",
    "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "whole_chrom_top1_accuracy",
    "chr_pair_oracle_top1_accuracy",
    "nearest_geometry_state_accuracy",
    "same_cross_accuracy",
    "truth_state_is_nearest_fraction",
    "posterior_top1_is_nearest_fraction",
    "posterior_top1_distance_rank_mean",
    "trans_force_per_k",
    "trans_mean_r",
    "trans_frac_k_tail_attractive",
    "trans_frac_k_zero_force",
    "trans_cancellation_ratio",
    "cis_force_per_k",
    "cis_mean_r",
    "status",
]


def read_kv(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    out: dict[str, str] = {}
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header != ["key", "value"]:
            return out
        for row in reader:
            if len(row) >= 2:
                out[row[0]] = row[1]
    return out


def read_single_row_tsv(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return dict(row)
    return {}


def read_force_summary(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    if not path.is_file():
        return rows
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cls = row.get("contact_class", "")
            if cls:
                rows[cls] = dict(row)
    return rows


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.9g}"
    s = str(value)
    return s if s != "" else "NA"


def pick(*values: str | None) -> str:
    for value in values:
        if value is not None and value != "":
            return value
    return "NA"


def summarize(root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    source_dir_env = os.environ.get("HK_FIXED_POSTERIOR_SOURCE_DIR", "")
    source_eval_env = os.environ.get("HK_FIXED_POSTERIOR_SOURCE_EVAL_DIR", "")
    if source_dir_env:
        source_dir = Path(source_dir_env)
        source_eval_dir = Path(source_eval_env) if source_eval_env else Path()
        manifest = read_kv(source_dir / "p9016_full.manifest.tsv")
        eval_summary = read_kv(source_eval_dir / "summary.tsv") if source_eval_env else {}
        sep = read_single_row_tsv(source_dir / "p9016_full.sep_diag.tsv")
        force = read_force_summary(source_dir / "p9016_full.force_class_diag.tsv")
        rows.append({
            "config_name": "source_016_before_018_relax",
            "graph_mode": "source_native_016",
            "dscale_mode": pick(manifest.get("dscale_mode"), manifest.get("d_scale_mode")),
            "dscale_gamma": pick(manifest.get("dscale_posterior_gamma"), manifest.get("d_scale_posterior_gamma")),
            "trans_dscale_multiplier": "1",
            "repulsion_multiplier": pick(manifest.get("repulsion_multiplier")),
            "native_softall_replay": "1",
            "min_sep_unit": pick(manifest.get("min_sep_unit")),
            "lambda_sep": pick(manifest.get("lambda_sep")),
            "n_fixed_edges": "NA",
            "final_contact_energy": pick(manifest.get("final_contact_energy")),
            "final_repulsion_energy": pick(manifest.get("final_repulsion_energy")),
            "final_backbone_energy": pick(manifest.get("final_backbone_energy")),
            "final_sep_energy": pick(manifest.get("final_sep_energy")),
            "final_min_sep": pick(manifest.get("final_min_sep"), sep.get("sep_min")),
            "sep_p05": pick(sep.get("sep_p05")),
            "final_mean_sep": pick(manifest.get("final_mean_sep"), sep.get("sep_mean")),
            "model_top1_accuracy_genome_all": pick(eval_summary.get("model_top1_accuracy_genome_all")),
            "model_top1_accuracy_genome_cis": pick(eval_summary.get("model_top1_accuracy_genome_cis")),
            "model_top1_accuracy_genome_trans": pick(eval_summary.get("model_top1_accuracy_genome_trans")),
            "model_same_cross_accuracy_genome_trans": pick(eval_summary.get("model_same_cross_accuracy_genome_trans")),
            "mean_per_chrom_cis_distance_spearman": pick(eval_summary.get("mean_per_chrom_cis_distance_spearman")),
            "whole_chrom_top1_accuracy": "NA",
            "chr_pair_oracle_top1_accuracy": "NA",
            "nearest_geometry_state_accuracy": "NA",
            "same_cross_accuracy": "NA",
            "truth_state_is_nearest_fraction": "NA",
            "posterior_top1_is_nearest_fraction": "NA",
            "posterior_top1_distance_rank_mean": "NA",
            "trans_force_per_k": "NA",
            "trans_mean_r": "NA",
            "trans_frac_k_tail_attractive": "NA",
            "trans_frac_k_zero_force": "NA",
            "trans_cancellation_ratio": "NA",
            "cis_force_per_k": "NA",
            "cis_mean_r": "NA",
            "status": pick(manifest.get("status"), "OK"),
        })
    outputs_root = root / "outputs"
    for manifest_path in sorted(outputs_root.glob("*/p9016_full.manifest.tsv")):
        config_dir = manifest_path.parent
        config = config_dir.name
        eval_dir = root / "eval" / config
        diag_dir = root / "diagnostics" / config
        manifest = read_kv(manifest_path)
        eval_summary = read_kv(eval_dir / "summary.tsv")
        sep = read_single_row_tsv(config_dir / "p9016_full.sep_diag.tsv")
        geom = read_single_row_tsv(diag_dir / "trans_gauge_geometry_summary.tsv")
        force = read_force_summary(diag_dir / "force_regime_summary.tsv")
        row = {
            "config_name": config,
            "graph_mode": pick(manifest.get("fixed_posterior_graph_mode"), manifest.get("mstep_graph_mode")),
            "dscale_mode": pick(manifest.get("dscale_mode"), manifest.get("d_scale_mode")),
            "dscale_gamma": pick(manifest.get("dscale_posterior_gamma"), manifest.get("d_scale_posterior_gamma")),
            "trans_dscale_multiplier": pick(manifest.get("trans_dscale_multiplier")),
            "repulsion_multiplier": pick(manifest.get("repulsion_multiplier")),
            "native_softall_replay": pick(manifest.get("native_softall_replay")),
            "min_sep_unit": pick(manifest.get("min_sep_unit")),
            "lambda_sep": pick(manifest.get("lambda_sep")),
            "n_fixed_edges": pick(manifest.get("n_fixed_edges")),
            "final_contact_energy": pick(manifest.get("final_contact_energy")),
            "final_repulsion_energy": pick(manifest.get("final_repulsion_energy")),
            "final_backbone_energy": pick(manifest.get("final_backbone_energy")),
            "final_sep_energy": pick(manifest.get("final_sep_energy")),
            "final_min_sep": pick(manifest.get("final_min_sep"), sep.get("sep_min")),
            "sep_p05": pick(sep.get("sep_p05")),
            "final_mean_sep": pick(manifest.get("final_mean_sep"), sep.get("sep_mean")),
            "model_top1_accuracy_genome_all": pick(eval_summary.get("model_top1_accuracy_genome_all")),
            "model_top1_accuracy_genome_cis": pick(eval_summary.get("model_top1_accuracy_genome_cis")),
            "model_top1_accuracy_genome_trans": pick(eval_summary.get("model_top1_accuracy_genome_trans")),
            "model_same_cross_accuracy_genome_trans": pick(eval_summary.get("model_same_cross_accuracy_genome_trans")),
            "mean_per_chrom_cis_distance_spearman": pick(eval_summary.get("mean_per_chrom_cis_distance_spearman")),
            "whole_chrom_top1_accuracy": pick(geom.get("whole_chrom_top1_accuracy")),
            "chr_pair_oracle_top1_accuracy": pick(geom.get("chr_pair_oracle_top1_accuracy")),
            "nearest_geometry_state_accuracy": pick(geom.get("nearest_geometry_state_accuracy")),
            "same_cross_accuracy": pick(geom.get("same_cross_accuracy")),
            "truth_state_is_nearest_fraction": pick(geom.get("truth_state_is_nearest_fraction")),
            "posterior_top1_is_nearest_fraction": pick(geom.get("posterior_top1_is_nearest_fraction")),
            "posterior_top1_distance_rank_mean": pick(geom.get("posterior_top1_distance_rank_mean")),
            "trans_force_per_k": pick(force.get("trans", {}).get("force_per_k")),
            "trans_mean_r": pick(force.get("trans", {}).get("mean_r")),
            "trans_frac_k_tail_attractive": pick(force.get("trans", {}).get("frac_k_tail_attractive")),
            "trans_frac_k_zero_force": pick(force.get("trans", {}).get("frac_k_zero_force")),
            "trans_cancellation_ratio": pick(force.get("trans", {}).get("cancellation_ratio")),
            "cis_force_per_k": pick(force.get("cis", {}).get("force_per_k")),
            "cis_mean_r": pick(force.get("cis", {}).get("mean_r")),
            "status": pick(manifest.get("status"), "OK"),
        }
        rows.append(row)
    return rows


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: scripts/summarize_p9016_018_fixed_posterior.py <run_root>", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=COLUMNS, lineterminator="\n")
    writer.writeheader()
    for row in summarize(root):
        writer.writerow({col: fmt(row.get(col)) for col in COLUMNS})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
