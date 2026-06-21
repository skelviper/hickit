#!/usr/bin/env python3
"""Summarize experiment 016 CHARM-contact model-condition sweep."""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path


COLUMNS = [
    "radius_slug",
    "radius_pr",
    "condition_name",
    "config_name",
    "status",
    "backend",
    "output_dir",
    "eval_dir",
    "input_pairs",
    "output_pairs_total",
    "output_pairs_cis",
    "output_pairs_trans",
    "candidate_pairs_total",
    "init_mode",
    "init_seed",
    "init_eps",
    "init_noise_scale",
    "init_scale",
    "init_eps_effective",
    "init_noise_scale_effective",
    "init_scale_effective",
    "d_scale_mode_input_string",
    "d_scale_mode",
    "dscale_mode",
    "d_scale_posterior_gamma",
    "dscale_posterior_gamma",
    "d_scale_eps_count",
    "dscale_effective_count_formula",
    "dscale_probability_weighted",
    "edge_k_probability_weighted",
    "legacy_expected_count_alias_used",
    "min_sep_unit",
    "lambda_sep",
    "prior_mode",
    "rho_train_mode",
    "temperature_start",
    "temperature_end",
    "uses_phase_labels",
    "uses_charm_or_reference",
    "uses_charm_for_training",
    "reference_derived_positive_control",
    "n_raw",
    "n_bpair",
    "n_raw_cis",
    "n_raw_trans",
    "n_raw_same_bin_excluded",
    "final_mean_entropy",
    "final_mean_pU",
    "final_min_sep",
    "final_mean_sep",
    "final_max_sep",
    "final_contact_energy",
    "final_repulsion_energy",
    "final_backbone_energy",
    "final_sep_energy",
    "final_sep_force_l1",
    "sep_p01",
    "sep_p05",
    "sep_p10",
    "sep_median",
    "sep_p90",
    "frac_sep_below_min",
    "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans",
    "model_same_cross_accuracy_genome_all",
    "model_same_cross_accuracy_genome_cis",
    "model_same_cross_accuracy_genome_trans",
    "model_pmax90_accuracy_genome_all",
    "model_pmax90_accuracy_genome_cis",
    "model_pmax90_accuracy_genome_trans",
    "model_pmax90_recall_genome_all",
    "model_pmax90_recall_genome_cis",
    "model_pmax90_recall_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "charm3dg_top1_accuracy_genome_trans",
    "charm3dg_same_cross_accuracy_genome_trans",
    "copytrack_sep_p05",
    "copytrack_sep_median",
    "copytrack_vector_cos_p05",
    "copytrack_vector_cos_median",
    "copytrack_frac_cos_lt_0",
    "copytrack_frac_projection_sign_switch",
]


def read_kv(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header != ["key", "value"]:
            return {}
        return {row[0]: row[1] for row in reader if len(row) >= 2}


def read_one(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        row = next(reader, None)
        return dict(row) if row else {}


def read_copytrack_all(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("chrom") == "ALL":
                return dict(row)
    return {}


def read_condition_matrix(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as fh:
        return {row["config_name"]: dict(row) for row in csv.DictReader(fh, delimiter="\t")}


def value(*sources: dict[str, str], key: str, default: str = "NA") -> str:
    for source in sources:
        v = source.get(key)
        if v not in (None, ""):
            return v
    return default


def discover_config_dirs(root: Path) -> list[Path]:
    outputs = root / "outputs"
    if not outputs.exists():
        return []
    dirs: set[Path] = set()
    for dirpath, _, filenames in os.walk(outputs, followlinks=True):
        if "p9016_full.manifest.tsv" in filenames:
            dirs.add(Path(dirpath))
    return sorted(dirs, key=lambda p: (p.parent.name, p.name))


def eval_dir_for(root: Path, radius_slug: str, config_name: str) -> Path:
    nested = root / "eval" / radius_slug / config_name
    if nested.exists():
        return nested
    flat = root / "eval" / config_name
    if flat.exists():
        return flat
    return nested


def status(manifest: dict[str, str], eval_dir: Path, eval_summary: dict[str, str]) -> str:
    eval_status = eval_dir / "eval_status.tsv"
    if eval_status.exists():
        with eval_status.open() as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 2 and parts[0] == "status":
                    return parts[1]
    train_status = value(manifest, key="status")
    if train_status not in ("NA", "OK"):
        return train_status
    if eval_summary:
        return "OK"
    return train_status


def summarize(root: Path, config_dir: Path, conditions: dict[str, dict[str, str]]) -> dict[str, str]:
    manifest = read_kv(config_dir / "p9016_full.manifest.tsv")
    sep = read_one(config_dir / "p9016_full.sep_diag.tsv")
    config_name = value(manifest, key="config_name", default=config_dir.name)
    cond = conditions.get(config_name, {})
    radius_slug = value(cond, key="radius_slug", default=config_dir.parent.name)
    eval_dir = eval_dir_for(root, radius_slug, config_name)
    eval_summary = read_kv(eval_dir / "summary.tsv")
    copytrack = read_copytrack_all(eval_dir / "copytrack_vector_diag.tsv")
    metadata = read_kv(root / "synthetic_pairs" / f"p9016_charm3dg20k_{radius_slug}_contacts.metadata.tsv")

    row = {col: "NA" for col in COLUMNS}
    row["radius_slug"] = radius_slug
    row["radius_pr"] = value(cond, metadata, key="radius_pr")
    row["condition_name"] = value(cond, key="condition_name")
    row["config_name"] = config_name
    row["status"] = status(manifest, eval_dir, eval_summary)
    row["backend"] = value(manifest, key="relax_backend")
    row["output_dir"] = str(config_dir)
    row["eval_dir"] = str(eval_dir)
    row["input_pairs"] = value(cond, key="input_pairs")
    for key in [
        "output_pairs_total",
        "output_pairs_cis",
        "output_pairs_trans",
        "candidate_pairs_total",
    ]:
        row[key] = value(metadata, key=key)
    for key in [
        "init_mode",
        "init_seed",
        "init_eps",
        "init_noise_scale",
        "init_scale",
        "init_eps_effective",
        "init_noise_scale_effective",
        "init_scale_effective",
        "d_scale_mode_input_string",
        "d_scale_mode",
        "dscale_mode",
        "d_scale_posterior_gamma",
        "dscale_posterior_gamma",
        "d_scale_eps_count",
        "dscale_effective_count_formula",
        "dscale_probability_weighted",
        "edge_k_probability_weighted",
        "legacy_expected_count_alias_used",
        "min_sep_unit",
        "lambda_sep",
        "prior_mode",
        "rho_train_mode",
        "temperature_start",
        "temperature_end",
        "uses_phase_labels",
        "uses_charm_or_reference",
        "uses_charm_for_training",
        "reference_derived_positive_control",
        "n_raw",
        "n_bpair",
        "n_raw_cis",
        "n_raw_trans",
        "n_raw_same_bin_excluded",
        "final_mean_entropy",
        "final_mean_pU",
        "final_min_sep",
        "final_mean_sep",
        "final_max_sep",
        "final_contact_energy",
        "final_repulsion_energy",
        "final_backbone_energy",
        "final_sep_energy",
        "final_sep_force_l1",
    ]:
        row[key] = value(manifest, key=key)
    for key in [
        "sep_p01",
        "sep_p05",
        "sep_p10",
        "sep_median",
        "sep_p90",
        "frac_sep_below_min",
    ]:
        row[key] = value(sep, key=key)
    for key in [
        "model_top1_accuracy_genome_all",
        "model_top1_accuracy_genome_cis",
        "model_top1_accuracy_genome_trans",
        "model_same_cross_accuracy_genome_all",
        "model_same_cross_accuracy_genome_cis",
        "model_same_cross_accuracy_genome_trans",
        "model_pmax90_accuracy_genome_all",
        "model_pmax90_accuracy_genome_cis",
        "model_pmax90_accuracy_genome_trans",
        "model_pmax90_recall_genome_all",
        "model_pmax90_recall_genome_cis",
        "model_pmax90_recall_genome_trans",
        "mean_per_chrom_cis_distance_spearman",
        "charm3dg_top1_accuracy_genome_trans",
        "charm3dg_same_cross_accuracy_genome_trans",
    ]:
        row[key] = value(eval_summary, key=key)
    row["copytrack_sep_p05"] = value(copytrack, key="sep_p05")
    row["copytrack_sep_median"] = value(copytrack, key="sep_median")
    row["copytrack_vector_cos_p05"] = value(copytrack, key="vector_cos_p05")
    row["copytrack_vector_cos_median"] = value(copytrack, key="vector_cos_median")
    row["copytrack_frac_cos_lt_0"] = value(copytrack, key="frac_cos_lt_0")
    row["copytrack_frac_projection_sign_switch"] = value(copytrack, key="frac_projection_sign_switch")
    return row


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: scripts/summarize_p9016_016.py <FULL_RUN_ROOT>", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    conditions = read_condition_matrix(root / "condition_matrix.tsv")
    writer = csv.DictWriter(sys.stdout, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for config_dir in discover_config_dirs(root):
        writer.writerow(summarize(root, config_dir, conditions))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
