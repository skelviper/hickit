#!/usr/bin/env python3
"""Merge P9016 blind-model training, eval, and copy-track diagnostics."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


COLUMNS = [
    "config_name",
    "status",
    "output_dir",
    "eval_dir",
    "git_commit",
    "binary_hash",
    "backend",
    "bin_size_bp",
    "init_mode",
    "init_seed",
    "init_eps",
    "init_noise_scale",
    "init_scale",
    "init_eps_effective",
    "init_noise_scale_effective",
    "init_scale_effective",
    "prior_mode",
    "rho_train_mode",
    "d_scale_mode",
    "d_scale_eps_count",
    "min_sep_unit",
    "lambda_sep",
    "temperature_start",
    "temperature_end",
    "uses_phase_labels",
    "uses_charm_or_reference",
    "final_mean_sep",
    "final_min_sep",
    "final_max_sep",
    "sep_p01",
    "sep_p05",
    "sep_p10",
    "sep_median",
    "sep_p90",
    "n_sep_below_min",
    "frac_sep_below_min",
    "final_mean_entropy",
    "final_mean_pU",
    "final_sum_wedge_k",
    "last_training_sum_wedge_k",
    "final_refreshed_sum_wedge_k",
    "final_refreshed_n_wedges",
    "final_mean_rho_train_bpair",
    "final_min_rho_train_bpair",
    "final_max_rho_train_bpair",
    "final_contact_energy",
    "final_repulsion_energy",
    "final_backbone_energy",
    "final_sep_energy",
    "final_sep_force_l1",
    "force_diag_total_contact_energy",
    "force_diag_total_contact_force_l1",
    "sep_force_l1_over_contact_force_l1",
    "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans",
    "model_same_cross_accuracy_genome_all",
    "model_same_cross_accuracy_genome_cis",
    "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "p90_accuracy_genome_all",
    "p90_recall_genome_all",
    "model_top1_cis_recon_mean_distance",
    "truth_or_ref_top1_cis_mean_distance",
    "charm3dg_posterior_top1_contact_distance_genome_cis_mean",
    "reconstruction_posterior_top1_contact_distance_genome_cis_mean",
    "charm3dg_posterior_top1_contact_distance_genome_trans_mean",
    "reconstruction_posterior_top1_contact_distance_genome_trans_mean",
    "copytrack_sep_p05",
    "copytrack_sep_median",
    "copytrack_vector_cos_p05",
    "copytrack_vector_cos_median",
    "copytrack_frac_cos_lt_0",
    "copytrack_frac_cos_lt_neg0p5",
    "copytrack_frac_projection_sign_switch",
    "copytrack_median_projection_run_length",
]


def read_key_value_tsv(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header != ["key", "value"]:
            return values
        for row in reader:
            if len(row) >= 2:
                values[row[0]] = row[1]
    return values


def read_single_row_tsv(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        row = next(reader, None)
        return dict(row) if row else {}


def read_force_total(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("class") == "total":
                return dict(row)
    return {}


def read_copytrack_all(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("chrom") == "ALL":
                return dict(row)
    return {}


def value(*sources: dict[str, str], key: str, default: str = "NA") -> str:
    for source in sources:
        v = source.get(key)
        if v not in (None, ""):
            return v
    return default


def ratio(num: str, den: str) -> str:
    try:
        n = float(num)
        d = float(den)
    except ValueError:
        return "NA"
    if d == 0.0:
        return "NA"
    return f"{n / d:.9g}"


def discover_config_dirs(root: Path) -> list[Path]:
    return sorted({path.parent for path in root.rglob("p9016_full.manifest.tsv")}, key=lambda p: str(p))


def find_eval_dir(root: Path, config_name: str) -> Path:
    candidates = [
        root / "eval" / config_name,
        root.parent / "eval" / config_name,
    ]
    for path in candidates:
        if path.exists():
            return path
    matches = sorted(root.rglob(f"eval/**/{config_name}"))
    if matches:
        return matches[0]
    return root / "eval" / config_name


def status_from_files(manifest: dict[str, str], eval_dir: Path, copytrack: dict[str, str]) -> str:
    train_status = value(manifest, key="status", default="NA")
    eval_status_path = eval_dir / "eval_status.tsv"
    if eval_status_path.exists():
        with eval_status_path.open() as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 2 and parts[0] == "status":
                    return parts[1]
    if train_status not in ("NA", "OK"):
        return train_status
    if (eval_dir / "summary.tsv").exists() or copytrack:
        return "OK"
    return train_status


def summarize_config(root: Path, config_dir: Path) -> dict[str, str]:
    manifest = read_key_value_tsv(config_dir / "p9016_full.manifest.tsv")
    loop_diag = read_single_row_tsv(config_dir / "p9016_full.loop_diag.tsv")
    force_total = read_force_total(config_dir / "p9016_full.force_class_diag.tsv")
    sep_diag = read_single_row_tsv(config_dir / "p9016_full.sep_diag.tsv")
    config_name = value(manifest, key="config_name", default=config_dir.name)
    eval_dir = find_eval_dir(root, config_name)
    eval_summary = read_key_value_tsv(eval_dir / "summary.tsv")
    copytrack = read_copytrack_all(eval_dir / "copytrack_vector_diag.tsv")

    row = {column: "NA" for column in COLUMNS}
    row["config_name"] = config_name
    row["output_dir"] = str(config_dir)
    row["eval_dir"] = str(eval_dir)
    row["status"] = status_from_files(manifest, eval_dir, copytrack)

    row["backend"] = value(manifest, key="relax_backend")
    row["bin_size_bp"] = value(manifest, eval_summary, key="bin_size_bp", default=value(manifest, key="resolution"))

    for key in [
        "git_commit",
        "binary_hash",
        "init_mode",
        "init_seed",
        "init_eps",
        "init_noise_scale",
        "init_scale",
        "init_eps_effective",
        "init_noise_scale_effective",
        "init_scale_effective",
        "prior_mode",
        "rho_train_mode",
        "d_scale_mode",
        "d_scale_eps_count",
        "min_sep_unit",
        "lambda_sep",
        "temperature_start",
        "temperature_end",
        "uses_phase_labels",
        "uses_charm_or_reference",
        "final_mean_sep",
        "final_min_sep",
        "final_max_sep",
        "final_mean_entropy",
        "final_mean_pU",
        "final_sum_wedge_k",
        "last_training_sum_wedge_k",
        "final_refreshed_sum_wedge_k",
        "final_refreshed_n_wedges",
        "final_mean_rho_train_bpair",
        "final_min_rho_train_bpair",
        "final_max_rho_train_bpair",
        "final_contact_energy",
        "final_repulsion_energy",
        "final_backbone_energy",
        "final_sep_energy",
        "final_sep_force_l1",
    ]:
        row[key] = value(manifest, loop_diag, key=key)

    for key in [
        "sep_p01",
        "sep_p05",
        "sep_p10",
        "sep_median",
        "sep_p90",
        "n_sep_below_min",
        "frac_sep_below_min",
    ]:
        row[key] = value(sep_diag, key=key)

    row["force_diag_total_contact_energy"] = value(force_total, key="contact_energy")
    row["force_diag_total_contact_force_l1"] = value(force_total, key="contact_force_l1")
    row["sep_force_l1_over_contact_force_l1"] = ratio(
        value(manifest, loop_diag, force_total, key="final_sep_force_l1", default=value(force_total, key="homolog_sep_force_l1")),
        row["force_diag_total_contact_force_l1"],
    )

    for key in [
        "model_top1_accuracy_genome_all",
        "model_top1_accuracy_genome_cis",
        "model_top1_accuracy_genome_trans",
        "model_same_cross_accuracy_genome_all",
        "model_same_cross_accuracy_genome_cis",
        "model_same_cross_accuracy_genome_trans",
        "mean_per_chrom_cis_distance_spearman",
        "charm3dg_posterior_top1_contact_distance_genome_cis_mean",
        "reconstruction_posterior_top1_contact_distance_genome_cis_mean",
        "charm3dg_posterior_top1_contact_distance_genome_trans_mean",
        "reconstruction_posterior_top1_contact_distance_genome_trans_mean",
    ]:
        row[key] = value(eval_summary, key=key)

    row["p90_accuracy_genome_all"] = value(eval_summary, key="model_pmax90_accuracy_genome_all")
    row["p90_recall_genome_all"] = value(eval_summary, key="model_pmax90_recall_genome_all")
    row["model_top1_cis_recon_mean_distance"] = value(
        eval_summary,
        key="reconstruction_posterior_top1_contact_distance_genome_cis_mean",
    )
    row["truth_or_ref_top1_cis_mean_distance"] = value(
        eval_summary,
        key="charm3dg_posterior_top1_contact_distance_genome_cis_mean",
    )

    row["copytrack_sep_p05"] = value(copytrack, key="sep_p05")
    row["copytrack_sep_median"] = value(copytrack, key="sep_median")
    row["copytrack_vector_cos_p05"] = value(copytrack, key="vector_cos_p05")
    row["copytrack_vector_cos_median"] = value(copytrack, key="vector_cos_median")
    row["copytrack_frac_cos_lt_0"] = value(copytrack, key="frac_cos_lt_0")
    row["copytrack_frac_cos_lt_neg0p5"] = value(copytrack, key="frac_cos_lt_neg0p5")
    row["copytrack_frac_projection_sign_switch"] = value(copytrack, key="frac_projection_sign_switch")
    row["copytrack_median_projection_run_length"] = value(copytrack, key="median_projection_run_length")
    return row


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: summarize_p9016_model_sweep.py FULL_RUN_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    if not root.is_dir():
        print(f"error: not a directory: {root}", file=sys.stderr)
        return 2
    config_dirs = discover_config_dirs(root)
    writer = csv.DictWriter(sys.stdout, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for config_dir in config_dirs:
        writer.writerow(summarize_config(root, config_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
