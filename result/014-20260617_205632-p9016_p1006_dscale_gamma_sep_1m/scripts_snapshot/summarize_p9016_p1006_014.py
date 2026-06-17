#!/usr/bin/env python3
"""Summarize the 014 P9016/P1006 gamma and separation sweep."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


COLUMNS = [
    "sample",
    "config_name",
    "status",
    "backend",
    "output_dir",
    "eval_dir",
    "git_commit",
    "git_dirty_count",
    "binary_hash",
    "bin_size_bp",
    "init_mode",
    "init_seed",
    "init_eps",
    "init_noise_scale",
    "prior_mode",
    "rho_train_mode",
    "d_scale_mode",
    "dscale_mode",
    "dscale_effective_count_formula",
    "d_scale_posterior_gamma",
    "dscale_posterior_gamma",
    "dscale_probability_weighted",
    "edge_k_probability_weighted",
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
    "final_contact_energy",
    "final_repulsion_energy",
    "final_backbone_energy",
    "final_sep_energy",
    "final_sep_force_l1",
    "force_diag_total_contact_force_l1",
    "sep_force_l1_over_contact_force_l1",
    "sep_p01",
    "sep_p05",
    "sep_p10",
    "sep_median",
    "sep_p90",
    "n_sep_below_min",
    "frac_sep_below_min",
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
    "copytrack_sep_p05",
    "copytrack_sep_median",
    "copytrack_vector_cos_p05",
    "copytrack_vector_cos_median",
    "copytrack_frac_cos_lt_0",
    "copytrack_frac_projection_sign_switch",
    "copytrack_median_projection_run_length",
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
	outputs = root / "outputs"
	if not outputs.exists():
		return []
	config_dirs: set[Path] = set()
	for sample_dir in outputs.iterdir():
		if not sample_dir.is_dir():
			continue
		for config_dir in sample_dir.iterdir():
			if config_dir.is_dir() and (config_dir / "p9016_full.manifest.tsv").exists():
				config_dirs.add(config_dir)
	config_dirs.update({path.parent for path in outputs.rglob("p9016_full.manifest.tsv")})
	return sorted(config_dirs, key=lambda p: (p.parent.name, p.name))


def eval_dir_for(root: Path, sample: str, config_name: str) -> Path:
    candidate = root / "eval" / sample / config_name
    if candidate.exists():
        return candidate
    flat = root / "eval" / config_name
    if flat.exists():
        return flat
    return candidate


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


def summarize(root: Path, config_dir: Path) -> dict[str, str]:
    manifest = read_kv(config_dir / "p9016_full.manifest.tsv")
    loop = read_one(config_dir / "p9016_full.loop_diag.tsv")
    sep = read_one(config_dir / "p9016_full.sep_diag.tsv")
    force = read_force_total(config_dir / "p9016_full.force_class_diag.tsv")
    sample = value(manifest, key="sample", default=config_dir.parent.name)
    config_name = value(manifest, key="config_name", default=config_dir.name)
    eval_dir = eval_dir_for(root, sample, config_name)
    eval_summary = read_kv(eval_dir / "summary.tsv")
    copytrack = read_copytrack_all(eval_dir / "copytrack_vector_diag.tsv")

    row = {column: "NA" for column in COLUMNS}
    row["sample"] = sample
    row["config_name"] = config_name
    row["status"] = status(manifest, eval_dir, eval_summary)
    row["backend"] = value(manifest, key="relax_backend")
    row["output_dir"] = str(config_dir)
    row["eval_dir"] = str(eval_dir)
    row["bin_size_bp"] = value(manifest, eval_summary, key="bin_size_bp", default=value(manifest, key="resolution"))
    row["dscale_mode"] = value(manifest, key="dscale_mode", default=value(manifest, key="d_scale_mode"))
    row["dscale_posterior_gamma"] = value(manifest, key="d_scale_posterior_gamma")

    for key in [
        "git_commit",
        "git_dirty_count",
        "binary_hash",
        "init_mode",
        "init_seed",
        "init_eps",
        "init_noise_scale",
        "prior_mode",
        "rho_train_mode",
        "d_scale_mode",
        "dscale_effective_count_formula",
        "d_scale_posterior_gamma",
        "dscale_probability_weighted",
        "edge_k_probability_weighted",
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
        "final_contact_energy",
        "final_repulsion_energy",
        "final_backbone_energy",
        "final_sep_energy",
        "final_sep_force_l1",
    ]:
        row[key] = value(manifest, loop, key=key)

    row["force_diag_total_contact_force_l1"] = value(force, key="contact_force_l1")
    row["sep_force_l1_over_contact_force_l1"] = ratio(row["final_sep_force_l1"], row["force_diag_total_contact_force_l1"])

    for key in ["sep_p01", "sep_p05", "sep_p10", "sep_median", "sep_p90", "n_sep_below_min", "frac_sep_below_min"]:
        row[key] = value(sep, key=key)

    for key in [
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
    ]:
        row[key] = value(eval_summary, key=key)

    copy_map = {
        "copytrack_sep_p05": "sep_p05",
        "copytrack_sep_median": "sep_median",
        "copytrack_vector_cos_p05": "vector_cos_p05",
        "copytrack_vector_cos_median": "vector_cos_median",
        "copytrack_frac_cos_lt_0": "frac_cos_lt_0",
        "copytrack_frac_projection_sign_switch": "frac_projection_sign_switch",
        "copytrack_median_projection_run_length": "median_projection_run_length",
    }
    for out_key, in_key in copy_map.items():
        row[out_key] = value(copytrack, key=in_key)
    return row


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: summarize_p9016_p1006_014.py FULL_RUN_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    writer = csv.DictWriter(sys.stdout, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for config_dir in discover_config_dirs(root):
        writer.writerow(summarize(root, config_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
