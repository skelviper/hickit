#!/usr/bin/env python3
"""Merge P9016 blind-model training, eval, and copy-track diagnostics."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


COLUMNS = [
    "config_name",
    "sample",
    "status",
    "output_dir",
    "eval_dir",
    "git_commit",
    "git_dirty_count",
    "binary_hash",
    "backend",
    "bin_size_bp",
    "resolution_label",
    "refinement_stage",
    "refinement_stage_index",
    "refinement_n_stages",
    "parent_bin_size_bp",
    "refinement_init_source",
    "n_raw",
    "n_raw_cis",
    "n_raw_trans",
    "n_bpair",
    "n_bpair_cis",
    "n_bpair_trans",
    "init_mode",
    "init_seed",
    "init_eps",
    "init_noise_scale",
    "init_scale",
    "init_eps_effective",
    "init_noise_scale_effective",
    "init_scale_effective",
    "prior_mode",
    "trans_chr_pair_prior_mode",
    "trans_chr_pair_prior_lambda",
    "trans_chr_pair_prior_eps",
    "trans_chr_pair_prior_power",
    "trans_chr_pair_prior_warmup_iter",
    "trans_chr_pair_mstep_mode",
    "trans_chr_pair_mstep_lambda",
    "trans_chr_pair_mstep_eps",
    "trans_chr_pair_mstep_power",
    "trans_chr_pair_mstep_warmup_iter",
    "trans_chr_pair_mstep_application_point",
    "rho_train_mode",
    "rho_train",
    "rho_train_floor",
    "rho_train_start",
    "rho_train_end",
    "temperature_start",
    "temperature_end",
    "min_sep_unit",
    "lambda_sep",
    "lambda_copytrack",
    "lambda_global_copytrack",
    "lambda_normdir_copytrack",
    "normdir_copytrack_eps_unit",
    "normdir_copytrack_prior_mode",
    "global_copytrack_prior_mode",
    "global_copytrack_weight_mode",
    "d_scale_mode",
    "d_scale_mode_input_string",
    "d_scale_eps_count",
    "d_scale_posterior_gamma",
    "dscale_effective_count_formula",
    "edge_k_probability_weighted",
    "dscale_probability_weighted",
    "legacy_expected_count_alias_used",
    "trans_contact_scaling_mode",
    "trans_k_multiplier",
    "trans_dscale_multiplier",
    "trans_top1_mstep_mode",
    "trans_top1_mstep_min_pmax",
    "trans_top1_mstep_min_margin",
    "trans_top1_mstep_mix_weight",
    "trans_gate_mode",
    "trans_gate_min_pmax",
    "trans_gate_min_margin",
    "trans_gate_min_neg_entropy",
    "trans_gate_scope",
    "trans_gate_source",
    "trans_gate_confidence_unit",
    "trans_gate_neg_entropy_formula",
    "trans_gate_active_criteria",
    "input_contact_source",
    "readgroup_mode",
    "readgroup_max_segments",
    "readgroup_eps",
    "readgroup_entries",
    "readgroup_groups",
    "readgroup_groups_used",
    "readgroup_raw_decoded",
    "readgroup_groups_skipped_too_large",
    "readgroup_groups_skipped_bad",
    "readgroup_uses_phase_labels",
    "readgroup_uses_charm_or_reference",
    "uses_phase_labels",
    "uses_charm_or_reference",
    "final_mean_sep",
    "final_min_sep",
    "final_max_sep",
    "sep_p05",
    "sep_median",
    "final_mean_entropy",
    "final_mean_pU",
    "final_n_softall_gate_skip_raw",
    "final_contact_energy",
    "final_repulsion_energy",
    "final_backbone_energy",
    "final_sep_energy",
    "final_sep_force_l1",
    "final_copytrack_energy",
    "final_copytrack_force_l1",
    "final_global_copytrack_energy",
    "final_global_copytrack_force_l1",
    "final_normdir_copytrack_energy",
    "final_normdir_copytrack_force_l1",
    "force_diag_total_contact_force_l1",
    "sep_force_l1_over_contact_force_l1",
    "copytrack_force_l1_over_contact_force_l1",
    "global_copytrack_force_l1_over_contact_force_l1",
    "normdir_copytrack_force_l1_over_contact_force_l1",
    "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans",
    "model_same_cross_accuracy_genome_all",
    "model_same_cross_accuracy_genome_cis",
    "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "copytrack_sep_p05",
    "copytrack_sep_median",
    "copytrack_vector_cos_p05",
    "copytrack_vector_cos_median",
    "copytrack_frac_cos_lt_0",
    "copytrack_frac_projection_sign_switch",
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


def parse_int(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def discover_config_dirs(root: Path) -> list[Path]:
    output_root = root / "outputs"
    if output_root.exists():
        config_dirs: set[Path] = set()
        for child in output_root.iterdir():
            if child.is_dir() and (child / "p9016_full.manifest.tsv").exists():
                config_dirs.add(child)
        config_dirs.update({path.parent for path in output_root.rglob("p9016_full.manifest.tsv")})
        return sorted(config_dirs, key=lambda p: p.name)
    return sorted({path.parent for path in root.rglob("p9016_full.manifest.tsv")}, key=lambda p: str(p))


def find_eval_dir(root: Path, config_name: str) -> Path:
    candidate = root / "eval" / config_name
    if candidate.exists():
        return candidate
    return candidate


def status_from_files(manifest: dict[str, str], eval_dir: Path, copytrack: dict[str, str]) -> str:
    eval_status = eval_dir / "eval_status.tsv"
    if eval_status.exists():
        with eval_status.open() as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 2 and parts[0] == "status":
                    return parts[1]
    train_status = value(manifest, key="status", default="NA")
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
    stage_index = parse_int(value(manifest, key="refinement_stage_index"))
    n_stages = parse_int(value(manifest, key="refinement_n_stages"))
    is_nonfinal_stage = (
        stage_index is not None
        and n_stages is not None
        and n_stages > 1
        and stage_index < n_stages
    )
    if is_nonfinal_stage:
        eval_dir = root / "eval" / f"{config_name}__stage{stage_index}_not_evaled"
        eval_summary = {}
        copytrack = {}
    else:
        eval_dir = find_eval_dir(root, config_name)
        eval_summary = read_key_value_tsv(eval_dir / "summary.tsv")
        copytrack = read_copytrack_all(eval_dir / "copytrack_vector_diag.tsv")

    row = {column: "NA" for column in COLUMNS}
    row["config_name"] = config_name
    row["sample"] = value(manifest, eval_summary, key="sample")
    row["output_dir"] = str(config_dir)
    row["eval_dir"] = "NA" if is_nonfinal_stage else str(eval_dir)
    row["status"] = status_from_files(manifest, eval_dir, copytrack)
    row["backend"] = value(manifest, key="relax_backend")
    row["bin_size_bp"] = value(manifest, eval_summary, key="bin_size_bp", default=value(manifest, key="resolution"))

    for key in [
        "git_commit",
        "git_dirty_count",
        "binary_hash",
        "n_raw",
        "resolution_label",
        "refinement_stage",
        "refinement_stage_index",
        "refinement_n_stages",
        "parent_bin_size_bp",
        "refinement_init_source",
        "n_raw",
        "n_raw_cis",
        "n_raw_trans",
        "n_bpair",
        "n_bpair_cis",
        "n_bpair_trans",
        "init_mode",
        "init_seed",
        "init_eps",
        "init_noise_scale",
        "init_scale",
        "init_eps_effective",
        "init_noise_scale_effective",
        "init_scale_effective",
        "prior_mode",
        "trans_chr_pair_prior_mode",
        "trans_chr_pair_prior_lambda",
        "trans_chr_pair_prior_eps",
        "trans_chr_pair_prior_power",
        "trans_chr_pair_prior_warmup_iter",
        "trans_chr_pair_mstep_mode",
        "trans_chr_pair_mstep_lambda",
        "trans_chr_pair_mstep_eps",
        "trans_chr_pair_mstep_power",
        "trans_chr_pair_mstep_warmup_iter",
        "trans_chr_pair_mstep_application_point",
        "rho_train_mode",
        "rho_train",
        "rho_train_floor",
        "rho_train_start",
        "rho_train_end",
        "temperature_start",
        "temperature_end",
        "min_sep_unit",
        "lambda_sep",
        "lambda_copytrack",
        "lambda_global_copytrack",
        "lambda_normdir_copytrack",
        "normdir_copytrack_eps_unit",
        "normdir_copytrack_prior_mode",
        "global_copytrack_prior_mode",
        "global_copytrack_weight_mode",
        "d_scale_mode",
        "d_scale_mode_input_string",
        "d_scale_eps_count",
        "d_scale_posterior_gamma",
        "dscale_effective_count_formula",
        "edge_k_probability_weighted",
        "dscale_probability_weighted",
        "legacy_expected_count_alias_used",
        "trans_contact_scaling_mode",
        "trans_k_multiplier",
        "trans_dscale_multiplier",
        "trans_top1_mstep_mode",
        "trans_top1_mstep_min_pmax",
        "trans_top1_mstep_min_margin",
        "trans_top1_mstep_mix_weight",
        "trans_gate_mode",
        "trans_gate_min_pmax",
        "trans_gate_min_margin",
        "trans_gate_min_neg_entropy",
        "trans_gate_scope",
        "trans_gate_source",
        "trans_gate_confidence_unit",
        "trans_gate_neg_entropy_formula",
        "trans_gate_active_criteria",
        "input_contact_source",
        "readgroup_mode",
        "readgroup_max_segments",
        "readgroup_eps",
        "readgroup_entries",
        "readgroup_groups",
        "readgroup_groups_used",
        "readgroup_raw_decoded",
        "readgroup_groups_skipped_too_large",
        "readgroup_groups_skipped_bad",
        "readgroup_uses_phase_labels",
        "readgroup_uses_charm_or_reference",
        "uses_phase_labels",
        "uses_charm_or_reference",
        "final_mean_sep",
        "final_min_sep",
        "final_max_sep",
        "final_mean_entropy",
        "final_mean_pU",
        "final_n_softall_gate_skip_raw",
        "final_contact_energy",
        "final_repulsion_energy",
        "final_backbone_energy",
        "final_sep_energy",
        "final_sep_force_l1",
        "final_copytrack_energy",
        "final_copytrack_force_l1",
        "final_global_copytrack_energy",
        "final_global_copytrack_force_l1",
        "final_normdir_copytrack_energy",
        "final_normdir_copytrack_force_l1",
    ]:
        row[key] = value(manifest, loop_diag, key=key)

    for key in ["sep_p05", "sep_median"]:
        row[key] = value(sep_diag, key=key)

    row["force_diag_total_contact_force_l1"] = value(force_total, key="contact_force_l1")
    row["sep_force_l1_over_contact_force_l1"] = ratio(
        row["final_sep_force_l1"],
        row["force_diag_total_contact_force_l1"],
    )
    row["copytrack_force_l1_over_contact_force_l1"] = ratio(
        row["final_copytrack_force_l1"],
        row["force_diag_total_contact_force_l1"],
    )
    row["global_copytrack_force_l1_over_contact_force_l1"] = ratio(
        row["final_global_copytrack_force_l1"],
        row["force_diag_total_contact_force_l1"],
    )
    row["normdir_copytrack_force_l1_over_contact_force_l1"] = ratio(
        row["final_normdir_copytrack_force_l1"],
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
    ]:
        row[key] = value(eval_summary, key=key)

    copy_map = {
        "copytrack_sep_p05": "sep_p05",
        "copytrack_sep_median": "sep_median",
        "copytrack_vector_cos_p05": "vector_cos_p05",
        "copytrack_vector_cos_median": "vector_cos_median",
        "copytrack_frac_cos_lt_0": "frac_cos_lt_0",
        "copytrack_frac_projection_sign_switch": "frac_projection_sign_switch",
    }
    for out_key, in_key in copy_map.items():
        row[out_key] = value(copytrack, key=in_key)
    return row


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: summarize_p9016_model_sweep.py FULL_RUN_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    if not root.is_dir():
        print(f"error: not a directory: {root}", file=sys.stderr)
        return 2
    writer = csv.DictWriter(sys.stdout, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for config_dir in discover_config_dirs(root):
        writer.writerow(summarize_config(root, config_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
