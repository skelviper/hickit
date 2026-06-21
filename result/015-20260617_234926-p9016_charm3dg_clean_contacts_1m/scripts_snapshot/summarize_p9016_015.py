#!/usr/bin/env python3
"""Summarize experiment 015 CHARM-derived clean-contact positive controls."""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path


COLUMNS = [
    "config_name",
    "status",
    "radius_pr",
    "particle_radius",
    "radius_distance",
    "target_mode",
    "target_count_requested",
    "output_pairs_total",
    "output_pairs_cis",
    "output_pairs_trans",
    "candidate_pairs_total",
    "sampling_with_replacement",
    "backend",
    "input_contact_source",
    "uses_charm_or_reference",
    "uses_charm_for_training",
    "reference_derived_positive_control",
    "n_raw",
    "n_bpair",
    "n_raw_cis",
    "n_raw_trans",
    "n_raw_same_bin_excluded",
    "d_scale_mode",
    "d_scale_posterior_gamma",
    "final_mean_entropy",
    "final_mean_pU",
    "final_min_sep",
    "final_mean_sep",
    "sep_p05",
    "sep_median",
    "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans",
    "model_same_cross_accuracy_genome_cis",
    "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "charm3dg_top1_accuracy_genome_trans",
    "charm3dg_same_cross_accuracy_genome_trans",
    "copytrack_frac_cos_lt_0",
    "copytrack_frac_projection_sign_switch",
    "output_dir",
    "eval_dir",
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


def value(*sources: dict[str, str], key: str, default: str = "NA") -> str:
    for source in sources:
        v = source.get(key)
        if v not in (None, ""):
            return v
    return default


def config_dirs(root: Path) -> list[Path]:
    outputs = root / "outputs"
    if not outputs.exists():
        return []
    dirs: set[Path] = set()
    for dirpath, _, filenames in os.walk(outputs, followlinks=True):
        if "p9016_full.manifest.tsv" in filenames:
            dirs.add(Path(dirpath))
    return sorted(dirs, key=lambda p: p.name)


def summarize(root: Path, config_dir: Path) -> dict[str, str]:
    manifest = read_kv(config_dir / "p9016_full.manifest.tsv")
    sep = read_one(config_dir / "p9016_full.sep_diag.tsv")
    config = value(manifest, key="config_name", default=config_dir.name)
    eval_dir = root / "eval" / config
    eval_summary = read_kv(eval_dir / "summary.tsv")
    copytrack = read_copytrack_all(eval_dir / "copytrack_vector_diag.tsv")
    metadata = read_kv(root / "synthetic_pairs" / f"{config}.metadata.tsv")

    row = {col: "NA" for col in COLUMNS}
    row["config_name"] = config
    row["status"] = value(manifest, key="status", default="MISSING_MANIFEST")
    row["output_dir"] = str(config_dir)
    row["eval_dir"] = str(eval_dir)
    for key in [
        "radius_pr",
        "particle_radius",
        "radius_distance",
        "target_mode",
        "target_count_requested",
        "output_pairs_total",
        "output_pairs_cis",
        "output_pairs_trans",
        "candidate_pairs_total",
        "sampling_with_replacement",
    ]:
        row[key] = value(metadata, key=key)
    for key in [
        "n_raw",
        "n_bpair",
        "n_raw_cis",
        "n_raw_trans",
        "n_raw_same_bin_excluded",
        "input_contact_source",
        "uses_charm_or_reference",
        "uses_charm_for_training",
        "reference_derived_positive_control",
        "d_scale_mode",
        "d_scale_posterior_gamma",
        "final_mean_entropy",
        "final_mean_pU",
        "final_min_sep",
        "final_mean_sep",
    ]:
        row[key] = value(manifest, key=key)
    row["backend"] = value(manifest, key="relax_backend")
    row["sep_p05"] = value(sep, key="sep_p05")
    row["sep_median"] = value(sep, key="sep_median")
    for key in [
        "model_top1_accuracy_genome_all",
        "model_top1_accuracy_genome_cis",
        "model_top1_accuracy_genome_trans",
        "model_same_cross_accuracy_genome_cis",
        "model_same_cross_accuracy_genome_trans",
        "mean_per_chrom_cis_distance_spearman",
        "charm3dg_top1_accuracy_genome_trans",
        "charm3dg_same_cross_accuracy_genome_trans",
    ]:
        row[key] = value(eval_summary, key=key)
    row["copytrack_frac_cos_lt_0"] = value(copytrack, key="frac_cos_lt_0")
    row["copytrack_frac_projection_sign_switch"] = value(copytrack, key="frac_projection_sign_switch")
    return row


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: scripts/summarize_p9016_015.py <FULL_RUN_ROOT>", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    writer = csv.DictWriter(sys.stdout, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for config_dir in config_dirs(root):
        writer.writerow(summarize(root, config_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
