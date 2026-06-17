#!/usr/bin/env python3
"""Compare current 011 legacy baselines against historical 1 Mb scaffold runs."""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


CONFIG_REFS = [
    ("legacy_1m_cpu_exact_ieps0p5_noise0_seed17", "004-20260616_123920-p9016_softall_1m_direct_cpu_baseline"),
    ("legacy_1m_gpu_exact_ieps0p5_noise0_seed17", "003-20260616_125639-p9016_softall_1m_direct_gpu_baseline"),
]

METRICS = [
    ("model_top1_accuracy_genome_all", 0.01),
    ("model_top1_accuracy_genome_cis", 0.01),
    ("model_top1_accuracy_genome_trans", 0.01),
    ("model_same_cross_accuracy_genome_cis", 0.02),
    ("mean_per_chrom_cis_distance_spearman", 0.03),
    ("final_mean_entropy", 0.05),
    ("final_mean_pU", 0.05),
    ("final_mean_sep", 0.5),
    ("final_min_sep", 0.5),
]

MANIFEST_CHECKS = [
    ("init_mode", "unphased_scaffold_split"),
    ("init_seed", "17"),
    ("init_eps_effective", "0.5"),
    ("init_noise_scale_effective", "0"),
    ("prior_mode", "uniform"),
    ("rho_train_mode", "constant"),
    ("d_scale_mode", "raw_count"),
    ("min_sep_unit", "0"),
    ("lambda_sep", "0"),
    ("temperature_start", "1"),
    ("temperature_end", "1"),
    ("uses_phase_labels", "0"),
    ("uses_charm_or_reference", "0"),
]

FALLBACK = {
    "model_top1_accuracy_genome_all": "0.471",
    "model_top1_accuracy_genome_cis": "0.569",
    "model_top1_accuracy_genome_trans": "0.330",
    "mean_per_chrom_cis_distance_spearman": "0.792",
    "final_mean_entropy": "0.608",
    "final_mean_pU": "0.439",
    "final_mean_sep": "4.977",
}


def read_kv(path: Path) -> dict[str, str]:
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


def read_table(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def find_historical_summary(repo_root: Path, ref_name: str) -> tuple[Path | None, str]:
    candidates = [
        repo_root / "result" / ref_name / "summary.tsv",
        repo_root / "hickit" / "result" / ref_name / "summary.tsv",
        repo_root / "test_res" / ref_name / "summary.tsv",
    ]
    for path in candidates:
        if path.exists():
            return path, "historical_summary"
    return None, "fallback_expected_values"


def as_float(value: str) -> float | None:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(x):
        return None
    return x


def compare_numeric(current: str, historical: str, tolerance: float) -> tuple[str, str]:
    c = as_float(current)
    h = as_float(historical)
    if c is None or h is None:
        return "NA", "NA"
    delta = c - h
    return f"{delta:.9g}", "PASS" if abs(delta) <= tolerance else "FAIL"


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: audit_p9016_baseline_regression.py FULL_RUN_ROOT", file=sys.stderr)
        return 2
    full_root = Path(sys.argv[1]).resolve()
    if not full_root.is_dir():
        print(f"error: not a directory: {full_root}", file=sys.stderr)
        return 2
    repo_root = Path(__file__).resolve().parents[1]
    rows = read_table(full_root / "summary.tsv")
    by_config = {row.get("config_name", ""): row for row in rows}

    out_fields = [
        "current_config",
        "historical_ref",
        "source",
        "metric",
        "current_value",
        "historical_value",
        "delta",
        "tolerance",
        "pass_fail",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=out_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()

    for config, ref_name in CONFIG_REFS:
        current = by_config.get(config)
        hist_path, source = find_historical_summary(repo_root, ref_name)
        historical = read_kv(hist_path) if hist_path else FALLBACK
        if current is None:
            writer.writerow({
                "current_config": config,
                "historical_ref": ref_name,
                "source": source,
                "metric": "config_present",
                "current_value": "MISSING",
                "historical_value": "present",
                "delta": "NA",
                "tolerance": "0",
                "pass_fail": "SKIP",
            })
            continue

        for metric, tolerance in METRICS:
            current_value = current.get(metric, "NA")
            historical_value = historical.get(metric, "NA")
            delta, pass_fail = compare_numeric(current_value, historical_value, tolerance)
            if historical_value == "NA":
                pass_fail = "NA"
            writer.writerow({
                "current_config": config,
                "historical_ref": ref_name,
                "source": source,
                "metric": metric,
                "current_value": current_value,
                "historical_value": historical_value,
                "delta": delta,
                "tolerance": f"{tolerance:.9g}",
                "pass_fail": pass_fail,
            })

        for metric, expected in MANIFEST_CHECKS:
            current_value = current.get(metric, "NA")
            writer.writerow({
                "current_config": config,
                "historical_ref": ref_name,
                "source": "expected_manifest_field",
                "metric": metric,
                "current_value": current_value,
                "historical_value": expected,
                "delta": "NA",
                "tolerance": "exact",
                "pass_fail": "PASS" if current_value == expected else "FAIL",
            })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
