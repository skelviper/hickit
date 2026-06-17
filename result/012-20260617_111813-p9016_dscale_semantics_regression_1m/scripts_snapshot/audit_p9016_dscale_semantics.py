#!/usr/bin/env python3
"""Audit whether posterior-count d_scale semantics explain the 1 Mb baseline regression."""

from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from pathlib import Path


BACKENDS = {
    "cpu": {
        "raw": "current_raw_cpu_exact_ieps0p5_noise0_seed17",
        "posterior": "posterior_count_cpu_exact_ieps0p5_noise0_seed17",
        "historical_ref": "004-20260616_123920-p9016_softall_1m_direct_cpu_baseline",
    },
    "gpu": {
        "raw": "current_raw_gpu_exact_ieps0p5_noise0_seed17",
        "posterior": "posterior_count_gpu_exact_ieps0p5_noise0_seed17",
        "historical_ref": "003-20260616_125639-p9016_softall_1m_direct_gpu_baseline",
    },
}

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

FALLBACK = {
    "model_top1_accuracy_genome_all": "0.471",
    "model_top1_accuracy_genome_cis": "0.569",
    "model_top1_accuracy_genome_trans": "0.330",
    "model_same_cross_accuracy_genome_cis": "0.969",
    "mean_per_chrom_cis_distance_spearman": "0.792",
    "final_mean_entropy": "0.608",
    "final_mean_pU": "0.439",
    "final_mean_sep": "4.977",
    "final_min_sep": "1.5",
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


def find_historical_summary(repo_root: Path, ref_name: str) -> tuple[Path | None, dict[str, str]]:
    candidates = [
        repo_root / "result" / ref_name / "summary.tsv",
        repo_root.parent / "result" / ref_name / "summary.tsv",
        repo_root.parent / "test_res" / ref_name / "summary.tsv",
    ]
    for path in candidates:
        if path.exists():
            return path, read_kv(path)
    return None, dict(FALLBACK)


def as_float(value: str | None) -> float | None:
    try:
        x = float(value) if value is not None else math.nan
    except ValueError:
        return None
    if not math.isfinite(x):
        return None
    return x


def classify(raw_value: str, posterior_value: str, historical_value: str, tol: float) -> tuple[str, str, str, str, str]:
    raw = as_float(raw_value)
    posterior = as_float(posterior_value)
    historical = as_float(historical_value)
    if raw is None or posterior is None or historical is None:
        return "NA", "NA", "NA", "NA", "INSUFFICIENT_DATA"
    raw_delta = raw - historical
    posterior_delta = posterior - historical
    raw_match = abs(raw_delta) <= tol
    posterior_match = abs(posterior_delta) <= tol
    posterior_closer = abs(posterior_delta) < abs(raw_delta)
    if posterior_match and not raw_match:
        conclusion = "HISTORICAL_D_SCALE_SEMANTICS_CONFIRMED"
    elif raw_match and posterior_match:
        conclusion = "BOTH_MATCH"
    elif (not raw_match) and (not posterior_match):
        conclusion = "BOTH_FAIL"
    else:
        conclusion = "HISTORICAL_D_SCALE_SEMANTICS_NOT_CONFIRMED"
    return (
        f"{raw_delta:.9g}",
        f"{posterior_delta:.9g}",
        "1" if posterior_closer else "0",
        "1" if posterior_match else "0",
        conclusion,
    )


def headline_from_rows(rows: list[dict[str, str]]) -> str:
    conclusions = [r["conclusion"] for r in rows if r["conclusion"] != "INSUFFICIENT_DATA"]
    if not conclusions:
        return "INSUFFICIENT_DATA"
    counts = Counter(conclusions)
    if counts["HISTORICAL_D_SCALE_SEMANTICS_CONFIRMED"] >= 3:
        return "HISTORICAL_D_SCALE_SEMANTICS_CONFIRMED"
    if counts["BOTH_MATCH"] >= 3:
        return "BOTH_MATCH"
    if counts["BOTH_FAIL"] >= 3:
        return "BOTH_FAIL"
    return "HISTORICAL_D_SCALE_SEMANTICS_NOT_CONFIRMED"


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: audit_p9016_dscale_semantics.py FULL_RUN_ROOT", file=sys.stderr)
        return 2
    full_root = Path(sys.argv[1]).resolve()
    if not full_root.is_dir():
        print(f"error: not a directory: {full_root}", file=sys.stderr)
        return 2
    repo_root = Path(__file__).resolve().parents[1]
    rows = read_table(full_root / "summary.tsv")
    by_config = {row.get("config_name", ""): row for row in rows}

    out_fields = [
        "comparison",
        "backend",
        "metric",
        "raw_value",
        "posterior_count_value",
        "historical_value",
        "raw_delta",
        "posterior_count_delta",
        "posterior_count_closer",
        "tolerance",
        "conclusion",
        "historical_ref",
        "historical_source",
    ]
    out_rows: list[dict[str, str]] = []
    for backend, spec in BACKENDS.items():
        raw_row = by_config.get(spec["raw"])
        posterior_row = by_config.get(spec["posterior"])
        hist_path, historical = find_historical_summary(repo_root, spec["historical_ref"])
        source = str(hist_path) if hist_path else "fallback_expected_values"
        for metric, tol in METRICS:
            raw_value = raw_row.get(metric, "NA") if raw_row else "NA"
            posterior_value = posterior_row.get(metric, "NA") if posterior_row else "NA"
            historical_value = historical.get(metric, "NA")
            raw_delta, posterior_delta, closer, _posterior_match, conclusion = classify(
                raw_value, posterior_value, historical_value, tol
            )
            out_rows.append({
                "comparison": f"{spec['raw']}__vs__{spec['posterior']}",
                "backend": backend,
                "metric": metric,
                "raw_value": raw_value,
                "posterior_count_value": posterior_value,
                "historical_value": historical_value,
                "raw_delta": raw_delta,
                "posterior_count_delta": posterior_delta,
                "posterior_count_closer": closer,
                "tolerance": f"{tol:.9g}",
                "conclusion": conclusion,
                "historical_ref": spec["historical_ref"],
                "historical_source": source,
            })

    writer = csv.DictWriter(sys.stdout, fieldnames=out_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(out_rows)

    headline = headline_from_rows(out_rows)
    (full_root / "dscale_semantics_headline.txt").write_text(headline + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
