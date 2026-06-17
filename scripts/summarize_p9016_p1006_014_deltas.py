#!/usr/bin/env python3
"""Write gamma, separation, and sample deltas for experiment 014."""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path


METRICS = {
    "top1_all": "model_top1_accuracy_genome_all",
    "top1_cis": "model_top1_accuracy_genome_cis",
    "top1_trans": "model_top1_accuracy_genome_trans",
    "same_cross_cis": "model_same_cross_accuracy_genome_cis",
    "cis_spearman": "mean_per_chrom_cis_distance_spearman",
    "entropy": "final_mean_entropy",
    "pU": "final_mean_pU",
    "min_sep": "final_min_sep",
    "sep_p05": "sep_p05",
    "mean_sep": "final_mean_sep",
    "copytrack_frac_cos_lt_0": "copytrack_frac_cos_lt_0",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f(row: dict[str, str], key: str) -> float | None:
    try:
        value = float(row.get(key, "NA"))
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def fmt(value: float | None) -> str:
    return "NA" if value is None else f"{value:.9g}"


def delta(row: dict[str, str], base: dict[str, str], key: str) -> str:
    a = f(row, key)
    b = f(base, key)
    if a is None or b is None:
        return "NA"
    return fmt(a - b)


def gamma_value(row: dict[str, str]) -> str:
    value = row.get("d_scale_posterior_gamma") or row.get("dscale_posterior_gamma") or "NA"
    try:
        return f"{float(value):.9g}"
    except ValueError:
        return value


def setting_key(row: dict[str, str]) -> str:
    return "gamma={}:msep={}:lsep={}".format(
        gamma_value(row),
        row.get("min_sep_unit", "NA"),
        row.get("lambda_sep", "NA"),
    )


def is_sep_off(row: dict[str, str]) -> bool:
    return f(row, "min_sep_unit") == 0.0 and f(row, "lambda_sep") == 0.0


def gamma_float(row: dict[str, str]) -> float:
    value = f(row, "d_scale_posterior_gamma")
    return value if value is not None else math.nan


def config_gamma_rank(rows: list[dict[str, str]], metric_key: str) -> dict[str, int]:
    ordered = sorted(
        [r for r in rows if f(r, metric_key) is not None],
        key=lambda r: f(r, metric_key),
        reverse=True,
    )
    return {r["config_name"]: i + 1 for i, r in enumerate(ordered)}


def write_gamma_delta(root: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "sample",
        "gamma",
        "config_name",
        "top1_all",
        "top1_cis",
        "top1_trans",
        "same_cross_cis",
        "cis_spearman",
        "entropy",
        "pU",
        "min_sep",
        "sep_p05",
        "mean_sep",
        "copytrack_frac_cos_lt_0",
        "delta_top1_all_vs_gamma1",
        "delta_top1_cis_vs_gamma1",
        "delta_top1_trans_vs_gamma1",
        "delta_same_cross_cis_vs_gamma1",
        "delta_spearman_vs_gamma1",
        "delta_entropy_vs_gamma1",
        "delta_pU_vs_gamma1",
        "delta_min_sep_vs_gamma1",
        "gamma_rank_by_top1_all",
        "gamma_rank_by_cis_spearman",
    ]
    out_rows: list[dict[str, str]] = []
    for sample in sorted({r.get("sample", "") for r in rows}):
        sample_rows = [r for r in rows if r.get("sample") == sample and is_sep_off(r)]
        base = next((r for r in sample_rows if abs(gamma_float(r) - 1.0) < 1e-8), None)
        top_rank = config_gamma_rank(sample_rows, METRICS["top1_all"])
        sp_rank = config_gamma_rank(sample_rows, METRICS["cis_spearman"])
        for row in sorted(sample_rows, key=gamma_float):
            out = {"sample": sample, "gamma": gamma_value(row), "config_name": row["config_name"]}
            for out_key, in_key in METRICS.items():
                out[out_key] = row.get(in_key, "NA")
            if base is None:
                for key in fields:
                    out.setdefault(key, "NA")
            else:
                out["delta_top1_all_vs_gamma1"] = delta(row, base, METRICS["top1_all"])
                out["delta_top1_cis_vs_gamma1"] = delta(row, base, METRICS["top1_cis"])
                out["delta_top1_trans_vs_gamma1"] = delta(row, base, METRICS["top1_trans"])
                out["delta_same_cross_cis_vs_gamma1"] = delta(row, base, METRICS["same_cross_cis"])
                out["delta_spearman_vs_gamma1"] = delta(row, base, METRICS["cis_spearman"])
                out["delta_entropy_vs_gamma1"] = delta(row, base, METRICS["entropy"])
                out["delta_pU_vs_gamma1"] = delta(row, base, METRICS["pU"])
                out["delta_min_sep_vs_gamma1"] = delta(row, base, METRICS["min_sep"])
            out["gamma_rank_by_top1_all"] = str(top_rank.get(row["config_name"], "NA"))
            out["gamma_rank_by_cis_spearman"] = str(sp_rank.get(row["config_name"], "NA"))
            out_rows.append(out)
    write_table(root / "gamma_delta.tsv", fields, out_rows)


def write_sep_delta(root: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "sample",
        "gamma",
        "sep_config_name",
        "base_config_name",
        "min_sep_unit",
        "lambda_sep",
        "top1_all",
        "top1_cis",
        "top1_trans",
        "same_cross_cis",
        "cis_spearman",
        "entropy",
        "pU",
        "min_sep",
        "sep_p05",
        "frac_sep_below_min",
        "delta_top1_all",
        "delta_top1_cis",
        "delta_top1_trans",
        "delta_same_cross_cis",
        "delta_spearman",
        "delta_entropy",
        "delta_pU",
        "delta_min_sep",
        "delta_sep_p05",
        "delta_frac_sep_below_min",
    ]
    base_by_sample_gamma = {
        (r.get("sample", ""), gamma_value(r)): r
        for r in rows
        if is_sep_off(r)
    }
    out_rows = []
    for row in rows:
        if is_sep_off(row):
            continue
        gamma = gamma_value(row)
        base = base_by_sample_gamma.get((row.get("sample", ""), gamma))
        out = {
            "sample": row.get("sample", "NA"),
            "gamma": gamma,
            "sep_config_name": row.get("config_name", "NA"),
            "base_config_name": base.get("config_name", "NA") if base else "NA",
            "min_sep_unit": row.get("min_sep_unit", "NA"),
            "lambda_sep": row.get("lambda_sep", "NA"),
            "frac_sep_below_min": row.get("frac_sep_below_min", "NA"),
        }
        for out_key, in_key in METRICS.items():
            if out_key == "copytrack_frac_cos_lt_0":
                continue
            out[out_key] = row.get(in_key, "NA")
        if base:
            out["delta_top1_all"] = delta(row, base, METRICS["top1_all"])
            out["delta_top1_cis"] = delta(row, base, METRICS["top1_cis"])
            out["delta_top1_trans"] = delta(row, base, METRICS["top1_trans"])
            out["delta_same_cross_cis"] = delta(row, base, METRICS["same_cross_cis"])
            out["delta_spearman"] = delta(row, base, METRICS["cis_spearman"])
            out["delta_entropy"] = delta(row, base, METRICS["entropy"])
            out["delta_pU"] = delta(row, base, METRICS["pU"])
            out["delta_min_sep"] = delta(row, base, METRICS["min_sep"])
            out["delta_sep_p05"] = delta(row, base, METRICS["sep_p05"])
            out["delta_frac_sep_below_min"] = delta(row, base, "frac_sep_below_min")
        for key in fields:
            out.setdefault(key, "NA")
        out_rows.append(out)
    write_table(root / "sep_delta.tsv", fields, out_rows)


def write_sample_delta(root: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "setting_key",
        "p9016_config",
        "p1006_config",
        "p9016_top1_all",
        "p1006_top1_all",
        "delta_top1_all_p1006_minus_p9016",
        "p9016_top1_cis",
        "p1006_top1_cis",
        "delta_top1_cis",
        "p9016_top1_trans",
        "p1006_top1_trans",
        "delta_top1_trans",
        "p9016_same_cross_cis",
        "p1006_same_cross_cis",
        "delta_same_cross_cis",
        "p9016_cis_spearman",
        "p1006_cis_spearman",
        "delta_cis_spearman",
        "p9016_entropy",
        "p1006_entropy",
        "delta_entropy",
        "p9016_min_sep",
        "p1006_min_sep",
        "delta_min_sep",
    ]
    by_setting_sample = {(setting_key(r), r.get("sample", "")): r for r in rows}
    out_rows = []
    for key in sorted({setting_key(r) for r in rows}, key=natural_setting_key):
        p9016 = by_setting_sample.get((key, "P9016"))
        p1006 = by_setting_sample.get((key, "P1006"))
        if not p9016 or not p1006:
            continue
        out = {
            "setting_key": key,
            "p9016_config": p9016.get("config_name", "NA"),
            "p1006_config": p1006.get("config_name", "NA"),
        }
        metric_pairs = [
            ("top1_all", METRICS["top1_all"]),
            ("top1_cis", METRICS["top1_cis"]),
            ("top1_trans", METRICS["top1_trans"]),
            ("same_cross_cis", METRICS["same_cross_cis"]),
            ("cis_spearman", METRICS["cis_spearman"]),
            ("entropy", METRICS["entropy"]),
            ("min_sep", METRICS["min_sep"]),
        ]
        for label, metric in metric_pairs:
            out[f"p9016_{label}"] = p9016.get(metric, "NA")
            out[f"p1006_{label}"] = p1006.get(metric, "NA")
            delta_key = f"delta_{label}"
            if label == "top1_all":
                delta_key = "delta_top1_all_p1006_minus_p9016"
            out[delta_key] = delta(p1006, p9016, metric)
        out_rows.append(out)
    write_table(root / "sample_delta.tsv", fields, out_rows)


def natural_setting_key(value: str) -> tuple[float, float, float, str]:
    match = re.search(r"gamma=([^:]+):msep=([^:]+):lsep=([^:]+)", value)
    if not match:
        return (math.inf, math.inf, math.inf, value)
    vals = []
    for item in match.groups():
        try:
            vals.append(float(item))
        except ValueError:
            vals.append(math.inf)
    return (vals[0], vals[1], vals[2], value)


def write_table(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fields})


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: summarize_p9016_p1006_014_deltas.py FULL_RUN_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    rows = read_rows(root / "summary.tsv")
    write_gamma_delta(root, rows)
    write_sep_delta(root, rows)
    write_sample_delta(root, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
