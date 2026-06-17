#!/usr/bin/env python3
"""Compute homolog-vector continuity diagnostics for a P9016 reconstruction."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


COLUMNS = [
    "chrom",
    "n_beads",
    "n_adjacent_pairs",
    "sep_min",
    "sep_p01",
    "sep_p05",
    "sep_p10",
    "sep_p25",
    "sep_median",
    "sep_mean",
    "sep_p75",
    "sep_p90",
    "sep_p95",
    "sep_p99",
    "sep_max",
    "vector_cos_min",
    "vector_cos_p01",
    "vector_cos_p05",
    "vector_cos_p10",
    "vector_cos_median",
    "vector_cos_mean",
    "vector_cos_p90",
    "vector_cos_max",
    "n_cos_lt_0",
    "frac_cos_lt_0",
    "n_cos_lt_neg0p5",
    "frac_cos_lt_neg0p5",
    "n_cos_lt_neg0p9",
    "frac_cos_lt_neg0p9",
    "mean_delta_v_norm",
    "median_delta_v_norm",
    "mean_delta_v_over_sep",
    "median_delta_v_over_sep",
    "n_projection_sign_switch",
    "frac_projection_sign_switch",
    "min_projection_run_length",
    "median_projection_run_length",
    "max_projection_run_length",
]


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    tail = chrom[3:] if chrom.startswith("chr") else chrom
    if tail.isdigit():
        return (int(tail), "")
    return (10_000, tail)


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.9g}"
    return str(value)


def percentiles(values: list[float], qs: list[tuple[str, float]]) -> dict[str, float | None]:
    if not values:
        return {name: None for name, _ in qs}
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {name: None for name, _ in qs}
    return {name: float(np.percentile(arr, q)) for name, q in qs}


def read_coords(path: Path) -> dict[str, dict[int, dict[int, np.ndarray]]]:
    coords: dict[str, dict[int, dict[int, np.ndarray]]] = defaultdict(lambda: defaultdict(dict))
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"chr", "start", "copy", "x", "y", "z"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} missing required columns {sorted(required)}")
        for row in reader:
            copy = int(row["copy"])
            if copy not in (0, 1):
                continue
            coords[row["chr"]][int(row["start"])][copy] = np.array(
                [float(row["x"]), float(row["y"]), float(row["z"])],
                dtype=float,
            )
    return coords


def run_lengths(signs: list[int]) -> list[int]:
    if not signs:
        return []
    lengths: list[int] = []
    cur = signs[0]
    n = 1
    for sign in signs[1:]:
        if sign == cur:
            n += 1
        else:
            lengths.append(n)
            cur = sign
            n = 1
    lengths.append(n)
    return lengths


def metric_row(
    chrom: str,
    seps: list[float],
    adj_cos: list[float],
    delta_norms: list[float],
    delta_over_sep: list[float],
    projection_switches: int | None,
    projection_denominator: int | None,
    projection_run_lengths: list[int] | None,
) -> dict[str, object]:
    row: dict[str, object] = {
        "chrom": chrom,
        "n_beads": len(seps),
        "n_adjacent_pairs": len(delta_norms),
    }
    row.update(percentiles(seps, [
        ("sep_min", 0), ("sep_p01", 1), ("sep_p05", 5), ("sep_p10", 10),
        ("sep_p25", 25), ("sep_median", 50), ("sep_p75", 75), ("sep_p90", 90),
        ("sep_p95", 95), ("sep_p99", 99), ("sep_max", 100),
    ]))
    row["sep_mean"] = float(np.mean(seps)) if seps else None

    row.update(percentiles(adj_cos, [
        ("vector_cos_min", 0), ("vector_cos_p01", 1), ("vector_cos_p05", 5),
        ("vector_cos_p10", 10), ("vector_cos_median", 50),
        ("vector_cos_p90", 90), ("vector_cos_max", 100),
    ]))
    row["vector_cos_mean"] = float(np.mean(adj_cos)) if adj_cos else None
    n_cos = len(adj_cos)
    for label, threshold in [("0", 0.0), ("neg0p5", -0.5), ("neg0p9", -0.9)]:
        n = sum(1 for x in adj_cos if x < threshold)
        row[f"n_cos_lt_{label}"] = n
        row[f"frac_cos_lt_{label}"] = (n / n_cos) if n_cos else None

    row["mean_delta_v_norm"] = float(np.mean(delta_norms)) if delta_norms else None
    row["median_delta_v_norm"] = float(np.median(delta_norms)) if delta_norms else None
    row["mean_delta_v_over_sep"] = float(np.mean(delta_over_sep)) if delta_over_sep else None
    row["median_delta_v_over_sep"] = float(np.median(delta_over_sep)) if delta_over_sep else None

    if projection_switches is None or projection_denominator is None or projection_run_lengths is None:
        row["n_projection_sign_switch"] = None
        row["frac_projection_sign_switch"] = None
        row["min_projection_run_length"] = None
        row["median_projection_run_length"] = None
        row["max_projection_run_length"] = None
    else:
        row["n_projection_sign_switch"] = projection_switches
        row["frac_projection_sign_switch"] = (
            projection_switches / projection_denominator if projection_denominator > 0 else 0.0
        )
        row["min_projection_run_length"] = min(projection_run_lengths) if projection_run_lengths else None
        row["median_projection_run_length"] = (
            float(np.median(projection_run_lengths)) if projection_run_lengths else None
        )
        row["max_projection_run_length"] = max(projection_run_lengths) if projection_run_lengths else None
    return row


def raw_metrics(items: list[tuple[int, np.ndarray]]) -> dict[str, object]:
    items = sorted(items, key=lambda x: x[0])
    vectors = [v for _, v in items if np.all(np.isfinite(v))]
    seps = [float(np.linalg.norm(v)) for v in vectors]
    adj_cos: list[float] = []
    delta_norms: list[float] = []
    delta_over_sep: list[float] = []
    for a, b in zip(vectors, vectors[1:]):
        na = float(np.linalg.norm(a))
        nb = float(np.linalg.norm(b))
        if na > 0.0 and nb > 0.0:
            adj_cos.append(float(np.dot(a, b) / (na * nb)))
        delta = float(np.linalg.norm(b - a))
        delta_norms.append(delta)
        denom = 0.5 * (na + nb)
        if denom > 0.0:
            delta_over_sep.append(delta / denom)

    vbar = np.mean(np.vstack(vectors), axis=0) if vectors else None
    if vbar is None or float(np.linalg.norm(vbar)) < 1e-9:
        projection_switches = None
        projection_denominator = None
        projection_run_lengths = None
    else:
        signs: list[int] = []
        for v in vectors:
            dot = float(np.dot(v, vbar))
            signs.append(1 if dot >= 0.0 else -1)
        projection_switches = sum(1 for a, b in zip(signs, signs[1:]) if a != b)
        projection_denominator = max(0, len(signs) - 1)
        projection_run_lengths = run_lengths(signs)
    return {
        "seps": seps,
        "adj_cos": adj_cos,
        "delta_norms": delta_norms,
        "delta_over_sep": delta_over_sep,
        "projection_switches": projection_switches,
        "projection_denominator": projection_denominator,
        "projection_run_lengths": projection_run_lengths,
    }


def summarize_vectors(chrom: str, items: list[tuple[int, np.ndarray]]) -> dict[str, object]:
    metrics = raw_metrics(items)
    return metric_row(
        chrom,
        metrics["seps"],  # type: ignore[arg-type]
        metrics["adj_cos"],  # type: ignore[arg-type]
        metrics["delta_norms"],  # type: ignore[arg-type]
        metrics["delta_over_sep"],  # type: ignore[arg-type]
        metrics["projection_switches"],  # type: ignore[arg-type]
        metrics["projection_denominator"],  # type: ignore[arg-type]
        metrics["projection_run_lengths"],  # type: ignore[arg-type]
    )


def summarize_all(chrom_metrics: list[dict[str, object]]) -> dict[str, object]:
    seps: list[float] = []
    adj_cos: list[float] = []
    delta_norms: list[float] = []
    delta_over_sep: list[float] = []
    projection_switches = 0
    projection_denominator = 0
    projection_run_lengths: list[int] = []
    projection_defined = False
    for metrics in chrom_metrics:
        seps.extend(metrics["seps"])  # type: ignore[arg-type]
        adj_cos.extend(metrics["adj_cos"])  # type: ignore[arg-type]
        delta_norms.extend(metrics["delta_norms"])  # type: ignore[arg-type]
        delta_over_sep.extend(metrics["delta_over_sep"])  # type: ignore[arg-type]
        if metrics["projection_switches"] is not None:
            projection_defined = True
            projection_switches += int(metrics["projection_switches"])
            projection_denominator += int(metrics["projection_denominator"])
            projection_run_lengths.extend(metrics["projection_run_lengths"])  # type: ignore[arg-type]
    return metric_row(
        "ALL",
        seps,
        adj_cos,
        delta_norms,
        delta_over_sep,
        projection_switches if projection_defined else None,
        projection_denominator if projection_defined else None,
        projection_run_lengths if projection_defined else None,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config-output-dir", type=Path, required=True)
    parser.add_argument("--eval-output-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    coords_path = args.config_output_dir / "p9016_full.coords.tsv"
    coords = read_coords(coords_path)
    rows: list[dict[str, object]] = []
    all_metrics: list[dict[str, object]] = []
    for chrom in sorted(coords, key=chrom_sort_key):
        items: list[tuple[int, np.ndarray]] = []
        for start in sorted(coords[chrom]):
            copies = coords[chrom][start]
            if 0 in copies and 1 in copies:
                v = copies[1] - copies[0]
                items.append((start, v))
        metrics = raw_metrics(items)
        all_metrics.append(metrics)
        rows.append(metric_row(
            chrom,
            metrics["seps"],  # type: ignore[arg-type]
            metrics["adj_cos"],  # type: ignore[arg-type]
            metrics["delta_norms"],  # type: ignore[arg-type]
            metrics["delta_over_sep"],  # type: ignore[arg-type]
            metrics["projection_switches"],  # type: ignore[arg-type]
            metrics["projection_denominator"],  # type: ignore[arg-type]
            metrics["projection_run_lengths"],  # type: ignore[arg-type]
        ))
    rows.append(summarize_all(all_metrics))

    args.eval_output_dir.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in COLUMNS})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
