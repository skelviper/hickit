#!/usr/bin/env python3
"""Audit whether raw-heldout diagnostics can select a high-trans P9016 basin."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


DEFAULT_BASELINE_CONFIG = "p9016_heldout_td1_seed17_noise0"

HELDOUT_METRICS = [
    ("heldout_mean_expected_energy", "heldout_raw_fit", "posterior-expected FDG-flat energy on heldout raw bpairs"),
    ("heldout_mean_min_energy", "heldout_raw_fit", "best-state FDG-flat energy on heldout raw bpairs"),
    ("heldout_mean_entropy", "heldout_raw_fit", "heldout posterior entropy"),
    ("heldout_mean_pU", "heldout_raw_fit", "heldout unknown posterior mass proxy"),
    (
        "heldout_mean_best_normalized_distance",
        "heldout_raw_fit",
        "distance/target distance for heldout best state",
    ),
    ("heldout_short_distance_frac", "heldout_raw_fit", "fraction of heldout bpairs with best normalized distance < 1"),
]

SUPPORTING_BLIND_METRICS = [
    ("final_mean_entropy", "training_loop", "training posterior entropy"),
    ("final_mean_pU", "training_loop", "training unknown posterior mass proxy"),
    ("final_contact_energy", "force", "final contact energy"),
    ("final_repulsion_energy", "force", "final repulsion energy"),
    ("final_backbone_energy", "force", "final backbone energy"),
    ("final_min_sep", "geometry", "minimum homolog copy separation"),
    ("sep_p05", "geometry", "5th percentile homolog copy separation"),
    ("sep_median", "geometry", "median homolog copy separation"),
    ("final_mean_sep", "geometry", "mean homolog copy separation"),
    ("copytrack_vector_cos_p05", "copytrack_geometry", "5th percentile adjacent homolog-vector cosine"),
    ("copytrack_vector_cos_median", "copytrack_geometry", "median adjacent homolog-vector cosine"),
    ("copytrack_frac_cos_lt_0", "copytrack_geometry", "fraction adjacent homolog-vector cosine < 0"),
    (
        "copytrack_frac_projection_sign_switch",
        "copytrack_geometry",
        "projection sign switch fraction",
    ),
]

OUTPUT_FIELDS = [
    "metric",
    "selection_direction",
    "metric_source",
    "metric_note",
    "n_ranked",
    "spearman_with_eval_trans",
    "selected_config",
    "selected_metric_value",
    "selected_trans",
    "selected_delta_vs_baseline",
    "selected_gap_to_target",
    "selected_is_eval_best",
    "selected_rank_by_trans",
    "target_met",
    "eval_best_config",
    "eval_best_trans",
    "baseline_config",
    "baseline_trans",
    "target_trans",
]


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def as_float(value: str | None) -> float | None:
    if value in (None, "", "NA"):
        return None
    try:
        x = float(value)
    except ValueError:
        return None
    if not math.isfinite(x):
        return None
    return x


def fmt(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.9g}"


def rankdata(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda idx: values[idx])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        avg_rank = (start + 1 + end) / 2.0
        for offset in range(start, end):
            ranks[order[offset]] = avg_rank
        start = end
    return ranks


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx <= 0.0 or vy <= 0.0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def spearman(xs: list[float | None], ys: list[float | None]) -> float | None:
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(pairs) < 3:
        return None
    xr = rankdata([x for x, _ in pairs])
    yr = rankdata([y for _, y in pairs])
    return pearson(xr, yr)


def find_summary(path: Path) -> Path:
    if path.is_dir():
        path = path / "summary.tsv"
    if not path.exists():
        raise FileNotFoundError(f"summary TSV not found: {path}")
    return path


def write_selector_rows(
    writer: csv.DictWriter[str],
    rows: list[dict[str, str]],
    metric_specs: list[tuple[str, str, str]],
    baseline_config: str,
    target_delta: float,
) -> int:
    baseline = next((row for row in rows if row.get("config_name") == baseline_config), None)
    if baseline is None:
        print(f"baseline config not found: {baseline_config}", file=sys.stderr)
        return 1
    baseline_trans = as_float(baseline.get("model_top1_accuracy_genome_trans"))
    if baseline_trans is None:
        print(f"baseline trans metric missing: {baseline_config}", file=sys.stderr)
        return 1
    target_trans = baseline_trans + target_delta
    evaluated_rows = [row for row in rows if as_float(row.get("model_top1_accuracy_genome_trans")) is not None]
    if not evaluated_rows:
        print("no evaluated rows found", file=sys.stderr)
        return 1
    best = max(evaluated_rows, key=lambda row: as_float(row["model_top1_accuracy_genome_trans"]) or -math.inf)
    best_trans = as_float(best.get("model_top1_accuracy_genome_trans"))
    ranked_by_trans = sorted(
        evaluated_rows, key=lambda row: as_float(row["model_top1_accuracy_genome_trans"]) or -math.inf, reverse=True
    )
    trans_rank = {row["config_name"]: idx + 1 for idx, row in enumerate(ranked_by_trans)}
    trans_values = [as_float(row.get("model_top1_accuracy_genome_trans")) for row in rows]

    for metric, source, note in metric_specs:
        values = [as_float(row.get(metric)) for row in rows]
        usable = [(row, value) for row, value in zip(rows, values) if value is not None]
        if not usable:
            continue
        rho = spearman(values, trans_values)
        for direction in ("low", "high"):
            selected, selected_value = (
                min(usable, key=lambda item: item[1]) if direction == "low" else max(usable, key=lambda item: item[1])
            )
            selected_trans = as_float(selected.get("model_top1_accuracy_genome_trans"))
            writer.writerow(
                {
                    "metric": metric,
                    "selection_direction": direction,
                    "metric_source": source,
                    "metric_note": note,
                    "n_ranked": len(usable),
                    "spearman_with_eval_trans": fmt(rho),
                    "selected_config": selected.get("config_name", "NA"),
                    "selected_metric_value": fmt(selected_value),
                    "selected_trans": fmt(selected_trans),
                    "selected_delta_vs_baseline": fmt(
                        None if selected_trans is None else selected_trans - baseline_trans
                    ),
                    "selected_gap_to_target": fmt(None if selected_trans is None else target_trans - selected_trans),
                    "selected_is_eval_best": "1" if selected.get("config_name") == best.get("config_name") else "0",
                    "selected_rank_by_trans": trans_rank.get(selected.get("config_name"), "NA"),
                    "target_met": "1" if selected_trans is not None and selected_trans >= target_trans else "0",
                    "eval_best_config": best.get("config_name", "NA"),
                    "eval_best_trans": fmt(best_trans),
                    "baseline_config": baseline_config,
                    "baseline_trans": fmt(baseline_trans),
                    "target_trans": fmt(target_trans),
                }
            )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root_or_summary", type=Path)
    parser.add_argument("--baseline-config", default=DEFAULT_BASELINE_CONFIG)
    parser.add_argument("--target-delta", type=float, default=0.1)
    args = parser.parse_args()

    summary_path = find_summary(args.run_root_or_summary)
    rows = [
        row
        for row in read_table(summary_path)
        if row.get("status", "OK") == "OK" and as_float(row.get("heldout_n_raw_heldout")) not in (None, 0.0)
    ]
    if not rows:
        print(f"no heldout-enabled evaluated OK rows found in {summary_path}", file=sys.stderr)
        return 1

    writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=OUTPUT_FIELDS, lineterminator="\n")
    writer.writeheader()
    return write_selector_rows(
        writer,
        rows,
        HELDOUT_METRICS + SUPPORTING_BLIND_METRICS,
        args.baseline_config,
        args.target_delta,
    )


if __name__ == "__main__":
    raise SystemExit(main())
