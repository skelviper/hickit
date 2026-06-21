#!/usr/bin/env python3
"""Audit whether raw-posterior-only summaries can select a high-trans P9016 basin."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


DEFAULT_BASELINE_CONFIG = "p9016_multiseed_td1_seed17_noise0"
STATE_COLUMNS = ("p00", "p01", "p10", "p11")
CLASSES = ("all", "cis", "trans")

SUMMARY_METRICS = [
    "n_bpair",
    "n_raw",
    "mean_entropy",
    "mean_pU",
    "mean_pmax",
    "mean_margin",
    "mean_psame_raw",
    "mean_pcross_raw",
    "mean_top1_same",
    "mean_top1_cross",
    "frac_top1_same",
    "frac_top1_cross",
    "frac_pmax_ge_0p5",
    "frac_pmax_ge_0p75",
    "frac_margin_ge_0p1",
    "frac_margin_ge_0p25",
    "frac_low_entropy_lt_0p5",
    "frac_low_pU_lt_0p5",
]

OUTPUT_FIELDS = [
    "metric",
    "selection_direction",
    "metric_source",
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


class Accumulator:
    def __init__(self) -> None:
        self.n_bpair = 0
        self.n_raw = 0.0
        self.sum_entropy = 0.0
        self.sum_pu = 0.0
        self.sum_pmax = 0.0
        self.sum_margin = 0.0
        self.sum_psame = 0.0
        self.sum_pcross = 0.0
        self.sum_top1_same = 0.0
        self.sum_top1_cross = 0.0
        self.count_top1_same = 0.0
        self.count_top1_cross = 0.0
        self.count_pmax_ge_0p5 = 0.0
        self.count_pmax_ge_0p75 = 0.0
        self.count_margin_ge_0p1 = 0.0
        self.count_margin_ge_0p25 = 0.0
        self.count_entropy_lt_0p5 = 0.0
        self.count_pu_lt_0p5 = 0.0

    def add(self, row: dict[str, str]) -> None:
        n_raw = as_float(row.get("n_raw")) or 0.0
        if n_raw <= 0.0:
            n_raw = 1.0
        p = [as_float(row.get(col)) or 0.0 for col in STATE_COLUMNS]
        pmax = as_float(row.get("pmax"))
        if pmax is None:
            pmax = max(p)
        margin = as_float(row.get("margin")) or 0.0
        entropy = as_float(row.get("entropy")) or 0.0
        pu = as_float(row.get("pU")) or 0.0
        psame = as_float(row.get("psame_raw"))
        if psame is None:
            psame = p[0] + p[3]
        pcross = as_float(row.get("pcross_raw"))
        if pcross is None:
            pcross = p[1] + p[2]
        top_state = max(range(4), key=lambda idx: p[idx])
        top1_same = 1.0 if top_state in (0, 3) else 0.0
        top1_cross = 1.0 - top1_same

        self.n_bpair += 1
        self.n_raw += n_raw
        self.sum_entropy += n_raw * entropy
        self.sum_pu += n_raw * pu
        self.sum_pmax += n_raw * pmax
        self.sum_margin += n_raw * margin
        self.sum_psame += n_raw * psame
        self.sum_pcross += n_raw * pcross
        self.sum_top1_same += n_raw * top1_same
        self.sum_top1_cross += n_raw * top1_cross
        self.count_top1_same += top1_same
        self.count_top1_cross += top1_cross
        self.count_pmax_ge_0p5 += n_raw if pmax >= 0.5 else 0.0
        self.count_pmax_ge_0p75 += n_raw if pmax >= 0.75 else 0.0
        self.count_margin_ge_0p1 += n_raw if margin >= 0.1 else 0.0
        self.count_margin_ge_0p25 += n_raw if margin >= 0.25 else 0.0
        self.count_entropy_lt_0p5 += n_raw if entropy < 0.5 else 0.0
        self.count_pu_lt_0p5 += n_raw if pu < 0.5 else 0.0

    def finish(self) -> dict[str, float]:
        if self.n_raw <= 0.0:
            return {name: math.nan for name in SUMMARY_METRICS}
        return {
            "n_bpair": float(self.n_bpair),
            "n_raw": self.n_raw,
            "mean_entropy": self.sum_entropy / self.n_raw,
            "mean_pU": self.sum_pu / self.n_raw,
            "mean_pmax": self.sum_pmax / self.n_raw,
            "mean_margin": self.sum_margin / self.n_raw,
            "mean_psame_raw": self.sum_psame / self.n_raw,
            "mean_pcross_raw": self.sum_pcross / self.n_raw,
            "mean_top1_same": self.sum_top1_same / self.n_raw,
            "mean_top1_cross": self.sum_top1_cross / self.n_raw,
            "frac_top1_same": self.count_top1_same / max(1.0, float(self.n_bpair)),
            "frac_top1_cross": self.count_top1_cross / max(1.0, float(self.n_bpair)),
            "frac_pmax_ge_0p5": self.count_pmax_ge_0p5 / self.n_raw,
            "frac_pmax_ge_0p75": self.count_pmax_ge_0p75 / self.n_raw,
            "frac_margin_ge_0p1": self.count_margin_ge_0p1 / self.n_raw,
            "frac_margin_ge_0p25": self.count_margin_ge_0p25 / self.n_raw,
            "frac_low_entropy_lt_0p5": self.count_entropy_lt_0p5 / self.n_raw,
            "frac_low_pU_lt_0p5": self.count_pu_lt_0p5 / self.n_raw,
        }


def posterior_path_for_row(row: dict[str, str]) -> Path | None:
    output_dir = row.get("output_dir")
    if not output_dir or output_dir == "NA":
        return None
    path = Path(output_dir) / "p9016_full.bpair_posterior.tsv"
    return path if path.exists() else None


def summarize_posterior(path: Path) -> dict[str, float]:
    acc = {name: Accumulator() for name in CLASSES}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cls = row.get("contact_class", "all")
            acc["all"].add(row)
            if cls in ("cis", "trans"):
                acc[cls].add(row)
    out: dict[str, float] = {}
    for cls in CLASSES:
        values = acc[cls].finish()
        for name, value in values.items():
            out[f"rawpost_{cls}_{name}"] = value
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--baseline-config", default=DEFAULT_BASELINE_CONFIG)
    parser.add_argument("--target-delta", type=float, default=0.1)
    parser.add_argument("--raw-summary-out", type=Path)
    args = parser.parse_args()

    summary_path = args.run_root / "summary.tsv"
    rows = read_table(summary_path)
    ok_rows = [row for row in rows if row.get("status", "OK") == "OK"]
    if not ok_rows:
        print(f"no OK rows found in {summary_path}", file=sys.stderr)
        return 1

    raw_rows: list[dict[str, str]] = []
    for row in ok_rows:
        config = row.get("config_name", "NA")
        path = posterior_path_for_row(row)
        raw_values: dict[str, float] = {}
        status = "OK"
        if path is None:
            status = "MISSING_POSTERIOR"
        else:
            raw_values = summarize_posterior(path)
        raw_row = {
            "config_name": config,
            "status": status,
            "posterior_path": "NA" if path is None else str(path),
        }
        for key in sorted(raw_values):
            raw_row[key] = fmt(raw_values[key])
        raw_rows.append(raw_row)

    raw_columns = ["config_name", "status", "posterior_path"] + [
        f"rawpost_{cls}_{metric}" for cls in CLASSES for metric in SUMMARY_METRICS
    ]
    if args.raw_summary_out is not None:
        with args.raw_summary_out.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=raw_columns, lineterminator="\n")
            writer.writeheader()
            for row in raw_rows:
                writer.writerow({col: row.get(col, "NA") for col in raw_columns})

    raw_by_config = {row["config_name"]: row for row in raw_rows if row["status"] == "OK"}
    evaluated_rows = [
        row
        for row in ok_rows
        if row.get("config_name") in raw_by_config and as_float(row.get("model_top1_accuracy_genome_trans")) is not None
    ]
    if not evaluated_rows:
        print("no evaluated rows with raw-posterior summaries found", file=sys.stderr)
        return 1
    baseline = next((row for row in evaluated_rows if row.get("config_name") == args.baseline_config), None)
    if baseline is None:
        print(f"baseline config not found: {args.baseline_config}", file=sys.stderr)
        return 1
    baseline_trans = as_float(baseline.get("model_top1_accuracy_genome_trans"))
    if baseline_trans is None:
        print(f"baseline trans metric missing: {args.baseline_config}", file=sys.stderr)
        return 1
    target_trans = baseline_trans + args.target_delta
    best = max(evaluated_rows, key=lambda row: as_float(row["model_top1_accuracy_genome_trans"]) or -math.inf)
    best_trans = as_float(best.get("model_top1_accuracy_genome_trans"))
    ranked_by_trans = sorted(
        evaluated_rows, key=lambda row: as_float(row["model_top1_accuracy_genome_trans"]) or -math.inf, reverse=True
    )
    trans_rank = {row["config_name"]: idx + 1 for idx, row in enumerate(ranked_by_trans)}
    trans_values = [as_float(row.get("model_top1_accuracy_genome_trans")) for row in evaluated_rows]

    metrics = raw_columns[3:]
    writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=OUTPUT_FIELDS, lineterminator="\n")
    writer.writeheader()
    for metric in metrics:
        values = [as_float(raw_by_config[row["config_name"]].get(metric)) for row in evaluated_rows]
        usable = [(row, value) for row, value in zip(evaluated_rows, values) if value is not None]
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
                    "metric_source": "final_raw_posterior",
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
                    "baseline_config": args.baseline_config,
                    "baseline_trans": fmt(baseline_trans),
                    "target_trans": fmt(target_trans),
                }
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
