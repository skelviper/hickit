#!/usr/bin/env python3
"""Diagnose softall contact-force regimes from final coordinates."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "eval"))

from evaluate_p9016_baseline import (  # noqa: E402
    FDG_C_C1,
    FDG_C_C2,
    FDG_D_C1,
    FDG_D_C2,
    FDG_D_C3,
    STATE_LABELS,
    distance_for_copies,
    read_posterior,
    read_reconstruction,
)


REGIME_COLUMNS = [
    "contact_class",
    "state",
    "regime",
    "n_edges",
    "sum_k",
    "sum_abs_force",
    "sum_endpoint_force_l2",
    "sum_energy",
    "mean_r",
    "median_r",
    "mean_distance",
    "mean_d_scale",
    "mean_prob",
]

SUMMARY_COLUMNS = [
    "contact_class",
    "n_edges",
    "sum_k",
    "sum_abs_force",
    "sum_energy",
    "mean_r",
    "median_r",
    "p10_r",
    "p90_r",
    "mean_distance",
    "mean_d_scale",
    "mean_prob",
    "force_per_k",
    "frac_k_repulsive",
    "frac_k_zero_force",
    "frac_k_quadratic_attractive",
    "frac_k_tail_attractive",
    "net_force_l2",
    "sum_endpoint_force_l2",
    "cancellation_ratio",
]

BPAIR_COLUMNS = [
    "chr1",
    "start1",
    "chr2",
    "start2",
    "contact_class",
    "n_raw",
    "state",
    "prob",
    "k",
    "d_scale",
    "distance",
    "r",
    "regime",
    "edge_abs_force",
    "energy",
]


def read_kv(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header != ["key", "value"]:
            raise ValueError(f"{path} is not a key/value TSV")
        for row in reader:
            if len(row) >= 2:
                values[row[0]] = row[1]
    return values


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.9g}"
    return str(value)


def regime_for_r(r: float) -> str:
    if r < FDG_D_C1:
        return "repulsive"
    if r <= FDG_D_C2:
        return "zero_force"
    if r <= FDG_D_C3:
        return "quadratic_attractive"
    return "tail_attractive"


def contact_energy_r(r: float, k: float) -> float:
    if r < FDG_D_C1:
        t = FDG_D_C1 - r
        return k * t * t
    if r <= FDG_D_C2:
        return 0.0
    if r <= FDG_D_C3:
        t = r - FDG_D_C2
        return k * t * t
    t = r - FDG_D_C2
    return k * (FDG_C_C1 * (r - FDG_D_C3) + FDG_C_C2 / t)


def force_mag_r(r: float, k: float) -> float:
    if r < FDG_D_C1:
        t = FDG_D_C1 - r
        return 2.0 * k * t
    if r <= FDG_D_C2:
        return 0.0
    if r <= FDG_D_C3:
        t = r - FDG_D_C2
        return -2.0 * k * t
    t = r - FDG_D_C2
    return -k * (FDG_C_C1 - FDG_C_C2 / (t * t))


def state_copies(state: int) -> tuple[int, int]:
    return state >> 1, state & 1


def state_d_scale(
    base_d_scale: float,
    n_raw: int,
    prob: float,
    dscale_mode: str,
    gamma: float,
    eps_count: float,
) -> float:
    if dscale_mode in {"posterior_count", "tempered_posterior_count"}:
        effective = n_raw * (max(prob, 0.0) ** gamma)
        if not math.isfinite(effective) or effective < eps_count:
            effective = eps_count
        return effective ** (-1.0 / 3.0)
    return base_d_scale


class Accumulator:
    def __init__(self) -> None:
        self.n_edges = 0
        self.sum_k = 0.0
        self.sum_abs_force = 0.0
        self.sum_endpoint_force_l2 = 0.0
        self.sum_energy = 0.0
        self.rs: list[float] = []
        self.r_weights: list[float] = []
        self.sum_distance_k = 0.0
        self.sum_d_scale_k = 0.0
        self.sum_prob_k = 0.0

    def update(self, k: float, abs_force: float, energy: float, r: float, distance: float, d_scale: float, prob: float) -> None:
        self.n_edges += 1
        self.sum_k += k
        self.sum_abs_force += abs_force
        self.sum_endpoint_force_l2 += 2.0 * abs_force
        self.sum_energy += energy
        self.rs.append(r)
        self.r_weights.append(k)
        self.sum_distance_k += distance * k
        self.sum_d_scale_k += d_scale * k
        self.sum_prob_k += prob * k

    def row(self, **labels: object) -> dict[str, object]:
        r_arr = np.asarray(self.rs, dtype=float)
        out = {
            **labels,
            "n_edges": self.n_edges,
            "sum_k": self.sum_k,
            "sum_abs_force": self.sum_abs_force,
            "sum_endpoint_force_l2": self.sum_endpoint_force_l2,
            "sum_energy": self.sum_energy,
            "mean_r": float(np.average(r_arr, weights=np.asarray(self.r_weights))) if self.rs and self.sum_k > 0 else float("nan"),
            "median_r": float(np.median(r_arr)) if self.rs else float("nan"),
            "mean_distance": self.sum_distance_k / self.sum_k if self.sum_k > 0 else float("nan"),
            "mean_d_scale": self.sum_d_scale_k / self.sum_k if self.sum_k > 0 else float("nan"),
            "mean_prob": self.sum_prob_k / self.sum_k if self.sum_k > 0 else float("nan"),
        }
        return out


def percentile(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    return float(np.percentile(np.asarray(values, dtype=float), q))


def build_rows(
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    coords: dict[tuple[str, int, int], np.ndarray],
    manifest: dict[str, str],
    write_bpair_rows: bool,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    dscale_mode = manifest.get("dscale_mode", manifest.get("d_scale_mode", "raw_count"))
    gamma = float(manifest.get("dscale_posterior_gamma", manifest.get("d_scale_posterior_gamma", "1")))
    eps_count = float(manifest.get("d_scale_eps_count", "1e-6"))
    unit = float(manifest.get("unit", "1"))
    by_regime: dict[tuple[str, str, str], Accumulator] = defaultdict(Accumulator)
    by_class: dict[str, Accumulator] = defaultdict(Accumulator)
    regime_k_by_class: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    net_force: dict[str, dict[tuple[str, int, int], np.ndarray]] = defaultdict(lambda: defaultdict(lambda: np.zeros(3, dtype=float)))
    bpair_rows: list[dict[str, object]] = []

    for item in posterior.values():
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        start1 = int(item["start1"])
        start2 = int(item["start2"])
        cls = str(item["contact_class"])
        n_raw = int(item["n_raw"])
        base_k = float(item["base_k"])
        base_d_scale = float(item["base_d_scale"])
        p4 = np.asarray(item["p4"], dtype=float)
        for state, prob in enumerate(p4):
            prob = float(prob)
            k = base_k * prob
            if k <= 0.0:
                continue
            copy1, copy2 = state_copies(state)
            distance = distance_for_copies(coords, chrom1, start1, copy1, chrom2, start2, copy2)
            if distance is None:
                continue
            d_scale = state_d_scale(base_d_scale, n_raw, prob, dscale_mode, gamma, eps_count)
            r = (distance / unit) / d_scale
            if not math.isfinite(r):
                continue
            regime = regime_for_r(r)
            energy = contact_energy_r(r, k)
            fmag = force_mag_r(r, k)
            abs_force = abs(fmag)
            key1 = (chrom1, start1, copy1)
            key2 = (chrom2, start2, copy2)
            delta = coords[key1] - coords[key2]
            norm = float(np.linalg.norm(delta))
            direction = delta / norm if norm > 0.0 else np.array([1.0, 0.0, 0.0], dtype=float)
            force_vec = fmag * direction
            net_force[cls][key1] += force_vec
            net_force[cls][key2] -= force_vec
            by_regime[(cls, STATE_LABELS[state], regime)].update(k, abs_force, energy, r, distance, d_scale, prob)
            by_class[cls].update(k, abs_force, energy, r, distance, d_scale, prob)
            regime_k_by_class[cls][regime] += k
            if write_bpair_rows:
                bpair_rows.append({
                    "chr1": chrom1,
                    "start1": start1,
                    "chr2": chrom2,
                    "start2": start2,
                    "contact_class": cls,
                    "n_raw": n_raw,
                    "state": STATE_LABELS[state],
                    "prob": prob,
                    "k": k,
                    "d_scale": d_scale,
                    "distance": distance,
                    "r": r,
                    "regime": regime,
                    "edge_abs_force": abs_force,
                    "energy": energy,
                })

    regime_rows = [
        acc.row(contact_class=cls, state=state, regime=regime)
        for (cls, state, regime), acc in sorted(by_regime.items())
    ]
    summary_rows: list[dict[str, object]] = []
    for cls, acc in sorted(by_class.items()):
        row = acc.row(contact_class=cls)
        row["p10_r"] = percentile(acc.rs, 10)
        row["p90_r"] = percentile(acc.rs, 90)
        row["force_per_k"] = acc.sum_abs_force / acc.sum_k if acc.sum_k > 0 else float("nan")
        for regime in ("repulsive", "zero_force", "quadratic_attractive", "tail_attractive"):
            row[f"frac_k_{regime}"] = regime_k_by_class[cls][regime] / acc.sum_k if acc.sum_k > 0 else float("nan")
        net_l2 = sum(float(np.linalg.norm(v)) for v in net_force[cls].values())
        sum_l2 = acc.sum_endpoint_force_l2
        row["net_force_l2"] = net_l2
        row["sum_endpoint_force_l2"] = sum_l2
        row["cancellation_ratio"] = 1.0 - (net_l2 / sum_l2) if sum_l2 > 0 else float("nan")
        summary_rows.append(row)
    return summary_rows, regime_rows, bpair_rows


def write_rows(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: fmt(row.get(col)) for col in columns})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config-output-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--write-bpair", action="store_true")
    args = parser.parse_args()

    config_dir = args.config_output_dir
    posterior = read_posterior(config_dir / "p9016_full.bpair_posterior.tsv")
    coords = read_reconstruction(config_dir / "p9016_full.coords.tsv")
    manifest = read_kv(config_dir / "p9016_full.manifest.tsv")
    summary_rows, regime_rows, bpair_rows = build_rows(posterior, coords, manifest, args.write_bpair)

    write_rows(args.outdir / "force_regime_summary.tsv", SUMMARY_COLUMNS, summary_rows)
    write_rows(args.outdir / "force_regime_diag.tsv", REGIME_COLUMNS, regime_rows)
    if args.write_bpair:
        write_rows(args.outdir / "force_regime_bpair.tsv", BPAIR_COLUMNS, bpair_rows)
    print(f"wrote\t{args.outdir / 'force_regime_summary.tsv'}")
    print(f"wrote\t{args.outdir / 'force_regime_diag.tsv'}")
    if args.write_bpair:
        print(f"wrote\t{args.outdir / 'force_regime_bpair.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
