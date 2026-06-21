#!/usr/bin/env python3
"""Decompose trans contact failures into gauge, posterior, and geometry terms."""

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
    STATE_LABELS,
    align_p4_to_truth_gauge,
    all_copy_pair_distances,
    chrom_sort_key,
    p4_from_coords_for_bpair,
    read_contact_truth_counts,
    read_posterior,
    read_reference_3dg,
    read_reconstruction,
    swap_state,
)


SUMMARY_COLUMNS = [
    "scope",
    "n_bpair",
    "n_contacts",
    "whole_chrom_top1_accuracy",
    "chr_pair_oracle_top1_accuracy",
    "nearest_geometry_state_accuracy",
    "same_cross_accuracy",
    "pmax90_accuracy",
    "pmax90_recall",
    "chr_pair_oracle_delta_vs_whole",
    "truth_state_distance_rank_mean",
    "truth_state_distance_rank_median",
    "posterior_top1_distance_rank_mean",
    "posterior_top1_distance_rank_median",
    "posterior_top1_distance_gap_mean",
    "posterior_top1_distance_gap_median",
    "truth_state_is_nearest_fraction",
    "posterior_top1_is_nearest_fraction",
]

CHR_PAIR_COLUMNS = [
    "chr1",
    "chr2",
    "n_bpair",
    "n_contacts",
    "whole_chrom_top1_accuracy",
    "chr_pair_oracle_top1_accuracy",
    "best_relative_flip_chr2",
    "nearest_geometry_state_accuracy",
    "same_cross_accuracy",
    "pmax90_accuracy",
    "pmax90_recall",
    "truth_state_distance_rank_mean",
    "truth_state_is_nearest_fraction",
    "posterior_top1_distance_rank_mean",
    "posterior_top1_is_nearest_fraction",
    "posterior_top1_distance_gap_mean",
]

BPAIR_COLUMNS = [
    "chr1",
    "start1",
    "chr2",
    "start2",
    "n_contacts",
    "truth_top_state",
    "truth_top_count",
    "posterior_top_state",
    "posterior_pmax",
    "posterior_margin",
    "whole_chrom_pred_state",
    "whole_chrom_correct_contacts",
    "chr_pair_oracle_correct_contacts",
    "nearest_geometry_state",
    "nearest_geometry_correct_contacts",
    "same_cross_correct_contacts",
    "truth_state_distance_rank",
    "posterior_top1_distance_rank",
    "posterior_top1_distance_gap",
    "nearest_geometry_state",
    "best_relative_flip_chr2",
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


def read_selected_swaps(path: Path) -> dict[str, int]:
    swaps: dict[str, int] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"chrom", "copy_swap", "selected_for_eval"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} missing required columns {sorted(required)}")
        for row in reader:
            if row["selected_for_eval"] == "1":
                swaps[row["chrom"]] = int(row["copy_swap"])
    return swaps


def state_same_cross(state: int) -> int:
    return 1 if state in (0, 3) else 0


def weighted_median(values: list[float], weights: list[int]) -> float:
    if not values:
        return float("nan")
    pairs = sorted(zip(values, weights), key=lambda x: x[0])
    total = sum(weights)
    if total <= 0:
        return float("nan")
    cutoff = 0.5 * total
    running = 0.0
    for value, weight in pairs:
        running += weight
        if running >= cutoff:
            return float(value)
    return float(pairs[-1][0])


def format_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.9g}"
    return str(value)


def state_distance_rank(distances: list[float], state: int) -> int:
    order = sorted(range(4), key=lambda idx: (distances[idx], idx))
    return order.index(state) + 1


def second_distance_gap(distances: list[float], state: int) -> float:
    sorted_dist = sorted(distances)
    if state != int(np.argmin(np.asarray(distances, dtype=float))):
        return sorted_dist[1] - sorted_dist[0]
    return sorted_dist[1] - sorted_dist[0]


class Accumulator:
    def __init__(self) -> None:
        self.n_bpair = 0
        self.n_contacts = 0
        self.whole_correct = 0
        self.oracle_correct = 0
        self.nearest_geometry_correct = 0
        self.same_cross_correct = 0
        self.called_contacts = 0
        self.called_correct = 0
        self.truth_rank_values: list[float] = []
        self.truth_rank_weights: list[int] = []
        self.posterior_rank_values: list[float] = []
        self.posterior_rank_weights: list[int] = []
        self.gap_values: list[float] = []
        self.gap_weights: list[int] = []
        self.truth_nearest_contacts = 0
        self.posterior_nearest_contacts = 0

    def update(
        self,
        weight: int,
        whole_correct: int,
        oracle_correct: int,
        nearest_geometry_correct: int,
        same_cross_correct: int,
        pmax: float,
        truth_rank: int,
        posterior_rank: int,
        gap: float,
    ) -> None:
        self.n_bpair += 1
        self.n_contacts += weight
        self.whole_correct += whole_correct
        self.oracle_correct += oracle_correct
        self.nearest_geometry_correct += nearest_geometry_correct
        self.same_cross_correct += same_cross_correct
        if pmax >= 0.9:
            self.called_contacts += weight
            self.called_correct += whole_correct
        self.truth_rank_values.append(float(truth_rank))
        self.truth_rank_weights.append(weight)
        self.posterior_rank_values.append(float(posterior_rank))
        self.posterior_rank_weights.append(weight)
        self.gap_values.append(float(gap))
        self.gap_weights.append(weight)
        if truth_rank == 1:
            self.truth_nearest_contacts += weight
        if posterior_rank == 1:
            self.posterior_nearest_contacts += weight

    def row(self, scope: str) -> dict[str, object]:
        n = self.n_contacts
        return {
            "scope": scope,
            "n_bpair": self.n_bpair,
            "n_contacts": n,
            "whole_chrom_top1_accuracy": self.whole_correct / n if n else float("nan"),
            "chr_pair_oracle_top1_accuracy": self.oracle_correct / n if n else float("nan"),
            "nearest_geometry_state_accuracy": self.nearest_geometry_correct / n if n else float("nan"),
            "same_cross_accuracy": self.same_cross_correct / n if n else float("nan"),
            "pmax90_accuracy": self.called_correct / self.called_contacts if self.called_contacts else float("nan"),
            "pmax90_recall": self.called_correct / n if n else float("nan"),
            "chr_pair_oracle_delta_vs_whole": (self.oracle_correct - self.whole_correct) / n if n else float("nan"),
            "truth_state_distance_rank_mean": np.average(self.truth_rank_values, weights=self.truth_rank_weights)
            if self.truth_rank_values else float("nan"),
            "truth_state_distance_rank_median": weighted_median(self.truth_rank_values, self.truth_rank_weights),
            "posterior_top1_distance_rank_mean": np.average(self.posterior_rank_values, weights=self.posterior_rank_weights)
            if self.posterior_rank_values else float("nan"),
            "posterior_top1_distance_rank_median": weighted_median(self.posterior_rank_values, self.posterior_rank_weights),
            "posterior_top1_distance_gap_mean": np.average(self.gap_values, weights=self.gap_weights)
            if self.gap_values else float("nan"),
            "posterior_top1_distance_gap_median": weighted_median(self.gap_values, self.gap_weights),
            "truth_state_is_nearest_fraction": self.truth_nearest_contacts / n if n else float("nan"),
            "posterior_top1_is_nearest_fraction": self.posterior_nearest_contacts / n if n else float("nan"),
        }


def choose_chr_pair_oracle(
    records: list[dict[str, object]],
    swaps: dict[str, int],
) -> dict[tuple[str, str], int]:
    correct_by_pair_flip: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0])
    for rec in records:
        chrom1 = str(rec["chrom1"])
        chrom2 = str(rec["chrom2"])
        p4 = rec["p4"]
        counts = rec["counts"]
        swap1 = swaps.get(chrom1, 0)
        base_swap2 = swaps.get(chrom2, 0)
        for rel_flip in (0, 1):
            aligned = align_p4_to_truth_gauge(p4, swap1, base_swap2 ^ rel_flip)
            pred = int(np.argmax(aligned))
            correct_by_pair_flip[(chrom1, chrom2)][rel_flip] += int(counts[pred])
    return {
        pair: 1 if values[1] > values[0] else 0
        for pair, values in correct_by_pair_flip.items()
    }


def load_records(
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    contact_counts: dict[tuple[str, int, str, int], np.ndarray],
    coords: dict[tuple[str, int, int], np.ndarray],
    reference: dict[tuple[str, int, int], np.ndarray] | None,
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for key, counts in contact_counts.items():
        item = posterior.get(key)
        if item is None:
            continue
        if reference is not None and p4_from_coords_for_bpair(reference, item) is None:
            continue
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        if chrom1 == chrom2:
            continue
        start1 = int(item["start1"])
        start2 = int(item["start2"])
        distances = all_copy_pair_distances(coords, chrom1, start1, chrom2, start2)
        if distances is None:
            continue
        total = int(counts.sum())
        if total <= 0:
            continue
        records.append({
            "key": key,
            "chrom1": chrom1,
            "start1": start1,
            "chrom2": chrom2,
            "start2": start2,
            "counts": counts,
            "p4": np.asarray(item["p4"], dtype=float),
            "pmax": float(item["pmax"]),
            "margin": float(item.get("margin", float("nan"))),
            "distances": distances,
            "n_contacts": total,
        })
    return records


def summarize_records(
    records: list[dict[str, object]],
    swaps: dict[str, int],
    pair_relative_flips: dict[tuple[str, str], int],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    all_acc = Accumulator()
    by_pair: dict[tuple[str, str], Accumulator] = defaultdict(Accumulator)
    by_pmax: dict[str, Accumulator] = defaultdict(Accumulator)
    bpair_rows: list[dict[str, object]] = []

    for rec in records:
        chrom1 = str(rec["chrom1"])
        chrom2 = str(rec["chrom2"])
        counts = rec["counts"]
        p4 = rec["p4"]
        distances = rec["distances"]
        weight = int(rec["n_contacts"])
        swap1 = swaps.get(chrom1, 0)
        swap2 = swaps.get(chrom2, 0)
        pair_flip = pair_relative_flips.get((chrom1, chrom2), 0)

        whole = align_p4_to_truth_gauge(p4, swap1, swap2)
        whole_pred = int(np.argmax(whole))
        oracle = align_p4_to_truth_gauge(p4, swap1, swap2 ^ pair_flip)
        oracle_pred = int(np.argmax(oracle))
        nearest_geometry_pred = swap_state(int(np.argmin(np.asarray(distances, dtype=float))), swap1, swap2)
        truth_top = int(np.argmax(counts))
        posterior_top_model = int(np.argmax(p4))
        posterior_top_truth = swap_state(posterior_top_model, swap1, swap2)
        truth_model_state = swap_state(truth_top, swap1, swap2)
        truth_rank = state_distance_rank(distances, truth_model_state)
        posterior_rank = state_distance_rank(distances, posterior_top_model)
        gap = second_distance_gap(distances, posterior_top_model)
        pred_same = state_same_cross(whole_pred)
        same_contacts = int(counts[0] + counts[3])
        cross_contacts = int(counts[1] + counts[2])
        same_cross_correct = same_contacts if pred_same else cross_contacts
        whole_correct = int(counts[whole_pred])
        oracle_correct = int(counts[oracle_pred])
        nearest_geometry_correct = int(counts[nearest_geometry_pred])

        for acc in (all_acc, by_pair[(chrom1, chrom2)]):
            acc.update(
                weight=weight,
                whole_correct=whole_correct,
                oracle_correct=oracle_correct,
                nearest_geometry_correct=nearest_geometry_correct,
                same_cross_correct=same_cross_correct,
                pmax=float(rec["pmax"]),
                truth_rank=truth_rank,
                posterior_rank=posterior_rank,
                gap=gap,
            )
        pmax = float(rec["pmax"])
        if pmax >= 0.9:
            pbin = "pmax_ge_0p9"
        elif pmax >= 0.75:
            pbin = "pmax_0p75_0p9"
        elif pmax >= 0.5:
            pbin = "pmax_0p5_0p75"
        else:
            pbin = "pmax_lt_0p5"
        by_pmax[pbin].update(
            weight=weight,
            whole_correct=whole_correct,
            oracle_correct=oracle_correct,
            nearest_geometry_correct=nearest_geometry_correct,
            same_cross_correct=same_cross_correct,
            pmax=pmax,
            truth_rank=truth_rank,
            posterior_rank=posterior_rank,
            gap=gap,
        )

        bpair_rows.append({
            "chr1": chrom1,
            "start1": int(rec["start1"]),
            "chr2": chrom2,
            "start2": int(rec["start2"]),
            "n_contacts": weight,
            "truth_top_state": STATE_LABELS[truth_top],
            "truth_top_count": int(counts[truth_top]),
            "posterior_top_state": STATE_LABELS[posterior_top_truth],
            "posterior_pmax": pmax,
            "posterior_margin": rec["margin"],
            "whole_chrom_pred_state": STATE_LABELS[whole_pred],
            "whole_chrom_correct_contacts": whole_correct,
            "chr_pair_oracle_correct_contacts": oracle_correct,
            "nearest_geometry_state": STATE_LABELS[nearest_geometry_pred],
            "nearest_geometry_correct_contacts": nearest_geometry_correct,
            "same_cross_correct_contacts": same_cross_correct,
            "truth_state_distance_rank": truth_rank,
            "posterior_top1_distance_rank": posterior_rank,
            "posterior_top1_distance_gap": gap,
            "nearest_geometry_state": STATE_LABELS[swap_state(int(np.argmin(np.asarray(distances, dtype=float))), swap1, swap2)],
            "best_relative_flip_chr2": pair_flip,
        })

    summary_rows = [all_acc.row("genome_trans")]
    for key in ("pmax_ge_0p9", "pmax_0p75_0p9", "pmax_0p5_0p75", "pmax_lt_0p5"):
        if key in by_pmax:
            summary_rows.append(by_pmax[key].row(key))

    pair_rows: list[dict[str, object]] = []
    for (chrom1, chrom2), acc in sorted(by_pair.items(), key=lambda item: (chrom_sort_key(item[0][0]), chrom_sort_key(item[0][1]))):
        row = acc.row("chr_pair")
        pair_rows.append({
            "chr1": chrom1,
            "chr2": chrom2,
            "n_bpair": row["n_bpair"],
            "n_contacts": row["n_contacts"],
            "whole_chrom_top1_accuracy": row["whole_chrom_top1_accuracy"],
            "chr_pair_oracle_top1_accuracy": row["chr_pair_oracle_top1_accuracy"],
            "best_relative_flip_chr2": pair_relative_flips.get((chrom1, chrom2), 0),
            "nearest_geometry_state_accuracy": row["nearest_geometry_state_accuracy"],
            "same_cross_accuracy": row["same_cross_accuracy"],
            "pmax90_accuracy": row["pmax90_accuracy"],
            "pmax90_recall": row["pmax90_recall"],
            "truth_state_distance_rank_mean": row["truth_state_distance_rank_mean"],
            "truth_state_is_nearest_fraction": row["truth_state_is_nearest_fraction"],
            "posterior_top1_distance_rank_mean": row["posterior_top1_distance_rank_mean"],
            "posterior_top1_is_nearest_fraction": row["posterior_top1_is_nearest_fraction"],
            "posterior_top1_distance_gap_mean": row["posterior_top1_distance_gap_mean"],
        })
    return summary_rows, pair_rows, bpair_rows


def write_rows(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: format_value(row.get(col)) for col in columns})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config-output-dir", type=Path, required=True)
    parser.add_argument("--eval-dir", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, default=None)
    parser.add_argument("--reference", type=Path, default=None, help="Optional CHARM/3DG reference used only to match the standard eval denominator.")
    parser.add_argument("--bin-size", type=int, default=None)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    config_dir = args.config_output_dir
    eval_dir = args.eval_dir
    summary = read_kv(eval_dir / "summary.tsv")
    pairs = args.pairs or Path(summary["pairs_path"])
    bin_size = args.bin_size or int(summary["bin_size_bp"])
    reference_path = args.reference or Path(summary["reference_3dg_path"])
    posterior = read_posterior(config_dir / "p9016_full.bpair_posterior.tsv")
    coords = read_reconstruction(config_dir / "p9016_full.coords.tsv")
    reference = read_reference_3dg(reference_path, bin_size) if reference_path else None
    swaps = read_selected_swaps(eval_dir / "whole_chrom_snp_oracle_swaps.tsv")
    contact_counts = read_contact_truth_counts(pairs, bin_size, posterior)
    records = load_records(posterior, contact_counts, coords, reference)
    pair_flips = choose_chr_pair_oracle(records, swaps)
    summary_rows, pair_rows, bpair_rows = summarize_records(records, swaps, pair_flips)

    write_rows(args.outdir / "trans_gauge_geometry_summary.tsv", SUMMARY_COLUMNS, summary_rows)
    write_rows(args.outdir / "trans_gauge_geometry_chr_pair.tsv", CHR_PAIR_COLUMNS, pair_rows)
    write_rows(args.outdir / "trans_gauge_geometry_bpair.tsv", BPAIR_COLUMNS, bpair_rows)
    print(f"wrote\t{args.outdir / 'trans_gauge_geometry_summary.tsv'}")
    print(f"wrote\t{args.outdir / 'trans_gauge_geometry_chr_pair.tsv'}")
    print(f"wrote\t{args.outdir / 'trans_gauge_geometry_bpair.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
