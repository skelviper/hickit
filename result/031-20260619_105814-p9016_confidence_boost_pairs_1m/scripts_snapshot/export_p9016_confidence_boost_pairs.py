#!/usr/bin/env python3
"""Create blind posterior-confidence boosted pairs from P9016 raw pairs.

The script keeps every original raw pair once and adds extra copies of raw
contacts whose 1 Mb trans bpair is selected by a blind posterior confidence
score. It never reads SNP truth or CHARM/3DG reference information.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


STATE_COLUMNS = ("p00", "p01", "p10", "p11")


@dataclass(frozen=True)
class PosteriorRow:
    key: tuple[str, int, str, int]
    chrom1: str
    start1: int
    chrom2: str
    start2: int
    n_raw: int
    pmax: float
    margin: float
    entropy: float
    contact_class: str
    score: float


def open_text(path: Path, mode: str = "rt") -> TextIO:
    if "b" in mode:
        raise ValueError("open_text expects text mode")
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode)


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    tail = chrom[3:] if chrom.startswith("chr") else chrom
    if tail.isdigit():
        return int(tail), ""
    return 10_000, tail


def canonical_bpair_key(chrom1: str, start1: int, chrom2: str, start2: int) -> tuple[str, int, str, int]:
    a = (chrom_sort_key(chrom1), start1, chrom1)
    b = (chrom_sort_key(chrom2), start2, chrom2)
    if a <= b:
        return chrom1, start1, chrom2, start2
    return chrom2, start2, chrom1, start1


def score_row(row: dict[str, str], mode: str) -> float:
    pmax = float(row["pmax"])
    margin = float(row.get("margin", "0") or "0")
    entropy = float(row.get("entropy", "nan") or "nan")
    if mode == "pmax":
        return pmax
    if mode == "margin":
        return margin
    if mode == "pmax_margin":
        return pmax * margin
    if mode == "pmax_margin_entropy":
        if not math.isfinite(entropy):
            entropy = math.log(4.0)
        return pmax * margin * (1.0 - min(max(entropy / math.log(4.0), 0.0), 1.0))
    raise ValueError(f"unknown score mode: {mode}")


def read_trans_posterior(path: Path, score_mode: str) -> list[PosteriorRow]:
    rows: list[PosteriorRow] = []
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"chr1", "start1", "chr2", "start2", "n_raw", "pmax", "contact_class", *STATE_COLUMNS}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"{path} is missing required columns: {sorted(required)}")
        for row in reader:
            chrom1 = row["chr1"]
            chrom2 = row["chr2"]
            if chrom1 == chrom2:
                continue
            if row.get("contact_class") != "trans":
                continue
            start1 = int(row["start1"])
            start2 = int(row["start2"])
            p4 = [float(row[col]) for col in STATE_COLUMNS]
            p_sorted = sorted(p4, reverse=True)
            margin = float(row.get("margin", "") or (p_sorted[0] - p_sorted[1]))
            enriched = dict(row)
            enriched["margin"] = str(margin)
            key = canonical_bpair_key(chrom1, start1, chrom2, start2)
            rows.append(
                PosteriorRow(
                    key=key,
                    chrom1=chrom1,
                    start1=start1,
                    chrom2=chrom2,
                    start2=start2,
                    n_raw=int(row["n_raw"]),
                    pmax=float(row["pmax"]),
                    margin=margin,
                    entropy=float(row.get("entropy", "nan") or "nan"),
                    contact_class=row["contact_class"],
                    score=score_row(enriched, score_mode),
                )
            )
    return rows


def select_keys(rows: list[PosteriorRow], top_frac: float, min_pmax: float, min_margin: float) -> tuple[set[tuple[str, int, str, int]], dict[str, object]]:
    eligible = [
        row for row in rows
        if math.isfinite(row.score)
        and row.pmax >= min_pmax
        and row.margin >= min_margin
    ]
    eligible.sort(key=lambda row: (row.score, row.pmax, row.margin, row.n_raw), reverse=True)
    n_select = int(math.ceil(len(eligible) * top_frac))
    if top_frac > 0.0 and eligible:
        n_select = max(1, n_select)
    selected = eligible[:n_select]
    threshold = selected[-1].score if selected else "NA"
    return {row.key for row in selected}, {
        "posterior_trans_bpairs": len(rows),
        "eligible_trans_bpairs": len(eligible),
        "selected_trans_bpairs": len(selected),
        "score_threshold": threshold,
        "selected_score_min": threshold,
        "selected_score_max": selected[0].score if selected else "NA",
        "selected_pmax_min": min((row.pmax for row in selected), default="NA"),
        "selected_margin_min": min((row.margin for row in selected), default="NA"),
        "selected_n_raw_sum": sum(row.n_raw for row in selected),
    }


def parse_columns(line: str) -> list[str]:
    if line.startswith("#columns:"):
        return line.split(":", 1)[1].strip().split()
    raise ValueError("not a #columns line")


def update_columns_line(columns: list[str], phase_mode: str) -> tuple[str, list[str]]:
    out_cols = list(columns)
    if phase_mode == "drop":
        out_cols = [col for col in out_cols if col not in {"phase0", "phase1"}]
    elif phase_mode == "dot":
        for col in ("phase0", "phase1"):
            if col not in out_cols:
                out_cols.append(col)
    elif phase_mode == "copy":
        pass
    else:
        raise ValueError(f"unknown phase mode: {phase_mode}")
    return "#columns: " + " ".join(out_cols), out_cols


def line_from_row(row: dict[str, str], columns: list[str], copy_index: int, phase_mode: str) -> str:
    if copy_index > 0:
        row = dict(row)
        row["readID"] = f"{row.get('readID', '.')}:hkboost{copy_index}"
    if phase_mode == "dot":
        row = dict(row)
        if "phase0" in columns:
            row["phase0"] = "."
        if "phase1" in columns:
            row["phase1"] = "."
    return "\t".join(row.get(col, ".") for col in columns)


def iter_pair_rows(path: Path) -> Iterable[tuple[str, list[str] | None]]:
    with open_text(path) as fh:
        for line in fh:
            yield line.rstrip("\n"), None


def write_boosted_pairs(
    pairs_path: Path,
    out_path: Path,
    selected_keys: set[tuple[str, int, str, int]],
    bin_size: int,
    boost_copies: int,
    phase_mode: str,
) -> dict[str, object]:
    if boost_copies < 1:
        raise ValueError("--boost-copies must be >= 1")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if out_path.suffix == ".gz" else open
    columns: list[str] | None = None
    out_columns: list[str] | None = None
    stats: dict[str, object] = {
        "input_pairs": 0,
        "input_cis_pairs": 0,
        "input_trans_pairs": 0,
        "selected_trans_raw_pairs": 0,
        "output_pairs": 0,
        "extra_copies_written": 0,
        "boost_copies": boost_copies,
        "phase_mode": phase_mode,
    }
    with open_text(pairs_path) as inp, opener(out_path, "wt") as out:
        for line in inp:
            line = line.rstrip("\n")
            if line.startswith("#columns:"):
                columns = parse_columns(line)
                columns_line, out_columns = update_columns_line(columns, phase_mode)
                out.write(columns_line + "\n")
                continue
            if line.startswith("#"):
                out.write(line + "\n")
                continue
            if not line:
                continue
            if columns is None:
                columns = ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]
                _, out_columns = update_columns_line(columns, phase_mode)
            assert out_columns is not None
            fields = line.split()
            row = dict(zip(columns, fields))
            chrom1 = row["chr1"]
            chrom2 = row["chr2"]
            start1 = (int(row["pos1"]) // bin_size) * bin_size
            start2 = (int(row["pos2"]) // bin_size) * bin_size
            key = canonical_bpair_key(chrom1, start1, chrom2, start2)
            is_trans = chrom1 != chrom2
            is_selected = is_trans and key in selected_keys
            stats["input_pairs"] = int(stats["input_pairs"]) + 1
            stats["input_trans_pairs" if is_trans else "input_cis_pairs"] = int(stats["input_trans_pairs" if is_trans else "input_cis_pairs"]) + 1
            copies = boost_copies if is_selected else 1
            if is_selected:
                stats["selected_trans_raw_pairs"] = int(stats["selected_trans_raw_pairs"]) + 1
            for copy_index in range(copies):
                out.write(line_from_row(row, out_columns, copy_index, phase_mode) + "\n")
            stats["output_pairs"] = int(stats["output_pairs"]) + copies
            stats["extra_copies_written"] = int(stats["extra_copies_written"]) + (copies - 1)
    return stats


def write_metadata(path: Path, values: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        for key in sorted(values):
            value = values[key]
            if isinstance(value, float):
                if math.isfinite(value):
                    value = f"{value:.12g}"
                else:
                    value = "NA"
            writer.writerow([key, value])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True, help="Original P9016 pairs.gz file.")
    parser.add_argument("--posterior", type=Path, required=True, help="Blind baseline bpair posterior TSV.")
    parser.add_argument("--out-pairs", type=Path, required=True, help="Output boosted pairs path, usually .pairs.gz.")
    parser.add_argument("--metadata", type=Path, required=True, help="Output key/value metadata TSV.")
    parser.add_argument("--selected-keys-json", type=Path, default=None, help="Optional JSON list of selected bpair keys.")
    parser.add_argument("--bin-size", type=int, default=1_000_000)
    parser.add_argument("--score", choices=["pmax", "margin", "pmax_margin", "pmax_margin_entropy"], default="pmax_margin")
    parser.add_argument("--trans-top-frac", type=float, required=True)
    parser.add_argument("--min-pmax", type=float, default=0.0)
    parser.add_argument("--min-margin", type=float, default=0.0)
    parser.add_argument("--boost-copies", type=int, required=True, help="Total copies for selected raw contacts, including the original.")
    parser.add_argument("--phase-mode", choices=["dot", "drop", "copy"], default="dot")
    args = parser.parse_args()

    if args.bin_size <= 0:
        raise SystemExit("--bin-size must be positive")
    if not (0.0 <= args.trans_top_frac <= 1.0):
        raise SystemExit("--trans-top-frac must be between 0 and 1")
    if args.min_pmax < 0.0 or args.min_margin < 0.0:
        raise SystemExit("--min-pmax and --min-margin must be non-negative")

    posterior_rows = read_trans_posterior(args.posterior, args.score)
    selected_keys, selection_stats = select_keys(
        posterior_rows,
        args.trans_top_frac,
        args.min_pmax,
        args.min_margin,
    )
    pair_stats = write_boosted_pairs(
        args.pairs,
        args.out_pairs,
        selected_keys,
        args.bin_size,
        args.boost_copies,
        args.phase_mode,
    )
    if args.selected_keys_json:
        args.selected_keys_json.parent.mkdir(parents=True, exist_ok=True)
        payload = [
            {"chr1": key[0], "start1": key[1], "chr2": key[2], "start2": key[3]}
            for key in sorted(selected_keys, key=lambda k: (chrom_sort_key(k[0]), k[1], chrom_sort_key(k[2]), k[3]))
        ]
        args.selected_keys_json.write_text(json.dumps(payload, indent=2) + "\n")

    metadata: dict[str, object] = {
        "method": "blind_posterior_confidence_trans_boost",
        "uses_phase_labels_for_training": 0,
        "uses_charm_or_reference_for_training": 0,
        "pairs": str(args.pairs),
        "posterior": str(args.posterior),
        "out_pairs": str(args.out_pairs),
        "bin_size_bp": args.bin_size,
        "score_mode": args.score,
        "trans_top_frac": args.trans_top_frac,
        "min_pmax": args.min_pmax,
        "min_margin": args.min_margin,
        "boost_copies": args.boost_copies,
        "phase_mode": args.phase_mode,
    }
    metadata.update(selection_stats)
    metadata.update(pair_stats)
    write_metadata(args.metadata, metadata)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
