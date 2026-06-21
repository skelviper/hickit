#!/usr/bin/env python3
"""Infer CHARM/3DG copy states for unphased synthetic pairs.

This is a QC-only helper. It does not alter the unphased training pairs; it
reconstructs the copy pair that most plausibly generated each synthetic contact
from the reference 3DG coordinates.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
from pathlib import Path

import numpy as np


COPY_RE = re.compile(r"^(?P<chrom>.+)\((?P<copy>mat|pat|copy0|copy1|0|1)\)$")


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def copy_to_int(value: str) -> int:
    if value in {"mat", "copy0", "0"}:
        return 0
    if value in {"pat", "copy1", "1"}:
        return 1
    raise ValueError(f"unsupported copy label: {value}")


def split_chrom_copy(raw: str) -> tuple[str, int]:
    match = COPY_RE.match(raw)
    if not match:
        return raw, 0
    return match.group("chrom"), copy_to_int(match.group("copy"))


def read_3dg(path: Path) -> dict[tuple[str, int, int], np.ndarray]:
    coords: dict[tuple[str, int, int], np.ndarray] = {}
    with open_text(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                raise ValueError(f"{path}:{line_no} expected at least 5 columns")
            chrom, copy = split_chrom_copy(fields[0])
            start = int(fields[1])
            xyz = np.array([float(fields[2]), float(fields[3]), float(fields[4])], dtype=np.float32)
            if not np.all(np.isfinite(xyz)):
                raise ValueError(f"{path}:{line_no} non-finite coordinate")
            coords[(chrom, start, copy)] = xyz
    if not coords:
        raise ValueError(f"no coordinates read from {path}")
    return coords


def infer_copy_pair(
    coords: dict[tuple[str, int, int], np.ndarray],
    chrom1: str,
    pos1: int,
    chrom2: str,
    pos2: int,
) -> tuple[int, int, float] | None:
    best: tuple[int, int, float] | None = None
    for copy1 in (0, 1):
        xyz1 = coords.get((chrom1, pos1, copy1))
        if xyz1 is None:
            continue
        for copy2 in (0, 1):
            xyz2 = coords.get((chrom2, pos2, copy2))
            if xyz2 is None:
                continue
            dist = float(np.linalg.norm(xyz1 - xyz2))
            if not math.isfinite(dist):
                continue
            if best is None or dist < best[2]:
                best = (copy1, copy2, dist)
    return best


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tdg", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    coords = read_3dg(args.tdg)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] | None = None
    used = missing = 0
    with open_text(args.pairs) as pairs_fh, gzip.open(args.out, "wt") as out:
        out.write("readID\tchr1\tpos1\tcopy1\tchr2\tpos2\tcopy2\tdistance\n")
        for line in pairs_fh:
            if not line.strip():
                continue
            if line.startswith("#columns:"):
                columns = line.split(":", 1)[1].strip().split()
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                fields = line.split()
            if columns is None:
                columns = ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2", "phase0", "phase1"]
            row = dict(zip(columns, fields))
            read_id = row.get("readID", f"row_{used + missing}")
            chrom1, chrom2 = row["chr1"], row["chr2"]
            pos1, pos2 = int(row["pos1"]), int(row["pos2"])
            inferred = infer_copy_pair(coords, chrom1, pos1, chrom2, pos2)
            if inferred is None:
                missing += 1
                continue
            copy1, copy2, dist = inferred
            out.write(f"{read_id}\t{chrom1}\t{pos1}\t{copy1}\t{chrom2}\t{pos2}\t{copy2}\t{dist:.9g}\n")
            used += 1
    print(f"wrote\t{args.out}")
    print(f"used\t{used}")
    print(f"missing\t{missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
