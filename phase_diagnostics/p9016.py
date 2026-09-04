"""P9016 truth readers used only by explicitly ORACLE evaluation commands."""

from __future__ import annotations

import csv
import gzip
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np


COPY_PATTERN = re.compile(r"^(?P<chrom>.+)\((?P<copy>mat|pat|copy0|copy1|0|1)\)$")
HICKIT_COPY_PATTERN = re.compile(r"^(?P<chrom>chr(?:[0-9]+|X|Y|M|MT))(?P<copy>a|b)$")


@dataclass
class PhaseContactData:
    chrom1: np.ndarray
    start1: np.ndarray
    chrom2: np.ndarray
    start2: np.ndarray
    truth_state: np.ndarray
    distances: np.ndarray
    multiplicity: np.ndarray
    genomic_separation: np.ndarray
    is_cis: np.ndarray

    @property
    def n_contacts(self) -> int:
        return int(self.truth_state.size)


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def read_3dg(path: Path, resolution: int) -> dict[tuple[str, int, int], np.ndarray]:
    grouped: dict[tuple[str, int, int], list[np.ndarray]] = {}
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                raise ValueError(f"{path}:{line_number}: expected five columns")
            match = COPY_PATTERN.match(fields[0])
            if not match:
                match = HICKIT_COPY_PATTERN.match(fields[0])
            if not match:
                raise ValueError(f"{path}:{line_number}: chromosome lacks a diploid copy suffix")
            copy_token = match.group("copy")
            copy = 0 if copy_token in {"mat", "copy0", "0", "a"} else 1
            start = int(fields[1]) // resolution * resolution
            value = np.asarray([float(item) for item in fields[2:5]], dtype=float)
            if not np.isfinite(value).all():
                raise ValueError(f"{path}:{line_number}: non-finite coordinate")
            grouped.setdefault((match.group("chrom"), start, copy), []).append(value)
    return {key: np.vstack(values).mean(axis=0) for key, values in grouped.items()}


def _pairs_columns(path: Path) -> tuple[list[str], list[str]]:
    headers: list[str] = []
    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            headers.append(line.rstrip("\n"))
    columns_lines = [line for line in headers if line.startswith("#columns:")]
    if len(columns_lines) != 1:
        raise ValueError(f"{path}: expected exactly one #columns header")
    columns = columns_lines[0].split(":", 1)[1].strip().split()
    required = {"chr1", "pos1", "chr2", "pos2", "phase0", "phase1"}
    if not required.issubset(columns):
        raise ValueError(f"{path}: phase-bearing ORACLE input lacks required columns")
    return headers, columns


def load_phase_contacts(
    pairs_path: Path,
    reference_path: Path,
    *,
    resolution: int = 1_000_000,
    max_contacts: int | None = None,
) -> PhaseContactData:
    """Load measured double-end phase and CHARM coordinates for post hoc evaluation.

    This function must never be called from a blind training entry point.
    """
    _, columns = _pairs_columns(pairs_path)
    by_name = {name: index for index, name in enumerate(columns)}
    reference = read_3dg(reference_path, resolution)
    multiplicity: Counter[tuple[str, int, str, int]] = Counter()
    candidates: list[tuple[str, int, str, int, int]] = []
    with _open_text(pairs_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom1, chrom2 = fields[by_name["chr1"]], fields[by_name["chr2"]]
            start1 = int(fields[by_name["pos1"]]) // resolution * resolution
            start2 = int(fields[by_name["pos2"]]) // resolution * resolution
            key = (chrom1, start1, chrom2, start2)
            multiplicity[key] += 1
            phase1, phase2 = fields[by_name["phase0"]], fields[by_name["phase1"]]
            if phase1 not in {"0", "1"} or phase2 not in {"0", "1"}:
                continue
            if start1 == start2 and chrom1 == chrom2:
                continue
            if not all((chrom, start, copy) in reference for chrom, start in ((chrom1, start1), (chrom2, start2)) for copy in (0, 1)):
                continue
            candidates.append((chrom1, start1, chrom2, start2, 2 * int(phase1) + int(phase2)))
            if max_contacts is not None and len(candidates) >= max_contacts:
                break
    if not candidates:
        raise ValueError("no double-end phased contacts overlap the reference coordinates")
    chrom1_values: list[str] = []
    chrom2_values: list[str] = []
    start1_values: list[int] = []
    start2_values: list[int] = []
    truth: list[int] = []
    distances: list[np.ndarray] = []
    counts: list[int] = []
    genomic: list[int] = []
    cis_values: list[bool] = []
    for chrom1, start1, chrom2, start2, state in candidates:
        distance = np.asarray([
            np.linalg.norm(reference[(chrom1, start1, copy1)] - reference[(chrom2, start2, copy2)])
            for copy1 in (0, 1) for copy2 in (0, 1)
        ])
        chrom1_values.append(chrom1)
        chrom2_values.append(chrom2)
        start1_values.append(start1)
        start2_values.append(start2)
        truth.append(state)
        distances.append(distance)
        counts.append(multiplicity[(chrom1, start1, chrom2, start2)])
        is_cis = chrom1 == chrom2
        cis_values.append(is_cis)
        genomic.append(abs(start2 - start1) if is_cis else -1)
    return PhaseContactData(
        chrom1=np.asarray(chrom1_values, dtype=object),
        start1=np.asarray(start1_values, dtype=np.int64),
        chrom2=np.asarray(chrom2_values, dtype=object),
        start2=np.asarray(start2_values, dtype=np.int64),
        truth_state=np.asarray(truth, dtype=np.int8),
        distances=np.vstack(distances),
        multiplicity=np.asarray(counts, dtype=np.int32),
        genomic_separation=np.asarray(genomic, dtype=np.int64),
        is_cis=np.asarray(cis_values, dtype=bool),
    )
