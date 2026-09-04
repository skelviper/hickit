"""Canonical four-state semantics and chromosome-gauge relabeling."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


STATE_LABELS = ("00", "01", "10", "11")
STATE_TO_INDEX = {label: index for index, label in enumerate(STATE_LABELS)}


def state_bits(state: int | str) -> tuple[int, int]:
    index = STATE_TO_INDEX[state] if isinstance(state, str) else int(state)
    if index < 0 or index >= 4:
        raise ValueError(f"invalid four-state index: {state!r}")
    return index >> 1, index & 1


def state_index(left_copy: int, right_copy: int) -> int:
    if left_copy not in (0, 1) or right_copy not in (0, 1):
        raise ValueError("copy labels must be 0 or 1")
    return 2 * left_copy + right_copy


def swap_endpoint_state(state: int | str) -> int:
    left, right = state_bits(state)
    return state_index(right, left)


def flip_state(state: int | str, flip_left: int, flip_right: int) -> int:
    left, right = state_bits(state)
    return state_index(left ^ int(bool(flip_left)), right ^ int(bool(flip_right)))


def relabel_probabilities(
    probabilities: Sequence[float], *, swap_endpoints: bool = False,
    flip_left: int = 0, flip_right: int = 0,
) -> np.ndarray:
    """Relabel a state distribution without changing its statistical mass."""
    values = np.asarray(probabilities, dtype=float)
    if values.shape != (4,) or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("probabilities must be four finite non-negative values")
    out = np.zeros(4, dtype=float)
    for source in range(4):
        left, right = state_bits(source)
        if swap_endpoints:
            left, right = right, left
        target = state_index(left ^ int(bool(flip_left)), right ^ int(bool(flip_right)))
        out[target] += values[source]
    if not np.isclose(out.sum(), values.sum(), atol=1e-12, rtol=1e-12):
        raise AssertionError("state relabeling changed total probability")
    return out


def canonical_endpoint_key(
    chrom1: str, start1: int, chrom2: str, start2: int,
) -> tuple[tuple[str, int, str, int], bool]:
    left = (str(chrom1), int(start1))
    right = (str(chrom2), int(start2))
    if right < left:
        return (right[0], right[1], left[0], left[1]), True
    return (left[0], left[1], right[0], right[1]), False


@dataclass
class CanonicalAccumulator:
    """Aggregate raw-order state probabilities into canonical binned pairs."""

    sums: dict[tuple[str, int, str, int], np.ndarray] = field(default_factory=dict)
    counts: dict[tuple[str, int, str, int], int] = field(default_factory=dict)

    def add(
        self,
        chrom1: str,
        start1: int,
        chrom2: str,
        start2: int,
        probabilities: Sequence[float],
    ) -> None:
        key, swapped = canonical_endpoint_key(chrom1, start1, chrom2, start2)
        canonical = relabel_probabilities(probabilities, swap_endpoints=swapped)
        self.sums.setdefault(key, np.zeros(4, dtype=float))
        self.sums[key] += canonical
        self.counts[key] = self.counts.get(key, 0) + 1

    def normalized(self) -> dict[tuple[str, int, str, int], np.ndarray]:
        result: dict[tuple[str, int, str, int], np.ndarray] = {}
        for key, values in self.sums.items():
            total = float(values.sum())
            if total <= 0 or not np.isfinite(total):
                raise ValueError(f"invalid aggregated probability mass for {key}")
            result[key] = values / total
        return result


def write_p4(path: Path, values: Mapping[tuple[str, int, str, int], Sequence[float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("chrom1", "start1", "chrom2", "start2", "p00", "p01", "p10", "p11"))
        for key in sorted(values):
            p4 = np.asarray(values[key], dtype=float)
            if p4.shape != (4,) or not np.isclose(p4.sum(), 1.0, atol=1e-7):
                raise ValueError(f"non-normalized p4 for {key}")
            writer.writerow((*key, *(format(float(item), ".17g") for item in p4)))


def read_p4(path: Path) -> dict[tuple[str, int, str, int], np.ndarray]:
    result: dict[tuple[str, int, str, int], np.ndarray] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom1", "start1", "chrom2", "start2", "p00", "p01", "p10", "p11"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: invalid p4 header")
        for line_number, row in enumerate(reader, start=2):
            key, swapped = canonical_endpoint_key(
                row["chrom1"], int(row["start1"]), row["chrom2"], int(row["start2"])
            )
            values = np.asarray([float(row[f"p{label}"]) for label in STATE_LABELS], dtype=float)
            values = relabel_probabilities(values, swap_endpoints=swapped)
            if not np.isclose(values.sum(), 1.0, atol=1e-7) or (values < 0).any():
                raise ValueError(f"{path}:{line_number}: invalid p4")
            if key in result:
                raise ValueError(f"{path}:{line_number}: duplicate canonical pair")
            result[key] = values
    return result


def apply_chromosome_gauge(
    probabilities: Sequence[float], chrom1: str, chrom2: str,
    flips: Mapping[str, int],
) -> np.ndarray:
    return relabel_probabilities(
        probabilities,
        flip_left=int(bool(flips.get(chrom1, 0))),
        flip_right=int(bool(flips.get(chrom2, 0))),
    )


def parity_probabilities(probabilities: np.ndarray) -> np.ndarray:
    values = np.asarray(probabilities, dtype=float)
    if values.shape[-1] != 4:
        raise ValueError("last probability dimension must contain 00,01,10,11")
    return np.stack((values[..., 0] + values[..., 3], values[..., 1] + values[..., 2]), axis=-1)
