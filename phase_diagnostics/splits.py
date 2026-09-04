"""Leakage-aware held-out split units for contact observations."""

from __future__ import annotations

import hashlib
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable


SPLIT_MODES = (
    "raw_observation",
    "molecule_or_readchain",
    "exact_haploid_bin_pair",
    "genomic_block",
    "chromosome_pair",
)


@dataclass(frozen=True)
class Observation:
    observation_id: str
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    molecule_id: str | None = None


class _UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[right] = left


def _ordered_endpoints(observation: Observation) -> tuple[tuple[str, int], tuple[str, int]]:
    left = (observation.chrom1, int(observation.pos1))
    right = (observation.chrom2, int(observation.pos2))
    return (right, left) if right < left else (left, right)


def split_key(observation: Observation, mode: str, *, bin_size: int, block_size: int) -> tuple[object, ...]:
    if mode not in SPLIT_MODES:
        raise ValueError(f"unsupported split mode: {mode}")
    left, right = _ordered_endpoints(observation)
    if mode == "raw_observation":
        return ("observation", observation.observation_id)
    if mode == "molecule_or_readchain":
        return ("molecule", observation.molecule_id or observation.observation_id)
    if mode == "exact_haploid_bin_pair":
        return ("bpair", left[0], left[1] // bin_size, right[0], right[1] // bin_size)
    if mode == "genomic_block":
        return ("block", left[0], left[1] // block_size, right[0], right[1] // block_size)
    return ("chromosome_pair", left[0], right[0])


def blocked_split(
    observations: Iterable[Observation], *, mode: str, heldout_fraction: float, seed: int,
    bin_size: int = 1_000_000, block_size: int = 5_000_000,
) -> dict[str, str]:
    records = list(observations)
    if not 0 < heldout_fraction < 1:
        raise ValueError("heldout_fraction must be in (0,1)")
    if not records:
        return {}
    union = _UnionFind(len(records))
    by_primary: dict[tuple[object, ...], int] = {}
    by_molecule: dict[str, int] = {}
    for index, record in enumerate(records):
        primary = split_key(record, mode, bin_size=bin_size, block_size=block_size)
        if primary in by_primary:
            union.union(index, by_primary[primary])
        else:
            by_primary[primary] = index
        # This invariant is applied for every split mode, including raw_observation.
        if record.molecule_id:
            if record.molecule_id in by_molecule:
                union.union(index, by_molecule[record.molecule_id])
            else:
                by_molecule[record.molecule_id] = index

    components: dict[int, list[int]] = defaultdict(list)
    for index in range(len(records)):
        components[union.find(index)].append(index)
    ranked = sorted(
        components.values(),
        key=lambda members: hashlib.sha256(
            f"{seed}|{mode}|".encode() + "|".join(sorted(records[i].observation_id for i in members)).encode()
        ).digest(),
    )
    target = heldout_fraction * len(records)
    heldout_count = 0
    heldout_components: set[int] = set()
    for component_index, members in enumerate(ranked):
        before = abs(heldout_count - target)
        after = abs(heldout_count + len(members) - target)
        if after < before or heldout_count == 0:
            heldout_components.add(component_index)
            heldout_count += len(members)
    if len(heldout_components) == len(ranked) and len(ranked) > 1:
        heldout_components.remove(max(heldout_components))
    assignment: dict[str, str] = {}
    for component_index, members in enumerate(ranked):
        partition = "heldout" if component_index in heldout_components else "train"
        for index in members:
            assignment[records[index].observation_id] = partition
    return assignment


def leakage_audit(
    observations: Iterable[Observation], assignment: dict[str, str], *, bin_size: int = 1_000_000,
    neighborhood_bp: int = 1_000_000,
) -> dict[str, int | bool]:
    records = list(observations)
    groups: dict[str, dict[tuple[object, ...], set[str]]] = {
        "duplicate": defaultdict(set),
        "exact_haploid_bin_pair": defaultdict(set),
        "molecule": defaultdict(set),
        "neighbor_block": defaultdict(set),
    }
    for record in records:
        partition = assignment[record.observation_id]
        left, right = _ordered_endpoints(record)
        groups["duplicate"][(left, right)].add(partition)
        groups["exact_haploid_bin_pair"][(left[0], left[1] // bin_size, right[0], right[1] // bin_size)].add(partition)
        if record.molecule_id:
            groups["molecule"][(record.molecule_id,)].add(partition)
        groups["neighbor_block"][(
            left[0], left[1] // neighborhood_bp, right[0], right[1] // neighborhood_bp
        )].add(partition)
    result: dict[str, int | bool] = {}
    any_leak = False
    for name, mapping in groups.items():
        count = sum(len(partitions) > 1 for partitions in mapping.values())
        result[f"{name}_leakage_groups"] = count
        any_leak |= count > 0
    result["pass"] = not any_leak
    return result
