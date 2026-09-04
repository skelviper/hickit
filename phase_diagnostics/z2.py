"""Blind chromosome-level Z2 synchronization with uncertainty diagnostics."""

from __future__ import annotations

import itertools
import math
from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Iterable, Mapping

import numpy as np


@dataclass(frozen=True)
class RelativeFlipEdge:
    chrom1: str
    chrom2: str
    q: float
    weight: float = 1.0

    def validate(self) -> None:
        if self.chrom1 == self.chrom2:
            raise ValueError("Z2 edges must connect distinct chromosomes")
        if not math.isfinite(self.q) or not -1 <= self.q <= 1:
            raise ValueError("q must be a calibrated signed estimate in [-1,1]")
        if not math.isfinite(self.weight) or self.weight < 0:
            raise ValueError("weight must be finite and non-negative")


@dataclass(frozen=True)
class Z2Result:
    spins: dict[str, int]
    component: dict[str, int]
    margins: dict[str, float]
    unresolved: frozenset[str]
    objective: float
    frustrated_edge_fraction: float
    weighted_frustrated_edge_fraction: float
    cycle_inconsistency: float


def _components(chromosomes: set[str], edges: list[RelativeFlipEdge]) -> list[list[str]]:
    adjacency: dict[str, set[str]] = defaultdict(set)
    for edge in edges:
        if edge.weight > 0 and edge.q != 0:
            adjacency[edge.chrom1].add(edge.chrom2)
            adjacency[edge.chrom2].add(edge.chrom1)
    result: list[list[str]] = []
    seen: set[str] = set()
    for root in sorted(chromosomes):
        if root in seen:
            continue
        queue = deque([root])
        seen.add(root)
        component: list[str] = []
        while queue:
            node = queue.popleft()
            component.append(node)
            for neighbor in sorted(adjacency[node]):
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        result.append(component)
    return result


def _objective(spins: Mapping[str, int], edges: Iterable[RelativeFlipEdge]) -> float:
    return float(sum(edge.weight * edge.q * spins[edge.chrom1] * spins[edge.chrom2] for edge in edges))


def _solve_component(nodes: list[str], edges: list[RelativeFlipEdge]) -> dict[str, int]:
    if len(nodes) == 1:
        return {nodes[0]: 1}
    if len(nodes) <= 16:
        best_score = -math.inf
        best: dict[str, int] | None = None
        root = nodes[0]
        for bits in itertools.product((-1, 1), repeat=len(nodes) - 1):
            candidate = {root: 1, **{node: bit for node, bit in zip(nodes[1:], bits, strict=True)}}
            score = _objective(candidate, edges)
            if score > best_score:
                best_score, best = score, candidate
        assert best is not None
        return best

    index = {node: idx for idx, node in enumerate(nodes)}
    matrix = np.zeros((len(nodes), len(nodes)), dtype=float)
    for edge in edges:
        i, j = index[edge.chrom1], index[edge.chrom2]
        matrix[i, j] += edge.weight * edge.q
        matrix[j, i] += edge.weight * edge.q
    _, vectors = np.linalg.eigh(matrix)
    leading = vectors[:, -1]
    spins = {node: (1 if leading[index[node]] >= 0 else -1) for node in nodes}
    if spins[nodes[0]] < 0:
        spins = {node: -value for node, value in spins.items()}
    changed = True
    while changed:
        changed = False
        for node in nodes[1:]:
            local_field = sum(
                edge.weight * edge.q * spins[edge.chrom2 if edge.chrom1 == node else edge.chrom1]
                for edge in edges if edge.chrom1 == node or edge.chrom2 == node
            )
            preferred = 1 if local_field >= 0 else -1
            if preferred != spins[node]:
                spins[node] = preferred
                changed = True
    return spins


def synchronize_z2(
    edges: Iterable[RelativeFlipEdge], *, chromosomes: Iterable[str] | None = None,
    minimum_weighted_degree: float = 0.25, minimum_margin: float = 0.05,
) -> Z2Result:
    edge_list = list(edges)
    for edge in edge_list:
        edge.validate()
    nodes = set(chromosomes or ())
    for edge in edge_list:
        nodes.update((edge.chrom1, edge.chrom2))
    components = _components(nodes, edge_list)
    spins: dict[str, int] = {}
    component_ids: dict[str, int] = {}
    for component_id, component_nodes in enumerate(components):
        component_edges = [
            edge for edge in edge_list if edge.chrom1 in component_nodes and edge.chrom2 in component_nodes
        ]
        solved = _solve_component(component_nodes, component_edges)
        spins.update(solved)
        component_ids.update({node: component_id for node in component_nodes})

    margins: dict[str, float] = {}
    degrees: dict[str, float] = defaultdict(float)
    for node in nodes:
        local_field = 0.0
        for edge in edge_list:
            if edge.chrom1 == node:
                local_field += edge.weight * edge.q * spins[edge.chrom2]
                degrees[node] += edge.weight * abs(edge.q)
            elif edge.chrom2 == node:
                local_field += edge.weight * edge.q * spins[edge.chrom1]
                degrees[node] += edge.weight * abs(edge.q)
        margins[node] = 2.0 * abs(local_field)
    unresolved = frozenset(
        node for node in nodes if degrees[node] < minimum_weighted_degree or margins[node] < minimum_margin
    )

    active = [edge for edge in edge_list if edge.weight > 0 and edge.q != 0]
    frustrated = [edge for edge in active if edge.q * spins[edge.chrom1] * spins[edge.chrom2] < 0]
    total_weight = sum(edge.weight * abs(edge.q) for edge in active)
    frustrated_weight = sum(edge.weight * abs(edge.q) for edge in frustrated)

    # Aggregate duplicate pair evidence before checking a deterministic
    # spanning-forest basis. Every non-tree edge closes one fundamental cycle,
    # so frustrated four-cycles and longer sparse cycles are visible too.
    pair_score: dict[tuple[str, str], float] = defaultdict(float)
    for edge in active:
        key = (min(edge.chrom1, edge.chrom2), max(edge.chrom1, edge.chrom2))
        pair_score[key] += edge.weight * edge.q
    signs = {key: (1 if value >= 0 else -1) for key, value in pair_score.items() if value != 0}
    adjacency: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for (left, right), sign in signs.items():
        adjacency[left].append((right, sign))
        adjacency[right].append((left, sign))
    tree_edges: set[tuple[str, str]] = set()
    path_spin: dict[str, int] = {}
    for root in sorted(nodes):
        if root in path_spin:
            continue
        path_spin[root] = 1
        queue = deque([root])
        while queue:
            node = queue.popleft()
            for neighbor, sign in sorted(adjacency[node]):
                key = (min(node, neighbor), max(node, neighbor))
                if neighbor not in path_spin:
                    path_spin[neighbor] = path_spin[node] * sign
                    tree_edges.add(key)
                    queue.append(neighbor)
    cycle_total = cycle_bad = 0
    for (left, right), sign in signs.items():
        if (left, right) in tree_edges:
            continue
        cycle_total += 1
        cycle_bad += path_spin[left] * path_spin[right] != sign
    return Z2Result(
        spins=spins,
        component=component_ids,
        margins=margins,
        unresolved=unresolved,
        objective=_objective(spins, edge_list),
        frustrated_edge_fraction=len(frustrated) / len(active) if active else math.nan,
        weighted_frustrated_edge_fraction=frustrated_weight / total_weight if total_weight else math.nan,
        cycle_inconsistency=cycle_bad / cycle_total if cycle_total else math.nan,
    )
