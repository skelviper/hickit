#!/usr/bin/env python3
"""Generate unphased synthetic pairs from CHARM/3DG geometry.

The output is a positive-control training input: coordinates from CHARM/3DG are
used to define contacts, but phase labels are deliberately omitted.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import json
import math
import random
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np


COPY_RE = re.compile(r"^(?P<chrom>.+)\((?P<copy>mat|pat|copy0|copy1|0|1)\)$")


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    tail = chrom[3:] if chrom.startswith("chr") else chrom
    if tail.isdigit():
        return (int(tail), "")
    if tail == "X":
        return (10_000, "")
    if tail == "Y":
        return (10_001, "")
    if tail in {"M", "MT"}:
        return (10_002, "")
    return (20_000, tail)


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


def read_3dg(path: Path) -> list[dict[str, object]]:
    beads: list[dict[str, object]] = []
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
            beads.append({"chrom": chrom, "copy": copy, "start": start, "xyz": xyz})
    if not beads:
        raise ValueError(f"no beads read from {path}")
    return beads


def chrom_sizes_from_beads(beads: list[dict[str, object]], tdg_resolution: int) -> dict[str, int]:
    sizes: dict[str, int] = {}
    for bead in beads:
        chrom = str(bead["chrom"])
        end = int(bead["start"]) + tdg_resolution
        sizes[chrom] = max(sizes.get(chrom, 0), end)
    return sizes


def estimate_particle_radius(beads: list[dict[str, object]]) -> float:
    by_track: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for bead in beads:
        by_track[(str(bead["chrom"]), int(bead["copy"]))].append(bead)
    distances: list[float] = []
    for track in by_track.values():
        track.sort(key=lambda b: int(b["start"]))
        for a, b in zip(track, track[1:]):
            if int(b["start"]) <= int(a["start"]):
                continue
            d = float(np.linalg.norm(np.asarray(b["xyz"]) - np.asarray(a["xyz"])))
            if math.isfinite(d) and d > 0.0:
                distances.append(d)
    if not distances:
        raise ValueError("cannot estimate particle radius: no positive adjacent distances")
    return float(np.median(np.asarray(distances, dtype=np.float64)))


def count_pairs(path: Path, bin_size: int) -> dict[str, int]:
    columns: list[str] | None = None
    total = cis = trans = same_bin = 0
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#columns:"):
                columns = line.split(":", 1)[1].strip().split()
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split()
            if columns is None:
                columns = ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]
            row = dict(zip(columns, fields))
            total += 1
            if row["chr1"] == row["chr2"]:
                cis += 1
                if int(row["pos1"]) // bin_size == int(row["pos2"]) // bin_size:
                    same_bin += 1
            else:
                trans += 1
    return {
        "observed_pairs_total": total,
        "observed_pairs_cis": cis,
        "observed_pairs_trans": trans,
        "observed_pairs_same_bin_at_training_resolution": same_bin,
        "observed_pairs_non_same_bin_at_training_resolution": total - same_bin,
    }


def reservoir_candidates(
    beads: list[dict[str, object]],
    radius: float,
    max_keep: int,
    seed: int,
    min_genomic_separation: int,
    training_bin_size: int,
    exclude_training_same_bin: bool,
) -> tuple[list[tuple[int, int, float]], dict[str, int]]:
    rng = random.Random(seed)
    radius2 = radius * radius
    cell_size = radius
    grid: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    candidates: list[tuple[int, int, float]] = []
    stats = {
        "candidate_pairs_total": 0,
        "candidate_pairs_cis": 0,
        "candidate_pairs_trans": 0,
        "candidate_pairs_skipped_same_copy_same_locus": 0,
        "candidate_pairs_skipped_opposite_copy_same_locus": 0,
        "candidate_pairs_skipped_min_genomic_separation": 0,
        "candidate_pairs_skipped_training_same_bin": 0,
    }

    def cell_of(xyz: np.ndarray) -> tuple[int, int, int]:
        return tuple(int(math.floor(float(v) / cell_size)) for v in xyz)

    for i, bead in enumerate(beads):
        xyz = np.asarray(bead["xyz"])
        cx, cy, cz = cell_of(xyz)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((cx + dx, cy + dy, cz + dz), []):
                        other = beads[j]
                        if bead["copy"] == other["copy"] and bead["chrom"] == other["chrom"] and bead["start"] == other["start"]:
                            stats["candidate_pairs_skipped_same_copy_same_locus"] += 1
                            continue
                        if bead["chrom"] == other["chrom"] and bead["start"] == other["start"]:
                            stats["candidate_pairs_skipped_opposite_copy_same_locus"] += 1
                            continue
                        if (
                            exclude_training_same_bin
                            and bead["chrom"] == other["chrom"]
                            and int(bead["start"]) // training_bin_size == int(other["start"]) // training_bin_size
                        ):
                            stats["candidate_pairs_skipped_training_same_bin"] += 1
                            continue
                        if bead["chrom"] == other["chrom"] and abs(int(bead["start"]) - int(other["start"])) < min_genomic_separation:
                            stats["candidate_pairs_skipped_min_genomic_separation"] += 1
                            continue
                        d2 = float(np.sum((xyz - np.asarray(other["xyz"])) ** 2))
                        if d2 > radius2:
                            continue
                        stats["candidate_pairs_total"] += 1
                        if bead["chrom"] == other["chrom"]:
                            stats["candidate_pairs_cis"] += 1
                        else:
                            stats["candidate_pairs_trans"] += 1
                        dist = math.sqrt(d2)
                        if len(candidates) < max_keep:
                            candidates.append((j, i, dist))
                        else:
                            k = stats["candidate_pairs_total"]
                            replace = rng.randrange(k)
                            if replace < max_keep:
                                candidates[replace] = (j, i, dist)
        grid[(cx, cy, cz)].append(i)
    return candidates, stats


def nearest_candidates(
    beads: list[dict[str, object]],
    max_keep: int,
    min_genomic_separation: int,
    training_bin_size: int,
    exclude_training_same_bin: bool,
) -> tuple[list[tuple[int, int, float]], dict[str, int]]:
    heap: list[tuple[float, int, int]] = []
    stats = {
        "candidate_pairs_total": 0,
        "candidate_pairs_cis": 0,
        "candidate_pairs_trans": 0,
        "candidate_pairs_skipped_same_copy_same_locus": 0,
        "candidate_pairs_skipped_opposite_copy_same_locus": 0,
        "candidate_pairs_skipped_min_genomic_separation": 0,
        "candidate_pairs_skipped_training_same_bin": 0,
    }
    n = len(beads)
    for i in range(n):
        xi = np.asarray(beads[i]["xyz"])
        for j in range(i + 1, n):
            if beads[i]["copy"] == beads[j]["copy"] and beads[i]["chrom"] == beads[j]["chrom"] and beads[i]["start"] == beads[j]["start"]:
                stats["candidate_pairs_skipped_same_copy_same_locus"] += 1
                continue
            if beads[i]["chrom"] == beads[j]["chrom"] and beads[i]["start"] == beads[j]["start"]:
                stats["candidate_pairs_skipped_opposite_copy_same_locus"] += 1
                continue
            if (
                exclude_training_same_bin
                and beads[i]["chrom"] == beads[j]["chrom"]
                and int(beads[i]["start"]) // training_bin_size == int(beads[j]["start"]) // training_bin_size
            ):
                stats["candidate_pairs_skipped_training_same_bin"] += 1
                continue
            if beads[i]["chrom"] == beads[j]["chrom"] and abs(int(beads[i]["start"]) - int(beads[j]["start"])) < min_genomic_separation:
                stats["candidate_pairs_skipped_min_genomic_separation"] += 1
                continue
            d = float(np.linalg.norm(xi - np.asarray(beads[j]["xyz"])))
            if len(heap) < max_keep:
                heapq.heappush(heap, (-d, i, j))
            elif d < -heap[0][0]:
                heapq.heapreplace(heap, (-d, i, j))
    candidates = [(i, j, -neg_d) for neg_d, i, j in heap]
    stats["candidate_pairs_total"] = len(candidates)
    stats["candidate_pairs_cis"] = sum(1 for i, j, _ in candidates if beads[i]["chrom"] == beads[j]["chrom"])
    stats["candidate_pairs_trans"] = len(candidates) - stats["candidate_pairs_cis"]
    return candidates, stats


def canonical_pair(beads: list[dict[str, object]], i: int, j: int) -> tuple[str, int, str, int]:
    a = beads[i]
    b = beads[j]
    key_a = (chrom_sort_key(str(a["chrom"])), int(a["start"]))
    key_b = (chrom_sort_key(str(b["chrom"])), int(b["start"]))
    if key_a <= key_b:
        return str(a["chrom"]), int(a["start"]), str(b["chrom"]), int(b["start"])
    return str(b["chrom"]), int(b["start"]), str(a["chrom"]), int(a["start"])


def sample_pairs(
    candidates: list[tuple[int, int, float]],
    target: int,
    seed: int,
    replacement: bool,
) -> list[tuple[int, int, float]]:
    rng = random.Random(seed)
    if not candidates or target <= 0:
        return []
    if replacement:
        return [candidates[rng.randrange(len(candidates))] for _ in range(target)]
    target = min(target, len(candidates))
    return rng.sample(candidates, target)


def write_pairs(
    path: Path,
    beads: list[dict[str, object]],
    sampled: Iterable[tuple[int, int, float]],
    chrom_sizes: dict[str, int],
) -> dict[str, int]:
    rows: list[tuple[str, int, str, int]] = [canonical_pair(beads, i, j) for i, j, _ in sampled]
    rows.sort(key=lambda r: (chrom_sort_key(r[0]), chrom_sort_key(r[2]), r[1], r[3]))
    path.parent.mkdir(parents=True, exist_ok=True)
    total = cis = trans = 0
    with gzip.open(path, "wt") as out:
        out.write("## pairs format v1.0\n")
        out.write("#sorted: chr1-chr2-pos1-pos2\n")
        out.write("#shape: upper triangle\n")
        for chrom in sorted(chrom_sizes, key=chrom_sort_key):
            out.write(f"#chromosome: {chrom} {chrom_sizes[chrom]}\n")
        out.write("#columns:\treadID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2\tphase0\tphase1\n")
        for idx, (chr1, pos1, chr2, pos2) in enumerate(rows):
            total += 1
            if chr1 == chr2:
                cis += 1
            else:
                trans += 1
            out.write(f"synthetic_charm3dg_{idx:012d}\t{chr1}\t{pos1}\t{chr2}\t{pos2}\t+\t+\t.\t.\n")
    return {"output_pairs_total": total, "output_pairs_cis": cis, "output_pairs_trans": trans}


def write_metadata(path: Path, values: dict[str, object]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        for key in sorted(values):
            writer.writerow([key, values[key]])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tdg", type=Path, required=True)
    parser.add_argument("--observed-pairs", type=Path, required=True)
    parser.add_argument("--out-pairs", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument("--training-bin-size", type=int, default=1_000_000)
    parser.add_argument("--tdg-resolution", type=int, default=20_000)
    parser.add_argument("--radius-pr", type=float, required=True)
    parser.add_argument("--particle-radius", type=float, default=0.0)
    parser.add_argument("--target-mode", choices=["observed_non_same_bin", "observed_total", "fixed"], default="observed_non_same_bin")
    parser.add_argument("--target-count", type=int, default=0)
    parser.add_argument("--candidate-multiplier", type=float, default=4.0)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--min-genomic-separation", type=int, default=0)
    parser.add_argument("--include-training-same-bin", action="store_true")
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--nearest-if-insufficient", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.radius_pr <= 0:
        raise SystemExit("--radius-pr must be > 0")
    if args.training_bin_size <= 0:
        raise SystemExit("--training-bin-size must be > 0")
    beads = read_3dg(args.tdg)
    observed = count_pairs(args.observed_pairs, args.training_bin_size)
    particle_radius = args.particle_radius if args.particle_radius > 0 else estimate_particle_radius(beads)
    radius = particle_radius * args.radius_pr
    if args.target_mode == "observed_total":
        target = observed["observed_pairs_total"]
    elif args.target_mode == "observed_non_same_bin":
        target = observed["observed_pairs_non_same_bin_at_training_resolution"]
    else:
        if args.target_count <= 0:
            raise SystemExit("--target-count must be > 0 for --target-mode fixed")
        target = args.target_count
    max_keep = max(target, int(math.ceil(target * args.candidate_multiplier)))
    candidates, candidate_stats = reservoir_candidates(
        beads,
        radius=radius,
        max_keep=max_keep,
        seed=args.seed,
        min_genomic_separation=args.min_genomic_separation,
        training_bin_size=args.training_bin_size,
        exclude_training_same_bin=not args.include_training_same_bin,
    )
    fallback_nearest = 0
    if len(candidates) < target and args.nearest_if_insufficient:
        candidates, candidate_stats = nearest_candidates(
            beads,
            target,
            args.min_genomic_separation,
            args.training_bin_size,
            not args.include_training_same_bin,
        )
        fallback_nearest = 1
    if len(candidates) < target and not args.replace:
        raise SystemExit(
            f"insufficient candidates after filtering: {len(candidates)} < target {target}; "
            "increase radius, lower target, set --replace, or explicitly enable --nearest-if-insufficient"
        )
    replacement = args.replace or len(candidates) < target
    sampled = sample_pairs(candidates, target, args.seed + 1, replacement)
    chrom_sizes = chrom_sizes_from_beads(beads, args.tdg_resolution)
    output_stats = write_pairs(args.out_pairs, beads, sampled, chrom_sizes)
    values: dict[str, object] = {
        "tdg": str(args.tdg),
        "observed_pairs": str(args.observed_pairs),
        "out_pairs": str(args.out_pairs),
        "training_bin_size": args.training_bin_size,
        "tdg_resolution": args.tdg_resolution,
        "radius_pr": args.radius_pr,
        "particle_radius": particle_radius,
        "radius_distance": radius,
        "target_mode": args.target_mode,
        "target_count_requested": target,
        "candidate_multiplier": args.candidate_multiplier,
        "seed": args.seed,
        "min_genomic_separation": args.min_genomic_separation,
        "exclude_training_same_bin": int(not args.include_training_same_bin),
        "sampling_with_replacement": int(replacement),
        "nearest_if_insufficient": int(args.nearest_if_insufficient),
        "fallback_nearest_used": fallback_nearest,
        "n_beads": len(beads),
        "n_chromosomes": len(chrom_sizes),
    }
    values.update(observed)
    values.update(candidate_stats)
    values.update(output_stats)
    write_metadata(args.metadata, values)
    if args.json:
        args.json.write_text(json.dumps(values, indent=2, sort_keys=True) + "\n")
    print(f"wrote\t{args.out_pairs}")
    print(f"wrote\t{args.metadata}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
