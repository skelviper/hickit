#!/usr/bin/env python3
"""Decode read-group-consistent copy assignments from contacts.seg.gz.

This is a blind post-training diagnostic. It uses a trained posterior and the
read-level grouping retained in contacts.seg.gz to enforce that all segment
pairs from the same read are compatible with a single binary copy assignment
per segment. SNP phase and CHARM/3DG are used only after decoding for standard
evaluation/scoring.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import itertools
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, TextIO

import numpy as np


def load_eval_module() -> Any:
    repo_root = Path(__file__).resolve().parents[1]
    eval_path = repo_root / "eval" / "evaluate_p9016_baseline.py"
    spec = importlib.util.spec_from_file_location("evaluate_p9016_baseline", eval_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load eval module from {eval_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def parse_segment(field: str) -> tuple[str, int, int | None] | None:
    parts = field.split("!")
    if len(parts) < 7:
        return None
    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError:
        return None
    try:
        mapq = int(parts[5])
    except ValueError:
        mapq = None
    return parts[0], (start + end) // 2, mapq


def start_bin(pos: int, bin_size: int) -> int:
    return (pos // bin_size) * bin_size


def posterior_key_for_segments(eval_mod: Any, a: tuple[str, int, int | None], b: tuple[str, int, int | None], bin_size: int) -> tuple[tuple[str, int, str, int], bool]:
    chrom1, pos1, _ = a
    chrom2, pos2, _ = b
    s1 = start_bin(pos1, bin_size)
    s2 = start_bin(pos2, bin_size)
    key = eval_mod.canonical_bpair_key(chrom1, s1, chrom2, s2)
    swapped = key != (chrom1, s1, chrom2, s2)
    return key, swapped


def state_from_assignment(copy_a: int, copy_b: int, swapped: bool) -> int:
    if swapped:
        copy_a, copy_b = copy_b, copy_a
    return copy_a * 2 + copy_b


def decode_read(
    eval_mod: Any,
    segments: list[tuple[str, int, int | None]],
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    bin_size: int,
    eps: float,
) -> tuple[list[tuple[tuple[str, int, str, int], int]], int, float] | None:
    edges: list[tuple[int, int, tuple[str, int, str, int], bool, np.ndarray]] = []
    for i, j in itertools.combinations(range(len(segments)), 2):
        key, swapped = posterior_key_for_segments(eval_mod, segments[i], segments[j], bin_size)
        item = posterior.get(key)
        if item is None:
            continue
        if key[0] == key[2] and key[1] == key[3]:
            continue
        p4 = np.asarray(item["p4"], dtype=float)
        if not np.all(np.isfinite(p4)):
            continue
        edges.append((i, j, key, swapped, np.maximum(p4, eps)))
    if not edges:
        return None
    best_score = float("-inf")
    best_mask = 0
    n = len(segments)
    for mask in range(1 << n):
        score = 0.0
        for i, j, _, swapped, p4 in edges:
            ci = (mask >> i) & 1
            cj = (mask >> j) & 1
            state = state_from_assignment(ci, cj, swapped)
            score += math.log(float(p4[state]))
        if score > best_score:
            best_score = score
            best_mask = mask
    decoded: list[tuple[tuple[str, int, str, int], int]] = []
    for i, j, key, swapped, _ in edges:
        ci = (best_mask >> i) & 1
        cj = (best_mask >> j) & 1
        decoded.append((key, state_from_assignment(ci, cj, swapped)))
    return decoded, len(edges), best_score


def read_selected_swaps(path: Path) -> dict[str, int]:
    swaps: dict[str, int] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if int(row.get("selected_for_eval", "0")) == 1:
                swaps[row["chrom"]] = int(row["copy_swap"])
    return swaps


def empty_stats() -> dict[str, float]:
    return {"contacts": 0.0, "correct": 0.0, "same_cross": 0.0, "decoded_contacts": 0.0}


def update_stats(eval_mod: Any, stats: dict[str, float], counts: np.ndarray, p4: np.ndarray, decoded: bool) -> None:
    n = int(counts.sum())
    if n <= 0:
        return
    pred = int(np.argmax(p4))
    stats["contacts"] += n
    stats["correct"] += int(counts[pred])
    pred_same = pred in (0, 3)
    same_counts = int(counts[0] + counts[3])
    cross_counts = int(counts[1] + counts[2])
    stats["same_cross"] += same_counts if pred_same else cross_counts
    if decoded:
        stats["decoded_contacts"] += n


def finalize(scope: str, policy: str, stats: dict[str, float], total_contacts: float) -> dict[str, object]:
    contacts = stats["contacts"]
    return {
        "policy": policy,
        "scope": scope,
        "n_eval_contacts": int(contacts),
        "top1_correct_contacts": int(stats["correct"]),
        "top1_accuracy": stats["correct"] / contacts if contacts else float("nan"),
        "same_cross_accuracy": stats["same_cross"] / contacts if contacts else float("nan"),
        "decoded_contact_fraction": stats["decoded_contacts"] / contacts if contacts else float("nan"),
        "scope_contact_fraction_vs_baseline": contacts / total_contacts if total_contacts else float("nan"),
        "construction_uses_phase_labels": 0,
        "construction_uses_charm_or_reference": 0,
        "eval_only_uses_phase_labels_for_scoring": 1,
    }


def evaluate_policies(
    eval_mod: Any,
    contact_counts: dict[tuple[str, int, str, int], np.ndarray],
    posterior: dict[tuple[str, int, str, int], dict[str, object]],
    decoded_counts: dict[tuple[str, int, str, int], np.ndarray],
    reference: dict[tuple[str, int, int], np.ndarray],
    swaps: dict[str, int],
    blend_lambdas: list[float],
) -> list[dict[str, object]]:
    stats_by_policy_scope: dict[tuple[str, str], dict[str, float]] = {}
    total_by_scope: Counter[str] = Counter()
    policies = ["baseline", "readgroup_only"] + [f"hybrid_lambda{lam:g}" for lam in blend_lambdas]
    for key, counts in contact_counts.items():
        item = posterior.get(key)
        if item is None:
            continue
        if eval_mod.p4_from_coords_for_bpair(reference, item) is None:
            continue
        chrom1 = str(item["chrom1"])
        chrom2 = str(item["chrom2"])
        scope = "genome_cis" if chrom1 == chrom2 else "genome_trans"
        scopes = ["genome_all", scope]
        base_p4 = np.asarray(item["p4"], dtype=float)
        decoded = decoded_counts.get(key)
        decoded_p4 = None
        if decoded is not None and int(decoded.sum()) > 0:
            decoded_p4 = decoded.astype(float) / float(decoded.sum())
        for s in scopes:
            total_by_scope[s] += int(counts.sum())
        for policy in policies:
            if policy == "baseline":
                p4 = base_p4
                has_decoded = False
            elif policy == "readgroup_only":
                if decoded_p4 is None:
                    continue
                p4 = decoded_p4
                has_decoded = True
            else:
                lam = float(policy.removeprefix("hybrid_lambda"))
                if decoded_p4 is None:
                    p4 = base_p4
                    has_decoded = False
                else:
                    p4 = (1.0 - lam) * base_p4 + lam * decoded_p4
                    p4 = p4 / p4.sum()
                    has_decoded = True
            aligned = eval_mod.align_p4_to_truth_gauge(p4, swaps.get(chrom1, 0), swaps.get(chrom2, 0))
            for s in scopes:
                stats = stats_by_policy_scope.setdefault((policy, s), empty_stats())
                update_stats(eval_mod, stats, counts, aligned, has_decoded)
    rows: list[dict[str, object]] = []
    for (policy, scope), stats in sorted(stats_by_policy_scope.items()):
        rows.append(finalize(scope, policy, stats, float(total_by_scope[scope])))
    return rows


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.9g}"
    return str(value)


def write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in fieldnames})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seg", required=True, type=Path)
    ap.add_argument("--config-output-dir", required=True, type=Path)
    ap.add_argument("--eval-output-dir", required=True, type=Path)
    ap.add_argument("--pairs", required=True, type=Path)
    ap.add_argument("--reference-3dg", required=True, type=Path)
    ap.add_argument("--bin-size", type=int, default=1_000_000)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--max-read-segments", type=int, default=8)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()

    eval_mod = load_eval_module()
    posterior = eval_mod.read_posterior(args.config_output_dir / "p9016_full.bpair_posterior.tsv")
    reference = eval_mod.read_reference_3dg(args.reference_3dg, args.bin_size)
    contact_counts = eval_mod.read_contact_truth_counts(args.pairs, args.bin_size, posterior)
    swaps = read_selected_swaps(args.eval_output_dir / "whole_chrom_snp_oracle_swaps.tsv")

    decoded_counts: dict[tuple[str, int, str, int], np.ndarray] = defaultdict(lambda: np.zeros(4, dtype=np.int64))
    diag = Counter()
    with open_text(args.seg) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            diag["reads_seen"] += 1
            segments = [parse_segment(field) for field in fields[1:]]
            usable = [seg for seg in segments if seg is not None and (seg[2] is None or seg[2] >= args.min_mapq)]
            if len(usable) < 2:
                continue
            if len(usable) > args.max_read_segments:
                diag["reads_skipped_too_many_segments"] += 1
                continue
            decoded = decode_read(eval_mod, usable, posterior, args.bin_size, 1e-9)
            if decoded is None:
                diag["reads_without_posterior_edges"] += 1
                continue
            pairs, n_edges, _ = decoded
            diag["reads_decoded"] += 1
            diag["posterior_edges_decoded"] += n_edges
            if len(usable) >= 3:
                diag["multiseg_reads_decoded"] += 1
            if len({seg[0] for seg in usable}) >= 2:
                diag["multichrom_reads_decoded"] += 1
            for key, state in pairs:
                decoded_counts[key][state] += 1

    rows = evaluate_policies(
        eval_mod,
        contact_counts,
        posterior,
        decoded_counts,
        reference,
        swaps,
        [0.25, 0.5, 0.75, 1.0],
    )
    fields = [
        "policy",
        "scope",
        "n_eval_contacts",
        "top1_correct_contacts",
        "top1_accuracy",
        "same_cross_accuracy",
        "decoded_contact_fraction",
        "scope_contact_fraction_vs_baseline",
        "construction_uses_phase_labels",
        "construction_uses_charm_or_reference",
        "eval_only_uses_phase_labels_for_scoring",
    ]
    write_rows(args.outdir / "readgroup_decoding_accuracy.tsv", rows, fields)

    args.outdir.mkdir(parents=True, exist_ok=True)
    with (args.outdir / "readgroup_decoding_diag.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        for key in sorted(diag):
            writer.writerow([key, diag[key]])
        writer.writerow(["decoded_bpairs", len(decoded_counts)])
        writer.writerow(["decoded_trans_bpairs", sum(1 for key in decoded_counts if key[0] != key[2])])
        writer.writerow(["construction_uses_phase_labels", 0])
        writer.writerow(["construction_uses_charm_or_reference", 0])
        writer.writerow(["eval_only_uses_phase_labels_for_scoring", 1])


if __name__ == "__main__":
    main()
