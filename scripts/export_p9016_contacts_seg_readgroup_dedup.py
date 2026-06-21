#!/usr/bin/env python3
"""Export blind-safe contacts.seg-derived pairs with hickit-style dedup metadata.

This script intentionally separates two concepts:

* hickit-like reference pairs: adjacent segment contacts, hickit boundary
  coordinates, hickit-like close-leg filtering, and hickit-style pair-level
  duplicate representatives.
* readgroup-preserving pairs: contacts carry read_group_id/read_seg metadata so
  downstream readgroup-aware blind models can consume the multi-segment source.
* readchain-deduplicated pairs: whole segment chains are deduplicated across
  reads before contact generation, then representative chains emit adjacent or
  all-combination contacts with readgroup metadata.

Phase labels from contacts.seg are never written to output pairs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


def open_text(path: Path, mode: str = "rt") -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def bucket_count(n: int) -> str:
    if n < 10:
        return str(n)
    if n < 100:
        lo = (n // 10) * 10
        return f"{lo}-{lo + 9}"
    if n < 1000:
        lo = (n // 100) * 100
        return f"{lo}-{lo + 99}"
    return "1000+"


def fmt_float(x: float | int | None) -> str:
    if x is None:
        return "NA"
    try:
        y = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not math.isfinite(y):
        return "NA"
    return f"{y:.9g}"


def percentile(values: list[int], q: float) -> float | None:
    if not values:
        return None
    vals = sorted(values)
    if len(vals) == 1:
        return float(vals[0])
    idx = q * (len(vals) - 1)
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    if lo == hi:
        return float(vals[lo])
    return float(vals[lo]) * (hi - idx) + float(vals[hi]) * (idx - lo)


@dataclass(frozen=True)
class Segment:
    chrom: str
    start: int
    end: int
    strand: str
    phase: str
    mapq: int
    count: str
    mate: str
    raw_index: int

    @property
    def midpoint(self) -> int:
        return (self.start + self.end) // 2


@dataclass
class PairRecord:
    read_id: str
    read_group_id: int
    read_seg0: int
    read_seg1: int
    read_n_segments: int
    read_n_chroms: int
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    strand1: str
    strand2: str
    source_mode: str
    cluster_id: int = -1
    cluster_size: int = 1
    is_representative: int = 1

    def sort_key(self, chrom_order: dict[str, int]) -> tuple[int, int, int, int, str, int, int, int]:
        return (
            chrom_order.get(self.chrom1, 10**9),
            chrom_order.get(self.chrom2, 10**9),
            self.pos1,
            self.pos2,
            self.read_id,
            self.read_group_id,
            self.read_seg0,
            self.read_seg1,
        )

    def duplicate_scan_key(self, chrom_order: dict[str, int]) -> tuple[int, int, int, int, str, str]:
        return (
            chrom_order.get(self.chrom1, 10**9),
            chrom_order.get(self.chrom2, 10**9),
            self.pos1,
            self.pos2,
            self.strand1,
            self.strand2,
	        )


@dataclass
class ReadChainRecord:
    read_id: str
    raw_segments: list[Segment]
    read_group_id: int = -1
    cluster_id: int = -1
    cluster_size: int = 1
    is_representative: int = 1

    @property
    def n_chroms(self) -> int:
        return len({seg.chrom for seg in self.raw_segments})

    def sort_key(self, chrom_order: dict[str, int]) -> tuple[int, int, int, int, str]:
        first = self.raw_segments[0]
        return (
            len(self.raw_segments),
            chrom_order.get(first.chrom, 10**9),
            first.start,
            first.end,
            self.read_id,
        )


def parse_segment(field: str, raw_index: int) -> Segment | None:
    parts = field.split("!")
    if len(parts) < 6:
        return None
    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError:
        return None
    if end < start:
        end = start
    strand = parts[3] if len(parts) > 3 and parts[3] in {"+", "-"} else "."
    phase = parts[4] if len(parts) > 4 else "."
    try:
        mapq = int(parts[5])
    except ValueError:
        mapq = 0
    count = parts[6] if len(parts) > 6 else "."
    mate = parts[7] if len(parts) > 7 else "."
    return Segment(parts[0], start, end, strand, phase, mapq, count, mate, raw_index)


def ordered_pair(
    read_id: str,
    read_group_id: int,
    seg_a: Segment,
    seg_b: Segment,
    read_n_segments: int,
    read_n_chroms: int,
    chrom_order: dict[str, int],
    source_mode: str,
    coordinate_mode: str,
) -> PairRecord:
    if coordinate_mode == "boundary":
        # Match hickit hk_seg2pair(): previous segment contributes end and next
        # segment contributes start before canonical ordering.
        raw_a_pos = seg_a.end
        raw_b_pos = seg_b.start
    elif coordinate_mode == "midpoint":
        raw_a_pos = seg_a.midpoint
        raw_b_pos = seg_b.midpoint
    else:
        raise ValueError(f"unknown coordinate_mode: {coordinate_mode}")
    ak = (chrom_order.get(seg_a.chrom, 10**9), raw_a_pos, seg_a.chrom)
    bk = (chrom_order.get(seg_b.chrom, 10**9), raw_b_pos, seg_b.chrom)
    if ak <= bk:
        return PairRecord(
            read_id=read_id,
            read_group_id=read_group_id,
            read_seg0=seg_a.raw_index,
            read_seg1=seg_b.raw_index,
            read_n_segments=read_n_segments,
            read_n_chroms=read_n_chroms,
            chrom1=seg_a.chrom,
            pos1=raw_a_pos,
            chrom2=seg_b.chrom,
            pos2=raw_b_pos,
            strand1=seg_a.strand,
            strand2=seg_b.strand,
            source_mode=source_mode,
        )
    return PairRecord(
        read_id=read_id,
        read_group_id=read_group_id,
        read_seg0=seg_b.raw_index,
        read_seg1=seg_a.raw_index,
        read_n_segments=read_n_segments,
        read_n_chroms=read_n_chroms,
        chrom1=seg_b.chrom,
        pos1=raw_b_pos,
        chrom2=seg_a.chrom,
        pos2=raw_a_pos,
        strand1=seg_b.strand,
        strand2=seg_a.strand,
        source_mode=source_mode,
    )


def close_leg_keep(pair: PairRecord, min_leg_dist: int, all_close_leg: bool) -> bool:
    if pair.chrom1 != pair.chrom2:
        return True
    if pair.pos2 - pair.pos1 >= min_leg_dist:
        return True
    same_or_unknown = pair.strand1 == pair.strand2 or pair.strand1 == "." or pair.strand2 == "."
    if same_or_unknown or all_close_leg:
        return False
    return True


def iter_pair_indices(n_segments: int, mode: str) -> Iterable[tuple[int, int]]:
    if mode == "adjacent":
        for i in range(1, n_segments):
            yield i - 1, i
    elif mode == "allcomb":
        yield from itertools.combinations(range(n_segments), 2)
    else:
        raise ValueError(f"unknown pair mode: {mode}")


def build_pairs(
    seg_path: Path,
    source_label: str,
    mode: str,
    chrom_order: dict[str, int],
    min_mapq: int,
    min_leg_dist: int,
    max_seg: int,
    all_close_leg: bool,
    coordinate_mode: str,
    keep_over_max_seg_for_readgroup: bool,
) -> tuple[list[PairRecord], dict[str, object], list[str], Counter[str]]:
    chrom_lines: list[str] = []
    pairs: list[PairRecord] = []
    read_segment_hist: Counter[str] = Counter()
    n_lines = 0
    n_readgroups_total = 0
    n_readgroups_ge2_raw = 0
    n_readgroups_ge2_mapq = 0
    n_readgroups_used = 0
    n_readgroups_skipped_maxseg = 0
    n_readgroups_multichrom_raw = 0
    n_readgroups_multichrom_used = 0
    n_readgroups_phase_labeled = 0
    n_segments_raw = 0
    n_segments_mapq = 0
    n_pairs_before_close = 0
    n_pairs_after_close = 0
    n_pairs_cis = 0
    n_pairs_trans = 0
    with open_text(seg_path) as fh:
        for line in fh:
            if line.startswith("#chromosome:") or line.startswith("#chromsize:"):
                chrom_lines.append(line.rstrip("\n"))
                parts = line.strip().split()
                if len(parts) >= 3 and parts[1] not in chrom_order:
                    chrom_order[parts[1]] = len(chrom_order)
                continue
            if not line.strip() or line.startswith("#"):
                continue
            n_lines += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            read_id = fields[0]
            n_readgroups_total += 1
            raw_segments = [parse_segment(field, idx) for idx, field in enumerate(fields[1:])]
            raw_segments = [seg for seg in raw_segments if seg is not None]
            n_segments_raw += len(raw_segments)
            read_segment_hist[bucket_count(len(raw_segments))] += 1
            if len(raw_segments) >= 2:
                n_readgroups_ge2_raw += 1
            if len({seg.chrom for seg in raw_segments}) > 1:
                n_readgroups_multichrom_raw += 1
            if any(seg.phase in {"0", "1"} for seg in raw_segments):
                n_readgroups_phase_labeled += 1
            if max_seg > 0 and len(raw_segments) > max_seg and not keep_over_max_seg_for_readgroup:
                n_readgroups_skipped_maxseg += 1
                continue
            n_segments_mapq += sum(1 for seg in raw_segments if seg.mapq >= min_mapq)
            if sum(1 for seg in raw_segments if seg.mapq >= min_mapq) < 2:
                continue
            n_readgroups_ge2_mapq += 1
            n_readgroups_used += 1
            read_group_id = n_readgroups_used - 1
            n_chroms = len({seg.chrom for seg in raw_segments})
            if n_chroms > 1:
                n_readgroups_multichrom_used += 1
            for ia, ib in iter_pair_indices(len(raw_segments), mode):
                if raw_segments[ia].mapq < min_mapq or raw_segments[ib].mapq < min_mapq:
                    continue
                pair = ordered_pair(
                    read_id,
                    read_group_id,
                    raw_segments[ia],
                    raw_segments[ib],
                    len(raw_segments),
                    n_chroms,
                    chrom_order,
                    mode,
                    coordinate_mode,
                )
                n_pairs_before_close += 1
                if not close_leg_keep(pair, min_leg_dist, all_close_leg):
                    continue
                pairs.append(pair)
                n_pairs_after_close += 1
                if pair.chrom1 == pair.chrom2:
                    n_pairs_cis += 1
                else:
                    n_pairs_trans += 1
    summary: dict[str, object] = {
        "source_label": source_label,
        "seg_path": str(seg_path),
        "export_mode": mode,
        "coordinate_mode": coordinate_mode,
        "min_mapq": min_mapq,
        "min_leg_dist": min_leg_dist,
        "max_seg": max_seg,
        "keep_over_max_seg_for_readgroup": int(keep_over_max_seg_for_readgroup),
        "all_close_leg": int(all_close_leg),
        "n_lines": n_lines,
        "n_readgroups_total": n_readgroups_total,
        "n_readgroups_ge2_raw": n_readgroups_ge2_raw,
        "n_readgroups_ge2_mapq": n_readgroups_ge2_mapq,
        "n_readgroups_used": n_readgroups_used,
        "n_readgroups_skipped_maxseg": n_readgroups_skipped_maxseg,
        "n_readgroups_multichrom_raw": n_readgroups_multichrom_raw,
        "n_readgroups_multichrom_used": n_readgroups_multichrom_used,
        "n_readgroups_phase_labeled": n_readgroups_phase_labeled,
        "n_segments_raw": n_segments_raw,
        "n_segments_mapq": n_segments_mapq,
        "n_pairs_before_close_filter": n_pairs_before_close,
        "n_pairs_after_close_filter": n_pairs_after_close,
        "n_pairs_cis": n_pairs_cis,
        "n_pairs_trans": n_pairs_trans,
        "phase_labels_written": 0,
        "training_uses_phase_labels": 0,
        "training_uses_charm_or_reference": 0,
    }
    return pairs, summary, chrom_lines, read_segment_hist


def same_readchain(a: ReadChainRecord, b: ReadChainRecord, dup_dist: int) -> bool:
    if len(a.raw_segments) != len(b.raw_segments):
        return False
    for sa, sb in zip(a.raw_segments, b.raw_segments):
        if sa.chrom != sb.chrom or sa.strand != sb.strand:
            return False
        if abs(sa.start - sb.start) >= dup_dist:
            return False
        if abs(sa.end - sb.end) >= dup_dist:
            return False
    return True


def annotate_readchain_clusters(
    records: list[ReadChainRecord],
    chrom_order: dict[str, int],
    dup_dist: int,
) -> dict[str, object]:
    records.sort(key=lambda r: r.sort_key(chrom_order))
    representatives: list[ReadChainRecord] = []
    cluster_members: dict[int, list[int]] = {}
    n_duplicates = 0
    if dup_dist <= 0:
        for idx, record in enumerate(records):
            record.cluster_id = idx
            record.cluster_size = 1
            record.is_representative = 1
        return {
            "n_readchains_before_dedup": len(records),
            "n_readchains_after_dedup": len(records),
            "n_readchains_removed_by_dedup": 0,
            "readchain_dedup_rate": 0.0,
            "n_readchain_dup_clusters": len(records),
            "readchain_cluster_size_p50": 1 if records else None,
            "readchain_cluster_size_p90": 1 if records else None,
            "readchain_cluster_size_p99": 1 if records else None,
            "readchain_cluster_size_max": 1 if records else None,
        }
    for idx, record in enumerate(records):
        first = record.raw_segments[0]
        first_chr = chrom_order.get(first.chrom, 10**9)
        matched_cluster = -1
        for prev in reversed(representatives):
            pfirst = prev.raw_segments[0]
            if len(prev.raw_segments) != len(record.raw_segments):
                if len(prev.raw_segments) < len(record.raw_segments):
                    break
                continue
            pfirst_chr = chrom_order.get(pfirst.chrom, 10**9)
            if pfirst_chr != first_chr:
                if pfirst_chr < first_chr:
                    break
                continue
            if first.start - pfirst.start >= dup_dist:
                break
            if same_readchain(prev, record, dup_dist):
                matched_cluster = prev.cluster_id
                break
        if matched_cluster >= 0:
            record.cluster_id = matched_cluster
            record.is_representative = 0
            cluster_members[matched_cluster].append(idx)
            n_duplicates += 1
        else:
            record.cluster_id = len(representatives)
            record.is_representative = 1
            representatives.append(record)
            cluster_members[record.cluster_id] = [idx]
    sizes = [len(v) for v in cluster_members.values()]
    for members in cluster_members.values():
        size = len(members)
        for idx in members:
            records[idx].cluster_size = size
    return {
        "n_readchains_before_dedup": len(records),
        "n_readchains_after_dedup": len(representatives),
        "n_readchains_removed_by_dedup": n_duplicates,
        "readchain_dedup_rate": n_duplicates / len(records) if records else 0.0,
        "n_readchain_dup_clusters": len(representatives),
        "readchain_cluster_size_p50": percentile(sizes, 0.50),
        "readchain_cluster_size_p90": percentile(sizes, 0.90),
        "readchain_cluster_size_p99": percentile(sizes, 0.99),
        "readchain_cluster_size_max": max(sizes) if sizes else None,
    }


def collect_readchain_representatives(
    seg_path: Path,
    chrom_order: dict[str, int],
    min_mapq: int,
    max_seg: int,
    dup_dist: int,
) -> tuple[list[ReadChainRecord], dict[str, object], list[str], Counter[str]]:
    chrom_lines: list[str] = []
    read_segment_hist: Counter[str] = Counter()
    records: list[ReadChainRecord] = []
    n_lines = 0
    n_readgroups_total = 0
    n_readgroups_ge2_raw = 0
    n_readgroups_ge2_mapq = 0
    n_readgroups_used = 0
    n_readgroups_skipped_maxseg = 0
    n_readgroups_multichrom_raw = 0
    n_readgroups_multichrom_used = 0
    n_readgroups_phase_labeled = 0
    n_segments_raw = 0
    n_segments_mapq = 0
    with open_text(seg_path) as fh:
        for line in fh:
            if line.startswith("#chromosome:") or line.startswith("#chromsize:"):
                chrom_lines.append(line.rstrip("\n"))
                parts = line.strip().split()
                if len(parts) >= 3 and parts[1] not in chrom_order:
                    chrom_order[parts[1]] = len(chrom_order)
                continue
            if not line.strip() or line.startswith("#"):
                continue
            n_lines += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            read_id = fields[0]
            n_readgroups_total += 1
            raw_segments = [parse_segment(field, idx) for idx, field in enumerate(fields[1:])]
            raw_segments = [seg for seg in raw_segments if seg is not None]
            n_segments_raw += len(raw_segments)
            read_segment_hist[bucket_count(len(raw_segments))] += 1
            if len(raw_segments) >= 2:
                n_readgroups_ge2_raw += 1
            if len({seg.chrom for seg in raw_segments}) > 1:
                n_readgroups_multichrom_raw += 1
            if any(seg.phase in {"0", "1"} for seg in raw_segments):
                n_readgroups_phase_labeled += 1
            usable_segments = [
                Segment(
                    seg.chrom,
                    seg.start,
                    seg.end,
                    seg.strand,
                    seg.phase,
                    seg.mapq,
                    seg.count,
                    seg.mate,
                    usable_idx,
                )
                for usable_idx, seg in enumerate(seg for seg in raw_segments if seg.mapq >= min_mapq)
            ]
            n_mapq = len(usable_segments)
            n_segments_mapq += n_mapq
            if n_mapq < 2:
                continue
            n_readgroups_ge2_mapq += 1
            if max_seg > 0 and n_mapq > max_seg:
                n_readgroups_skipped_maxseg += 1
                continue
            n_readgroups_used += 1
            if len({seg.chrom for seg in usable_segments}) > 1:
                n_readgroups_multichrom_used += 1
            records.append(ReadChainRecord(read_id=read_id, raw_segments=usable_segments))

    readchain_summary = annotate_readchain_clusters(records, chrom_order, dup_dist)
    representatives = [record for record in records if record.is_representative]
    representatives.sort(key=lambda r: r.sort_key(chrom_order))
    n_representative_multichrom = sum(1 for record in representatives if record.n_chroms > 1)
    summary: dict[str, object] = {
        "dedup_unit": "readchain",
        "n_lines": n_lines,
        "n_readgroups_total": n_readgroups_total,
        "n_readgroups_ge2_raw": n_readgroups_ge2_raw,
        "n_readgroups_ge2_mapq": n_readgroups_ge2_mapq,
        "n_readgroups_used": n_readgroups_used,
        "n_readgroups_after_readchain_dedup": len(representatives),
        "n_readgroups_skipped_maxseg": n_readgroups_skipped_maxseg,
        "n_readgroups_multichrom_raw": n_readgroups_multichrom_raw,
        "n_readgroups_multichrom_used": n_readgroups_multichrom_used,
        "n_readgroups_multichrom_after_readchain_dedup": n_representative_multichrom,
        "n_readgroups_phase_labeled": n_readgroups_phase_labeled,
        "n_segments_raw": n_segments_raw,
        "n_segments_mapq": n_segments_mapq,
        "phase_labels_written": 0,
        "training_uses_phase_labels": 0,
        "training_uses_charm_or_reference": 0,
        **readchain_summary,
    }
    return representatives, summary, chrom_lines, read_segment_hist


def pairs_from_readchain_representatives(
    representatives: list[ReadChainRecord],
    mode: str,
    chrom_order: dict[str, int],
    min_leg_dist: int,
    all_close_leg: bool,
    coordinate_mode: str,
) -> tuple[list[PairRecord], dict[str, object]]:
    pairs: list[PairRecord] = []
    n_pairs_before_close = 0
    n_pairs_after_close = 0
    n_pairs_cis = 0
    n_pairs_trans = 0
    for read_group_id, record in enumerate(representatives):
        record.read_group_id = read_group_id
        for ia, ib in iter_pair_indices(len(record.raw_segments), mode):
            pair = ordered_pair(
                record.read_id,
                record.read_group_id,
                record.raw_segments[ia],
                record.raw_segments[ib],
                len(record.raw_segments),
                record.n_chroms,
                chrom_order,
                mode,
                coordinate_mode,
            )
            n_pairs_before_close += 1
            if not close_leg_keep(pair, min_leg_dist, all_close_leg):
                continue
            pairs.append(pair)
            n_pairs_after_close += 1
            if pair.chrom1 == pair.chrom2:
                n_pairs_cis += 1
            else:
                n_pairs_trans += 1
    summary = {
        "n_pairs_before_close_filter": n_pairs_before_close,
        "n_pairs_after_close_filter": n_pairs_after_close,
        "n_pairs_cis": n_pairs_cis,
        "n_pairs_trans": n_pairs_trans,
    }
    return pairs, summary


def build_pairs_readchain_dedup(
    seg_path: Path,
    source_label: str,
    mode: str,
    chrom_order: dict[str, int],
    min_mapq: int,
    min_leg_dist: int,
    max_seg: int,
    all_close_leg: bool,
    coordinate_mode: str,
    dup_dist: int,
) -> tuple[list[PairRecord], dict[str, object], list[str], Counter[str]]:
    representatives, readchain_summary, chrom_lines, read_segment_hist = collect_readchain_representatives(
        seg_path,
        chrom_order,
        min_mapq,
        max_seg,
        dup_dist,
    )
    pairs, pair_summary = pairs_from_readchain_representatives(
        representatives,
        mode,
        chrom_order,
        min_leg_dist,
        all_close_leg,
        coordinate_mode,
    )
    summary: dict[str, object] = {
        "source_label": source_label,
        "seg_path": str(seg_path),
        "export_mode": mode,
        "coordinate_mode": coordinate_mode,
        "dedup_unit": "readchain",
        "min_mapq": min_mapq,
        "min_leg_dist": min_leg_dist,
        "max_seg": max_seg,
        "keep_over_max_seg_for_readgroup": 0,
        "all_close_leg": int(all_close_leg),
        **readchain_summary,
        **pair_summary,
    }
    return pairs, summary, chrom_lines, read_segment_hist


def annotate_dup_clusters(pairs: list[PairRecord], chrom_order: dict[str, int], dup_dist: int) -> dict[str, object]:
    pairs.sort(key=lambda p: p.sort_key(chrom_order))
    representatives: list[PairRecord] = []
    cluster_members: dict[int, list[int]] = {}
    n_duplicates = 0
    if dup_dist <= 0:
        for idx, pair in enumerate(pairs):
            pair.cluster_id = idx
            pair.cluster_size = 1
            pair.is_representative = 1
        return {
            "dup_dist": dup_dist,
            "n_pairs_before_dedup": len(pairs),
            "n_pairs_after_dedup": len(pairs),
            "n_pairs_removed_by_dedup": 0,
            "dedup_rate": 0.0,
            "n_dup_clusters": len(pairs),
            "cluster_size_p50": 1 if pairs else None,
            "cluster_size_p90": 1 if pairs else None,
            "cluster_size_p99": 1 if pairs else None,
            "cluster_size_max": 1 if pairs else None,
        }
    for idx, pair in enumerate(pairs):
        matched_cluster = -1
        for prev in reversed(representatives):
            if prev.chrom1 != pair.chrom1 or prev.chrom2 != pair.chrom2:
                break
            if pair.pos1 - prev.pos1 >= dup_dist:
                break
            if abs(pair.pos2 - prev.pos2) < dup_dist and pair.strand1 == prev.strand1 and pair.strand2 == prev.strand2:
                matched_cluster = prev.cluster_id
                break
        if matched_cluster >= 0:
            pair.cluster_id = matched_cluster
            pair.is_representative = 0
            cluster_members[matched_cluster].append(idx)
            n_duplicates += 1
        else:
            pair.cluster_id = len(representatives)
            pair.is_representative = 1
            representatives.append(pair)
            cluster_members[pair.cluster_id] = [idx]
    sizes = [len(v) for v in cluster_members.values()]
    for members in cluster_members.values():
        size = len(members)
        for idx in members:
            pairs[idx].cluster_size = size
    return {
        "dup_dist": dup_dist,
        "n_pairs_before_dedup": len(pairs),
        "n_pairs_after_dedup": len(representatives),
        "n_pairs_removed_by_dedup": n_duplicates,
        "dedup_rate": n_duplicates / len(pairs) if pairs else 0.0,
        "n_dup_clusters": len(representatives),
        "cluster_size_p50": percentile(sizes, 0.50),
        "cluster_size_p90": percentile(sizes, 0.90),
        "cluster_size_p99": percentile(sizes, 0.99),
        "cluster_size_max": max(sizes) if sizes else None,
    }


def write_pairs(path: Path, chrom_lines: list[str], pairs: list[PairRecord], readgroup_cols: bool, representatives_only: bool) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open_text(path, "wt") as out:
        out.write("## pairs format v1.0\n")
        out.write("#sorted: chr1-chr2-pos1-pos2\n")
        out.write("#shape: upper triangle\n")
        for line in chrom_lines:
            out.write(line.rstrip("\n") + "\n")
        if readgroup_cols:
            out.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1 read_group_id read_seg0 read_seg1 read_n_segments\n")
        else:
            out.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1\n")
        for pair in pairs:
            if representatives_only and not bool(pair.is_representative):
                continue
            if readgroup_cols:
                out.write(
                    f"{pair.read_id}\t{pair.chrom1}\t{pair.pos1}\t{pair.chrom2}\t{pair.pos2}\t"
                    f"{pair.strand1}\t{pair.strand2}\t.\t.\t{pair.read_group_id}\t"
                    f"{pair.read_seg0}\t{pair.read_seg1}\t{pair.read_n_segments}\n"
                )
            else:
                out.write(
                    f".\t{pair.chrom1}\t{pair.pos1}\t{pair.chrom2}\t{pair.pos2}\t"
                    f"{pair.strand1}\t{pair.strand2}\t.\t.\n"
                )
            n += 1
    return n


def write_annotation(path: Path, pairs: list[PairRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "readID",
        "read_group_id",
        "read_seg0",
        "read_seg1",
        "read_n_segments",
        "read_n_chroms",
        "chr1",
        "pos1",
        "chr2",
        "pos2",
        "strand1",
        "strand2",
        "source_mode",
        "hickit_dup_cluster_id",
        "hickit_dup_cluster_size",
        "hickit_dup_is_representative",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        for pair in pairs:
            writer.writerow(
                {
                    "readID": pair.read_id,
                    "read_group_id": pair.read_group_id,
                    "read_seg0": pair.read_seg0,
                    "read_seg1": pair.read_seg1,
                    "read_n_segments": pair.read_n_segments,
                    "read_n_chroms": pair.read_n_chroms,
                    "chr1": pair.chrom1,
                    "pos1": pair.pos1,
                    "chr2": pair.chrom2,
                    "pos2": pair.pos2,
                    "strand1": pair.strand1,
                    "strand2": pair.strand2,
                    "source_mode": pair.source_mode,
                    "hickit_dup_cluster_id": pair.cluster_id,
                    "hickit_dup_cluster_size": pair.cluster_size,
                    "hickit_dup_is_representative": pair.is_representative,
                }
            )


def write_kv(path: Path, rows: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        for key, value in rows.items():
            writer.writerow([key, fmt_float(value) if isinstance(value, float) else ("NA" if value is None else value)])


def write_hist(path: Path, hist: Counter[str], value_name: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["bucket", value_name])
        for key in sorted(hist):
            writer.writerow([key, hist[key]])


def write_export_artifacts(
    outdir: Path,
    export_name: str,
    pairs: list[PairRecord],
    summary: dict[str, object],
    chrom_lines: list[str],
    seg_hist: Counter[str],
    chrom_order: dict[str, int],
    readgroup_cols: bool,
    representatives_only: bool,
    write_annotation_flag: bool,
    pair_mode: str,
    coordinate_mode: str,
    min_mapq: int,
    min_leg_dist: int,
    max_seg: int,
    dup_dist: int,
    dedup_unit: str,
) -> None:
    dedup_summary = annotate_dup_clusters(pairs, chrom_order=chrom_order, dup_dist=dup_dist)
    out_pairs = outdir / f"{export_name}.pairs.gz"
    n_written = write_pairs(out_pairs, chrom_lines, pairs, readgroup_cols, representatives_only)
    metadata = {
        **summary,
        **dedup_summary,
        "export_name": export_name,
        "out_pairs": str(out_pairs),
        "readgroup_cols_written": int(readgroup_cols),
        "representatives_only": int(representatives_only),
        "dedup_unit": dedup_unit,
        "n_pairs_written": n_written,
        "annotation_written": int(write_annotation_flag),
        "hickit_like_semantics": int(
            pair_mode == "adjacent"
            and coordinate_mode == "boundary"
            and min_mapq == 20
            and min_leg_dist == 1000
            and max_seg == 3
            and dup_dist == 100
            and representatives_only
            and not readgroup_cols
        ),
    }
    write_kv(outdir / f"{export_name}.metadata.tsv", metadata)
    write_hist(outdir / f"{export_name}.read_segment_count_hist.tsv", seg_hist, "n_readgroups")
    if write_annotation_flag:
        write_annotation(outdir / f"{export_name}.dedup_annotation.tsv", pairs)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seg", required=True, type=Path)
    ap.add_argument("--source-label", required=True)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--export-name", required=True)
    ap.add_argument("--pair-mode", choices=["adjacent", "allcomb", "both"], required=True)
    ap.add_argument("--dedup-unit", choices=["pair", "readchain"], default="pair")
    ap.add_argument("--coordinate-mode", choices=["boundary", "midpoint"], default="boundary")
    ap.add_argument("--readgroup-cols", action="store_true")
    ap.add_argument("--representatives-only", action="store_true")
    ap.add_argument("--write-annotation", action="store_true")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-leg-dist", type=int, default=1000)
    ap.add_argument("--max-seg", type=int, default=3)
    ap.add_argument("--keep-over-max-seg-for-readgroup", action="store_true")
    ap.add_argument("--dup-dist", type=int, default=100)
    ap.add_argument("--all-close-leg", action="store_true")
    args = ap.parse_args()

    if args.representatives_only and args.readgroup_cols:
        # This is allowed, but representative-only readgroup files no longer
        # preserve every observed readgroup pair. Keep the mode explicit.
        pass

    args.outdir.mkdir(parents=True, exist_ok=True)
    chrom_order: dict[str, int] = {}
    if args.dedup_unit == "readchain":
        if args.keep_over_max_seg_for_readgroup:
            raise SystemExit("--keep-over-max-seg-for-readgroup is not supported with --dedup-unit readchain")
        representatives, readchain_summary, chrom_lines, seg_hist = collect_readchain_representatives(
            args.seg,
            chrom_order,
            args.min_mapq,
            args.max_seg,
            args.dup_dist,
        )
        modes = ["adjacent", "allcomb"] if args.pair_mode == "both" else [args.pair_mode]
        for mode in modes:
            pairs, pair_summary = pairs_from_readchain_representatives(
                representatives,
                mode,
                chrom_order,
                args.min_leg_dist,
                args.all_close_leg,
                args.coordinate_mode,
            )
            summary: dict[str, object] = {
                "source_label": args.source_label,
                "seg_path": str(args.seg),
                "export_mode": mode,
                "coordinate_mode": args.coordinate_mode,
                "dedup_unit": "readchain",
                "min_mapq": args.min_mapq,
                "min_leg_dist": args.min_leg_dist,
                "max_seg": args.max_seg,
                "keep_over_max_seg_for_readgroup": 0,
                "all_close_leg": int(args.all_close_leg),
                **readchain_summary,
                **pair_summary,
            }
            export_name = args.export_name if len(modes) == 1 else f"{args.export_name}_{mode}"
            write_export_artifacts(
                args.outdir,
                export_name,
                pairs,
                summary,
                chrom_lines,
                seg_hist,
                chrom_order,
                args.readgroup_cols,
                args.representatives_only,
                args.write_annotation,
                mode,
                args.coordinate_mode,
                args.min_mapq,
                args.min_leg_dist,
                args.max_seg,
                args.dup_dist,
                args.dedup_unit,
            )
        return 0
    else:
        if args.pair_mode == "both":
            raise SystemExit("--pair-mode both is only supported with --dedup-unit readchain")
        pairs, summary, chrom_lines, seg_hist = build_pairs(
            args.seg,
            args.source_label,
            args.pair_mode,
            chrom_order,
            args.min_mapq,
            args.min_leg_dist,
            args.max_seg,
            args.all_close_leg,
            args.coordinate_mode,
            args.keep_over_max_seg_for_readgroup,
        )
    write_export_artifacts(
        args.outdir,
        args.export_name,
        pairs,
        summary,
        chrom_lines,
        seg_hist,
        chrom_order,
        args.readgroup_cols,
        args.representatives_only,
        args.write_annotation,
        args.pair_mode,
        args.coordinate_mode,
        args.min_mapq,
        args.min_leg_dist,
        args.max_seg,
        args.dup_dist,
        args.dedup_unit,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
