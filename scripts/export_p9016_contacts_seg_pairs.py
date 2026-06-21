#!/usr/bin/env python3
"""Export blind-safe pairs from a CHARM contacts.seg.gz file.

The exporter preserves readID but drops phase labels by writing phase0/phase1
as ".". This keeps the output usable as a raw-contact training source without
leaking SNP phase into blind training.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import random
from pathlib import Path
from typing import TextIO


def open_text(path: Path, mode: str = "rt") -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    rest = chrom[3:] if chrom.startswith("chr") else chrom
    if rest.isdigit():
        return (int(rest), "")
    if rest == "X":
        return (1000, "")
    if rest == "Y":
        return (1001, "")
    return (2000, rest)


def parse_segment(field: str) -> tuple[str, int, str, int | None] | None:
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
    strand = parts[3] if parts[3] in {"+", "-"} else "+"
    return parts[0], (start + end) // 2, strand, mapq


def canonical_pair(a: tuple[str, int, str, int | None], b: tuple[str, int, str, int | None]) -> tuple[tuple[str, int, str, int | None], tuple[str, int, str, int | None]]:
    ak = (chrom_sort_key(a[0]), a[1], a[0])
    bk = (chrom_sort_key(b[0]), b[1], b[0])
    if ak <= bk:
        return a, b
    return b, a


def should_keep_pair(
    a: tuple[str, int, str, int | None],
    b: tuple[str, int, str, int | None],
    mode: str,
    n_segments: int,
    n_chroms: int,
    bin_size: int,
    drop_same_bin: bool,
) -> bool:
    if mode == "all":
        pass
    elif mode == "multiseg":
        if n_segments < 3:
            return False
    elif mode == "multichrom":
        if n_chroms < 2:
            return False
        if a[0] == b[0]:
            return False
    elif mode == "multiseg_multichrom":
        if n_segments < 3 or n_chroms < 2:
            return False
        if a[0] == b[0]:
            return False
    else:
        raise ValueError(f"unknown mode: {mode}")
    if drop_same_bin and a[0] == b[0] and a[1] // bin_size == b[1] // bin_size:
        return False
    return True


def write_header(out: TextIO, chrom_lines: list[str], write_readgroup_cols: bool = False) -> None:
	out.write("## pairs format v1.0\n")
	out.write("#sorted: chr1-chr2-pos1-pos2\n")
	out.write("#shape: upper triangle\n")
	for line in chrom_lines:
		out.write(line.rstrip("\n") + "\n")
	if write_readgroup_cols:
		out.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1 read_group_id read_seg0 read_seg1 read_n_segments\n")
	else:
		out.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seg", required=True, type=Path)
    ap.add_argument("--out-pairs", required=True, type=Path)
    ap.add_argument("--metadata", required=True, type=Path)
    ap.add_argument("--mode", choices=["all", "multiseg", "multichrom", "multiseg_multichrom"], default="all")
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--bin-size", type=int, default=1_000_000)
    ap.add_argument("--drop-same-bin", action="store_true")
    ap.add_argument("--max-pairs", type=int, default=0, help="Reservoir-sample this many pairs after filtering; 0 keeps all.")
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--write-readgroup-cols", action="store_true", help="Add numeric read-group metadata columns for blind readgroup-aware training.")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    chrom_lines: list[str] = []
    candidates: list[tuple[str, int, tuple[str, int, str, int | None], tuple[str, int, str, int | None], int, int, int]] = []
    n_reads = n_usable_reads = n_pairs_seen = n_pairs_kept = 0
    n_pairs_cis = n_pairs_trans = 0
    read_group_id = 0

    with open_text(args.seg) as fh:
        for line in fh:
            if line.startswith("#chromosome:") or line.startswith("#chromsize:"):
                chrom_lines.append(line)
                continue
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            n_reads += 1
            read_id = fields[0]
            raw_segments = [parse_segment(field) for field in fields[1:]]
            segments = [seg for seg in raw_segments if seg is not None and (seg[3] is None or seg[3] >= args.min_mapq)]
            if len(segments) < 2:
                continue
            n_usable_reads += 1
            current_group_id = read_group_id
            read_group_id += 1
            n_chroms = len({seg[0] for seg in segments})
            for ia, ib in itertools.combinations(range(len(segments)), 2):
                a = segments[ia]
                b = segments[ib]
                if not should_keep_pair(a, b, args.mode, len(segments), n_chroms, args.bin_size, args.drop_same_bin):
                    continue
                n_pairs_seen += 1
                ca, cb = canonical_pair(a, b)
                if ca == a and cb == b:
                    seg0, seg1 = ia, ib
                else:
                    seg0, seg1 = ib, ia
                if ca[0] == cb[0]:
                    n_pairs_cis += 1
                else:
                    n_pairs_trans += 1
                record = (read_id, current_group_id, ca, cb, seg0, seg1, len(segments))
                if args.max_pairs <= 0:
                    candidates.append(record)
                elif len(candidates) < args.max_pairs:
                    candidates.append(record)
                else:
                    j = rng.randrange(n_pairs_seen)
                    if j < args.max_pairs:
                        candidates[j] = record
    candidates.sort(key=lambda r: (chrom_sort_key(r[2][0]), chrom_sort_key(r[3][0]), r[2][1], r[3][1], r[0], r[1], r[4], r[5]))
    n_pairs_kept = len(candidates)

    args.out_pairs.parent.mkdir(parents=True, exist_ok=True)
    args.metadata.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.out_pairs, "wt") as out:
        write_header(out, chrom_lines, args.write_readgroup_cols)
        for read_id, group_id, a, b, seg0, seg1, n_segments in candidates:
            if args.write_readgroup_cols:
                out.write(f"{read_id}\t{a[0]}\t{a[1]}\t{b[0]}\t{b[1]}\t{a[2]}\t{b[2]}\t.\t.\t{group_id}\t{seg0}\t{seg1}\t{n_segments}\n")
            else:
                out.write(f"{read_id}\t{a[0]}\t{a[1]}\t{b[0]}\t{b[1]}\t{a[2]}\t{b[2]}\t.\t.\n")

    with args.metadata.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        writer.writerow(["seg_path", str(args.seg)])
        writer.writerow(["out_pairs", str(args.out_pairs)])
        writer.writerow(["mode", args.mode])
        writer.writerow(["min_mapq", args.min_mapq])
        writer.writerow(["bin_size_bp", args.bin_size])
        writer.writerow(["drop_same_bin", int(args.drop_same_bin)])
        writer.writerow(["max_pairs", args.max_pairs])
        writer.writerow(["seed", args.seed])
        writer.writerow(["write_readgroup_cols", int(args.write_readgroup_cols)])
        writer.writerow(["n_reads", n_reads])
        writer.writerow(["n_usable_reads_ge2_segments", n_usable_reads])
        writer.writerow(["n_pairs_after_filter_before_sampling", n_pairs_seen])
        writer.writerow(["n_pairs_written", n_pairs_kept])
        writer.writerow(["n_pairs_after_filter_cis_before_sampling", n_pairs_cis])
        writer.writerow(["n_pairs_after_filter_trans_before_sampling", n_pairs_trans])
        writer.writerow(["phase_labels_written", 0])
        writer.writerow(["training_uses_phase_labels", 0])
        writer.writerow(["training_uses_charm_or_reference", 0])
        writer.writerow(["readgroup_columns_written", int(args.write_readgroup_cols)])


if __name__ == "__main__":
    main()
