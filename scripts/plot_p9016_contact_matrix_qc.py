#!/usr/bin/env python3
"""Plot observed-vs-experiment contact matrices at a shared bin size."""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_CHROM_ORDER = [f"chr{i}" for i in range(1, 20)] + ["chrX"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", action="append", required=True, help="Matrix label/path as LABEL=PAIRS.gz")
    parser.add_argument("--observed-copy-matrix", action="append", default=[], help="Copy-resolved observed pairs as LABEL=PAIRS.gz")
    parser.add_argument(
        "--synthetic-copy-matrix",
        action="append",
        default=[],
        help="Copy-resolved synthetic sidecar as LABEL=COPY_STATE.tsv.gz",
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--bin-size", type=int, default=1_000_000)
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--exclude-same-bin", action="store_true")
    return parser.parse_args()


def parse_matrix_specs(specs: list[str]) -> list[tuple[str, Path]]:
    parsed: list[tuple[str, Path]] = []
    for spec in specs:
        if "=" not in spec:
            raise ValueError(f"--matrix must be LABEL=PATH, got {spec!r}")
        label, path = spec.rsplit("=", 1)
        label = label.strip()
        if not label:
            raise ValueError(f"empty matrix label in {spec!r}")
        pairs_path = Path(path)
        if not pairs_path.is_file():
            raise FileNotFoundError(pairs_path)
        parsed.append((label, pairs_path))
    return parsed


def parse_optional_specs(specs: list[str]) -> list[tuple[str, Path]]:
    return parse_matrix_specs(specs) if specs else []


def read_chrom_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("#chromosome:"):
                parts = line.strip().split()
                if len(parts) >= 3:
                    sizes[parts[1]] = int(parts[2])
    return sizes


def build_layout(chrom_sizes: dict[str, int], bin_size: int) -> tuple[list[str], dict[str, int], dict[str, int], list[int], list[float]]:
    chroms = [chrom for chrom in DEFAULT_CHROM_ORDER if chrom in chrom_sizes]
    if not chroms:
        raise ValueError("no recognized chromosomes in pairs header")
    bins_per_chrom = {chrom: int(math.ceil(chrom_sizes[chrom] / bin_size)) for chrom in chroms}
    offsets: dict[str, int] = {}
    offset = 0
    boundaries = [0]
    for chrom in chroms:
        offsets[chrom] = offset
        offset += bins_per_chrom[chrom]
        boundaries.append(offset)
    centers = [offsets[chrom] + bins_per_chrom[chrom] / 2 for chrom in chroms]
    return chroms, bins_per_chrom, offsets, boundaries, centers


def build_matrix(
    path: Path,
    offsets: dict[str, int],
    bins_per_chrom: dict[str, int],
    n_bins: int,
    bin_size: int,
    exclude_same_bin: bool,
) -> tuple[np.ndarray, dict[str, int]]:
    mat = np.zeros((n_bins, n_bins), dtype=np.uint32)
    stats = {
        "rows": 0,
        "used_contacts": 0,
        "skipped_same_1mb_bin": 0,
        "skipped_chrom": 0,
        "skipped_bounds": 0,
    }
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                parts = line.strip().split()
            if len(parts) < 5:
                continue
            stats["rows"] += 1
            chrom1, chrom2 = parts[1], parts[3]
            if chrom1 not in offsets or chrom2 not in offsets:
                stats["skipped_chrom"] += 1
                continue
            try:
                pos1, pos2 = int(parts[2]), int(parts[4])
            except ValueError:
                continue
            bin1 = offsets[chrom1] + pos1 // bin_size
            bin2 = offsets[chrom2] + pos2 // bin_size
            if (
                bin1 < offsets[chrom1]
                or bin1 >= offsets[chrom1] + bins_per_chrom[chrom1]
                or bin2 < offsets[chrom2]
                or bin2 >= offsets[chrom2] + bins_per_chrom[chrom2]
            ):
                stats["skipped_bounds"] += 1
                continue
            if exclude_same_bin and bin1 == bin2:
                stats["skipped_same_1mb_bin"] += 1
                continue
            mat[bin1, bin2] += 1
            if bin1 != bin2:
                mat[bin2, bin1] += 1
            stats["used_contacts"] += 1
    stats["matrix_sum_symmetric"] = int(mat.sum())
    stats["nonzero_bins"] = int(np.count_nonzero(mat))
    stats["max_count"] = int(mat.max())
    return mat, stats


def add_copy_contact(
    mats: np.ndarray,
    bin1: int,
    copy1: int,
    bin2: int,
    copy2: int,
) -> None:
    mats[copy1 * 2 + copy2, bin1, bin2] += 1
    if bin1 != bin2 or copy1 != copy2:
        mats[copy2 * 2 + copy1, bin2, bin1] += 1


def build_observed_copy_matrices(
    path: Path,
    offsets: dict[str, int],
    bins_per_chrom: dict[str, int],
    n_bins: int,
    bin_size: int,
    exclude_same_bin: bool,
) -> tuple[np.ndarray, dict[str, int]]:
    mats = np.zeros((4, n_bins, n_bins), dtype=np.uint32)
    columns: list[str] | None = None
    stats = {
        "rows": 0,
        "used_contacts": 0,
        "skipped_missing_phase": 0,
        "skipped_same_1mb_bin": 0,
        "skipped_chrom": 0,
        "skipped_bounds": 0,
    }
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#columns:"):
                columns = line.split(":", 1)[1].strip().split()
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                parts = line.strip().split()
            if columns is None:
                columns = ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2", "phase0", "phase1"]
            row = dict(zip(columns, parts))
            stats["rows"] += 1
            phase1 = row.get("phase0", ".")
            phase2 = row.get("phase1", ".")
            if phase1 not in {"0", "1"} or phase2 not in {"0", "1"}:
                stats["skipped_missing_phase"] += 1
                continue
            chrom1, chrom2 = row["chr1"], row["chr2"]
            if chrom1 not in offsets or chrom2 not in offsets:
                stats["skipped_chrom"] += 1
                continue
            try:
                pos1 = int(row["pos1"])
                pos2 = int(row["pos2"])
            except ValueError:
                continue
            bin1 = offsets[chrom1] + pos1 // bin_size
            bin2 = offsets[chrom2] + pos2 // bin_size
            if (
                bin1 < offsets[chrom1]
                or bin1 >= offsets[chrom1] + bins_per_chrom[chrom1]
                or bin2 < offsets[chrom2]
                or bin2 >= offsets[chrom2] + bins_per_chrom[chrom2]
            ):
                stats["skipped_bounds"] += 1
                continue
            if exclude_same_bin and bin1 == bin2:
                stats["skipped_same_1mb_bin"] += 1
                continue
            add_copy_contact(mats, bin1, int(phase1), bin2, int(phase2))
            stats["used_contacts"] += 1
    update_copy_stats(stats, mats)
    return mats, stats


def build_synthetic_copy_matrices(
    path: Path,
    offsets: dict[str, int],
    bins_per_chrom: dict[str, int],
    n_bins: int,
    bin_size: int,
    exclude_same_bin: bool,
) -> tuple[np.ndarray, dict[str, int]]:
    mats = np.zeros((4, n_bins, n_bins), dtype=np.uint32)
    stats = {
        "rows": 0,
        "used_contacts": 0,
        "skipped_missing_phase": 0,
        "skipped_same_1mb_bin": 0,
        "skipped_chrom": 0,
        "skipped_bounds": 0,
    }
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            stats["rows"] += 1
            chrom1, chrom2 = row["chr1"], row["chr2"]
            if chrom1 not in offsets or chrom2 not in offsets:
                stats["skipped_chrom"] += 1
                continue
            try:
                pos1 = int(row["pos1"])
                pos2 = int(row["pos2"])
                copy1 = int(row["copy1"])
                copy2 = int(row["copy2"])
            except ValueError:
                stats["skipped_missing_phase"] += 1
                continue
            if copy1 not in (0, 1) or copy2 not in (0, 1):
                stats["skipped_missing_phase"] += 1
                continue
            bin1 = offsets[chrom1] + pos1 // bin_size
            bin2 = offsets[chrom2] + pos2 // bin_size
            if (
                bin1 < offsets[chrom1]
                or bin1 >= offsets[chrom1] + bins_per_chrom[chrom1]
                or bin2 < offsets[chrom2]
                or bin2 >= offsets[chrom2] + bins_per_chrom[chrom2]
            ):
                stats["skipped_bounds"] += 1
                continue
            if exclude_same_bin and bin1 == bin2:
                stats["skipped_same_1mb_bin"] += 1
                continue
            add_copy_contact(mats, bin1, copy1, bin2, copy2)
            stats["used_contacts"] += 1
    update_copy_stats(stats, mats)
    return mats, stats


def update_copy_stats(stats: dict[str, int], mats: np.ndarray) -> None:
    states = ("00", "01", "10", "11")
    for idx, state in enumerate(states):
        stats[f"state_{state}_matrix_sum_symmetric"] = int(mats[idx].sum())
        stats[f"state_{state}_nonzero_bins"] = int(np.count_nonzero(mats[idx]))
        stats[f"state_{state}_max_count"] = int(mats[idx].max())


def log_image(mat: np.ndarray) -> np.ndarray:
    return np.log10(mat.astype(np.float32) + 1.0)


def draw_scope(
    *,
    matrices: list[np.ndarray],
    labels: list[str],
    matrix_slice: slice,
    scope_label: str,
    out_path: Path,
    tick_kind: str,
    bin_size: int,
    chroms: list[str],
    boundaries: list[int],
    centers: list[float],
) -> None:
    imgs = [log_image(mat[matrix_slice, matrix_slice]) for mat in matrices]
    nonzero = np.concatenate([img[img > 0] for img in imgs if np.any(img > 0)])
    vmax = float(np.quantile(nonzero, 0.995)) if nonzero.size else 1.0
    fig, axes = plt.subplots(1, len(labels), figsize=(3.4 * len(labels) + 0.8, 3.35), dpi=320, constrained_layout=True)
    if len(labels) == 1:
        axes = [axes]
    im = None
    for ax, label, img in zip(axes, labels, imgs):
        im = ax.imshow(img, origin="lower", interpolation="nearest", cmap="magma", vmin=0.0, vmax=vmax, rasterized=True)
        ax.set_title(label, fontsize=7)
        ax.tick_params(labelsize=7, length=2, width=0.5)
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
        if tick_kind == "chrom":
            n = img.shape[0]
            ticks = np.linspace(0, n - 1, 5)
            tick_labels = [f"{int(round(t))}" for t in np.linspace(0, (n - 1) * bin_size / 1e6, 5)]
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            ax.set_xticklabels(tick_labels)
            ax.set_yticklabels(tick_labels)
            ax.set_xlabel(f"{scope_label} position (Mb)", fontsize=7)
            ax.set_ylabel(f"{scope_label} position (Mb)", fontsize=7)
        else:
            local_boundaries = [b - matrix_slice.start for b in boundaries if matrix_slice.start <= b <= matrix_slice.stop]
            for boundary in local_boundaries:
                ax.axhline(boundary - 0.5, color="white", lw=0.22, alpha=0.55)
                ax.axvline(boundary - 0.5, color="white", lw=0.22, alpha=0.55)
            tick_idx = np.linspace(0, len(chroms) - 1, 5).round().astype(int)
            ticks = [centers[i] - matrix_slice.start for i in tick_idx]
            tick_labels = [chroms[i].replace("chr", "") for i in tick_idx]
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            ax.set_xticklabels(tick_labels)
            ax.set_yticklabels(tick_labels)
            ax.set_xlabel("chromosome", fontsize=7)
            ax.set_ylabel("chromosome", fontsize=7)
    assert im is not None
    cbar = fig.colorbar(im, ax=axes, fraction=0.028, pad=0.012)
    cbar.set_label("log10(contact count + 1)", fontsize=7)
    cbar.ax.tick_params(labelsize=7, length=2, width=0.5)
    suffix = "same-bin contacts excluded" if tick_kind == "chrom" else "same-bin contacts excluded"
    fig.suptitle(f"1 Mb contact matrices, {scope_label}; {suffix}", fontsize=7)
    fig.savefig(out_path, dpi=320)
    plt.close(fig)


def draw_copy_scope(
    *,
    copy_matrices: list[np.ndarray],
    labels: list[str],
    matrix_slice: slice,
    scope_label: str,
    out_path: Path,
    tick_kind: str,
    bin_size: int,
    chroms: list[str],
    boundaries: list[int],
    centers: list[float],
) -> None:
    states = ("00", "01", "10", "11")
    imgs_by_source = [[log_image(mat[state_idx, matrix_slice, matrix_slice]) for state_idx in range(4)] for mat in copy_matrices]
    nonzero_arrays = [img[img > 0] for imgs in imgs_by_source for img in imgs if np.any(img > 0)]
    nonzero = np.concatenate(nonzero_arrays) if nonzero_arrays else np.array([], dtype=float)
    vmax = float(np.quantile(nonzero, 0.995)) if nonzero.size else 1.0
    fig, axes = plt.subplots(
        4,
        len(labels),
        figsize=(3.2 * len(labels) + 0.9, 3.0 * 4),
        dpi=320,
        constrained_layout=True,
        squeeze=False,
    )
    im = None
    for col, label in enumerate(labels):
        for row, state in enumerate(states):
            ax = axes[row, col]
            img = imgs_by_source[col][row]
            im = ax.imshow(img, origin="lower", interpolation="nearest", cmap="magma", vmin=0.0, vmax=vmax, rasterized=True)
            if row == 0:
                ax.set_title(label, fontsize=7)
            if col == 0:
                ax.set_ylabel(f"copy {state[0]}-{state[1]}", fontsize=7)
            ax.tick_params(labelsize=7, length=2, width=0.5)
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)
            if tick_kind == "chrom":
                n = img.shape[0]
                ticks = np.linspace(0, n - 1, 5)
                tick_labels = [f"{int(round(t))}" for t in np.linspace(0, (n - 1) * bin_size / 1e6, 5)]
                ax.set_xticks(ticks)
                ax.set_yticks(ticks)
                ax.set_xticklabels(tick_labels if row == 3 else [])
                ax.set_yticklabels(tick_labels)
                if row == 3:
                    ax.set_xlabel(f"{scope_label} position (Mb)", fontsize=7)
            else:
                local_boundaries = [b - matrix_slice.start for b in boundaries if matrix_slice.start <= b <= matrix_slice.stop]
                for boundary in local_boundaries:
                    ax.axhline(boundary - 0.5, color="white", lw=0.22, alpha=0.55)
                    ax.axvline(boundary - 0.5, color="white", lw=0.22, alpha=0.55)
                tick_idx = np.linspace(0, len(chroms) - 1, 5).round().astype(int)
                ticks = [centers[i] - matrix_slice.start for i in tick_idx]
                tick_labels = [chroms[i].replace("chr", "") for i in tick_idx]
                ax.set_xticks(ticks)
                ax.set_yticks(ticks)
                ax.set_xticklabels(tick_labels if row == 3 else [])
                ax.set_yticklabels(tick_labels)
                if row == 3:
                    ax.set_xlabel("chromosome", fontsize=7)
    assert im is not None
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.018, pad=0.01)
    cbar.set_label("log10(contact count + 1)", fontsize=7)
    cbar.ax.tick_params(labelsize=7, length=2, width=0.5)
    fig.suptitle(f"1 Mb copy-resolved contact matrices, {scope_label}; same-bin contacts excluded", fontsize=7)
    fig.savefig(out_path, dpi=320)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    specs = parse_matrix_specs(args.matrix)
    observed_copy_specs = parse_optional_specs(args.observed_copy_matrix)
    synthetic_copy_specs = parse_optional_specs(args.synthetic_copy_matrix)
    chrom_sizes = read_chrom_sizes(specs[0][1])
    chroms, bins_per_chrom, offsets, boundaries, centers = build_layout(chrom_sizes, args.bin_size)
    if args.chrom not in offsets:
        raise ValueError(f"{args.chrom} not found in pairs header")
    n_bins = boundaries[-1]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = args.out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    matrices = []
    stats_rows = []
    labels = []
    for label, path in specs:
        mat, stats = build_matrix(path, offsets, bins_per_chrom, n_bins, args.bin_size, args.exclude_same_bin)
        matrices.append(mat)
        labels.append(label)
        stats_rows.append((label, path, stats))

    chrom_slice = slice(offsets[args.chrom], offsets[args.chrom] + bins_per_chrom[args.chrom])
    draw_scope(
        matrices=matrices,
        labels=labels,
        matrix_slice=chrom_slice,
        scope_label=args.chrom,
        out_path=plot_dir / f"p9016_{args.chrom}_contact_matrix_1mb_observed_vs_charm_synthetic.png",
        tick_kind="chrom",
        bin_size=args.bin_size,
        chroms=chroms,
        boundaries=boundaries,
        centers=centers,
    )
    draw_scope(
        matrices=matrices,
        labels=labels,
        matrix_slice=slice(0, n_bins),
        scope_label="whole genome",
        out_path=plot_dir / "p9016_genome_contact_matrix_1mb_observed_vs_charm_synthetic.png",
        tick_kind="genome",
        bin_size=args.bin_size,
        chroms=chroms,
        boundaries=boundaries,
        centers=centers,
    )

    with (args.out_dir / "matrix_plot_summary.tsv").open("w") as out:
        out.write(
            "matrix_label\tpairs_path\tbin_size_bp\tfilter_same_1mb_bin\trows\tused_contacts\t"
            "skipped_same_1mb_bin\tskipped_chrom\tskipped_bounds\tmatrix_sum_symmetric\tnonzero_bins\tmax_count\n"
        )
        for label, path, stats in stats_rows:
            out.write(
                "\t".join(
                    [
                        label,
                        str(path),
                        str(args.bin_size),
                        "1" if args.exclude_same_bin else "0",
                        str(stats["rows"]),
                        str(stats["used_contacts"]),
                        str(stats["skipped_same_1mb_bin"]),
                        str(stats["skipped_chrom"]),
                        str(stats["skipped_bounds"]),
                        str(stats["matrix_sum_symmetric"]),
                        str(stats["nonzero_bins"]),
                        str(stats["max_count"]),
                    ]
                )
                + "\n"
            )
    copy_specs = [("observed", label, path) for label, path in observed_copy_specs]
    copy_specs.extend(("synthetic", label, path) for label, path in synthetic_copy_specs)
    if copy_specs:
        copy_matrices = []
        copy_labels = []
        copy_stats_rows = []
        for source_type, label, path in copy_specs:
            if source_type == "observed":
                mat, stats = build_observed_copy_matrices(path, offsets, bins_per_chrom, n_bins, args.bin_size, args.exclude_same_bin)
            else:
                mat, stats = build_synthetic_copy_matrices(path, offsets, bins_per_chrom, n_bins, args.bin_size, args.exclude_same_bin)
            copy_matrices.append(mat)
            copy_labels.append(label)
            copy_stats_rows.append((source_type, label, path, stats))
        draw_copy_scope(
            copy_matrices=copy_matrices,
            labels=copy_labels,
            matrix_slice=chrom_slice,
            scope_label=args.chrom,
            out_path=plot_dir / f"p9016_{args.chrom}_copy_resolved_contact_matrix_1mb_observed_vs_charm_synthetic.png",
            tick_kind="chrom",
            bin_size=args.bin_size,
            chroms=chroms,
            boundaries=boundaries,
            centers=centers,
        )
        draw_copy_scope(
            copy_matrices=copy_matrices,
            labels=copy_labels,
            matrix_slice=slice(0, n_bins),
            scope_label="whole genome",
            out_path=plot_dir / "p9016_genome_copy_resolved_contact_matrix_1mb_observed_vs_charm_synthetic.png",
            tick_kind="genome",
            bin_size=args.bin_size,
            chroms=chroms,
            boundaries=boundaries,
            centers=centers,
        )
        with (args.out_dir / "copy_resolved_matrix_plot_summary.tsv").open("w") as out:
            state_cols = []
            for state in ("00", "01", "10", "11"):
                state_cols.extend(
                    [
                        f"state_{state}_matrix_sum_symmetric",
                        f"state_{state}_nonzero_bins",
                        f"state_{state}_max_count",
                    ]
                )
            out.write(
                "source_type\tmatrix_label\tpath\tbin_size_bp\tfilter_same_1mb_bin\trows\tused_contacts\t"
                "skipped_missing_phase\tskipped_same_1mb_bin\tskipped_chrom\tskipped_bounds\t"
                + "\t".join(state_cols)
                + "\n"
            )
            for source_type, label, path, stats in copy_stats_rows:
                out.write(
                    "\t".join(
                        [
                            source_type,
                            label,
                            str(path),
                            str(args.bin_size),
                            "1" if args.exclude_same_bin else "0",
                            str(stats["rows"]),
                            str(stats["used_contacts"]),
                            str(stats["skipped_missing_phase"]),
                            str(stats["skipped_same_1mb_bin"]),
                            str(stats["skipped_chrom"]),
                            str(stats["skipped_bounds"]),
                        ]
                        + [str(stats[col]) for col in state_cols]
                    )
                    + "\n"
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
