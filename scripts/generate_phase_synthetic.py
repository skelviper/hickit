#!/usr/bin/env python3
"""Generate one reproducible diploid synthetic dataset with blind/truth separation."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
from dataclasses import asdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from phase_diagnostics.provenance import collect_provenance, file_sha256, write_provenance
from phase_diagnostics.synthetic import SyntheticCondition, generate_synthetic, observation_masses


def require_analysis() -> None:
    if os.environ.get("CONDA_DEFAULT_ENV") != "analysis":
        raise RuntimeError("Synthetic generation must run inside conda environment 'analysis'")


def write_gzip_tsv(path: Path, fields: list[str], rows) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--generator-family", choices=("well_specified", "model_mismatched"), default="well_specified")
    parser.add_argument("--homolog-separation", type=float, default=1.0)
    parser.add_argument("--contacts-per-cell", type=int, default=4000)
    parser.add_argument("--cis-fraction", type=float, default=0.7)
    parser.add_argument("--background-fraction", type=float, default=0.05)
    parser.add_argument("--resolution-bp", type=int, default=1_000_000)
    parser.add_argument("--chromosome-count", type=int, default=3)
    parser.add_argument("--chromosome-length-bins", type=int, default=12)
    parser.add_argument("--readchain-length", type=int, default=2)
    parser.add_argument("--molecule-count", type=int, default=0)
    parser.add_argument("--count-overdispersion", type=float, default=0.0)
    parser.add_argument("--anchor-fraction", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    require_analysis()
    condition = SyntheticCondition(
        generator_family=args.generator_family,
        homolog_separation=args.homolog_separation,
        contacts_per_cell=args.contacts_per_cell,
        cis_fraction=args.cis_fraction,
        background_fraction=args.background_fraction,
        resolution_bp=args.resolution_bp,
        chromosome_count=args.chromosome_count,
        chromosome_length_bins=args.chromosome_length_bins,
        readchain_length=args.readchain_length,
        molecule_count=args.molecule_count,
        count_overdispersion=args.count_overdispersion,
        anchor_fraction=args.anchor_fraction,
        seed=args.seed,
    )
    dataset = generate_synthetic(condition)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "condition.json").write_text(
        json.dumps(asdict(condition), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    blind_rows = []
    truth_rows = []
    anchor_rows = []
    masses = observation_masses(dataset)
    for index, (left, right) in enumerate(zip(dataset.pair_i, dataset.pair_j, strict=True)):
        blind = {
            "observation_id": index,
            "chrom1": f"chr{int(dataset.chrom_index[left]) + 1}",
            "bin1": int(dataset.bin_offset[left]),
            "chrom2": f"chr{int(dataset.chrom_index[right]) + 1}",
            "bin2": int(dataset.bin_offset[right]),
            "molecule_id": (int(dataset.molecule_id[index]) if dataset.molecule_id[index] >= 0 else "NA"),
            "mass": float(masses[index]),
        }
        blind_rows.append(blind)
        truth_rows.append({
            "observation_id": index,
            "state": int(dataset.truth_state[index]),
            "is_background": int(dataset.is_background[index]),
        })
        if dataset.anchor_mask[index]:
            anchor_rows.append({"observation_id": index, "state": int(dataset.truth_state[index])})
    write_gzip_tsv(
        args.output_dir / "blind_contacts.tsv.gz",
        ["observation_id", "chrom1", "bin1", "chrom2", "bin2", "molecule_id", "mass"], blind_rows,
    )
    write_gzip_tsv(args.output_dir / "truth_contacts.eval_only.tsv.gz", ["observation_id", "state", "is_background"], truth_rows)
    write_gzip_tsv(args.output_dir / "phased_anchors.explicit_external.tsv.gz", ["observation_id", "state"], anchor_rows)
    coordinate_rows = []
    for bin_index in range(dataset.n_bins):
        for copy in (0, 1):
            value = dataset.truth_coordinates[bin_index, copy]
            coordinate_rows.append({
                "chrom": f"chr{int(dataset.chrom_index[bin_index]) + 1}",
                "bin": int(dataset.bin_offset[bin_index]), "copy": copy,
                "x": value[0], "y": value[1], "z": value[2],
            })
    write_gzip_tsv(
        args.output_dir / "truth_coordinates.eval_only.tsv.gz",
        ["chrom", "bin", "copy", "x", "y", "z"], coordinate_rows,
    )
    artifacts = sorted(args.output_dir.glob("*.gz"))
    manifest = {
        "boundary": {
            "blind_input": "blind_contacts.tsv.gz",
            "external_anchor_input": "phased_anchors.explicit_external.tsv.gz",
            "eval_only": ["truth_contacts.eval_only.tsv.gz", "truth_coordinates.eval_only.tsv.gz"],
        },
        "sha256": {path.name: file_sha256(path) for path in artifacts},
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    provenance = collect_provenance(
        Path(__file__).resolve().parents[1], sys.argv, args.seed, asdict(condition), ()
    )
    write_provenance(args.output_dir / "provenance.json", provenance)


if __name__ == "__main__":
    main()
