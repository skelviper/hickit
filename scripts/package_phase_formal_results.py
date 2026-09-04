#!/usr/bin/env python3
"""Package logical phase diagnostics 057-059 into new formal test_res folders.

The default mode is read-only. Pass --execute explicitly to create the three
folders. Packaging never reuses a sequence number or overwrites a destination.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import socket
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Mapping, Sequence


REPO = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS_ROOT = REPO / "results/phase_diagnostics"
DEFAULT_TEST_RES_ROOT = REPO / "test_res"
DEFAULT_057_STANDARD = Path(
    "/work/phase3/hickit/test_res/"
    "083-20260713_111134-p9016_phased_hickit_oracle_gate_1m/plots"
)
DEFAULT_058_STANDARD = Path(
    "/work/phase3/hickit/test_res/"
    "078-20260712_211825-p9016_corrected_combined_candidate_1m/plots/K1_final"
)
STANDARD_PLOT_NAMES = (
    "all_chrom_3d_scatter.png",
    "chr1_copy_3d_scatter.png",
    "chr1_distance_maps.png",
)
TIMESTAMP_PATTERN = re.compile(r"^[0-9]{8}_[0-9]{6}$")
SEQUENCE_PATTERN = re.compile(r"^([0-9]+)-")
MACHINE_SUFFIXES = {".tsv", ".csv", ".json"}


@dataclass(frozen=True)
class ExperimentSpec:
    logical_id: str
    slug: str
    required_files: tuple[str, ...]
    required_pdfs: tuple[str, ...]


@dataclass(frozen=True)
class Artifact:
    artifact: str
    role: str
    source: str
    source_sha256: str
    packaged_sha256: str
    truth_boundary: str
    copy_swap_policy: str


@dataclass(frozen=True)
class Exp057BStatus:
    declared_status: str
    complete: bool
    seeds: tuple[int, ...]
    rows: tuple[Mapping[str, str], ...]
    problems: tuple[str, ...]


@dataclass(frozen=True)
class ProvenanceIssues:
    fatal: tuple[str, ...]
    dirty: tuple[str, ...]


SPECS = (
    ExperimentSpec(
        logical_id="057",
        slug="phase_exp057_oracle_decomposition_1m",
        required_files=(
            "run_manifest.tsv",
            "metrics_long.tsv",
            "per_chromosome.tsv",
            "calibration.tsv",
            "selective_calling.tsv",
            "near_truth_initialization.tsv",
            "provenance.json",
            "configs/experiment_057.json",
        ),
        required_pdfs=(
            "oracle_coordinate_estep_accuracy.pdf",
            "selective_calling.pdf",
            "near_truth_initialization.pdf",
        ),
    ),
    ExperimentSpec(
        logical_id="058",
        slug="phase_exp058_heldout_invariants_partial_1m",
        required_files=(
            "split_manifest.tsv",
            "heldout_metrics.tsv",
            "seed_stability.tsv",
            "model_selection.tsv",
            "provenance.json",
        ),
        required_pdfs=("blocked_split_accuracy.pdf", "seed_stability.pdf"),
    ),
    ExperimentSpec(
        logical_id="059",
        slug="phase_exp059_identifiability_synthetic_1m",
        required_files=(
            "conditions.tsv",
            "metrics_long.tsv",
            "replicate_summary.tsv",
            "config.json",
            "provenance.json",
        ),
        required_pdfs=(
            "contacts_x_homolog_separation.pdf",
            "trans_fraction_x_background.pdf",
            "anchor_fraction_x_contacts.pdf",
            "readchain_length_x_molecule_count.pdf",
        ),
    ),
)

EXPECTED_057B_ARMS = (
    "native_hickit_cis_only",
    "native_hickit_cis_plus_trans",
    "experimental_fixed_cis_rigid_trans_placement",
)
LEGACY_UNSUPPORTED_057B_CONFIGS = {
    "oracle_phase_cis_only_reconstruction",
    "oracle_phase_cis_plus_trans_reconstruction",
    "oracle_phase_trans_placement_fixed_cis_shapes",
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require_analysis() -> None:
    if os.environ.get("CONDA_DEFAULT_ENV") != "analysis":
        raise RuntimeError(
            "package_phase_formal_results.py must run inside conda environment 'analysis'"
        )


def git_output(*args: str) -> str:
    return subprocess.check_output(
        ("git", *args), cwd=REPO, text=True, stderr=subprocess.DEVNULL
    ).rstrip("\r\n")


def read_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def manifest_fields(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    if not rows or "field" not in rows[0] or "value" not in rows[0]:
        raise ValueError(f"{path}: expected field/value columns")
    fields: dict[str, str] = {}
    for row in rows:
        field = row.get("field", "")
        if not field:
            raise ValueError(f"{path}: empty manifest field")
        if field in fields:
            raise ValueError(f"{path}: duplicate manifest field {field!r}")
        fields[field] = row.get("value", "")
    return fields


def parse_seed_csv(value: str) -> tuple[int, ...]:
    try:
        seeds = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as error:
        raise ValueError(f"invalid comma-separated seed list {value!r}") from error
    if not seeds or len(seeds) != len(set(seeds)) or any(seed < 0 for seed in seeds):
        raise ValueError(f"invalid comma-separated seed list {value!r}")
    return seeds


def assess_057b(source_root: Path) -> Exp057BStatus:
    """Cross-check the declared 057B status against all current result tables."""
    manifest = manifest_fields(source_root / "run_manifest.tsv")
    declared = manifest.get("057B_status", "UNDECLARED")
    run_path = source_root / "oracle_mstep_runs.tsv"
    long_rows = read_tsv(source_root / "metrics_long.tsv")
    chromosome_rows = read_tsv(source_root / "per_chromosome.tsv")
    problems: list[str] = []

    legacy = sorted(
        {
            row.get("configuration", "")
            for row in long_rows
            if row.get("configuration", "") in LEGACY_UNSUPPORTED_057B_CONFIGS
        }
    )
    if legacy:
        problems.append(
            "metrics_long.tsv contains unsupported historical pseudo-arms: "
            + ", ".join(legacy)
        )

    if not run_path.is_file():
        if declared == "COMPLETE_THREE_EXPLICIT_ARMS":
            problems.append("run_manifest declares complete 057B but oracle_mstep_runs.tsv is missing")
        return Exp057BStatus(declared, False, (), (), tuple(problems))

    rows = read_tsv(run_path)
    required_columns = {
        "arm",
        "seed",
        "n_chromosomes",
        "mean_copy_swap_invariant_distance_spearman",
        "runtime_seconds",
        "truth_used_in_inference",
    }
    columns = set(rows[0]) if rows else set()
    missing_columns = sorted(required_columns - columns)
    if missing_columns:
        problems.append(
            f"oracle_mstep_runs.tsv lacks columns: {', '.join(missing_columns)}"
        )

    try:
        seeds = parse_seed_csv(manifest.get("057B_seeds", ""))
    except ValueError as error:
        seeds = ()
        problems.append(str(error))
    expected_keys = {(arm, str(seed)) for arm in EXPECTED_057B_ARMS for seed in seeds}
    observed_keys = {(row.get("arm", ""), row.get("seed", "")) for row in rows}
    missing_keys = sorted(expected_keys - observed_keys)
    extra_keys = sorted(observed_keys - expected_keys)
    unexpected_arms = sorted({arm for arm, _ in observed_keys if arm not in EXPECTED_057B_ARMS})
    if missing_keys:
        problems.append(
            "oracle_mstep_runs.tsv lacks arm/seed rows: "
            + ", ".join(f"{arm}/seed{seed}" for arm, seed in missing_keys)
        )
    if unexpected_arms:
        problems.append(
            "oracle_mstep_runs.tsv contains unexpected arms: " + ", ".join(unexpected_arms)
        )
    if extra_keys:
        problems.append(
            "oracle_mstep_runs.tsv contains undeclared arm/seed rows: "
            + ", ".join(f"{arm}/seed{seed}" for arm, seed in extra_keys)
        )
    if len(observed_keys) != len(rows):
        problems.append("oracle_mstep_runs.tsv contains duplicate arm/seed rows")
    bad_truth = [
        f"{row.get('arm', '')}/seed{row.get('seed', '')}"
        for row in rows
        if row.get("truth_used_in_inference") != "ORACLE_CONTACT_PHASE"
    ]
    if bad_truth:
        problems.append(
            "057B rows lack ORACLE_CONTACT_PHASE boundary: " + ", ".join(bad_truth)
        )

    explicit_long_rows = [
        row for row in long_rows
        if row.get("component") == "057B_explicit_oracle_phase_mstep"
    ]
    explicit_chromosome_rows = [
        row for row in chromosome_rows
        if row.get("component") == "057B_oracle_phase_mstep"
    ]
    bad_long_truth = [
        f"{row.get('configuration', '')}/{row.get('metric', '')}"
        for row in explicit_long_rows
        if row.get("truth_used_in_inference") != "ORACLE_CONTACT_PHASE"
    ]
    bad_chromosome_truth = [
        f"{row.get('configuration', '')}/{row.get('chromosome', '')}"
        for row in explicit_chromosome_rows
        if row.get("truth_used_in_inference") != "ORACLE_CONTACT_PHASE"
    ]
    if bad_long_truth:
        problems.append(
            "metrics_long.tsv has non-ORACLE explicit 057B rows: "
            + ", ".join(bad_long_truth)
        )
    if bad_chromosome_truth:
        problems.append(
            "per_chromosome.tsv has non-ORACLE explicit 057B rows: "
            + ", ".join(bad_chromosome_truth)
        )
    explicit_long = {
        row.get("configuration", "")
        for row in explicit_long_rows
    }
    explicit_chromosome = {
        row.get("configuration", "")
        for row in explicit_chromosome_rows
    }
    expected_configurations = {f"{arm}_seed{seed}" for arm, seed in expected_keys}
    missing_long = sorted(expected_configurations - explicit_long)
    missing_chromosome = sorted(expected_configurations - explicit_chromosome)
    extra_long = sorted(explicit_long - expected_configurations)
    extra_chromosome = sorted(explicit_chromosome - expected_configurations)
    if missing_long:
        problems.append(
            "metrics_long.tsv lacks explicit 057B configurations: " + ", ".join(missing_long)
        )
    if missing_chromosome:
        problems.append(
            "per_chromosome.tsv lacks explicit 057B configurations: "
            + ", ".join(missing_chromosome)
        )
    if extra_long:
        problems.append(
            "metrics_long.tsv contains undeclared explicit 057B configurations: "
            + ", ".join(extra_long)
        )
    if extra_chromosome:
        problems.append(
            "per_chromosome.tsv contains undeclared explicit 057B configurations: "
            + ", ".join(extra_chromosome)
        )

    long_by_configuration: dict[str, dict[str, list[str]]] = {}
    for row in long_rows:
        configuration = row.get("configuration", "")
        if configuration not in expected_configurations or row.get("component") != "057B_explicit_oracle_phase_mstep":
            continue
        long_by_configuration.setdefault(configuration, {}).setdefault(
            row.get("metric", ""), []
        ).append(row.get("value", ""))
    chromosome_by_configuration: dict[str, list[Mapping[str, str]]] = {
        configuration: [
            row for row in chromosome_rows
            if row.get("component") == "057B_oracle_phase_mstep"
            and row.get("configuration") == configuration
        ]
        for configuration in expected_configurations
    }
    metadata_fields = {"arm", "seed", "truth_used_in_inference", "placement_optimizer"}
    required_run_metrics = {
        "n_chromosomes",
        "mean_copy_swap_invariant_distance_spearman",
        "mean_homolog_separation_spearman",
        "chromosome_centroid_distance_spearman",
        "runtime_seconds",
    }
    for run_row in rows:
        arm, seed = run_row.get("arm", ""), run_row.get("seed", "")
        configuration = f"{arm}_seed{seed}"
        if configuration not in expected_configurations:
            continue
        for metric in required_run_metrics:
            try:
                value = float(run_row.get(metric, "nan"))
            except ValueError:
                value = math.nan
            if not math.isfinite(value):
                problems.append(f"oracle_mstep_runs.tsv has non-finite {configuration}/{metric}")
        try:
            expected_chromosomes = int(run_row.get("n_chromosomes", "0"))
        except ValueError:
            expected_chromosomes = 0
        chromosome_subset = chromosome_by_configuration[configuration]
        chromosome_names = [row.get("chromosome", "") for row in chromosome_subset]
        if expected_chromosomes <= 0 or len(chromosome_subset) != expected_chromosomes:
            problems.append(
                f"per_chromosome.tsv has {len(chromosome_subset)} rows for {configuration}, "
                f"expected {expected_chromosomes}"
            )
        if len(chromosome_names) != len(set(chromosome_names)) or any(not name for name in chromosome_names):
            problems.append(f"per_chromosome.tsv has duplicate/empty chromosomes for {configuration}")
        chromosome_metrics = (
            "procrustes_rmsd",
            "rigid_rmsd",
            "distance_pearson",
            "distance_spearman",
            "chromosome_shape_preservation",
            "homolog_separation_spearman",
        )
        for chromosome_row in chromosome_subset:
            for metric in chromosome_metrics:
                try:
                    value = float(chromosome_row.get(metric, "nan"))
                except ValueError:
                    value = math.nan
                if not math.isfinite(value):
                    problems.append(
                        f"per_chromosome.tsv has non-finite {configuration}/"
                        f"{chromosome_row.get('chromosome', '')}/{metric}"
                    )
        for metric, raw_value in run_row.items():
            if metric in metadata_fields or raw_value == "":
                continue
            try:
                expected_value = float(raw_value)
            except ValueError:
                continue
            observed_values = long_by_configuration.get(configuration, {}).get(metric, [])
            if len(observed_values) != 1:
                problems.append(
                    f"metrics_long.tsv has {len(observed_values)} rows for {configuration}/{metric}"
                )
                continue
            try:
                observed_value = float(observed_values[0])
            except ValueError:
                observed_value = math.nan
            if not math.isfinite(expected_value) or not math.isfinite(observed_value) or not math.isclose(
                expected_value, observed_value, rel_tol=1e-12, abs_tol=1e-12
            ):
                problems.append(f"057B table mismatch for {configuration}/{metric}")
        if chromosome_subset:
            summary_mapping = {
                "procrustes_rmsd": "mean_copy_swap_invariant_procrustes_rmsd",
                "rigid_rmsd": "mean_copy_swap_invariant_rigid_rmsd",
                "distance_pearson": "mean_copy_swap_invariant_distance_pearson",
                "distance_spearman": "mean_copy_swap_invariant_distance_spearman",
                "chromosome_shape_preservation": "mean_chromosome_shape_preservation",
                "homolog_separation_spearman": "mean_homolog_separation_spearman",
            }
            for chromosome_metric, summary_metric in summary_mapping.items():
                per_chromosome_mean = sum(
                    float(row[chromosome_metric]) for row in chromosome_subset
                ) / len(chromosome_subset)
                try:
                    run_mean = float(run_row[summary_metric])
                except (KeyError, ValueError):
                    run_mean = math.nan
                if not math.isfinite(run_mean) or not math.isclose(
                    per_chromosome_mean, run_mean, rel_tol=1e-12, abs_tol=1e-12
                ):
                    problems.append(
                        f"057B chromosome-summary mismatch for {configuration}/{summary_metric}"
                    )
    if declared != "COMPLETE_THREE_EXPLICIT_ARMS":
        problems.append(
            f"run_manifest 057B_status={declared!r}, expected 'COMPLETE_THREE_EXPLICIT_ARMS'"
        )
    complete = bool(rows and seeds and not problems)
    return Exp057BStatus(declared, complete, seeds, tuple(rows), tuple(problems))


def non_output_git_status(lines: Iterable[str]) -> list[str]:
    """Exclude generated result roots when deciding whether source code was dirty."""
    retained: list[str] = []
    for line in lines:
        path = line[3:].strip() if len(line) >= 4 else line.strip()
        paths = [part.strip() for part in path.split(" -> ")]
        if paths and all(
            item == "results"
            or item.startswith("results/")
            or item == "test_res"
            or item.startswith("test_res/")
            for item in paths
        ):
            continue
        retained.append(line)
    return retained


def format_float(value: object, digits: int = 6) -> str:
    try:
        number = float(str(value))
    except (TypeError, ValueError):
        return "NA"
    if not (number == number and abs(number) != float("inf")):
        return "NA"
    return f"{number:.{digits}f}"


def markdown_escape(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def shell_command(argv: Sequence[str]) -> str:
    return shlex.join([str(item) for item in argv])


def discover_sequence_start(test_res_root: Path) -> int:
    maximum = 0
    if test_res_root.exists():
        for child in test_res_root.iterdir():
            match = SEQUENCE_PATTERN.match(child.name)
            if match:
                maximum = max(maximum, int(match.group(1)))
    return maximum + 1


def occupied_sequences(test_res_root: Path) -> dict[int, list[Path]]:
    occupied: dict[int, list[Path]] = {}
    if not test_res_root.exists():
        return occupied
    for child in test_res_root.iterdir():
        match = SEQUENCE_PATTERN.match(child.name)
        if match:
            occupied.setdefault(int(match.group(1)), []).append(child)
    return occupied


def validate_timestamp(value: str) -> str:
    if not TIMESTAMP_PATTERN.fullmatch(value):
        raise ValueError("--timestamp must use YYYYMMDD_HHMMSS")
    datetime.strptime(value, "%Y%m%d_%H%M%S")
    return value


def target_paths(
    test_res_root: Path, sequence_start: int, timestamp: str
) -> dict[str, Path]:
    if sequence_start < 1:
        raise ValueError("--sequence-start must be positive")
    occupied = occupied_sequences(test_res_root)
    targets: dict[str, Path] = {}
    for offset, spec in enumerate(SPECS):
        sequence = sequence_start + offset
        if sequence in occupied:
            names = ", ".join(str(path) for path in occupied[sequence])
            raise FileExistsError(f"formal sequence {sequence} is already occupied by: {names}")
        destination = test_res_root / f"{sequence:03d}-{timestamp}-{spec.slug}"
        if destination.exists():
            raise FileExistsError(f"refusing to overwrite existing destination: {destination}")
        targets[spec.logical_id] = destination
    return targets


def validate_059_standard(
    root: Path, seed: int, *, require_closure_manifest: bool = False
) -> list[str]:
    """Validate that reusable/generated standard plots and tables form one run."""
    problems: list[str] = []
    plot_paths = [root / "plots" / name for name in STANDARD_PLOT_NAMES]
    table_paths = [
        root / "summary.tsv",
        root / "per_chromosome.tsv",
        root / "standard_plot_config.json",
    ]
    for path in (*plot_paths, *table_paths):
        if not path.is_file():
            problems.append(f"missing representative 059 artifact: {path}")
    if problems:
        return problems

    config_path = root / "standard_plot_config.json"
    config = read_json(config_path)
    condition = config.get("condition", {})
    inference = config.get("inference", {})
    if not isinstance(condition, dict):
        problems.append(f"{config_path}: condition must be a JSON object")
        condition = {}
    if not isinstance(inference, dict):
        problems.append(f"{config_path}: inference must be a JSON object")
        inference = {}
    try:
        if int(condition.get("resolution_bp", -1)) != 1_000_000:
            problems.append(f"{config_path}: representative resolution is not 1 Mb")
        if int(condition.get("seed", -1)) != seed:
            problems.append(
                f"{config_path}: representative seed {condition.get('seed')!r} != requested {seed}"
            )
        if int(inference.get("seed", -1)) != seed:
            problems.append(
                f"{config_path}: inference seed {inference.get('seed')!r} != requested {seed}"
            )
        if float(condition.get("anchor_fraction", -1.0)) != 0.0:
            problems.append(f"{config_path}: representative condition must have zero anchors")
        chromosome_count = int(condition.get("chromosome_count", -1))
    except (TypeError, ValueError) as error:
        problems.append(f"{config_path}: invalid numeric condition/inference field ({error})")
        chromosome_count = -1
    if str(config.get("truth_used_in_inference", "")) != "none":
        problems.append(f"{config_path}: truth_used_in_inference must be 'none'")
    contract = {
        "condition": condition,
        "inference": inference,
        "truth_used_in_inference": config.get("truth_used_in_inference"),
    }
    serialized_contract = json.dumps(contract, sort_keys=True, separators=(",", ":"))
    calculated_contract_hash = hashlib.sha256(serialized_contract.encode("utf-8")).hexdigest()
    run_contract_sha256 = str(config.get("run_contract_sha256", ""))
    run_id = str(config.get("run_id", ""))
    if run_contract_sha256 != calculated_contract_hash:
        problems.append(f"{config_path}: run_contract_sha256 does not match condition/inference")
    if not run_id or run_id != f"phase059_standard_{run_contract_sha256[:16]}":
        problems.append(f"{config_path}: run_id does not match run_contract_sha256")

    summary_path = root / "summary.tsv"
    summary_rows = read_tsv(summary_path)
    if len(summary_rows) != 1:
        problems.append(f"{summary_path}: expected exactly one summary row")
    else:
        summary = summary_rows[0]
        if summary.get("truth_used_in_inference") != "none":
            problems.append(f"{summary_path}: truth_used_in_inference must be 'none'")
        if summary.get("copy_swap_policy") != "ORACLE_GEOMETRY_ALIGNMENT_EVAL_ONLY":
            problems.append(f"{summary_path}: copy-swap policy is missing or unexpected")
        for field, expected in (
            ("run_id", run_id),
            ("run_contract_sha256", run_contract_sha256),
            ("seed", str(seed)),
            ("generator_family", str(condition.get("generator_family", ""))),
        ):
            if summary.get(field) != expected:
                problems.append(f"{summary_path}: {field} does not match standard_plot_config.json")
        try:
            if int(summary.get("resolution_bp", -1)) != 1_000_000:
                problems.append(f"{summary_path}: resolution is not 1 Mb")
            if int(summary.get("n_contacts", 0)) <= 0:
                problems.append(f"{summary_path}: n_contacts must be positive")
            if int(summary.get("converged", 0)) != 1:
                problems.append(f"{summary_path}: representative inference did not converge")
        except (TypeError, ValueError) as error:
            problems.append(f"{summary_path}: invalid numeric summary field ({error})")

    chromosome_path = root / "per_chromosome.tsv"
    chromosome_rows = read_tsv(chromosome_path)
    chromosomes = [row.get("chromosome", "") for row in chromosome_rows]
    if chromosome_count >= 0 and len(chromosome_rows) != chromosome_count:
        problems.append(
            f"{chromosome_path}: {len(chromosome_rows)} rows != chromosome_count {chromosome_count}"
        )
    if len(chromosomes) != len(set(chromosomes)) or any(not value for value in chromosomes):
        problems.append(f"{chromosome_path}: chromosome labels must be nonempty and unique")
    bad_flips = [
        row.get("chromosome", "")
        for row in chromosome_rows
        if row.get("geometry_selected_copy_flip") not in {"0", "1"}
    ]
    if bad_flips:
        problems.append(
            f"{chromosome_path}: invalid geometry-selected copy flip for {', '.join(bad_flips)}"
        )
    for row in chromosome_rows:
        for field, expected in (
            ("run_id", run_id),
            ("run_contract_sha256", run_contract_sha256),
            ("seed", str(seed)),
            ("generator_family", str(condition.get("generator_family", ""))),
        ):
            if row.get(field) != expected:
                problems.append(
                    f"{chromosome_path}: {row.get('chromosome', 'unknown')} {field} "
                    "does not match standard_plot_config.json"
                )

    artifact_manifest = root / "artifact_manifest.tsv"
    if require_closure_manifest and not artifact_manifest.is_file():
        problems.append(f"{artifact_manifest}: required for reusable-run hash closure")
    if artifact_manifest.is_file():
        manifest_rows = read_tsv(artifact_manifest)
        by_artifact = {row.get("artifact", ""): row for row in manifest_rows}
        for path in (*plot_paths, *table_paths):
            relative = path.relative_to(root).as_posix()
            row = by_artifact.get(relative)
            if row is None:
                problems.append(f"{artifact_manifest}: missing closure row for {relative}")
                continue
            expected_hash = row.get("packaged_sha256") or row.get("sha256")
            if not expected_hash:
                problems.append(f"{artifact_manifest}: missing hash for {relative}")
            elif expected_hash != file_sha256(path):
                problems.append(f"{artifact_manifest}: hash mismatch for {relative}")
    return problems


def validate_sources(
    results_root: Path,
    standard_057: Path,
    standard_058: Path,
    reuse_059: Path | None,
    seed: int,
) -> list[str]:
    problems: list[str] = []
    for spec in SPECS:
        source_root = results_root / spec.logical_id
        for relative in spec.required_files:
            path = source_root / relative
            if not path.is_file():
                problems.append(f"missing required result: {path}")
        for name in spec.required_pdfs:
            path = source_root / "figures" / name
            if not path.is_file():
                problems.append(f"missing required PDF: {path}")
    for label, root in (("057", standard_057), ("058", standard_058)):
        for name in STANDARD_PLOT_NAMES:
            path = root / name
            if not path.is_file():
                problems.append(f"missing historical standard plot for {label}: {path}")
    generator = REPO / "scripts/generate_phase_standard_plots.py"
    if not generator.is_file():
        problems.append(f"missing 059 standard-plot generator: {generator}")
    source_057 = results_root / "057"
    if (source_057 / "run_manifest.tsv").is_file():
        status_057b = assess_057b(source_057)
        problems.extend(f"057B source inconsistency: {item}" for item in status_057b.problems)
    if reuse_059 is not None:
        problems.extend(
            validate_059_standard(reuse_059, seed, require_closure_manifest=True)
        )
    return problems


def provenance_issues(results_root: Path) -> ProvenanceIssues:
    fatal: list[str] = []
    dirty: list[str] = []
    commits: set[str] = set()
    for spec in SPECS:
        path = results_root / spec.logical_id / "provenance.json"
        if not path.is_file():
            continue
        payload = read_json(path)
        environment = str(payload.get("conda_default_env", ""))
        if environment != "analysis":
            fatal.append(f"{path}: conda_default_env={environment!r}, expected 'analysis'")
        recorded_status = payload.get("git_status", [])
        if not isinstance(recorded_status, list):
            recorded_status = []
        source_status = non_output_git_status(str(line) for line in recorded_status)
        if bool(payload.get("git_dirty", True)) and (not recorded_status or source_status):
            detail = "; ".join(source_status[:5]) or "dirty paths were not recorded"
            dirty.append(f"{path}: source-code provenance is dirty ({detail})")
        commit = str(payload.get("git_commit", ""))
        if not commit:
            fatal.append(f"{path}: missing git_commit")
        else:
            commits.add(commit)
    if len(commits) > 1:
        fatal.append("logical results were generated from different git commits: " + ", ".join(sorted(commits)))
    current_commit = git_output("rev-parse", "HEAD")
    if commits and commits != {current_commit}:
        fatal.append(
            f"result commit(s) {sorted(commits)} do not match current HEAD {current_commit}"
        )
    current_status = git_output("status", "--porcelain")
    current_source_status = non_output_git_status(current_status.splitlines())
    if current_source_status:
        dirty.append("current packaging source tree is dirty before staging")
    return ProvenanceIssues(tuple(fatal), tuple(dirty))


def copy_artifact(
    source: Path,
    destination: Path,
    package_root: Path,
    *,
    role: str,
    truth_boundary: str,
    copy_swap_policy: str,
) -> Artifact:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    source_hash = file_sha256(source)
    packaged_hash = file_sha256(destination)
    if source_hash != packaged_hash:
        raise RuntimeError(f"copy hash mismatch: {source} -> {destination}")
    return Artifact(
        artifact=destination.relative_to(package_root).as_posix(),
        role=role,
        source=str(source.resolve()),
        source_sha256=source_hash,
        packaged_sha256=packaged_hash,
        truth_boundary=truth_boundary,
        copy_swap_policy=copy_swap_policy,
    )


def copy_result_files(spec: ExperimentSpec, source_root: Path, package_root: Path) -> list[Artifact]:
    truth_boundary = {
        "057": "mixed diagnostic; every truth-dependent row is explicitly ORACLE",
        "058": "P9016 split uses raw-safe fields; repeated fit/selection metrics are synthetic",
        "059": "synthetic truth for generation/evaluation; nonzero anchors are ORACLE",
    }[spec.logical_id]
    policy = {
        "057": "row-specific; copy0/copy1 are gauge labels",
        "058": "row-specific; seed agreement may use eval-only chromosome flips",
        "059": "row-specific; truth alignment is evaluation/visualization only",
    }[spec.logical_id]
    records: list[Artifact] = []
    machine_files = sorted(
        path
        for path in source_root.rglob("*")
        if path.is_file()
        and path.suffix.lower() in MACHINE_SUFFIXES
        and "figures" not in path.relative_to(source_root).parts
    )
    for source in machine_files:
        relative = source.relative_to(source_root)
        destination = package_root / "tables" / relative
        records.append(
            copy_artifact(
                source,
                destination,
                package_root,
                role=f"logical experiment {spec.logical_id} machine-readable result",
                truth_boundary=truth_boundary,
                copy_swap_policy=policy,
            )
        )
    figure_root = source_root / "figures"
    for source in sorted(figure_root.rglob("*.pdf")):
        destination = package_root / "plots" / source.relative_to(figure_root)
        records.append(
            copy_artifact(
                source,
                destination,
                package_root,
                role=f"logical experiment {spec.logical_id} diagnostic PDF",
                truth_boundary=truth_boundary,
                copy_swap_policy=policy,
            )
        )
    return records


def copy_historical_standard_plots(
    logical_id: str, source_root: Path, package_root: Path
) -> list[Artifact]:
    if logical_id == "057":
        role = "historical experiment 083 native Hickit joint ORACLE-phase geometry plot"
        truth_boundary = "ORACLE_CONTACT_PHASE; historical joint cis+trans run only"
        policy = (
            "per-chromosome geometry-selected copy swap for evaluation; "
            "copy labels are not maternal/paternal"
        )
    elif logical_id == "058":
        role = "historical experiment 078 corrected blind K1 geometry plot"
        truth_boundary = "blind K1 training; CHARM/3DG and phase used post-closure for evaluation only"
        policy = (
            "geometry plots use post-closure geometry-selected swaps; phase metrics use the "
            "declared post-closure SNP-cis swap; neither is blind seed selection"
        )
    else:
        raise ValueError(logical_id)
    return [
        copy_artifact(
            source_root / name,
            package_root / "plots" / name,
            package_root,
            role=role,
            truth_boundary=truth_boundary,
            copy_swap_policy=policy,
        )
        for name in STANDARD_PLOT_NAMES
    ]


def register_generated_artifact(
    path: Path,
    package_root: Path,
    *,
    role: str,
    source: Path,
    truth_boundary: str,
    copy_swap_policy: str,
) -> Artifact:
    if not path.is_file():
        raise FileNotFoundError(f"expected generated artifact is missing: {path}")
    return Artifact(
        artifact=path.relative_to(package_root).as_posix(),
        role=role,
        source=str(source.resolve()),
        source_sha256=file_sha256(source),
        packaged_sha256=file_sha256(path),
        truth_boundary=truth_boundary,
        copy_swap_policy=copy_swap_policy,
    )


def generate_or_reuse_059_standard(
    package_root: Path, seed: int, reuse_root: Path | None
) -> tuple[list[Artifact], str]:
    records: list[Artifact] = []
    truth_boundary = "synthetic truth used for generation/evaluation; no phase anchors in representative run"
    policy = "per-chromosome ORACLE geometry swap, then global Procrustes, evaluation/visualization only"
    if reuse_root is not None:
        for name in STANDARD_PLOT_NAMES:
            records.append(
                copy_artifact(
                    reuse_root / "plots" / name,
                    package_root / "plots" / name,
                    package_root,
                    role="reused representative synthetic standard plot",
                    truth_boundary=truth_boundary,
                    copy_swap_policy=policy,
                )
            )
        for name in ("summary.tsv", "per_chromosome.tsv", "standard_plot_config.json"):
            records.append(
                copy_artifact(
                    reuse_root / name,
                    package_root / name,
                    package_root,
                    role="reused representative synthetic standard evaluation",
                    truth_boundary=truth_boundary,
                    copy_swap_policy=policy,
                )
            )
        problems = validate_059_standard(package_root, seed)
        if problems:
            raise RuntimeError("invalid copied 059 standard run: " + "; ".join(problems))
        return records, f"reused from {reuse_root.resolve()}"

    generator = REPO / "scripts/generate_phase_standard_plots.py"
    command = (
        sys.executable,
        str(generator),
        "--output-dir",
        str(package_root / "plots"),
        "--seed",
        str(seed),
    )
    environment = os.environ.copy()
    environment.setdefault("MPLCONFIGDIR", "/tmp/hickit-phase-formal-matplotlib")
    subprocess.run(command, cwd=REPO, env=environment, check=True)
    for name in STANDARD_PLOT_NAMES:
        records.append(
            register_generated_artifact(
                package_root / "plots" / name,
                package_root,
                role="generated representative synthetic standard plot",
                source=generator,
                truth_boundary=truth_boundary,
                copy_swap_policy=policy,
            )
        )
    for name in ("summary.tsv", "per_chromosome.tsv", "standard_plot_config.json"):
        records.append(
            register_generated_artifact(
                package_root / name,
                package_root,
                role="generated representative synthetic standard evaluation",
                source=generator,
                truth_boundary=truth_boundary,
                copy_swap_policy=policy,
            )
        )
    problems = validate_059_standard(package_root, seed)
    if problems:
        raise RuntimeError("invalid generated 059 standard run: " + "; ".join(problems))
    return records, shell_command(command)


def metric_lookup(
    rows: Iterable[Mapping[str, str]], configuration: str, scope: str, metric: str
) -> str:
    for row in rows:
        if (
            row.get("configuration") == configuration
            and row.get("scope") == scope
            and row.get("metric") == metric
        ):
            return format_float(row.get("value"))
    return "NA"


def provenance_inputs_markdown(source_root: Path) -> str:
    payload = read_json(source_root / "provenance.json")
    hashes = payload.get("input_sha256", {})
    if not isinstance(hashes, dict) or not hashes:
        return "- No external input file; conditions are generated from the packaged configuration."
    return "\n".join(
        f"- `{markdown_escape(path)}`; SHA256 `{markdown_escape(digest)}`"
        for path, digest in sorted(hashes.items())
    )


def mean_metric(rows: Sequence[Mapping[str, str]], field: str) -> str:
    values: list[float] = []
    for row in rows:
        try:
            value = float(row.get(field, "nan"))
        except (TypeError, ValueError):
            continue
        if value == value and abs(value) != float("inf"):
            values.append(value)
    return format_float(sum(values) / len(values)) if values else "NA"


def oracle_mstep_markdown(status: Exp057BStatus) -> str:
    if not status.complete:
        details = "; ".join(status.problems) if status.problems else "explicit full-mode table absent"
        return (
            f"057B execution status: **{markdown_escape(status.declared_status)} (PARTIAL)**. "
            f"The packager found: {markdown_escape(details)}."
        )
    lines = [
        "057B execution status: **COMPLETE_THREE_EXPLICIT_ARMS**. The packager ",
        "cross-validated `run_manifest.tsv`, `oracle_mstep_runs.tsv`, ",
        "`metrics_long.tsv`, and `per_chromosome.tsv` for every declared seed.",
        "",
        "| explicit ORACLE-phase arm | seeds | unordered distance Spearman | homolog-separation Spearman | chromosome-placement Spearman | runtime (s) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    labels = {
        "native_hickit_cis_only": "current CPU Hickit, cis only",
        "native_hickit_cis_plus_trans": "current CPU Hickit, cis + trans",
        "experimental_fixed_cis_rigid_trans_placement": "fixed cis shapes, rigid trans placement",
    }
    for arm in EXPECTED_057B_ARMS:
        rows = [row for row in status.rows if row.get("arm") == arm]
        lines.append(
            "| {label} | {seeds} | {structure} | {separation} | {placement} | {runtime} |".format(
                label=labels[arm],
                seeds=len(rows),
                structure=mean_metric(rows, "mean_copy_swap_invariant_distance_spearman"),
                separation=mean_metric(rows, "mean_homolog_separation_spearman"),
                placement=mean_metric(rows, "chromosome_centroid_distance_spearman"),
                runtime=mean_metric(rows, "runtime_seconds"),
            )
        )
    lines.extend([
        "",
        "Runtime semantics: native Hickit rows report the complete arm run; the ",
        "fixed-cis rigid-trans row reports only its trans placement optimizer and ",
        "excludes the upstream cis-shape reconstruction time.",
    ])
    return "\n".join(lines)


def standard_source_markdown(records: Sequence[Artifact]) -> str:
    lines = []
    for record in records:
        lines.append(
            f"- `{markdown_escape(record.artifact)}` <- `{markdown_escape(record.source)}`; "
            f"SHA256 `{record.source_sha256}`"
        )
    return "\n".join(lines)


def readme_057(
    sequence: int,
    source_root: Path,
    source_hashes: str,
    package_command: str,
) -> str:
    metrics = read_tsv(source_root / "metrics_long.tsv")
    status_057b = assess_057b(source_root)
    oracle_rows = [
        row
        for row in metrics
        if row.get("component") == "057B_historical_joint_oracle_phase_hickit_partial"
        and row.get("metric") == "strict_unordered_pair_spearman_mean"
    ]
    oracle_values = [float(row["value"]) for row in oracle_rows]
    oracle_summary = (
        f"{min(oracle_values):.6f}-{max(oracle_values):.6f} across {len(oracle_values)} historical seeds"
        if oracle_values
        else "NA"
    )
    mstep_summary = oracle_mstep_markdown(status_057b)
    inputs = provenance_inputs_markdown(source_root)
    package_status = "COMPLETE" if status_057b.complete else "PARTIAL"
    return f"""# Experiment {sequence:03d}: logical 057 ORACLE decomposition at 1 Mb

Status: **{package_status}**. This package contains fixed-coordinate E-step,
prior, negative-control, near-truth synthetic, and ORACLE-phase M-step evidence.
The 057B execution status below is derived from machine-readable source tables,
not assigned by the packager.

## Scope and boundary

- Logical experiment: `057`; formal sequence: `{sequence:03d}`.
- Resolution: **1,000,000 bp (1 Mb)**.
- P9016 contacts are used by fixed-coordinate, prior, and negative-control
  diagnostics. SNP phase and CHARM/3DG truth enter only explicitly ORACLE or
  final-evaluation paths.
- `057A` is `ORACLE_COORDINATES`; the imported historical 083 run is
  `ORACLE_CONTACT_PHASE`; near-truth initialization is synthetic and ORACLE.
- Prior and randomized controls do not use SNP phase in inference.
- `copy0/copy1` are gauge labels, never maternal/paternal labels.

Exact source inputs:

{inputs}

## Quantitative snapshot

| scorer | scope | exact top1 | same/cross | exact ECE |
|---|---|---:|---:|---:|
| FDG-flat | all | {metric_lookup(metrics, 'oracle_coords_fdg_flat', 'all', 'exact_top1_accuracy')} | {metric_lookup(metrics, 'oracle_coords_fdg_flat', 'all', 'same_cross_accuracy')} | {metric_lookup(metrics, 'oracle_coords_fdg_flat', 'all', 'exact_ece')} |
| FDG-flat | trans | {metric_lookup(metrics, 'oracle_coords_fdg_flat', 'trans', 'exact_top1_accuracy')} | {metric_lookup(metrics, 'oracle_coords_fdg_flat', 'trans', 'same_cross_accuracy')} | {metric_lookup(metrics, 'oracle_coords_fdg_flat', 'trans', 'exact_ece')} |
| monotone-log | all | {metric_lookup(metrics, 'oracle_coords_monotone_log', 'all', 'exact_top1_accuracy')} | {metric_lookup(metrics, 'oracle_coords_monotone_log', 'all', 'same_cross_accuracy')} | {metric_lookup(metrics, 'oracle_coords_monotone_log', 'all', 'exact_ece')} |
| monotone-log | trans | {metric_lookup(metrics, 'oracle_coords_monotone_log', 'trans', 'exact_top1_accuracy')} | {metric_lookup(metrics, 'oracle_coords_monotone_log', 'trans', 'same_cross_accuracy')} | {metric_lookup(metrics, 'oracle_coords_monotone_log', 'trans', 'exact_ece')} |

Historical 083 joint cis+trans ORACLE-phase unordered structure Spearman range:
`{oracle_summary}`. This is a positive control for one joint native Hickit setup.
It is distinct from the current explicit 057B arms and is not used to infer
their status.

{mstep_summary}

## Standard plot provenance

The standard P9016 geometry plots come from historical experiment 083, not the
new explicit arm table. They use per-chromosome geometry-selected copy swaps for
evaluation. This is ORACLE-phase evidence and does not establish
maternal/paternal identity.

{source_hashes}

## Artifacts and reproduction

Machine-readable logical results are copied under `tables/`; diagnostic PDFs and
standard plots are under `plots/`. `artifact_manifest.tsv` records source and
packaged SHA256 values for every packaged artifact.

Packaging command:

```bash
conda activate analysis
{package_command}
```
"""


def p9016_split_markdown(rows: Sequence[Mapping[str, str]]) -> str:
    selected = [row for row in rows if row.get("source", "").startswith("P9016")]
    lines = [
        "| split mode | N | train | held-out | exact-pair leakage groups | pass |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in selected:
        lines.append(
            "| {mode} | {n} | {train} | {heldout} | {leak} | {passed} |".format(
                mode=markdown_escape(row.get("split_mode", "NA")),
                n=markdown_escape(row.get("n_observations", "NA")),
                train=markdown_escape(row.get("n_train", "NA")),
                heldout=markdown_escape(row.get("n_heldout", "NA")),
                leak=markdown_escape(row.get("exact_haploid_bin_pair_leakage_groups", "NA")),
                passed=markdown_escape(row.get("pass", "NA")),
            )
        )
    return "\n".join(lines)


def readme_058(
    sequence: int,
    source_root: Path,
    source_hashes: str,
    package_command: str,
) -> str:
    split_rows = read_tsv(source_root / "split_manifest.tsv")
    inputs = provenance_inputs_markdown(source_root)
    return f"""# Experiment {sequence:03d}: logical 058 held-out and invariants at 1 Mb

Status: **PARTIAL**. The P9016 raw-safe split audit is packaged. Repeated
held-out fitting, seed stability, and blind model-selection analysis for the
P9016 approved/current baseline are **NOT_RUN**.

## Scope and boundary

- Logical experiment: `058`; formal sequence: `{sequence:03d}`.
- Resolution: **1,000,000 bp (1 Mb)**.
- The P9016 scan reads chromosome/position fields only; phase columns do not
  enter split construction.
- P9016 `readID` is `.` throughout the scanned source, so molecule atomicity is
  `NOT_AUDITABLE` from this pairs file and molecule/readchain split degenerates
  to observation-level grouping.
- New repeated held-out score, seed-selection, and stability rows are synthetic.
  Historical 078 is frozen post-closure context, not a new selected seed.
- `copy0/copy1` remain chromosome-gauge labels.

Exact source inputs:

{inputs}

## Quantitative P9016 split snapshot

{p9016_split_markdown(split_rows)}

These counts audit split membership only. They are not evidence that a P9016
blind-visible held-out metric selects a better reconstruction basin.

## Standard plot provenance

The standard P9016 plots are frozen corrected blind K1 artifacts from historical
experiment 078. Training was blind. Geometry swaps and any phase-based swaps are
post-closure evaluation policies and were not used for seed selection.

{source_hashes}

## Artifacts and reproduction

Machine-readable logical results are copied under `tables/`; diagnostic PDFs and
standard plots are under `plots/`. `artifact_manifest.tsv` records source and
packaged SHA256 values for every packaged artifact.

Packaging command:

```bash
conda activate analysis
{package_command}
```
"""


def readme_059(
    sequence: int,
    source_root: Path,
    package_root: Path,
    standard_method: str,
    package_command: str,
) -> str:
    conditions = read_tsv(source_root / "conditions.tsv")
    families = sorted({row.get("generator_family", "NA") for row in conditions})
    contacts = [int(row["contacts_per_cell"]) for row in conditions if row.get("contacts_per_cell", "").isdigit()]
    anchor_conditions = sum(float(row.get("anchor_fraction", "0")) > 0 for row in conditions)
    summary_rows = read_tsv(package_root / "summary.tsv")
    summary = summary_rows[0] if summary_rows else {}
    inputs = provenance_inputs_markdown(source_root)
    return f"""# Experiment {sequence:03d}: logical 059 synthetic identifiability at 1 Mb

Status: **SYNTHETIC_DIAGNOSTIC_ONLY**. This package is not a P9016 production
training result and does not establish full SNP-free diploid reconstruction.

## Scope and boundary

- Logical experiment: `059`; formal sequence: `{sequence:03d}`.
- Nominal resolution: **1,000,000 bp (1 Mb)**.
- Conditions: `{len(conditions)}`; generator families: `{', '.join(families)}`;
  contact range: `{min(contacts) if contacts else 'NA'}-{max(contacts) if contacts else 'NA'}`.
- Conditions with nonzero phased-anchor fraction: `{anchor_conditions}`. Every
  such condition must be labeled `ORACLE_EXTERNAL_PHASE_ANCHORS` in the machine
  tables; zero-anchor conditions are blind synthetic inference.
- Synthetic truth is used for generation and final evaluation. `copy0/copy1`
  are gauge labels; truth-based geometry swaps are evaluation/visualization only.

Input contract:

{inputs}

## Representative standard evaluation

| metric | value |
|---|---:|
| evaluated non-background, non-anchor contacts | {markdown_escape(summary.get('n_contacts', 'NA'))} |
| exact four-state top1, geometry-gauge aligned (evaluation only) | {format_float(summary.get('exact_top1_accuracy'))} |
| same/cross, geometry-gauge aligned (evaluation only) | {format_float(summary.get('same_cross_accuracy'))} |
| exact ECE, geometry-gauge aligned (evaluation only) | {format_float(summary.get('exact_ece'))} |
| same/cross ECE, geometry-gauge aligned (evaluation only) | {format_float(summary.get('same_cross_ece'))} |
| pmax >= 0.9 coverage, geometry-gauge aligned (evaluation only) | {format_float(summary.get('pmax_ge_0p9_coverage'))} |
| complete synthetic inference runtime seconds | {format_float(summary.get('runtime_seconds'), 3)} |

The representative standard figures were produced as follows: `{markdown_escape(standard_method)}`.
Their per-chromosome copy swaps and global Procrustes alignment are explicitly
evaluation-only and do not demonstrate biological maternal/paternal identity.
Raw arbitrary-gauge contact metrics remain available in the machine-readable
phase-diagram tables and are not interchangeable with the aligned values above.

## Artifacts and reproduction

Machine-readable phase-diagram results are copied under `tables/`; phase-diagram
PDFs and the three standard plots are under `plots/`. The representative
`summary.tsv`, `per_chromosome.tsv`, and `standard_plot_config.json` are at this
folder's top level. `artifact_manifest.tsv` records source and packaged hashes.

Packaging command:

```bash
conda activate analysis
{package_command}
```
"""


def write_package_receipt(
    path: Path,
    *,
    logical_id: str,
    sequence: int,
    timestamp: str,
    source_root: Path,
    command: str,
    source_git_commit: str,
    source_git_status: str,
) -> None:
    payload = {
        "logical_experiment_id": logical_id,
        "formal_sequence": sequence,
        "timestamp_local": timestamp,
        "packaged_at": datetime.now().astimezone().isoformat(),
        "hostname": socket.gethostname(),
        "python": sys.version.replace("\n", " "),
        "conda_default_env": os.environ.get("CONDA_DEFAULT_ENV", ""),
        "git_commit_before_packaging": source_git_commit,
        "git_dirty_before_packaging": bool(source_git_status),
        "git_source_code_dirty_before_packaging": bool(
            non_output_git_status(source_git_status.splitlines())
        ),
        "git_status_before_packaging": source_git_status.splitlines(),
        "packaging_script": str(Path(__file__).resolve()),
        "packaging_script_sha256": file_sha256(Path(__file__).resolve()),
        "source_results": str(source_root.resolve()),
        "source_provenance_sha256": file_sha256(source_root / "provenance.json"),
        "command": command,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_artifact_manifest(path: Path, artifacts: Sequence[Artifact]) -> None:
    fields = (
        "artifact",
        "role",
        "source",
        "source_sha256",
        "packaged_sha256",
        "truth_boundary",
        "copy_swap_policy",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for artifact in sorted(artifacts, key=lambda row: row.artifact):
            writer.writerow(artifact.__dict__)


def register_generated_file(
    path: Path,
    package_root: Path,
    *,
    role: str,
) -> Artifact:
    script = Path(__file__).resolve()
    return Artifact(
        artifact=path.relative_to(package_root).as_posix(),
        role=role,
        source=str(script),
        source_sha256=file_sha256(script),
        packaged_sha256=file_sha256(path),
        truth_boundary="descriptive package metadata only",
        copy_swap_policy="not applicable",
    )


def validate_package(spec: ExperimentSpec, package_root: Path) -> None:
    required = [
        package_root / "README.md",
        package_root / "package_manifest.json",
        package_root / "artifact_manifest.tsv",
        package_root / "plots",
        package_root / "tables",
    ]
    required.extend(package_root / "plots" / name for name in STANDARD_PLOT_NAMES)
    if spec.logical_id == "059":
        required.extend(
            package_root / name
            for name in ("summary.tsv", "per_chromosome.tsv", "standard_plot_config.json")
        )
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise RuntimeError("incomplete staged formal package: " + ", ".join(missing))


def package_one(
    spec: ExperimentSpec,
    sequence: int,
    stage_root: Path,
    results_root: Path,
    standard_057: Path,
    standard_058: Path,
    reuse_059: Path | None,
    seed: int,
    timestamp: str,
    package_command: str,
    source_git_commit: str,
    source_git_status: str,
) -> None:
    source_root = results_root / spec.logical_id
    package_root = stage_root / f"{sequence:03d}-{timestamp}-{spec.slug}"
    (package_root / "plots").mkdir(parents=True)
    (package_root / "tables").mkdir()
    artifacts = copy_result_files(spec, source_root, package_root)

    if spec.logical_id == "057":
        standard_records = copy_historical_standard_plots("057", standard_057, package_root)
        artifacts.extend(standard_records)
        readme = readme_057(
            sequence,
            source_root,
            standard_source_markdown(standard_records),
            package_command,
        )
    elif spec.logical_id == "058":
        standard_records = copy_historical_standard_plots("058", standard_058, package_root)
        artifacts.extend(standard_records)
        readme = readme_058(
            sequence,
            source_root,
            standard_source_markdown(standard_records),
            package_command,
        )
    else:
        standard_records, standard_method = generate_or_reuse_059_standard(
            package_root, seed, reuse_059
        )
        artifacts.extend(standard_records)
        readme = readme_059(
            sequence, source_root, package_root, standard_method, package_command
        )

    readme_path = package_root / "README.md"
    readme_path.write_text(readme, encoding="utf-8")
    receipt_path = package_root / "package_manifest.json"
    write_package_receipt(
        receipt_path,
        logical_id=spec.logical_id,
        sequence=sequence,
        timestamp=timestamp,
        source_root=source_root,
        command=package_command,
        source_git_commit=source_git_commit,
        source_git_status=source_git_status,
    )
    artifacts.append(register_generated_file(readme_path, package_root, role="formal experiment README"))
    artifacts.append(register_generated_file(receipt_path, package_root, role="formal packaging provenance"))
    write_artifact_manifest(package_root / "artifact_manifest.tsv", artifacts)
    validate_package(spec, package_root)


def commit_staged(stage_root: Path, destinations: Mapping[str, Path]) -> None:
    moved: list[tuple[Path, Path]] = []
    try:
        for spec in SPECS:
            destination = destinations[spec.logical_id]
            source = stage_root / destination.name
            if destination.exists():
                raise FileExistsError(f"refusing to overwrite destination created during staging: {destination}")
            source.rename(destination)
            moved.append((source, destination))
    except Exception:
        for source, destination in reversed(moved):
            if destination.exists() and not source.exists():
                destination.rename(source)
        raise


def canonical_execute_command(
    *,
    results_root: Path,
    test_res_root: Path,
    sequence_start: int,
    timestamp: str,
    seed: int,
    standard_057: Path,
    standard_058: Path,
    reuse_059: Path | None,
) -> str:
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--results-root",
        str(results_root),
        "--test-res-root",
        str(test_res_root),
        "--sequence-start",
        str(sequence_start),
        "--timestamp",
        timestamp,
        "--seed",
        str(seed),
        "--exp057-standard-plots",
        str(standard_057),
        "--exp058-standard-plots",
        str(standard_058),
    ]
    if reuse_059 is not None:
        command.extend(("--reuse-059-standard-dir", str(reuse_059)))
    command.append("--execute")
    return shell_command(command)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--test-res-root", type=Path, default=DEFAULT_TEST_RES_ROOT)
    parser.add_argument(
        "--sequence-start",
        type=int,
        help="First of three consecutive formal sequence numbers; default is max existing + 1",
    )
    parser.add_argument(
        "--timestamp",
        help="Shared folder timestamp in YYYYMMDD_HHMMSS; default is current local time",
    )
    parser.add_argument("--seed", type=int, default=17, help="Seed for the representative 059 standard plots")
    parser.add_argument("--exp057-standard-plots", type=Path, default=DEFAULT_057_STANDARD)
    parser.add_argument("--exp058-standard-plots", type=Path, default=DEFAULT_058_STANDARD)
    parser.add_argument(
        "--reuse-059-standard-dir",
        type=Path,
        help="Reuse a validated folder containing plots/ plus summary/per_chromosome/config instead of regenerating",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--execute", action="store_true", help="Create the formal folders")
    mode.add_argument("--dry-run", action="store_true", help="Show and validate the plan without writing (default)")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    require_analysis()
    if args.seed < 0:
        raise ValueError("--seed must be non-negative")
    results_root = args.results_root.resolve()
    test_res_root = args.test_res_root.resolve()
    standard_057 = args.exp057_standard_plots.resolve()
    standard_058 = args.exp058_standard_plots.resolve()
    reuse_059 = args.reuse_059_standard_dir.resolve() if args.reuse_059_standard_dir else None
    timestamp = validate_timestamp(args.timestamp) if args.timestamp else datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
    sequence_start = (
        args.sequence_start
        if args.sequence_start is not None
        else discover_sequence_start(test_res_root)
    )
    destinations = target_paths(test_res_root, sequence_start, timestamp)

    source_problems = validate_sources(
        results_root, standard_057, standard_058, reuse_059, args.seed
    )
    if source_problems:
        raise FileNotFoundError("\n".join(source_problems))
    provenance = provenance_issues(results_root)

    print("Formal phase-result packaging plan")
    print(f"  mode: {'EXECUTE' if args.execute else 'DRY_RUN'}")
    print(f"  results root: {results_root}")
    print(f"  test_res root: {test_res_root}")
    print(f"  timestamp: {timestamp}")
    for spec in SPECS:
        print(f"  logical {spec.logical_id} -> {destinations[spec.logical_id]}")
    print(f"  059 standard plots: {'reuse ' + str(reuse_059) if reuse_059 else 'regenerate with seed ' + str(args.seed)}")
    for problem in provenance.fatal:
        print(f"  FATAL PROVENANCE ERROR: {problem}")
    for problem in provenance.dirty:
        print(f"  NON-OVERRIDABLE DIRTY PROVENANCE: {problem}")

    if not args.execute:
        print("Dry run complete; no files or directories were created.")
        if provenance.fatal:
            print("Execution is blocked by non-overridable provenance errors.")
        if provenance.dirty:
            print("Execution requires clean, committed source-result and packaging provenance.")
        return 0
    if provenance.fatal or provenance.dirty:
        raise RuntimeError(
            "refusing to package non-clean or inconsistent source provenance; "
            "rerun all sources from clean committed code"
        )

    source_git_commit = git_output("rev-parse", "HEAD")
    source_git_status = git_output("status", "--porcelain")
    package_command = canonical_execute_command(
        results_root=results_root,
        test_res_root=test_res_root,
        sequence_start=sequence_start,
        timestamp=timestamp,
        seed=args.seed,
        standard_057=standard_057,
        standard_058=standard_058,
        reuse_059=reuse_059,
    )
    test_res_root.mkdir(parents=True, exist_ok=True)
    lock_path = test_res_root / ".phase-formal-package.lock"
    lock_fd: int | None = None
    owns_lock = False
    stage_root: Path | None = None
    try:
        lock_fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        owns_lock = True
        os.write(lock_fd, f"pid={os.getpid()} host={socket.gethostname()}\n".encode())
        stage_root = Path(tempfile.mkdtemp(prefix=".phase-formal-stage-", dir=test_res_root))
        for offset, spec in enumerate(SPECS):
            package_one(
                spec,
                sequence_start + offset,
                stage_root,
                results_root,
                standard_057,
                standard_058,
                reuse_059,
                args.seed,
                timestamp,
                package_command,
                source_git_commit,
                source_git_status,
            )
        # Recheck sequence occupancy immediately before committing staged folders.
        occupied = occupied_sequences(test_res_root)
        conflicts = [sequence_start + offset for offset in range(len(SPECS)) if sequence_start + offset in occupied]
        if conflicts:
            raise FileExistsError(f"formal sequences became occupied during staging: {conflicts}")
        commit_staged(stage_root, destinations)
        print("Created formal result folders:")
        for spec in SPECS:
            print(f"  {destinations[spec.logical_id]}")
        return 0
    finally:
        if lock_fd is not None:
            os.close(lock_fd)
        if owns_lock and lock_path.exists():
            lock_path.unlink()
        if stage_root is not None and stage_root.exists():
            shutil.rmtree(stage_root)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
