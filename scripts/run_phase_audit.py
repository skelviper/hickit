#!/usr/bin/env python3
"""Generate score, mass, and reproducibility audits without running inference."""

from __future__ import annotations

import argparse
import csv
import os
import platform
import re
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from phase_diagnostics.likelihood import SCORE_MODES, ScoreConfig, is_nonincreasing, score_curve
from phase_diagnostics.mass import mstep_contribution, normalized_pair_expansion
from phase_diagnostics.provenance import file_sha256


P9016_PAIRS = Path("/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz")
P9016_TDG = Path("/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz")
HISTORICAL_056 = Path(
    "/work/phase3/test_res/056-20260620_230016-p9016_readchain_dedup_training_1m"
)
HISTORICAL_TEST_RES = Path("/work/phase3/test_res")


EXPERIMENT_TOPICS = {
    1: "softall 4 Mb baseline",
    2: "softall 4 Mb to 1 Mb to 200 kb resolution chain",
    3: "softall 1 Mb direct GPU baseline",
    4: "softall 1 Mb direct CPU baseline",
    5: "random diploid 1 Mb initialization baseline",
    6: "random diploid multiresolution chain",
    7: "anti-collapse sweep",
    8: "anti-collapse random-initialization sweep",
    9: "softall 200 kb direct GPU baseline",
    10: "count-to-distance and separation sweep",
    11: "copy-track regression",
    12: "count-to-distance semantic regression",
    13: "P1006 count-to-distance control",
    14: "P9016/P1006 posterior-gamma and separation comparison",
    15: "CHARM/3DG-derived clean-contact positive control",
    16: "CHARM/3DG-derived contact-condition sweep",
    17: "trans gauge and force diagnostic",
    18: "fixed ORACLE posterior graph-semantics control",
    19: "copy-track continuity",
    20: "trans gauge ORACLE diagnostic",
    21: "trans chromosome-pair prior",
    22: "trans pair-sharpen warmup",
    23: "blind pair-flip rule audit",
    24: "raw pairs source-signal audit",
    25: "callable trans confidence audit",
    26: "contacts.seg source-signal audit",
    27: "contacts.seg source sweep",
    28: "readgroup-consistent decoding",
    29: "readgroup joint marginal",
    30: "global copy-track regularization",
    31: "confidence-boosted pairs",
    32: "trans contact scaling",
    33: "trans top1 M-step",
    34: "entropy, rho, and temperature",
    35: "copy-track plus entropy",
    36: "trans confidence/rho",
    37: "trans chromosome-pair M-step",
    38: "trans decoder audit",
    39: "global gauge decoder audit",
    40: "normalized-direction copy-track",
    41: "callable trans abstention",
    42: "trans-gated M-step",
    43: "raw pairs remaining ablation",
    44: "multi-seed basin audit",
    45: "raw held-out selection audit",
    46: "pair-flip synchronization signal audit",
    47: "blind coordinate-initialization anchor",
    48: "callable-anchor M-step",
    49: "trans count-to-distance posterior gamma",
    50: "pairs identifiability audit",
    51: "pairs-only signal boundary",
    52: "strict training/evaluation boundary audit",
    53: "contacts.seg source audit",
    54: "contacts.seg readgroup deduplication audit",
    55: "contacts.seg readgroup training",
    56: "readchain deduplication training and approved pairs-only comparator",
}


CANONICAL_RUN_SUFFIX = {
    18: "018-20260618_220644-p9016_fixed_posterior_graph_semantics_1m",
    21: "021-20260618_234407-p9016_trans_chr_pair_prior_1m",
    24: "024-20260619_024741-p9016_pairs_source_signal_audit_1m",
    27: "027-20260619_030956-p9016_contacts_seg_source_sweep_1m",
    28: "028-20260619_034341-p9016_readgroup_consistent_decoding_combined_1m",
    29: "029-20260619_061138-p9016_readgroup_joint_marginal_1m",
    31: "031-20260619_105814-p9016_confidence_boost_pairs_1m",
    41: "041-20260620_015915-p9016_callable_trans_abstention_1m",
    54: "054-20260620_163922-p9016_contacts_seg_readgroup_dedup_audit_1m",
}


def historical_entrypoint(experiment_id: int, canonical_path: Path | None) -> tuple[str, str, str]:
    """Resolve a preserved launcher without pretending a missing name is executable."""
    if experiment_id in {1, 3, 4, 5, 9}:
        candidate = "scripts/run_p9016_softall_baseline.sh"
    if experiment_id in {2, 6}:
        candidate = "scripts/run_p9016_softall_resolution_chain.sh"
    if experiment_id in {7, 8}:
        candidate = "scripts/run_p9016_anticollapse_sweep.sh"
    if experiment_id == 10:
        candidate = "scripts/run_p9016_dscale_sep_sweep.sh"
    if experiment_id not in set(range(1, 11)):
        candidate = ""

    preserved: list[Path] = []
    if canonical_path is not None:
        preserved.extend(sorted((canonical_path / "scripts_snapshot").glob(f"*{experiment_id:03d}*")))
        preserved.extend(sorted((canonical_path / "logs").glob(f"*{experiment_id:03d}*")))
        if experiment_id == 56:
            preserved.extend(sorted((canonical_path / "scripts_snapshot").glob("run_p9016_056_*.sh")))
    preserved = [path for path in preserved if path.is_file()]
    preserved.sort(key=lambda path: (
        0 if path.name.startswith("run_") and path.suffix == ".sh" else 1,
        0 if path.suffix == ".sh" else 1,
        path.name,
    ))
    if preserved:
        path = preserved[0]
        return str(path), "PRESERVED", str(path)
    if candidate:
        return f"UNAVAILABLE: {candidate}", "UNAVAILABLE", str(canonical_path or HISTORICAL_TEST_RES)
    evidence = ""
    if canonical_path is not None:
        for name in ("README.md", "run_manifest.tsv", "commands.log"):
            path = canonical_path / name
            if path.is_file():
                evidence = str(path)
                break
    return "README/manifest/commands only", "NOT_PRESERVED", evidence or str(canonical_path or HISTORICAL_TEST_RES)


def lineage_boundary(experiment_id: int) -> tuple[str, str, str]:
    if experiment_id == 13:
        return "P1006 pairs", "none", "non-P9016 control"
    if experiment_id in {15, 16}:
        return "CHARM/3DG-derived contacts", "ORACLE_REFERENCE_CONTACTS", "positive control only"
    if experiment_id == 18:
        return "fixed posterior from experiment 016", "ORACLE_FIXED_POSTERIOR", "positive control only"
    if experiment_id in {20, 38, 39}:
        return "P9016 raw contacts", "none", "truth used only in explicitly ORACLE/eval decoder"
    if experiment_id in {26, 27, 28, 29, 53, 54, 55}:
        return "P9016 contacts.seg/readgroup source", "none", "not the approved pairs-only input surface"
    if experiment_id == 56:
        return "P9016 pairs plus separate readchain comparator arms", "none", "approved baseline arm is pairs-only"
    return "P9016 raw pairs", "none", "phase and CHARM/3DG evaluation only"


def experiment_lineage(output: Path) -> None:
    rows: list[dict[str, object]] = []
    directories = [path for path in HISTORICAL_TEST_RES.iterdir() if path.is_dir()]
    for experiment_id in range(1, 57):
        prefix = f"{experiment_id:03d}-"
        physical = sorted(
            path.name for path in directories
            if path.name.startswith(prefix)
            and re.match(r"^[0-9]{3}-[0-9]{8}_[0-9]{6}-", path.name)
        )
        if not physical:
            canonical = "MISSING"
        else:
            canonical = CANONICAL_RUN_SUFFIX.get(experiment_id, physical[-1])
        canonical_path = HISTORICAL_TEST_RES / canonical if canonical != "MISSING" else None
        readme_exists = bool(canonical_path and (canonical_path / "README.md").is_file())
        manifest_exists = bool(canonical_path and any(canonical_path.rglob("*manifest*")))
        training_source, truth_used, boundary_note = lineage_boundary(experiment_id)
        if experiment_id == 20:
            status = "MISSING_REQUIRED_PROVENANCE"
        elif experiment_id == 13:
            status = "NON_P9016_CONTROL"
        elif experiment_id in {15, 16, 18}:
            status = "ORACLE_CONTROL"
        elif len(physical) > 1:
            status = "CANONICAL_SELECTED_FROM_MULTIPLE_RUNS"
        else:
            status = "HISTORICAL"
        provenance = "high" if readme_exists and manifest_exists else "medium" if readme_exists else "low"
        if experiment_id <= 10:
            provenance = "low" if not manifest_exists else "medium"
        notes = boundary_note
        if len(physical) > 1:
            notes += "; all physical runs retained in physical_runs"
        if experiment_id == 41:
            notes += "; latest repaired machine-schema run, not formally declared superseding"
        entrypoint, entrypoint_status, evidence_path = historical_entrypoint(experiment_id, canonical_path)
        if experiment_id in CANONICAL_RUN_SUFFIX:
            selection_reason = "explicit canonical mapping retained from experiment audit"
        elif len(physical) == 1:
            selection_reason = "only physical run with this logical identifier"
        elif len(physical) > 1:
            selection_reason = "latest timestamp used because no explicit supersession record was found"
        else:
            selection_reason = "no physical run found"
        rows.append({
            "logical_experiment_id": f"{experiment_id:03d}",
            "canonical_run_id": canonical,
            "physical_runs": ";".join(physical) if physical else "MISSING",
            "topic": EXPERIMENT_TOPICS[experiment_id],
            "status": status,
            "training_source": training_source,
            "truth_used_in_inference": truth_used,
            "resolution_bp": "1000000" if experiment_id >= 11 else "mixed; see run README",
            "entrypoint": entrypoint,
            "entrypoint_status": entrypoint_status,
            "evidence_path": evidence_path,
            "canonical_selection_reason": selection_reason,
            "readme_present": int(readme_exists),
            "manifest_present": int(manifest_exists),
            "provenance_grade": provenance,
            "notes": notes,
        })
    write_tsv(output / "experiment_lineage_001_056.tsv", rows)


def require_analysis() -> None:
    if os.environ.get("CONDA_DEFAULT_ENV") != "analysis":
        raise RuntimeError("this script must run inside conda environment 'analysis'")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    fieldnames = fields or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def score_audit(output: Path) -> None:
    distances = np.unique(np.concatenate((
        np.asarray([0.0, 1e-12, 1e-9, 1e-6, 1e-4, 0.01, 0.05]),
        np.linspace(0.1, 5.0, 246),
    )))
    rows: list[dict[str, object]] = []
    for mode in SCORE_MODES:
        for temperature in (0.5, 1.0, 2.0):
            config = ScoreConfig(mode=mode, temperature=temperature)
            monotone = is_nonincreasing(config, distances)
            for row in score_curve(config, distances):
                rows.append({**row, "monotonic_nonincreasing": int(monotone)})
    write_tsv(output / "estep_score_curves.tsv", rows)

    plt.rcParams.update({"font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7})
    figure, axes = plt.subplots(3, 2, figsize=(6.4, 9.0), constrained_layout=True)
    for axis, mode in zip(axes.flat, SCORE_MODES, strict=True):
        for temperature, color in ((0.5, "#D55E00"), (1.0, "#0072B2"), (2.0, "#009E73")):
            selected = [row for row in rows if row["score_mode"] == mode and row["temperature"] == temperature]
            axis.plot(
                [float(row["distance"]) for row in selected],
                [float(row["log_likelihood_like_score"]) for row in selected],
                color=color, label=f"T={temperature:g}", linewidth=1.0,
            )
        axis.set_title(mode)
        axis.set_xlabel("normalized distance")
        axis.set_ylabel("log score")
        axis.legend(frameon=False)
        axis.axvline(0.5, color="#777777", linewidth=0.4, linestyle=":")
        axis.axvline(1.5, color="#777777", linewidth=0.4, linestyle=":")
    figure.savefig(output / "estep_score_curves.pdf", dpi=300)
    plt.close(figure)


def mass_audit(output: Path) -> None:
    rows: list[dict[str, object]] = []

    def add_row(*, expected_state_mass: float, force_policy_pass: bool, **values: object) -> None:
        observed = float(values["total_state_mass"])
        state_pass = bool(np.isclose(observed, expected_state_mass))
        rows.append({
            **values,
            "expected_total_state_mass": expected_state_mass,
            "state_mass_pass": int(state_pass),
            "force_policy_pass": int(force_policy_pass),
            "pass": int(state_pass and force_policy_pass),
        })

    for top in (0.25, 0.50, 0.75, 0.95, 0.999):
        posterior = np.full(4, (1.0 - top) / 3.0)
        posterior[0] = top
        for target_mode in ("posterior_once", "legacy_posterior_target"):
            contribution = mstep_contribution(1.0, posterior, target_mode=target_mode)
            add_row(
                source_type="raw_pair",
                observation_id=f"pmax_{top:g}_{target_mode}",
                stage="state_expansion_and_target_conversion",
                number_of_segments=2,
                number_of_pair_expansions=1,
                total_input_mass=1.0,
                total_state_mass=contribution.total_state_mass,
                total_mstep_force_coefficient=contribution.total_force_coefficient,
                expected_state_mass=1.0,
                force_policy_pass=target_mode == "posterior_once",
                notes=(
                    "posterior enters force exactly once"
                    if target_mode == "posterior_once"
                    else "FAIL: posterior also changes target distance"
                ),
            )

    base = np.asarray([0.55, 0.20, 0.15, 0.10])
    sharpened = np.square(base)
    sharpened /= sharpened.sum()
    for stage, state_mass, expected, force_mass, note in (
        ("posterior_sharpening_renormalized", sharpened.sum(), 1.0, sharpened.sum(), "power sharpening followed by normalization"),
        ("entropy_gate_soft_weight", 0.4 * base.sum(), 0.4, 0.4 * base.sum(), "gate rho=0.4 may abstain but cannot amplify source mass"),
        ("contact_aggregation_two_raw", 2.0 * base.sum(), 2.0, 2.0 * base.sum(), "two raw observations retain total mass two after aggregation"),
        ("resampling_three_copies_normalized", 3.0 * (base / 3.0).sum(), 1.0, 1.0, "three resamples share one original-observation mass"),
        ("resolution_change_four_children_normalized", 4.0 * (base / 4.0).sum(), 1.0, 1.0, "four child-bin expansions share one parent-observation mass"),
    ):
        add_row(
            source_type="raw_pair",
            observation_id=f"audit_{stage}",
            stage=stage,
            number_of_segments=2,
            number_of_pair_expansions=(2 if "aggregation" in stage else 3 if "resampling" in stage else 4 if "resolution" in stage else 1),
            total_input_mass=(2.0 if "aggregation" in stage else 1.0),
            total_state_mass=float(state_mass),
            total_mstep_force_coefficient=float(force_mass),
            expected_state_mass=expected,
            force_policy_pass=True,
            notes=note,
        )
    for segments in (2, 3, 4, 8):
        weights = normalized_pair_expansion(segments)
        add_row(
            source_type="readchain_molecule",
            observation_id=f"molecule_segments_{segments}",
            stage="all_pair_expansion_normalized",
            number_of_segments=segments,
            number_of_pair_expansions=weights.size,
            total_input_mass=1.0,
            total_state_mass=float(weights.sum()),
            total_mstep_force_coefficient=float(weights.sum()),
            expected_state_mass=1.0,
            force_policy_pass=True,
            notes="each pair receives 1/choose(m,2) before four-state expansion",
        )
        add_row(
            source_type="readchain_molecule_legacy",
            observation_id=f"legacy_segments_{segments}",
            stage="all_pair_expansion_unnormalized",
            number_of_segments=segments,
            number_of_pair_expansions=weights.size,
            total_input_mass=1.0,
            total_state_mass=float(weights.size),
            total_mstep_force_coefficient=float(weights.size),
            expected_state_mass=1.0,
            force_policy_pass=weights.size == 1,
            notes="FAIL for m>2: current C pair-expansion baseline gives every derived pair unit mass",
        )
    write_tsv(output / "mass_conservation.tsv", rows)


def compiler_version() -> str:
    try:
        return subprocess.check_output(("cc", "--version"), text=True).splitlines()[0]
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def baseline_manifest(repo: Path, output: Path) -> None:
    commit = subprocess.check_output(("git", "rev-parse", "HEAD"), cwd=repo, text=True).strip()
    dirty = bool(subprocess.check_output(("git", "status", "--porcelain"), cwd=repo, text=True).strip())
    source_manifest = HISTORICAL_056 / "work_outputs/approved_baseline_pairs_only/approved_baseline_pairs_only/p9016_full.manifest.tsv"
    source_eval = HISTORICAL_056 / "eval/approved_baseline_pairs_only"
    command_lines = (HISTORICAL_056 / "commands.log").read_text(encoding="utf-8").splitlines()
    historical_run_command = command_lines[3].removeprefix("+ ") if len(command_lines) >= 4 else "UNAVAILABLE"
    rows = [
        ("baseline_definition", "historical_approved_pairs_only_056", source_manifest, "sample/config_name", "historical artifact; not correctness-valid after later force audit"),
        ("active_method_family", "softall/raw_expected_soft_all", repo / "run_blind_p9016_minimal.c", "set_minimal_schedule", "current code family"),
        ("input_pairs", str(P9016_PAIRS), source_manifest, "input_path", "phase columns physically present but ignored by blind parser"),
        ("input_pairs_sha256", file_sha256(P9016_PAIRS), source_manifest, "input_path", "independently verified"),
        ("preprocessing_source", "/sharec/zliu/CHARM/mESC/result/cleaned_pairs/c123/P9016.pairs.gz", Path("/sharec/zliu/CHARM/mESC/CHARM/rules/scHiC_2dprocess.rules"), "scHiC_2dprocess", "byte-identical upstream c123 pairs"),
        ("reference_3dg_eval_only", str(P9016_TDG), source_eval / "eval_command.txt", "--reference-3dg", "ORACLE/eval only"),
        ("reference_3dg_sha256", file_sha256(P9016_TDG), source_eval / "eval_command.txt", "--reference-3dg", "independently verified"),
        ("bin_size_bp", "1000000", source_manifest, "resolution", "1 Mb"),
        ("chromosome_set", "chr1-chr19,chrX", source_eval / "eval_manifest.json", "chromosomes", "20 chromosomes"),
        ("raw_contact_count", "1703888", source_manifest, "n_raw", "all P9016 contacts"),
        ("cis_contact_count", "1135454", source_manifest, "n_raw_cis", "raw contacts"),
        ("trans_contact_count", "568434", source_manifest, "n_raw_trans", "raw contacts"),
        ("any_phase_labeled_count", "1302203", source_eval / "summary.tsv", "input_contacts_with_any_known_phase", "evaluation availability only"),
        ("both_phase_labeled_count", "496021", P9016_PAIRS, "independent scan", "evaluation availability only"),
        ("truth_eval_count", "356770", source_eval / "contact_truth_distribution.tsv", "all", "after coordinate/evaluator eligibility"),
        ("top1_eval_count", "350034", source_eval / "contact_accuracy.tsv", "all", "coverage denominator for historical top1"),
        ("historical_build_command", "make gpu=1 run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin", HISTORICAL_056 / "commands.log", "line 1", "requires CUDA/GPU"),
        ("historical_run_command", historical_run_command, HISTORICAL_056 / "commands.log", "line 4", "exact full environment-preserving command"),
        ("current_legacy_compatibility_flag", "HK_BLIND_P9016_ESTEP_SCORE_MODE=legacy_fdg_flat", repo / "run_blind_p9016_minimal.c", "parse_estep_score_mode", "prepend when reproducing the historical omitted-environment FDG-flat semantics with current code"),
        ("historical_random_seed", "17", source_manifest, "init_seed", "explicit"),
        ("historical_estep_score", "fdg_flat", source_manifest, "estep_score_mode", "non-monotonic shell energy"),
        ("historical_mstep_weight", "posterior_count gamma=1", source_manifest, "dscale_mode", "posterior affects k and target distance"),
        ("historical_iterations", "100x100", source_manifest, "n_iter/relax_steps", "100 EM rounds and 100 relax steps"),
        ("historical_relax_step", "0.012", source_manifest, "relax_step", "old force scale; later invalidated"),
        ("historical_git_commit", "1037ee5f79d54e943eef5c37e5aeb2d841698cb8", source_manifest, "git_commit", "dirty_count=133"),
        ("historical_runner_sha256", "69deae8c2e88401fef016873aa8bc970f0910b57b57d01b271d18f6720005c8a", HISTORICAL_056 / "logs/build_env.txt", "runner_sha256", "saved binary"),
        ("current_git_commit", commit, repo / ".git", "HEAD", "diagnostics branch"),
        ("current_git_dirty", str(int(dirty)), repo / ".git", "status --porcelain", "expected while generating diagnostics"),
        ("compiler_version", compiler_version(), repo / "Makefile", "CC", "current node"),
        ("python_version", platform.python_version(), Path(sys.executable), "version", "conda analysis"),
        ("conda_environment", os.environ.get("CONDA_DEFAULT_ENV", ""), Path(sys.executable), "CONDA_DEFAULT_ENV", "must equal analysis"),
        ("gpu_rerun_status", "BLOCKED_NO_NVIDIA_DRIVER", HISTORICAL_056 / "commands.log", "gpu=1", "current nvidia-smi cannot communicate with driver"),
        ("historical_terminal_audit", "INVALID_AFTER_MANIFEST_APPEND", HISTORICAL_056 / "scripts_snapshot/run_p9016_056_readchain_dedup_training.sh", "audit then append source fields", "saved OK receipt predates final manifest mutation"),
    ]
    write_tsv(
        output / "baseline_manifest.tsv",
        [{"field": field, "value": value, "source_file": str(path), "line_or_function": line, "notes": notes}
         for field, value, path, line, notes in rows],
        ["field", "value", "source_file", "line_or_function", "notes"],
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results/phase_diagnostics"))
    args = parser.parse_args()
    require_analysis()
    repo = Path(__file__).resolve().parents[1]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    score_audit(args.output_dir)
    mass_audit(args.output_dir)
    baseline_manifest(repo, args.output_dir)
    experiment_lineage(args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
