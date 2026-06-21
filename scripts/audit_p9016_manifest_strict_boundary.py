#!/usr/bin/env python3
"""Manifest-aware audit for strict P9016 pairs-only trans results."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable


APPROVED_P9016_PAIRS = Path("/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz")
DEFAULT_BASELINE_TRANS = 0.358585


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_kv(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return {}
    start = 1 if rows[0][:2] == ["key", "value"] else 0
    out: dict[str, str] = {}
    for row in rows[start:]:
        if len(row) >= 2:
            out[row[0]] = row[1]
    return out


def safe_float(value: str | None) -> float | None:
    if value is None or value in ("", "NA", "nan", "NaN"):
        return None
    try:
        x = float(value)
    except ValueError:
        return None
    return x if math.isfinite(x) else None


def norm(value: str | None) -> str:
    return "" if value is None else str(value).strip()


def is_zero_or_missing(value: str | None) -> bool:
    v = norm(value)
    if v in ("", "NA", "na", "none", "None", "off", "false", "False"):
        return True
    try:
        return float(v) == 0.0
    except ValueError:
        return False


def approved_pairs_path(manifest: dict[str, str]) -> bool:
    if norm(manifest.get("approved_p9016_raw_pairs_realpath")) == "1":
        return True
    input_path = norm(manifest.get("input_path"))
    if not input_path:
        return False
    try:
        return Path(input_path).resolve() == APPROVED_P9016_PAIRS.resolve()
    except OSError:
        return Path(input_path) == APPROVED_P9016_PAIRS


def iter_summary_rows(test_res_root: Path) -> Iterable[tuple[Path, dict[str, str]]]:
    for summary in sorted(test_res_root.glob("0*-*/summary.tsv")):
        try:
            rows = read_rows(summary)
        except Exception:
            continue
        if not rows or "config_name" not in rows[0]:
            continue
        for row in rows:
            if safe_float(row.get("model_top1_accuracy_genome_trans")) is not None:
                yield summary, row


def find_manifest(summary: Path, row: dict[str, str]) -> Path | None:
    output_dir = norm(row.get("output_dir"))
    config = norm(row.get("config_name"))
    candidates: list[Path] = []
    if output_dir:
        out = Path(output_dir)
        candidates.extend([out / "p9016_full.manifest.tsv", out / config / "p9016_full.manifest.tsv"])
        if out.exists():
            candidates.extend(sorted(out.glob("**/p9016_full.manifest.tsv")))
    run_root = summary.parent
    if config:
        candidates.extend(sorted(run_root.glob(f"outputs/**/{config}/p9016_full.manifest.tsv")))
        candidates.extend(sorted(run_root.glob(f"outputs/{config}/**/p9016_full.manifest.tsv")))
    candidates.extend(sorted(run_root.glob("outputs/**/p9016_full.manifest.tsv")))

    seen: set[Path] = set()
    for path in candidates:
        if not path or path in seen or not path.exists():
            continue
        seen.add(path)
        kv = read_kv(path)
        if config and norm(kv.get("config_name")) not in ("", config):
            continue
        return path
    return None


TRAINING_REFERENCE_FLAG_KEYS = [
    "uses_phase_labels",
    "uses_charm_or_reference",
    "uses_charm_for_training",
    "trans_chr_pair_prior_uses_phase_labels",
    "trans_chr_pair_prior_uses_charm_or_reference",
    "trans_chr_pair_mstep_uses_phase_labels",
    "trans_chr_pair_mstep_uses_charm_or_reference",
    "trans_callable_anchor_uses_phase_labels",
    "trans_callable_anchor_uses_charm_or_reference",
    "trans_contact_scaling_uses_phase_labels",
    "trans_contact_scaling_uses_charm_or_reference",
    "trans_top1_mstep_uses_phase_labels",
    "trans_top1_mstep_uses_charm_or_reference",
    "trans_gate_uses_phase_labels",
    "trans_gate_uses_charm_or_reference",
    "readgroup_uses_phase_labels",
    "readgroup_uses_charm_or_reference",
    "init_coord_anchor_uses_phase_labels",
    "init_coord_anchor_uses_charm_or_reference",
]


def classify(summary: Path, row: dict[str, str]) -> dict[str, str]:
    manifest_path = find_manifest(summary, row)
    manifest = read_kv(manifest_path) if manifest_path else {}
    reasons: list[str] = []

    if not manifest_path:
        reasons.append("manifest_missing")
    if norm(row.get("status")) not in ("", "OK"):
        reasons.append(f"status={norm(row.get('status'))}")
    sample = norm(manifest.get("sample") or row.get("sample"))
    if sample != "P9016":
        reasons.append(f"sample={sample or 'missing'}")
    source = norm(manifest.get("input_contact_source") or row.get("input_contact_source"))
    if source != "raw_pairs":
        reasons.append(f"input_contact_source={source or 'missing'}")
    if source == "raw_pairs" and not approved_pairs_path(manifest):
        reasons.append("approved_p9016_pairs_unverified")
    for key in TRAINING_REFERENCE_FLAG_KEYS:
        if not is_zero_or_missing(manifest.get(key)):
            reasons.append(f"{key}={manifest.get(key)}")
    readgroup = norm(manifest.get("readgroup_mode"))
    if readgroup not in ("", "off", "NA"):
        reasons.append(f"readgroup_mode={readgroup}")
    if not is_zero_or_missing(manifest.get("fixed_posterior")):
        reasons.append(f"fixed_posterior={manifest.get('fixed_posterior')}")

    strict = len(reasons) == 0
    return {
        "run_id": summary.parent.name,
        "config_name": norm(row.get("config_name")),
        "status": norm(row.get("status")) or "NA",
        "strict_pairs_only_candidate": "1" if strict else "0",
        "boundary_reasons": "OK" if strict else ";".join(reasons),
        "model_top1_accuracy_genome_trans": norm(row.get("model_top1_accuracy_genome_trans")),
        "model_top1_accuracy_genome_cis": norm(row.get("model_top1_accuracy_genome_cis")),
        "model_top1_accuracy_genome_all": norm(row.get("model_top1_accuracy_genome_all")),
        "model_same_cross_accuracy_genome_trans": norm(row.get("model_same_cross_accuracy_genome_trans")),
        "mean_per_chrom_cis_distance_spearman": norm(row.get("mean_per_chrom_cis_distance_spearman")),
        "input_contact_source": source or "NA",
        "sample": sample or "NA",
        "uses_phase_labels": norm(manifest.get("uses_phase_labels")) or "NA",
        "uses_charm_or_reference": norm(manifest.get("uses_charm_or_reference")) or "NA",
        "uses_charm_for_training": norm(manifest.get("uses_charm_for_training")) or "NA",
        "readgroup_mode": readgroup or "NA",
        "fixed_posterior": norm(manifest.get("fixed_posterior")) or "NA",
        "manifest_path": str(manifest_path) if manifest_path else "NA",
        "summary_path": str(summary),
    }


AUDIT_FIELDS = [
    "run_id",
    "config_name",
    "status",
    "strict_pairs_only_candidate",
    "boundary_reasons",
    "model_top1_accuracy_genome_trans",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_all",
    "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "input_contact_source",
    "sample",
    "uses_phase_labels",
    "uses_charm_or_reference",
    "uses_charm_for_training",
    "readgroup_mode",
    "fixed_posterior",
    "manifest_path",
    "summary_path",
]


def sort_by_trans(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        rows,
        key=lambda r: safe_float(r.get("model_top1_accuracy_genome_trans")) or float("-inf"),
        reverse=True,
    )


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_plot(outdir: Path, strict: list[dict[str, str]], nonstrict: list[dict[str, str]], baseline: float) -> None:
    plots = outdir / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        (plots / "PLOT_SKIPPED.txt").write_text(f"matplotlib unavailable: {exc}\n")
        return
    rows = sort_by_trans(strict)[:15] + sort_by_trans(nonstrict)[:8]
    if not rows:
        (plots / "PLOT_SKIPPED.txt").write_text("no rows to plot\n")
        return
    labels = [f"{r['run_id']}\n{r['config_name'][:36]}" for r in rows]
    vals = [safe_float(r["model_top1_accuracy_genome_trans"]) or 0.0 for r in rows]
    colors = ["#2563eb" if r["strict_pairs_only_candidate"] == "1" else "#dc2626" for r in rows]
    fig, ax = plt.subplots(figsize=(7.5, max(3.0, 0.28 * len(rows))), dpi=300)
    ax.barh(range(len(rows)), vals, color=colors)
    ax.axvline(baseline, color="#111827", lw=0.8, ls="--", label="baseline")
    ax.axvline(baseline + 0.1, color="#7f1d1d", lw=0.8, ls=":", label="+0.1 target")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_xlabel("trans top1 accuracy", fontsize=7)
    ax.set_title("Manifest-aware trans leaderboard", fontsize=7)
    ax.invert_yaxis()
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(plots / "manifest_strict_trans_leaderboard.png")
    plt.close(fig)


def write_readme(outdir: Path, baseline: float, strict: list[dict[str, str]], nonstrict: list[dict[str, str]], all_rows: list[dict[str, str]]) -> None:
    best_strict = sort_by_trans(strict)[0] if strict else None
    best_control = sort_by_trans(nonstrict)[0] if nonstrict else None
    target = baseline + 0.1
    lines = [
        "# 052 P9016 Manifest-Strict Boundary Audit 1Mb",
        "",
        "This diagnostic audits every numbered `test_res` summary row with a manifest-aware strict P9016 pairs-only boundary.",
        "",
        "## Problem",
        "",
        "A summary-only leaderboard can misclassify positive controls as valid blind candidates when the summary table does not carry all boundary fields. This matters because several reference-derived controls have high trans accuracy and can look close to the +0.1 goal unless their manifests are checked.",
        "",
        "Small conceptual example: a row can report `model_top1_accuracy_genome_trans=0.455`, but its manifest may say `input_contact_source=charm3dg_derived_pairs` and `uses_charm_for_training=1`. That row is useful as a positive control, but it is not a strict P9016.pairs.gz-only training result.",
        "",
        "Validation: this audit joins each summary row to `outputs/**/p9016_full.manifest.tsv`, then requires raw P9016 pairs, sample P9016, no training phase/reference flags, no fixed posterior, and an approved P9016 pairs path.",
        "",
        "## Headline",
        "",
        f"- baseline trans used for +0.1 target: `{baseline:.9g}`",
        f"- target trans: `{target:.9g}`",
        f"- audited rows: `{len(all_rows)}`",
        f"- strict pairs-only rows: `{len(strict)}`",
        f"- non-strict/control rows: `{len(nonstrict)}`",
    ]
    if best_strict:
        trans = safe_float(best_strict["model_top1_accuracy_genome_trans"]) or float("nan")
        lines.extend(
            [
                f"- best strict config: `{best_strict['config_name']}`",
                f"- best strict run: `{best_strict['run_id']}`",
                f"- best strict trans: `{trans:.9g}`",
                f"- best strict delta vs baseline: `{trans - baseline:.9g}`",
                f"- reaches +0.1 target: `{int(trans >= target)}`",
            ]
        )
    if best_control:
        trans = safe_float(best_control["model_top1_accuracy_genome_trans"]) or float("nan")
        lines.extend(
            [
                f"- best excluded control config: `{best_control['config_name']}`",
                f"- best excluded control trans: `{trans:.9g}`",
                f"- exclusion reason: `{best_control['boundary_reasons']}`",
            ]
        )
    lines.extend(
        [
            "",
            "## Top Strict Rows",
            "",
            "| trans | config | run | cis | same/cross trans | cis Spearman |",
            "|---:|---|---|---:|---:|---:|",
        ]
    )
    for row in sort_by_trans(strict)[:12]:
        lines.append(
            f"| {row['model_top1_accuracy_genome_trans']} | `{row['config_name']}` | `{row['run_id']}` | {row['model_top1_accuracy_genome_cis']} | {row['model_same_cross_accuracy_genome_trans']} | {row['mean_per_chrom_cis_distance_spearman']} |"
        )
    lines.extend(
        [
            "",
            "## Top Excluded Controls",
            "",
            "| trans | config | run | reason |",
            "|---:|---|---|---|",
        ]
    )
    for row in sort_by_trans(nonstrict)[:12]:
        lines.append(f"| {row['model_top1_accuracy_genome_trans']} | `{row['config_name']}` | `{row['run_id']}` | {row['boundary_reasons']} |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The high 018 rows are positive controls, not blind candidates: their manifests trace back to CHARM/3DG-derived contacts and fixed posteriors. Under the strict P9016.pairs.gz-only boundary, the best verified full-denominator trans result remains below the requested +0.1 target.",
            "",
            "The useful scientific lesson is that reference-derived clean contacts can nearly reach the target, but the approved pairs-only training surface has not exposed enough blind signal to reproduce that behavior. The next validation step is either a boundary change to an upstream read-level source, or a different objective that reports callable/abstention accuracy and recall.",
        ]
    )
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("test_res_root", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--baseline-trans", type=float, default=DEFAULT_BASELINE_TRANS)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    all_rows = [classify(summary, row) for summary, row in iter_summary_rows(args.test_res_root)]
    strict = [r for r in all_rows if r["strict_pairs_only_candidate"] == "1"]
    nonstrict = [r for r in all_rows if r["strict_pairs_only_candidate"] != "1"]
    write_tsv(args.outdir / "manifest_boundary_audit.tsv", sort_by_trans(all_rows), AUDIT_FIELDS)
    write_tsv(args.outdir / "strict_pairs_only_leaderboard.tsv", sort_by_trans(strict), AUDIT_FIELDS)
    write_tsv(args.outdir / "non_strict_positive_control_leaderboard.tsv", sort_by_trans(nonstrict), AUDIT_FIELDS)

    best_strict = sort_by_trans(strict)[0] if strict else None
    best_control = sort_by_trans(nonstrict)[0] if nonstrict else None
    best_strict_trans = safe_float(best_strict["model_top1_accuracy_genome_trans"]) if best_strict else None
    best_control_trans = safe_float(best_control["model_top1_accuracy_genome_trans"]) if best_control else None
    target = args.baseline_trans + 0.1
    decision = {
        "baseline_trans": f"{args.baseline_trans:.9g}",
        "target_trans_plus_0p1": f"{target:.9g}",
        "n_audited_rows": str(len(all_rows)),
        "n_strict_pairs_only_rows": str(len(strict)),
        "n_non_strict_rows": str(len(nonstrict)),
        "best_strict_config": best_strict["config_name"] if best_strict else "NA",
        "best_strict_run_id": best_strict["run_id"] if best_strict else "NA",
        "best_strict_trans": "NA" if best_strict_trans is None else f"{best_strict_trans:.9g}",
        "best_strict_delta_vs_baseline": "NA" if best_strict_trans is None else f"{best_strict_trans - args.baseline_trans:.9g}",
        "best_strict_meets_plus_0p1": "0" if best_strict_trans is None else str(int(best_strict_trans >= target)),
        "best_non_strict_config": best_control["config_name"] if best_control else "NA",
        "best_non_strict_run_id": best_control["run_id"] if best_control else "NA",
        "best_non_strict_trans": "NA" if best_control_trans is None else f"{best_control_trans:.9g}",
        "best_non_strict_reason": best_control["boundary_reasons"] if best_control else "NA",
        "recommendation": "do_not_continue_pairs_only_knob_sweeps_without_new_blind_signal",
    }
    for name in ("decision.tsv", "summary.tsv"):
        with (args.outdir / name).open("w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(["key", "value"])
            for key, value in decision.items():
                writer.writerow([key, value])

    write_plot(args.outdir, strict, nonstrict, args.baseline_trans)
    write_readme(args.outdir, args.baseline_trans, strict, nonstrict, all_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
