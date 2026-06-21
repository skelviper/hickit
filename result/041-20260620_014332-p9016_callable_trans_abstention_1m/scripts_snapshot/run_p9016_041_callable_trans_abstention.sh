#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PHASE3_ROOT=$(cd "$REPO_ROOT/.." && pwd)

TEST_RES_ROOT=${HK_BLIND_TEST_RES_ROOT:-"$PHASE3_ROOT/test_res"}
RUN_ID=${HK_BLIND_041_RUN_ID:-"041-$(date +%Y%m%d_%H%M%S)-p9016_callable_trans_abstention_1m"}
FULL_RUN_ROOT="$TEST_RES_ROOT/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
DIAG_ROOT="$FULL_RUN_ROOT/diagnostics"
LOG_ROOT="$FULL_RUN_ROOT/logs"
MAX_PARALLEL=${HK_BLIND_041_MAX_PARALLEL:-8}
PRIMARY_SCORE=${HK_BLIND_041_PRIMARY_SCORE:-neg_entropy}
FIXED_COVERAGES=${HK_BLIND_041_FIXED_COVERAGES:-"0.01 0.02 0.05 0.1 0.2 0.3"}

SOURCE_RUN_ROOT=${HK_BLIND_041_SOURCE_RUN_ROOT:-}
if [ -z "$SOURCE_RUN_ROOT" ]; then
	SOURCE_RUN_ROOT=$(find "$TEST_RES_ROOT" -maxdepth 1 -type d -name '040-*-p9016_normdir_copytrack_1m' | sort | tail -n 1)
fi
if [ -z "$SOURCE_RUN_ROOT" ] || [ ! -d "$SOURCE_RUN_ROOT" ]; then
	echo "cannot find source run root; set HK_BLIND_041_SOURCE_RUN_ROOT" >&2
	exit 2
fi
SOURCE_RUN_ROOT=$(cd "$SOURCE_RUN_ROOT" && pwd)

PAIRS=${HK_BLIND_P9016_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}

if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "HK_BLIND_041_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi

mkdir -p "$DIAG_ROOT" "$LOG_ROOT" "$FULL_RUN_ROOT/scripts_snapshot"
cp "$0" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/diagnose_p9016_callable_trans.py" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/p9016_publish_light_result.sh" "$FULL_RUN_ROOT/scripts_snapshot/"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SOURCE_RUN_ROOT=$SOURCE_RUN_ROOT"
echo "PRIMARY_SCORE=$PRIMARY_SCORE"
echo "FIXED_COVERAGES=$FIXED_COVERAGES"

log_cmd() {
	printf '%q ' "$@" >> "$FULL_RUN_ROOT/commands.log"
	printf '\n' >> "$FULL_RUN_ROOT/commands.log"
}

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "source_run_root	$SOURCE_RUN_ROOT"
	echo "pairs	$PAIRS"
	echo "eval_tdg	$EVAL_TDG"
	echo "bin_size_bp	$BIN_SIZE"
	echo "purpose	callable_trans_abstention_audit"
	echo "training_experiment	0"
	echo "primary_score	$PRIMARY_SCORE"
	echo "fixed_coverages	$FIXED_COVERAGES"
	echo "call_scores_use_phase_labels	0"
	echo "call_scores_use_charm_or_reference	0"
	echo "phase_labels_used_eval_only_for_scoring	1"
	echo "charm_used_eval_only_for_shared_denominator	1"
	echo "git_commit	$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo NA)"
	echo "git_dirty_count	$(git -C "$REPO_ROOT" status --short 2>/dev/null | wc -l | awk '{print $1}')"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

git -C "$REPO_ROOT" status --short > "$LOG_ROOT/git_status.txt" || true
git -C "$REPO_ROOT" diff --stat > "$LOG_ROOT/git_diff_stat.txt" || true
git -C "$REPO_ROOT" diff > "$LOG_ROOT/git_diff.patch" || true
{
	echo "timestamp	$(date -Is)"
	echo "shell	$SHELL"
	echo "python	$(command -v python || true)"
	echo "conda	$(command -v conda || true)"
} > "$LOG_ROOT/build_env.txt"

cd "$REPO_ROOT"
if command -v conda >/dev/null 2>&1; then
	eval "$(conda shell.bash hook)"
	conda activate analysis
fi

log_cmd python -m py_compile scripts/diagnose_p9016_callable_trans.py
python -m py_compile scripts/diagnose_p9016_callable_trans.py

CONFIGS=()
while IFS= read -r -d '' eval_dir; do
	config=$(basename "$eval_dir")
	out_dir="$SOURCE_RUN_ROOT/work_outputs/$config/$config"
	if [ ! -d "$out_dir" ]; then
		out_dir="$SOURCE_RUN_ROOT/outputs/$config"
	fi
	if [ -s "$out_dir/p9016_full.manifest.tsv" ] \
		&& [ -s "$out_dir/p9016_full.coords.tsv" ] \
		&& [ -s "$out_dir/p9016_full.bpair_posterior.tsv" ] \
		&& [ -s "$eval_dir/whole_chrom_snp_oracle_swaps.tsv" ]; then
		CONFIGS+=("$config|$out_dir|$eval_dir")
	fi
done < <(find "$SOURCE_RUN_ROOT/eval" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

if [ "${#CONFIGS[@]}" -eq 0 ]; then
	echo "no eligible configs under $SOURCE_RUN_ROOT" >&2
	exit 2
fi

run_one() {
	local spec=$1
	IFS='|' read -r config out_dir eval_dir <<< "$spec"
	local diag_dir="$DIAG_ROOT/$config"
	local log_file="$LOG_ROOT/$config.log"
	rm -rf "$diag_dir"
	mkdir -p "$diag_dir"
	log_cmd python scripts/diagnose_p9016_callable_trans.py --config-output-dir "$out_dir" --eval-output-dir "$eval_dir" --pairs "$PAIRS" --reference-3dg "$EVAL_TDG" --bin-size "$BIN_SIZE" --label "$config" --outdir "$diag_dir"
	python scripts/diagnose_p9016_callable_trans.py \
		--config-output-dir "$out_dir" \
		--eval-output-dir "$eval_dir" \
		--pairs "$PAIRS" \
		--reference-3dg "$EVAL_TDG" \
		--bin-size "$BIN_SIZE" \
		--label "$config" \
		--outdir "$diag_dir" > "$log_file" 2>&1
}

failures=0
running=0
for spec in "${CONFIGS[@]}"; do
	run_one "$spec" &
	running=$((running + 1))
	if [ "$running" -ge "$MAX_PARALLEL" ]; then
		if ! wait -n; then
			failures=$((failures + 1))
		fi
		running=$((running - 1))
	fi
done
while [ "$running" -gt 0 ]; do
	if ! wait -n; then
		failures=$((failures + 1))
	fi
	running=$((running - 1))
done
if [ "$failures" -ne 0 ]; then
	echo "one or more callable diagnostics failed: $failures" >&2
	exit 1
fi

PRIMARY_SCORE="$PRIMARY_SCORE" FIXED_COVERAGES="$FIXED_COVERAGES" python - "$FULL_RUN_ROOT" <<'PY'
from __future__ import annotations

import csv
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

root = Path(sys.argv[1])
diag_root = root / "diagnostics"
primary_score = os.environ.get("PRIMARY_SCORE", "neg_entropy")
fixed_coverages = [float(x) for x in os.environ.get("FIXED_COVERAGES", "0.01 0.02 0.05 0.1 0.2 0.3").split()]


def read_kv(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                out[row[0]] = row[1]
    return out


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def choose(rows: list[dict[str, str]], fraction: str) -> dict[str, str] | None:
    candidates = [r for r in rows if r.get("target_called_fraction") == fraction]
    if not candidates:
        return None
    return max(candidates, key=lambda r: fnum(r.get("top1_accuracy", "nan")))


def fmt(value: object) -> str:
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.9g}"
    return str(value)


def fixed_policy_rows(
    config_name: str,
    records: list[dict[str, str]],
    summary: dict[str, str],
    score_name: str,
    coverages: list[float],
) -> list[dict[str, str]]:
    total_contacts = sum(int(float(r.get("n_contacts", "0"))) for r in records)
    total_bin_pairs = len(records)
    baseline = fnum(summary.get("baseline_top1_accuracy", "nan"))
    ordered = sorted(records, key=lambda r: fnum(r.get(score_name, "nan")), reverse=True)
    rows: list[dict[str, str]] = []
    for coverage in coverages:
        target_contacts = max(1, int(math.ceil(total_contacts * coverage)))
        called: list[dict[str, str]] = []
        called_contacts = 0
        for rec in ordered:
            if called_contacts >= target_contacts:
                break
            called.append(rec)
            called_contacts += int(float(rec.get("n_contacts", "0")))
        top1_correct = sum(int(float(r.get("top1_correct", "0"))) for r in called)
        same_cross_correct = sum(int(float(r.get("same_cross_correct", "0"))) for r in called)
        by_pair: dict[str, int] = defaultdict(int)
        for rec in called:
            by_pair[rec.get("chrom_pair", "NA")] += int(float(rec.get("n_contacts", "0")))
        top_pair = "NA"
        top_pair_contacts = 0
        if by_pair:
            top_pair, top_pair_contacts = max(by_pair.items(), key=lambda item: item[1])
        min_pair_contacts = min(by_pair.values()) if by_pair else 0
        threshold = fnum(called[-1].get(score_name, "nan")) if called else float("nan")
        top1 = top1_correct / called_contacts if called_contacts else float("nan")
        rows.append(
            {
                "config_name": config_name,
                "predeclared_primary_score": score_name,
                "predeclared_coverage_level": fmt(coverage),
                "score_threshold": fmt(threshold),
                "called_contacts": str(called_contacts),
                "total_contacts": str(total_contacts),
                "actual_called_fraction": fmt(called_contacts / total_contacts if total_contacts else float("nan")),
                "called_bin_pairs": str(len(called)),
                "total_bin_pairs": str(total_bin_pairs),
                "called_bin_pair_fraction": fmt(len(called) / total_bin_pairs if total_bin_pairs else float("nan")),
                "top1_accuracy": fmt(top1),
                "delta_top1_vs_config_all_trans": fmt(top1 - baseline if math.isfinite(top1) and math.isfinite(baseline) else float("nan")),
                "same_cross_accuracy": fmt(same_cross_correct / called_contacts if called_contacts else float("nan")),
                "top1_recall": fmt(top1_correct / total_contacts if total_contacts else float("nan")),
                "n_chrom_pairs_called": str(len(by_pair)),
                "top_chrom_pair": top_pair,
                "top_chrom_pair_contacts": str(top_pair_contacts),
                "top_chrom_pair_fraction_of_called_contacts": fmt(top_pair_contacts / called_contacts if called_contacts else float("nan")),
                "min_chrom_pair_contacts_called": str(min_pair_contacts),
                "input_contact_source": summary.get("train_manifest_input_contact_source", "NA"),
                "train_manifest_uses_phase_labels": summary.get("train_manifest_uses_phase_labels", "NA"),
                "train_manifest_uses_charm_or_reference": summary.get("train_manifest_uses_charm_or_reference", "NA"),
                "train_manifest_uses_charm_for_training": summary.get("train_manifest_uses_charm_for_training", "NA"),
                "train_manifest_reference_derived_positive_control": summary.get("train_manifest_reference_derived_positive_control", "NA"),
                "train_manifest_approved_p9016_raw_pairs_realpath": summary.get("train_manifest_approved_p9016_raw_pairs_realpath", "NA"),
                "best_selected_after_eval": "0",
                "selection_uses_truth": "0",
                "selection_uses_charm": "0",
                "scoring_uses_truth": "1",
                "denominator_uses_charm_filter": "1",
            }
        )
    return rows


summary_rows: list[dict[str, str]] = []
curve_rows: list[dict[str, str]] = []
fixed_rows: list[dict[str, str]] = []
for diag in sorted(p for p in diag_root.iterdir() if p.is_dir()):
    s_path = diag / "summary.tsv"
    c_path = diag / "callable_trans_curve.tsv"
    b_path = diag / "callable_trans_best.tsv"
    r_path = diag / "callable_trans_records.tsv"
    if not s_path.exists() or not c_path.exists() or not b_path.exists() or not r_path.exists():
        continue
    summary = read_kv(s_path)
    curves = read_rows(c_path)
    best_rows = read_rows(b_path)
    records = read_rows(r_path)
    curve_rows.extend(curves)
    fixed_rows.extend(fixed_policy_rows(diag.name, records, summary, primary_score, fixed_coverages))
    best = max(best_rows, key=lambda r: fnum(r.get("top1_accuracy", "nan"))) if best_rows else {}
    at10 = choose(curves, "0.1") or {}
    at20 = choose(curves, "0.2") or {}
    at30 = choose(curves, "0.3") or {}
    baseline = fnum(summary.get("baseline_top1_accuracy", "nan"))
    row = {
        "config_name": diag.name,
        "baseline_trans_top1": summary.get("baseline_top1_accuracy", "NA"),
        "baseline_trans_same_cross": summary.get("baseline_same_cross_accuracy", "NA"),
        "best_score_le50": best.get("score_name", "NA"),
        "best_called_fraction_le50": best.get("actual_called_fraction", "NA"),
        "best_trans_top1_le50": best.get("top1_accuracy", "NA"),
        "best_delta_top1_le50": best.get("delta_top1_vs_all_trans", "NA"),
        "best_same_cross_le50": best.get("same_cross_accuracy", "NA"),
        "best_meets_plus0p1_le50": "1" if fnum(best.get("top1_accuracy", "nan")) >= baseline + 0.1 else "0",
        "best_score_at10": at10.get("score_name", "NA"),
        "top1_at10": at10.get("top1_accuracy", "NA"),
        "delta_at10": at10.get("delta_top1_vs_all_trans", "NA"),
        "same_cross_at10": at10.get("same_cross_accuracy", "NA"),
        "best_score_at20": at20.get("score_name", "NA"),
        "top1_at20": at20.get("top1_accuracy", "NA"),
        "delta_at20": at20.get("delta_top1_vs_all_trans", "NA"),
        "same_cross_at20": at20.get("same_cross_accuracy", "NA"),
        "best_score_at30": at30.get("score_name", "NA"),
        "top1_at30": at30.get("top1_accuracy", "NA"),
        "delta_at30": at30.get("delta_top1_vs_all_trans", "NA"),
        "same_cross_at30": at30.get("same_cross_accuracy", "NA"),
        "call_score_uses_phase_labels": summary.get("call_score_uses_phase_labels", "NA"),
        "call_score_uses_charm_or_reference": summary.get("call_score_uses_charm_or_reference", "NA"),
        "eval_only_uses_phase_labels_for_scoring": summary.get("eval_only_uses_phase_labels_for_scoring", "NA"),
        "input_contact_source": summary.get("train_manifest_input_contact_source", "NA"),
        "train_manifest_uses_phase_labels": summary.get("train_manifest_uses_phase_labels", "NA"),
        "train_manifest_uses_charm_or_reference": summary.get("train_manifest_uses_charm_or_reference", "NA"),
    }
    summary_rows.append(row)

fields = [
    "config_name",
    "baseline_trans_top1",
    "baseline_trans_same_cross",
    "best_score_le50",
    "best_called_fraction_le50",
    "best_trans_top1_le50",
    "best_delta_top1_le50",
    "best_same_cross_le50",
    "best_meets_plus0p1_le50",
    "best_score_at10",
    "top1_at10",
    "delta_at10",
    "same_cross_at10",
    "best_score_at20",
    "top1_at20",
    "delta_at20",
    "same_cross_at20",
    "best_score_at30",
    "top1_at30",
    "delta_at30",
    "same_cross_at30",
    "call_score_uses_phase_labels",
    "call_score_uses_charm_or_reference",
    "eval_only_uses_phase_labels_for_scoring",
    "input_contact_source",
    "train_manifest_uses_phase_labels",
    "train_manifest_uses_charm_or_reference",
]
with (root / "summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in sorted(summary_rows, key=lambda r: fnum(r["best_delta_top1_le50"]), reverse=True):
        writer.writerow(row)

all_curve_fields = list(curve_rows[0].keys()) if curve_rows else []
if all_curve_fields:
    with (root / "callable_trans_all_curves.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=all_curve_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(curve_rows)

fixed_fields = [
    "config_name",
    "predeclared_primary_score",
    "predeclared_coverage_level",
    "score_threshold",
    "called_contacts",
    "total_contacts",
    "actual_called_fraction",
    "called_bin_pairs",
    "total_bin_pairs",
    "called_bin_pair_fraction",
    "top1_accuracy",
    "delta_top1_vs_config_all_trans",
    "delta_top1_vs_best035_replay_fixed_policy",
    "same_cross_accuracy",
    "top1_recall",
    "n_chrom_pairs_called",
    "top_chrom_pair",
    "top_chrom_pair_contacts",
    "top_chrom_pair_fraction_of_called_contacts",
    "min_chrom_pair_contacts_called",
    "input_contact_source",
    "train_manifest_uses_phase_labels",
    "train_manifest_uses_charm_or_reference",
    "train_manifest_uses_charm_for_training",
    "train_manifest_reference_derived_positive_control",
    "train_manifest_approved_p9016_raw_pairs_realpath",
    "best_selected_after_eval",
    "selection_uses_truth",
    "selection_uses_charm",
    "scoring_uses_truth",
    "denominator_uses_charm_filter",
]
base_fixed = {
    row["predeclared_coverage_level"]: fnum(row["top1_accuracy"])
    for row in fixed_rows
    if row["config_name"].endswith("best035_replay")
}
for row in fixed_rows:
    base = base_fixed.get(row["predeclared_coverage_level"], float("nan"))
    cur = fnum(row["top1_accuracy"])
    row["delta_top1_vs_best035_replay_fixed_policy"] = fmt(cur - base if math.isfinite(cur) and math.isfinite(base) else float("nan"))
if fixed_rows:
    with (root / "predeclared_callable_trans_summary.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fixed_fields, lineterminator="\n")
        writer.writeheader()
        for row in sorted(
            fixed_rows,
            key=lambda r: (fnum(r["predeclared_coverage_level"]), fnum(r["top1_accuracy"])),
        ):
            writer.writerow({key: row.get(key, "NA") for key in fixed_fields})

best_overall = max(summary_rows, key=lambda r: fnum(r["best_delta_top1_le50"])) if summary_rows else {}
headline = "NO_CALLABLE_PLUS_0P1"
if best_overall and best_overall.get("best_meets_plus0p1_le50") == "1":
    headline = "CALLABLE_SUBSET_PLUS_0P1"
with (root / "headline.txt").open("w") as fh:
    fh.write(headline + "\n")

with (root / "README.md").open("w") as out:
    out.write("# 041 P9016 Callable Trans Abstention Audit\n\n")
    out.write("This is a diagnostic run, not a new training experiment. It asks whether raw-blind confidence scores can select a trans subset whose four-state top1 accuracy exceeds the all-contact trans denominator by at least 0.1.\n\n")
    out.write("## Boundary\n\n")
    out.write("- Call scores use only posterior probabilities, reconstruction geometry, and raw contact count.\n")
    out.write("- Call scores do not use SNP labels, phase labels, CHARM/3DG, or reference structure.\n")
    out.write("- SNP labels are used only after scoring to compute accuracy under the standard whole-chromosome cis SNP gauge.\n")
    out.write("- CHARM/3DG is used only as the same shared-denominator filter used by the standard evaluator.\n\n")
    out.write("## Headline\n\n")
    out.write(f"- `{headline}`\n")
    if best_overall:
        out.write(f"- best config: `{best_overall['config_name']}`\n")
        out.write(f"- all-trans top1: `{best_overall['baseline_trans_top1']}`\n")
        out.write(f"- best called top1 <=50% coverage: `{best_overall['best_trans_top1_le50']}`\n")
        out.write(f"- best delta: `{best_overall['best_delta_top1_le50']}`\n")
        out.write(f"- best called fraction: `{best_overall['best_called_fraction_le50']}`\n\n")
    primary_rows = [
        row for row in fixed_rows
        if row["config_name"].endswith("best035_replay") and row["predeclared_primary_score"] == primary_score
    ]
    if primary_rows:
        out.write("## Predeclared Policy\n\n")
        out.write(f"- primary score: `{primary_score}`\n")
        out.write("- fixed coverage levels: `" + ", ".join(fmt(c) for c in fixed_coverages) + "`\n\n")
        out.write("| coverage | top1 | delta vs all-trans | same/cross | called contacts | called bin pairs | chrom pairs | top chrom-pair fraction |\n")
        out.write("|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in sorted(primary_rows, key=lambda r: fnum(r["predeclared_coverage_level"])):
            out.write(
                f"| {row['predeclared_coverage_level']} | {row['top1_accuracy']} | "
                f"{row['delta_top1_vs_config_all_trans']} | {row['same_cross_accuracy']} | "
                f"{row['called_contacts']} | {row['called_bin_pairs']} | {row['n_chrom_pairs_called']} | "
                f"{row['top_chrom_pair_fraction_of_called_contacts']} |\n"
            )
        out.write("\n")
    out.write("## Top Rows\n\n")
    out.write("| config | all trans top1 | best score | called frac | called top1 | delta | top1 at 10% | top1 at 20% | top1 at 30% |\n")
    out.write("|---|---:|---|---:|---:|---:|---:|---:|---:|\n")
    for row in sorted(summary_rows, key=lambda r: fnum(r["best_delta_top1_le50"]), reverse=True)[:12]:
        out.write(
            f"| `{row['config_name']}` | {row['baseline_trans_top1']} | {row['best_score_le50']} | "
            f"{row['best_called_fraction_le50']} | {row['best_trans_top1_le50']} | {row['best_delta_top1_le50']} | "
            f"{row['top1_at10']} | {row['top1_at20']} | {row['top1_at30']} |\n"
        )
    out.write("\n## Interpretation\n\n")
    if headline == "CALLABLE_SUBSET_PLUS_0P1":
        out.write("- A blind call/no-call objective can produce a trans subset above all-contact trans top1 by more than 0.1.\n")
        out.write("- This does not solve the original full-denominator trans identity problem; it identifies a reliable subset.\n")
        out.write("- The next trainable/reporting step should expose a calibrated trans callability score and predeclare coverage levels.\n")
    else:
        out.write("- These raw-blind confidence scores do not produce a +0.1 callable subset, so abstention alone does not currently rescue trans.\n")
    out.write("- Because 024 found no usable readID/molecule-group anchor in the current pairs file, full-denominator +0.1 likely needs a new raw data source or external/reference information.\n")
PY

log_cmd scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$LOG_ROOT/publish.log"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt")"
