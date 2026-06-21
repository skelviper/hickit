#!/usr/bin/env bash
set -euo pipefail

if [ "${CONDA_DEFAULT_ENV:-}" != "analysis" ]; then
	echo "error: activate conda env 'analysis' before running this script" >&2
	exit 2
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PHASE_ROOT=$(cd "$REPO_ROOT/.." && pwd)
cd "$REPO_ROOT"

RUN_ID=${HK_BLIND_050_RUN_ID:-"050-$(date +%Y%m%d_%H%M%S)-p9016_pairs_identifiability_audit_1m"}
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE_ROOT/test_res}/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
LOG_ROOT="$FULL_RUN_ROOT/logs"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"

mkdir -p "$LOG_ROOT" "$FULL_RUN_ROOT/scripts_snapshot" "$FULL_RUN_ROOT/diagnostics"
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)
exec > >(tee -a "$LOG_ROOT/runner.stdout.log") 2> >(tee -a "$LOG_ROOT/runner.stderr.log" >&2)

APPROVED_PAIRS=/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz
PAIRS=${HK_BLIND_P9016_PAIRS:-$APPROVED_PAIRS}
REFERENCE_3DG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
CONFIG_OUTPUT_DIR=${HK_BLIND_050_CONFIG_OUTPUT_DIR:-/mnt/ssd/zliu/phase3/test_res/049-20260620_103248-p9016_trans_dscale_gamma_1m/outputs/p9016_transgamma_g1_td0p5_seed71_best035_replay/p9016_transgamma_g1_td0p5_seed71_best035_replay}
EVAL_OUTPUT_DIR=${HK_BLIND_050_EVAL_OUTPUT_DIR:-/mnt/ssd/zliu/phase3/test_res/049-20260620_103248-p9016_trans_dscale_gamma_1m/eval/p9016_transgamma_g1_td0p5_seed71_best035_replay}

log_command_line() {
	local line="+"
	local arg
	for arg in "$@"; do
		line+=" $(printf '%q' "$arg")"
	done
	printf '%s\n' "$line" >> "$COMMANDS_LOG"
}

log_cmd() {
	log_command_line "$@"
	"$@"
}

same_realpath() {
	local a=$1 b=$2
	[ "$(realpath "$a")" = "$(realpath "$b")" ]
}

if [ "$BIN_SIZE" != "1000000" ]; then
	echo "error: 050 is fixed at 1Mb" >&2
	exit 2
fi
if ! same_realpath "$PAIRS" "$APPROVED_PAIRS"; then
	echo "error: 050 audits approved P9016 pairs only: $APPROVED_PAIRS" >&2
	echo "got: $PAIRS" >&2
	exit 2
fi
if [ ! -s "$CONFIG_OUTPUT_DIR/p9016_full.bpair_posterior.tsv" ] ||
   [ ! -s "$CONFIG_OUTPUT_DIR/p9016_full.coords.tsv" ] ||
   [ ! -s "$CONFIG_OUTPUT_DIR/p9016_full.manifest.tsv" ]; then
	echo "error: missing source config output files in $CONFIG_OUTPUT_DIR" >&2
	exit 2
fi
if [ ! -s "$EVAL_OUTPUT_DIR/whole_chrom_snp_oracle_swaps.tsv" ]; then
	echo "error: missing eval swaps in $EVAL_OUTPUT_DIR" >&2
	exit 2
fi

GIT_COMMIT=$(git rev-parse HEAD 2>/dev/null || echo NA)
GIT_DIRTY_COUNT=$(git status --porcelain | wc -l | awk '{print $1}')
git status --porcelain=v1 > "$LOG_ROOT/git_status.txt" || true
git diff --stat > "$LOG_ROOT/git_diff_stat.txt" || true
git diff > "$LOG_ROOT/git_diff.patch" || true
{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "pairs	$PAIRS"
	echo "reference_3dg	$REFERENCE_3DG"
	echo "bin_size_bp	$BIN_SIZE"
	echo "config_output_dir	$CONFIG_OUTPUT_DIR"
	echo "eval_output_dir	$EVAL_OUTPUT_DIR"
	echo "experiment_type	posthoc_pairs_only_identifiability_audit"
	echo "training_uses_phase_labels	0"
	echo "training_uses_charm_or_reference	0"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

cp "$0" \
   scripts/diagnose_p9016_pairs_identifiability.py \
   scripts/p9016_publish_light_result.sh \
   "$FULL_RUN_ROOT/scripts_snapshot/"

DIAG_DIR="$FULL_RUN_ROOT/diagnostics/pairs_identifiability"
log_cmd python scripts/diagnose_p9016_pairs_identifiability.py \
	--pairs "$PAIRS" \
	--reference-3dg "$REFERENCE_3DG" \
	--config-output-dir "$CONFIG_OUTPUT_DIR" \
	--eval-output-dir "$EVAL_OUTPUT_DIR" \
	--bin-size "$BIN_SIZE" \
	--outdir "$DIAG_DIR"

cp "$DIAG_DIR/summary.tsv" "$FULL_RUN_ROOT/summary.tsv"

python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
light = Path(sys.argv[2])
summary = {}
with (root / "summary.tsv").open(newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        summary[row["key"]] = row["value"]

curve_rows = list(csv.DictReader((root / "diagnostics/pairs_identifiability/coverage_accuracy_curve.tsv").open(), delimiter="\t"))
blind_rows = [
    r for r in curve_rows
    if r["score_uses_phase_labels"] == "0" and r["score_uses_charm_or_reference"] == "0"
]
blind_rows.sort(key=lambda r: float(r["top1_accuracy"]), reverse=True)
top_blind = blind_rows[:10]

with (root / "README.md").open("w") as out:
    out.write("# 050 P9016 Pairs Identifiability Audit 1Mb\n\n")
    out.write("This is a post-hoc diagnostic, not a new model. It audits whether fields visible in the approved P9016 pairs file and blind-model confidence scores contain enough signal to support a +0.1 full-denominator trans top1 improvement.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- source config output: `{summary.get('config_output_dir','NA')}`\n")
    out.write(f"- baseline trans top1: `{summary.get('baseline_trans_top1_accuracy','NA')}`\n")
    out.write(f"- pairs-visible read linkage available: `{summary.get('pairs_visible_read_linkage_available','NA')}`\n")
    out.write(f"- pairs-visible strand information available: `{summary.get('pairs_visible_strand_information_available','NA')}`\n")
    out.write(f"- best blind score <=50% coverage: `{summary.get('best_blind_score_le_50pct','NA')}`\n")
    out.write(f"- best blind called fraction <=50% coverage: `{summary.get('best_blind_called_fraction_le_50pct','NA')}`\n")
    out.write(f"- best blind top1 accuracy <=50% coverage: `{summary.get('best_blind_top1_accuracy_le_50pct','NA')}`\n")
    out.write(f"- best blind recall <=50% coverage: `{summary.get('best_blind_recall_le_50pct','NA')}`\n\n")
    out.write("## Top Blind Coverage Rows\n\n")
    fields = ["score", "called_fraction", "top1_accuracy", "top1_recall", "same_cross_accuracy", "delta_top1_vs_baseline"]
    out.write("| " + " | ".join(fields) + " |\n")
    out.write("|" + "|".join(["---"] * len(fields)) + "|\n")
    for row in top_blind:
        out.write("| " + " | ".join(row.get(f, "NA") for f in fields) + " |\n")
    out.write("\n## Interpretation\n\n")
    out.write("The approved P9016 pairs file has phase columns, but those are eval-only and cannot enter blind training. The pairs-visible readID/strand fields provide no useful read-level linkage for training. Therefore, any remaining blind selection signal must come from raw count or model confidence. If high accuracy appears only at small coverage and recall stays low, that supports a callable/abstention interpretation rather than a full-denominator +0.1 model improvement.\n")
PY

log_cmd scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"
{
	echo "end_time	$(date -Is)"
} >> "$FULL_RUN_ROOT/run_manifest.tsv"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "DIAG_DIR=$DIAG_DIR"
