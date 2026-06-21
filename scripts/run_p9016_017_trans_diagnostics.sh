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

SOURCE_RUN_ROOT=${HK_BLIND_017_SOURCE_RUN_ROOT:-"$PHASE_ROOT/test_res/016-20260618_150851-p9016_charm_contacts_condition_sweep_1m"}
SOURCE_RADIUS=${HK_BLIND_017_SOURCE_RADIUS:-1pr}
SOURCE_CONFIG=${HK_BLIND_017_SOURCE_CONFIG:-p9016_charm3dg20k_1pr_pcgamma1_common_msep1p5_lsep1}
SOURCE_OUTPUT_DIR=${HK_BLIND_017_SOURCE_OUTPUT_DIR:-"$SOURCE_RUN_ROOT/work_outputs/$SOURCE_RADIUS/$SOURCE_CONFIG/$SOURCE_CONFIG"}
SOURCE_EVAL_DIR=${HK_BLIND_017_SOURCE_EVAL_DIR:-"$SOURCE_RUN_ROOT/eval/$SOURCE_RADIUS/$SOURCE_CONFIG"}

RUN_ID="017-$(date +%Y%m%d_%H%M%S)-p9016_trans_gauge_force_diagnostics"
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE_ROOT/test_res}/$RUN_ID"
mkdir -p "$FULL_RUN_ROOT"/{logs,diagnostics,scripts_snapshot,outputs,eval}
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"

exec > >(tee -a "$FULL_RUN_ROOT/logs/runner.stdout.log") 2> >(tee -a "$FULL_RUN_ROOT/logs/runner.stderr.log" >&2)

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

if [ ! -s "$SOURCE_OUTPUT_DIR/p9016_full.coords.tsv" ]; then
	echo "error: source coords missing: $SOURCE_OUTPUT_DIR/p9016_full.coords.tsv" >&2
	exit 1
fi
if [ ! -s "$SOURCE_OUTPUT_DIR/p9016_full.bpair_posterior.tsv" ]; then
	echo "error: source posterior missing: $SOURCE_OUTPUT_DIR/p9016_full.bpair_posterior.tsv" >&2
	exit 1
fi
if [ ! -s "$SOURCE_EVAL_DIR/summary.tsv" ]; then
	echo "error: source eval summary missing: $SOURCE_EVAL_DIR/summary.tsv" >&2
	exit 1
fi

GIT_COMMIT=$(git rev-parse HEAD)
GIT_DIRTY_COUNT=$(git status --porcelain | wc -l | awk '{print $1}')
git status --porcelain=v1 > "$FULL_RUN_ROOT/logs/git_status.txt"
git diff --stat > "$FULL_RUN_ROOT/logs/git_diff_stat.txt"
git diff > "$FULL_RUN_ROOT/logs/git_diff.patch"
{
	echo "cc	$(${CC:-cc} --version | head -1)"
	echo "make	$(make --version | head -1)"
	echo "conda_env	${CONDA_DEFAULT_ENV:-NA}"
	echo "path	$PATH"
} > "$FULL_RUN_ROOT/logs/build_env.txt"

cp scripts/diagnose_p9016_trans_gauge_geometry.py \
   scripts/diagnose_p9016_force_regime.py \
   scripts/run_p9016_017_trans_diagnostics.sh \
   scripts/p9016_publish_light_result.sh \
   "$FULL_RUN_ROOT/scripts_snapshot/"

log_cmd() {
	local line="+"
	local arg
	for arg in "$@"; do
		line+=" $(printf '%q' "$arg")"
	done
	printf '%s\n' "$line" >> "$COMMANDS_LOG"
	"$@"
}

PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-$(awk -F '\t' '$1=="pairs_path"{print $2}' "$SOURCE_EVAL_DIR/summary.tsv")}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-$(awk -F '\t' '$1=="bin_size_bp"{print $2}' "$SOURCE_EVAL_DIR/summary.tsv")}
if [ ! -s "$PAIRS" ]; then
	echo "error: eval pairs missing: $PAIRS" >&2
	exit 1
fi

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "purpose	trans gauge/geometry and force-regime diagnostics only; no retraining"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "source_run_root	$SOURCE_RUN_ROOT"
	echo "source_radius	$SOURCE_RADIUS"
	echo "source_config	$SOURCE_CONFIG"
	echo "source_output_dir	$SOURCE_OUTPUT_DIR"
	echo "source_eval_dir	$SOURCE_EVAL_DIR"
	echo "pairs_path	$PAIRS"
	echo "bin_size_bp	$BIN_SIZE"
	echo "training_boundary	no new training; reads existing posterior/coords; phase labels used eval-only for diagnostic accuracy"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

ln -s "$SOURCE_OUTPUT_DIR" "$FULL_RUN_ROOT/outputs/$SOURCE_CONFIG"
ln -s "$SOURCE_EVAL_DIR" "$FULL_RUN_ROOT/eval/$SOURCE_CONFIG"

DIAG_DIR="$FULL_RUN_ROOT/diagnostics/$SOURCE_CONFIG"
mkdir -p "$DIAG_DIR"

log_cmd python scripts/diagnose_p9016_trans_gauge_geometry.py \
	--config-output-dir "$SOURCE_OUTPUT_DIR" \
	--eval-dir "$SOURCE_EVAL_DIR" \
	--pairs "$PAIRS" \
	--bin-size "$BIN_SIZE" \
	--outdir "$DIAG_DIR"

log_cmd python scripts/diagnose_p9016_force_regime.py \
	--config-output-dir "$SOURCE_OUTPUT_DIR" \
	--outdir "$DIAG_DIR"

cp "$DIAG_DIR/trans_gauge_geometry_summary.tsv" "$FULL_RUN_ROOT/summary.tsv"

python - "$FULL_RUN_ROOT" "$SOURCE_CONFIG" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
config = sys.argv[2]
diag = root / "diagnostics" / config

def read_rows(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

summary = read_rows(diag / "trans_gauge_geometry_summary.tsv")
force = read_rows(diag / "force_regime_summary.tsv")
genome = next((r for r in summary if r["scope"] == "genome_trans"), {})
force_by_class = {r["contact_class"]: r for r in force}
trans_force = force_by_class.get("trans", {})
cis_force = force_by_class.get("cis", {})

def get(row, key):
    return row.get(key, "NA")

with (root / "README.md").open("w") as fh:
    fh.write("# 017 P9016 Trans Gauge and Force Diagnostics\n\n")
    fh.write("This is a diagnostic-only run. It does not retrain the model. It reads the existing 016 best positive-control posterior/coordinates and uses phase labels only for eval-only trans diagnostics.\n\n")
    fh.write("## Paths\n\n")
    fh.write(f"- full result root: `{root}`\n")
	    fh.write(f"- source config: `{config}`\n")
	    fh.write(f"- diagnostics: `{diag}`\n\n")
	    fh.write("## Main Trans Gauge/Geometry Result\n\n")
	    fh.write("`chr_pair_oracle_*` is an eval-only upper-bound diagnostic: within each trans chromosome pair it chooses the relative copy flip using SNP-truth contacts, so it is not a blind selection metric. `nearest_geometry_state_accuracy` uses the nearest copy-pair state in the reconstruction under the whole-chrom SNP gauge; it is not CHARM/3DG nearest-state accuracy and not posterior-top1 accuracy.\n\n")
	    fh.write("| metric | value |\n|---|---:|\n")
    for key in [
        "n_contacts",
        "whole_chrom_top1_accuracy",
        "chr_pair_oracle_top1_accuracy",
        "chr_pair_oracle_delta_vs_whole",
        "nearest_geometry_state_accuracy",
        "same_cross_accuracy",
        "pmax90_accuracy",
        "pmax90_recall",
        "truth_state_distance_rank_mean",
        "truth_state_is_nearest_fraction",
        "posterior_top1_distance_rank_mean",
        "posterior_top1_is_nearest_fraction",
        "posterior_top1_distance_gap_mean",
    ]:
        fh.write(f"| {key} | {get(genome, key)} |\n")
	    fh.write("\n## Force Regime Summary\n\n")
	    fh.write("These force-regime rows are reconstructed from final bpair posterior state edges. `sum_abs_force` and `force_per_k` use posterior-edge `|fmag|` for regime comparison; they are not expected to equal the runner's exact `contact_force_l1`, which is reported separately in the source `p9016_full.force_class_diag.tsv` with direction-L1 components on the exact training graph.\n\n")
	    fh.write("| class | sum_k | force_per_k | mean_r | frac_k_zero_force | frac_k_tail_attractive | cancellation_ratio |\n")
    fh.write("|---|---:|---:|---:|---:|---:|---:|\n")
    for cls, row in [("cis", cis_force), ("trans", trans_force)]:
        fh.write(
            f"| {cls} | {get(row, 'sum_k')} | {get(row, 'force_per_k')} | "
            f"{get(row, 'mean_r')} | {get(row, 'frac_k_zero_force')} | "
            f"{get(row, 'frac_k_tail_attractive')} | {get(row, 'cancellation_ratio')} |\n"
        )
    fh.write("\n## Interpretation Rules\n\n")
    fh.write("- If chr_pair_oracle_top1_accuracy is much higher than whole_chrom_top1_accuracy, trans is partly a relative chromosome-pair gauge problem.\n")
    fh.write("- If truth_state_is_nearest_fraction remains low, the reconstruction geometry itself does not put SNP-supported trans states nearest.\n")
    fh.write("- Force-regime rows are reconstructed from final bpair posterior state edges, not the unsaved raw softall graph. Use them as a target-distance and cancellation diagnostic, while the source p9016_full.force_class_diag.tsv remains the exact runner graph summary.\n")
    fh.write("- If trans has high force_per_k but high cancellation_ratio or high tail-attractive fraction, the problem is likely long-range force conflict or dscale target semantics rather than absence of trans force.\n")
    fh.write("- This run is not a blind model improvement claim; it is a mechanism diagnostic for choosing the next controlled experiment.\n")
PY

log_cmd scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "DIAGNOSTICS_DIR=$DIAG_DIR"
