#!/usr/bin/env bash
set -euo pipefail

if [ "${CONDA_DEFAULT_ENV:-}" != "analysis" ]; then
	echo "error: activate conda env 'analysis' before running this script" >&2
	exit 2
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PHASE3_ROOT=$(cd "$REPO_ROOT/.." && pwd)
cd "$REPO_ROOT"

RUN_ID=${HK_BLIND_032_RUN_ID:-"032-$(date +%Y%m%d_%H%M%S)-p9016_trans_contact_scaling_1m"}
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE3_ROOT/test_res}/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
OUTPUT_ROOT="$FULL_RUN_ROOT/outputs"
EVAL_ROOT="$FULL_RUN_ROOT/eval"
LOG_ROOT="$FULL_RUN_ROOT/logs"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"

if [ -e "$FULL_RUN_ROOT" ] &&
   [ "${HK_BLIND_ALLOW_NONEMPTY_RUN_ROOT:-0}" != "1" ] &&
   find "$FULL_RUN_ROOT" -mindepth 1 -print -quit | grep -q .; then
	echo "error: run root already exists and is not empty: $FULL_RUN_ROOT" >&2
	exit 1
fi
mkdir -p "$OUTPUT_ROOT" "$EVAL_ROOT" "$LOG_ROOT" "$FULL_RUN_ROOT/scripts_snapshot"
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)

exec > >(tee -a "$LOG_ROOT/runner.stdout.log") 2> >(tee -a "$LOG_ROOT/runner.stderr.log" >&2)

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

PAIR_PATH=${HK_BLIND_P9016_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
EVAL_PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
APPROVED_P9016_PAIRS=/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-${HK_BLIND_P9016_RELAX_STEP:-0.012}}
BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-gpu}
MAX_PARALLEL=${HK_BLIND_032_MAX_PARALLEL:-${HK_BLIND_MAX_PARALLEL:-10}}
RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-1}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-1}
CONFIG_SET=${HK_BLIND_032_CONFIG_SET:-full}
ALLOW_CUSTOM=${HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS:-0}
MAKE_GPU=${HK_BLIND_MAKE_GPU:-}
if [ -z "$MAKE_GPU" ]; then
	if [ "$BACKEND" = "gpu" ]; then
		MAKE_GPU=1
	else
		MAKE_GPU=0
	fi
fi

case "$MAX_PARALLEL" in
	''|*[!0-9]*)
		echo "error: HK_BLIND_032_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_032_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi

same_realpath() {
	local a=$1
	local b=$2
	local ra rb
	ra=$(realpath "$a")
	rb=$(realpath "$b")
	[ "$ra" = "$rb" ]
}

if [ "$CONFIG_SET" != "smoke" ]; then
	if ! same_realpath "$PAIR_PATH" "$APPROVED_P9016_PAIRS" &&
	   [ "${HK_BLIND_032_ALLOW_CUSTOM_TRAINING_PAIRS:-0}" != "1" ]; then
		echo "error: formal 032 requires approved P9016 training pairs: $APPROVED_P9016_PAIRS" >&2
		echo "got: $PAIR_PATH" >&2
		exit 1
	fi
	if ! same_realpath "$EVAL_PAIRS" "$APPROVED_P9016_PAIRS" &&
	   [ "${HK_BLIND_032_ALLOW_CUSTOM_EVAL_PAIRS:-0}" != "1" ]; then
		echo "error: formal 032 requires approved P9016 eval pairs: $APPROVED_P9016_PAIRS" >&2
		echo "got: $EVAL_PAIRS" >&2
		exit 1
	fi
fi

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

hash_file() {
	if command -v sha256sum >/dev/null 2>&1; then
		sha256sum "$1" | awk '{print $1}'
	else
		echo NA
	fi
}

clean_env_prefix() {
	local env_cmd=(
		env -i
		PATH="$PATH"
		LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
		HOME="${HOME:-}"
		CONDA_DEFAULT_ENV="${CONDA_DEFAULT_ENV:-}"
		CONDA_PREFIX="${CONDA_PREFIX:-}"
	)
	"${env_cmd[@]}" "$@"
}

GIT_COMMIT=$(git rev-parse HEAD 2>/dev/null || echo NA)
GIT_DIRTY_COUNT=$(git status --porcelain | wc -l | awk '{print $1}')
if [ "${HK_BLIND_REQUIRE_CLEAN_TREE:-0}" = "1" ] && [ "$GIT_DIRTY_COUNT" != "0" ]; then
	echo "error: hickit tree is dirty and HK_BLIND_REQUIRE_CLEAN_TREE=1" >&2
	exit 1
fi
git status --porcelain=v1 > "$LOG_ROOT/git_status.txt" || true
git diff --stat > "$LOG_ROOT/git_diff_stat.txt" || true
git diff > "$LOG_ROOT/git_diff.patch" || true
{
	echo "cc	$(${CC:-cc} --version | head -1)"
	echo "make	$(make --version | head -1)"
	echo "cflags	${CFLAGS:-default}"
	echo "conda_env	${CONDA_DEFAULT_ENV:-NA}"
	echo "path	$PATH"
} > "$LOG_ROOT/build_env.txt"

cp "$0" "$FULL_RUN_ROOT/scripts_snapshot/"
cp scripts/p9016_common_eval.sh \
   scripts/p9016_publish_light_result.sh \
   scripts/compute_p9016_copytrack_diag.py \
   scripts/summarize_p9016_model_sweep.py \
   "$FULL_RUN_ROOT/scripts_snapshot/"

log_cmd make smoke_blind_p9016_minimal
if [ "$MAKE_GPU" = "1" ]; then
	log_cmd make gpu=1 run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
else
	log_cmd make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
fi

RUN_BIN="$FULL_RUN_ROOT/logs/run_blind_p9016_minimal.bin"
cp run_blind_p9016_minimal.bin "$RUN_BIN"
chmod +x "$RUN_BIN"
RUN_HASH=$(hash_file "$RUN_BIN")
AUDIT_HASH=$(hash_file audit_blind_p9016_full_cpu_output.bin)
{
	echo "run_blind_p9016_minimal_sha256	$RUN_HASH"
	echo "audit_blind_p9016_full_cpu_output_sha256	$AUDIT_HASH"
} >> "$LOG_ROOT/build_env.txt"

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "run_blind_p9016_minimal_sha256	$RUN_HASH"
	echo "audit_blind_p9016_full_cpu_output_sha256	$AUDIT_HASH"
	echo "pairs	$PAIR_PATH"
	echo "eval_pairs	$EVAL_PAIRS"
	echo "eval_tdg	$EVAL_TDG"
	echo "backend	$BACKEND"
	echo "make_gpu	$MAKE_GPU"
	echo "bin_size_bp	$BIN_SIZE"
	echo "n_iter	$N_ITER"
	echo "relax_steps	$RELAX_STEPS"
	echo "relax_step	$RELAX_STEP"
	echo "max_parallel	$MAX_PARALLEL"
	echo "config_set	$CONFIG_SET"
	echo "training_boundary	P9016 pairs only; phase labels and CHARM/3DG are eval-only"
	echo "experiment_control	posterior_count gamma1 baseline plus trans-only edge k/dscale scaling"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

run_config() {
	local config=$1
	local trans_k=$2
	local trans_dscale=$3
	local min_sep=$4
	local lambda_sep=$5
	local out_dir="$OUTPUT_ROOT/$config"
	local eval_dir="$EVAL_ROOT/$config"
	local log_prefix="$LOG_ROOT/$config"
	local train_cmd audit_cmd eval_cmd copytrack_cmd

	echo "CONFIG_START=$config trans_k=$trans_k trans_dscale=$trans_dscale min_sep=$min_sep lambda_sep=$lambda_sep"
	rm -rf "$out_dir" "$eval_dir"
	mkdir -p "$eval_dir"

	train_cmd=(
		clean_env_prefix
		HK_BLIND_SAMPLE=P9016
		HK_BLIND_P9016_SAMPLE=P9016
		HK_BLIND_GIT_COMMIT="$GIT_COMMIT"
		HK_BLIND_GIT_DIRTY_COUNT="$GIT_DIRTY_COUNT"
		HK_BLIND_BINARY_HASH="$RUN_HASH"
		HK_BLIND_P9016_PAIRS="$PAIR_PATH"
		HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS="$ALLOW_CUSTOM"
		HK_BLIND_P9016_OUTPUT_ROOT="$OUTPUT_ROOT"
		HK_BLIND_P9016_CONFIG_NAME="$config"
		HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE"
		HK_BLIND_P9016_MINIMAL_N_ITER="$N_ITER"
		HK_BLIND_P9016_MINIMAL_RELAX_STEPS="$RELAX_STEPS"
		HK_BLIND_P9016_RELAX_STEP="$RELAX_STEP"
		HK_BLIND_P9016_RELAX_BACKEND="$BACKEND"
		HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split
		HK_BLIND_P9016_INIT_SEED=17
		HK_BLIND_P9016_INIT_EPS=0.5
		HK_BLIND_P9016_INIT_NOISE_SCALE=0
		HK_BLIND_P9016_D_SCALE_MODE=posterior_count
		HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA=1
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6
		HK_BLIND_P9016_MIN_SEP_UNIT="$min_sep"
		HK_BLIND_P9016_LAMBDA_SEP="$lambda_sep"
		HK_BLIND_P9016_LAMBDA_COPYTRACK=0
		HK_BLIND_P9016_LAMBDA_GLOBAL_COPYTRACK=0
		HK_BLIND_P9016_TRANS_CHR_PAIR_PRIOR_LAMBDA=0
		HK_BLIND_P9016_TRANS_K_MULTIPLIER="$trans_k"
		HK_BLIND_P9016_TRANS_DSCALE_MULTIPLIER="$trans_dscale"
		HK_BLIND_WRITE_RAW=0
		"$RUN_BIN"
	)
	log_command_line "${train_cmd[@]}"
	"${train_cmd[@]}" > "$log_prefix.train.log" 2>&1

	audit_cmd=(./audit_blind_p9016_full_cpu_output.bin "$out_dir")
	log_command_line "${audit_cmd[@]}"
	"${audit_cmd[@]}" > "$log_prefix.audit.log" 2>&1

	if [ "$RUN_EVAL" = "1" ]; then
		eval_cmd=(
			clean_env_prefix
			HK_BLIND_SAMPLE=P9016
			HK_BLIND_P9016_SAMPLE=P9016
			HK_BLIND_P9016_EVAL_PAIRS="$EVAL_PAIRS"
			HK_BLIND_P9016_EVAL_TDG="$EVAL_TDG"
			HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE"
			HK_BLIND_REQUIRE_EVAL="$REQUIRE_EVAL"
			scripts/p9016_common_eval.sh "$out_dir" "$eval_dir"
		)
		log_command_line "${eval_cmd[@]}"
		"${eval_cmd[@]}" > "$log_prefix.eval.log" 2>&1
	else
		printf 'status\tSKIPPED_EVAL_DISABLED\n' > "$eval_dir/eval_status.tsv"
	fi

	if [ -s "$out_dir/p9016_full.coords.tsv" ]; then
		copytrack_cmd=(
			python scripts/compute_p9016_copytrack_diag.py
			--config-output-dir "$out_dir"
			--eval-output-dir "$eval_dir"
			--out "$eval_dir/copytrack_vector_diag.tsv"
		)
		log_command_line "${copytrack_cmd[@]}"
		"${copytrack_cmd[@]}" > "$log_prefix.copytrack.log" 2>&1
	fi
	echo "CONFIG_DONE=$config"
}

if [ "$CONFIG_SET" = "smoke" ]; then
	CONFIGS=(
		"p9016_transscale_smoke_tk1_td1|1|1|0|0"
		"p9016_transscale_smoke_tk2_td0p75|2|0.75|0|0"
	)
else
	CONFIGS=(
		"p9016_pcgamma1_tk1_td1_sep0_ieps0p5_noise0_seed17|1|1|0|0"
		"p9016_pcgamma1_tk1p25_td1_sep0_ieps0p5_noise0_seed17|1.25|1|0|0"
		"p9016_pcgamma1_tk1p5_td1_sep0_ieps0p5_noise0_seed17|1.5|1|0|0"
		"p9016_pcgamma1_tk2_td1_sep0_ieps0p5_noise0_seed17|2|1|0|0"
		"p9016_pcgamma1_tk3_td1_sep0_ieps0p5_noise0_seed17|3|1|0|0"
		"p9016_pcgamma1_tk1_td0p9_sep0_ieps0p5_noise0_seed17|1|0.9|0|0"
		"p9016_pcgamma1_tk1_td0p75_sep0_ieps0p5_noise0_seed17|1|0.75|0|0"
		"p9016_pcgamma1_tk1_td0p5_sep0_ieps0p5_noise0_seed17|1|0.5|0|0"
		"p9016_pcgamma1_tk1_td1p25_sep0_ieps0p5_noise0_seed17|1|1.25|0|0"
		"p9016_pcgamma1_tk1p5_td0p75_sep0_ieps0p5_noise0_seed17|1.5|0.75|0|0"
		"p9016_pcgamma1_tk2_td0p75_sep0_ieps0p5_noise0_seed17|2|0.75|0|0"
		"p9016_pcgamma1_tk2_td0p5_sep0_ieps0p5_noise0_seed17|2|0.5|0|0"
		"p9016_pcgamma1_tk3_td0p75_sep0_ieps0p5_noise0_seed17|3|0.75|0|0"
		"p9016_pcgamma1_tk1_td1_msep1_lsep1_ieps0p5_noise0_seed17|1|1|1|1"
		"p9016_pcgamma1_tk2_td0p75_msep1_lsep1_ieps0p5_noise0_seed17|2|0.75|1|1"
	)
fi

failures=0
running=0
for spec in "${CONFIGS[@]}"; do
	IFS='|' read -r config trans_k trans_dscale min_sep lambda_sep <<< "$spec"
	run_config "$config" "$trans_k" "$trans_dscale" "$min_sep" "$lambda_sep" &
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
	echo "one or more configs failed: $failures" >&2
	exit 1
fi

scripts/summarize_p9016_model_sweep.py "$FULL_RUN_ROOT" > "$FULL_RUN_ROOT/summary.tsv"

python - "$FULL_RUN_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))
base = None
for row in rows:
    if row.get("trans_k_multiplier") in ("1", "1.0") and row.get("trans_dscale_multiplier") in ("1", "1.0") and row.get("min_sep_unit") in ("0", "0.0"):
        base = row
        break

def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None

if base is not None:
    b = fnum(base.get("model_top1_accuracy_genome_trans"))
    if b is not None:
        for row in rows:
            t = fnum(row.get("model_top1_accuracy_genome_trans"))
            if t is None:
                continue
            row["delta_trans_vs_baseline"] = f"{t - b:.9g}"
            row["target_plus_0p1_met"] = "1" if t >= b + 0.1 else "0"
for row in rows:
    row["eval_copy_swap_policy"] = "whole_chrom_snp_eval_only"

fields = [
    "config_name", "trans_contact_scaling_mode", "trans_k_multiplier",
    "trans_dscale_multiplier", "min_sep_unit", "lambda_sep",
    "eval_copy_swap_policy",
    "model_top1_accuracy_genome_all", "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans", "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman", "final_mean_entropy", "final_mean_pU",
    "final_min_sep", "sep_p05", "final_mean_sep", "force_diag_total_contact_force_l1",
    "delta_trans_vs_baseline", "target_plus_0p1_met",
]
with (root / "trans_delta_summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: row.get(key, "NA") for key in fields})

headline = "PLUS_0P1_MET" if any(row.get("target_plus_0p1_met") == "1" for row in rows) else "NO_PLUS_0P1"
(root / "headline.txt").write_text(headline + "\n")
PY

python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
light = Path(sys.argv[2])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))

def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("-inf")

rows_sorted = sorted(rows, key=lambda r: fnum(r.get("model_top1_accuracy_genome_trans")), reverse=True)
headline = (root / "headline.txt").read_text().strip()
fields = [
    "config_name", "trans_k_multiplier", "trans_dscale_multiplier",
    "min_sep_unit", "lambda_sep", "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis", "model_top1_accuracy_genome_trans",
    "model_same_cross_accuracy_genome_trans", "mean_per_chrom_cis_distance_spearman",
    "final_mean_entropy", "final_mean_pU", "final_min_sep", "sep_p05", "final_mean_sep",
]
with (root / "README.md").open("w") as out:
    out.write("# 032 P9016 Trans Contact Scaling 1Mb\n\n")
    out.write("This controlled blind-training experiment tests whether trans-specific M-step edge strength or target-distance scaling can improve four-state trans accuracy.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- headline: `{headline}`\n")
    out.write("- training input: raw P9016 pairs only\n")
    out.write("- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper\n")
    out.write("- reported top1 metrics use the standard whole-chromosome SNP copy-swap policy during eval; they are not reference-free model-selection metrics\n")
    out.write("- unchanged model settings: uniform prior, constant rho_train, temperature 1, posterior_count gamma 1, no chromosome-pair prior, no copytrack force\n\n")
    out.write("## Main Results\n\n")
    out.write("| " + " | ".join(fields) + " |\n")
    out.write("|" + "|".join(["---"] * len(fields)) + "|\n")
    for row in rows_sorted:
        out.write("| " + " | ".join(row.get(key, "NA") for key in fields) + " |\n")
    out.write("\n## Interpretation Boundary\n\n")
    out.write("A real success requires full-denominator trans top1 accuracy to increase by at least 0.1 over the no-scaling baseline without a major cis accuracy or cis distance-Spearman collapse. These knobs change only trans M-step contact edges; they do not solve copy-gauge synchronization directly.\n")
PY

scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$FULL_RUN_ROOT/publish.log"

{
	echo "end_time	$(date -Is)"
	echo "headline	$(cat "$FULL_RUN_ROOT/headline.txt")"
} >> "$FULL_RUN_ROOT/run_manifest.tsv"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "TRANS_DELTA_SUMMARY_TSV=$FULL_RUN_ROOT/trans_delta_summary.tsv"
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt")"
