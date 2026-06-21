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

RUN_ID=${HK_BLIND_043_RUN_ID:-"043-$(date +%Y%m%d_%H%M%S)-p9016_raw_pairs_remaining_ablation_1m"}
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

APPROVED_P9016_PAIRS=/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz
PAIR_PATH=${HK_BLIND_P9016_PAIRS:-$APPROVED_P9016_PAIRS}
EVAL_PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-$APPROVED_P9016_PAIRS}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_RELAX_STEP:-${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-0.012}}
BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-gpu}
MAX_PARALLEL=${HK_BLIND_043_MAX_PARALLEL:-${HK_BLIND_MAX_PARALLEL:-10}}
RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-1}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-1}
CONFIG_SET=${HK_BLIND_043_CONFIG_SET:-full}
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
		echo "error: HK_BLIND_043_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_043_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi
if [ "$BIN_SIZE" != "1000000" ]; then
	echo "error: 043 is fixed at final 1Mb eval; HK_BLIND_P9016_BIN_SIZE_BP=$BIN_SIZE" >&2
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
	if ! same_realpath "$PAIR_PATH" "$APPROVED_P9016_PAIRS"; then
		echo "error: formal 043 requires approved P9016 training pairs: $APPROVED_P9016_PAIRS" >&2
		echo "got: $PAIR_PATH" >&2
		exit 1
	fi
	if ! same_realpath "$EVAL_PAIRS" "$APPROVED_P9016_PAIRS"; then
		echo "error: formal 043 requires approved P9016 eval pairs: $APPROVED_P9016_PAIRS" >&2
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

find_final_output_dir() {
	local config=$1
	local root="$OUTPUT_ROOT/$config"
	if [ -s "$root/p9016_full.manifest.tsv" ]; then
		printf '%s\n' "$root"
		return 0
	fi
	if [ -s "$root/1m/$config/p9016_full.manifest.tsv" ]; then
		printf '%s\n' "$root/1m/$config"
		return 0
	fi
	local found
	found=$(find "$root" -type f -name p9016_full.manifest.tsv 2>/dev/null | sort | tail -n 1 || true)
	if [ -n "$found" ]; then
		dirname "$found"
		return 0
	fi
	printf '%s\n' "$root"
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

log_cmd make smoke_blind_p9016_minimal
if [ "$MAKE_GPU" = "1" ]; then
	log_cmd make gpu=1 run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
else
	log_cmd make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
fi
log_cmd python -m py_compile scripts/summarize_p9016_model_sweep.py scripts/compute_p9016_copytrack_diag.py

cp "$0" "$FULL_RUN_ROOT/scripts_snapshot/"
cp scripts/p9016_common_eval.sh \
   scripts/p9016_publish_light_result.sh \
   scripts/compute_p9016_copytrack_diag.py \
   scripts/summarize_p9016_model_sweep.py \
   "$FULL_RUN_ROOT/scripts_snapshot/"

RUN_BIN="$FULL_RUN_ROOT/logs/run_blind_p9016_minimal.bin"
AUDIT_BIN="$FULL_RUN_ROOT/logs/audit_blind_p9016_full_cpu_output.bin"
cp run_blind_p9016_minimal.bin "$RUN_BIN"
cp audit_blind_p9016_full_cpu_output.bin "$AUDIT_BIN"
chmod +x "$RUN_BIN" "$AUDIT_BIN"
RUN_HASH=$(hash_file "$RUN_BIN")
AUDIT_HASH=$(hash_file "$AUDIT_BIN")
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
	echo "training_boundary	P9016 raw pairs only; phase labels and CHARM/3DG are eval-only"
	echo "experiment_control	remaining raw-pairs-only ablation after readID/molecule signal was found absent"
	echo "success_criterion	full-denominator trans top1 accuracy +0.1 over td1 baseline"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

run_config() {
	local config=$1
	local chain=$2
	local prior_lambda=$3
	local prior_warmup=$4
	local gate_mode=$5
	local min_pmax=$6
	local min_margin=$7
	local trans_dscale=$8
	local lambda_copytrack=$9
	local lambda_global=${10}
	local out_root="$OUTPUT_ROOT/$config"
	local eval_dir="$EVAL_ROOT/$config"
	local log_prefix="$LOG_ROOT/$config"
	local train_cmd audit_cmd eval_cmd copytrack_cmd final_out

	echo "CONFIG_START=$config chain=$chain prior_lambda=$prior_lambda prior_warmup=$prior_warmup gate_mode=$gate_mode trans_dscale=$trans_dscale lambda_copytrack=$lambda_copytrack lambda_global=$lambda_global"
	rm -rf "$out_root" "$eval_dir"
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
		HK_BLIND_P9016_OUTPUT_ROOT="$out_root"
		HK_BLIND_P9016_CONFIG_NAME="$config"
		HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE"
		HK_BLIND_P9016_MINIMAL_N_ITER="$N_ITER"
		HK_BLIND_P9016_MINIMAL_RELAX_STEPS="$RELAX_STEPS"
		HK_BLIND_P9016_RELAX_STEP="$RELAX_STEP"
		HK_BLIND_P9016_RELAX_BACKEND="$BACKEND"
		HK_BLIND_P9016_CHAIN="$chain"
		HK_BLIND_P9016_RESOLUTION_CHAIN=4000000,1000000
		HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split
		HK_BLIND_P9016_INIT_SEED=17
		HK_BLIND_P9016_INIT_EPS=0.5
		HK_BLIND_P9016_INIT_NOISE_SCALE=0
		HK_BLIND_P9016_D_SCALE_MODE=posterior_count
		HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA=1
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6
		HK_BLIND_P9016_RHO_TRAIN_MODE=constant
		HK_BLIND_P9016_RHO_TRAIN_FLOOR=0
		HK_BLIND_P9016_TEMPERATURE_START=1
		HK_BLIND_P9016_TEMPERATURE_END=1
		HK_BLIND_P9016_TRANS_DSCALE_MULTIPLIER="$trans_dscale"
		HK_BLIND_P9016_TRANS_K_MULTIPLIER=1
		HK_BLIND_P9016_LAMBDA_COPYTRACK="$lambda_copytrack"
		HK_BLIND_P9016_LAMBDA_GLOBAL_COPYTRACK="$lambda_global"
		HK_BLIND_P9016_LAMBDA_NORMDIR_COPYTRACK=0
		HK_BLIND_P9016_MIN_SEP_UNIT=0
		HK_BLIND_P9016_LAMBDA_SEP=0
		HK_BLIND_P9016_TRANS_TOP1_MODE=off
		HK_BLIND_P9016_TRANS_GATE_MODE="$gate_mode"
		HK_BLIND_P9016_TRANS_GATE_MIN_PMAX="$min_pmax"
		HK_BLIND_P9016_TRANS_GATE_MIN_MARGIN="$min_margin"
		HK_BLIND_P9016_TRANS_GATE_MIN_NEG_ENTROPY=-1.38629436
		HK_BLIND_P9016_TRANS_CHR_PAIR_PRIOR_LAMBDA="$prior_lambda"
		HK_BLIND_P9016_TRANS_CHR_PAIR_PRIOR_EPS=0.001
		HK_BLIND_P9016_TRANS_CHR_PAIR_PRIOR_POWER=1
		HK_BLIND_P9016_TRANS_CHR_PAIR_PRIOR_WARMUP_ITER="$prior_warmup"
		HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_LAMBDA=0
		HK_BLIND_WRITE_RAW=0
		"$RUN_BIN"
	)
	log_command_line "${train_cmd[@]}"
	"${train_cmd[@]}" > "$log_prefix.train.log" 2>&1

	final_out=$(find_final_output_dir "$config")
	audit_cmd=(
		clean_env_prefix
		HK_BLIND_AUDIT_ALLOW_CUSTOM_RAW_PAIRS="$ALLOW_CUSTOM"
		"$AUDIT_BIN" "$final_out"
	)
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
			scripts/p9016_common_eval.sh "$final_out" "$eval_dir"
		)
		log_command_line "${eval_cmd[@]}"
		"${eval_cmd[@]}" > "$log_prefix.eval.log" 2>&1
	else
		printf 'status\tSKIPPED_EVAL_DISABLED\n' > "$eval_dir/eval_status.tsv"
	fi

	if [ -s "$final_out/p9016_full.coords.tsv" ]; then
		copytrack_cmd=(
			python scripts/compute_p9016_copytrack_diag.py
			--config-output-dir "$final_out"
			--eval-output-dir "$eval_dir"
			--out "$eval_dir/copytrack_vector_diag.tsv"
		)
		log_command_line "${copytrack_cmd[@]}"
		"${copytrack_cmd[@]}" > "$log_prefix.copytrack.log" 2>&1
	fi
	echo "CONFIG_DONE=$config final_out=$final_out"
}

if [ "$CONFIG_SET" = "smoke" ]; then
	CONFIGS=(
		"p9016_rawremain_smoke_td1_baseline|0|0|0|off|1|1|1|0|0"
		"p9016_rawremain_smoke_best035|0|0|0|off|1|1|0.5|0.03|0.003"
	)
else
	CONFIGS=(
		"p9016_rawremain_td1_baseline|0|0|0|off|1|1|1|0|0"
		"p9016_rawremain_best035_td0p5_lct0p03_gct0p003|0|0|0|off|1|1|0.5|0.03|0.003"
		"p9016_rawremain_chain4m1m_best035|1|0|0|off|1|1|0.5|0.03|0.003"
		"p9016_rawremain_chrpair_estep_l0p10_w20_best035|0|0.10|20|off|1|1|0.5|0.03|0.003"
		"p9016_rawremain_pmaxmargin_p045_m002_best035|0|0|0|pmax_margin|0.45|0.02|0.5|0.03|0.003"
		"p9016_rawremain_chain4m1m_chrpair_l0p10_w20_best035|1|0.10|20|off|1|1|0.5|0.03|0.003"
	)
fi

failures=0
running=0
for spec in "${CONFIGS[@]}"; do
	IFS='|' read -r config chain prior_lambda prior_warmup gate_mode min_pmax min_margin trans_dscale lambda_copytrack lambda_global <<< "$spec"
	run_config "$config" "$chain" "$prior_lambda" "$prior_warmup" "$gate_mode" "$min_pmax" "$min_margin" "$trans_dscale" "$lambda_copytrack" "$lambda_global" &
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

log_command_line scripts/summarize_p9016_model_sweep.py "$FULL_RUN_ROOT" ">" "$FULL_RUN_ROOT/summary.tsv"
scripts/summarize_p9016_model_sweep.py "$FULL_RUN_ROOT" > "$FULL_RUN_ROOT/summary.tsv"

log_command_line python - "$FULL_RUN_ROOT" "[inline raw_pairs_remaining_delta_summary]"
python - "$FULL_RUN_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))
eval_rows = [
    row for row in rows
    if row.get("model_top1_accuracy_genome_trans", "NA") not in ("", "NA")
]

def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None

base = next((r for r in eval_rows if r.get("config_name") == "p9016_rawremain_td1_baseline"), None)
best = next((r for r in eval_rows if r.get("config_name") == "p9016_rawremain_best035_td0p5_lct0p03_gct0p003"), None)
base_t = fnum(base.get("model_top1_accuracy_genome_trans")) if base else None
best_t = fnum(best.get("model_top1_accuracy_genome_trans")) if best else None
for row in eval_rows:
    t = fnum(row.get("model_top1_accuracy_genome_trans"))
    if t is not None and base_t is not None:
        row["delta_trans_vs_td1_baseline"] = f"{t - base_t:.9g}"
        row["target_plus_0p1_vs_td1_met"] = "1" if t >= base_t + 0.1 else "0"
    else:
        row["delta_trans_vs_td1_baseline"] = "NA"
        row["target_plus_0p1_vs_td1_met"] = "0"
    if t is not None and best_t is not None:
        row["delta_trans_vs_best035"] = f"{t - best_t:.9g}"
    else:
        row["delta_trans_vs_best035"] = "NA"

fields = [
    "config_name",
    "status",
    "bin_size_bp",
    "backend",
    "refinement_stage",
    "refinement_stage_index",
    "refinement_n_stages",
    "trans_chr_pair_prior_lambda",
    "trans_chr_pair_prior_warmup_iter",
    "trans_gate_mode",
    "trans_gate_min_pmax",
    "trans_gate_min_margin",
    "trans_dscale_multiplier",
    "lambda_copytrack",
    "lambda_global_copytrack",
    "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans",
    "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "final_mean_entropy",
    "final_mean_pU",
    "final_min_sep",
    "sep_p05",
    "final_mean_sep",
    "copytrack_frac_cos_lt_0",
    "copytrack_frac_projection_sign_switch",
    "delta_trans_vs_td1_baseline",
    "delta_trans_vs_best035",
    "target_plus_0p1_vs_td1_met",
]
rows_sorted = sorted(eval_rows, key=lambda r: fnum(r.get("model_top1_accuracy_genome_trans")) or -1, reverse=True)
with (root / "raw_pairs_remaining_delta_summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows_sorted:
        writer.writerow({key: row.get(key, "NA") for key in fields})

headline = "PLUS_0P1_MET" if any(row.get("target_plus_0p1_vs_td1_met") == "1" for row in eval_rows) else "NO_PLUS_0P1"
(root / "headline.txt").write_text(headline + "\n")
PY

log_command_line python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" "[inline README]"
python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
light = Path(sys.argv[2])
rows = list(csv.DictReader((root / "raw_pairs_remaining_delta_summary.tsv").open(), delimiter="\t"))
headline = (root / "headline.txt").read_text().strip()
fields = [
    "config_name",
    "trans_chr_pair_prior_lambda",
    "trans_chr_pair_prior_warmup_iter",
    "trans_gate_mode",
    "trans_gate_min_pmax",
    "trans_gate_min_margin",
    "trans_dscale_multiplier",
    "lambda_copytrack",
    "lambda_global_copytrack",
    "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans",
    "delta_trans_vs_td1_baseline",
    "delta_trans_vs_best035",
    "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "final_mean_entropy",
    "final_mean_pU",
    "final_min_sep",
    "sep_p05",
    "copytrack_frac_cos_lt_0",
    "copytrack_frac_projection_sign_switch",
]
with (root / "README.md").open("w") as out:
    out.write("# 043 P9016 Raw-Pairs Remaining Ablation 1Mb\n\n")
    out.write("This controlled blind-training experiment tests the last small set of raw-pairs-only mechanisms that remained after prior P9016 trans experiments: coarse-to-fine 4Mb-to-1Mb initialization path, weak trans chromosome-pair E-step prior, and a mild pmax+margin trans M-step gate.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- headline: `{headline}`\n")
    out.write("- training input: approved P9016 raw pairs only\n")
    out.write("- eval-only inputs: SNP phase labels and CHARM/3DG reference are used only after training through the standard eval wrapper\n")
    out.write("- copy0/copy1 are gauge labels; this run does not use maternal/paternal labels in training\n")
    out.write("- success criterion: full-denominator trans top1 accuracy must improve by at least 0.1 over the td1 baseline\n\n")
    out.write("## Why These Configs\n\n")
    out.write("- `chain4m1m`: tests whether a coarse-to-fine basin changes trans copy-gauge identity when local copytrack alone does not.\n")
    out.write("- `chrpair_estep_l0p10_w20`: tests a weak posterior-derived chromosome-pair E-step prior without changing M-step graph topology.\n")
    out.write("- `pmaxmargin_p045_m002`: fills the mild pmax+margin gate hole left after stronger confidence gates failed.\n\n")
    out.write("## Main Results\n\n")
    out.write("| " + " | ".join(fields) + " |\n")
    out.write("|" + "|".join(["---"] * len(fields)) + "|\n")
    for row in rows:
        out.write("| " + " | ".join(row.get(key, "NA") for key in fields) + " |\n")
    out.write("\n## Interpretation Boundary\n\n")
    out.write("This is not a new baseline unless the full-denominator trans top1 criterion is met without degrading cis accuracy or cis distance Spearman. If these configs do not exceed the best035 family, the remaining bottleneck is unlikely to be solvable by reusing the current flattened P9016 pairs as only pairwise contact counts.\n")
PY

if [ "$CONFIG_SET" = "smoke" ]; then
	echo "SKIPPED_LIGHT_PUBLISH_FOR_SMOKE" > "$FULL_RUN_ROOT/publish.log"
	if [ "${HK_BLIND_KEEP_SMOKE_RESULT:-0}" != "1" ]; then
		echo "SMOKE_RESULT_REMOVED=$FULL_RUN_ROOT"
		rm -rf "$FULL_RUN_ROOT"
	fi
else
	log_command_line scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" ">" "$FULL_RUN_ROOT/publish.log"
	scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$FULL_RUN_ROOT/publish.log"
fi

if [ -d "$FULL_RUN_ROOT" ]; then
	{
		echo "end_time	$(date -Is)"
		echo "headline	$(cat "$FULL_RUN_ROOT/headline.txt" 2>/dev/null || echo SMOKE_RESULT_REMOVED)"
	} >> "$FULL_RUN_ROOT/run_manifest.tsv"
fi

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
if [ -s "$FULL_RUN_ROOT/summary.tsv" ]; then
	echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
fi
if [ -s "$FULL_RUN_ROOT/raw_pairs_remaining_delta_summary.tsv" ]; then
	echo "RAW_PAIRS_REMAINING_DELTA_SUMMARY_TSV=$FULL_RUN_ROOT/raw_pairs_remaining_delta_summary.tsv"
fi
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt" 2>/dev/null || echo SMOKE_RESULT_REMOVED)"
