#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
HICKIT_DIR="$REPO_ROOT/hickit"
APPROVED_PAIRS="/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz"
DEFAULT_EVAL_TDG="/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz"
RUN_ID="011-$(date +%Y%m%d_%H%M%S)-p9016_regression_copytrack_1m"
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$REPO_ROOT/test_res}/$RUN_ID"
FULL_RUN_ROOT=$(realpath -m "$FULL_RUN_ROOT")
OUTPUT_ROOT="$FULL_RUN_ROOT/outputs"
LOG_DIR="$FULL_RUN_ROOT/logs"
SUMMARY_TSV="$FULL_RUN_ROOT/summary.tsv"
REGRESSION_AUDIT_TSV="$FULL_RUN_ROOT/regression_audit.tsv"
COPYTRACK_SUMMARY_TSV="$FULL_RUN_ROOT/copytrack_summary.tsv"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"
LIGHT_RESULT_ROOT="$HICKIT_DIR/result/$RUN_ID"

if [[ "${CONDA_DEFAULT_ENV:-}" != "analysis" ]]; then
	printf 'error: activate the analysis Conda environment before running this script\n' >&2
	exit 2
fi

mkdir -p "$OUTPUT_ROOT" "$LOG_DIR" "$FULL_RUN_ROOT/eval"
: > "$COMMANDS_LOG"

log() {
	printf '%s\n' "$*" | tee -a "$COMMANDS_LOG"
}

run_logged() {
	log "+ $*"
	"$@"
}

run_logged_capture() {
	local log_path="$1"
	shift
	log "+ $* > $log_path 2>&1"
	"$@" >"$log_path" 2>&1
}

run_logged_shell() {
	local description="$1"
	shift
	log "+ $description"
	bash -c "$*"
}

write_status_eval() {
	local config_name="$1"
	local status="$2"
	local eval_dir="$FULL_RUN_ROOT/eval/$config_name"
	mkdir -p "$eval_dir"
	{
		printf 'status\t%s\n' "$status"
		printf 'config_name\t%s\n' "$config_name"
		printf 'created_at\t%s\n' "$(date -Is)"
	} >"$eval_dir/eval_status.tsv"
	printf '%s\n' "$status" >"$eval_dir/eval.log"
}

write_skip_manifest() {
	local config_name="$1"
	local status="$2"
	local config_dir="$OUTPUT_ROOT/$config_name"
	mkdir -p "$config_dir"
	{
		printf 'key\tvalue\n'
		printf 'sample\tP9016\n'
		printf 'runner_family\tp9016_minimal\n'
		printf 'config_name\t%s\n' "$config_name"
		printf 'status\t%s\n' "$status"
		printf 'output_dir\t%s\n' "$config_dir"
		printf 'relax_backend\tgpu\n'
		printf 'bin_size_bp\t%s\n' "$BIN_SIZE"
		printf 'resolution\t%s\n' "$BIN_SIZE"
		printf 'init_mode\tunphased_scaffold_split\n'
		printf 'init_seed\t17\n'
		printf 'init_eps\t0.5\n'
		printf 'init_noise_scale\t0\n'
		printf 'init_scale\t10\n'
		printf 'init_eps_effective\t0.5\n'
		printf 'init_noise_scale_effective\t0\n'
		printf 'init_scale_effective\t0\n'
		printf 'prior_mode\tuniform\n'
		printf 'rho_train_mode\tconstant\n'
		printf 'd_scale_mode\traw_count\n'
		printf 'd_scale_eps_count\t1e-6\n'
		printf 'min_sep_unit\t0\n'
		printf 'lambda_sep\t0\n'
		printf 'temperature_start\t1\n'
		printf 'temperature_end\t1\n'
		printf 'uses_phase_labels\t0\n'
		printf 'uses_charm_or_reference\t0\n'
	} >"$config_dir/p9016_full.manifest.tsv"
}

metric_value() {
	local table="$1"
	local config="$2"
	local key="$3"
	awk -F'\t' -v cfg="$config" -v key="$key" '
		NR == 1 {
			for (i = 1; i <= NF; ++i) {
				idx[$i] = i
			}
			next
		}
		$(idx["config_name"]) == cfg {
			if (key in idx) {
				print $(idx[key])
			} else {
				print "NA"
			}
			found = 1
			exit
		}
		END {
			if (!found) print "NA"
		}
	' "$table"
}

write_copytrack_summary() {
	python - "$SUMMARY_TSV" "$COPYTRACK_SUMMARY_TSV" <<'PY'
import csv
import sys

summary_path, out_path = sys.argv[1:3]
cols = [
    "config_name",
    "backend",
    "copytrack_sep_p05",
    "copytrack_sep_median",
    "copytrack_vector_cos_p05",
    "copytrack_vector_cos_median",
    "copytrack_frac_cos_lt_0",
    "copytrack_frac_cos_lt_neg0p5",
    "copytrack_frac_projection_sign_switch",
    "copytrack_median_projection_run_length",
]
with open(summary_path, newline="") as fh, open(out_path, "w", newline="") as out:
    reader = csv.DictReader(fh, delimiter="\t")
    writer = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in reader:
        writer.writerow({col: row.get(col, "NA") for col in cols})
PY
}

regression_headline() {
	if [[ ! -s "$REGRESSION_AUDIT_TSV" ]]; then
		printf 'NOT_RUN'
		return
	fi
	if awk -F'\t' 'NR > 1 && $9 == "FAIL" { fail=1 } END { exit fail?0:1 }' "$REGRESSION_AUDIT_TSV"; then
		printf 'FAIL'
	else
		printf 'PASS'
	fi
}

write_readme() {
	local regression_status="$1"
	local light_root="$2"
	{
		printf '# %s\n\n' "$RUN_ID"
		printf 'This is a controlled P9016 blind-diploid infrastructure and regression audit. It does not introduce a new model.\n\n'
		printf '## Paths\n\n'
		printf -- '- full result root: `%s`\n' "$FULL_RUN_ROOT"
		printf -- '- lightweight hickit result: `%s`\n' "$light_root"
		printf -- '- summary TSV: `%s`\n' "$SUMMARY_TSV"
		printf -- '- regression audit TSV: `%s`\n' "$REGRESSION_AUDIT_TSV"
		printf -- '- copytrack summary TSV: `%s`\n\n' "$COPYTRACK_SUMMARY_TSV"
		printf '## Training Boundary\n\n'
		printf 'Training reads only the P9016 raw pairs file. Phase labels and CHARM/3DG are passed only to post-training eval. Copy labels are gauge labels, not parental labels.\n\n'
		printf '## Configs\n\n'
		printf '| config | backend | init_eps | init_noise_scale | min_sep_unit | lambda_sep | d_scale_mode |\n'
		printf '| --- | --- | ---: | ---: | ---: | ---: | --- |\n'
		printf '| legacy_1m_cpu_exact_ieps0p5_noise0_seed17 | cpu | 0.5 | 0.0 | 0.0 | 0.0 | raw_count |\n'
		printf '| legacy_1m_gpu_exact_ieps0p5_noise0_seed17 | gpu | 0.5 | 0.0 | 0.0 | 0.0 | raw_count |\n'
		printf '| common_raw_sep_off_ieps1_noise0p05_seed17 | %s | 1.0 | 0.05 | 0.0 | 0.0 | raw_count |\n' "$PREFERRED_BACKEND"
		printf '| common_raw_msep1p5_lsep1_ieps1_noise0p05_seed17 | %s | 1.0 | 0.05 | 1.5 | 1.0 | raw_count |\n\n' "$PREFERRED_BACKEND"
		printf '## Regression Audit\n\n'
		if [[ "$regression_status" == "FAIL" ]]; then
			printf 'REGRESSION_FAIL: at least one exact legacy baseline metric or manifest field differs from the historical 003/004-style baseline beyond tolerance. Debug this before interpreting later rows as model changes.\n\n'
		else
			printf 'REGRESSION_PASS: exact legacy baseline rows match the available historical 003/004-style reference within the configured tolerances, or skipped rows were explicitly marked.\n\n'
		fi
		printf '## Headline Metrics\n\n'
		printf '| config | backend | top1 all | top1 cis | top1 trans | same/cross cis | mean cis Spearman | entropy | pU | min sep | mean sep | frac cos<0 | projection switch frac |\n'
		printf '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n'
		for config in \
			legacy_1m_cpu_exact_ieps0p5_noise0_seed17 \
			legacy_1m_gpu_exact_ieps0p5_noise0_seed17 \
			common_raw_sep_off_ieps1_noise0p05_seed17 \
			common_raw_msep1p5_lsep1_ieps1_noise0p05_seed17
		do
			printf '| `%s` | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n' \
				"$config" \
				"$(metric_value "$SUMMARY_TSV" "$config" backend)" \
				"$(metric_value "$SUMMARY_TSV" "$config" model_top1_accuracy_genome_all)" \
				"$(metric_value "$SUMMARY_TSV" "$config" model_top1_accuracy_genome_cis)" \
				"$(metric_value "$SUMMARY_TSV" "$config" model_top1_accuracy_genome_trans)" \
				"$(metric_value "$SUMMARY_TSV" "$config" model_same_cross_accuracy_genome_cis)" \
				"$(metric_value "$SUMMARY_TSV" "$config" mean_per_chrom_cis_distance_spearman)" \
				"$(metric_value "$SUMMARY_TSV" "$config" final_mean_entropy)" \
				"$(metric_value "$SUMMARY_TSV" "$config" final_mean_pU)" \
				"$(metric_value "$SUMMARY_TSV" "$config" final_min_sep)" \
				"$(metric_value "$SUMMARY_TSV" "$config" final_mean_sep)" \
				"$(metric_value "$SUMMARY_TSV" "$config" copytrack_frac_cos_lt_0)" \
				"$(metric_value "$SUMMARY_TSV" "$config" copytrack_frac_projection_sign_switch)"
		done
		printf '\n## Interpretation\n\n'
		printf 'This run is a diagnostic and regression audit, not a final model. Stronger separation is only a guardrail: it should increase min/p05 separation, but it is not expected to solve 00-vs-11 identity alone. High `copytrack_frac_cos_lt_0` or high projection sign-switch rate supports local homolog-vector flipping as a bottleneck. If four-state top1 accuracy is low while same/cross cis is high, and copy-track discontinuity is high, the next controlled experiment should test a small homolog-vector continuity prior. If copy-track discontinuity is low, the next diagnostic should target oracle posterior or E-step prior behavior.\n'
	} >"$FULL_RUN_ROOT/README.md"
}

append_run_manifest() {
	local key="$1"
	local value="$2"
	printf '%s\t%s\n' "$key" "$value" >> "$FULL_RUN_ROOT/run_manifest.tsv"
}

gpu_build_available() {
	command -v nvidia-smi >/dev/null 2>&1 &&
		command -v /usr/local/cuda/bin/nvcc >/dev/null 2>&1 &&
		nvidia-smi >/dev/null 2>&1
}

resolve_existing_path() {
	local path="$1"
	if [[ "$path" == /* && -e "$path" ]]; then
		realpath -e "$path"
		return
	fi
	if [[ -e "$path" ]]; then
		realpath -e "$path"
		return
	fi
	if [[ -e "$REPO_ROOT/$path" ]]; then
		realpath -e "$REPO_ROOT/$path"
		return
	fi
	if [[ -e "$HICKIT_DIR/$path" ]]; then
		realpath -e "$HICKIT_DIR/$path"
		return
	fi
	printf 'error: path does not exist: %s\n' "$path" >&2
	return 1
}

run_config() {
	local config_name="$1"
	local backend="$2"
	local init_eps="$3"
	local init_noise="$4"
	local min_sep="$5"
	local lambda_sep="$6"
	local config_dir="$OUTPUT_ROOT/$config_name"
	local eval_dir="$FULL_RUN_ROOT/eval/$config_name"
	local log_prefix="$LOG_DIR/$config_name"

	log "CONFIG_START=$config_name"
	mkdir -p "$eval_dir"
	run_logged_capture "${log_prefix}.train.log" \
		env \
		HK_BLIND_P9016_PAIRS="$PAIRS_REAL" \
		HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS="${HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS:-0}" \
		HK_BLIND_P9016_OUTPUT_ROOT="$OUTPUT_ROOT" \
		HK_BLIND_P9016_CONFIG_NAME="$config_name" \
		HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE" \
		HK_BLIND_P9016_MINIMAL_N_ITER="$N_ITER" \
		HK_BLIND_P9016_MINIMAL_RELAX_STEPS="$RELAX_STEPS" \
		HK_BLIND_P9016_RELAX_STEP="$RELAX_STEP" \
		HK_BLIND_P9016_RELAX_BACKEND="$backend" \
		HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split \
		HK_BLIND_P9016_INIT_SCALE=10 \
		HK_BLIND_P9016_INIT_SEED=17 \
		HK_BLIND_P9016_INIT_EPS="$init_eps" \
		HK_BLIND_P9016_INIT_NOISE_SCALE="$init_noise" \
		HK_BLIND_P9016_MIN_SEP_UNIT="$min_sep" \
		HK_BLIND_P9016_LAMBDA_SEP="$lambda_sep" \
		HK_BLIND_P9016_D_SCALE_MODE=raw_count \
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6 \
		"$HICKIT_DIR/run_blind_p9016_minimal.bin"

	if [[ -s "$config_dir/p9016_full.manifest.tsv" ]]; then
		{
			printf 'git_commit\t%s\n' "$GIT_COMMIT"
			printf 'binary_hash\t%s\n' "$BINARY_HASH"
		} >> "$config_dir/p9016_full.manifest.tsv"
	fi

	run_logged_capture "${log_prefix}.audit.log" \
		"$HICKIT_DIR/audit_blind_p9016_full_cpu_output.bin" "$config_dir"

	if [[ "$RUN_EVAL" == "1" ]]; then
		run_logged_capture "${log_prefix}.eval_wrapper.log" \
			env \
			HK_BLIND_P9016_EVAL_PAIRS="$EVAL_PAIRS" \
			HK_BLIND_P9016_EVAL_TDG="$EVAL_TDG" \
			HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE" \
			"$SCRIPT_DIR/p9016_common_eval.sh" "$config_dir" "$eval_dir"
	else
		write_status_eval "$config_name" "SKIPPED_EVAL_CUSTOM_PAIRS"
	fi

	run_logged_capture "${log_prefix}.copytrack.log" \
		python "$SCRIPT_DIR/compute_p9016_copytrack_diag.py" \
		--config-output-dir "$config_dir" \
		--eval-output-dir "$eval_dir" \
		--out "$eval_dir/copytrack_vector_diag.tsv"
	log "CONFIG_DONE=$config_name"
}

BIN_SIZE="${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}"
N_ITER="${HK_BLIND_P9016_MINIMAL_N_ITER:-100}"
RELAX_STEPS="${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}"
RELAX_STEP="${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-${HK_BLIND_P9016_RELAX_STEP:-0.012}}"
PAIRS="${HK_BLIND_P9016_PAIRS:-$APPROVED_PAIRS}"
EVAL_PAIRS="${HK_BLIND_P9016_EVAL_PAIRS:-$APPROVED_PAIRS}"
EVAL_TDG="${HK_BLIND_P9016_EVAL_TDG:-$DEFAULT_EVAL_TDG}"

if [[ "$BIN_SIZE" != "1000000" ]]; then
	printf 'error: experiment 011 is fixed at 1 Mb; got HK_BLIND_P9016_BIN_SIZE_BP=%s\n' "$BIN_SIZE" >&2
	exit 2
fi

APPROVED_PAIRS_REAL=$(realpath -e "$APPROVED_PAIRS")
PAIRS_REAL=$(resolve_existing_path "$PAIRS")
if [[ "$PAIRS_REAL" != "$APPROVED_PAIRS_REAL" && "${HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS:-0}" != "1" ]]; then
	printf 'error: training pairs must be approved P9016 raw pairs realpath: %s\n' "$APPROVED_PAIRS_REAL" >&2
	printf '       got: %s\n' "$PAIRS_REAL" >&2
	printf '       set HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS=1 only for explicit non-baseline diagnostics\n' >&2
	exit 2
fi

RUN_EVAL="${HK_BLIND_P9016_RUN_EVAL:-1}"
if [[ "$PAIRS_REAL" != "$APPROVED_PAIRS_REAL" && "${HK_BLIND_P9016_RUN_EVAL+x}" != "x" ]]; then
	RUN_EVAL=0
fi
PREFERRED_BACKEND="${HK_BLIND_P9016_RELAX_BACKEND:-auto}"
case "$PREFERRED_BACKEND" in
	cpu|gpu|auto) ;;
	*)
		printf 'error: HK_BLIND_P9016_RELAX_BACKEND must be cpu, gpu, or auto; got %s\n' "$PREFERRED_BACKEND" >&2
		exit 2
		;;
esac

GPU_AVAILABLE=0
MAKE_GPU=0
if gpu_build_available; then
	GPU_AVAILABLE=1
	MAKE_GPU=1
	if [[ "$PREFERRED_BACKEND" == "auto" ]]; then
		PREFERRED_BACKEND=gpu
	fi
elif [[ "$PREFERRED_BACKEND" == "auto" ]]; then
	PREFERRED_BACKEND=cpu
fi

printf 'FULL_RESULT_ROOT=%s\n' "$FULL_RUN_ROOT"
log "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
log "OUTPUT_ROOT=$OUTPUT_ROOT"
log "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

run_logged cp "$0" "$LOG_DIR/run_p9016_011_regression_copytrack.sh"
mkdir -p "$FULL_RUN_ROOT/scripts_snapshot"
run_logged cp "$SCRIPT_DIR/p9016_common_eval.sh" "$SCRIPT_DIR/compute_p9016_copytrack_diag.py" \
	"$SCRIPT_DIR/summarize_p9016_model_sweep.py" "$SCRIPT_DIR/audit_p9016_baseline_regression.py" \
	"$SCRIPT_DIR/p9016_publish_light_result.sh" "$FULL_RUN_ROOT/scripts_snapshot/"

cd "$HICKIT_DIR"
run_logged make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
run_logged make smoke_blind_p9016_minimal
if (( MAKE_GPU == 1 )); then
	run_logged make gpu=1 run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
fi

GIT_COMMIT=$(git -C "$HICKIT_DIR" rev-parse HEAD)
BINARY_HASH=$(sha256sum "$HICKIT_DIR/run_blind_p9016_minimal.bin" | awk '{print $1}')
{
	printf 'key\tvalue\n'
} > "$FULL_RUN_ROOT/run_manifest.tsv"
append_run_manifest run_id "$RUN_ID"
append_run_manifest full_run_root "$FULL_RUN_ROOT"
append_run_manifest light_result_root "$LIGHT_RESULT_ROOT"
append_run_manifest created_at "$(date -Is)"
append_run_manifest hickit_dir "$HICKIT_DIR"
append_run_manifest hickit_git_commit "$GIT_COMMIT"
append_run_manifest hickit_git_dirty "$(git -C "$HICKIT_DIR" status --short | wc -l | awk '{print $1}')"
append_run_manifest binary_hash "$BINARY_HASH"
append_run_manifest pairs "$PAIRS_REAL"
append_run_manifest eval_pairs "$EVAL_PAIRS"
append_run_manifest eval_tdg "$EVAL_TDG"
append_run_manifest bin_size_bp "$BIN_SIZE"
append_run_manifest n_iter "$N_ITER"
append_run_manifest relax_steps "$RELAX_STEPS"
append_run_manifest relax_step "$RELAX_STEP"
append_run_manifest preferred_backend "$PREFERRED_BACKEND"
append_run_manifest gpu_available "$GPU_AVAILABLE"
append_run_manifest run_eval "$RUN_EVAL"

run_config legacy_1m_cpu_exact_ieps0p5_noise0_seed17 cpu 0.5 0.0 0.0 0.0
if (( GPU_AVAILABLE == 1 )); then
	run_config legacy_1m_gpu_exact_ieps0p5_noise0_seed17 gpu 0.5 0.0 0.0 0.0
else
	write_skip_manifest legacy_1m_gpu_exact_ieps0p5_noise0_seed17 SKIPPED_GPU_UNAVAILABLE
	write_status_eval legacy_1m_gpu_exact_ieps0p5_noise0_seed17 SKIPPED_GPU_UNAVAILABLE
fi
run_config common_raw_sep_off_ieps1_noise0p05_seed17 "$PREFERRED_BACKEND" 1.0 0.05 0.0 0.0
run_config common_raw_msep1p5_lsep1_ieps1_noise0p05_seed17 "$PREFERRED_BACKEND" 1.0 0.05 1.5 1.0

run_logged_shell "scripts/summarize_p9016_model_sweep.py \"$FULL_RUN_ROOT\" > \"$SUMMARY_TSV\" 2> \"$LOG_DIR/summarize.log\"" \
	"\"$SCRIPT_DIR/summarize_p9016_model_sweep.py\" \"$FULL_RUN_ROOT\" > \"$SUMMARY_TSV\" 2> \"$LOG_DIR/summarize.log\""

run_logged_shell "scripts/audit_p9016_baseline_regression.py \"$FULL_RUN_ROOT\" > \"$REGRESSION_AUDIT_TSV\" 2> \"$LOG_DIR/regression_audit.log\"" \
	"\"$SCRIPT_DIR/audit_p9016_baseline_regression.py\" \"$FULL_RUN_ROOT\" > \"$REGRESSION_AUDIT_TSV\" 2> \"$LOG_DIR/regression_audit.log\""

write_copytrack_summary
REGRESSION_STATUS=$(regression_headline)
write_readme "$REGRESSION_STATUS" "$LIGHT_RESULT_ROOT"
run_logged "$SCRIPT_DIR/p9016_publish_light_result.sh" "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"

printf 'FULL_RESULT_ROOT=%s\n' "$FULL_RUN_ROOT"
printf 'LIGHT_RESULT_ROOT=%s\n' "$LIGHT_RESULT_ROOT"
printf 'SUMMARY_TSV=%s\n' "$SUMMARY_TSV"
printf 'REGRESSION_AUDIT_TSV=%s\n' "$REGRESSION_AUDIT_TSV"
printf 'COPYTRACK_SUMMARY_TSV=%s\n' "$COPYTRACK_SUMMARY_TSV"
