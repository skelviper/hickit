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

RUN_ID="016-$(date +%Y%m%d_%H%M%S)-p9016_charm_contacts_condition_sweep_1m"
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE_ROOT/test_res}/$RUN_ID"
mkdir -p "$FULL_RUN_ROOT"/{logs,outputs,eval,work_outputs,synthetic_pairs,scripts_snapshot,contact_matrix_qc}
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"
COMMANDS_LOCK_DIR="$FULL_RUN_ROOT/.commands.lock"

exec > >(tee -a "$FULL_RUN_ROOT/logs/runner.stdout.log") 2> >(tee -a "$FULL_RUN_ROOT/logs/runner.stderr.log" >&2)

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

GIT_COMMIT=$(git rev-parse HEAD)
GIT_DIRTY_COUNT=$(git status --porcelain | wc -l | awk '{print $1}')
git status --porcelain=v1 > "$FULL_RUN_ROOT/logs/git_status.txt"
git diff --stat > "$FULL_RUN_ROOT/logs/git_diff_stat.txt"
git diff > "$FULL_RUN_ROOT/logs/git_diff.patch"
{
	echo "cc	$(${CC:-cc} --version | head -1)"
	echo "make	$(make --version | head -1)"
	echo "cflags	${CFLAGS:-default}"
	echo "conda_env	${CONDA_DEFAULT_ENV:-NA}"
	echo "path	$PATH"
} > "$FULL_RUN_ROOT/logs/build_env.txt"

if [ "${HK_BLIND_REQUIRE_CLEAN_TREE:-0}" = "1" ] && [ "$GIT_DIRTY_COUNT" != "0" ]; then
	echo "error: hickit tree is dirty and HK_BLIND_REQUIRE_CLEAN_TREE=1" >&2
	exit 1
fi

cp "$0" "$FULL_RUN_ROOT/scripts_snapshot/"
cp scripts/generate_charm3dg_synthetic_pairs.py \
   scripts/infer_charm3dg_pair_copy_state.py \
   scripts/plot_p9016_contact_matrix_qc.py \
   scripts/summarize_p9016_016.py \
   scripts/p9016_common_eval.sh \
   scripts/p9016_publish_light_result.sh \
   scripts/compute_p9016_copytrack_diag.py \
   "$FULL_RUN_ROOT/scripts_snapshot/"

with_commands_lock() {
	while ! mkdir "$COMMANDS_LOCK_DIR" 2>/dev/null; do
		sleep 0.05
	done
	"$@" >> "$COMMANDS_LOG"
	rmdir "$COMMANDS_LOCK_DIR"
}

log_line() {
	with_commands_lock printf '%s\n' "$*"
}

log_command_line() {
	local line="+"
	local arg
	for arg in "$@"; do
		line+=" $(printf '%q' "$arg")"
	done
	log_line "$line"
}

log_cmd() {
	log_command_line "$@"
	"$@"
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
	if [ "${CUDA_VISIBLE_DEVICES+x}" = "x" ]; then
		env_cmd+=(CUDA_VISIBLE_DEVICES="$CUDA_VISIBLE_DEVICES")
	fi
	"${env_cmd[@]}" "$@"
}

hash_file() {
	if command -v sha256sum >/dev/null 2>&1; then
		sha256sum "$1" | awk '{print $1}'
	else
		echo NA
	fi
}

backend_available() {
	if [ "$1" = "gpu" ]; then
		command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1
	else
		return 0
	fi
}

radius_slug() {
	case "$1" in
	1|1.0) echo "1pr" ;;
	2|2.0) echo "2pr" ;;
	*) echo "${1//./p}pr" ;;
	esac
}

gamma_slug() {
	case "$1" in
	0|0.0) echo "0" ;;
	0.25) echo "0p25" ;;
	0.5) echo "0p5" ;;
	0.75) echo "0p75" ;;
	1|1.0) echo "1" ;;
	*) echo "$1" | tr '.' 'p' ;;
	esac
}

sep_slug() {
	local min_sep=$1
	local lambda_sep=$2
	if [ "$min_sep" = "0" ] || [ "$min_sep" = "0.0" ]; then
		echo "sep_off"
	else
		echo "msep${min_sep//./p}_lsep${lambda_sep//./p}"
	fi
}

PAIRS=${HK_BLIND_P9016_OBSERVED_PAIRS:-${HK_BLIND_P9016_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}}
EVAL_PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-$PAIRS}
TDG20=${HK_BLIND_P9016_20K_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.20k.3dg.gz}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
TDG_RES=${HK_BLIND_016_TDG_RESOLUTION_BP:-20000}
RADII=${HK_BLIND_016_RADII_PR:-"1 2"}
GAMMAS=${HK_BLIND_016_GAMMAS:-"0 0.25 0.5 0.75 1"}
INCLUDE_RANDOM=${HK_BLIND_016_INCLUDE_RANDOM:-1}
RUN_SEP=${HK_BLIND_016_RUN_SEP:-1}
RUN_COMMON_INIT=${HK_BLIND_016_RUN_COMMON_INIT:-1}
TARGET_MODE=${HK_BLIND_016_TARGET_MODE:-observed_non_same_bin}
TARGET_COUNT=${HK_BLIND_016_TARGET_COUNT:-0}
CANDIDATE_MULTIPLIER=${HK_BLIND_016_CANDIDATE_MULTIPLIER:-4}
PARTICLE_RADIUS=${HK_BLIND_016_PARTICLE_RADIUS:-0}
MIN_GENOMIC_SEP=${HK_BLIND_016_MIN_GENOMIC_SEPARATION:-0}
SEED=${HK_BLIND_016_SEED:-17}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-0.012}
MAX_PARALLEL=${HK_BLIND_016_MAX_PARALLEL:-10}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-1}
NEAREST_IF_INSUFFICIENT=${HK_BLIND_016_NEAREST_IF_INSUFFICIENT:-0}
REQUESTED_BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-gpu}
if [ "$REQUESTED_BACKEND" = "auto" ]; then
	if backend_available gpu; then
		REQUESTED_BACKEND=gpu
	else
		REQUESTED_BACKEND=cpu
	fi
fi
if [ "$REQUESTED_BACKEND" = "gpu" ] && ! backend_available gpu; then
	echo "error: GPU backend requested but unavailable" >&2
	exit 1
fi

case "$MAX_PARALLEL" in
	''|*[!0-9]*)
		echo "error: HK_BLIND_016_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_016_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi

for required in "$PAIRS" "$EVAL_PAIRS" "$TDG20" "$EVAL_TDG"; do
	if [ ! -s "$required" ]; then
		echo "error: required input missing: $required" >&2
		exit 1
	fi
done

log_cmd make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
log_cmd make smoke_blind_p9016_minimal
log_cmd scripts/smoke_p9016_dscale_gamma.sh
if [ "$REQUESTED_BACKEND" = "gpu" ]; then
	log_cmd make gpu=1 run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
else
	log_cmd make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
fi

RUN_BIN="$FULL_RUN_ROOT/logs/run_blind_p9016_minimal.$REQUESTED_BACKEND.bin"
cp run_blind_p9016_minimal.bin "$RUN_BIN"
chmod +x "$RUN_BIN"
RUN_HASH=$(hash_file "$RUN_BIN")
AUDIT_HASH=$(hash_file audit_blind_p9016_full_cpu_output.bin)

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "run_blind_p9016_minimal_sha256	$RUN_HASH"
	echo "audit_blind_p9016_full_cpu_output_sha256	$AUDIT_HASH"
	echo "sample	P9016"
	echo "observed_pairs	$PAIRS"
	echo "training_source_20k_tdg	$TDG20"
	echo "eval_pairs	$EVAL_PAIRS"
	echo "eval_tdg	$EVAL_TDG"
	echo "radii_pr	$RADII"
	echo "gammas	$GAMMAS"
	echo "run_sep	$RUN_SEP"
	echo "run_common_init	$RUN_COMMON_INIT"
	echo "include_random	$INCLUDE_RANDOM"
	echo "target_mode	$TARGET_MODE"
	echo "backend	$REQUESTED_BACKEND"
	echo "bin_size_bp	$BIN_SIZE"
	echo "n_iter	$N_ITER"
	echo "relax_steps	$RELAX_STEPS"
	echo "relax_step	$RELAX_STEP"
	echo "max_parallel	$MAX_PARALLEL"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"
{
	echo "run_blind_p9016_minimal_sha256	$RUN_HASH"
	echo "audit_blind_p9016_full_cpu_output_sha256	$AUDIT_HASH"
} >> "$FULL_RUN_ROOT/logs/build_env.txt"

generate_pairs_for_radius() {
	local radius=$1
	local slug=$2
	local out_pairs="$FULL_RUN_ROOT/synthetic_pairs/p9016_charm3dg20k_${slug}_contacts.pairs.gz"
	local copy_state="$FULL_RUN_ROOT/synthetic_pairs/p9016_charm3dg20k_${slug}_contacts.copy_state.tsv.gz"
	local metadata="$FULL_RUN_ROOT/synthetic_pairs/p9016_charm3dg20k_${slug}_contacts.metadata.tsv"
	local metadata_json="$FULL_RUN_ROOT/synthetic_pairs/p9016_charm3dg20k_${slug}_contacts.json"
	local cmd=(
		python scripts/generate_charm3dg_synthetic_pairs.py
		--tdg "$TDG20"
		--observed-pairs "$PAIRS"
		--out-pairs "$out_pairs"
		--copy-state-out "$copy_state"
		--metadata "$metadata"
		--json "$metadata_json"
		--training-bin-size "$BIN_SIZE"
		--tdg-resolution "$TDG_RES"
		--radius-pr "$radius"
		--particle-radius "$PARTICLE_RADIUS"
		--target-mode "$TARGET_MODE"
		--target-count "$TARGET_COUNT"
		--candidate-multiplier "$CANDIDATE_MULTIPLIER"
		--seed "$SEED"
		--min-genomic-separation "$MIN_GENOMIC_SEP"
	)
	if [ "$NEAREST_IF_INSUFFICIENT" = "1" ]; then
		cmd+=(--nearest-if-insufficient)
	fi
	log_command_line "${cmd[@]}"
	"${cmd[@]}" > "$FULL_RUN_ROOT/logs/generate_${slug}.log" 2>&1
}

add_condition() {
	local radius=$1 radius_slug=$2 condition=$3 init_mode=$4 init_eps=$5 init_noise=$6 init_scale=$7 dscale=$8 gamma=$9 eps_count=${10} min_sep=${11} lambda_sep=${12}
	local gslug
	gslug=$(gamma_slug "$gamma")
	local sslug
	sslug=$(sep_slug "$min_sep" "$lambda_sep")
	local config="p9016_charm3dg20k_${radius_slug}_${condition}"
	local input_pairs="$FULL_RUN_ROOT/synthetic_pairs/p9016_charm3dg20k_${radius_slug}_contacts.pairs.gz"
	{
		printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
			"$radius_slug" "$radius" "$condition" "$config" "$input_pairs" "$init_mode" "$init_eps" "$init_noise" "$init_scale" "$dscale" "$gamma" "$eps_count" "$min_sep" "$lambda_sep"
	} >> "$FULL_RUN_ROOT/condition_matrix.tsv"
}

write_condition_matrix() {
	printf '%s\n' "radius_slug	radius_pr	condition_name	config_name	input_pairs	init_mode	init_eps	init_noise_scale	init_scale	d_scale_mode	d_scale_posterior_gamma	d_scale_eps_count	min_sep_unit	lambda_sep" > "$FULL_RUN_ROOT/condition_matrix.tsv"
	local radius slug gamma gslug
	for radius in $RADII; do
		slug=$(radius_slug "$radius")
		add_condition "$radius" "$slug" "legacy_raw_ieps0p5_noise0_${REQUESTED_BACKEND}" unphased_scaffold_split 0.5 0.0 10 raw_count 1 1e-6 0 0
		add_condition "$radius" "$slug" "expected_alias_gamma1_ieps0p5_noise0" unphased_scaffold_split 0.5 0.0 10 expected_count 1 0.001 0 0
		for gamma in $GAMMAS; do
			gslug=$(gamma_slug "$gamma")
			add_condition "$radius" "$slug" "pcgamma${gslug}_sep_off_ieps0p5_noise0" unphased_scaffold_split 0.5 0.0 10 posterior_count "$gamma" 1e-6 0 0
			add_condition "$radius" "$slug" "tempered_pcgamma${gslug}_sep_off_ieps0p5_noise0" unphased_scaffold_split 0.5 0.0 10 tempered_posterior_count "$gamma" 1e-6 0 0
		done
		if [ "$RUN_SEP" = "1" ]; then
			add_condition "$radius" "$slug" "pcgamma0p5_msep1_lsep0p5_ieps0p5_noise0" unphased_scaffold_split 0.5 0.0 10 posterior_count 0.5 1e-6 1.0 0.5
			add_condition "$radius" "$slug" "pcgamma0p5_msep1_lsep1_ieps0p5_noise0" unphased_scaffold_split 0.5 0.0 10 posterior_count 0.5 1e-6 1.0 1.0
			add_condition "$radius" "$slug" "pcgamma1_msep1_lsep0p5_ieps0p5_noise0" unphased_scaffold_split 0.5 0.0 10 posterior_count 1 1e-6 1.0 0.5
			add_condition "$radius" "$slug" "pcgamma1_msep1_lsep1_ieps0p5_noise0" unphased_scaffold_split 0.5 0.0 10 posterior_count 1 1e-6 1.0 1.0
			add_condition "$radius" "$slug" "pcgamma1_msep1p5_lsep1_ieps0p5_noise0" unphased_scaffold_split 0.5 0.0 10 posterior_count 1 1e-6 1.5 1.0
			add_condition "$radius" "$slug" "raw_msep1p5_lsep1_ieps0p5_noise0" unphased_scaffold_split 0.5 0.0 10 raw_count 1 1e-6 1.5 1.0
		fi
		if [ "$RUN_COMMON_INIT" = "1" ]; then
			add_condition "$radius" "$slug" "raw_common_ieps1_noise0p05" unphased_scaffold_split 1.0 0.05 10 raw_count 1 1e-6 0 0
			add_condition "$radius" "$slug" "pcgamma1_common_ieps1_noise0p05" unphased_scaffold_split 1.0 0.05 10 posterior_count 1 1e-6 0 0
			add_condition "$radius" "$slug" "pcgamma1_common_msep1p5_lsep1" unphased_scaffold_split 1.0 0.05 10 posterior_count 1 1e-6 1.5 1.0
			add_condition "$radius" "$slug" "expected_common_msep1p5_lsep0p5" unphased_scaffold_split 1.0 0.05 10 expected_count 1 0.001 1.5 0.5
		fi
		if [ "$INCLUDE_RANDOM" = "1" ]; then
			add_condition "$radius" "$slug" "random_diploid_pcgamma1_sep_off" random_diploid 0.5 0.0 10 posterior_count 1 1e-6 0 0
			add_condition "$radius" "$slug" "random_haploid_pcgamma1_sep_off" random_haploid_split 0.5 0.0 10 posterior_count 1 1e-6 0 0
		fi
	done
}

CONFIG_PIDS=()
CONFIG_LABELS=()

wait_one_if_needed() {
	while [ "${#CONFIG_PIDS[@]}" -ge "$MAX_PARALLEL" ]; do
		local pid=${CONFIG_PIDS[0]}
		local label=${CONFIG_LABELS[0]}
		if wait "$pid"; then
			echo "CONFIG_STATUS=$label OK"
		else
			echo "CONFIG_STATUS=$label FAILED" >&2
			return 1
		fi
		CONFIG_PIDS=("${CONFIG_PIDS[@]:1}")
		CONFIG_LABELS=("${CONFIG_LABELS[@]:1}")
	done
}

wait_all_configs() {
	local failed=0
	local i
	for i in "${!CONFIG_PIDS[@]}"; do
		if wait "${CONFIG_PIDS[$i]}"; then
			echo "CONFIG_STATUS=${CONFIG_LABELS[$i]} OK"
		else
			echo "CONFIG_STATUS=${CONFIG_LABELS[$i]} FAILED" >&2
			failed=1
		fi
	done
	return "$failed"
}

run_config() {
	local radius_slug=$1 condition=$2 config=$3 input_pairs=$4 init_mode=$5 init_eps=$6 init_noise=$7 init_scale=$8 dscale=$9 gamma=${10} eps_count=${11} min_sep=${12} lambda_sep=${13}
	local sample_out_root="$FULL_RUN_ROOT/work_outputs/$radius_slug/$config"
	local work_out_dir="$sample_out_root/$config"
	local out_dir="$FULL_RUN_ROOT/outputs/$radius_slug/$config"
	local eval_dir="$FULL_RUN_ROOT/eval/$radius_slug/$config"
	local log_prefix="$FULL_RUN_ROOT/logs/${radius_slug}_${condition}"

	echo "CONFIG_START=$config"
	rm -rf "$sample_out_root" "$out_dir" "$eval_dir"
	mkdir -p "$sample_out_root" "$(dirname "$out_dir")" "$eval_dir"
	local train_cmd=(
		clean_env_prefix
		HK_BLIND_SAMPLE=P9016
		HK_BLIND_P9016_SAMPLE=P9016
		HK_BLIND_GIT_COMMIT="$GIT_COMMIT"
		HK_BLIND_GIT_DIRTY_COUNT="$GIT_DIRTY_COUNT"
		HK_BLIND_BINARY_HASH="$RUN_HASH"
		HK_BLIND_P9016_PAIRS="$input_pairs"
		HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS=1
		HK_BLIND_P9016_ALLOW_REFERENCE_DERIVED_PAIRS=1
		HK_BLIND_P9016_INPUT_CONTACT_SOURCE=charm3dg_derived_pairs
		HK_BLIND_P9016_REFERENCE_SOURCE_TDG="$TDG20"
		HK_BLIND_P9016_OUTPUT_ROOT="$sample_out_root"
		HK_BLIND_P9016_CONFIG_NAME="$config"
		HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE"
		HK_BLIND_P9016_MINIMAL_N_ITER="$N_ITER"
		HK_BLIND_P9016_MINIMAL_RELAX_STEPS="$RELAX_STEPS"
		HK_BLIND_P9016_RELAX_STEP="$RELAX_STEP"
		HK_BLIND_P9016_RELAX_BACKEND="$REQUESTED_BACKEND"
		HK_BLIND_P9016_INIT_MODE="$init_mode"
		HK_BLIND_P9016_INIT_SCALE="$init_scale"
		HK_BLIND_P9016_INIT_SEED="$SEED"
		HK_BLIND_P9016_INIT_EPS="$init_eps"
		HK_BLIND_P9016_INIT_NOISE_SCALE="$init_noise"
		HK_BLIND_P9016_MIN_SEP_UNIT="$min_sep"
		HK_BLIND_P9016_LAMBDA_SEP="$lambda_sep"
		HK_BLIND_P9016_D_SCALE_MODE="$dscale"
		HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA="$gamma"
		HK_BLIND_P9016_D_SCALE_EPS_COUNT="$eps_count"
		"$RUN_BIN"
	)
	log_command_line "${train_cmd[@]}"
	"${train_cmd[@]}" > "$log_prefix.train.log" 2>&1
	ln -s "$work_out_dir" "$out_dir"

	local audit_cmd=(clean_env_prefix HK_BLIND_SAMPLE=P9016 HK_BLIND_P9016_SAMPLE=P9016 ./audit_blind_p9016_full_cpu_output.bin "$out_dir")
	log_command_line "${audit_cmd[@]}"
	"${audit_cmd[@]}" > "$log_prefix.audit.log" 2>&1

	local eval_cmd=(
		clean_env_prefix
		HK_BLIND_SAMPLE=P9016
		HK_BLIND_P9016_SAMPLE=P9016
		HK_BLIND_ALLOW_REFERENCE_DERIVED_TRAINING=1
		HK_BLIND_P9016_EVAL_PAIRS="$EVAL_PAIRS"
		HK_BLIND_P9016_EVAL_TDG="$EVAL_TDG"
		HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE"
		scripts/p9016_common_eval.sh "$out_dir" "$eval_dir"
	)
	log_command_line "${eval_cmd[@]}"
	"${eval_cmd[@]}" > "$log_prefix.eval_wrapper.log" 2>&1
	if [ "$REQUIRE_EVAL" = "1" ] && [ ! -s "$eval_dir/summary.tsv" ]; then
		echo "required eval missing for $config" >&2
		return 1
	fi

	if [ -s "$out_dir/p9016_full.coords.tsv" ]; then
		local copytrack_cmd=(
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

launch_config() {
	local radius_slug=$1 condition=$2 config=$3 input_pairs=$4 init_mode=$5 init_eps=$6 init_noise=$7 init_scale=$8 dscale=$9 gamma=${10} eps_count=${11} min_sep=${12} lambda_sep=${13}
	wait_one_if_needed
	run_config "$radius_slug" "$condition" "$config" "$input_pairs" "$init_mode" "$init_eps" "$init_noise" "$init_scale" "$dscale" "$gamma" "$eps_count" "$min_sep" "$lambda_sep" &
	CONFIG_PIDS+=("$!")
	CONFIG_LABELS+=("$config")
}

for radius in $RADII; do
	slug=$(radius_slug "$radius")
	generate_pairs_for_radius "$radius" "$slug"
done

write_condition_matrix

while IFS=$'\t' read -r radius_slug radius condition config input_pairs init_mode init_eps init_noise init_scale dscale gamma eps_count min_sep lambda_sep; do
	launch_config "$radius_slug" "$condition" "$config" "$input_pairs" "$init_mode" "$init_eps" "$init_noise" "$init_scale" "$dscale" "$gamma" "$eps_count" "$min_sep" "$lambda_sep"
done < <(tail -n +2 "$FULL_RUN_ROOT/condition_matrix.tsv")
wait_all_configs

matrix_qc_cmd=(
	python scripts/plot_p9016_contact_matrix_qc.py
	--matrix "Observed P9016 pairs=$PAIRS"
	--out-dir "$FULL_RUN_ROOT/contact_matrix_qc"
	--bin-size "$BIN_SIZE"
	--chrom chr1
	--exclude-same-bin
	--observed-copy-matrix "Observed SNP-labeled P9016 pairs=$PAIRS"
)
for radius in $RADII; do
	slug=$(radius_slug "$radius")
	matrix_qc_cmd+=(--matrix "CHARM20k <=${radius} PR synthetic=$FULL_RUN_ROOT/synthetic_pairs/p9016_charm3dg20k_${slug}_contacts.pairs.gz")
	matrix_qc_cmd+=(--synthetic-copy-matrix "CHARM20k <=${radius} PR synthetic=$FULL_RUN_ROOT/synthetic_pairs/p9016_charm3dg20k_${slug}_contacts.copy_state.tsv.gz")
done
log_command_line "${matrix_qc_cmd[@]}"
"${matrix_qc_cmd[@]}" > "$FULL_RUN_ROOT/logs/contact_matrix_qc.log" 2>&1

log_line "+ scripts/summarize_p9016_016.py $FULL_RUN_ROOT > $FULL_RUN_ROOT/summary.tsv"
scripts/summarize_p9016_016.py "$FULL_RUN_ROOT" > "$FULL_RUN_ROOT/summary.tsv" 2> "$FULL_RUN_ROOT/logs/summarize.log"

python - <<'PY' "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" "$GIT_COMMIT" "$GIT_DIRTY_COUNT" "$RUN_HASH" "$REQUESTED_BACKEND"
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
light, git_commit, dirty, run_hash, backend = sys.argv[2:7]
rows = list(csv.DictReader(open(root / "summary.tsv"), delimiter="\t"))

def f(row, key):
    v = row.get(key, "NA")
    try:
        return float(v)
    except ValueError:
        return float("-inf")

best_all = sorted(rows, key=lambda r: f(r, "model_top1_accuracy_genome_all"), reverse=True)[:8]
best_trans = sorted(rows, key=lambda r: f(r, "model_top1_accuracy_genome_trans"), reverse=True)[:8]
best_spearman = sorted(rows, key=lambda r: f(r, "mean_per_chrom_cis_distance_spearman"), reverse=True)[:8]

with open(root / "README.md", "w") as out:
    out.write(f"# {root.name}\n\n")
    out.write("Controlled 016 sweep: reuse CHARM20k-derived 1PR/2PR unphased synthetic contacts, then test prior P9016 model-condition knobs that are already implemented and blind-safe with respect to SNP labels. Training remains reference-derived positive control because contacts are generated from CHARM/3DG.\n\n")
    out.write("## Paths\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- summary.tsv: `{root / 'summary.tsv'}`\n")
    out.write(f"- condition_matrix.tsv: `{root / 'condition_matrix.tsv'}`\n")
    out.write(f"- contact matrix QC: `{root / 'contact_matrix_qc'}`\n")
    out.write(f"- commands.log: `{root / 'commands.log'}`\n\n")
    out.write("## Boundary\n\n")
    out.write("This is not a blind baseline: training contacts are generated from CHARM/3DG. Synthetic pairs remain unphased, `uses_phase_labels=0`, and SNP labels plus P9016.1m.3DG are used only by post-training evaluation.\n\n")
    out.write("## Git And Build\n\n")
    out.write(f"- git commit at launch: `{git_commit}`\n")
    out.write(f"- git dirty count at launch: `{dirty}`\n")
    if dirty != "0":
        out.write("- WARNING: dirty tree at launch; inspect `logs/git_status.txt` and `logs/git_diff.patch` in full result root.\n")
    out.write(f"- run binary sha256: `{run_hash}`\n")
    out.write(f"- backend: `{backend}`\n\n")
    out.write("## Best By Top1 All\n\n")
    out.write("| radius | condition | top1 all | top1 cis | top1 trans | same/cross trans | cis Spearman | entropy | pU | min sep | mean sep |\n")
    out.write("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
    for r in best_all:
        out.write(f"| {r['radius_slug']} | `{r['condition_name']}` | {r['model_top1_accuracy_genome_all']} | {r['model_top1_accuracy_genome_cis']} | {r['model_top1_accuracy_genome_trans']} | {r['model_same_cross_accuracy_genome_trans']} | {r['mean_per_chrom_cis_distance_spearman']} | {r['final_mean_entropy']} | {r['final_mean_pU']} | {r['final_min_sep']} | {r['final_mean_sep']} |\n")
    out.write("\n## Best By Top1 Trans\n\n")
    out.write("| radius | condition | top1 trans | top1 all | top1 cis | same/cross trans | cis Spearman | entropy | pU | copytrack cos<0 |\n")
    out.write("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
    for r in best_trans:
        out.write(f"| {r['radius_slug']} | `{r['condition_name']}` | {r['model_top1_accuracy_genome_trans']} | {r['model_top1_accuracy_genome_all']} | {r['model_top1_accuracy_genome_cis']} | {r['model_same_cross_accuracy_genome_trans']} | {r['mean_per_chrom_cis_distance_spearman']} | {r['final_mean_entropy']} | {r['final_mean_pU']} | {r['copytrack_frac_cos_lt_0']} |\n")
    out.write("\n## Best By Cis Spearman\n\n")
    out.write("| radius | condition | cis Spearman | top1 all | top1 cis | top1 trans | entropy | pU |\n")
    out.write("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |\n")
    for r in best_spearman:
        out.write(f"| {r['radius_slug']} | `{r['condition_name']}` | {r['mean_per_chrom_cis_distance_spearman']} | {r['model_top1_accuracy_genome_all']} | {r['model_top1_accuracy_genome_cis']} | {r['model_top1_accuracy_genome_trans']} | {r['final_mean_entropy']} | {r['final_mean_pU']} |\n")
    out.write("\n## Interpretation Rules\n\n")
    out.write("- A condition is useful only if it improves trans without sacrificing cis Spearman, coverage, entropy/pU, or copytrack continuity.\n")
    out.write("- Any positive signal here remains a reference-derived positive control, not a blind production baseline.\n")
    out.write("- If no condition beats 015 robustly, the bottleneck is not fixed by previously implemented dscale/init/sep knobs on clean contacts.\n")
PY

log_line "+ scripts/p9016_publish_light_result.sh $FULL_RUN_ROOT $LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "CONDITION_MATRIX_TSV=$FULL_RUN_ROOT/condition_matrix.tsv"
