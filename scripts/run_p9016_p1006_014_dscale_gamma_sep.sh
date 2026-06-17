#!/usr/bin/env bash
set -euo pipefail

if [ "${CONDA_DEFAULT_ENV:-}" != "analysis" ]; then
	echo "error: activate conda env 'analysis' before running this script" >&2
	exit 2
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

RUN_ID="014-$(date +%Y%m%d_%H%M%S)-p9016_p1006_dscale_gamma_sep_1m"
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$REPO_ROOT/test_res}/$RUN_ID"
mkdir -p "$FULL_RUN_ROOT"/{logs,outputs,eval,work_outputs,scripts_snapshot}
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)
OUTPUT_ROOT="$FULL_RUN_ROOT/outputs"
WORK_OUTPUT_ROOT="$FULL_RUN_ROOT/work_outputs"
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
	echo "cppflags	${CPPFLAGS:-default}"
	echo "ldflags	${LDFLAGS:-default}"
	echo "conda_env	${CONDA_DEFAULT_ENV:-NA}"
	echo "path	$PATH"
} > "$FULL_RUN_ROOT/logs/build_env.txt"

if [ "${HK_BLIND_REQUIRE_CLEAN_TREE:-0}" = "1" ] && [ "$GIT_DIRTY_COUNT" != "0" ]; then
	echo "error: hickit tree is dirty and HK_BLIND_REQUIRE_CLEAN_TREE=1" >&2
	exit 1
fi

cp "$0" "$FULL_RUN_ROOT/scripts_snapshot/"
cp scripts/p9016_common_eval.sh \
   scripts/p9016_publish_light_result.sh \
   scripts/compute_p9016_copytrack_diag.py \
   scripts/summarize_p9016_p1006_014.py \
   scripts/summarize_p9016_p1006_014_deltas.py \
   scripts/smoke_p9016_dscale_gamma.sh \
   "$FULL_RUN_ROOT/scripts_snapshot/"

with_commands_lock() {
	local delay=0.05
	while ! mkdir "$COMMANDS_LOCK_DIR" 2>/dev/null; do
		sleep "$delay"
	done
	"$@" >> "$COMMANDS_LOG"
	rmdir "$COMMANDS_LOCK_DIR"
}

log_cmd() {
	local line="+"
	local arg
	for arg in "$@"; do
		line+=" $(printf '%q' "$arg")"
	done
	with_commands_lock printf '%s\n' "$line"
	"$@"
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

clean_env_prefix() {
	local env_cmd=(
		env -i
		PATH="$PATH"
		LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
		HOME="${HOME:-}"
		CONDA_DEFAULT_ENV="${CONDA_DEFAULT_ENV:-}"
		CONDA_PREFIX="${CONDA_PREFIX:-}"
	)
	# Passing CUDA_VISIBLE_DEVICES as an empty string hides all GPUs.
	if [ "${CUDA_VISIBLE_DEVICES+x}" = "x" ]; then
		env_cmd+=(CUDA_VISIBLE_DEVICES="$CUDA_VISIBLE_DEVICES")
	fi
	if [ "${CUDA_DEVICE_ORDER+x}" = "x" ]; then
		env_cmd+=(CUDA_DEVICE_ORDER="$CUDA_DEVICE_ORDER")
	fi
	if [ "${NVIDIA_VISIBLE_DEVICES+x}" = "x" ]; then
		env_cmd+=(NVIDIA_VISIBLE_DEVICES="$NVIDIA_VISIBLE_DEVICES")
	fi
	if [ "${NVIDIA_DRIVER_CAPABILITIES+x}" = "x" ]; then
		env_cmd+=(NVIDIA_DRIVER_CAPABILITIES="$NVIDIA_DRIVER_CAPABILITIES")
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

sample_pairs_default() {
	case "$1" in
	P9016) echo "/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz" ;;
	P1006) echo "/shared/zliu/CHARM/CHARM_mesc/data/pairs/P1006.pairs.gz" ;;
	*) echo "/shared/zliu/CHARM/CHARM_mesc/data/pairs/$1.pairs.gz" ;;
	esac
}

sample_tdg_default() {
	case "$1" in
	P9016) echo "/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz" ;;
	P1006) echo "/shared/zliu/CHARM/CHARM_mesc/data/tdg/P1006.1m.3dg.gz" ;;
	*) echo "/shared/zliu/CHARM/CHARM_mesc/data/tdg/$1.1m.3dg.gz" ;;
	esac
}

resolve_existing_path() {
	local p=$1
	if [ -s "$p" ]; then
		cd "$(dirname "$p")" && printf '%s/%s\n' "$(pwd)" "$(basename "$p")"
	elif [ -s "$REPO_ROOT/$p" ]; then
		cd "$(dirname "$REPO_ROOT/$p")" && printf '%s/%s\n' "$(pwd)" "$(basename "$p")"
	else
		echo "$p"
	fi
}

gamma_slug() {
	case "$1" in
	0) echo "0" ;;
	0.25) echo "0p25" ;;
	0.5) echo "0p5" ;;
	0.75) echo "0p75" ;;
	1|1.0) echo "1" ;;
	*) echo "$1" | tr '.' 'p' ;;
	esac
}

backend_available() {
	if [ "$1" = "gpu" ]; then
		command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1
	else
		return 0
	fi
}

SAMPLES=${HK_BLIND_014_SAMPLES:-"P9016 P1006"}
GAMMAS=${HK_BLIND_014_GAMMAS:-"0 0.25 0.5 0.75 1"}
RUN_SEP=${HK_BLIND_014_RUN_SEP:-1}
MAX_PARALLEL=${HK_BLIND_014_MAX_PARALLEL:-2}
case "$MAX_PARALLEL" in
	''|*[!0-9]*)
		echo "error: HK_BLIND_014_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_014_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-0.012}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-0}
ALLOW_CUSTOM=${HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS:-0}
REQUESTED_BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-auto}
if [ "$REQUESTED_BACKEND" = "auto" ]; then
	if backend_available gpu; then
		REQUESTED_BACKEND=gpu
	else
		REQUESTED_BACKEND=cpu
	fi
fi
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
	echo "samples	$SAMPLES"
	echo "gammas	$GAMMAS"
	echo "run_sep	$RUN_SEP"
	echo "max_parallel	$MAX_PARALLEL"
	echo "backend	$REQUESTED_BACKEND"
	echo "bin_size_bp	$BIN_SIZE"
	echo "n_iter	$N_ITER"
	echo "relax_steps	$RELAX_STEPS"
	echo "relax_step	$RELAX_STEP"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"
{
	echo "run_blind_p9016_minimal_sha256	$RUN_HASH"
	echo "audit_blind_p9016_full_cpu_output_sha256	$AUDIT_HASH"
} >> "$FULL_RUN_ROOT/logs/build_env.txt"

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
	local sample=$1
	local config=$2
	local gamma=$3
	local min_sep=$4
	local lambda_sep=$5
	local pairs=$6
	local eval_pairs=$7
	local eval_tdg=$8
	local allow_custom=$9
	local sample_out_root="$WORK_OUTPUT_ROOT/$sample/$config"
	local work_out_dir="$sample_out_root/$config"
	local out_dir="$OUTPUT_ROOT/$sample/$config"
	local eval_dir="$FULL_RUN_ROOT/eval/$sample/$config"
	local log_prefix="$FULL_RUN_ROOT/logs/${sample}_${config}"
	local run_eval=1
	local train_cmd
	local audit_cmd
	local eval_cmd
	local copytrack_cmd

	echo "CONFIG_START=$sample/$config"
	rm -rf "$sample_out_root" "$out_dir" "$eval_dir"
	mkdir -p "$sample_out_root" "$(dirname "$out_dir")" "$eval_dir"
	train_cmd=(
		clean_env_prefix
		HK_BLIND_SAMPLE="$sample" \
		HK_BLIND_P9016_SAMPLE="$sample" \
		HK_BLIND_GIT_COMMIT="$GIT_COMMIT" \
		HK_BLIND_GIT_DIRTY_COUNT="$GIT_DIRTY_COUNT" \
		HK_BLIND_BINARY_HASH="$RUN_HASH" \
		HK_BLIND_P9016_PAIRS="$pairs" \
		HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS="$allow_custom" \
		HK_BLIND_P9016_OUTPUT_ROOT="$sample_out_root" \
		HK_BLIND_P9016_CONFIG_NAME="$config" \
		HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE" \
		HK_BLIND_P9016_MINIMAL_N_ITER="$N_ITER" \
		HK_BLIND_P9016_MINIMAL_RELAX_STEPS="$RELAX_STEPS" \
		HK_BLIND_P9016_RELAX_STEP="$RELAX_STEP" \
		HK_BLIND_P9016_RELAX_BACKEND="$REQUESTED_BACKEND" \
		HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split \
		HK_BLIND_P9016_INIT_SCALE=10 \
		HK_BLIND_P9016_INIT_SEED=17 \
		HK_BLIND_P9016_INIT_EPS=0.5 \
		HK_BLIND_P9016_INIT_NOISE_SCALE=0.0 \
		HK_BLIND_P9016_MIN_SEP_UNIT="$min_sep" \
		HK_BLIND_P9016_LAMBDA_SEP="$lambda_sep" \
		HK_BLIND_P9016_D_SCALE_MODE=posterior_count \
		HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA="$gamma" \
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6 \
		"$RUN_BIN"
	)
	log_command_line "${train_cmd[@]}"
	"${train_cmd[@]}" > "$log_prefix.train.log" 2>&1

	ln -s "$work_out_dir" "$out_dir"
	audit_cmd=(
		clean_env_prefix
		HK_BLIND_SAMPLE="$sample"
		HK_BLIND_P9016_SAMPLE="$sample"
		./audit_blind_p9016_full_cpu_output.bin "$out_dir"
	)
	log_command_line "${audit_cmd[@]}"
	"${audit_cmd[@]}" > "$log_prefix.audit.log" 2>&1

	if [ "${HK_BLIND_P9016_RUN_EVAL:-1}" = "0" ]; then
		run_eval=0
	fi
	if [ "$run_eval" = "1" ]; then
		eval_cmd=(
			clean_env_prefix
			HK_BLIND_SAMPLE="$sample" \
			HK_BLIND_P9016_SAMPLE="$sample" \
			HK_BLIND_P9016_EVAL_PAIRS="$eval_pairs" \
			HK_BLIND_P9016_EVAL_TDG="$eval_tdg" \
			HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE" \
			scripts/p9016_common_eval.sh "$out_dir" "$eval_dir"
		)
		log_command_line "${eval_cmd[@]}"
		"${eval_cmd[@]}" > "$log_prefix.eval_wrapper.log" 2>&1
		if [ "$REQUIRE_EVAL" = "1" ] && [ ! -s "$eval_dir/summary.tsv" ]; then
			echo "required eval missing for $sample/$config" >&2
			return 1
		fi
	else
		{
			echo "status	SKIPPED_EVAL_SMOKE_MODE"
			echo "sample	$sample"
			echo "reason	HK_BLIND_P9016_RUN_EVAL=0"
		} > "$eval_dir/eval_status.tsv"
		echo "SKIPPED_EVAL_SMOKE_MODE" > "$eval_dir/eval.log"
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
	echo "CONFIG_DONE=$sample/$config"
}

launch_config() {
	local sample=$1
	local config=$2
	shift 2
	wait_one_if_needed
	run_config "$sample" "$config" "$@" &
	CONFIG_PIDS+=("$!")
	CONFIG_LABELS+=("$sample/$config")
}

for sample in $SAMPLES; do
	pairs_var="HK_BLIND_${sample}_PAIRS"
	eval_pairs_var="HK_BLIND_${sample}_EVAL_PAIRS"
	eval_tdg_var="HK_BLIND_${sample}_EVAL_TDG"
	pairs=$(resolve_existing_path "${!pairs_var:-${HK_BLIND_P9016_PAIRS:-$(sample_pairs_default "$sample")}}")
	eval_pairs=$(resolve_existing_path "${!eval_pairs_var:-${HK_BLIND_P9016_EVAL_PAIRS:-$(sample_pairs_default "$sample")}}")
	eval_tdg=$(resolve_existing_path "${!eval_tdg_var:-${HK_BLIND_P9016_EVAL_TDG:-$(sample_tdg_default "$sample")}}")
	allow_custom=$ALLOW_CUSTOM
	if [ "$pairs" != "$(sample_pairs_default P9016)" ]; then
		allow_custom=1
	fi
	for gamma in $GAMMAS; do
		slug=$(gamma_slug "$gamma")
		config="${sample,,}_pcgamma${slug}_sep_off_ieps0p5_noise0_seed17"
		launch_config "$sample" "$config" "$gamma" 0 0 "$pairs" "$eval_pairs" "$eval_tdg" "$allow_custom"
	done
	if [ "$RUN_SEP" = "1" ]; then
		for spec in \
			"0.5 1 0.5" \
			"0.5 1 1" \
			"1 1 0.5" \
			"1 1 1" \
			"1 1.5 1"
		do
			set -- $spec
			gamma=$1
			min_sep=$2
			lambda_sep=$3
			gslug=$(gamma_slug "$gamma")
			mslug=$(gamma_slug "$min_sep")
			lslug=$(gamma_slug "$lambda_sep")
			config="${sample,,}_pcgamma${gslug}_msep${mslug}_lsep${lslug}_ieps0p5_noise0_seed17"
			launch_config "$sample" "$config" "$gamma" "$min_sep" "$lambda_sep" "$pairs" "$eval_pairs" "$eval_tdg" "$allow_custom"
		done
	fi
done

wait_all_configs

log_line "+ scripts/summarize_p9016_p1006_014.py $FULL_RUN_ROOT > $FULL_RUN_ROOT/summary.tsv"
scripts/summarize_p9016_p1006_014.py "$FULL_RUN_ROOT" > "$FULL_RUN_ROOT/summary.tsv" 2> "$FULL_RUN_ROOT/logs/summarize.log"
log_line "+ scripts/summarize_p9016_p1006_014_deltas.py $FULL_RUN_ROOT"
scripts/summarize_p9016_p1006_014_deltas.py "$FULL_RUN_ROOT" 2> "$FULL_RUN_ROOT/logs/deltas.log"

python - <<'PY' "$FULL_RUN_ROOT/summary.tsv" "$FULL_RUN_ROOT/copytrack_summary.tsv"
import csv, sys
summary, out = sys.argv[1], sys.argv[2]
cols = ["sample","config_name","copytrack_sep_p05","copytrack_sep_median","copytrack_vector_cos_p05","copytrack_vector_cos_median","copytrack_frac_cos_lt_0","copytrack_frac_projection_sign_switch"]
with open(summary, newline="") as f, open(out, "w", newline="") as g:
    rows = list(csv.DictReader(f, delimiter="\t"))
    w = csv.DictWriter(g, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "NA") for c in cols})
PY

python - <<'PY' "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" "$GIT_COMMIT" "$GIT_DIRTY_COUNT" "$RUN_HASH" "$REQUESTED_BACKEND"
import csv, math, sys
from pathlib import Path
root = Path(sys.argv[1])
light = sys.argv[2]
git_commit = sys.argv[3]
dirty = sys.argv[4]
run_hash = sys.argv[5]
backend = sys.argv[6]
rows = list(csv.DictReader(open(root / "summary.tsv"), delimiter="\t"))
gamma_rows = list(csv.DictReader(open(root / "gamma_delta.tsv"), delimiter="\t")) if (root / "gamma_delta.tsv").exists() else []
sep_rows = list(csv.DictReader(open(root / "sep_delta.tsv"), delimiter="\t")) if (root / "sep_delta.tsv").exists() else []

def flt(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None

def best(rows, sample, metric):
    candidates = [r for r in rows if r.get("sample") == sample and flt(r.get(metric)) is not None]
    if not candidates:
        return "NA"
    return max(candidates, key=lambda r: flt(r[metric])).get("gamma", "NA")

with open(root / "README.md", "w") as out:
    out.write(f"# {root.name}\n\n")
    out.write("Posterior-count d_scale gamma sweep across P9016 and P1006, plus modest homolog-separation guardrail tests at 1 Mb.\n\n")
    out.write("## Paths\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- lightweight hickit result: `{light}`\n")
    for name in ["summary.tsv", "gamma_delta.tsv", "sep_delta.tsv", "sample_delta.tsv", "commands.log"]:
        out.write(f"- {name}: `{root / name}`\n")
    out.write("\n## Training Boundary\n\n")
    out.write("Training uses only raw sample pairs. Phase labels and CHARM/3DG are eval-only through `scripts/p9016_common_eval.sh`. Manifests should report `uses_phase_labels=0` and `uses_charm_or_reference=0`.\n\n")
    out.write("## Git And Build\n\n")
    out.write(f"- git commit: `{git_commit}`\n")
    out.write(f"- git dirty count at start: `{dirty}`\n")
    if dirty != "0":
        out.write("- WARNING: dirty tree at launch; inspect `logs/git_status.txt` and `logs/git_diff.patch` in full test_res.\n")
    out.write(f"- run binary sha256: `{run_hash}`\n")
    out.write(f"- backend: `{backend}`\n\n")
    out.write("## Configs\n\n")
    out.write("| sample | config | gamma | min_sep_unit | lambda_sep | backend |\n")
    out.write("| --- | --- | ---: | ---: | ---: | --- |\n")
    for r in rows:
        out.write(f"| {r['sample']} | `{r['config_name']}` | {r['d_scale_posterior_gamma']} | {r['min_sep_unit']} | {r['lambda_sep']} | {r['backend']} |\n")
    out.write("\n## Main Results\n\n")
    out.write("| sample | config | gamma | min_sep | lambda_sep | top1 all | top1 cis | top1 trans | same/cross cis | cis Spearman | entropy | pU | min sep | sep_p05 | mean sep | cos<0 |\n")
    out.write("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
    for r in rows:
        out.write(f"| {r['sample']} | `{r['config_name']}` | {r['d_scale_posterior_gamma']} | {r['min_sep_unit']} | {r['lambda_sep']} | {r['model_top1_accuracy_genome_all']} | {r['model_top1_accuracy_genome_cis']} | {r['model_top1_accuracy_genome_trans']} | {r['model_same_cross_accuracy_genome_cis']} | {r['mean_per_chrom_cis_distance_spearman']} | {r['final_mean_entropy']} | {r['final_mean_pU']} | {r['final_min_sep']} | {r['sep_p05']} | {r['final_mean_sep']} | {r['copytrack_frac_cos_lt_0']} |\n")
    out.write("\n## Interpretation\n\n")
    out.write(f"- Best P9016 gamma by top1_all: `{best(gamma_rows, 'P9016', 'top1_all')}`.\n")
    out.write(f"- Best P1006 gamma by top1_all: `{best(gamma_rows, 'P1006', 'top1_all')}`.\n")
    out.write(f"- Best P9016 gamma by cis Spearman: `{best(gamma_rows, 'P9016', 'cis_spearman')}`.\n")
    out.write(f"- Best P1006 gamma by cis Spearman: `{best(gamma_rows, 'P1006', 'cis_spearman')}`.\n")
    out.write("- Treat gamma=1 as the historical posterior-count baseline. If another gamma is not better by more than about 0.005 top1_all, prefer gamma=1 for continuity.\n")
    out.write("- Separation guardrails should be interpreted as geometry controls: accept them only if sep_p05/min separation improves without a meaningful top1 or Spearman cost.\n")
    out.write("- This run does not introduce priors, entropy-aware rho, annealing, sharpening, random init, copytrack smoothness, oracle training, 200 kb, or multiresolution changes.\n")
PY

log_line "+ scripts/p9016_publish_light_result.sh $FULL_RUN_ROOT $LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "GAMMA_DELTA_TSV=$FULL_RUN_ROOT/gamma_delta.tsv"
echo "SEP_DELTA_TSV=$FULL_RUN_ROOT/sep_delta.tsv"
echo "SAMPLE_DELTA_TSV=$FULL_RUN_ROOT/sample_delta.tsv"
