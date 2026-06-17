#!/usr/bin/env bash
set -euo pipefail

if [ "${CONDA_DEFAULT_ENV:-}" != "analysis" ]; then
	echo "error: activate conda env 'analysis' before running this script" >&2
	exit 2
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

RUN_ID="012-$(date +%Y%m%d_%H%M%S)-p9016_dscale_semantics_regression_1m"
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$REPO_ROOT/test_res}/$RUN_ID"
mkdir -p "$FULL_RUN_ROOT"/{logs,outputs,eval,plots,scripts_snapshot}
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)
OUTPUT_ROOT="$FULL_RUN_ROOT/outputs"
WORK_OUTPUT_ROOT="$FULL_RUN_ROOT/work_outputs"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"
COMMANDS_LOCK_DIR="$FULL_RUN_ROOT/.commands.lock"

exec > >(tee -a "$FULL_RUN_ROOT/logs/runner.stdout.log") 2> >(tee -a "$FULL_RUN_ROOT/logs/runner.stderr.log" >&2)

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "OUTPUT_ROOT=$OUTPUT_ROOT"
echo "WORK_OUTPUT_ROOT=$WORK_OUTPUT_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

GIT_COMMIT=$(git rev-parse HEAD)
GIT_DIRTY_COUNT=$(git status --porcelain | wc -l | awk '{print $1}')
{
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"
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

cp "$0" "$FULL_RUN_ROOT/logs/run_p9016_012_dscale_semantics_regression.sh"
cp scripts/p9016_common_eval.sh \
   scripts/p9016_publish_light_result.sh \
   scripts/compute_p9016_copytrack_diag.py \
   scripts/summarize_p9016_model_sweep.py \
   scripts/audit_p9016_dscale_semantics.py \
   scripts/smoke_p9016_dscale_modes.sh \
   "$FULL_RUN_ROOT/scripts_snapshot/"

log_cmd() {
	with_commands_lock printf '+ '
	with_commands_lock printf '%q ' "$@"
	with_commands_lock printf '\n'
	"$@"
}

with_commands_lock() {
	local delay=0.05
	while ! mkdir "$COMMANDS_LOCK_DIR" 2>/dev/null; do
		sleep "$delay"
	done
	"$@" >> "$COMMANDS_LOG"
	rmdir "$COMMANDS_LOCK_DIR"
}

log_line() {
	with_commands_lock printf '%s\n' "$*"
}

log_cmd make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
log_cmd make smoke_blind_p9016_minimal
log_cmd scripts/smoke_p9016_dscale_modes.sh

hash_file() {
	if command -v sha256sum >/dev/null 2>&1; then
		sha256sum "$1" | awk '{print $1}'
	else
		echo NA
	fi
}

CPU_RUN_BIN="$FULL_RUN_ROOT/logs/run_blind_p9016_minimal.cpu.bin"
cp run_blind_p9016_minimal.bin "$CPU_RUN_BIN"
chmod +x "$CPU_RUN_BIN"
CPU_RUN_HASH=$(hash_file "$CPU_RUN_BIN")
AUDIT_HASH=$(hash_file audit_blind_p9016_full_cpu_output.bin)

resolve_pairs_path() {
	local p=${HK_BLIND_P9016_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
	if [ -s "$p" ]; then
		cd "$(dirname "$p")" && printf '%s/%s\n' "$(pwd)" "$(basename "$p")"
	elif [ -s "$REPO_ROOT/$p" ]; then
		cd "$(dirname "$REPO_ROOT/$p")" && printf '%s/%s\n' "$(pwd)" "$(basename "$p")"
	else
		echo "$p"
	fi
}

PAIRS_PATH=$(resolve_pairs_path)
ALLOW_CUSTOM=${HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS:-0}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-0.012}
RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-1}
if [ "$ALLOW_CUSTOM" = "1" ] || [[ "$PAIRS_PATH" == *p9016_blind_smoke.pairs ]]; then
	RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-0}
fi

GPU_AVAILABLE=0
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
	GPU_AVAILABLE=1
fi
PREFERRED_BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-auto}
if [ "$PREFERRED_BACKEND" = "auto" ]; then
	if [ "$GPU_AVAILABLE" = "1" ]; then
		PREFERRED_BACKEND=gpu
	else
		PREFERRED_BACKEND=cpu
	fi
fi

if [ "$GPU_AVAILABLE" = "1" ]; then
	log_cmd make gpu=1 run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
	GPU_RUN_BIN="$FULL_RUN_ROOT/logs/run_blind_p9016_minimal.gpu.bin"
	cp run_blind_p9016_minimal.bin "$GPU_RUN_BIN"
	chmod +x "$GPU_RUN_BIN"
	GPU_RUN_HASH=$(hash_file "$GPU_RUN_BIN")
	AUDIT_HASH=$(hash_file audit_blind_p9016_full_cpu_output.bin)
else
	GPU_RUN_BIN="$CPU_RUN_BIN"
	GPU_RUN_HASH="$CPU_RUN_HASH"
fi

RUN_HASH="$GPU_RUN_HASH"
if [ "$GPU_AVAILABLE" != "1" ]; then
	RUN_HASH="$CPU_RUN_HASH"
fi
{
	echo "run_blind_p9016_minimal_sha256	$RUN_HASH"
	echo "run_blind_p9016_minimal_cpu_sha256	$CPU_RUN_HASH"
	echo "run_blind_p9016_minimal_gpu_sha256	$GPU_RUN_HASH"
	echo "audit_blind_p9016_full_cpu_output_sha256	$AUDIT_HASH"
} >> "$FULL_RUN_ROOT/run_manifest.tsv"
{
	echo "run_blind_p9016_minimal_sha256	$RUN_HASH"
	echo "run_blind_p9016_minimal_cpu_sha256	$CPU_RUN_HASH"
	echo "run_blind_p9016_minimal_gpu_sha256	$GPU_RUN_HASH"
	echo "audit_blind_p9016_full_cpu_output_sha256	$AUDIT_HASH"
} >> "$FULL_RUN_ROOT/logs/build_env.txt"

{
	echo "pairs_path	$PAIRS_PATH"
	echo "allow_custom_pairs	$ALLOW_CUSTOM"
	echo "bin_size_bp	$BIN_SIZE"
	echo "n_iter	$N_ITER"
	echo "relax_steps	$RELAX_STEPS"
	echo "relax_step	$RELAX_STEP"
	echo "run_eval	$RUN_EVAL"
	echo "gpu_available	$GPU_AVAILABLE"
	echo "preferred_backend	$PREFERRED_BACKEND"
} >> "$FULL_RUN_ROOT/run_manifest.tsv"

append_config_summary() {
	local config=$1
	local out_dir="$OUTPUT_ROOT/$config"
	if [ -s "$out_dir/p9016_full.manifest.tsv" ]; then
		:
	else
		return 0
	fi
}

run_config() {
	local config=$1
	local backend=$2
	local init_eps=$3
	local init_noise=$4
	local min_sep=$5
	local lambda_sep=$6
	local dscale_mode=$7
	local run_bin=$8
	local run_hash=$9
	local log_prefix="$FULL_RUN_ROOT/logs/$config"
	local work_root="$WORK_OUTPUT_ROOT/$config"
	local work_out_dir="$work_root/$config"
	local out_dir="$OUTPUT_ROOT/$config"
	local eval_dir="$FULL_RUN_ROOT/eval/$config"

	echo "CONFIG_START=$config"
	rm -rf "$work_root" "$out_dir"
	mkdir -p "$work_root"
	log_line "+ env HK_BLIND_P9016_CONFIG_NAME=$config HK_BLIND_P9016_RELAX_BACKEND=$backend HK_BLIND_P9016_D_SCALE_MODE=$dscale_mode ... run_blind"
	env \
		HK_BLIND_GIT_COMMIT="$GIT_COMMIT" \
		HK_BLIND_GIT_DIRTY_COUNT="$GIT_DIRTY_COUNT" \
		HK_BLIND_BINARY_HASH="$run_hash" \
		HK_BLIND_P9016_PAIRS="$PAIRS_PATH" \
		HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS="$ALLOW_CUSTOM" \
		HK_BLIND_P9016_OUTPUT_ROOT="$work_root" \
		HK_BLIND_P9016_CONFIG_NAME="$config" \
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
		HK_BLIND_P9016_D_SCALE_MODE="$dscale_mode" \
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6 \
		HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA=1.0 \
		"$run_bin" > "$log_prefix.train.log" 2>&1

	ln -s "$work_out_dir" "$out_dir"
	log_line "+ ./audit_blind_p9016_full_cpu_output.bin $out_dir"
	./audit_blind_p9016_full_cpu_output.bin "$out_dir" > "$log_prefix.audit.log" 2>&1

	mkdir -p "$eval_dir"
	if [ "$RUN_EVAL" = "1" ]; then
		log_line "+ scripts/p9016_common_eval.sh $out_dir $eval_dir"
		env \
			HK_BLIND_P9016_EVAL_PAIRS="${HK_BLIND_P9016_EVAL_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}" \
			HK_BLIND_P9016_EVAL_TDG="${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}" \
			HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE" \
			scripts/p9016_common_eval.sh "$out_dir" "$eval_dir" > "$log_prefix.eval_wrapper.log" 2>&1
	else
		{
			echo "status	SKIPPED_EVAL_SMOKE_MODE"
			echo "reason	custom_or_smoke_pairs"
		} > "$eval_dir/eval_status.tsv"
		echo "SKIPPED_EVAL_SMOKE_MODE" > "$eval_dir/eval.log"
	fi
	if [ -s "$out_dir/p9016_full.coords.tsv" ]; then
		log_line "+ python scripts/compute_p9016_copytrack_diag.py --config-output-dir $out_dir --eval-output-dir $eval_dir --out $eval_dir/copytrack_vector_diag.tsv"
		python scripts/compute_p9016_copytrack_diag.py \
			--config-output-dir "$out_dir" \
			--eval-output-dir "$eval_dir" \
			--out "$eval_dir/copytrack_vector_diag.tsv" > "$log_prefix.copytrack.log" 2>&1
	fi
	append_config_summary "$config"
	echo "CONFIG_DONE=$config"
}

CONFIG_PIDS=()
CONFIG_NAMES=()

launch_config() {
	local config=$1
	shift
	run_config "$config" "$@" &
	CONFIG_PIDS+=("$!")
	CONFIG_NAMES+=("$config")
}

wait_for_configs() {
	local failed=0
	local i pid config
	for i in "${!CONFIG_PIDS[@]}"; do
		pid=${CONFIG_PIDS[$i]}
		config=${CONFIG_NAMES[$i]}
		if wait "$pid"; then
			echo "CONFIG_STATUS=$config OK"
		else
			echo "CONFIG_STATUS=$config FAILED" >&2
			failed=1
		fi
	done
	return "$failed"
}

launch_config current_raw_cpu_exact_ieps0p5_noise0_seed17 cpu 0.5 0.0 0.0 0.0 raw_count "$CPU_RUN_BIN" "$CPU_RUN_HASH"
launch_config posterior_count_cpu_exact_ieps0p5_noise0_seed17 cpu 0.5 0.0 0.0 0.0 posterior_count "$CPU_RUN_BIN" "$CPU_RUN_HASH"

if [ "$GPU_AVAILABLE" = "1" ]; then
	launch_config current_raw_gpu_exact_ieps0p5_noise0_seed17 gpu 0.5 0.0 0.0 0.0 raw_count "$GPU_RUN_BIN" "$GPU_RUN_HASH"
	launch_config posterior_count_gpu_exact_ieps0p5_noise0_seed17 gpu 0.5 0.0 0.0 0.0 posterior_count "$GPU_RUN_BIN" "$GPU_RUN_HASH"
else
	for config in current_raw_gpu_exact_ieps0p5_noise0_seed17 posterior_count_gpu_exact_ieps0p5_noise0_seed17; do
		mkdir -p "$OUTPUT_ROOT/$config" "$FULL_RUN_ROOT/eval/$config"
		echo "status	SKIPPED_GPU_UNAVAILABLE" > "$FULL_RUN_ROOT/eval/$config/eval_status.tsv"
	done
fi

if [ "${HK_BLIND_P9016_012_CORE_ONLY:-0}" != "1" ]; then
	PREFERRED_RUN_BIN="$CPU_RUN_BIN"
	PREFERRED_RUN_HASH="$CPU_RUN_HASH"
	if [ "$PREFERRED_BACKEND" = "gpu" ]; then
		PREFERRED_RUN_BIN="$GPU_RUN_BIN"
		PREFERRED_RUN_HASH="$GPU_RUN_HASH"
	fi
	launch_config posterior_count_common_sep_off_ieps1_noise0p05_seed17 "$PREFERRED_BACKEND" 1.0 0.05 0.0 0.0 posterior_count "$PREFERRED_RUN_BIN" "$PREFERRED_RUN_HASH"
	launch_config posterior_count_common_msep1p5_lsep1_ieps1_noise0p05_seed17 "$PREFERRED_BACKEND" 1.0 0.05 1.5 1.0 posterior_count "$PREFERRED_RUN_BIN" "$PREFERRED_RUN_HASH"
fi

wait_for_configs

log_line "+ scripts/summarize_p9016_model_sweep.py $FULL_RUN_ROOT > $FULL_RUN_ROOT/summary.tsv"
scripts/summarize_p9016_model_sweep.py "$FULL_RUN_ROOT" > "$FULL_RUN_ROOT/summary.tsv" 2> "$FULL_RUN_ROOT/logs/summarize.log"
log_line "+ scripts/audit_p9016_dscale_semantics.py $FULL_RUN_ROOT > $FULL_RUN_ROOT/dscale_semantics_audit.tsv"
scripts/audit_p9016_dscale_semantics.py "$FULL_RUN_ROOT" > "$FULL_RUN_ROOT/dscale_semantics_audit.tsv" 2> "$FULL_RUN_ROOT/logs/dscale_semantics_audit.log"

python - <<'PY' "$FULL_RUN_ROOT/summary.tsv" "$FULL_RUN_ROOT/copytrack_summary.tsv"
import csv, sys
summary, out = sys.argv[1], sys.argv[2]
cols = ["config_name","backend","copytrack_sep_p05","copytrack_sep_median","copytrack_vector_cos_p05","copytrack_vector_cos_median","copytrack_frac_cos_lt_0","copytrack_frac_projection_sign_switch"]
with open(summary, newline="") as f, open(out, "w", newline="") as g:
    rows = list(csv.DictReader(f, delimiter="\t"))
    w = csv.DictWriter(g, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "NA") for c in cols})
PY

HEADLINE=$(cat "$FULL_RUN_ROOT/dscale_semantics_headline.txt")

python - <<'PY' "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" "$HEADLINE" "$GIT_DIRTY_COUNT"
import csv, sys
from pathlib import Path
root = Path(sys.argv[1])
light = sys.argv[2]
headline = sys.argv[3]
dirty = sys.argv[4]
rows = list(csv.DictReader(open(root / "summary.tsv"), delimiter="\t"))
with open(root / "README.md", "w") as out:
    out.write(f"# {root.name}\n\n")
    out.write("This controlled run tests whether historical 003/004 1 Mb scaffold softall behavior corresponds to posterior-count d_scale semantics.\n\n")
    out.write("## Paths\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- lightweight hickit result: `{light}`\n")
    out.write(f"- summary TSV: `{root / 'summary.tsv'}`\n")
    out.write(f"- d_scale semantics audit TSV: `{root / 'dscale_semantics_audit.tsv'}`\n")
    out.write(f"- d_scale semantics headline: `{headline}`\n\n")
    out.write("## Tree State\n\n")
    out.write(f"- hickit git dirty entries at start: `{dirty}`\n")
    if dirty != "0":
        out.write("- WARNING: this run started from a dirty tree; inspect `logs/git_status.txt` and `logs/git_diff.patch`.\n")
    out.write("\n## Training Boundary\n\n")
    out.write("Training reads only raw P9016 pairs. Phase labels and CHARM/3DG are used only by post-training eval.\n\n")
    out.write("## Configs\n\n")
    out.write("| config | backend | init_eps | init_noise | d_scale_mode | formula | min_sep_unit | lambda_sep |\n")
    out.write("| --- | --- | ---: | ---: | --- | --- | ---: | ---: |\n")
    for r in rows:
        out.write(f"| `{r['config_name']}` | {r['backend']} | {r['init_eps']} | {r['init_noise_scale']} | {r['d_scale_mode']} | {r['dscale_effective_count_formula']} | {r['min_sep_unit']} | {r['lambda_sep']} |\n")
    out.write("\n## Main Results\n\n")
    out.write("| config | backend | top1 all | top1 cis | top1 trans | same/cross cis | cis Spearman | entropy | pU | min sep | mean sep | cos<0 | projection switch |\n")
    out.write("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
    for r in rows:
        out.write(f"| `{r['config_name']}` | {r['backend']} | {r['model_top1_accuracy_genome_all']} | {r['model_top1_accuracy_genome_cis']} | {r['model_top1_accuracy_genome_trans']} | {r['model_same_cross_accuracy_genome_cis']} | {r['mean_per_chrom_cis_distance_spearman']} | {r['final_mean_entropy']} | {r['final_mean_pU']} | {r['final_min_sep']} | {r['final_mean_sep']} | {r['copytrack_frac_cos_lt_0']} | {r['copytrack_frac_projection_sign_switch']} |\n")
    out.write("\n## d_scale Semantic Audit\n\n")
    out.write(f"`{headline}`\n\n")
    out.write("Interpretation rules: if posterior_count exact matches 003/004 while raw_count exact fails, the 011 regression is explained by d_scale semantic drift. If both fail, d_scale semantics alone is insufficient and inputs/init/backend/compiler/dirty-tree differences should be checked. This run is a regression audit, not a model improvement claim.\n")
PY

log_line "+ scripts/p9016_publish_light_result.sh $FULL_RUN_ROOT $LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "DSCALE_SEMANTICS_AUDIT_TSV=$FULL_RUN_ROOT/dscale_semantics_audit.tsv"
echo "DSCALE_SEMANTICS_HEADLINE=$HEADLINE"
