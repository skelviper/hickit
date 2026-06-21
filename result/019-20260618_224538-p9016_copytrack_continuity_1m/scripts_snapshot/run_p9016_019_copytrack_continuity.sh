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

RUN_ID="019-$(date +%Y%m%d_%H%M%S)-p9016_copytrack_continuity_1m"
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE3_ROOT/test_res}/$RUN_ID"
mkdir -p "$FULL_RUN_ROOT"/{logs,outputs,eval,scripts_snapshot}
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
OUTPUT_ROOT="$FULL_RUN_ROOT/outputs"
EVAL_ROOT="$FULL_RUN_ROOT/eval"
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
cp scripts/p9016_common_eval.sh \
   scripts/p9016_publish_light_result.sh \
   scripts/compute_p9016_copytrack_diag.py \
   scripts/summarize_p9016_model_sweep.py \
   scripts/smoke_p9016_copytrack.sh \
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

slug_float() {
	case "$1" in
	0) echo "0" ;;
	0.001) echo "0p001" ;;
	0.003) echo "0p003" ;;
	0.01) echo "0p01" ;;
	0.03) echo "0p03" ;;
	0.1) echo "0p1" ;;
	0.3) echo "0p3" ;;
	0.5) echo "0p5" ;;
	1|1.0) echo "1" ;;
	*) echo "$1" | tr '.' 'p' ;;
	esac
}

PAIR_PATH=${HK_BLIND_P9016_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
EVAL_PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-$PAIR_PATH}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-0.012}
MAX_PARALLEL=${HK_BLIND_019_MAX_PARALLEL:-${HK_BLIND_MAX_PARALLEL:-6}}
REQUESTED_BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-auto}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-0}
RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-1}
ALLOW_CUSTOM=${HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS:-0}

case "$MAX_PARALLEL" in
	''|*[!0-9]*)
		echo "error: HK_BLIND_019_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_019_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi

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
log_cmd scripts/smoke_p9016_copytrack.sh
if [ "$REQUESTED_BACKEND" = "gpu" ]; then
	log_cmd make gpu=1 run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
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
	echo "pairs	$PAIR_PATH"
	echo "eval_pairs	$EVAL_PAIRS"
	echo "eval_tdg	$EVAL_TDG"
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
	local config=$1
	local gamma=$2
	local min_sep=$3
	local lambda_sep=$4
	local lambda_copytrack=$5
	local out_dir="$OUTPUT_ROOT/$config"
	local eval_dir="$EVAL_ROOT/$config"
	local log_prefix="$FULL_RUN_ROOT/logs/$config"
	local train_cmd audit_cmd eval_cmd copytrack_cmd

	echo "CONFIG_START=$config"
	rm -rf "$out_dir" "$eval_dir"
	mkdir -p "$OUTPUT_ROOT" "$EVAL_ROOT" "$eval_dir"
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
		HK_BLIND_P9016_RELAX_BACKEND="$REQUESTED_BACKEND"
		HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split
		HK_BLIND_P9016_INIT_SCALE=10
		HK_BLIND_P9016_INIT_SEED=17
		HK_BLIND_P9016_INIT_EPS=0.5
		HK_BLIND_P9016_INIT_NOISE_SCALE=0.0
		HK_BLIND_P9016_MIN_SEP_UNIT="$min_sep"
		HK_BLIND_P9016_LAMBDA_SEP="$lambda_sep"
		HK_BLIND_P9016_LAMBDA_COPYTRACK="$lambda_copytrack"
		HK_BLIND_P9016_D_SCALE_MODE=posterior_count
		HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA="$gamma"
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6
		"$RUN_BIN"
	)
	log_command_line "${train_cmd[@]}"
	"${train_cmd[@]}" > "$log_prefix.train.log" 2>&1

	audit_cmd=(
		clean_env_prefix
		HK_BLIND_SAMPLE=P9016
		HK_BLIND_P9016_SAMPLE=P9016
		./audit_blind_p9016_full_cpu_output.bin "$out_dir"
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
			scripts/p9016_common_eval.sh "$out_dir" "$eval_dir"
		)
		log_command_line "${eval_cmd[@]}"
		"${eval_cmd[@]}" > "$log_prefix.eval_wrapper.log" 2>&1
		if [ "$REQUIRE_EVAL" = "1" ] && [ ! -s "$eval_dir/summary.tsv" ]; then
			echo "required eval missing for $config" >&2
			return 1
		fi
	else
		{
			echo "status	SKIPPED_EVAL_SMOKE_MODE"
			echo "sample	P9016"
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
	echo "CONFIG_DONE=$config"
}

launch_config() {
	local config=$1
	shift
	wait_one_if_needed
	run_config "$config" "$@" &
	CONFIG_PIDS+=("$!")
	CONFIG_LABELS+=("$config")
}

launch_config p9016_pcgamma1_copytrack0_sep_off_ieps0p5_noise0_seed17 1 0 0 0

COPYTRACK_LAMBDAS=${HK_BLIND_019_COPYTRACK_LAMBDAS-"0.001 0.003 0.01 0.03 0.1"}
COPYTRACK_LAMBDAS_GAMMA05=${HK_BLIND_019_COPYTRACK_LAMBDAS_GAMMA05-"0.003 0.01 0.03"}

for lambda_copytrack in $COPYTRACK_LAMBDAS; do
	cslug=$(slug_float "$lambda_copytrack")
	launch_config "p9016_pcgamma1_copytrack${cslug}_sep_off_ieps0p5_noise0_seed17" 1 0 0 "$lambda_copytrack"
done

for lambda_copytrack in $COPYTRACK_LAMBDAS_GAMMA05; do
	cslug=$(slug_float "$lambda_copytrack")
	launch_config "p9016_pcgamma0p5_copytrack${cslug}_sep_off_ieps0p5_noise0_seed17" 0.5 0 0 "$lambda_copytrack"
done

if [ "${HK_BLIND_019_RUN_SEP:-1}" = "1" ]; then
	for spec in \
		"1 1 0.5 0.01" \
		"1 1 1 0.01" \
		"1 1.5 1 0.01" \
		"1 1 0.5 0.03" \
		"1 1 1 0.03" \
		"0.5 1 0.5 0.01" \
		"0.5 1 1 0.01"
	do
		set -- $spec
		gamma=$1
		min_sep=$2
		lambda_sep=$3
		lambda_copytrack=$4
		gslug=$(slug_float "$gamma")
		mslug=$(slug_float "$min_sep")
		lslug=$(slug_float "$lambda_sep")
		cslug=$(slug_float "$lambda_copytrack")
		launch_config "p9016_pcgamma${gslug}_msep${mslug}_lsep${lslug}_copytrack${cslug}_ieps0p5_noise0_seed17" \
			"$gamma" "$min_sep" "$lambda_sep" "$lambda_copytrack"
	done
fi

wait_all_configs

log_line "+ scripts/summarize_p9016_model_sweep.py $FULL_RUN_ROOT > $FULL_RUN_ROOT/summary.tsv"
scripts/summarize_p9016_model_sweep.py "$FULL_RUN_ROOT" > "$FULL_RUN_ROOT/summary.tsv" 2> "$FULL_RUN_ROOT/logs/summarize.log"

python - <<'PY' "$FULL_RUN_ROOT/summary.tsv" "$FULL_RUN_ROOT/trans_delta.tsv"
import csv, math, sys
summary, out = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(summary), delimiter="\t"))
base = next((r for r in rows if r.get("config_name") == "p9016_pcgamma1_copytrack0_sep_off_ieps0p5_noise0_seed17"), None)
fields = [
    "config_name", "d_scale_posterior_gamma", "min_sep_unit", "lambda_sep", "lambda_copytrack",
    "top1_all", "top1_cis", "top1_trans", "same_cross_cis", "cis_spearman",
    "entropy", "pU", "final_min_sep", "sep_p05", "final_mean_sep",
    "copytrack_frac_cos_lt_0", "copytrack_force_l1_over_contact_force_l1",
    "delta_top1_trans_vs_baseline", "delta_top1_all_vs_baseline", "target_plus_0p1_met"
]
def f(row, key):
    try:
        x = float(row.get(key, "NA"))
    except ValueError:
        return None
    return x if math.isfinite(x) else None
def fmt(x):
    return "NA" if x is None else f"{x:.9g}"
with open(out, "w", newline="") as g:
    w = csv.DictWriter(g, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    base_trans = f(base, "model_top1_accuracy_genome_trans") if base else None
    base_all = f(base, "model_top1_accuracy_genome_all") if base else None
    for r in rows:
        trans = f(r, "model_top1_accuracy_genome_trans")
        allv = f(r, "model_top1_accuracy_genome_all")
        dtrans = None if trans is None or base_trans is None else trans - base_trans
        dall = None if allv is None or base_all is None else allv - base_all
        w.writerow({
            "config_name": r.get("config_name", "NA"),
            "d_scale_posterior_gamma": r.get("d_scale_posterior_gamma", "NA"),
            "min_sep_unit": r.get("min_sep_unit", "NA"),
            "lambda_sep": r.get("lambda_sep", "NA"),
            "lambda_copytrack": r.get("lambda_copytrack", "NA"),
            "top1_all": r.get("model_top1_accuracy_genome_all", "NA"),
            "top1_cis": r.get("model_top1_accuracy_genome_cis", "NA"),
            "top1_trans": r.get("model_top1_accuracy_genome_trans", "NA"),
            "same_cross_cis": r.get("model_same_cross_accuracy_genome_cis", "NA"),
            "cis_spearman": r.get("mean_per_chrom_cis_distance_spearman", "NA"),
            "entropy": r.get("final_mean_entropy", "NA"),
            "pU": r.get("final_mean_pU", "NA"),
            "final_min_sep": r.get("final_min_sep", "NA"),
            "sep_p05": r.get("sep_p05", "NA"),
            "final_mean_sep": r.get("final_mean_sep", "NA"),
            "copytrack_frac_cos_lt_0": r.get("copytrack_frac_cos_lt_0", "NA"),
            "copytrack_force_l1_over_contact_force_l1": r.get("copytrack_force_l1_over_contact_force_l1", "NA"),
            "delta_top1_trans_vs_baseline": fmt(dtrans),
            "delta_top1_all_vs_baseline": fmt(dall),
            "target_plus_0p1_met": "1" if dtrans is not None and dtrans >= 0.1 else "0",
        })
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
summary = list(csv.DictReader(open(root / "summary.tsv"), delimiter="\t"))
delta = list(csv.DictReader(open(root / "trans_delta.tsv"), delimiter="\t"))
def f(row, key):
    try:
        x = float(row.get(key, "NA"))
    except ValueError:
        return None
    return x if math.isfinite(x) else None
valid = [r for r in delta if f(r, "top1_trans") is not None]
best = max(valid, key=lambda r: f(r, "top1_trans")) if valid else None
target = [r for r in delta if r.get("target_plus_0p1_met") == "1"]
with open(root / "README.md", "w") as out:
    out.write(f"# {root.name}\n\n")
    out.write("Controlled blind P9016 copy-track continuity experiment at 1 Mb. The experiment tests whether a homolog-vector continuity force can reduce local copy-track flipping and improve trans phase accuracy.\n\n")
    out.write("## Paths\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- lightweight hickit result: `{light}`\n")
    out.write(f"- summary.tsv: `{root / 'summary.tsv'}`\n")
    out.write(f"- trans_delta.tsv: `{root / 'trans_delta.tsv'}`\n")
    out.write(f"- commands.log: `{root / 'commands.log'}`\n\n")
    out.write("## Training Boundary\n\n")
    out.write("Training uses only P9016 pairs. SNP/phase labels and CHARM/3DG are used only in eval through `scripts/p9016_common_eval.sh`. The copy-track force is blind geometry: adjacent same-chromosome homolog vectors are encouraged to be continuous.\n\n")
    out.write("## Git And Build\n\n")
    out.write(f"- git commit: `{git_commit}`\n")
    out.write(f"- git dirty count at start: `{dirty}`\n")
    if dirty != "0":
        out.write("- WARNING: dirty tree at launch; inspect `logs/git_status.txt` and `logs/git_diff.patch`.\n")
    out.write(f"- run binary sha256: `{run_hash}`\n")
    out.write(f"- backend: `{backend}`\n\n")
    out.write("## Main Results\n\n")
    out.write("| config | gamma | min_sep | lambda_sep | lambda_copytrack | top1 all | top1 cis | top1 trans | delta trans | cis Spearman | entropy | pU | sep_p05 | cos<0 |\n")
    out.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
    for r in delta:
        out.write(f"| `{r['config_name']}` | {r['d_scale_posterior_gamma']} | {r['min_sep_unit']} | {r['lambda_sep']} | {r['lambda_copytrack']} | {r['top1_all']} | {r['top1_cis']} | {r['top1_trans']} | {r['delta_top1_trans_vs_baseline']} | {r['cis_spearman']} | {r['entropy']} | {r['pU']} | {r['sep_p05']} | {r['copytrack_frac_cos_lt_0']} |\n")
    out.write("\n## Interpretation\n\n")
    if best:
        out.write(f"- Best observed trans config: `{best['config_name']}` with trans top1 `{best['top1_trans']}` and delta `{best['delta_top1_trans_vs_baseline']}` versus baseline.\n")
    if target:
        out.write("- The +0.1 trans target was met by at least one config in this sweep.\n")
    else:
        out.write("- The +0.1 trans target was not met in this sweep; continue with chromosome-level gauge synchronization or E-step prior diagnostics.\n")
    out.write("- Accept copy-track continuity only if trans improves without a large cis/Spearman/entropy penalty. A drop in copytrack cos<0 without trans improvement means internal chromosome copy tracks stabilized but cross-chromosome gauge remains unresolved.\n")
PY

log_line "+ scripts/p9016_publish_light_result.sh $FULL_RUN_ROOT $LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "TRANS_DELTA_TSV=$FULL_RUN_ROOT/trans_delta.tsv"
