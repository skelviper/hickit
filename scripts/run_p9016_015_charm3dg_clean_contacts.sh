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

RUN_ID="015-$(date +%Y%m%d_%H%M%S)-p9016_charm3dg_clean_contacts_1m"
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE_ROOT/test_res}/$RUN_ID"
mkdir -p "$FULL_RUN_ROOT"/{logs,outputs,eval,work_outputs,synthetic_pairs,scripts_snapshot}
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
   scripts/plot_p9016_contact_matrix_qc.py \
   scripts/summarize_p9016_015.py \
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

PAIRS=${HK_BLIND_P9016_OBSERVED_PAIRS:-${HK_BLIND_P9016_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}}
EVAL_PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-$PAIRS}
TDG20=${HK_BLIND_P9016_20K_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.20k.3dg.gz}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
TDG_RES=${HK_BLIND_015_TDG_RESOLUTION_BP:-20000}
RADII=${HK_BLIND_015_RADII_PR:-"1 2"}
TARGET_MODE=${HK_BLIND_015_TARGET_MODE:-observed_non_same_bin}
TARGET_COUNT=${HK_BLIND_015_TARGET_COUNT:-0}
CANDIDATE_MULTIPLIER=${HK_BLIND_015_CANDIDATE_MULTIPLIER:-4}
PARTICLE_RADIUS=${HK_BLIND_015_PARTICLE_RADIUS:-0}
MIN_GENOMIC_SEP=${HK_BLIND_015_MIN_GENOMIC_SEPARATION:-0}
SEED=${HK_BLIND_015_SEED:-17}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-0.012}
MAX_PARALLEL=${HK_BLIND_015_MAX_PARALLEL:-2}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-1}
NEAREST_IF_INSUFFICIENT=${HK_BLIND_015_NEAREST_IF_INSUFFICIENT:-0}
REQUESTED_BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-auto}
if [ "$REQUESTED_BACKEND" = "auto" ]; then
	if backend_available gpu; then
		REQUESTED_BACKEND=gpu
	else
		REQUESTED_BACKEND=cpu
	fi
fi

case "$MAX_PARALLEL" in
	''|*[!0-9]*)
		echo "error: HK_BLIND_015_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_015_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi

for required in "$PAIRS" "$EVAL_PAIRS" "$TDG20" "$EVAL_TDG"; do
	if [ ! -s "$required" ]; then
		echo "error: required input missing: $required" >&2
		exit 1
	fi
done

log_cmd make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
log_line "+ HK_BLIND_TEST_RES_ROOT=$FULL_RUN_ROOT make smoke_blind_p9016_minimal"
HK_BLIND_TEST_RES_ROOT="$FULL_RUN_ROOT" make smoke_blind_p9016_minimal
rm -rf "$FULL_RUN_ROOT/smoke"
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
	echo "target_mode	$TARGET_MODE"
	echo "backend	$REQUESTED_BACKEND"
	echo "bin_size_bp	$BIN_SIZE"
	echo "n_iter	$N_ITER"
	echo "relax_steps	$RELAX_STEPS"
	echo "relax_step	$RELAX_STEP"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

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

generate_pairs_for_radius() {
	local radius=$1
	local config=$2
	local out_pairs="$FULL_RUN_ROOT/synthetic_pairs/${config}.pairs.gz"
	local copy_state="$FULL_RUN_ROOT/synthetic_pairs/${config}.copy_state.tsv.gz"
	local metadata="$FULL_RUN_ROOT/synthetic_pairs/${config}.metadata.tsv"
	local metadata_json="$FULL_RUN_ROOT/synthetic_pairs/${config}.json"
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
	"${cmd[@]}" > "$FULL_RUN_ROOT/logs/${config}.generate.log" 2>&1
}

run_config() {
	local radius=$1
	local config=$2
	local synthetic_pairs="$FULL_RUN_ROOT/synthetic_pairs/${config}.pairs.gz"
	local sample_out_root="$FULL_RUN_ROOT/work_outputs/$config"
	local work_out_dir="$sample_out_root/$config"
	local out_dir="$FULL_RUN_ROOT/outputs/$config"
	local eval_dir="$FULL_RUN_ROOT/eval/$config"
	local log_prefix="$FULL_RUN_ROOT/logs/$config"

	echo "CONFIG_START=$config"
	rm -rf "$sample_out_root" "$out_dir" "$eval_dir"
	mkdir -p "$sample_out_root" "$eval_dir"
	local train_cmd=(
		clean_env_prefix
		HK_BLIND_SAMPLE=P9016
		HK_BLIND_P9016_SAMPLE=P9016
		HK_BLIND_GIT_COMMIT="$GIT_COMMIT"
		HK_BLIND_GIT_DIRTY_COUNT="$GIT_DIRTY_COUNT"
		HK_BLIND_BINARY_HASH="$RUN_HASH"
		HK_BLIND_P9016_PAIRS="$synthetic_pairs"
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
		HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split
		HK_BLIND_P9016_INIT_SCALE=10
		HK_BLIND_P9016_INIT_SEED=17
		HK_BLIND_P9016_INIT_EPS=0.5
		HK_BLIND_P9016_INIT_NOISE_SCALE=0.0
		HK_BLIND_P9016_MIN_SEP_UNIT=0
		HK_BLIND_P9016_LAMBDA_SEP=0
		HK_BLIND_P9016_D_SCALE_MODE=posterior_count
		HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA=1
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6
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
	local radius=$1
	local config=$2
	wait_one_if_needed
	run_config "$radius" "$config" &
	CONFIG_PIDS+=("$!")
	CONFIG_LABELS+=("$config")
}

for radius in $RADII; do
	slug=$(radius_slug "$radius")
	config="p9016_charm3dg20k_${slug}_clean_contacts_pcgamma1_ieps0p5_noise0_seed17"
	generate_pairs_for_radius "$radius" "$config"
	launch_config "$radius" "$config"
done

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
	config="p9016_charm3dg20k_${slug}_clean_contacts_pcgamma1_ieps0p5_noise0_seed17"
	matrix_qc_cmd+=(--matrix "CHARM20k <=${radius} PR synthetic=$FULL_RUN_ROOT/synthetic_pairs/${config}.pairs.gz")
	matrix_qc_cmd+=(--synthetic-copy-matrix "CHARM20k <=${radius} PR synthetic=$FULL_RUN_ROOT/synthetic_pairs/${config}.copy_state.tsv.gz")
done
log_command_line "${matrix_qc_cmd[@]}"
"${matrix_qc_cmd[@]}" > "$FULL_RUN_ROOT/logs/contact_matrix_qc.log" 2>&1

log_line "+ scripts/summarize_p9016_015.py $FULL_RUN_ROOT > $FULL_RUN_ROOT/summary.tsv"
scripts/summarize_p9016_015.py "$FULL_RUN_ROOT" > "$FULL_RUN_ROOT/summary.tsv" 2> "$FULL_RUN_ROOT/logs/summarize.log"
cp "$FULL_RUN_ROOT/summary.tsv" "$FULL_RUN_ROOT/synthetic_contact_summary.tsv"

python - <<'PY' "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" "$GIT_COMMIT" "$GIT_DIRTY_COUNT" "$RUN_HASH" "$REQUESTED_BACKEND"
import csv
import sys
from pathlib import Path
root = Path(sys.argv[1])
light, git_commit, dirty, run_hash, backend = sys.argv[2:7]
rows = list(csv.DictReader(open(root / "summary.tsv"), delimiter="\t"))
with open(root / "README.md", "w") as out:
    out.write(f"# {root.name}\n\n")
    out.write("Positive-control experiment: contacts were generated from P9016 CHARM 20kb 3DG geometry at radius thresholds, downsampled to the observed contact level, then used as unphased input for the existing blind softall runner.\n\n")
    out.write("## Paths\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- summary.tsv: `{root / 'summary.tsv'}`\n")
    out.write(f"- synthetic pairs: `{root / 'synthetic_pairs'}`\n")
    out.write(f"- commands.log: `{root / 'commands.log'}`\n\n")
    out.write("## Boundary\n\n")
    out.write("This is not a blind baseline. Training contacts are reference-derived from CHARM/3DG, but no SNP phase labels are written into the synthetic pairs. Manifests must report `input_contact_source=charm3dg_derived_pairs`, `uses_phase_labels=0`, and `uses_charm_for_training=1`.\n\n")
    out.write("## Git And Build\n\n")
    out.write(f"- git commit at launch: `{git_commit}`\n")
    out.write(f"- git dirty count at launch: `{dirty}`\n")
    if dirty != "0":
        out.write("- WARNING: dirty tree at launch; inspect `logs/git_status.txt` and `logs/git_diff.patch` in the full result root.\n")
    out.write(f"- run binary sha256: `{run_hash}`\n")
    out.write(f"- backend: `{backend}`\n\n")
    out.write("## Results\n\n")
    out.write("| config | radius_pr | contacts | trans contacts | top1 all | top1 cis | top1 trans | same/cross trans | cis Spearman | entropy | pU | min sep | mean sep | cos<0 |\n")
    out.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
    for r in rows:
        out.write(f"| `{r['config_name']}` | {r['radius_pr']} | {r['output_pairs_total']} | {r['output_pairs_trans']} | {r['model_top1_accuracy_genome_all']} | {r['model_top1_accuracy_genome_cis']} | {r['model_top1_accuracy_genome_trans']} | {r['model_same_cross_accuracy_genome_trans']} | {r['mean_per_chrom_cis_distance_spearman']} | {r['final_mean_entropy']} | {r['final_mean_pU']} | {r['final_min_sep']} | {r['final_mean_sep']} | {r['copytrack_frac_cos_lt_0']} |\n")
    out.write("\n## Interpretation Rules\n\n")
    out.write("- If trans improves strongly here, the raw-pairs failure is consistent with noise/sparsity or misleading contacts being a major bottleneck.\n")
    out.write("- If trans still fails with CHARM-derived clean contacts, the bottleneck is likely in gauge/copy-track identifiability or EM/FDG dynamics rather than contact noise alone.\n")
    out.write("- Do not treat this as a production model improvement because training used CHARM/3DG-derived contacts.\n")
PY

log_line "+ scripts/p9016_publish_light_result.sh $FULL_RUN_ROOT $LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "SYNTHETIC_CONTACT_SUMMARY_TSV=$FULL_RUN_ROOT/synthetic_contact_summary.tsv"
