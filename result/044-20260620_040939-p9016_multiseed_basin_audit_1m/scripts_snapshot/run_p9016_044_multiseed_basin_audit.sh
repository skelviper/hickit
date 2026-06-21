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

RUN_ID=${HK_BLIND_044_RUN_ID:-"044-$(date +%Y%m%d_%H%M%S)-p9016_multiseed_basin_audit_1m"}
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
MAX_PARALLEL=${HK_BLIND_044_MAX_PARALLEL:-${HK_BLIND_MAX_PARALLEL:-10}}
RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-1}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-1}
CONFIG_SET=${HK_BLIND_044_CONFIG_SET:-full}
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
		echo "error: HK_BLIND_044_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_044_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi
if [ "$BIN_SIZE" != "1000000" ]; then
	echo "error: 044 is fixed at 1Mb; HK_BLIND_P9016_BIN_SIZE_BP=$BIN_SIZE" >&2
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
		echo "error: formal 044 requires approved P9016 training pairs: $APPROVED_P9016_PAIRS" >&2
		echo "got: $PAIR_PATH" >&2
		exit 1
	fi
	if ! same_realpath "$EVAL_PAIRS" "$APPROVED_P9016_PAIRS"; then
		echo "error: formal 044 requires approved P9016 eval pairs: $APPROVED_P9016_PAIRS" >&2
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
	echo "experiment_control	multiseed basin audit; no new model mechanism"
	echo "success_criterion	full-denominator trans top1 accuracy +0.1 over td1 seed17 baseline"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

run_config() {
	local config=$1
	local seed=$2
	local init_noise=$3
	local trans_dscale=$4
	local lambda_copytrack=$5
	local lambda_global=$6
	local mstep_lambda=$7
	local mstep_power=$8
	local mstep_warmup=$9
	local out_root="$OUTPUT_ROOT/$config"
	local eval_dir="$EVAL_ROOT/$config"
	local log_prefix="$LOG_ROOT/$config"
	local train_cmd audit_cmd eval_cmd copytrack_cmd

	echo "CONFIG_START=$config seed=$seed noise=$init_noise trans_dscale=$trans_dscale lambda_copytrack=$lambda_copytrack lambda_global=$lambda_global mstep_lambda=$mstep_lambda mstep_power=$mstep_power mstep_warmup=$mstep_warmup"
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
		HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split
		HK_BLIND_P9016_INIT_SEED="$seed"
		HK_BLIND_P9016_INIT_EPS=0.5
		HK_BLIND_P9016_INIT_NOISE_SCALE="$init_noise"
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
		HK_BLIND_P9016_TRANS_GATE_MODE=off
		HK_BLIND_P9016_TRANS_CHR_PAIR_PRIOR_LAMBDA=0
		HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_LAMBDA="$mstep_lambda"
		HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_EPS=0.001
		HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_POWER="$mstep_power"
		HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_WARMUP_ITER="$mstep_warmup"
		HK_BLIND_WRITE_RAW=0
		"$RUN_BIN"
	)
	log_command_line "${train_cmd[@]}"
	"${train_cmd[@]}" > "$log_prefix.train.log" 2>&1

	audit_cmd=(
		clean_env_prefix
		HK_BLIND_AUDIT_ALLOW_CUSTOM_RAW_PAIRS="$ALLOW_CUSTOM"
		"$AUDIT_BIN" "$out_root/$config"
	)
	log_command_line "${audit_cmd[@]}"
	"${audit_cmd[@]}" > "$log_prefix.audit.log" 2>&1

	manifest="$out_root/$config/p9016_full.manifest.tsv"
	grep -q $'input_contact_source\traw_pairs' "$manifest"
	grep -q $'uses_phase_labels\t0' "$manifest"
	grep -q $'uses_charm_or_reference\t0' "$manifest"
	grep -q $'copy_labels_are_gauge_only\t1' "$manifest"
	if [ "$CONFIG_SET" != "smoke" ]; then
		grep -q $'approved_p9016_raw_pairs_realpath\t1' "$manifest"
	fi

	if [ "$RUN_EVAL" = "1" ]; then
		eval_cmd=(
			clean_env_prefix
			HK_BLIND_SAMPLE=P9016
			HK_BLIND_P9016_SAMPLE=P9016
			HK_BLIND_P9016_EVAL_PAIRS="$EVAL_PAIRS"
			HK_BLIND_P9016_EVAL_TDG="$EVAL_TDG"
			HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE"
			HK_BLIND_REQUIRE_EVAL="$REQUIRE_EVAL"
			scripts/p9016_common_eval.sh "$out_root/$config" "$eval_dir"
		)
		log_command_line "${eval_cmd[@]}"
		"${eval_cmd[@]}" > "$log_prefix.eval.log" 2>&1
	else
		printf 'status\tSKIPPED_EVAL_DISABLED\n' > "$eval_dir/eval_status.tsv"
	fi

	if [ -s "$out_root/$config/p9016_full.coords.tsv" ]; then
		copytrack_cmd=(
			python scripts/compute_p9016_copytrack_diag.py
			--config-output-dir "$out_root/$config"
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
		"p9016_multiseed_smoke_best035_seed17|17|0|0.5|0.03|0.003|0|1|0"
		"p9016_multiseed_smoke_best035_seed23_noise0p05|23|0.05|0.5|0.03|0.003|0|1|0"
	)
else
	SEEDS=(${HK_BLIND_044_SEEDS:-17 23 31 47 59 71 83 97})
	CONFIGS=(
		"p9016_multiseed_td1_seed17_noise0|17|0|1|0|0|0|1|0"
		"p9016_multiseed_best035_seed17_noise0|17|0|0.5|0.03|0.003|0|1|0"
	)
	for seed in "${SEEDS[@]}"; do
		CONFIGS+=("p9016_multiseed_best035_seed${seed}_noise0p05|$seed|0.05|0.5|0.03|0.003|0|1|0")
		CONFIGS+=("p9016_multiseed_chrpair_m0p5_p4_w10_seed${seed}_noise0p05|$seed|0.05|0.5|0|0|0.5|4|10")
		CONFIGS+=("p9016_multiseed_chrpair_m1_p1_w0_seed${seed}_noise0p05|$seed|0.05|0.5|0|0|1|1|0")
	done
fi

failures=0
running=0
for spec in "${CONFIGS[@]}"; do
	IFS='|' read -r config seed init_noise trans_dscale lambda_copytrack lambda_global mstep_lambda mstep_power mstep_warmup <<< "$spec"
	run_config "$config" "$seed" "$init_noise" "$trans_dscale" "$lambda_copytrack" "$lambda_global" "$mstep_lambda" "$mstep_power" "$mstep_warmup" &
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

log_command_line python - "$FULL_RUN_ROOT" "[inline multiseed_basin_delta_summary]"
python - "$FULL_RUN_ROOT" <<'PY'
import csv
import re
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = [
    row for row in csv.DictReader((root / "summary.tsv").open(), delimiter="\t")
    if row.get("model_top1_accuracy_genome_trans", "NA") not in ("", "NA")
]

def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None

def family(name):
    if "_td1_" in name:
        return "td1"
    if "_best035_" in name:
        return "best035"
    if "_chrpair_m0p5_p4_w10_" in name:
        return "chrpair_m0p5_p4_w10"
    if "_chrpair_m1_p1_w0_" in name:
        return "chrpair_m1_p1_w0"
    return "other"

def seed_from(name):
    m = re.search(r"_seed([0-9]+)_", name)
    return m.group(1) if m else "NA"

def noise_from(name):
    if "noise0p05" in name:
        return "0.05"
    if "noise0" in name:
        return "0"
    return "NA"

base = next((r for r in rows if r.get("config_name") == "p9016_multiseed_td1_seed17_noise0"), None)
anchor = next((r for r in rows if r.get("config_name") == "p9016_multiseed_best035_seed17_noise0"), None)
base_t = fnum(base.get("model_top1_accuracy_genome_trans")) if base else None
anchor_t = fnum(anchor.get("model_top1_accuracy_genome_trans")) if anchor else None
for row in rows:
    name = row.get("config_name", "")
    row["basin_family"] = family(name)
    row["seed"] = seed_from(name)
    row["init_noise_short"] = noise_from(name)
    t = fnum(row.get("model_top1_accuracy_genome_trans"))
    if t is not None and base_t is not None:
        row["delta_trans_vs_td1_seed17"] = f"{t - base_t:.9g}"
        row["target_plus_0p1_vs_td1_met"] = "1" if t >= base_t + 0.1 else "0"
    else:
        row["delta_trans_vs_td1_seed17"] = "NA"
        row["target_plus_0p1_vs_td1_met"] = "0"
    if t is not None and anchor_t is not None:
        row["delta_trans_vs_best035_seed17"] = f"{t - anchor_t:.9g}"
    else:
        row["delta_trans_vs_best035_seed17"] = "NA"

fields = [
    "config_name",
    "basin_family",
    "seed",
    "init_noise_short",
    "status",
    "backend",
    "init_seed",
    "init_noise_scale",
    "trans_dscale_multiplier",
    "lambda_copytrack",
    "lambda_global_copytrack",
    "trans_chr_pair_mstep_lambda",
    "trans_chr_pair_mstep_power",
    "trans_chr_pair_mstep_warmup_iter",
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
    "delta_trans_vs_td1_seed17",
    "delta_trans_vs_best035_seed17",
    "target_plus_0p1_vs_td1_met",
]
rows_sorted = sorted(rows, key=lambda r: fnum(r.get("model_top1_accuracy_genome_trans")) or -1, reverse=True)
with (root / "multiseed_basin_delta_summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows_sorted:
        writer.writerow({key: row.get(key, "NA") for key in fields})

headline = "PLUS_0P1_MET" if any(row.get("target_plus_0p1_vs_td1_met") == "1" for row in rows) else "NO_PLUS_0P1"
(root / "headline.txt").write_text(headline + "\n")
PY

log_command_line python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" "[inline README]"
python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
light = Path(sys.argv[2])
rows = list(csv.DictReader((root / "multiseed_basin_delta_summary.tsv").open(), delimiter="\t"))
headline = (root / "headline.txt").read_text().strip()
fields = [
    "config_name",
    "basin_family",
    "seed",
    "init_noise_short",
    "trans_chr_pair_mstep_lambda",
    "trans_chr_pair_mstep_power",
    "trans_chr_pair_mstep_warmup_iter",
    "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans",
    "delta_trans_vs_td1_seed17",
    "delta_trans_vs_best035_seed17",
    "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "final_mean_entropy",
    "final_mean_pU",
    "final_min_sep",
    "sep_p05",
    "copytrack_frac_cos_lt_0",
    "target_plus_0p1_vs_td1_met",
]
with (root / "README.md").open("w") as out:
    out.write("# 044 P9016 Multiseed Basin Audit 1Mb\n\n")
    out.write("This controlled blind-training experiment tests whether the remaining trans gap is merely a seed/basin issue. It does not add a new model mechanism. Training uses only the approved P9016 raw pairs file; phase labels and CHARM/3DG are eval-only.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- headline: `{headline}`\n")
    out.write("- success criterion: full-denominator trans top1 accuracy must improve by at least 0.1 over the td1 seed17/noise0 baseline\n")
    out.write("- tested families: td1 anchor, best035 anchor/multiseed, and two chr-pair M-step semantics whose 038 eval-only pair oracle had +0.1 ceiling\n")
    out.write("- copy0/copy1 are gauge labels; no maternal/paternal labels are used in training\n\n")
    out.write("## Main Results\n\n")
    out.write("| " + " | ".join(fields) + " |\n")
    out.write("|" + "|".join(["---"] * len(fields)) + "|\n")
    for row in rows:
        out.write("| " + " | ".join(row.get(key, "NA") for key in fields) + " |\n")
    out.write("\n## Interpretation Boundary\n\n")
    out.write("If no seed reaches the +0.1 full-denominator trans target, this argues against a simple basin-selection explanation under the current raw-pairs-only 1Mb training boundary. A positive result here would still need a blind selection rule, because eval labels are used only after training.\n")
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
if [ -s "$FULL_RUN_ROOT/multiseed_basin_delta_summary.tsv" ]; then
	echo "MULTISEED_BASIN_DELTA_SUMMARY_TSV=$FULL_RUN_ROOT/multiseed_basin_delta_summary.tsv"
fi
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt" 2>/dev/null || echo SMOKE_RESULT_REMOVED)"
