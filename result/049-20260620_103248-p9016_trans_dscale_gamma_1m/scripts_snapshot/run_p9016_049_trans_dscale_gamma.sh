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

RUN_ID=${HK_BLIND_049_RUN_ID:-"049-$(date +%Y%m%d_%H%M%S)-p9016_trans_dscale_gamma_1m"}
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE_ROOT/test_res}/$RUN_ID"
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
MAX_PARALLEL=${HK_BLIND_049_MAX_PARALLEL:-${HK_BLIND_MAX_PARALLEL:-10}}
RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-1}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-1}
CONFIG_SET=${HK_BLIND_049_CONFIG_SET:-full}
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
		echo "error: HK_BLIND_049_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_049_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi
if [ "$BIN_SIZE" != "1000000" ]; then
	echo "error: 049 is fixed at 1Mb; HK_BLIND_P9016_BIN_SIZE_BP=$BIN_SIZE" >&2
	exit 2
fi
if [ "$CONFIG_SET" != "smoke" ]; then
	if [ "$N_ITER" != "100" ] || [ "$RELAX_STEPS" != "100" ]; then
		echo "error: formal 049 requires n_iter=100 and relax_steps=100" >&2
		exit 2
	fi
	if [ "$RUN_EVAL" != "1" ] || [ "$REQUIRE_EVAL" != "1" ]; then
		echo "error: formal 049 requires eval enabled and required" >&2
		exit 2
	fi
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
		echo "error: formal 049 requires approved P9016 training pairs: $APPROVED_P9016_PAIRS" >&2
		echo "got: $PAIR_PATH" >&2
		exit 1
	fi
	if ! same_realpath "$EVAL_PAIRS" "$APPROVED_P9016_PAIRS"; then
		echo "error: formal 049 requires approved P9016 eval pairs: $APPROVED_P9016_PAIRS" >&2
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
	echo "experiment_control	trans-only posterior dscale gamma and trans dscale multiplier; no posterior/prior/rho/temperature changes"
	echo "success_criterion	full-denominator trans top1 accuracy +0.1 over current best blind-safe standard trans"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

run_config() {
	local config=$1
	local seed=$2
	local init_noise=$3
	local trans_gamma=$4
	local trans_dscale=$5
	local out_root="$OUTPUT_ROOT/$config"
	local eval_dir="$EVAL_ROOT/$config"
	local log_prefix="$LOG_ROOT/$config"
	local manifest
	local train_cmd audit_cmd eval_cmd copytrack_cmd

	echo "CONFIG_START=$config seed=$seed noise=$init_noise trans_gamma=$trans_gamma trans_dscale=$trans_dscale"
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
		HK_BLIND_P9016_TRANS_D_SCALE_POSTERIOR_GAMMA="$trans_gamma"
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6
		HK_BLIND_P9016_RHO_TRAIN_MODE=constant
		HK_BLIND_P9016_RHO_TRAIN_FLOOR=0
		HK_BLIND_P9016_TEMPERATURE_START=1
		HK_BLIND_P9016_TEMPERATURE_END=1
		HK_BLIND_P9016_TRANS_DSCALE_MULTIPLIER="$trans_dscale"
		HK_BLIND_P9016_TRANS_K_MULTIPLIER=1
		HK_BLIND_P9016_LAMBDA_COPYTRACK=0.03
		HK_BLIND_P9016_LAMBDA_GLOBAL_COPYTRACK=0.003
		HK_BLIND_P9016_LAMBDA_NORMDIR_COPYTRACK=0
		HK_BLIND_P9016_INIT_COORD_ANCHOR_K=0
		HK_BLIND_P9016_MIN_SEP_UNIT=0
		HK_BLIND_P9016_LAMBDA_SEP=0
		HK_BLIND_P9016_TRANS_TOP1_MODE=off
		HK_BLIND_P9016_TRANS_GATE_MODE=off
		HK_BLIND_P9016_TRANS_CHR_PAIR_PRIOR_LAMBDA=0
		HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_LAMBDA=0
		HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_MODE=state4
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
	grep -q $'dscale_mode\tposterior_count' "$manifest"
	grep -q "trans_d_scale_posterior_gamma	$trans_gamma" "$manifest"
	grep -q $'trans_contact_scaling_uses_phase_labels\t0' "$manifest"
	grep -q $'trans_contact_scaling_uses_charm_or_reference\t0' "$manifest"
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
		"p9016_transgamma_smoke_g1_td1_seed17|17|0|1|1"
		"p9016_transgamma_smoke_g0p5_td0p75_seed17|17|0|0.5|0.75"
	)
else
	CONFIGS=(
		"p9016_transgamma_g0_td0p75_seed17_best035|17|0|0|0.75"
		"p9016_transgamma_g0p25_td0p75_seed17_best035|17|0|0.25|0.75"
		"p9016_transgamma_g0p5_td0p75_seed17_best035|17|0|0.5|0.75"
		"p9016_transgamma_g0p75_td0p75_seed17_best035|17|0|0.75|0.75"
		"p9016_transgamma_g1_td0p75_seed17_best035|17|0|1|0.75"
		"p9016_transgamma_g0_td1_seed17_best035|17|0|0|1"
		"p9016_transgamma_g0p25_td1_seed17_best035|17|0|0.25|1"
		"p9016_transgamma_g0p5_td1_seed17_best035|17|0|0.5|1"
		"p9016_transgamma_g0p75_td1_seed17_best035|17|0|0.75|1"
		"p9016_transgamma_g1_td1_seed17_best035|17|0|1|1"
		"p9016_transgamma_g1_td0p5_seed17_best035_replay|17|0|1|0.5"
		"p9016_transgamma_g0_td0p75_seed71_best035|71|0.05|0|0.75"
		"p9016_transgamma_g0p25_td0p75_seed71_best035|71|0.05|0.25|0.75"
		"p9016_transgamma_g0p5_td0p75_seed71_best035|71|0.05|0.5|0.75"
		"p9016_transgamma_g0p75_td0p75_seed71_best035|71|0.05|0.75|0.75"
		"p9016_transgamma_g1_td0p75_seed71_best035|71|0.05|1|0.75"
		"p9016_transgamma_g0_td1_seed71_best035|71|0.05|0|1"
		"p9016_transgamma_g0p25_td1_seed71_best035|71|0.05|0.25|1"
		"p9016_transgamma_g0p5_td1_seed71_best035|71|0.05|0.5|1"
		"p9016_transgamma_g0p75_td1_seed71_best035|71|0.05|0.75|1"
		"p9016_transgamma_g1_td1_seed71_best035|71|0.05|1|1"
		"p9016_transgamma_g1_td0p5_seed71_best035_replay|71|0.05|1|0.5"
	)
fi

failures=0
running=0
for spec in "${CONFIGS[@]}"; do
	IFS='|' read -r config seed init_noise trans_gamma trans_dscale <<< "$spec"
	run_config "$config" "$seed" "$init_noise" "$trans_gamma" "$trans_dscale" &
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

log_command_line python - "$FULL_RUN_ROOT" "$CONFIG_SET" "${CONFIGS[@]}" "[inline expected matrix audit]"
python - "$FULL_RUN_ROOT" "$CONFIG_SET" "${CONFIGS[@]}" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
config_set = sys.argv[2]
expected = [spec.split("|", 1)[0] for spec in sys.argv[3:]]
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))
seen = [row.get("config_name", "") for row in rows]
allowed_status = {"OK"}
if config_set == "smoke":
    allowed_status.add("SKIPPED_EVAL_DISABLED")
missing = sorted(set(expected) - set(seen))
extra = sorted(set(seen) - set(expected))
bad_status = sorted(row.get("config_name", "") for row in rows if row.get("status") not in allowed_status)
missing_eval = sorted(
    row.get("config_name", "")
    for row in rows
    if config_set != "smoke" and row.get("model_top1_accuracy_genome_trans", "NA") in ("", "NA")
)
with (root / "matrix_completeness_audit.tsv").open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["field", "value"])
    writer.writerow(["expected_n", len(expected)])
    writer.writerow(["observed_n", len(rows)])
    writer.writerow(["missing", ",".join(missing) or "NA"])
    writer.writerow(["extra", ",".join(extra) or "NA"])
    writer.writerow(["bad_status", ",".join(bad_status) or "NA"])
    writer.writerow(["missing_eval", ",".join(missing_eval) or "NA"])
if missing or extra or bad_status or missing_eval or len(rows) != len(expected):
    raise SystemExit("incomplete 049 matrix")
PY

log_command_line python - "$FULL_RUN_ROOT" "[inline trans_dscale_gamma_delta_summary]"
python - "$FULL_RUN_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = [
    row for row in csv.DictReader((root / "summary.tsv").open(), delimiter="\t")
    if row.get("model_top1_accuracy_genome_trans", "NA") not in ("", "NA")
]
CURRENT_BEST_TRANS = 0.358585

def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None

def seed_family(row):
    return "seed71_noise0p05" if row.get("init_seed") == "71" else "seed17_noise0"

base = {}
for row in rows:
    if row.get("trans_d_scale_posterior_gamma") == "1" and row.get("trans_dscale_multiplier") == "0.5":
        base[seed_family(row)] = fnum(row.get("model_top1_accuracy_genome_trans"))

for row in rows:
    row["seed_family"] = seed_family(row)
    t = fnum(row.get("model_top1_accuracy_genome_trans"))
    local_base = base.get(row["seed_family"])
    row["delta_trans_vs_same_seed_g1_td0p5"] = (
        f"{t - local_base:.9g}" if t is not None and local_base is not None else "NA"
    )
    if t is not None:
        row["delta_trans_vs_current_best047"] = f"{t - CURRENT_BEST_TRANS:.9g}"
        row["target_plus_0p1_vs_current_best047_met"] = "1" if t >= CURRENT_BEST_TRANS + 0.1 else "0"
    else:
        row["delta_trans_vs_current_best047"] = "NA"
        row["target_plus_0p1_vs_current_best047_met"] = "0"

fields = [
    "config_name",
    "seed_family",
    "status",
    "backend",
    "init_seed",
    "init_noise_scale",
    "trans_d_scale_posterior_gamma",
    "trans_dscale_multiplier",
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
    "force_diag_total_contact_force_l1",
    "copytrack_force_l1_over_contact_force_l1",
    "global_copytrack_force_l1_over_contact_force_l1",
    "copytrack_frac_cos_lt_0",
    "delta_trans_vs_same_seed_g1_td0p5",
    "delta_trans_vs_current_best047",
    "target_plus_0p1_vs_current_best047_met",
]
rows_sorted = sorted(rows, key=lambda r: fnum(r.get("model_top1_accuracy_genome_trans")) or -1, reverse=True)
with (root / "trans_dscale_gamma_delta_summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows_sorted:
        writer.writerow({key: row.get(key, "NA") for key in fields})
headline = "PLUS_0P1_MET" if any(row.get("target_plus_0p1_vs_current_best047_met") == "1" for row in rows) else "NO_PLUS_0P1"
(root / "headline.txt").write_text(headline + "\n")
PY

log_command_line python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" "[inline README]"
python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
light = Path(sys.argv[2])
rows = list(csv.DictReader((root / "trans_dscale_gamma_delta_summary.tsv").open(), delimiter="\t"))
headline = (root / "headline.txt").read_text().strip()
fields = [
    "config_name",
    "seed_family",
    "trans_d_scale_posterior_gamma",
    "trans_dscale_multiplier",
    "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans",
    "delta_trans_vs_same_seed_g1_td0p5",
    "delta_trans_vs_current_best047",
    "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "final_mean_entropy",
    "final_mean_pU",
    "force_diag_total_contact_force_l1",
    "copytrack_frac_cos_lt_0",
    "target_plus_0p1_vs_current_best047_met",
]
with (root / "README.md").open("w") as out:
    out.write("# 049 P9016 Trans Dscale Gamma 1Mb\n\n")
    out.write("This controlled blind-training experiment tests whether trans-specific posterior-count dscale gamma can fix the softall trans geometry bottleneck without changing posterior calculation, priors, rho, temperature, or graph topology.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- headline: `{headline}`\n")
    out.write("- training inputs: approved P9016 raw pairs only\n")
    out.write("- eval inputs: P9016 SNP labels and CHARM/3DG are used only after training\n")
    out.write("- fixed knobs: posterior_count global gamma=1, lambda_copytrack=0.03, lambda_global_copytrack=0.003, no sep force, no callable anchor, no chr-pair M-step\n")
    out.write("- controlled knobs: trans_d_scale_posterior_gamma and trans_dscale_multiplier only\n")
    out.write("- success criterion: +0.1 full-denominator trans top1 over current best blind-safe standard trans, recorded as 0.358585\n\n")
    out.write("## Main Results\n\n")
    out.write("| " + " | ".join(fields) + " |\n")
    out.write("|" + "|".join(["---"] * len(fields)) + "|\n")
    for row in rows:
        out.write("| " + " | ".join(row.get(key, "NA") for key in fields) + " |\n")
    out.write("\n## Interpretation Boundary\n\n")
    out.write("A positive result would mean the current trans gap is sensitive to count-to-distance semantics in softall, not merely posterior assignment. A negative result means the largest known dscale lever is insufficient; then the remaining problem is more likely raw-pair identifiability or missing data linkage rather than another posterior self-training rule.\n")
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
if [ -s "$FULL_RUN_ROOT/trans_dscale_gamma_delta_summary.tsv" ]; then
	echo "TRANS_DSCALE_GAMMA_DELTA_SUMMARY_TSV=$FULL_RUN_ROOT/trans_dscale_gamma_delta_summary.tsv"
fi
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt" 2>/dev/null || echo SMOKE_RESULT_REMOVED)"
