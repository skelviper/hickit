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

RUN_ID=${HK_BLIND_030_RUN_ID:-"030-$(date +%Y%m%d_%H%M%S)-p9016_global_copytrack_1m"}
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE3_ROOT/test_res}/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
OUTPUT_ROOT="$FULL_RUN_ROOT/outputs"
EVAL_ROOT="$FULL_RUN_ROOT/eval"
LOG_ROOT="$FULL_RUN_ROOT/logs"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"
mkdir -p "$OUTPUT_ROOT" "$EVAL_ROOT" "$LOG_ROOT" "$FULL_RUN_ROOT/scripts_snapshot"
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)

exec > >(tee -a "$LOG_ROOT/runner.stdout.log") 2> >(tee -a "$LOG_ROOT/runner.stderr.log" >&2)

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

PAIR_PATH=${HK_BLIND_P9016_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
EVAL_PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-$PAIR_PATH}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-${HK_BLIND_P9016_RELAX_STEP:-0.012}}
MAX_PARALLEL=${HK_BLIND_030_MAX_PARALLEL:-${HK_BLIND_MAX_PARALLEL:-8}}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-1}
RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-1}
ALLOW_CUSTOM=${HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS:-0}

case "$MAX_PARALLEL" in
	''|*[!0-9]*)
		echo "error: HK_BLIND_030_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_030_MAX_PARALLEL must be >= 1" >&2
	exit 2
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

slug_float() {
	case "$1" in
		0|0.0) echo "0" ;;
		0.0003) echo "0p0003" ;;
		0.001) echo "0p001" ;;
		0.003) echo "0p003" ;;
		0.01) echo "0p01" ;;
		0.03) echo "0p03" ;;
		0.1) echo "0p1" ;;
		0.3) echo "0p3" ;;
		0.5) echo "0p5" ;;
		1|1.0) echo "1" ;;
		1.5) echo "1p5" ;;
		*) echo "$1" | tr '.' 'p' ;;
	esac
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

log_cmd make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
log_cmd make smoke_blind_p9016_minimal

RUN_BIN="$FULL_RUN_ROOT/logs/run_blind_p9016_minimal.cpu.bin"
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
	echo "backend	cpu"
	echo "bin_size_bp	$BIN_SIZE"
	echo "n_iter	$N_ITER"
	echo "relax_steps	$RELAX_STEPS"
	echo "relax_step	$RELAX_STEP"
	echo "max_parallel	$MAX_PARALLEL"
	echo "training_boundary	raw P9016 pairs only; phase labels and CHARM/3DG are eval-only"
	echo "global_copytrack_semantics	chromosome_mean homolog-vector gauge alignment, blind geometry force"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

run_config() {
	local config=$1
	local gamma=$2
	local min_sep=$3
	local lambda_sep=$4
	local lambda_copytrack=$5
	local lambda_global=$6
	local out_dir="$OUTPUT_ROOT/$config"
	local eval_dir="$EVAL_ROOT/$config"
	local log_prefix="$LOG_ROOT/$config"
	local train_cmd audit_cmd eval_cmd copytrack_cmd

	echo "CONFIG_START=$config"
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
		HK_BLIND_P9016_RELAX_BACKEND=cpu
		HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split
		HK_BLIND_P9016_INIT_SCALE=10
		HK_BLIND_P9016_INIT_SEED=17
		HK_BLIND_P9016_INIT_EPS=0.5
		HK_BLIND_P9016_INIT_NOISE_SCALE=0.0
		HK_BLIND_P9016_MIN_SEP_UNIT="$min_sep"
		HK_BLIND_P9016_LAMBDA_SEP="$lambda_sep"
		HK_BLIND_P9016_LAMBDA_COPYTRACK="$lambda_copytrack"
		HK_BLIND_P9016_LAMBDA_GLOBAL_COPYTRACK="$lambda_global"
		HK_BLIND_P9016_D_SCALE_MODE=posterior_count
		HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA="$gamma"
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6
		HK_BLIND_WRITE_RAW=0
		"$RUN_BIN"
	)
	log_command_line "${train_cmd[@]}"
	"${train_cmd[@]}" > "$log_prefix.train.log" 2>&1

	audit_cmd=(clean_env_prefix ./audit_blind_p9016_full_cpu_output.bin "$out_dir")
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

PIDS=()
LABELS=()

wait_for_slot() {
	while [ "${#PIDS[@]}" -ge "$MAX_PARALLEL" ]; do
		local pid=${PIDS[0]}
		local label=${LABELS[0]}
		if wait "$pid"; then
			echo "CONFIG_STATUS=$label OK"
		else
			echo "CONFIG_STATUS=$label FAILED" >&2
			return 1
		fi
		PIDS=("${PIDS[@]:1}")
		LABELS=("${LABELS[@]:1}")
	done
}

launch_config() {
	local config=$1
	shift
	wait_for_slot
	run_config "$config" "$@" &
	PIDS+=("$!")
	LABELS+=("$config")
}

launch_config p9016_pcgamma1_sep0_lct0_gct0_ieps0p5_noise0_seed17 1 0 0 0 0

for lambda_copytrack in 0.001 0.003 0.01 0.03 0.1; do
	cslug=$(slug_float "$lambda_copytrack")
	launch_config "p9016_pcgamma1_sep0_lct${cslug}_gct0_ieps0p5_noise0_seed17" 1 0 0 "$lambda_copytrack" 0
done

for lambda_global in 0.001 0.003 0.01 0.03 0.1; do
	gslug=$(slug_float "$lambda_global")
	launch_config "p9016_pcgamma1_sep0_lct0_gct${gslug}_ieps0p5_noise0_seed17" 1 0 0 0 "$lambda_global"
done

for spec in \
	"0.003 0.003" \
	"0.01 0.003" \
	"0.003 0.01" \
	"0.01 0.01" \
	"0.03 0.01" \
	"0.01 0.03"
do
	set -- $spec
	lct=$1
	gct=$2
	launch_config "p9016_pcgamma1_sep0_lct$(slug_float "$lct")_gct$(slug_float "$gct")_ieps0p5_noise0_seed17" 1 0 0 "$lct" "$gct"
done

if [ "${HK_BLIND_030_RUN_SEP:-1}" = "1" ]; then
	for spec in \
		"1 0.5 0 0.01" \
		"1 1 0 0.01" \
		"1 0.5 0.01 0.01" \
		"1 1 0.01 0.01" \
		"1 0.5 0.03 0.01" \
		"1 1 0.03 0.01" \
		"1.5 1 0.01 0.01"
	do
		set -- $spec
		min_sep=$1
		lsep=$2
		lct=$3
		gct=$4
		launch_config "p9016_pcgamma1_msep$(slug_float "$min_sep")_lsep$(slug_float "$lsep")_lct$(slug_float "$lct")_gct$(slug_float "$gct")_ieps0p5_noise0_seed17" 1 "$min_sep" "$lsep" "$lct" "$gct"
	done
fi

failed=0
for i in "${!PIDS[@]}"; do
	if wait "${PIDS[$i]}"; then
		echo "CONFIG_STATUS=${LABELS[$i]} OK"
	else
		echo "CONFIG_STATUS=${LABELS[$i]} FAILED" >&2
		failed=1
	fi
done
if [ "$failed" -ne 0 ]; then
	echo "one or more configs failed" >&2
	exit 1
fi

scripts/summarize_p9016_model_sweep.py "$FULL_RUN_ROOT" > "$FULL_RUN_ROOT/summary.tsv"

python - "$FULL_RUN_ROOT" <<'PY'
from __future__ import annotations
import csv
import math
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))
base = next((r for r in rows if r.get("config_name") == "p9016_pcgamma1_sep0_lct0_gct0_ieps0p5_noise0_seed17"), None)

def num(value: str | None) -> float:
    try:
        x = float(value if value not in (None, "", "NA") else "nan")
    except ValueError:
        return float("nan")
    return x

base_trans = num(base.get("model_top1_accuracy_genome_trans")) if base else float("nan")
for row in rows:
    trans = num(row.get("model_top1_accuracy_genome_trans"))
    if math.isfinite(trans) and math.isfinite(base_trans):
        row["delta_trans_vs_baseline"] = f"{trans - base_trans:.9g}"
        row["target_plus_0p1_met"] = "1" if trans >= base_trans + 0.1 else "0"
    else:
        row["delta_trans_vs_baseline"] = "NA"
        row["target_plus_0p1_met"] = "0"

fields = [
    "config_name",
    "lambda_copytrack",
    "lambda_global_copytrack",
    "min_sep_unit",
    "lambda_sep",
    "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans",
    "delta_trans_vs_baseline",
    "target_plus_0p1_met",
    "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "final_mean_entropy",
    "final_mean_pU",
    "final_min_sep",
    "sep_p05",
    "final_mean_sep",
    "copytrack_frac_cos_lt_0",
    "copytrack_frac_projection_sign_switch",
    "copytrack_force_l1_over_contact_force_l1",
    "global_copytrack_force_l1_over_contact_force_l1",
]
rows_sorted = sorted(rows, key=lambda r: num(r.get("model_top1_accuracy_genome_trans")), reverse=True)
with (root / "trans_delta_summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    for row in rows_sorted:
        writer.writerow({field: row.get(field, "NA") for field in fields})
headline = "PLUS_0P1_MET" if any(r.get("target_plus_0p1_met") == "1" for r in rows) else "NO_PLUS_0P1"
(root / "headline.txt").write_text(headline + "\n")
PY

python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" <<'PY'
from __future__ import annotations
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
light = Path(sys.argv[2])
rows = list(csv.DictReader((root / "trans_delta_summary.tsv").open(), delimiter="\t"))
headline = (root / "headline.txt").read_text().strip()
fields = [
    "config_name",
    "lambda_copytrack",
    "lambda_global_copytrack",
    "min_sep_unit",
    "lambda_sep",
    "model_top1_accuracy_genome_trans",
    "delta_trans_vs_baseline",
    "model_top1_accuracy_genome_cis",
    "mean_per_chrom_cis_distance_spearman",
    "final_mean_entropy",
    "final_mean_pU",
    "final_min_sep",
    "sep_p05",
    "final_mean_sep",
    "copytrack_frac_cos_lt_0",
    "global_copytrack_force_l1_over_contact_force_l1",
]
with (root / "README.md").open("w") as out:
    out.write("# 030 P9016 Global Copytrack 1Mb\n\n")
    out.write("This controlled blind-training experiment tests whether local and chromosome-mean homolog-vector gauge alignment can improve P9016 trans four-state accuracy.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- summary: `{root / 'summary.tsv'}`\n")
    out.write(f"- trans delta summary: `{root / 'trans_delta_summary.tsv'}`\n")
    out.write(f"- headline: `{headline}`\n")
    out.write("- training input: approved raw P9016 pairs only\n")
    out.write("- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper\n")
    out.write("- backend: CPU for this experiment, because GPU extra-force semantics use operator splitting and are not used for conclusions here\n\n")
    out.write("## Main Results\n\n")
    out.write("| " + " | ".join(fields) + " |\n")
    out.write("|" + "|".join(["---"] * len(fields)) + "|\n")
    for row in rows:
        out.write("| " + " | ".join(row.get(field, "NA") for field in fields) + " |\n")
    out.write("\n## Interpretation Boundary\n\n")
    out.write("Global copytrack is an internal gauge-alignment force on chromosome-mean copy1-minus-copy0 vectors. It does not assign maternal/paternal identity and does not use phase labels or CHARM/3DG during training. A real success requires trans top1 to increase by at least 0.1 without a major cis/Spearman collapse.\n")
PY

scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$FULL_RUN_ROOT/publish.log"

{
	echo "end_time	$(date -Is)"
	echo "headline	$(cat "$FULL_RUN_ROOT/headline.txt")"
} >> "$FULL_RUN_ROOT/run_manifest.tsv"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "TRANS_DELTA_TSV=$FULL_RUN_ROOT/trans_delta_summary.tsv"
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt")"
