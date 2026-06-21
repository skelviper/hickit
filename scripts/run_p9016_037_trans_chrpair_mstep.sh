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

RUN_ID=${HK_BLIND_037_RUN_ID:-"037-$(date +%Y%m%d_%H%M%S)-p9016_trans_chrpair_mstep_1m"}
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

PAIR_PATH=${HK_BLIND_P9016_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
EVAL_PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
APPROVED_P9016_PAIRS=/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_RELAX_STEP:-${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-0.012}}
BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-gpu}
MAX_PARALLEL=${HK_BLIND_037_MAX_PARALLEL:-${HK_BLIND_MAX_PARALLEL:-10}}
RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-1}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-1}
CONFIG_SET=${HK_BLIND_037_CONFIG_SET:-full}
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
		echo "error: HK_BLIND_037_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_037_MAX_PARALLEL must be >= 1" >&2
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
	if ! same_realpath "$PAIR_PATH" "$APPROVED_P9016_PAIRS" &&
	   [ "${HK_BLIND_037_ALLOW_CUSTOM_TRAINING_PAIRS:-0}" != "1" ]; then
		echo "error: formal 037 requires approved P9016 training pairs: $APPROVED_P9016_PAIRS" >&2
		echo "got: $PAIR_PATH" >&2
		exit 1
	fi
	if ! same_realpath "$EVAL_PAIRS" "$APPROVED_P9016_PAIRS" &&
	   [ "${HK_BLIND_037_ALLOW_CUSTOM_EVAL_PAIRS:-0}" != "1" ]; then
		echo "error: formal 037 requires approved P9016 eval pairs: $APPROVED_P9016_PAIRS" >&2
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

cp "$0" "$FULL_RUN_ROOT/scripts_snapshot/"
cp scripts/p9016_common_eval.sh \
   scripts/p9016_publish_light_result.sh \
   scripts/compute_p9016_copytrack_diag.py \
   scripts/summarize_p9016_model_sweep.py \
   "$FULL_RUN_ROOT/scripts_snapshot/"

log_cmd make smoke_blind_p9016_minimal
if [ "$MAKE_GPU" = "1" ]; then
	log_cmd make gpu=1 run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
else
	log_cmd make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
fi

SMOKE_ROOT=/tmp/hk_blind_p9016_037_chrpair_mstep_smoke
rm -rf "$SMOKE_ROOT"
log_command_line env HK_BLIND_P9016_CONFIG_NAME=chrpair_mstep_smoke ./run_blind_p9016_minimal.bin
env \
	HK_BLIND_P9016_PAIRS=testdata/p9016_blind_smoke.pairs \
	HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS=1 \
	HK_BLIND_P9016_OUTPUT_ROOT="$SMOKE_ROOT" \
	HK_BLIND_P9016_CONFIG_NAME=chrpair_mstep_smoke \
	HK_BLIND_P9016_MINIMAL_N_ITER=2 \
	HK_BLIND_P9016_MINIMAL_RELAX_STEPS=1 \
	HK_BLIND_P9016_D_SCALE_MODE=posterior_count \
	HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA=1 \
	HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_LAMBDA=0.5 \
	HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_POWER=2 \
	HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_WARMUP_ITER=1 \
	./run_blind_p9016_minimal.bin > "$LOG_ROOT/chrpair_mstep_smoke.log" 2>&1
log_cmd ./audit_blind_p9016_full_cpu_output.bin "$SMOKE_ROOT/chrpair_mstep_smoke"
grep -q $'trans_chr_pair_mstep_mode\tblind_posterior_chrom_pair_mstep' "$SMOKE_ROOT/chrpair_mstep_smoke/p9016_full.manifest.tsv"
grep -q $'trans_chr_pair_mstep_lambda\t0.5' "$SMOKE_ROOT/chrpair_mstep_smoke/p9016_full.manifest.tsv"
grep -q $'trans_chr_pair_mstep_power\t2' "$SMOKE_ROOT/chrpair_mstep_smoke/p9016_full.manifest.tsv"
grep -q $'trans_chr_pair_mstep_uses_phase_labels\t0' "$SMOKE_ROOT/chrpair_mstep_smoke/p9016_full.manifest.tsv"
grep -q $'trans_chr_pair_mstep_uses_charm_or_reference\t0' "$SMOKE_ROOT/chrpair_mstep_smoke/p9016_full.manifest.tsv"
rm -rf "$SMOKE_ROOT"

RUN_BIN="$LOG_ROOT/run_blind_p9016_minimal.bin"
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
	echo "backend	$BACKEND"
	echo "make_gpu	$MAKE_GPU"
	echo "bin_size_bp	$BIN_SIZE"
	echo "n_iter	$N_ITER"
	echo "relax_steps	$RELAX_STEPS"
	echo "relax_step	$RELAX_STEP"
	echo "max_parallel	$MAX_PARALLEL"
	echo "config_set	$CONFIG_SET"
	echo "training_boundary	P9016 pairs only; phase labels and CHARM/3DG are eval-only"
	echo "experiment_control	trans chromosome-pair posterior synchronization applied only to softall M-step training probabilities"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

run_config() {
	local config=$1
	local trans_dscale=$2
	local lambda_copytrack=$3
	local lambda_global=$4
	local mstep_lambda=$5
	local mstep_power=$6
	local mstep_warmup=$7
	local out_dir="$OUTPUT_ROOT/$config"
	local eval_dir="$EVAL_ROOT/$config"
	local work_root="$FULL_RUN_ROOT/work_outputs/$config"
	local produced_dir="$work_root/$config"
	local log_prefix="$LOG_ROOT/$config"
	local train_cmd
	local manifest

	echo "CONFIG_START=$config trans_dscale=$trans_dscale lambda_copytrack=$lambda_copytrack lambda_global=$lambda_global mstep_lambda=$mstep_lambda mstep_power=$mstep_power mstep_warmup=$mstep_warmup"
	rm -rf "$out_dir" "$eval_dir" "$work_root"
	mkdir -p "$eval_dir" "$work_root"

	train_cmd=(
		clean_env_prefix
		HK_BLIND_SAMPLE=P9016
		HK_BLIND_P9016_SAMPLE=P9016
		HK_BLIND_GIT_COMMIT="$GIT_COMMIT"
		HK_BLIND_GIT_DIRTY_COUNT="$GIT_DIRTY_COUNT"
		HK_BLIND_BINARY_HASH="$RUN_HASH"
		HK_BLIND_P9016_PAIRS="$PAIR_PATH"
		HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS="$ALLOW_CUSTOM"
		HK_BLIND_P9016_OUTPUT_ROOT="$work_root"
		HK_BLIND_P9016_CONFIG_NAME="$config"
		HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE"
		HK_BLIND_P9016_MINIMAL_N_ITER="$N_ITER"
		HK_BLIND_P9016_MINIMAL_RELAX_STEPS="$RELAX_STEPS"
		HK_BLIND_P9016_RELAX_STEP="$RELAX_STEP"
		HK_BLIND_P9016_RELAX_BACKEND="$BACKEND"
		HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split
		HK_BLIND_P9016_INIT_SEED=17
		HK_BLIND_P9016_INIT_EPS=0.5
		HK_BLIND_P9016_INIT_NOISE_SCALE=0
		HK_BLIND_P9016_D_SCALE_MODE=posterior_count
		HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA=1
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6
		HK_BLIND_P9016_MIN_SEP_UNIT=0
		HK_BLIND_P9016_LAMBDA_SEP=0
		HK_BLIND_P9016_LAMBDA_COPYTRACK="$lambda_copytrack"
		HK_BLIND_P9016_LAMBDA_GLOBAL_COPYTRACK="$lambda_global"
		HK_BLIND_P9016_TRANS_DSCALE_MULTIPLIER="$trans_dscale"
		HK_BLIND_P9016_TRANS_CHR_PAIR_PRIOR_LAMBDA=0
		HK_BLIND_P9016_TRANS_TOP1_MODE=off
		HK_BLIND_P9016_RHO_TRAIN_MODE=constant
		HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_LAMBDA="$mstep_lambda"
		HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_EPS=0.001
		HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_POWER="$mstep_power"
		HK_BLIND_P9016_TRANS_CHR_PAIR_MSTEP_WARMUP_ITER="$mstep_warmup"
		HK_BLIND_WRITE_RAW=0
		"$RUN_BIN"
	)
	log_command_line "${train_cmd[@]}"
	"${train_cmd[@]}" > "${log_prefix}.train.log" 2>&1
	if [ ! -d "$produced_dir" ]; then
		echo "missing config output directory after training: $produced_dir" >&2
		return 1
	fi
	mkdir -p "$OUTPUT_ROOT"
	rm -rf "$out_dir"
	mv "$produced_dir" "$out_dir"
	if [ -s "$work_root/matrix_summary.tsv" ]; then
		mv "$work_root/matrix_summary.tsv" "$out_dir/matrix_summary.tsv"
	fi
	if [ -s "$work_root/bmap_summary.tsv" ]; then
		mv "$work_root/bmap_summary.tsv" "$out_dir/bmap_summary.tsv"
	fi
	rmdir "$work_root" 2>/dev/null || true
	log_cmd ./audit_blind_p9016_full_cpu_output.bin "$out_dir" > "${log_prefix}.audit.log" 2>&1
	manifest="$out_dir/p9016_full.manifest.tsv"
	grep -q $'input_contact_source\traw_pairs' "$manifest"
	grep -q $'uses_phase_labels\t0' "$manifest"
	grep -q $'uses_charm_or_reference\t0' "$manifest"
	grep -q $'trans_chr_pair_prior_lambda\t0' "$manifest"
	grep -q $'trans_top1_mstep_mode\toff' "$manifest"
	grep -q $'trans_chr_pair_mstep_scope\ttrans_edges_only' "$manifest"
	grep -q $'trans_chr_pair_mstep_application_point\tpost_estep_pre_softall_mstep' "$manifest"
	grep -q $'trans_chr_pair_mstep_uses_phase_labels\t0' "$manifest"
	grep -q $'trans_chr_pair_mstep_uses_charm_or_reference\t0' "$manifest"
	if [ "$CONFIG_SET" != "smoke" ]; then
		grep -q $'approved_p9016_raw_pairs_realpath\t1' "$manifest"
	fi
	if [ "$RUN_EVAL" = "1" ]; then
		if [ ! -s "$EVAL_PAIRS" ] || [ ! -s "$EVAL_TDG" ]; then
			if [ "$REQUIRE_EVAL" = "1" ]; then
				echo "missing eval inputs for $config" >&2
				return 1
			fi
			printf 'key\tvalue\nstatus\tSKIPPED_EVAL_INPUT_MISSING\n' > "$eval_dir/eval_status.tsv"
		else
			HK_BLIND_SAMPLE=P9016 \
			HK_BLIND_P9016_SAMPLE=P9016 \
			HK_BLIND_P9016_EVAL_PAIRS="$EVAL_PAIRS" \
			HK_BLIND_P9016_EVAL_TDG="$EVAL_TDG" \
			HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE" \
			scripts/p9016_common_eval.sh "$out_dir" "$eval_dir" > "${log_prefix}.eval.log" 2>&1
			if [ -s "$out_dir/p9016_full.coords.tsv" ]; then
				scripts/compute_p9016_copytrack_diag.py \
					--config-output-dir "$out_dir" \
					--eval-output-dir "$eval_dir" \
					--out "$eval_dir/copytrack_vector_diag.tsv" > "${log_prefix}.copytrack.log" 2>&1
			fi
		fi
	fi
	echo "CONFIG_DONE=$config"
}

if [ "$CONFIG_SET" = "smoke" ]; then
	CONFIGS=(
		"p9016_chrpair_mstep_smoke_off|0.5|0|0|0|1|0"
		"p9016_chrpair_mstep_smoke_l0p5_p2|0.5|0|0|0.5|2|1"
	)
else
	CONFIGS=(
		"p9016_pcgamma1_td1_mstep0_sep0|1|0|0|0|1|0"
		"p9016_pcgamma1_td0p5_mstep0_sep0|0.5|0|0|0|1|0"
		"p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035|0.5|0.03|0.003|0|1|0"
		"p9016_pcgamma1_td0p5_mstep0p25_p1_w0_sep0|0.5|0|0|0.25|1|0"
		"p9016_pcgamma1_td0p5_mstep0p5_p1_w0_sep0|0.5|0|0|0.5|1|0"
		"p9016_pcgamma1_td0p5_mstep0p75_p1_w0_sep0|0.5|0|0|0.75|1|0"
		"p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0|0.5|0|0|1|1|0"
		"p9016_pcgamma1_td0p5_mstep0p5_p2_w10_sep0|0.5|0|0|0.5|2|10"
		"p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0|0.5|0|0|0.5|4|10"
		"p9016_pcgamma1_td0p5_mstep0p25_p4_w10_sep0|0.5|0|0|0.25|4|10"
		"p9016_pcgamma1_td0p5_mstep0p5_p2_w10_lct0p03_gct0p003|0.5|0.03|0.003|0.5|2|10"
		"p9016_pcgamma1_td0p5_mstep0p25_p2_w10_lct0p03_gct0p003|0.5|0.03|0.003|0.25|2|10"
	)
fi

failures=0
running=0
for spec in "${CONFIGS[@]}"; do
	IFS='|' read -r config trans_dscale lambda_copytrack lambda_global mstep_lambda mstep_power mstep_warmup <<< "$spec"
	run_config "$config" "$trans_dscale" "$lambda_copytrack" "$lambda_global" "$mstep_lambda" "$mstep_power" "$mstep_warmup" &
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

scripts/summarize_p9016_model_sweep.py "$FULL_RUN_ROOT" > "$FULL_RUN_ROOT/summary.tsv"

python - "$FULL_RUN_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))

def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None

base = next((r for r in rows if r.get("config_name") == "p9016_pcgamma1_td0p5_mstep0_sep0"), None)
base_t = fnum(base.get("model_top1_accuracy_genome_trans")) if base else None
fields = [
    "config_name", "trans_dscale_multiplier", "lambda_copytrack",
    "lambda_global_copytrack", "trans_chr_pair_mstep_lambda",
    "trans_chr_pair_mstep_power", "trans_chr_pair_mstep_warmup_iter",
    "model_top1_accuracy_genome_all", "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans", "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman", "final_mean_entropy",
    "final_mean_pU", "final_min_sep", "sep_p05", "final_mean_sep",
    "copytrack_frac_cos_lt_0", "delta_trans_vs_td0p5_mstep0",
    "target_plus_0p1_met",
]
for r in rows:
    t = fnum(r.get("model_top1_accuracy_genome_trans"))
    if base_t is not None and t is not None:
        r["delta_trans_vs_td0p5_mstep0"] = f"{t - base_t:.9g}"
        r["target_plus_0p1_met"] = "1" if t >= base_t + 0.1 else "0"
    else:
        r["delta_trans_vs_td0p5_mstep0"] = "NA"
        r["target_plus_0p1_met"] = "NA"
with (root / "trans_chrpair_mstep_delta_summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for r in rows:
        writer.writerow({k: r.get(k, "NA") for k in fields})
headline = "NO_PLUS_0P1"
if any(r.get("target_plus_0p1_met") == "1" for r in rows):
    headline = "PLUS_0P1_MET"
(root / "headline.txt").write_text(headline + "\n")
PY

python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
light = Path(sys.argv[2])
rows = list(csv.DictReader((root / "trans_chrpair_mstep_delta_summary.tsv").open(), delimiter="\t"))
def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None
rows_sorted = sorted(rows, key=lambda r: (fnum(r.get("model_top1_accuracy_genome_trans")) is None,
                                         -(fnum(r.get("model_top1_accuracy_genome_trans")) or -1)))
headline = (root / "headline.txt").read_text().strip()
with (root / "README.md").open("w") as out:
    out.write("# 037 P9016 trans chromosome-pair M-step synchronization\\n\\n")
    out.write("This controlled blind-training run tests whether chromosome-pair posterior aggregation can help trans identity when it is applied only to the softall M-step training graph. E-step posterior, prior mode, temperature, rho mode, and training input stay unchanged.\\n\\n")
    out.write(f"- full result root: `{root}`\\n")
    out.write(f"- light result root: `{light}`\\n")
    out.write(f"- headline: `{headline}`\\n")
    out.write("- training boundary: P9016 raw pairs only; SNP/phase and CHARM/3DG are eval-only.\\n\\n")
    out.write("## Top Results\\n\\n")
    cols = [
        "config_name", "trans_chr_pair_mstep_lambda", "trans_chr_pair_mstep_power",
        "model_top1_accuracy_genome_all", "model_top1_accuracy_genome_cis",
        "model_top1_accuracy_genome_trans", "model_same_cross_accuracy_genome_trans",
        "mean_per_chrom_cis_distance_spearman", "delta_trans_vs_td0p5_mstep0",
    ]
    out.write("| " + " | ".join(cols) + " |\\n")
    out.write("|" + "|".join(["---"] * len(cols)) + "|\\n")
    for r in rows_sorted[:12]:
        out.write("| " + " | ".join(r.get(c, "NA") for c in cols) + " |\\n")
    out.write("\\n## Interpretation Boundary\\n\\n")
    out.write("A real success requires full-denominator trans top1 to increase by at least 0.1 over the td0.5 M-step-off baseline without a major cis or cis-Spearman collapse. This experiment tests blind posterior self-synchronization; it does not use truth labels or reference structure during training.\\n")
PY

if [ "$CONFIG_SET" = "smoke" ]; then
	echo "SMOKE_RESULT_ROOT=$FULL_RUN_ROOT"
	rm -rf "$FULL_RUN_ROOT"
	echo "SMOKE_CLEANED=1"
	exit 0
fi

scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "TRANS_CHRPAIR_MSTEP_DELTA_TSV=$FULL_RUN_ROOT/trans_chrpair_mstep_delta_summary.tsv"
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt")"
