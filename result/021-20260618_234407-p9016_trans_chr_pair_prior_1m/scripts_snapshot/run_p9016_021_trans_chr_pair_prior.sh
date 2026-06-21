#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PHASE3_ROOT=$(cd "$REPO_ROOT/.." && pwd)

TEST_RES_ROOT=${HK_BLIND_TEST_RES_ROOT:-"$PHASE3_ROOT/test_res"}
RUN_ID=${HK_BLIND_021_RUN_ID:-"021-$(date +%Y%m%d_%H%M%S)-p9016_trans_chr_pair_prior_1m"}
FULL_RUN_ROOT="$TEST_RES_ROOT/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
OUTPUT_ROOT="$FULL_RUN_ROOT/outputs"
WORK_OUTPUT_ROOT="$FULL_RUN_ROOT/work_outputs"
EVAL_ROOT="$FULL_RUN_ROOT/eval"
LOG_ROOT="$FULL_RUN_ROOT/logs"
MAX_PARALLEL=${HK_BLIND_021_MAX_PARALLEL:-10}

PAIRS=${HK_BLIND_P9016_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
EVAL_PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-$PAIRS}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_RELAX_STEP:-${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-0.012}}
BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-auto}

if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "HK_BLIND_021_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi

mkdir -p "$OUTPUT_ROOT" "$WORK_OUTPUT_ROOT" "$EVAL_ROOT" "$LOG_ROOT" "$FULL_RUN_ROOT/scripts_snapshot"
cp "$0" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/p9016_common_eval.sh" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/p9016_publish_light_result.sh" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/summarize_p9016_model_sweep.py" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/compute_p9016_copytrack_diag.py" "$FULL_RUN_ROOT/scripts_snapshot/"

echo "FULL_RUN_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "pairs	$PAIRS"
	echo "eval_pairs	$EVAL_PAIRS"
	echo "eval_tdg	$EVAL_TDG"
	echo "bin_size_bp	$BIN_SIZE"
	echo "n_iter	$N_ITER"
	echo "relax_steps	$RELAX_STEPS"
	echo "relax_step	$RELAX_STEP"
	echo "backend	$BACKEND"
	echo "max_parallel	$MAX_PARALLEL"
	echo "git_commit	$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo NA)"
	echo "git_dirty_count	$(git -C "$REPO_ROOT" status --short 2>/dev/null | wc -l | awk '{print $1}')"
} > "$FULL_RUN_ROOT/run_manifest.tsv"
git -C "$REPO_ROOT" status --short > "$LOG_ROOT/git_status.txt" || true
git -C "$REPO_ROOT" diff --stat > "$LOG_ROOT/git_diff_stat.txt" || true
git -C "$REPO_ROOT" diff > "$LOG_ROOT/git_diff.patch" || true

log_cmd() {
	printf '%q ' "$@" >> "$FULL_RUN_ROOT/commands.log"
	printf '\n' >> "$FULL_RUN_ROOT/commands.log"
}

cd "$REPO_ROOT"
log_cmd make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
log_cmd make smoke_blind_p9016_minimal
make smoke_blind_p9016_minimal

CONFIGS=(
	"p9016_pcgamma1_baseline_sep_off_copytrack0_prior0|1|0|0|0|0"
	"p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0|1|1|0.5|0.03|0"
	"p9016_pcgamma1_prior0p25_sep_off|1|0|0|0|0.25"
	"p9016_pcgamma1_prior0p5_sep_off|1|0|0|0|0.5"
	"p9016_pcgamma1_prior0p75_sep_off|1|0|0|0|0.75"
	"p9016_pcgamma1_prior1_sep_off|1|0|0|0|1"
	"p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p25|1|1|0.5|0.03|0.25"
	"p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p5|1|1|0.5|0.03|0.5"
	"p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p75|1|1|0.5|0.03|0.75"
	"p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior1|1|1|0.5|0.03|1"
	"p9016_pcgamma0p5_prior0p5_sep_off|0.5|0|0|0|0.5"
	"p9016_pcgamma0p5_msep1_lsep0p5_copytrack0p01_prior0p5|0.5|1|0.5|0.01|0.5"
)

run_one() {
	local spec=$1
	IFS='|' read -r config gamma min_sep lambda_sep lambda_copytrack prior_lambda <<< "$spec"
	local run_root="$WORK_OUTPUT_ROOT/$config"
	local out_dir="$run_root/$config"
	local link_dir="$OUTPUT_ROOT/$config"
	local eval_dir="$EVAL_ROOT/$config"
	local log_file="$LOG_ROOT/${config}.log"
	rm -rf "$run_root" "$link_dir" "$eval_dir"
	mkdir -p "$run_root" "$eval_dir"
	{
		echo "CONFIG=$config"
		echo "gamma=$gamma min_sep=$min_sep lambda_sep=$lambda_sep lambda_copytrack=$lambda_copytrack prior_lambda=$prior_lambda"
	} > "$log_file"
	log_cmd env HK_BLIND_P9016_CONFIG_NAME="$config" HK_BLIND_P9016_OUTPUT_ROOT="$run_root" ./run_blind_p9016_minimal.bin
	(
		export HK_BLIND_SAMPLE=P9016
		export HK_BLIND_P9016_SAMPLE=P9016
		export HK_BLIND_P9016_PAIRS="$PAIRS"
		export HK_BLIND_P9016_EVAL_PAIRS="$EVAL_PAIRS"
		export HK_BLIND_P9016_EVAL_TDG="$EVAL_TDG"
		export HK_BLIND_P9016_OUTPUT_ROOT="$run_root"
		export HK_BLIND_P9016_CONFIG_NAME="$config"
		export HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE"
		export HK_BLIND_P9016_MINIMAL_N_ITER="$N_ITER"
		export HK_BLIND_P9016_MINIMAL_RELAX_STEPS="$RELAX_STEPS"
		export HK_BLIND_P9016_RELAX_STEP="$RELAX_STEP"
		export HK_BLIND_P9016_RELAX_BACKEND="$BACKEND"
		export HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split
		export HK_BLIND_P9016_INIT_SEED=17
		export HK_BLIND_P9016_INIT_EPS=0.5
		export HK_BLIND_P9016_INIT_NOISE_SCALE=0
		export HK_BLIND_P9016_D_SCALE_MODE=posterior_count
		export HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA="$gamma"
		export HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6
		export HK_BLIND_P9016_MIN_SEP_UNIT="$min_sep"
		export HK_BLIND_P9016_LAMBDA_SEP="$lambda_sep"
		export HK_BLIND_P9016_LAMBDA_COPYTRACK="$lambda_copytrack"
		export HK_BLIND_P9016_TRANS_CHR_PAIR_PRIOR_LAMBDA="$prior_lambda"
		export HK_BLIND_P9016_TRANS_CHR_PAIR_PRIOR_EPS=0.001
		./run_blind_p9016_minimal.bin >> "$log_file" 2>&1
		ln -sfn "../work_outputs/$config/$config" "$link_dir"
		./audit_blind_p9016_full_cpu_output.bin "$out_dir" >> "$log_file" 2>&1
		scripts/p9016_common_eval.sh "$out_dir" "$eval_dir" >> "$log_file" 2>&1
		if [ -s "$out_dir/p9016_full.coords.tsv" ]; then
			scripts/compute_p9016_copytrack_diag.py \
				--config-output-dir "$out_dir" \
				--eval-output-dir "$eval_dir" \
				--out "$eval_dir/copytrack_vector_diag.tsv" >> "$log_file" 2>&1
		fi
	) || {
		echo "FAILED	$config" > "$eval_dir/eval_status.tsv"
		return 1
	}
}

failures=0
running=0
for spec in "${CONFIGS[@]}"; do
	run_one "$spec" &
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
import csv, sys
from pathlib import Path
root = Path(sys.argv[1])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))
base = next((r for r in rows if r["config_name"] == "p9016_pcgamma1_baseline_sep_off_copytrack0_prior0"), None)
best = None
def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None
if base:
    b = fnum(base.get("model_top1_accuracy_genome_trans"))
    if b is not None:
        for r in rows:
            t = fnum(r.get("model_top1_accuracy_genome_trans"))
            if t is None:
                continue
            r["delta_trans_vs_baseline"] = f"{t-b:.9g}"
            r["target_plus_0p1_met"] = "1" if t >= b + 0.1 else "0"
            if best is None or t > (fnum(best.get("model_top1_accuracy_genome_trans")) or float("-inf")):
                best = r
with (root / "trans_delta.tsv").open("w", newline="") as fh:
    fields = [
        "config_name", "d_scale_posterior_gamma", "min_sep_unit", "lambda_sep",
        "lambda_copytrack", "trans_chr_pair_prior_lambda",
        "model_top1_accuracy_genome_all", "model_top1_accuracy_genome_cis",
        "model_top1_accuracy_genome_trans", "model_same_cross_accuracy_genome_trans",
        "mean_per_chrom_cis_distance_spearman", "final_mean_entropy", "final_mean_pU",
        "final_min_sep", "sep_p05", "final_mean_sep", "delta_trans_vs_baseline",
        "target_plus_0p1_met",
    ]
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for r in rows:
        out = {k: r.get(k, "NA") for k in fields}
        writer.writerow(out)
headline = "NO_PLUS_0P1"
if any(r.get("target_plus_0p1_met") == "1" for r in rows):
    headline = "PLUS_0P1_MET"
(root / "headline.txt").write_text(headline + "\n")
PY

python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" <<'PY'
import csv, sys
from pathlib import Path
root = Path(sys.argv[1])
light = Path(sys.argv[2])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))
def fnum(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return float("-inf")
    return x
rows_sorted = sorted(rows, key=lambda r: fnum(r.get("model_top1_accuracy_genome_trans")), reverse=True)
headline = (root / "headline.txt").read_text().strip()
with (root / "README.md").open("w") as out:
    out.write("# 021 P9016 Trans Chromosome-Pair Prior 1Mb\n\n")
    out.write("This controlled blind-training experiment tests whether a posterior-derived trans chromosome-pair prior can recover the pair-specific gauge space identified by the eval-only oracle diagnostic.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- headline: `{headline}`\n")
    out.write("- training input: raw P9016 pairs only\n")
    out.write("- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper\n\n")
    out.write("## Main Results\n\n")
    fields = [
        "config_name", "d_scale_posterior_gamma", "min_sep_unit", "lambda_sep",
        "lambda_copytrack", "trans_chr_pair_prior_lambda",
        "model_top1_accuracy_genome_all", "model_top1_accuracy_genome_cis",
        "model_top1_accuracy_genome_trans", "model_same_cross_accuracy_genome_trans",
        "mean_per_chrom_cis_distance_spearman", "final_mean_entropy", "final_mean_pU",
    ]
    out.write("| " + " | ".join(fields) + " |\n")
    out.write("|" + "|".join(["---"] * len(fields)) + "|\n")
    for r in rows_sorted:
        out.write("| " + " | ".join(r.get(k, "NA") for k in fields) + " |\n")
    out.write("\n## Interpretation Boundary\n\n")
    out.write("The trans chromosome-pair prior is blind because it is estimated from previous posterior probabilities and raw contact counts only. It does not use phase labels or CHARM/3DG during training. If the +0.1 target is not met, the eval-only pair oracle is not self-bootstrappable from the current posterior family.\n")
PY

scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$FULL_RUN_ROOT/publish.log"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "TRANS_DELTA_TSV=$FULL_RUN_ROOT/trans_delta.tsv"
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt")"
