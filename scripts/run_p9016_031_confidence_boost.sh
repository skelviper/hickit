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

RUN_ID=${HK_BLIND_031_RUN_ID:-"031-$(date +%Y%m%d_%H%M%S)-p9016_confidence_boost_pairs_1m"}
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE3_ROOT/test_res}/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
OUTPUT_ROOT="$FULL_RUN_ROOT/outputs"
EVAL_ROOT="$FULL_RUN_ROOT/eval"
LOG_ROOT="$FULL_RUN_ROOT/logs"
GENERATED_ROOT="$FULL_RUN_ROOT/generated_pairs"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"
if [ -e "$FULL_RUN_ROOT" ] &&
   [ "${HK_BLIND_ALLOW_NONEMPTY_RUN_ROOT:-0}" != "1" ] &&
   find "$FULL_RUN_ROOT" -mindepth 1 -print -quit | grep -q .; then
	echo "error: run root already exists and is not empty: $FULL_RUN_ROOT" >&2
	echo "set HK_BLIND_ALLOW_NONEMPTY_RUN_ROOT=1 only for an explicit resume/debug run" >&2
	exit 1
fi
mkdir -p "$OUTPUT_ROOT" "$EVAL_ROOT" "$LOG_ROOT" "$GENERATED_ROOT" "$FULL_RUN_ROOT/scripts_snapshot"
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
RELAX_STEP=${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-${HK_BLIND_P9016_RELAX_STEP:-0.012}}
BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-gpu}
MAX_PARALLEL=${HK_BLIND_031_MAX_PARALLEL:-${HK_BLIND_MAX_PARALLEL:-10}}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-1}
RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-1}
ALLOW_CUSTOM=${HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS:-0}
CONFIG_SET=${HK_BLIND_031_CONFIG_SET:-full}
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
		echo "error: HK_BLIND_031_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_031_MAX_PARALLEL must be >= 1" >&2
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
	   [ "${HK_BLIND_031_ALLOW_CUSTOM_TRAINING_PAIRS:-0}" != "1" ]; then
		echo "error: formal 031 requires approved P9016 training pairs: $APPROVED_P9016_PAIRS" >&2
		echo "got: $PAIR_PATH" >&2
		echo "set HK_BLIND_031_ALLOW_CUSTOM_TRAINING_PAIRS=1 only for an explicit diagnostic run" >&2
		exit 1
	fi
	if ! same_realpath "$EVAL_PAIRS" "$APPROVED_P9016_PAIRS" &&
	   [ "${HK_BLIND_031_ALLOW_CUSTOM_EVAL_PAIRS:-0}" != "1" ]; then
		echo "error: formal 031 requires approved P9016 eval pairs: $APPROVED_P9016_PAIRS" >&2
		echo "got: $EVAL_PAIRS" >&2
		echo "set HK_BLIND_031_ALLOW_CUSTOM_EVAL_PAIRS=1 only for an explicit diagnostic run" >&2
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

slug_float() {
	case "$1" in
		0|0.0) echo "0" ;;
		0.001) echo "0p001" ;;
		0.002) echo "0p002" ;;
		0.005) echo "0p005" ;;
		0.01) echo "0p01" ;;
		0.02) echo "0p02" ;;
		0.05) echo "0p05" ;;
		0.1) echo "0p1" ;;
		0.2) echo "0p2" ;;
		0.5) echo "0p5" ;;
		1|1.0) echo "1" ;;
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
cp scripts/export_p9016_confidence_boost_pairs.py \
   scripts/p9016_common_eval.sh \
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

RUN_BIN="$FULL_RUN_ROOT/logs/run_blind_p9016_minimal.bin"
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
	echo "training_boundary	P9016 pairs only; boosted pairs are derived from blind baseline posterior; phase labels and CHARM/3DG are eval-only"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

run_config() {
	local config=$1
	local pairs_path=$2
	local allow_custom=$3
	local min_sep=$4
	local lambda_sep=$5
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
		HK_BLIND_P9016_PAIRS="$pairs_path"
		HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS="$allow_custom"
		HK_BLIND_P9016_INPUT_CONTACT_SOURCE=raw_pairs
		HK_BLIND_P9016_OUTPUT_ROOT="$OUTPUT_ROOT"
		HK_BLIND_P9016_CONFIG_NAME="$config"
		HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE"
		HK_BLIND_P9016_MINIMAL_N_ITER="$N_ITER"
		HK_BLIND_P9016_MINIMAL_RELAX_STEPS="$RELAX_STEPS"
		HK_BLIND_P9016_RELAX_STEP="$RELAX_STEP"
		HK_BLIND_P9016_RELAX_BACKEND="$BACKEND"
		HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split
		HK_BLIND_P9016_INIT_SCALE=10
		HK_BLIND_P9016_INIT_SEED=17
		HK_BLIND_P9016_INIT_EPS=0.5
		HK_BLIND_P9016_INIT_NOISE_SCALE=0.0
		HK_BLIND_P9016_MIN_SEP_UNIT="$min_sep"
		HK_BLIND_P9016_LAMBDA_SEP="$lambda_sep"
		HK_BLIND_P9016_LAMBDA_COPYTRACK=0
		HK_BLIND_P9016_LAMBDA_GLOBAL_COPYTRACK=0
		HK_BLIND_P9016_D_SCALE_MODE=posterior_count
		HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA=1
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

derive_pairs() {
	local score=$1
	local frac=$2
	local boost=$3
	local name="p9016_confboost_${score}_top$(slug_float "$frac")_b${boost}"
	local pairs_out="$GENERATED_ROOT/${name}.pairs.gz"
	local meta_out="$GENERATED_ROOT/${name}.metadata.tsv"
	local keys_out="$GENERATED_ROOT/${name}.selected_keys.json"
	local cmd

	if [ ! -s "$pairs_out" ] || [ ! -s "$meta_out" ]; then
		cmd=(
			python scripts/export_p9016_confidence_boost_pairs.py
			--pairs "$PAIR_PATH"
			--posterior "$BASELINE_POSTERIOR"
			--out-pairs "$pairs_out"
			--metadata "$meta_out"
			--selected-keys-json "$keys_out"
			--bin-size "$BIN_SIZE"
			--score "$score"
			--trans-top-frac "$frac"
			--boost-copies "$boost"
			--phase-mode dot
		)
		log_command_line "${cmd[@]}"
		"${cmd[@]}" > "$LOG_ROOT/${name}.export.log" 2>&1
	fi
	printf '%s\n' "$pairs_out"
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

BASELINE_CONFIG=p9016_pcgamma1_approved_baseline_ieps0p5_noise0_seed17
run_config "$BASELINE_CONFIG" "$PAIR_PATH" "$ALLOW_CUSTOM" 0 0
BASELINE_POSTERIOR="$OUTPUT_ROOT/$BASELINE_CONFIG/p9016_full.bpair_posterior.tsv"
if [ ! -s "$BASELINE_POSTERIOR" ]; then
	echo "missing baseline posterior: $BASELINE_POSTERIOR" >&2
	exit 1
fi

launch_config p9016_pcgamma1_approved_msep1_lsep1_ieps0p5_noise0_seed17 "$PAIR_PATH" "$ALLOW_CUSTOM" 1 1

if [ "$CONFIG_SET" = "smoke" ]; then
	SPECS=(
		"pmax_margin 0.01 2 0 0"
	)
else
	SPECS=(
		"pmax_margin 0.005 2 0 0"
		"pmax_margin 0.01 2 0 0"
		"pmax_margin 0.02 2 0 0"
		"pmax_margin 0.05 2 0 0"
		"pmax_margin 0.1 2 0 0"
		"pmax_margin 0.2 2 0 0"
		"pmax_margin 0.01 5 0 0"
		"pmax_margin 0.02 5 0 0"
		"pmax_margin 0.05 5 0 0"
		"pmax_margin 0.1 5 0 0"
		"pmax 0.02 2 0 0"
		"pmax 0.05 2 0 0"
		"pmax 0.02 5 0 0"
		"pmax 0.05 5 0 0"
		"pmax_margin_entropy 0.02 2 0 0"
		"pmax_margin_entropy 0.05 2 0 0"
		"pmax_margin_entropy 0.02 5 0 0"
		"pmax_margin_entropy 0.05 5 0 0"
		"pmax_margin 0.02 2 1 1"
		"pmax_margin 0.05 2 1 1"
		"pmax_margin 0.1 2 1 1"
		"pmax_margin 0.02 5 1 1"
		"pmax_margin 0.05 5 1 1"
	)
fi

for spec in "${SPECS[@]}"; do
	set -- $spec
	score=$1
	frac=$2
	boost=$3
	min_sep=$4
	lsep=$5
	frac_slug=$(slug_float "$frac")
	pairs_path=$(derive_pairs "$score" "$frac" "$boost")
	if [ "$min_sep" = "0" ] && [ "$lsep" = "0" ]; then
		config="p9016_pcgamma1_confboost_${score}_top${frac_slug}_b${boost}_sep0_ieps0p5_noise0_seed17"
	else
		config="p9016_pcgamma1_confboost_${score}_top${frac_slug}_b${boost}_msep$(slug_float "$min_sep")_lsep$(slug_float "$lsep")_ieps0p5_noise0_seed17"
	fi
	launch_config "$config" "$pairs_path" 1 "$min_sep" "$lsep"
done

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

python - "$FULL_RUN_ROOT" "$BASELINE_CONFIG" <<'PY'
from __future__ import annotations
import csv
import math
import sys
from pathlib import Path

root = Path(sys.argv[1])
baseline_config = sys.argv[2]
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))

def num(value: str | None) -> float:
    try:
        return float(value if value not in (None, "", "NA") else "nan")
    except ValueError:
        return float("nan")

def read_kv(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header != ["key", "value"]:
            return {}
        return {row[0]: row[1] for row in reader if len(row) >= 2}

meta_by_stem = {p.name.replace(".metadata.tsv", ""): read_kv(p) for p in (root / "generated_pairs").glob("*.metadata.tsv")}
base = next((r for r in rows if r.get("config_name") == baseline_config), None)
base_trans = num(base.get("model_top1_accuracy_genome_trans")) if base else float("nan")
base_all = num(base.get("model_top1_accuracy_genome_all")) if base else float("nan")

out_rows = []
for row in rows:
    config = row.get("config_name", "")
    trans = num(row.get("model_top1_accuracy_genome_trans"))
    all_acc = num(row.get("model_top1_accuracy_genome_all"))
    row = dict(row)
    row["delta_trans_vs_baseline"] = f"{trans - base_trans:.9g}" if math.isfinite(trans) and math.isfinite(base_trans) else "NA"
    row["delta_all_vs_baseline"] = f"{all_acc - base_all:.9g}" if math.isfinite(all_acc) and math.isfinite(base_all) else "NA"
    row["target_plus_0p1_met"] = "1" if math.isfinite(trans) and math.isfinite(base_trans) and trans >= base_trans + 0.1 else "0"
    matched = {}
    for stem, meta in meta_by_stem.items():
        token = stem.replace("p9016_confboost_", "confboost_")
        if token in config:
            matched = meta
            break
    for key in [
        "score_mode",
        "trans_top_frac",
        "boost_copies",
        "selected_trans_bpairs",
        "selected_trans_raw_pairs",
        "extra_copies_written",
        "output_pairs",
        "score_threshold",
        "phase_mode",
    ]:
        row[key] = matched.get(key, "NA")
    out_rows.append(row)

fields = [
    "config_name",
    "score_mode",
    "trans_top_frac",
    "boost_copies",
    "min_sep_unit",
    "lambda_sep",
    "selected_trans_bpairs",
    "selected_trans_raw_pairs",
    "extra_copies_written",
    "output_pairs",
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
]
rows_sorted = sorted(out_rows, key=lambda r: num(r.get("model_top1_accuracy_genome_trans")), reverse=True)
with (root / "trans_delta_summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    for row in rows_sorted:
        writer.writerow({field: row.get(field, "NA") for field in fields})
headline = "PLUS_0P1_MET" if any(r.get("target_plus_0p1_met") == "1" for r in out_rows) else "NO_PLUS_0P1"
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
    "score_mode",
    "trans_top_frac",
    "boost_copies",
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
]
with (root / "README.md").open("w") as out:
    out.write("# 031 P9016 Confidence-Boosted Trans Pairs 1Mb\n\n")
    out.write("This controlled blind-training experiment tests whether high-confidence trans contacts from a baseline posterior can act as a self-training scaffold.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- summary: `{root / 'summary.tsv'}`\n")
    out.write(f"- trans delta summary: `{root / 'trans_delta_summary.tsv'}`\n")
    out.write(f"- headline: `{headline}`\n")
    out.write("- resolution: 1 Mb (`bin_size_bp=1000000` unless explicitly overridden in `run_manifest.tsv`)\n")
    out.write("- training input: P9016 pairs only; boosted pairs duplicate original raw contacts selected by blind baseline posterior confidence\n")
    out.write("- training boundary: SNP phase labels and CHARM/3DG are not used for training; generated training pairs write phase columns as `.`\n")
    out.write("- eval-only inputs: original P9016 pairs with SNP labels and CHARM/3DG are used only by the standard eval wrapper after training\n\n")
    out.write("## Main Results\n\n")
    out.write("| " + " | ".join(fields) + " |\n")
    out.write("|" + "|".join(["---"] * len(fields)) + "|\n")
    for row in rows:
        out.write("| " + " | ".join(row.get(field, "NA") for field in fields) + " |\n")
    out.write("\n## Interpretation Boundary\n\n")
    out.write("A real success requires full-denominator trans top1 accuracy to increase by at least 0.1 over the approved baseline without collapsing cis accuracy or cis distance Spearman. Copy0/copy1 are gauge labels. Eval applies the standard whole-chromosome SNP gauge policy for accuracy and uses CHARM/3DG only as post-training reference. This experiment tests blind self-training geometry only; it does not use oracle posterior labels, CHARM/3DG, or SNP truth during training.\n")
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
