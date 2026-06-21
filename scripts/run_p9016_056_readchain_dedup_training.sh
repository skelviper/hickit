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

RUN_ID=${HK_BLIND_056_RUN_ID:-"056-$(date +%Y%m%d_%H%M%S)-p9016_readchain_dedup_training_1m"}
EXPECTED_TEST_RES_ROOT="$PHASE_ROOT/test_res"
REQUESTED_TEST_RES_ROOT="${HK_BLIND_TEST_RES_ROOT:-$EXPECTED_TEST_RES_ROOT}"
mkdir -p "$EXPECTED_TEST_RES_ROOT"
EXPECTED_TEST_RES_ROOT=$(cd "$EXPECTED_TEST_RES_ROOT" && pwd)
if [ "$REQUESTED_TEST_RES_ROOT" != "$EXPECTED_TEST_RES_ROOT" ]; then
	echo "error: 056 full outputs must be written under $EXPECTED_TEST_RES_ROOT" >&2
	echo "error: got HK_BLIND_TEST_RES_ROOT=$REQUESTED_TEST_RES_ROOT" >&2
	exit 2
fi
FULL_RUN_ROOT="$EXPECTED_TEST_RES_ROOT/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
EXPORT_ROOT="$FULL_RUN_ROOT/exports/sharec"
OUTPUT_ROOT="$FULL_RUN_ROOT/work_outputs"
LINK_OUTPUT_ROOT="$FULL_RUN_ROOT/outputs"
EVAL_ROOT="$FULL_RUN_ROOT/eval"
LOG_ROOT="$FULL_RUN_ROOT/logs"
SOURCE_MANIFEST_ROOT="$FULL_RUN_ROOT/source_pairs_manifest"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"

if [ -e "$FULL_RUN_ROOT" ] &&
   [ "${HK_BLIND_ALLOW_NONEMPTY_RUN_ROOT:-0}" != "1" ] &&
   find "$FULL_RUN_ROOT" -mindepth 1 -print -quit | grep -q .; then
	echo "error: run root already exists and is not empty: $FULL_RUN_ROOT" >&2
	exit 1
fi
mkdir -p "$EXPORT_ROOT" "$OUTPUT_ROOT" "$LINK_OUTPUT_ROOT" "$EVAL_ROOT" "$LOG_ROOT" "$SOURCE_MANIFEST_ROOT" "$FULL_RUN_ROOT/scripts_snapshot"
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)

exec > >(tee -a "$LOG_ROOT/runner.stdout.log") 2> >(tee -a "$LOG_ROOT/runner.stderr.log" >&2)

log_command_line() {
	local line="+"
	local arg
	for arg in "$@"; do
		line+=" $(printf '%q' "$arg")"
	done
	printf '%s\n' "$line" >> "$COMMANDS_LOG"
}

run_cmd() {
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

latest_054_root() {
	find "$PHASE_ROOT/test_res" -maxdepth 1 -type d -name '054-*-p9016_contacts_seg_readgroup_dedup_audit_1m' -printf '%T@ %p\n' 2>/dev/null |
		sort -nr | awk 'NR==1{print $2}'
}

SEG=${HK_BLIND_P9016_CONTACTS_SEG_SHAREC:-/sharec/zliu/CHARM/mESC/processed/P9016/2d_info/contacts.seg.gz}
SOURCE_054_ROOT=${HK_BLIND_054_ROOT:-$(latest_054_root)}
APPROVED_P9016_PAIRS=/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz
EVAL_PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-$APPROVED_P9016_PAIRS}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_RELAX_STEP:-${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-0.012}}
BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-gpu}
MAX_PARALLEL=${HK_BLIND_056_MAX_PARALLEL:-${HK_BLIND_MAX_PARALLEL:-6}}
RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-1}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-1}
CONFIG_SET=${HK_BLIND_056_CONFIG_SET:-full}
MIN_MAPQ=${HK_BLIND_056_MIN_MAPQ:-20}
MIN_LEG_DIST=${HK_BLIND_056_MIN_LEG_DIST:-1000}
MAX_SEG=${HK_BLIND_056_MAX_SEG:-8}
DUP_DIST=${HK_BLIND_056_DUP_DIST:-100}
MAX_SEG_SLUG="maxseg${MAX_SEG}"
DUP_DIST_SLUG="dup${DUP_DIST}"
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
		echo "error: HK_BLIND_056_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_056_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi
if [ ! -s "$SEG" ]; then
	echo "error: contacts.seg missing: $SEG" >&2
	exit 2
fi
if [ -z "$SOURCE_054_ROOT" ] || [ ! -d "$SOURCE_054_ROOT" ]; then
	echo "error: HK_BLIND_054_ROOT missing or no latest 054 root found" >&2
	exit 2
fi

HICKIT_LIKE="$SOURCE_054_ROOT/exports/sharec/hickit_like_dup100.pairs.gz"
OLD_READGROUP_ADJ="$SOURCE_054_ROOT/exports/sharec/readgroup_adjacent_maxseg10_annotated.pairs.gz"
for p in "$HICKIT_LIKE" "$OLD_READGROUP_ADJ"; do
	if [ ! -s "$p" ]; then
		echo "error: source pairs missing: $p" >&2
		exit 2
	fi
done

GIT_COMMIT=$(git rev-parse HEAD 2>/dev/null || echo NA)
GIT_DIRTY_COUNT=$(git status --porcelain | wc -l | awk '{print $1}')
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

if [ "$MAKE_GPU" = "1" ]; then
	run_cmd make gpu=1 run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
else
	run_cmd make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
fi
run_cmd python -m py_compile scripts/export_p9016_contacts_seg_readgroup_dedup.py scripts/summarize_p9016_model_sweep.py scripts/compute_p9016_copytrack_diag.py

cp "$0" \
   scripts/export_p9016_contacts_seg_readgroup_dedup.py \
   scripts/p9016_common_eval.sh \
   scripts/p9016_publish_light_result.sh \
   scripts/summarize_p9016_model_sweep.py \
   scripts/compute_p9016_copytrack_diag.py \
   "$FULL_RUN_ROOT/scripts_snapshot/"

RUN_BIN="$LOG_ROOT/run_blind_p9016_minimal.bin"
AUDIT_BIN="$LOG_ROOT/audit_blind_p9016_full_cpu_output.bin"
cp run_blind_p9016_minimal.bin "$RUN_BIN"
cp audit_blind_p9016_full_cpu_output.bin "$AUDIT_BIN"
chmod +x "$RUN_BIN" "$AUDIT_BIN"
RUN_HASH=$(hash_file "$RUN_BIN")
AUDIT_HASH=$(hash_file "$AUDIT_BIN")
{
	echo "run_blind_p9016_minimal_sha256	$RUN_HASH"
	echo "audit_blind_p9016_full_cpu_output_sha256	$AUDIT_HASH"
} >> "$LOG_ROOT/build_env.txt"

export_readchain_both() {
	local name=$1
	run_cmd python scripts/export_p9016_contacts_seg_readgroup_dedup.py \
		--seg "$SEG" \
		--source-label sharec \
		--outdir "$EXPORT_ROOT" \
		--export-name "$name" \
		--pair-mode both \
		--dedup-unit readchain \
		--coordinate-mode boundary \
		--readgroup-cols \
		--min-mapq "$MIN_MAPQ" \
		--min-leg-dist "$MIN_LEG_DIST" \
		--max-seg "$MAX_SEG" \
		--dup-dist "$DUP_DIST"
}
READCHAIN_BASE="readchain_${DUP_DIST_SLUG}_${MAX_SEG_SLUG}"
export_readchain_both "$READCHAIN_BASE"

READCHAIN_ADJ="$EXPORT_ROOT/${READCHAIN_BASE}_adjacent.pairs.gz"
READCHAIN_ALL="$EXPORT_ROOT/${READCHAIN_BASE}_allcomb.pairs.gz"

copy_source_metadata() {
	local name=$1
	local pairs=$2
	local meta="${pairs%.pairs.gz}.metadata.tsv"
	local out="$SOURCE_MANIFEST_ROOT/$name.source.tsv"
	{
		echo "key	value"
		echo "name	$name"
		echo "pairs	$pairs"
		echo "pairs_size_bytes	$(stat -c '%s' "$pairs")"
		echo "pairs_sha256	$(hash_file "$pairs")"
		echo "source_metadata	$meta"
		if [ -s "$meta" ]; then
			tail -n +2 "$meta"
		fi
	} > "$out"
}
copy_source_metadata hickit_like_dup100 "$HICKIT_LIKE"
copy_source_metadata old_readgroup_adjacent "$OLD_READGROUP_ADJ"
copy_source_metadata readchain_adjacent "$READCHAIN_ADJ"
copy_source_metadata readchain_allcomb "$READCHAIN_ALL"

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "source_054_root	$SOURCE_054_ROOT"
	echo "contacts_seg	$SEG"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "run_blind_p9016_minimal_sha256	$RUN_HASH"
	echo "audit_blind_p9016_full_cpu_output_sha256	$AUDIT_HASH"
	echo "approved_pairs	$APPROVED_P9016_PAIRS"
	echo "eval_pairs	$EVAL_PAIRS"
	echo "eval_tdg	$EVAL_TDG"
	echo "hickit_like_pairs	$HICKIT_LIKE"
	echo "old_readgroup_adjacent_pairs	$OLD_READGROUP_ADJ"
	echo "readchain_adjacent_pairs	$READCHAIN_ADJ"
	echo "readchain_allcomb_pairs	$READCHAIN_ALL"
	echo "backend	$BACKEND"
	echo "make_gpu	$MAKE_GPU"
	echo "bin_size_bp	$BIN_SIZE"
	echo "n_iter	$N_ITER"
	echo "relax_steps	$RELAX_STEPS"
	echo "relax_step	$RELAX_STEP"
	echo "max_parallel	$MAX_PARALLEL"
	echo "config_set	$CONFIG_SET"
	echo "min_mapq	$MIN_MAPQ"
	echo "min_leg_dist	$MIN_LEG_DIST"
	echo "max_seg	$MAX_SEG"
	echo "dup_dist	$DUP_DIST"
	echo "training_boundary	custom contacts.seg-derived readchain-dedup pairs; phase labels and CHARM/3DG are eval-only"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "READCHAIN_ADJ=$READCHAIN_ADJ"
echo "READCHAIN_ALL=$READCHAIN_ALL"

run_config() {
	local config=$1
	local pairs=$2
	local allow_custom=$3
	local readgroup_mode=$4
	local source_note=$5
	local input_contact_source="contacts_seg_derived_pairs"
	local allow_nonstandard=1
	local out_root="$OUTPUT_ROOT/$config"
	local eval_dir="$EVAL_ROOT/$config"
	local log_prefix="$LOG_ROOT/$config"
	local manifest
	if [ "$source_note" = "approved_pairs" ]; then
		input_contact_source="raw_pairs"
		allow_nonstandard=0
	fi

	echo "CONFIG_START=$config pairs=$pairs readgroup=$readgroup_mode"
	rm -rf "$out_root" "$LINK_OUTPUT_ROOT/$config" "$eval_dir"
	mkdir -p "$eval_dir"

	local train_cmd=(
		clean_env_prefix
		HK_BLIND_SAMPLE=P9016
		HK_BLIND_P9016_SAMPLE=P9016
		HK_BLIND_GIT_COMMIT="$GIT_COMMIT"
		HK_BLIND_GIT_DIRTY_COUNT="$GIT_DIRTY_COUNT"
		HK_BLIND_BINARY_HASH="$RUN_HASH"
		HK_BLIND_P9016_PAIRS="$pairs"
		HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS="$allow_custom"
		HK_BLIND_P9016_INPUT_CONTACT_SOURCE="$input_contact_source"
		HK_BLIND_P9016_ALLOW_NONSTANDARD_PAIRS="$allow_nonstandard"
		HK_BLIND_P9016_OUTPUT_ROOT="$out_root"
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
		HK_BLIND_P9016_READGROUP_MODE="$readgroup_mode"
		HK_BLIND_P9016_READGROUP_MAX_SEGMENTS="$MAX_SEG"
		HK_BLIND_P9016_READGROUP_EPS=1e-6
		HK_BLIND_WRITE_RAW=0
		"$RUN_BIN"
	)
	log_command_line "${train_cmd[@]}"
	"${train_cmd[@]}" > "$log_prefix.train.log" 2>&1

	ln -s "../work_outputs/$config/$config" "$LINK_OUTPUT_ROOT/$config"

	local audit_cmd=(
		clean_env_prefix
		HK_BLIND_AUDIT_ALLOW_CUSTOM_RAW_PAIRS=1
		HK_BLIND_AUDIT_ALLOW_NONSTANDARD_PAIRS=1
		"$AUDIT_BIN" "$LINK_OUTPUT_ROOT/$config"
	)
	log_command_line "${audit_cmd[@]}"
	"${audit_cmd[@]}" > "$log_prefix.audit.log" 2>&1

	manifest="$LINK_OUTPUT_ROOT/$config/p9016_full.manifest.tsv"
	{
		echo "source_note	$source_note"
		echo "source_pairs	$pairs"
		echo "source_pairs_sha256	$(hash_file "$pairs")"
	} >> "$manifest"

	if [ "$RUN_EVAL" = "1" ]; then
		local eval_cmd=(
			clean_env_prefix
			HK_BLIND_SAMPLE=P9016
			HK_BLIND_P9016_SAMPLE=P9016
			HK_BLIND_P9016_EVAL_PAIRS="$EVAL_PAIRS"
			HK_BLIND_P9016_EVAL_TDG="$EVAL_TDG"
			HK_BLIND_P9016_BIN_SIZE_BP="$BIN_SIZE"
			HK_BLIND_REQUIRE_EVAL="$REQUIRE_EVAL"
			scripts/p9016_common_eval.sh "$LINK_OUTPUT_ROOT/$config" "$eval_dir"
		)
		log_command_line "${eval_cmd[@]}"
		"${eval_cmd[@]}" > "$log_prefix.eval.log" 2>&1
	else
		printf 'status\tSKIPPED_EVAL_DISABLED\n' > "$eval_dir/eval_status.tsv"
	fi

	if [ -s "$LINK_OUTPUT_ROOT/$config/p9016_full.coords.tsv" ]; then
		local copytrack_cmd=(
			python scripts/compute_p9016_copytrack_diag.py
			--config-output-dir "$LINK_OUTPUT_ROOT/$config"
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
		"seg056_smoke_readchain_adj_joint|$READCHAIN_ADJ|1|joint_marginal|readchain_adjacent"
	)
else
	CONFIGS=(
		"approved_baseline_pairs_only|$APPROVED_P9016_PAIRS|0|off|approved_pairs"
		"hickit_like_dup100_pairwise|$HICKIT_LIKE|1|off|hickit_like_dup100"
		"old_readgroup_adjacent_joint|$OLD_READGROUP_ADJ|1|joint_marginal|old_readgroup_adjacent"
		"readchain_adjacent_pairwise|$READCHAIN_ADJ|1|off|readchain_adjacent"
		"readchain_adjacent_joint|$READCHAIN_ADJ|1|joint_marginal|readchain_adjacent"
		"readchain_allcomb_pairwise|$READCHAIN_ALL|1|off|readchain_allcomb"
		"readchain_allcomb_joint|$READCHAIN_ALL|1|joint_marginal|readchain_allcomb"
	)
fi

failures=0
running=0
for spec in "${CONFIGS[@]}"; do
	IFS='|' read -r config pairs allow_custom readgroup source_note <<< "$spec"
	run_config "$config" "$pairs" "$allow_custom" "$readgroup" "$source_note" &
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

log_command_line python - "$FULL_RUN_ROOT" "[inline 056 deltas]"
python - "$FULL_RUN_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))
by_name = {row.get("config_name", ""): row for row in rows}

def f(row, key):
    if not row:
        return None
    try:
        v = row.get(key, "NA")
        if v in (None, "", "NA", "nan", "NaN"):
            return None
        return float(v)
    except (TypeError, ValueError):
        return None

def fmt(x):
    return "NA" if x is None else f"{x:.9g}"

base = by_name.get("approved_baseline_pairs_only")
fields = [
    "config_name","source_note","source_pairs","readgroup_mode","backend",
    "n_raw","n_bpair","n_raw_trans","n_bpair_trans",
    "readgroup_entries","readgroup_groups_used","readgroup_groups_skipped_too_large","readgroup_groups_skipped_bad",
    "top1_all","top1_cis","top1_trans","same_cross_cis","same_cross_trans","cis_spearman",
    "entropy","pU","final_min_sep","final_mean_sep","delta_top1_trans_vs_approved","delta_top1_all_vs_approved","delta_cis_spearman_vs_approved",
]
with (root / "source_delta_summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", lineterminator="\n", fieldnames=fields)
    writer.writeheader()
    for row in rows:
        btrans = f(base, "model_top1_accuracy_genome_trans")
        ball = f(base, "model_top1_accuracy_genome_all")
        bspear = f(base, "mean_per_chrom_cis_distance_spearman")
        rtrans = f(row, "model_top1_accuracy_genome_trans")
        rall = f(row, "model_top1_accuracy_genome_all")
        rspear = f(row, "mean_per_chrom_cis_distance_spearman")
        writer.writerow({
            "config_name": row.get("config_name", "NA"),
            "source_note": row.get("source_note", "NA"),
            "source_pairs": row.get("source_pairs", "NA"),
            "readgroup_mode": row.get("readgroup_mode", "NA"),
            "backend": row.get("backend", "NA"),
            "n_raw": row.get("n_raw", "NA"),
            "n_bpair": row.get("n_bpair", "NA"),
            "n_raw_trans": row.get("n_raw_trans", "NA"),
            "n_bpair_trans": row.get("n_bpair_trans", "NA"),
            "readgroup_entries": row.get("readgroup_entries", "NA"),
            "readgroup_groups_used": row.get("readgroup_groups_used", "NA"),
            "readgroup_groups_skipped_too_large": row.get("readgroup_groups_skipped_too_large", "NA"),
            "readgroup_groups_skipped_bad": row.get("readgroup_groups_skipped_bad", "NA"),
            "top1_all": row.get("model_top1_accuracy_genome_all", "NA"),
            "top1_cis": row.get("model_top1_accuracy_genome_cis", "NA"),
            "top1_trans": row.get("model_top1_accuracy_genome_trans", "NA"),
            "same_cross_cis": row.get("model_same_cross_accuracy_genome_cis", "NA"),
            "same_cross_trans": row.get("model_same_cross_accuracy_genome_trans", "NA"),
            "cis_spearman": row.get("mean_per_chrom_cis_distance_spearman", "NA"),
            "entropy": row.get("final_mean_entropy", "NA"),
            "pU": row.get("final_mean_pU", "NA"),
            "final_min_sep": row.get("final_min_sep", "NA"),
            "final_mean_sep": row.get("final_mean_sep", "NA"),
            "delta_top1_trans_vs_approved": fmt(rtrans - btrans) if rtrans is not None and btrans is not None else "NA",
            "delta_top1_all_vs_approved": fmt(rall - ball) if rall is not None and ball is not None else "NA",
            "delta_cis_spearman_vs_approved": fmt(rspear - bspear) if rspear is not None and bspear is not None else "NA",
        })

pairs = [
    ("readchain_adjacent_joint", "readchain_adjacent_pairwise"),
    ("readchain_allcomb_joint", "readchain_allcomb_pairwise"),
]
with (root / "readgroup_matched_delta_summary.tsv").open("w", newline="") as fh:
    fields2 = ["joint_config","pairwise_config","delta_top1_all","delta_top1_cis","delta_top1_trans","delta_same_cross_cis","delta_same_cross_trans","delta_cis_spearman","delta_entropy","delta_pU"]
    writer = csv.DictWriter(fh, delimiter="\t", lineterminator="\n", fieldnames=fields2)
    writer.writeheader()
    for joint, pairwise in pairs:
        jr = by_name.get(joint)
        pr = by_name.get(pairwise)
        def d(key):
            a = f(jr, key); b = f(pr, key)
            return fmt(a - b) if a is not None and b is not None else "NA"
        writer.writerow({
            "joint_config": joint,
            "pairwise_config": pairwise,
            "delta_top1_all": d("model_top1_accuracy_genome_all"),
            "delta_top1_cis": d("model_top1_accuracy_genome_cis"),
            "delta_top1_trans": d("model_top1_accuracy_genome_trans"),
            "delta_same_cross_cis": d("model_same_cross_accuracy_genome_cis"),
            "delta_same_cross_trans": d("model_same_cross_accuracy_genome_trans"),
            "delta_cis_spearman": d("mean_per_chrom_cis_distance_spearman"),
            "delta_entropy": d("final_mean_entropy"),
            "delta_pU": d("final_mean_pU"),
        })

copy_fields = ["config_name","source_note","readgroup_mode","backend","copytrack_sep_p05","copytrack_sep_median","copytrack_vector_cos_p05","copytrack_vector_cos_median","copytrack_frac_cos_lt_0","copytrack_frac_projection_sign_switch"]
with (root / "copytrack_summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", lineterminator="\n", fieldnames=copy_fields)
    writer.writeheader()
    for row in rows:
        writer.writerow({k: row.get(k, "NA") for k in copy_fields})

source_rows = list(csv.DictReader((root / "source_delta_summary.tsv").open(), delimiter="\t"))
matched_rows = list(csv.DictReader((root / "readgroup_matched_delta_summary.tsv").open(), delimiter="\t"))
lines = [
    "# 056 P9016 Readchain-Dedup Readgroup Training 1Mb",
    "",
    "This experiment tests whole-read/segment-chain deduplication before pair generation. A duplicate chain such as A1-B1-C1 and A2-B2-C2 is collapsed into one representative A-B-C, then adjacent or all-combination contacts are emitted with readgroup metadata.",
    "",
    "## Paths",
    "",
    f"- full result root: `{root}`",
    f"- light result root: `{Path('/mnt/ssd/zliu/phase3/hickit/result') / root.name}`",
    f"- summary: `{root / 'summary.tsv'}`",
    f"- source delta summary: `{root / 'source_delta_summary.tsv'}`",
    f"- readgroup matched delta summary: `{root / 'readgroup_matched_delta_summary.tsv'}`",
    "",
    "## Boundary",
    "",
    "- `approved_baseline_pairs_only` is the only strict approved `/shared/.../P9016.pairs.gz` baseline row.",
    "- All readchain rows are contacts.seg-derived custom-pairs diagnostics. Phase labels are dropped before training.",
    "- SNP phase and CHARM/3DG are used only by post-training eval through `scripts/p9016_common_eval.sh`.",
    "- copy0/copy1 remain gauge labels.",
    "",
    "## Main Results",
    "",
    "| config | source | readgroup | n raw | trans raw | top1 all | top1 cis | top1 trans | cis Spearman | entropy | pU | delta trans vs approved |",
    "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
]
for row in source_rows:
    lines.append("| {config} | {source} | {rg} | {nraw} | {ntr} | {all} | {cis} | {trans} | {spear} | {entropy} | {pU} | {dtrans} |".format(
        config=row.get("config_name", "NA"), source=row.get("source_note", "NA"), rg=row.get("readgroup_mode", "NA"),
        nraw=row.get("n_raw", "NA"), ntr=row.get("n_raw_trans", "NA"), all=row.get("top1_all", "NA"),
        cis=row.get("top1_cis", "NA"), trans=row.get("top1_trans", "NA"), spear=row.get("cis_spearman", "NA"),
        entropy=row.get("entropy", "NA"), pU=row.get("pU", "NA"), dtrans=row.get("delta_top1_trans_vs_approved", "NA")))
lines.extend([
    "",
    "## Matched Readgroup Deltas",
    "",
    "| joint config | pairwise config | delta top1 all | delta top1 cis | delta top1 trans | delta cis Spearman |",
    "|---|---|---:|---:|---:|---:|",
])
for row in matched_rows:
    lines.append("| {joint} | {pairwise} | {dall} | {dcis} | {dtrans} | {dspear} |".format(
        joint=row.get("joint_config", "NA"), pairwise=row.get("pairwise_config", "NA"),
        dall=row.get("delta_top1_all", "NA"), dcis=row.get("delta_top1_cis", "NA"),
        dtrans=row.get("delta_top1_trans", "NA"), dspear=row.get("delta_cis_spearman", "NA")))
lines.extend([
    "",
    "## Interpretation",
    "",
    "This is a controlled source/interface test. A positive readgroup claim requires the joint row to beat the matched same-source pairwise row without damaging cis accuracy or cis distance Spearman, and ideally to approach or exceed the approved baseline.",
])
(root / "README.md").write_text("\n".join(lines) + "\n")
PY

log_command_line scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$LOG_ROOT/publish.log" 2>&1

echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "SOURCE_DELTA_TSV=$FULL_RUN_ROOT/source_delta_summary.tsv"
echo "READGROUP_MATCHED_DELTA_TSV=$FULL_RUN_ROOT/readgroup_matched_delta_summary.tsv"
echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
