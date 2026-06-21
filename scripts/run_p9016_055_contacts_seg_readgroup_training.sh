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

RUN_ID=${HK_BLIND_055_RUN_ID:-"055-$(date +%Y%m%d_%H%M%S)-p9016_contacts_seg_readgroup_training_1m"}
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE_ROOT/test_res}/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
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
mkdir -p "$OUTPUT_ROOT" "$LINK_OUTPUT_ROOT" "$EVAL_ROOT" "$LOG_ROOT" "$SOURCE_MANIFEST_ROOT" "$FULL_RUN_ROOT/scripts_snapshot"
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

APPROVED_P9016_PAIRS=/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz
EVAL_PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-$APPROVED_P9016_PAIRS}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
SOURCE_054_ROOT=${HK_BLIND_054_ROOT:-$(latest_054_root)}
SOURCE_LABEL=${HK_BLIND_055_SOURCE_LABEL:-sharec}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_RELAX_STEP:-${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-0.012}}
BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-gpu}
MAX_PARALLEL=${HK_BLIND_055_MAX_PARALLEL:-${HK_BLIND_MAX_PARALLEL:-6}}
RUN_EVAL=${HK_BLIND_P9016_RUN_EVAL:-1}
REQUIRE_EVAL=${HK_BLIND_REQUIRE_EVAL:-1}
CONFIG_SET=${HK_BLIND_055_CONFIG_SET:-full}
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
		echo "error: HK_BLIND_055_MAX_PARALLEL must be a positive integer" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_055_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi
if [ -z "$SOURCE_054_ROOT" ] || [ ! -d "$SOURCE_054_ROOT" ]; then
	echo "error: HK_BLIND_054_ROOT missing or no latest 054 root found" >&2
	exit 2
fi

HICKIT_LIKE="$SOURCE_054_ROOT/exports/$SOURCE_LABEL/hickit_like_dup100.pairs.gz"
READGROUP_ADJ="$SOURCE_054_ROOT/exports/$SOURCE_LABEL/readgroup_adjacent_maxseg10_annotated.pairs.gz"
READGROUP_ALL="$SOURCE_054_ROOT/exports/$SOURCE_LABEL/readgroup_allcomb_maxseg10_annotated.pairs.gz"
for p in "$HICKIT_LIKE" "$READGROUP_ADJ" "$READGROUP_ALL"; do
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
run_cmd python -m py_compile scripts/summarize_p9016_model_sweep.py scripts/compute_p9016_copytrack_diag.py

cp "$0" \
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
copy_source_metadata readgroup_adjacent "$READGROUP_ADJ"
copy_source_metadata readgroup_allcomb "$READGROUP_ALL"

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "source_054_root	$SOURCE_054_ROOT"
	echo "source_label	$SOURCE_LABEL"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "run_blind_p9016_minimal_sha256	$RUN_HASH"
	echo "audit_blind_p9016_full_cpu_output_sha256	$AUDIT_HASH"
	echo "approved_pairs	$APPROVED_P9016_PAIRS"
	echo "eval_pairs	$EVAL_PAIRS"
	echo "eval_tdg	$EVAL_TDG"
	echo "hickit_like_pairs	$HICKIT_LIKE"
	echo "readgroup_adjacent_pairs	$READGROUP_ADJ"
	echo "readgroup_allcomb_pairs	$READGROUP_ALL"
	echo "backend	$BACKEND"
	echo "make_gpu	$MAKE_GPU"
	echo "bin_size_bp	$BIN_SIZE"
	echo "n_iter	$N_ITER"
	echo "relax_steps	$RELAX_STEPS"
	echo "relax_step	$RELAX_STEP"
	echo "max_parallel	$MAX_PARALLEL"
	echo "config_set	$CONFIG_SET"
	echo "training_boundary	custom 054 contacts.seg-derived blind-safe raw-like pairs; phase labels and CHARM/3DG are eval-only"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SOURCE_054_ROOT=$SOURCE_054_ROOT"

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
		HK_BLIND_P9016_READGROUP_MAX_SEGMENTS=8
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
		"seg055_smoke_hickit_like|$HICKIT_LIKE|1|off|hickit_like_dup100"
		"seg055_smoke_readgroup_adj_joint|$READGROUP_ADJ|1|joint_marginal|readgroup_adjacent"
	)
else
	CONFIGS=(
		"approved_baseline_pairs_only|$APPROVED_P9016_PAIRS|0|off|approved_pairs"
		"hickit_like_dup100_pairwise|$HICKIT_LIKE|1|off|hickit_like_dup100"
		"readgroup_adjacent_pairwise|$READGROUP_ADJ|1|off|readgroup_adjacent"
		"readgroup_adjacent_joint|$READGROUP_ADJ|1|joint_marginal|readgroup_adjacent"
		"readgroup_allcomb_pairwise|$READGROUP_ALL|1|off|readgroup_allcomb"
		"readgroup_allcomb_joint|$READGROUP_ALL|1|joint_marginal|readgroup_allcomb"
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

log_command_line python - "$FULL_RUN_ROOT" "[inline 055 deltas]"
python - "$FULL_RUN_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))
by_name = {row.get("config_name", ""): row for row in rows}

def f(row, key):
    try:
        return float(row.get(key, "NA"))
    except (TypeError, ValueError):
        return None

def fmt(x):
    return "NA" if x is None else f"{x:.9g}"

base = by_name.get("approved_baseline_pairs_only")
fields = [
    "config_name","source_note","readgroup_mode","top1_all","top1_cis","top1_trans",
    "same_cross_cis","cis_spearman","entropy","pU","delta_top1_trans_vs_approved",
    "delta_top1_all_vs_approved","delta_cis_spearman_vs_approved",
]
with (root / "source_delta_summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", lineterminator="\n", fieldnames=fields)
    writer.writeheader()
    for row in rows:
        btrans = f(base, "model_top1_accuracy_genome_trans") if base else None
        ball = f(base, "model_top1_accuracy_genome_all") if base else None
        bspear = f(base, "mean_per_chrom_cis_distance_spearman") if base else None
        rtrans = f(row, "model_top1_accuracy_genome_trans")
        rall = f(row, "model_top1_accuracy_genome_all")
        rspear = f(row, "mean_per_chrom_cis_distance_spearman")
        writer.writerow({
            "config_name": row.get("config_name", "NA"),
            "source_note": row.get("source_note", "NA"),
            "readgroup_mode": row.get("readgroup_mode", "NA"),
            "top1_all": row.get("model_top1_accuracy_genome_all", "NA"),
            "top1_cis": row.get("model_top1_accuracy_genome_cis", "NA"),
            "top1_trans": row.get("model_top1_accuracy_genome_trans", "NA"),
            "same_cross_cis": row.get("model_same_cross_accuracy_genome_cis", "NA"),
            "cis_spearman": row.get("mean_per_chrom_cis_distance_spearman", "NA"),
            "entropy": row.get("final_mean_entropy", "NA"),
            "pU": row.get("final_mean_pU", "NA"),
            "delta_top1_trans_vs_approved": fmt(rtrans - btrans) if rtrans is not None and btrans is not None else "NA",
            "delta_top1_all_vs_approved": fmt(rall - ball) if rall is not None and ball is not None else "NA",
            "delta_cis_spearman_vs_approved": fmt(rspear - bspear) if rspear is not None and bspear is not None else "NA",
        })

pairs = [
    ("readgroup_adjacent_joint", "readgroup_adjacent_pairwise"),
    ("readgroup_allcomb_joint", "readgroup_allcomb_pairwise"),
]
with (root / "readgroup_matched_delta_summary.tsv").open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["joint_config","pairwise_config","delta_top1_all","delta_top1_cis","delta_top1_trans","delta_cis_spearman"])
    for joint, pairwise in pairs:
        jr = by_name.get(joint)
        pr = by_name.get(pairwise)
        writer.writerow([
            joint, pairwise,
            fmt(f(jr, "model_top1_accuracy_genome_all") - f(pr, "model_top1_accuracy_genome_all")) if jr and pr and f(jr, "model_top1_accuracy_genome_all") is not None and f(pr, "model_top1_accuracy_genome_all") is not None else "NA",
            fmt(f(jr, "model_top1_accuracy_genome_cis") - f(pr, "model_top1_accuracy_genome_cis")) if jr and pr and f(jr, "model_top1_accuracy_genome_cis") is not None and f(pr, "model_top1_accuracy_genome_cis") is not None else "NA",
            fmt(f(jr, "model_top1_accuracy_genome_trans") - f(pr, "model_top1_accuracy_genome_trans")) if jr and pr and f(jr, "model_top1_accuracy_genome_trans") is not None and f(pr, "model_top1_accuracy_genome_trans") is not None else "NA",
            fmt(f(jr, "mean_per_chrom_cis_distance_spearman") - f(pr, "mean_per_chrom_cis_distance_spearman")) if jr and pr and f(jr, "mean_per_chrom_cis_distance_spearman") is not None and f(pr, "mean_per_chrom_cis_distance_spearman") is not None else "NA",
        ])
PY

log_command_line python - "$FULL_RUN_ROOT" "[inline README 055]"
python - "$FULL_RUN_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = list(csv.DictReader((root / "source_delta_summary.tsv").open(), delimiter="\t"))
lines = [
    "# 055 P9016 contacts.seg Readgroup Training 1Mb",
    "",
    "This experiment trains blind diploid reconstructions on 054 contacts.seg-derived blind-safe pairs and evaluates with the standard P9016 eval wrapper.",
    "",
    "## Paths",
    "",
    f"- full result root: `{root}`",
    f"- summary: `{root / 'summary.tsv'}`",
    f"- source delta summary: `{root / 'source_delta_summary.tsv'}`",
    f"- readgroup matched delta summary: `{root / 'readgroup_matched_delta_summary.tsv'}`",
    "",
    "## Boundary",
    "",
    "- Training uses approved P9016 pairs or 054 contacts.seg-derived raw-like pairs with phase labels dropped.",
    "- SNP phase and CHARM/3DG are used only by post-training eval.",
    "- copy0/copy1 remain gauge labels.",
    "",
    "## Main Results",
    "",
    "| config | source | readgroup | top1 all | top1 cis | top1 trans | cis Spearman | delta trans vs approved |",
    "|---|---|---|---:|---:|---:|---:|---:|",
]
for row in rows:
    lines.append(
        "| {config} | {source} | {rg} | {all} | {cis} | {trans} | {spear} | {dtrans} |".format(
            config=row.get("config_name", "NA"),
            source=row.get("source_note", "NA"),
            rg=row.get("readgroup_mode", "NA"),
            all=row.get("top1_all", "NA"),
            cis=row.get("top1_cis", "NA"),
            trans=row.get("top1_trans", "NA"),
            spear=row.get("cis_spearman", "NA"),
            dtrans=row.get("delta_top1_trans_vs_approved", "NA"),
        )
    )
lines.extend([
    "",
    "## Interpretation",
    "",
    "Read this as a source/interface experiment, not a solved model claim. hickit_like_dup100 tests hickit-style upstream contacts; readgroup_adjacent/allcomb test whether numeric readgroup metadata helps beyond matched pairwise contacts.",
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
