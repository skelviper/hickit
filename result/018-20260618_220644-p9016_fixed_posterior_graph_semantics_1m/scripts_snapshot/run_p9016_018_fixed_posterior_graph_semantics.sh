#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PHASE3_ROOT=$(cd "$REPO_ROOT/.." && pwd)

if [ "${CONDA_DEFAULT_ENV:-}" != "analysis" ]; then
	echo "This repository must be run inside Conda env 'analysis'." >&2
	echo "Run: conda activate analysis" >&2
	exit 2
fi

RUN_ID=${HK_BLIND_018_RUN_ID:-018-$(date +%Y%m%d_%H%M%S)-p9016_fixed_posterior_graph_semantics_1m}
FULL_RUN_ROOT=${HK_BLIND_TEST_RES_ROOT:-"$PHASE3_ROOT/test_res"}/$RUN_ID
LIGHT_RESULT_ROOT=$REPO_ROOT/result/$RUN_ID
COMMANDS_LOG=$FULL_RUN_ROOT/commands.log

SOURCE_DIR=${HK_FIXED_POSTERIOR_SOURCE_DIR:-$PHASE3_ROOT/test_res/016-20260618_150851-p9016_charm_contacts_condition_sweep_1m/work_outputs/1pr/p9016_charm3dg20k_1pr_pcgamma1_common_msep1p5_lsep1/p9016_charm3dg20k_1pr_pcgamma1_common_msep1p5_lsep1}
SOURCE_EVAL_DIR=${HK_FIXED_POSTERIOR_SOURCE_EVAL_DIR:-$PHASE3_ROOT/test_res/016-20260618_150851-p9016_charm_contacts_condition_sweep_1m/eval/1pr/p9016_charm3dg20k_1pr_pcgamma1_common_msep1p5_lsep1}

mkdir -p "$FULL_RUN_ROOT"/{logs,outputs,eval,diagnostics,scripts_snapshot}

log_cmd() {
	printf '%s\n' "$*" >> "$COMMANDS_LOG"
}

run_logged() {
	log_cmd "$*"
	"$@"
}

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SOURCE_DIR=$SOURCE_DIR"

if [ ! -s "$SOURCE_DIR/p9016_full.coords.tsv" ] || [ ! -s "$SOURCE_DIR/p9016_full.bpair_posterior.tsv" ]; then
	echo "missing fixed posterior source outputs: $SOURCE_DIR" >&2
	exit 2
fi

cd "$REPO_ROOT"
git rev-parse HEAD > "$FULL_RUN_ROOT/logs/git_commit.txt"
git status --porcelain > "$FULL_RUN_ROOT/logs/git_status.txt"
git diff --stat > "$FULL_RUN_ROOT/logs/git_diff_stat.txt" || true
git diff > "$FULL_RUN_ROOT/logs/git_diff.patch" || true
{
	echo "timestamp	$(date -Is)"
	echo "repo_root	$REPO_ROOT"
	echo "phase3_root	$PHASE3_ROOT"
	echo "run_id	$RUN_ID"
	echo "source_dir	$SOURCE_DIR"
	echo "source_eval_dir	$SOURCE_EVAL_DIR"
	echo "cc	${CC:-cc}"
	echo "cflags	${CFLAGS:-default}"
} > "$FULL_RUN_ROOT/logs/build_env.txt"

cp "$0" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/summarize_p9016_018_fixed_posterior.py" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/diagnose_p9016_trans_gauge_geometry.py" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/diagnose_p9016_force_regime.py" "$FULL_RUN_ROOT/scripts_snapshot/"

run_logged make run_blind_p9016_fixed_posterior_relax.bin
BINARY_HASH=$(sha256sum run_blind_p9016_fixed_posterior_relax.bin | awk '{print $1}')
GIT_COMMIT=$(git rev-parse HEAD)
GIT_DIRTY_COUNT=$(git status --porcelain | wc -l | awk '{print $1}')

cat > "$FULL_RUN_ROOT/run_manifest.tsv" <<EOF
key	value
run_id	$RUN_ID
full_run_root	$FULL_RUN_ROOT
light_result_root	$LIGHT_RESULT_ROOT
source_dir	$SOURCE_DIR
source_eval_dir	$SOURCE_EVAL_DIR
git_commit	$GIT_COMMIT
git_dirty_count	$GIT_DIRTY_COUNT
binary_hash	$BINARY_HASH
experiment_type	fixed_posterior_graph_semantics
training_boundary	reference_derived_positive_control_not_blind
EOF

run_config() {
	local config=$1
	local graph_mode=$2
	local trans_mult=$3
	local min_sep=$4
	local lambda_sep=$5
	local out_dir="$FULL_RUN_ROOT/outputs/$config"
	local eval_dir="$FULL_RUN_ROOT/eval/$config"
	local diag_dir="$FULL_RUN_ROOT/diagnostics/$config"
	mkdir -p "$out_dir" "$eval_dir" "$diag_dir"
	echo "[018] running $config"
	log_cmd "HK_FIXED_POSTERIOR_SOURCE_DIR=$SOURCE_DIR HK_FIXED_POSTERIOR_OUTPUT_DIR=$out_dir HK_FIXED_POSTERIOR_CONFIG_NAME=$config HK_FIXED_POSTERIOR_GRAPH_MODE=$graph_mode HK_FIXED_POSTERIOR_TRANS_DSCALE_MULTIPLIER=$trans_mult HK_FIXED_POSTERIOR_MIN_SEP_UNIT=$min_sep HK_FIXED_POSTERIOR_LAMBDA_SEP=$lambda_sep HK_FIXED_POSTERIOR_REPULSION_MULTIPLIER=${HK_FIXED_POSTERIOR_REPULSION_MULTIPLIER:-1} ./run_blind_p9016_fixed_posterior_relax.bin"
	HK_BLIND_GIT_COMMIT="$GIT_COMMIT" \
	HK_BLIND_GIT_DIRTY_COUNT="$GIT_DIRTY_COUNT" \
	HK_BLIND_BINARY_HASH="$BINARY_HASH" \
	HK_BLIND_SAMPLE=P9016 \
	HK_FIXED_POSTERIOR_SOURCE_DIR="$SOURCE_DIR" \
	HK_FIXED_POSTERIOR_OUTPUT_DIR="$out_dir" \
	HK_FIXED_POSTERIOR_CONFIG_NAME="$config" \
	HK_FIXED_POSTERIOR_GRAPH_MODE="$graph_mode" \
	HK_FIXED_POSTERIOR_D_SCALE_MODE=posterior_count \
	HK_FIXED_POSTERIOR_D_SCALE_POSTERIOR_GAMMA=1 \
	HK_FIXED_POSTERIOR_D_SCALE_EPS_COUNT=1e-6 \
	HK_FIXED_POSTERIOR_TRANS_DSCALE_MULTIPLIER="$trans_mult" \
	HK_FIXED_POSTERIOR_MIN_SEP_UNIT="$min_sep" \
	HK_FIXED_POSTERIOR_LAMBDA_SEP="$lambda_sep" \
	HK_FIXED_POSTERIOR_REPULSION_MULTIPLIER=${HK_FIXED_POSTERIOR_REPULSION_MULTIPLIER:-1} \
	HK_FIXED_POSTERIOR_RELAX_STEP=${HK_FIXED_POSTERIOR_RELAX_STEP:-0.012} \
	HK_FIXED_POSTERIOR_RELAX_STEPS=${HK_FIXED_POSTERIOR_RELAX_STEPS:-100} \
	./run_blind_p9016_fixed_posterior_relax.bin > "$FULL_RUN_ROOT/logs/$config.train.log" 2>&1

	log_cmd "HK_BLIND_ALLOW_REFERENCE_DERIVED_TRAINING=1 scripts/p9016_common_eval.sh $out_dir $eval_dir"
	HK_BLIND_ALLOW_REFERENCE_DERIVED_TRAINING=1 \
	HK_BLIND_SAMPLE=P9016 \
	HK_BLIND_P9016_SAMPLE=P9016 \
	"$SCRIPT_DIR/p9016_common_eval.sh" "$out_dir" "$eval_dir"

	log_cmd "python scripts/diagnose_p9016_trans_gauge_geometry.py --config-output-dir $out_dir --eval-dir $eval_dir --outdir $diag_dir"
	python "$SCRIPT_DIR/diagnose_p9016_trans_gauge_geometry.py" \
		--config-output-dir "$out_dir" \
		--eval-dir "$eval_dir" \
		--outdir "$diag_dir" > "$FULL_RUN_ROOT/logs/$config.trans_gauge.log" 2>&1

	log_cmd "python scripts/diagnose_p9016_force_regime.py --config-output-dir $out_dir --outdir $diag_dir"
	python "$SCRIPT_DIR/diagnose_p9016_force_regime.py" \
		--config-output-dir "$out_dir" \
		--outdir "$diag_dir" > "$FULL_RUN_ROOT/logs/$config.force_regime.log" 2>&1
}

run_config p9016_fixed_state_softall softall 1 1.5 1
run_config p9016_fixed_state_top1_only top1_only 1 1.5 1
run_config p9016_fixed_state_top2_only top2_only 1 1.5 1
run_config p9016_fixed_state_same_cross_gated same_cross_gated 1 1.5 1
run_config p9016_fixed_state_softall_trans_dscale0p75 softall 0.75 1.5 1
run_config p9016_fixed_state_softall_trans_dscale0p5 softall 0.5 1.5 1

HK_FIXED_POSTERIOR_SOURCE_DIR="$SOURCE_DIR" \
HK_FIXED_POSTERIOR_SOURCE_EVAL_DIR="$SOURCE_EVAL_DIR" \
python "$SCRIPT_DIR/summarize_p9016_018_fixed_posterior.py" "$FULL_RUN_ROOT" > "$FULL_RUN_ROOT/summary.tsv"

python - <<'PY' "$FULL_RUN_ROOT"
from pathlib import Path
import csv
import sys

root = Path(sys.argv[1])
summary = root / "summary.tsv"
rows = list(csv.DictReader(summary.open(), delimiter="\t")) if summary.exists() else []
manifest = {}
with (root / "run_manifest.tsv").open() as fh:
    reader = csv.reader(fh, delimiter="\t")
    next(reader, None)
    for row in reader:
        if len(row) >= 2:
            manifest[row[0]] = row[1]

def val(row, key):
    v = row.get(key, "NA")
    return v if v not in ("", "NA") else "NA"

best_nearest = max(rows, key=lambda r: float(r.get("nearest_geometry_state_accuracy", "nan")) if r.get("nearest_geometry_state_accuracy", "NA") != "NA" else float("-inf"), default=None)
best_spearman = max(rows, key=lambda r: float(r.get("mean_per_chrom_cis_distance_spearman", "nan")) if r.get("mean_per_chrom_cis_distance_spearman", "NA") != "NA" else float("-inf"), default=None)

lines = [
    "# 018 P9016 Fixed-Posterior Graph Semantics",
    "",
    "This is a reference-derived positive-control diagnostic, not a blind baseline. It starts from the best 016 CHARM/3DG-derived 1-pr condition, fixes its binned posterior and coordinates, then tests direct fixed-state-edge M-step variants for one relaxation.",
    "",
    f"- full result root: `{root}`",
    f"- source fixed posterior: `{manifest.get('source_dir', 'NA')}`",
    f"- summary: `{summary}`",
    "",
    "## Boundary",
    "",
    "- Training for this diagnostic uses CHARM/3DG-derived fixed posterior and coordinates from the 016 positive control.",
    "- Additional SNP labels and reference reads inside eval scripts are post-relaxation diagnostics only.",
    "- These fixed-state-edge graphs are not native Hickit raw_expected_soft_all replay; the source_016_before_018_relax row is the native-source baseline.",
    "- Because posterior is fixed, posterior-top1 accuracy is expected to be constant across fixed-state configs; geometry metrics are the main signal.",
    "",
    "## Configs",
    "",
    "| config | graph | native replay | trans dscale | repulsion | min sep | lambda sep |",
    "|---|---:|---:|---:|---:|---:|---:|",
]
for row in rows:
    lines.append(f"| {row['config_name']} | {row['graph_mode']} | {row.get('native_softall_replay','NA')} | {row['trans_dscale_multiplier']} | {row.get('repulsion_multiplier','NA')} | {row['min_sep_unit']} | {row['lambda_sep']} |")
lines += [
    "",
    "## Main Results",
    "",
    "| config | top1 trans | nearest geometry | truth nearest frac | cis Spearman | trans mean r | trans tail k frac | trans force/k |",
    "|---|---:|---:|---:|---:|---:|---:|---:|",
]
for row in rows:
    lines.append(
        f"| {row['config_name']} | {val(row,'model_top1_accuracy_genome_trans')} | "
        f"{val(row,'nearest_geometry_state_accuracy')} | {val(row,'truth_state_is_nearest_fraction')} | "
        f"{val(row,'mean_per_chrom_cis_distance_spearman')} | {val(row,'trans_mean_r')} | "
        f"{val(row,'trans_frac_k_tail_attractive')} | {val(row,'trans_force_per_k')} |"
    )
lines += [
    "",
    "## Headline",
    "",
    f"- best nearest-geometry config: `{best_nearest['config_name'] if best_nearest else 'NA'}` with nearest_geometry_state_accuracy={best_nearest.get('nearest_geometry_state_accuracy','NA') if best_nearest else 'NA'}",
    f"- best cis-Spearman config: `{best_spearman['config_name'] if best_spearman else 'NA'}` with mean_per_chrom_cis_distance_spearman={best_spearman.get('mean_per_chrom_cis_distance_spearman','NA') if best_spearman else 'NA'}",
    "",
    "## Interpretation Rules",
    "",
    "- If fixed-state-edge variants improve nearest-geometry accuracy while posterior-top1 accuracy stays fixed, the bottleneck includes M-step graph/force geometry, not only E-step posterior labels.",
    "- If top1/top2/same-cross pruning improves trans geometry but hurts cis Spearman or separation, it is a diagnostic lead, not a baseline candidate.",
    "- If trans dscale multipliers below 1.0 make normalized trans r larger and increase tail-attractive k fraction, that is not evidence for a useful rescue.",
    "- Do not claim model improvement from this run; this is not blind training.",
]
(root / "README.md").write_text("\n".join(lines) + "\n")
PY

"$SCRIPT_DIR/p9016_publish_light_result.sh" "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$FULL_RUN_ROOT/logs/publish.log"
cat "$FULL_RUN_ROOT/logs/publish.log"
echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
