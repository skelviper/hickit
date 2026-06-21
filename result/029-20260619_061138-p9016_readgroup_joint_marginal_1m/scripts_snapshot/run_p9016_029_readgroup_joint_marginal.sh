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

TEST_RES_ROOT=${HK_BLIND_TEST_RES_ROOT:-"$PHASE3_ROOT/test_res"}
RUN_ID=${HK_BLIND_029_RUN_ID:-"029-$(date +%Y%m%d_%H%M%S)-p9016_readgroup_joint_marginal_1m"}
FULL_RUN_ROOT="$TEST_RES_ROOT/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
SEG=${HK_BLIND_P9016_CONTACTS_SEG:-/sharec/zliu/CHARM/mESC/processed/P9016/2d_info/contacts.seg.gz}
APPROVED_PAIRS=${HK_BLIND_P9016_APPROVED_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
EVAL_PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
EVAL_TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
N_ITER=${HK_BLIND_P9016_MINIMAL_N_ITER:-100}
RELAX_STEPS=${HK_BLIND_P9016_MINIMAL_RELAX_STEPS:-100}
RELAX_STEP=${HK_BLIND_P9016_RELAX_STEP:-${HK_BLIND_P9016_MINIMAL_RELAX_STEP:-0.012}}
BACKEND=${HK_BLIND_P9016_RELAX_BACKEND:-auto}
MAX_PARALLEL=${HK_BLIND_029_MAX_PARALLEL:-3}
MAX_PAIRS=${HK_BLIND_029_MAX_PAIRS:-1000000}

PAIRS_ROOT="$FULL_RUN_ROOT/generated_pairs"
OUTPUT_ROOT="$FULL_RUN_ROOT/outputs"
WORK_OUTPUT_ROOT="$FULL_RUN_ROOT/work_outputs"
EVAL_ROOT="$FULL_RUN_ROOT/eval"
LOG_ROOT="$FULL_RUN_ROOT/logs"

mkdir -p "$PAIRS_ROOT" "$OUTPUT_ROOT" "$WORK_OUTPUT_ROOT" "$EVAL_ROOT" "$LOG_ROOT" "$FULL_RUN_ROOT/scripts_snapshot"
cp "$0" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/export_p9016_contacts_seg_pairs.py" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/p9016_common_eval.sh" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/p9016_publish_light_result.sh" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/summarize_p9016_model_sweep.py" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/compute_p9016_copytrack_diag.py" "$FULL_RUN_ROOT/scripts_snapshot/"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

log_cmd() {
	printf '%q ' "$@" >> "$FULL_RUN_ROOT/commands.log"
	printf '\n' >> "$FULL_RUN_ROOT/commands.log"
}

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "contacts_seg	$SEG"
	echo "approved_pairs	$APPROVED_PAIRS"
	echo "eval_pairs	$EVAL_PAIRS"
	echo "eval_tdg	$EVAL_TDG"
	echo "bin_size_bp	$BIN_SIZE"
	echo "n_iter	$N_ITER"
	echo "relax_steps	$RELAX_STEPS"
	echo "relax_step	$RELAX_STEP"
	echo "backend	$BACKEND"
	echo "max_parallel	$MAX_PARALLEL"
	echo "max_pairs_per_export	$MAX_PAIRS"
	echo "training_boundary	readgroup joint marginal uses contacts.seg read grouping only; phase labels and CHARM/3DG are eval-only"
	echo "git_commit	$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo NA)"
	echo "git_dirty_count	$(git -C "$REPO_ROOT" status --short 2>/dev/null | wc -l | awk '{print $1}')"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

git -C "$REPO_ROOT" status --short > "$LOG_ROOT/git_status.txt" || true
git -C "$REPO_ROOT" diff --stat > "$LOG_ROOT/git_diff_stat.txt" || true
git -C "$REPO_ROOT" diff > "$LOG_ROOT/git_diff.patch" || true

cd "$REPO_ROOT"
log_cmd make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
log_cmd python -m py_compile scripts/export_p9016_contacts_seg_pairs.py
python -m py_compile scripts/export_p9016_contacts_seg_pairs.py

export_pairs() {
	local mode=$1
	local out="$PAIRS_ROOT/P9016.seg_${mode}.readgroup.pairs.gz"
	local meta="$PAIRS_ROOT/P9016.seg_${mode}.readgroup.metadata.tsv"
	if [ -s "$out" ]; then
		return
	fi
	log_cmd python scripts/export_p9016_contacts_seg_pairs.py --seg "$SEG" --out-pairs "$out" --metadata "$meta" --mode "$mode" --bin-size "$BIN_SIZE" --drop-same-bin --max-pairs "$MAX_PAIRS" --seed 17 --write-readgroup-cols
	python scripts/export_p9016_contacts_seg_pairs.py \
		--seg "$SEG" \
		--out-pairs "$out" \
		--metadata "$meta" \
		--mode "$mode" \
		--bin-size "$BIN_SIZE" \
		--drop-same-bin \
		--max-pairs "$MAX_PAIRS" \
		--seed 17 \
		--write-readgroup-cols > "$LOG_ROOT/export_${mode}.log" 2>&1
}

export_pairs all
export_pairs multiseg_multichrom

CONFIGS=(
	"approved_baseline|$APPROVED_PAIRS|0|raw_pairs|off"
	"seg_all_pairwise|$PAIRS_ROOT/P9016.seg_all.readgroup.pairs.gz|1|contacts_seg_derived_pairs|off"
	"seg_all_readgroup_joint|$PAIRS_ROOT/P9016.seg_all.readgroup.pairs.gz|1|contacts_seg_derived_pairs|joint_marginal"
	"seg_multiseg_multichrom_pairwise|$PAIRS_ROOT/P9016.seg_multiseg_multichrom.readgroup.pairs.gz|1|contacts_seg_derived_pairs|off"
	"seg_multiseg_multichrom_readgroup_joint|$PAIRS_ROOT/P9016.seg_multiseg_multichrom.readgroup.pairs.gz|1|contacts_seg_derived_pairs|joint_marginal"
)

run_one() {
	local spec=$1
	IFS='|' read -r label pairs allow_custom source readgroup_mode <<< "$spec"
	local config="p9016_${label}_pcgamma1_ieps0p5_noise0_seed17"
	local run_root="$WORK_OUTPUT_ROOT/$config"
	local out_dir="$run_root/$config"
	local link_dir="$OUTPUT_ROOT/$config"
	local eval_dir="$EVAL_ROOT/$config"
	local log_file="$LOG_ROOT/$config.log"
	rm -rf "$run_root" "$link_dir" "$eval_dir"
	mkdir -p "$run_root" "$eval_dir"
	{
		echo "CONFIG=$config"
		echo "pairs=$pairs"
		echo "allow_custom=$allow_custom"
		echo "source=$source"
		echo "readgroup_mode=$readgroup_mode"
	} > "$log_file"
	(
		export HK_BLIND_SAMPLE=P9016
		export HK_BLIND_P9016_SAMPLE=P9016
		export HK_BLIND_P9016_PAIRS="$pairs"
		export HK_BLIND_P9016_EVAL_PAIRS="$EVAL_PAIRS"
		export HK_BLIND_P9016_EVAL_TDG="$EVAL_TDG"
		export HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS="$allow_custom"
		export HK_BLIND_P9016_INPUT_CONTACT_SOURCE="$source"
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
		export HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA=1
		export HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6
		export HK_BLIND_P9016_READGROUP_MODE="$readgroup_mode"
		export HK_BLIND_P9016_READGROUP_MAX_SEGMENTS=8
		export HK_BLIND_P9016_READGROUP_EPS=1e-6
		export HK_BLIND_WRITE_RAW=0
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
		echo -e "status\tFAILED" > "$eval_dir/eval_status.tsv"
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
from __future__ import annotations
import csv, math, sys
from pathlib import Path

root = Path(sys.argv[1])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))
base = next((r for r in rows if "approved_baseline" in r.get("config_name", "")), None)
base_trans = float(base["model_top1_accuracy_genome_trans"]) if base and base.get("model_top1_accuracy_genome_trans") not in {"", "NA"} else float("nan")
for r in rows:
    try:
        t = float(r.get("model_top1_accuracy_genome_trans", "nan"))
        r["delta_trans_vs_approved_baseline"] = f"{t - base_trans:.9g}"
        r["target_plus0p1_full_trans_met"] = "1" if t >= base_trans + 0.1 else "0"
    except ValueError:
        r["delta_trans_vs_approved_baseline"] = "NA"
        r["target_plus0p1_full_trans_met"] = "0"

def numeric(value: str | None, default: float = -999.0) -> float:
    try:
        x = float(value if value not in (None, "") else default)
    except ValueError:
        return default
    return x if math.isfinite(x) else default

fields = [
    "config_name", "input_contact_source", "readgroup_mode", "readgroup_entries",
    "readgroup_groups_used", "readgroup_raw_decoded",
    "model_top1_accuracy_genome_all", "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans", "delta_trans_vs_approved_baseline",
    "target_plus0p1_full_trans_met", "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman", "final_mean_entropy", "final_mean_pU",
    "final_min_sep", "final_mean_sep", "copytrack_frac_cos_lt_0",
]
with (root / "readgroup_delta_summary.tsv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    for r in sorted(rows, key=lambda x: numeric(x.get("delta_trans_vs_approved_baseline")), reverse=True):
        w.writerow(r)
headline = "FULL_TRANS_PLUS_0P1" if any(r.get("target_plus0p1_full_trans_met") == "1" for r in rows) else "NO_FULL_TRANS_PLUS_0P1"
(root / "headline.txt").write_text(headline + "\n")
with (root / "README.md").open("w") as out:
    out.write("# 029 P9016 Readgroup Joint-Marginal Training\n\n")
    out.write("This experiment tests a blind readgroup-aware E-step prototype. The runner uses `contacts.seg.gz` read grouping only; phase labels and CHARM/3DG are used only by eval after training.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write("- light result root: under `hickit/result` after publish\n")
    out.write(f"- headline: `{headline}`\n\n")
    out.write("| config | readgroup | top1 all | top1 cis | top1 trans | delta trans | same/cross trans | cis Spearman | entropy | pU |\n")
    out.write("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    for r in sorted(rows, key=lambda x: numeric(x.get("delta_trans_vs_approved_baseline")), reverse=True):
        out.write(f"| `{r.get('config_name','NA')}` | {r.get('readgroup_mode','NA')} | {r.get('model_top1_accuracy_genome_all','NA')} | {r.get('model_top1_accuracy_genome_cis','NA')} | {r.get('model_top1_accuracy_genome_trans','NA')} | {r.get('delta_trans_vs_approved_baseline','NA')} | {r.get('model_same_cross_accuracy_genome_trans','NA')} | {r.get('mean_per_chrom_cis_distance_spearman','NA')} | {r.get('final_mean_entropy','NA')} | {r.get('final_mean_pU','NA')} |\n")
    out.write("\n## Interpretation\n\n")
    if headline == "FULL_TRANS_PLUS_0P1":
        out.write("- The readgroup-aware prototype reached full-denominator trans +0.1 over approved-pairs baseline.\n")
        out.write("- Next step is independent code review plus seed/source replication before promoting it.\n")
    else:
        out.write("- The readgroup-aware prototype did not reach full-denominator trans +0.1 in this run.\n")
        out.write("- If it improves only slightly, read grouping is not sufficient by itself and the remaining problem is likely global trans copy identity/gauge propagation.\n")
PY

scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$LOG_ROOT/publish.log"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "READGROUP_DELTA_TSV=$FULL_RUN_ROOT/readgroup_delta_summary.tsv"
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt")"
