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

RUN_ID=${HK_BLIND_038_RUN_ID:-"038-$(date +%Y%m%d_%H%M%S)-p9016_trans_decoder_audit_1m"}
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE_ROOT/test_res}/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
LOG_ROOT="$FULL_RUN_ROOT/logs"
DIAG_ROOT="$FULL_RUN_ROOT/diagnostics"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"

mkdir -p "$LOG_ROOT" "$DIAG_ROOT" "$FULL_RUN_ROOT/scripts_snapshot" "$FULL_RUN_ROOT/outputs" "$FULL_RUN_ROOT/eval"
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)

exec > >(tee -a "$LOG_ROOT/runner.stdout.log") 2> >(tee -a "$LOG_ROOT/runner.stderr.log" >&2)

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

SOURCE_RUN_ROOT=${HK_BLIND_038_SOURCE_RUN_ROOT:-"$PHASE_ROOT/test_res/037-20260619_230802-p9016_trans_chrpair_mstep_1m"}
PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}

CONFIGS=(
	"p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035"
	"p9016_pcgamma1_td0p5_mstep0_sep0"
	"p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0"
	"p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0"
)

log_cmd() {
	local line="+"
	local arg
	for arg in "$@"; do
		line+=" $(printf '%q' "$arg")"
	done
	printf '%s\n' "$line" >> "$COMMANDS_LOG"
	"$@"
}

GIT_COMMIT=$(git rev-parse HEAD 2>/dev/null || echo NA)
GIT_DIRTY_COUNT=$(git status --porcelain 2>/dev/null | wc -l | awk '{print $1}')
git status --porcelain=v1 > "$LOG_ROOT/git_status.txt" || true
git diff --stat > "$LOG_ROOT/git_diff_stat.txt" || true
git diff > "$LOG_ROOT/git_diff.patch" || true
{
	echo "cc	$(${CC:-cc} --version | head -1)"
	echo "make	$(make --version | head -1)"
	echo "conda_env	${CONDA_DEFAULT_ENV:-NA}"
	echo "path	$PATH"
} > "$LOG_ROOT/build_env.txt"

cp scripts/diagnose_p9016_trans_decoder_rules.py \
   scripts/p9016_publish_light_result.sh \
   scripts/run_p9016_038_trans_decoder_audit.sh \
   "$FULL_RUN_ROOT/scripts_snapshot/"

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "purpose	post-training blind trans decoder audit; no retraining"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "source_run_root	$SOURCE_RUN_ROOT"
	echo "pairs_path	$PAIRS"
	echo "reference_3dg_path	$TDG"
	echo "bin_size_bp	$BIN_SIZE"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "training_boundary	no training; SNP labels and CHARM/3DG are eval-only for scoring and denominator matching"
	echo "decoder_boundary	blind_* decoder policies select chromosome-pair swaps from raw posterior bpair weights and reconstruction geometry only; SNP labels and CHARM/3DG enter only after selection for scoring/denominator matching"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

log_cmd python -m py_compile scripts/diagnose_p9016_trans_decoder_rules.py

if [ ! -s "$PAIRS" ]; then
	echo "error: eval pairs missing: $PAIRS" >&2
	exit 1
fi
if [ ! -s "$TDG" ]; then
	echo "error: eval 3DG missing: $TDG" >&2
	exit 1
fi

kv_get() {
	local file=$1
	local key=$2
	awk -F'\t' -v k="$key" 'NR > 1 && $1 == k {print $2; found=1; exit} END {if (!found) exit 1}' "$file"
}

validate_source_manifest() {
	local manifest=$1
	local config=$2
	local source bin uses_phase uses_ref
	if [ ! -s "$manifest" ]; then
		echo "error: missing source manifest for $config: $manifest" >&2
		exit 1
	fi
	source=$(kv_get "$manifest" input_contact_source || echo NA)
	bin=$(kv_get "$manifest" bin_size_bp || kv_get "$manifest" resolution || echo NA)
	uses_phase=$(kv_get "$manifest" uses_phase_labels || echo NA)
	uses_ref=$(kv_get "$manifest" uses_charm_or_reference || echo NA)
	if [ "$source" != "raw_pairs" ]; then
		echo "error: $config input_contact_source=$source, expected raw_pairs" >&2
		exit 1
	fi
	if [ "$bin" != "$BIN_SIZE" ]; then
		echo "error: $config bin_size/resolution=$bin, expected $BIN_SIZE" >&2
		exit 1
	fi
	if [ "$uses_phase" != "0" ]; then
		echo "error: $config uses_phase_labels=$uses_phase, expected 0" >&2
		exit 1
	fi
	if [ "$uses_ref" != "0" ]; then
		echo "error: $config uses_charm_or_reference=$uses_ref, expected 0" >&2
		exit 1
	fi
	printf '%s\t%s\t%s\t%s\t%s\n' "$config" "$source" "$bin" "$uses_phase" "$uses_ref" >> "$FULL_RUN_ROOT/source_manifest_audit.tsv"
}

{
	echo -e "config_name\tinput_contact_source\tbin_size_bp\tuses_phase_labels\tuses_charm_or_reference"
} > "$FULL_RUN_ROOT/source_manifest_audit.tsv"

for config in "${CONFIGS[@]}"; do
	out_dir="$SOURCE_RUN_ROOT/outputs/$config"
	eval_dir="$SOURCE_RUN_ROOT/eval/$config"
	diag_dir="$DIAG_ROOT/$config"
	log_file="$LOG_ROOT/$config.log"
	if [ ! -s "$out_dir/p9016_full.bpair_posterior.tsv" ]; then
		echo "error: missing posterior for $config: $out_dir" >&2
		exit 1
	fi
	if [ ! -s "$out_dir/p9016_full.coords.tsv" ]; then
		echo "error: missing coords for $config: $out_dir" >&2
		exit 1
	fi
	if [ ! -s "$eval_dir/whole_chrom_snp_oracle_swaps.tsv" ]; then
		echo "error: missing eval swaps for $config: $eval_dir" >&2
		exit 1
	fi
	validate_source_manifest "$out_dir/p9016_full.manifest.tsv" "$config"
	mkdir -p "$diag_dir"
	ln -s "$out_dir" "$FULL_RUN_ROOT/outputs/$config"
	ln -s "$eval_dir" "$FULL_RUN_ROOT/eval/$config"
	log_cmd python scripts/diagnose_p9016_trans_decoder_rules.py \
		--config-output-dir "$out_dir" \
		--eval-dir "$eval_dir" \
		--pairs "$PAIRS" \
		--reference-3dg "$TDG" \
		--bin-size "$BIN_SIZE" \
		--config-name "$config" \
		--outdir "$diag_dir" > "$log_file" 2>&1
done

python - "$FULL_RUN_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = []
for path in sorted((root / "diagnostics").glob("*/trans_decoder_rule_summary.tsv")):
    with path.open(newline="") as fh:
        rows.extend(csv.DictReader(fh, delimiter="\t"))

fields = [
    "config_name",
    "policy",
    "flip_source",
    "n_eval_contacts",
    "top1_accuracy",
    "same_cross_accuracy",
    "pmax_threshold_accuracy",
    "pmax_threshold_recall",
    "delta_top1_vs_cis_selected",
    "target_plus_0p1_met",
    "uses_snp_for_base_gauge",
    "uses_snp_for_pair_decoder_selection",
    "uses_snp_labeled_denominator_for_selection_weight",
    "uses_charm_for_denominator_filter",
    "decoder_selection_denominator",
    "decoder_selection_n_records",
    "decoder_selection_total_weight",
]
with (root / "summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: row.get(key, "NA") for key in fields})

headline = "NO_RAW_BLIND_DECODER_PLUS_0P1"
if any(
    row.get("target_plus_0p1_met") == "1"
    and row.get("uses_snp_for_pair_decoder_selection") == "0"
    and row.get("uses_snp_labeled_denominator_for_selection_weight") == "0"
    for row in rows
):
    headline = "RAW_BLIND_DECODER_PLUS_0P1_MET"
(root / "headline.txt").write_text(headline + "\n")
PY

python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
light = Path(sys.argv[2])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))

def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None

by_config = {}
for row in rows:
    by_config.setdefault(row["config_name"], []).append(row)

with (root / "README.md").open("w") as out:
    out.write("# 038 P9016 Trans Decoder Audit\n\n")
    out.write("This is a post-training diagnostic. It does not retrain Hickit. It asks whether current reconstruction/posterior quantities contain a blind-selectable trans copy-state decoder that can close the +0.1 trans top1 gap.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- headline: `{(root / 'headline.txt').read_text().strip()}`\n")
    out.write("- training boundary: no training is run here; source runs are required to have `input_contact_source=raw_pairs`, `uses_phase_labels=0`, and `uses_charm_or_reference=0`.\n")
    out.write("- base gauge: `cis_selected_whole_chrom` uses the standard eval-only whole-chromosome SNP cis oracle gauge, matching the existing evaluator headline.\n")
    out.write("- blind decoder boundary: blind policies select chromosome-pair swaps using raw posterior trans bpair weights (`n_raw`) and reconstruction geometry/posterior fields only.\n")
    out.write("- scoring boundary: SNP labels and CHARM/3DG are used only after decoder selection for scoring and shared denominator matching.\n")
    out.write("- oracle decoder policy uses SNP truth for pair selection and is an eval-only ceiling.\n\n")
    out.write("## Main Table\n\n")
    cols = [
        "config_name",
        "policy",
        "flip_source",
        "top1_accuracy",
        "same_cross_accuracy",
        "delta_top1_vs_cis_selected",
        "target_plus_0p1_met",
        "uses_snp_for_pair_decoder_selection",
        "uses_snp_labeled_denominator_for_selection_weight",
        "decoder_selection_denominator",
    ]
    out.write("| " + " | ".join(cols) + " |\n")
    out.write("|" + "|".join(["---"] * len(cols)) + "|\n")
    for row in sorted(rows, key=lambda r: (r["config_name"], -(fnum(r.get("top1_accuracy")) or -1))):
        out.write("| " + " | ".join(row.get(c, "NA") for c in cols) + " |\n")
    out.write("\n## Interpretation\n\n")
    out.write("If a blind decoder reaches +0.1, the next step is to turn that decoder into a controlled training or post-processing candidate. If only the SNP oracle improves, the current posterior/geometry does not contain enough blind information to select the needed trans copy-state gauge.\n")
PY

log_cmd scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt")"
