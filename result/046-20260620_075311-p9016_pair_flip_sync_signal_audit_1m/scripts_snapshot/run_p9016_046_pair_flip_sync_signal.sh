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

RUN_ID=${HK_BLIND_046_RUN_ID:-"046-$(date +%Y%m%d_%H%M%S)-p9016_pair_flip_sync_signal_audit_1m"}
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE_ROOT/test_res}/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
LOG_ROOT="$FULL_RUN_ROOT/logs"
DIAG_ROOT="$FULL_RUN_ROOT/diagnostics"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"

mkdir -p "$LOG_ROOT" "$DIAG_ROOT" "$FULL_RUN_ROOT/scripts_snapshot" "$FULL_RUN_ROOT/outputs" "$FULL_RUN_ROOT/eval" "$FULL_RUN_ROOT/plots"
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)

exec > >(tee -a "$LOG_ROOT/runner.stdout.log") 2> >(tee -a "$LOG_ROOT/runner.stderr.log" >&2)

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
TDG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
SPLIT_SEED=${HK_BLIND_046_SPLIT_SEED:-17}
MAX_PARALLEL=${HK_BLIND_046_MAX_PARALLEL:-4}

if [ "$BIN_SIZE" != "1000000" ]; then
	echo "error: experiment 046 is fixed at 1Mb; HK_BLIND_P9016_BIN_SIZE_BP=$BIN_SIZE" >&2
	exit 2
fi
case "$SPLIT_SEED" in
	''|*[!0-9]*)
		echo "error: HK_BLIND_046_SPLIT_SEED must be an unsigned integer; got '$SPLIT_SEED'" >&2
		exit 2
		;;
esac
case "$MAX_PARALLEL" in
	''|*[!0-9]*)
		echo "error: HK_BLIND_046_MAX_PARALLEL must be a positive integer; got '$MAX_PARALLEL'" >&2
		exit 2
		;;
esac
if [ "$MAX_PARALLEL" -lt 1 ]; then
	echo "error: HK_BLIND_046_MAX_PARALLEL must be >= 1" >&2
	exit 2
fi

if [ -n "${HK_BLIND_046_SOURCE_CONFIGS:-}" ]; then
	mapfile -t SOURCE_CONFIGS <<< "$HK_BLIND_046_SOURCE_CONFIGS"
else
	SOURCE_CONFIGS=(
		"$PHASE_ROOT/test_res/037-20260619_230802-p9016_trans_chrpair_mstep_1m::p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035"
		"$PHASE_ROOT/test_res/037-20260619_230802-p9016_trans_chrpair_mstep_1m::p9016_pcgamma1_td0p5_mstep0_sep0"
		"$PHASE_ROOT/test_res/037-20260619_230802-p9016_trans_chrpair_mstep_1m::p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0"
		"$PHASE_ROOT/test_res/037-20260619_230802-p9016_trans_chrpair_mstep_1m::p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0"
		"$PHASE_ROOT/test_res/044-20260620_040939-p9016_multiseed_basin_audit_1m::p9016_multiseed_best035_seed71_noise0p05"
		"$PHASE_ROOT/test_res/044-20260620_040939-p9016_multiseed_basin_audit_1m::p9016_multiseed_best035_seed17_noise0"
		"$PHASE_ROOT/test_res/044-20260620_040939-p9016_multiseed_basin_audit_1m::p9016_multiseed_chrpair_m0p5_p4_w10_seed71_noise0p05"
		"$PHASE_ROOT/test_res/045-20260620_055343-p9016_raw_heldout_selection_1m::p9016_heldout_chrpair_m1_p1_w0_seed23_noise0p05"
		"$PHASE_ROOT/test_res/045-20260620_055343-p9016_raw_heldout_selection_1m::p9016_heldout_best035_seed31_noise0p05"
		"$PHASE_ROOT/test_res/045-20260620_055343-p9016_raw_heldout_selection_1m::p9016_heldout_td1_seed17_noise0"
	)
fi

log_cmd() {
	local line="+"
	local arg
	for arg in "$@"; do
		line+=" $(printf '%q' "$arg")"
	done
	printf '%s\n' "$line" >> "$COMMANDS_LOG"
	"$@"
}

kv_get() {
	local file=$1
	local key=$2
	awk -F'\t' -v k="$key" 'NR > 1 && $1 == k {print $2; found=1; exit} END {if (!found) exit 1}' "$file"
}

resolve_config_output_dir() {
	local source_root=$1
	local config=$2
	local base="$source_root/outputs/$config"
	if [ -s "$base/p9016_full.bpair_posterior.tsv" ] && [ -s "$base/p9016_full.coords.tsv" ]; then
		printf '%s\n' "$base"
		return 0
	fi
	if [ -s "$base/$config/p9016_full.bpair_posterior.tsv" ] && [ -s "$base/$config/p9016_full.coords.tsv" ]; then
		printf '%s\n' "$base/$config"
		return 0
	fi
	return 1
}

validate_source_manifest() {
	local manifest=$1
	local config=$2
	local source_root=$3
	local status source bin res_label uses_phase uses_ref uses_charm_train gauge_only baseline graph_mode
	if [ ! -s "$manifest" ]; then
		echo "error: missing source manifest for $config: $manifest" >&2
		exit 1
	fi
	status=$(kv_get "$manifest" status || echo NA)
	source=$(kv_get "$manifest" input_contact_source || echo NA)
	bin=$(kv_get "$manifest" bin_size_bp || kv_get "$manifest" resolution || echo NA)
	res_label=$(kv_get "$manifest" resolution_label || echo NA)
	uses_phase=$(kv_get "$manifest" uses_phase_labels || echo NA)
	uses_ref=$(kv_get "$manifest" uses_charm_or_reference || echo NA)
	uses_charm_train=$(kv_get "$manifest" uses_charm_for_training || echo NA)
	gauge_only=$(kv_get "$manifest" copy_labels_are_gauge_only || echo NA)
	baseline=$(kv_get "$manifest" baseline || echo NA)
	graph_mode=$(kv_get "$manifest" mstep_graph_mode || echo NA)
	if [ "$status" != "OK" ]; then
		echo "error: $config status=$status, expected OK" >&2
		exit 1
	fi
	if [ "$source" != "raw_pairs" ]; then
		echo "error: $config input_contact_source=$source, expected raw_pairs" >&2
		exit 1
	fi
	if [ "$bin" != "$BIN_SIZE" ]; then
		echo "error: $config bin_size/resolution=$bin, expected $BIN_SIZE" >&2
		exit 1
	fi
	if [ "$res_label" != "1Mb" ]; then
		echo "error: $config resolution_label=$res_label, expected 1Mb" >&2
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
	if [ "$uses_charm_train" != "0" ]; then
		echo "error: $config uses_charm_for_training=$uses_charm_train, expected 0" >&2
		exit 1
	fi
	if [ "$gauge_only" != "1" ]; then
		echo "error: $config copy_labels_are_gauge_only=$gauge_only, expected 1" >&2
		exit 1
	fi
	if [ "$baseline" != "softall" ]; then
		echo "error: $config baseline=$baseline, expected softall" >&2
		exit 1
	fi
	if [ "$graph_mode" != "raw_expected_soft_all" ]; then
		echo "error: $config mstep_graph_mode=$graph_mode, expected raw_expected_soft_all" >&2
		exit 1
	fi
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$config" "$source_root" "$status" "$source" "$bin" "$res_label" "$uses_phase" "$uses_ref" "$uses_charm_train" "$gauge_only" "$graph_mode" >> "$FULL_RUN_ROOT/source_manifest_audit.tsv"
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

cp scripts/diagnose_p9016_pair_flip_sync_signal.py \
   scripts/p9016_publish_light_result.sh \
   scripts/run_p9016_046_pair_flip_sync_signal.sh \
   "$FULL_RUN_ROOT/scripts_snapshot/"

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "purpose	post-training raw-only chromosome-pair copy-flip sync signal audit; no retraining"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "pairs_path	$PAIRS"
	echo "reference_3dg_path	$TDG"
	echo "bin_size_bp	$BIN_SIZE"
	echo "split_seed	$SPLIT_SEED"
	echo "max_parallel	$MAX_PARALLEL"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "training_boundary	no training is run here; source runs must be raw_pairs blind training"
	echo "selector_boundary	raw_* policies select chromosome-pair swaps from posterior bpair weights and reconstruction geometry only; SNP labels and CHARM/3DG enter only after selection for scoring and denominator matching"
	echo "copy_swap_policy	copy0/copy1 are gauge labels; cis SNP gauge is eval-only and pair flips are diagnostic"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

{
	echo -e "config_name\tsource_run_root\tstatus\tinput_contact_source\tbin_size_bp\tresolution_label\tuses_phase_labels\tuses_charm_or_reference\tuses_charm_for_training\tcopy_labels_are_gauge_only\tmstep_graph_mode"
} > "$FULL_RUN_ROOT/source_manifest_audit.tsv"

log_cmd python -m py_compile scripts/diagnose_p9016_pair_flip_sync_signal.py

if [ ! -s "$PAIRS" ]; then
	echo "error: eval pairs missing: $PAIRS" >&2
	exit 1
fi
if [ ! -s "$TDG" ]; then
	echo "error: eval 3DG missing: $TDG" >&2
	exit 1
fi

RESOLVED_ENTRIES=()
for entry in "${SOURCE_CONFIGS[@]}"; do
	[ -n "$entry" ] || continue
	source_root=${entry%%::*}
	config=${entry##*::}
	if [ "$source_root" = "$entry" ] || [ -z "$source_root" ] || [ -z "$config" ]; then
		echo "error: malformed source config entry '$entry'; expected <source_run_root>::<config_name>" >&2
		exit 2
	fi
	out_dir=$(resolve_config_output_dir "$source_root" "$config") || {
		echo "error: missing posterior/coords for $config under $source_root" >&2
		exit 1
	}
	eval_dir="$source_root/eval/$config"
	if [ ! -s "$eval_dir/summary.tsv" ]; then
		echo "error: missing eval summary for $config: $eval_dir" >&2
		exit 1
	fi
	if [ ! -s "$eval_dir/whole_chrom_snp_oracle_swaps.tsv" ]; then
		echo "error: missing eval whole-chrom swaps for $config: $eval_dir" >&2
		exit 1
	fi
	validate_source_manifest "$out_dir/p9016_full.manifest.tsv" "$config" "$source_root"
	mkdir -p "$DIAG_ROOT/$config"
	ln -sfn "$out_dir" "$FULL_RUN_ROOT/outputs/$config"
	ln -sfn "$eval_dir" "$FULL_RUN_ROOT/eval/$config"
	RESOLVED_ENTRIES+=("$source_root::$config::$out_dir::$eval_dir::$DIAG_ROOT/$config::$LOG_ROOT/$config.log")
done

run_one_config() {
	local entry=$1
	local rest source_root config out_dir eval_dir diag_dir log_file
	source_root=${entry%%::*}
	rest=${entry#*::}
	config=${rest%%::*}
	rest=${rest#*::}
	out_dir=${rest%%::*}
	rest=${rest#*::}
	eval_dir=${rest%%::*}
	rest=${rest#*::}
	diag_dir=${rest%%::*}
	log_file=${rest#*::}
	{
		printf '+ python scripts/diagnose_p9016_pair_flip_sync_signal.py --config-output-dir %q --eval-dir %q --pairs %q --reference-3dg %q --bin-size %q --config-name %q --split-seed %q --outdir %q\n' \
			"$out_dir" "$eval_dir" "$PAIRS" "$TDG" "$BIN_SIZE" "$config" "$SPLIT_SEED" "$diag_dir" >> "$COMMANDS_LOG"
		python scripts/diagnose_p9016_pair_flip_sync_signal.py \
			--config-output-dir "$out_dir" \
			--eval-dir "$eval_dir" \
			--pairs "$PAIRS" \
			--reference-3dg "$TDG" \
			--bin-size "$BIN_SIZE" \
			--config-name "$config" \
			--split-seed "$SPLIT_SEED" \
			--outdir "$diag_dir"
	} > "$log_file" 2>&1
}

pids=()
failed=0
for entry in "${RESOLVED_ENTRIES[@]}"; do
	run_one_config "$entry" &
	pids+=("$!")
	if [ "${#pids[@]}" -ge "$MAX_PARALLEL" ]; then
		pid=${pids[0]}
		if ! wait "$pid"; then
			failed=1
		fi
		pids=("${pids[@]:1}")
	fi
done
for pid in "${pids[@]}"; do
	if ! wait "$pid"; then
		failed=1
	fi
done
if [ "$failed" -ne 0 ]; then
	echo "error: at least one pair-flip sync diagnostic failed; see $LOG_ROOT" >&2
	exit 1
fi

python - "$FULL_RUN_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])

def read_rows(name: str) -> list[dict[str, str]]:
    rows = []
    for path in sorted((root / "diagnostics").glob(f"*/{name}")):
        with path.open(newline="") as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    return rows

summary_rows = read_rows("pair_flip_sync_summary.tsv")
audit_rows = read_rows("pair_flip_sync_audit.tsv")
pair_rows = read_rows("pair_flip_sync_chr_pair.tsv")
score_rows = read_rows("pair_flip_signal_scores.tsv")

def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fields})

write_table(root / "summary.tsv", summary_rows)
write_table(root / "pair_flip_sync_signal_audit.tsv", audit_rows)
write_table(root / "pair_flip_sync_chr_pair.tsv", pair_rows)
write_table(root / "pair_flip_signal_scores.tsv", score_rows)

headline = "INSUFFICIENT_DATA"
blind_rows = [
    row for row in summary_rows
    if row.get("selection_source") == "blind_raw_posterior_reconstruction_geometry"
    and row.get("uses_snp_for_pair_selection") == "0"
]
def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")

if blind_rows:
    best = max(blind_rows, key=lambda r: fnum(r.get("delta_top1_vs_cis_selected")))
    delta = fnum(best.get("delta_top1_vs_cis_selected"))
    if any(row.get("target_plus_0p1_met") == "1" for row in blind_rows):
        headline = "RAW_PAIR_FLIP_SYNC_PLUS_0P1_MET"
    elif delta >= 0.05:
        headline = "RAW_PAIR_FLIP_SYNC_PLUS_0P05_ONLY"
    else:
        headline = "NO_RAW_PAIR_FLIP_SYNC_PLUS_0P05"
    (root / "best_blind_pair_flip_policy.tsv").write_text(
        "\t".join(best.keys()) + "\n" + "\t".join(best.get(k, "NA") for k in best.keys()) + "\n"
    )
(root / "headline.txt").write_text(headline + "\n")
PY

python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
light = Path(sys.argv[2])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t")) if (root / "summary.tsv").exists() else []

def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")

headline = (root / "headline.txt").read_text().strip() if (root / "headline.txt").exists() else "NA"
blind = [r for r in rows if r.get("selection_source") == "blind_raw_posterior_reconstruction_geometry"]
top = sorted(blind, key=lambda r: fnum(r.get("delta_top1_vs_cis_selected")), reverse=True)[:20]

with (root / "plots" / "README.md").open("w") as out:
    out.write("No plot is generated for 046. This audit is a tabular post-training pair-flip signal diagnostic.\n")

with (root / "README.md").open("w") as out:
    out.write("# 046 P9016 Pair-Flip Sync Signal Audit\n\n")
    out.write("This is a post-training diagnostic. It does not retrain Hickit. It asks whether raw posterior/reconstruction geometry contains a stable blind signal for chromosome-pair copy flips that could explain the remaining trans deficit.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- headline: `{headline}`\n")
    out.write("- resolution: 1Mb\n")
    out.write("- training inputs: inherited source runs, all audited as `input_contact_source=raw_pairs`, `uses_phase_labels=0`, and `uses_charm_or_reference=0`.\n")
    out.write("- eval inputs: P9016 SNP labels from pairs and P9016 1Mb CHARM/3DG are used only after the blind pair-flip choices are fixed.\n")
    out.write("- copy policy: copy0/copy1 are gauge labels. The cis whole-chromosome SNP gauge is eval-only alignment; raw pair-flip policies must not use SNP for selection.\n\n")
    out.write("Important wording: `blind raw policy` in this README means the extra chromosome-pair flip is selected without SNP or CHARM/3DG. The reporting frame still uses the standard eval-only whole-chromosome SNP cis gauge and eval-only SNP/CHARM denominator matching.\n\n")
    out.write("## Main Blind Rows\n\n")
    cols = [
        "config_name",
        "policy",
        "base_policy",
        "policy_family",
        "top1_accuracy",
        "delta_top1_vs_cis_selected",
        "same_cross_accuracy",
        "split_agree_pair_fraction",
        "oracle_match_pair_fraction",
        "oracle_match_eval_contact_fraction",
        "target_plus_0p1_met",
    ]
    out.write("| " + " | ".join(cols) + " |\n")
    out.write("|" + "|".join(["---"] * len(cols)) + "|\n")
    for row in top:
        out.write("| " + " | ".join(row.get(col, "NA") for col in cols) + " |\n")
    out.write("\n## Interpretation Rule\n\n")
    out.write("If a blind raw policy reaches +0.1 full-denominator trans top1, it is a candidate to turn into a controlled training or post-processing experiment. A +0.05 headline is only an exploratory warning sign, not the predeclared success target. If only the SNP oracle improves, the data still contain pair-level eval-only space but the current raw posterior/geometry does not expose a blind selector. Split-agreement rows test whether the raw pair-flip signal is stable enough to trust without SNP labels.\n")
PY

log_cmd scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt")"
