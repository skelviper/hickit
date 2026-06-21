#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PHASE3_ROOT=$(cd "$REPO_ROOT/.." && pwd)

TEST_RES_ROOT=${HK_BLIND_TEST_RES_ROOT:-"$PHASE3_ROOT/test_res"}
RUN_ID=${HK_BLIND_024_RUN_ID:-"024-$(date +%Y%m%d_%H%M%S)-p9016_pairs_source_signal_audit_1m"}
FULL_RUN_ROOT="$TEST_RES_ROOT/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
DIAG_DIR="$FULL_RUN_ROOT/diagnostics/pairs_source_signal"
LOG_ROOT="$FULL_RUN_ROOT/logs"

PAIRS=${HK_BLIND_P9016_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
SAMPLE=${HK_BLIND_SAMPLE:-P9016}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}

mkdir -p "$DIAG_DIR" "$LOG_ROOT" "$FULL_RUN_ROOT/scripts_snapshot"
cp "$0" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/diagnose_p9016_pairs_source_signal.py" "$FULL_RUN_ROOT/scripts_snapshot/"
cp "$SCRIPT_DIR/p9016_publish_light_result.sh" "$FULL_RUN_ROOT/scripts_snapshot/"

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
	echo "sample	$SAMPLE"
	echo "pairs	$PAIRS"
	echo "bin_size_bp	$BIN_SIZE"
	echo "purpose	pairs_source_signal_audit"
	echo "training_experiment	0"
	echo "eval_only_phase_labels_used_for_truth_stability	1"
	echo "training_uses_phase_labels	0"
	echo "training_uses_charm_or_reference	0"
	echo "git_commit	$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo NA)"
	echo "git_dirty_count	$(git -C "$REPO_ROOT" status --short 2>/dev/null | wc -l | awk '{print $1}')"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

git -C "$REPO_ROOT" status --short > "$LOG_ROOT/git_status.txt" || true
git -C "$REPO_ROOT" diff --stat > "$LOG_ROOT/git_diff_stat.txt" || true
git -C "$REPO_ROOT" diff > "$LOG_ROOT/git_diff.patch" || true
{
	echo "timestamp	$(date -Is)"
	echo "shell	$SHELL"
	echo "python	$(command -v python || true)"
	echo "conda	$(command -v conda || true)"
} > "$LOG_ROOT/build_env.txt"

cd "$REPO_ROOT"

log_cmd python -m py_compile scripts/diagnose_p9016_pairs_source_signal.py
if command -v conda >/dev/null 2>&1; then
	eval "$(conda shell.bash hook)"
	conda activate analysis
fi
python -m py_compile scripts/diagnose_p9016_pairs_source_signal.py

log_cmd python scripts/diagnose_p9016_pairs_source_signal.py --pairs "$PAIRS" --bin-size "$BIN_SIZE" --sample "$SAMPLE" --outdir "$DIAG_DIR"
python scripts/diagnose_p9016_pairs_source_signal.py \
	--pairs "$PAIRS" \
	--bin-size "$BIN_SIZE" \
	--sample "$SAMPLE" \
	--outdir "$DIAG_DIR" > "$LOG_ROOT/pairs_source_signal.log" 2>&1

python - "$FULL_RUN_ROOT" <<'PY'
from __future__ import annotations

import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
diag = root / "diagnostics" / "pairs_source_signal"


def read_kv(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                out[row[0]] = row[1]
    return out


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def val(d: dict[str, str], key: str) -> str:
    return d.get(key, "NA")


summary = read_kv(diag / "pairs_source_signal_summary.tsv")
split = rows(diag / "split_half_truth_reproducibility.tsv")

trans_split = [r for r in split if r.get("scope") == "trans"]
trans_split_focus = [r for r in trans_split if r.get("truth_n_bucket") in {"2", "3", "4", "5", "6-10", "11-20", "21-50"}]

with (root / "summary.tsv").open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "value"])
    for key in [
        "sample",
        "pairs_path",
        "bin_size_bp",
        "n_rows",
        "n_phase_both_rows",
        "n_phase_both_trans_rows",
        "frac_phase_both_rows",
        "readid_signal_status",
        "n_readid_missing_rows",
        "n_readid_nonmissing_rows",
        "repeated_nonmissing_readid_groups",
        "rows_in_repeated_nonmissing_readids",
        "trans_truth_bpairs",
        "trans_truth_contact_weight",
        "trans_truth_singleton_bpairs",
        "trans_truth_singleton_contact_fraction",
    ]:
        writer.writerow([key, val(summary, key)])

readid_status = val(summary, "readid_signal_status")
trans_singleton = float(val(summary, "trans_truth_singleton_contact_fraction") or "nan")

headline = "NO_READID_GROUP_ANCHOR"
if readid_status == "POTENTIALLY_USABLE_REPEATED_READID_GROUPS":
    headline = "POTENTIAL_READID_GROUP_ANCHOR"

with (root / "README.md").open("w") as out:
    out.write("# 024 P9016 Pairs Source Signal Audit\n\n")
    out.write("This is a diagnostic run, not a training experiment. It asks whether the current P9016 `.pairs.gz` file retains blind read/molecule-group information that could anchor trans copy identity beyond single binned pairwise contacts.\n\n")
    out.write("## Paths\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- diagnostics: `{diag}`\n")
    out.write("- light result root: written under `hickit/result` by `scripts/p9016_publish_light_result.sh`\n\n")
    out.write("## Boundary\n\n")
    out.write("- Training uses no SNP, phase, CHARM/3DG, or reference information here because no training is run.\n")
    out.write("- Phase labels are used only for the eval-only split-half truth-stability audit.\n")
    out.write("- CHARM/3DG is not used.\n\n")
    out.write("## Headline\n\n")
    out.write(f"- `{headline}`\n")
    out.write(f"- readID status: `{readid_status}`\n")
    out.write(f"- total rows: `{val(summary, 'n_rows')}`\n")
    out.write(f"- non-missing readID rows: `{val(summary, 'n_readid_nonmissing_rows')}`\n")
    out.write(f"- repeated readID groups: `{val(summary, 'repeated_nonmissing_readid_groups')}`\n")
    out.write(f"- trans singleton contact fraction at 1 Mb: `{val(summary, 'trans_truth_singleton_contact_fraction')}`\n\n")
    out.write("## Split-Half Truth Stability\n\n")
    out.write("| scope | truth_n bucket | contact weight | 4-state agreement | same/cross agreement |\n")
    out.write("|---|---:|---:|---:|---:|\n")
    for r in trans_split_focus:
        out.write(
            f"| {r.get('scope','NA')} | {r.get('truth_n_bucket','NA')} | {r.get('contact_weight','NA')} | "
            f"{r.get('dominant_state_agreement_contact_weighted','NA')} | {r.get('same_cross_agreement_contact_weighted','NA')} |\n"
        )
    out.write("\n## Interpretation\n\n")
    if readid_status.startswith("UNUSABLE"):
        out.write("- The current P9016 pairs file does not preserve usable readID or molecule-group identifiers. A blind multi-contact/read-group trans anchor cannot be built from this file alone.\n")
    else:
        out.write("- Repeated readID groups exist and should be inspected as a possible blind multi-contact anchor before adding another local force knob.\n")
    if trans_singleton == trans_singleton and trans_singleton > 0.5:
        out.write("- Most SNP-labeled trans binned pairs are singleton at 1 Mb, so full-contact four-state trans accuracy is heavily affected by low-count bpair strata.\n")
    out.write("- The next trainable route needs a new blind information source or a changed objective such as high-confidence call/no-call; more sep/dscale tuning alone is unlikely to give a stable +0.1 all-contact trans gain.\n")

with (root / "headline.txt").open("w") as fh:
    fh.write(headline + "\n")
PY

log_cmd scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$LOG_ROOT/publish.log"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "HEADLINE=$(cat "$FULL_RUN_ROOT/headline.txt")"
