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

RUN_ID=${HK_BLIND_051_RUN_ID:-"051-$(date +%Y%m%d_%H%M%S)-p9016_pairs_only_signal_boundary_1m"}
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE_ROOT/test_res}/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
LOG_ROOT="$FULL_RUN_ROOT/logs"
mkdir -p "$LOG_ROOT" "$FULL_RUN_ROOT/scripts_snapshot"
FULL_RUN_ROOT=$(cd "$FULL_RUN_ROOT" && pwd)
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"

exec > >(tee -a "$LOG_ROOT/runner.stdout.log") 2> >(tee -a "$LOG_ROOT/runner.stderr.log" >&2)

log_command_line() {
	local line="+"
	local arg
	for arg in "$@"; do
		line+=" $(printf '%q' "$arg")"
	done
	printf '%s\n' "$line" >> "$COMMANDS_LOG"
}

GIT_COMMIT=$(git rev-parse HEAD 2>/dev/null || echo NA)
GIT_DIRTY_COUNT=$(git status --porcelain | wc -l | awk '{print $1}')
git status --porcelain=v1 > "$LOG_ROOT/git_status.txt" || true
git diff --stat > "$LOG_ROOT/git_diff_stat.txt" || true
git diff > "$LOG_ROOT/git_diff.patch" || true

cp "$0" scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT/scripts_snapshot/"

python - "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
light = Path(sys.argv[2])
result = Path("result")

def read_rows(path: Path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def kv(path: Path):
    rows = read_rows(path)
    return {r["key"]: r["value"] for r in rows}

def get_top1(path: Path, config=None, policy=None):
    rows = read_rows(path)
    best = None
    for r in rows:
        if config is not None and r.get("config_name") != config:
            continue
        if policy is not None and r.get("policy") != policy:
            continue
        v = r.get("model_top1_accuracy_genome_trans") or r.get("top1_accuracy")
        if v in (None, "", "NA", "nan"):
            continue
        try:
            fv = float(v)
        except ValueError:
            continue
        if best is None or fv > best:
            best = fv
    return best

evidence = []
def add(name, category, candidate, trans, target_base, evidence_path, note):
    target = target_base + 0.1 if target_base is not None else None
    evidence.append({
        "evidence_name": name,
        "category": category,
        "strict_pairs_only_candidate": candidate,
        "trans_top1": "NA" if trans is None else f"{trans:.9g}",
        "target_plus_0p1": "NA" if target is None else f"{target:.9g}",
        "meets_plus_0p1": "NA" if trans is None or target is None else str(int(trans >= target)),
        "evidence_path": str(evidence_path),
        "note": note,
    })

current_best = 0.358585
approved_029 = 0.315370

add("current_best_047_init_anchor", "strict_pairs_only_model", "1", current_best, current_best,
    result / "047-20260620_082138-p9016_init_coord_anchor_1m/summary.tsv",
    "Current best strict pairs-only standard trans result used as active target reference.")
add("049_trans_dscale_gamma_best", "strict_pairs_only_model", "1",
    get_top1(result / "049-20260620_103248-p9016_trans_dscale_gamma_1m/summary.tsv"),
    current_best,
    result / "049-20260620_103248-p9016_trans_dscale_gamma_1m/summary.tsv",
    "Trans-specific dscale gamma failed to exceed current best.")
summary050 = kv(result / "050-20260620_113123-p9016_pairs_identifiability_audit_1m/summary.tsv")
add("050_pairs_only_identifiability_baseline", "strict_pairs_only_audit", "1",
    float(summary050["baseline_trans_top1_accuracy"]), current_best,
    result / "050-20260620_113123-p9016_pairs_identifiability_audit_1m/summary.tsv",
    "Approved P9016 pairs have no read linkage or strand variation; high callable accuracy has low recall.")
add("050_best_blind_callable_subset_le50pct", "strict_pairs_only_callable_subset", "0",
    float(summary050["best_blind_top1_accuracy_le_50pct"]), current_best,
    result / "050-20260620_113123-p9016_pairs_identifiability_audit_1m/diagnostics/pairs_identifiability/coverage_accuracy_curve.tsv",
    "High top1 on selected subset only; recall is not full denominator model performance.")
add("038_pair_independent_snp_oracle", "eval_only_oracle", "0",
    get_top1(result / "038-20260620_000849-p9016_trans_decoder_audit_1m/summary.tsv", policy="pair_independent_4state_oracle"),
    current_best,
    result / "038-20260620_000849-p9016_trans_decoder_audit_1m/summary.tsv",
    "SNP oracle decoder improves only about 0.065 over its cis-selected base and still misses +0.1 over current best.")
add("039_global_chrom_snp_oracle", "eval_only_oracle", "0",
    get_top1(result / "039-20260620_002902-p9016_global_gauge_decoder_audit_1m/summary.tsv", policy="global_chrom_snp_oracle"),
    current_best,
    result / "039-20260620_002902-p9016_global_gauge_decoder_audit_1m/summary.tsv",
    "Whole chromosome gauge oracle gives only small improvement.")
add("029_contacts_seg_pairwise", "source_boundary_violation", "0",
    get_top1(result / "029-20260619_061138-p9016_readgroup_joint_marginal_1m/summary.tsv", config="p9016_seg_all_pairwise_pcgamma1_ieps0p5_noise0_seed17"),
    approved_029,
    result / "029-20260619_061138-p9016_readgroup_joint_marginal_1m/summary.tsv",
    "Uses contacts.seg-derived pairs, not approved P9016 pairs; improves but not +0.1.")
add("029_contacts_seg_readgroup_joint", "source_boundary_violation", "0",
    get_top1(result / "029-20260619_061138-p9016_readgroup_joint_marginal_1m/summary.tsv", config="p9016_seg_all_readgroup_joint_pcgamma1_ieps0p5_noise0_seed17"),
    approved_029,
    result / "029-20260619_061138-p9016_readgroup_joint_marginal_1m/summary.tsv",
    "Uses contacts.seg read grouping and is worse than seg pairwise.")

fields = ["evidence_name","category","strict_pairs_only_candidate","trans_top1","target_plus_0p1","meets_plus_0p1","evidence_path","note"]
with (root / "signal_boundary_summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(evidence)

decision = {
    "goal_current_best_trans": f"{current_best:.9g}",
    "goal_target_trans_plus_0p1": f"{current_best + 0.1:.9g}",
    "pairs_only_read_linkage_available": summary050.get("pairs_visible_read_linkage_available", "NA"),
    "pairs_only_strand_information_available": summary050.get("pairs_visible_strand_information_available", "NA"),
    "best_pairs_only_full_denominator_trans_seen": f"{current_best:.9g}",
    "best_pairs_only_model_meets_plus_0p1": "0",
    "best_callable_subset_trans": summary050.get("best_blind_top1_accuracy_le_50pct", "NA"),
    "best_callable_subset_called_fraction": summary050.get("best_blind_called_fraction_le_50pct", "NA"),
    "best_callable_subset_recall": summary050.get("best_blind_recall_le_50pct", "NA"),
    "recommended_next_step": "stop_pairs_only_knob_sweeps; request boundary change to upstream read-level source or report callable/abstention only",
}
with (root / "decision.tsv").open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["key","value"])
    for k,v in decision.items():
        writer.writerow([k,v])

with (root / "summary.tsv").open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["key","value"])
    for k,v in decision.items():
        writer.writerow([k,v])

with (root / "README.md").open("w") as out:
    out.write("# 051 P9016 Pairs-Only Signal Boundary 1Mb\n\n")
    out.write("This report consolidates 024/028/029/038/039/046/049/050 evidence to decide whether another strict P9016.pairs.gz-only knob sweep is likely to deliver +0.1 full-denominator trans top1.\n\n")
    out.write(f"- full result root: `{root}`\n")
    out.write(f"- light result root: `{light}`\n")
    out.write(f"- current best strict pairs-only trans: `{decision['goal_current_best_trans']}`\n")
    out.write(f"- +0.1 target: `{decision['goal_target_trans_plus_0p1']}`\n")
    out.write(f"- pairs-only read linkage available: `{decision['pairs_only_read_linkage_available']}`\n")
    out.write(f"- pairs-only strand information available: `{decision['pairs_only_strand_information_available']}`\n")
    out.write(f"- recommended next step: `{decision['recommended_next_step']}`\n\n")
    out.write("## Evidence Table\n\n")
    out.write("| " + " | ".join(fields) + " |\n")
    out.write("|" + "|".join(["---"] * len(fields)) + "|\n")
    for row in evidence:
        out.write("| " + " | ".join(row[f] for f in fields) + " |\n")
    out.write("\n## Interpretation\n\n")
    out.write("The strict candidate space has not produced a +0.1 full-denominator trans gain. The approved P9016 pairs file has no usable readID/read-group linkage and no strand variation; phase columns are eval-only. High-accuracy trans calls exist only as low-coverage callable subsets. Contacts.seg-derived experiments are useful source-boundary controls but are not valid strict pairs-only candidates. Continuing EM knob sweeps inside the same pairs-only signal surface is therefore unlikely to meet the requested +0.1 goal.\n")
PY

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "experiment_type	signal_boundary_report"
	echo "training_uses_phase_labels	0"
	echo "training_uses_charm_or_reference	0"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

log_command_line scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$FULL_RUN_ROOT/publish.log"
echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "SIGNAL_BOUNDARY_SUMMARY_TSV=$FULL_RUN_ROOT/signal_boundary_summary.tsv"
