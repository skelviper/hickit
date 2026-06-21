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

RUN_ID=${HK_BLIND_052_RUN_ID:-"052-$(date +%Y%m%d_%H%M%S)-p9016_manifest_strict_boundary_audit_1m"}
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE_ROOT/test_res}/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
LOG_ROOT="$FULL_RUN_ROOT/logs"
mkdir -p "$LOG_ROOT" "$FULL_RUN_ROOT/scripts_snapshot" "$FULL_RUN_ROOT/plots"
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

cp "$0" "$SCRIPT_DIR/audit_p9016_manifest_strict_boundary.py" "$SCRIPT_DIR/p9016_publish_light_result.sh" "$FULL_RUN_ROOT/scripts_snapshot/"

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "test_res_root	$PHASE_ROOT/test_res"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "experiment_type	manifest_strict_boundary_audit"
	echo "resolution	1000000"
	echo "training_run	0"
	echo "training_uses_phase_labels	0"
	echo "training_uses_charm_or_reference	0"
	echo "baseline_trans	${HK_BLIND_052_BASELINE_TRANS:-0.358585}"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

log_command_line python scripts/audit_p9016_manifest_strict_boundary.py "$PHASE_ROOT/test_res" --outdir "$FULL_RUN_ROOT" --baseline-trans "${HK_BLIND_052_BASELINE_TRANS:-0.358585}"
python scripts/audit_p9016_manifest_strict_boundary.py \
	"$PHASE_ROOT/test_res" \
	--outdir "$FULL_RUN_ROOT" \
	--baseline-trans "${HK_BLIND_052_BASELINE_TRANS:-0.358585}" \
	> "$LOG_ROOT/audit.log" 2>&1

log_command_line scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$FULL_RUN_ROOT/publish.log"

cat "$FULL_RUN_ROOT/summary.tsv"
echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "STRICT_LEADERBOARD_TSV=$FULL_RUN_ROOT/strict_pairs_only_leaderboard.tsv"
echo "MANIFEST_BOUNDARY_AUDIT_TSV=$FULL_RUN_ROOT/manifest_boundary_audit.tsv"
