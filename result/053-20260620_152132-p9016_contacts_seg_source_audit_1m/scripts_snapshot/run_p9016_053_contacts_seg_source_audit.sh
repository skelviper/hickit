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

RUN_ID=${HK_BLIND_053_RUN_ID:-"053-$(date +%Y%m%d_%H%M%S)-p9016_contacts_seg_source_audit_1m"}
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

SEG_SHAREC=${HK_BLIND_P9016_CONTACTS_SEG_SHAREC:-/sharec/zliu/CHARM/mESC/processed/P9016/2d_info/contacts.seg.gz}
SEG_ARCHIVE=${HK_BLIND_P9016_CONTACTS_SEG_ARCHIVE:-/shared/zliu/CHARM/CHARM_mesc/3_cellcycle/archieve/coverage/segs/P9016.sort.contacts_seg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
MIN_MAPQ=${HK_BLIND_P9016_SEG_MIN_MAPQ:-30}

GIT_COMMIT=$(git rev-parse HEAD 2>/dev/null || echo NA)
GIT_DIRTY_COUNT=$(git status --porcelain | wc -l | awk '{print $1}')
git status --porcelain=v1 > "$LOG_ROOT/git_status.txt" || true
git diff --stat > "$LOG_ROOT/git_diff_stat.txt" || true
git diff > "$LOG_ROOT/git_diff.patch" || true

cp "$0" "$SCRIPT_DIR/audit_p9016_contacts_seg_sources.py" "$SCRIPT_DIR/p9016_publish_light_result.sh" "$FULL_RUN_ROOT/scripts_snapshot/"

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "experiment_type	contacts_seg_source_audit"
	echo "training_run	0"
	echo "resolution	$BIN_SIZE"
	echo "min_mapq	$MIN_MAPQ"
	echo "sharec_seg	$SEG_SHAREC"
	echo "archive_seg	$SEG_ARCHIVE"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
echo "SEG_SHAREC=$SEG_SHAREC"
echo "SEG_ARCHIVE=$SEG_ARCHIVE"

if [ ! -s "$SEG_SHAREC" ]; then
	echo "missing SEG_SHAREC: $SEG_SHAREC" >&2
	exit 2
fi
if [ ! -s "$SEG_ARCHIVE" ]; then
	echo "missing SEG_ARCHIVE: $SEG_ARCHIVE" >&2
	exit 2
fi

log_command_line python scripts/audit_p9016_contacts_seg_sources.py --source "sharec=$SEG_SHAREC" --source "archive=$SEG_ARCHIVE" --outdir "$FULL_RUN_ROOT" --bin-size "$BIN_SIZE" --min-mapq "$MIN_MAPQ"
python scripts/audit_p9016_contacts_seg_sources.py \
	--source "sharec=$SEG_SHAREC" \
	--source "archive=$SEG_ARCHIVE" \
	--outdir "$FULL_RUN_ROOT" \
	--bin-size "$BIN_SIZE" \
	--min-mapq "$MIN_MAPQ" \
	> "$LOG_ROOT/audit.log" 2>&1

log_command_line scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$FULL_RUN_ROOT/publish.log"

echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "SEG_SOURCE_COMPARISON_TSV=$FULL_RUN_ROOT/seg_source_comparison.tsv"
echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
