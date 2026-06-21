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

RUN_ID=${HK_BLIND_054_RUN_ID:-"054-$(date +%Y%m%d_%H%M%S)-p9016_contacts_seg_readgroup_dedup_audit_1m"}
FULL_RUN_ROOT="${HK_BLIND_TEST_RES_ROOT:-$PHASE_ROOT/test_res}/$RUN_ID"
LIGHT_RESULT_ROOT="$REPO_ROOT/result/$RUN_ID"
LOG_ROOT="$FULL_RUN_ROOT/logs"
EXPORT_ROOT="$FULL_RUN_ROOT/exports"
COMMANDS_LOG="$FULL_RUN_ROOT/commands.log"

if [ -e "$FULL_RUN_ROOT" ] &&
   [ "${HK_BLIND_ALLOW_NONEMPTY_RUN_ROOT:-0}" != "1" ] &&
   find "$FULL_RUN_ROOT" -mindepth 1 -print -quit | grep -q .; then
	echo "error: run root already exists and is not empty: $FULL_RUN_ROOT" >&2
	exit 1
fi
mkdir -p "$LOG_ROOT" "$EXPORT_ROOT" "$FULL_RUN_ROOT/scripts_snapshot"
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

SEG_SHAREC=${HK_BLIND_P9016_CONTACTS_SEG_SHAREC:-/sharec/zliu/CHARM/mESC/processed/P9016/2d_info/contacts.seg.gz}
SEG_ARCHIVE=${HK_BLIND_P9016_CONTACTS_SEG_ARCHIVE:-/shared/zliu/CHARM/CHARM_mesc/3_cellcycle/archieve/coverage/segs/P9016.sort.contacts_seg.gz}
SOURCES=${HK_BLIND_054_SOURCES:-sharec}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
MIN_MAPQ=${HK_BLIND_054_MIN_MAPQ:-20}
MIN_LEG_DIST=${HK_BLIND_054_MIN_LEG_DIST:-1000}
DUP_DIST=${HK_BLIND_054_DUP_DIST:-100}
READGROUP_MAX_SEG=${HK_BLIND_054_READGROUP_MAX_SEG:-10}

GIT_COMMIT=$(git rev-parse HEAD 2>/dev/null || echo NA)
GIT_DIRTY_COUNT=$(git status --porcelain | wc -l | awk '{print $1}')
git status --porcelain=v1 > "$LOG_ROOT/git_status.txt" || true
git diff --stat > "$LOG_ROOT/git_diff_stat.txt" || true
git diff > "$LOG_ROOT/git_diff.patch" || true
{
	echo "python	$(python --version 2>&1)"
	echo "conda_env	${CONDA_DEFAULT_ENV:-NA}"
	echo "path	$PATH"
} > "$LOG_ROOT/build_env.txt"

run_cmd python -m py_compile scripts/export_p9016_contacts_seg_readgroup_dedup.py
cp "$0" \
   scripts/export_p9016_contacts_seg_readgroup_dedup.py \
   scripts/p9016_publish_light_result.sh \
   "$FULL_RUN_ROOT/scripts_snapshot/"

{
	echo "key	value"
	echo "run_id	$RUN_ID"
	echo "full_run_root	$FULL_RUN_ROOT"
	echo "light_result_root	$LIGHT_RESULT_ROOT"
	echo "git_commit	$GIT_COMMIT"
	echo "git_dirty_count	$GIT_DIRTY_COUNT"
	echo "experiment_type	contacts_seg_readgroup_dedup_audit"
	echo "training_run	0"
	echo "bin_size_bp	$BIN_SIZE"
	echo "sources	$SOURCES"
	echo "sharec_seg	$SEG_SHAREC"
	echo "archive_seg	$SEG_ARCHIVE"
	echo "min_mapq	$MIN_MAPQ"
	echo "min_leg_dist	$MIN_LEG_DIST"
	echo "dup_dist	$DUP_DIST"
	echo "readgroup_max_seg	$READGROUP_MAX_SEG"
	echo "start_time	$(date -Is)"
} > "$FULL_RUN_ROOT/run_manifest.tsv"

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"

source_path() {
	case "$1" in
		sharec) echo "$SEG_SHAREC" ;;
		archive) echo "$SEG_ARCHIVE" ;;
		*) echo "unknown source: $1" >&2; return 2 ;;
	esac
}

export_one() {
	local source=$1
	local seg_path=$2
	local outdir="$EXPORT_ROOT/$source"
	mkdir -p "$outdir"
	if [ ! -s "$seg_path" ]; then
		echo "missing seg source for $source: $seg_path" >&2
		exit 2
	fi

	run_cmd python scripts/export_p9016_contacts_seg_readgroup_dedup.py \
		--seg "$seg_path" \
		--source-label "$source" \
		--outdir "$outdir" \
		--export-name hickit_like_dup100 \
		--pair-mode adjacent \
		--coordinate-mode boundary \
		--min-mapq "$MIN_MAPQ" \
		--min-leg-dist "$MIN_LEG_DIST" \
		--max-seg 3 \
		--dup-dist "$DUP_DIST" \
		--representatives-only \
		--write-annotation \
		> "$LOG_ROOT/$source.hickit_like_dup100.log" 2>&1

	run_cmd python scripts/export_p9016_contacts_seg_readgroup_dedup.py \
		--seg "$seg_path" \
		--source-label "$source" \
		--outdir "$outdir" \
		--export-name hickit_like_dup0 \
		--pair-mode adjacent \
		--coordinate-mode boundary \
		--min-mapq "$MIN_MAPQ" \
		--min-leg-dist "$MIN_LEG_DIST" \
		--max-seg 3 \
		--dup-dist 0 \
		--representatives-only \
		> "$LOG_ROOT/$source.hickit_like_dup0.log" 2>&1

	run_cmd python scripts/export_p9016_contacts_seg_readgroup_dedup.py \
		--seg "$seg_path" \
		--source-label "$source" \
		--outdir "$outdir" \
		--export-name readgroup_adjacent_maxseg${READGROUP_MAX_SEG}_annotated \
		--pair-mode adjacent \
		--coordinate-mode boundary \
		--min-mapq "$MIN_MAPQ" \
		--min-leg-dist "$MIN_LEG_DIST" \
		--max-seg "$READGROUP_MAX_SEG" \
		--dup-dist "$DUP_DIST" \
		--readgroup-cols \
		> "$LOG_ROOT/$source.readgroup_adjacent.log" 2>&1

	run_cmd python scripts/export_p9016_contacts_seg_readgroup_dedup.py \
		--seg "$seg_path" \
		--source-label "$source" \
		--outdir "$outdir" \
		--export-name readgroup_allcomb_maxseg${READGROUP_MAX_SEG}_annotated \
		--pair-mode allcomb \
		--coordinate-mode boundary \
		--min-mapq "$MIN_MAPQ" \
		--min-leg-dist "$MIN_LEG_DIST" \
		--max-seg "$READGROUP_MAX_SEG" \
		--dup-dist "$DUP_DIST" \
		--readgroup-cols \
		> "$LOG_ROOT/$source.readgroup_allcomb.log" 2>&1
}

for source in $SOURCES; do
	seg_path=$(source_path "$source")
	export_one "$source" "$seg_path"
done

log_command_line python - "$FULL_RUN_ROOT" "[inline summarize 054]"
python - "$FULL_RUN_ROOT" <<'PY'
import csv
import gzip
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])

def read_kv(path):
    data = {}
    if not path.exists():
        return data
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            if len(row) >= 2:
                data[row[0]] = row[1]
    return data

rows = []
for meta in sorted((root / "exports").glob("*/*.metadata.tsv")):
    data = read_kv(meta)
    pair_path = Path(data.get("out_pairs", ""))
    data["metadata_path"] = str(meta)
    data["pairs_exists"] = int(pair_path.exists())
    data["pairs_size_bytes"] = pair_path.stat().st_size if pair_path.exists() else "NA"
    data["pairs_sha256"] = "NA"
    rows.append(data)
fields = [
    "source_label","export_name","export_mode","hickit_like_semantics","readgroup_cols_written",
    "representatives_only","out_pairs","pairs_size_bytes","metadata_path","min_mapq","min_leg_dist",
    "max_seg","keep_over_max_seg_for_readgroup","dup_dist","n_readgroups_total","n_readgroups_ge2_raw",
    "n_readgroups_ge2_mapq","n_readgroups_used","n_readgroups_skipped_maxseg",
    "n_readgroups_multichrom_raw","n_readgroups_multichrom_used","n_readgroups_phase_labeled",
    "n_pairs_before_close_filter","n_pairs_after_close_filter","n_pairs_cis","n_pairs_trans",
    "n_pairs_before_dedup","n_pairs_after_dedup","n_pairs_removed_by_dedup","dedup_rate",
    "n_dup_clusters","cluster_size_p90","cluster_size_p99","cluster_size_max","n_pairs_written",
    "phase_labels_written","training_uses_phase_labels","training_uses_charm_or_reference",
]
with (root / "summary.tsv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", lineterminator="\n", fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: row.get(k, "NA") for k in fields})
with (root / "export_manifest.json").open("w") as fh:
    json.dump(rows, fh, indent=2)
PY

log_command_line python - "$FULL_RUN_ROOT" "[inline README 054]"
python - "$FULL_RUN_ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = list(csv.DictReader((root / "summary.tsv").open(), delimiter="\t"))
lines = [
    "# 054 P9016 contacts.seg Readgroup/Dedup Audit 1Mb",
    "",
    "This run exports blind-safe pairs from upstream contacts.seg sources and audits hickit-like pair-level duplicate semantics. No reconstruction is trained in 054.",
    "",
    "## Paths",
    "",
    f"- full result root: `{root}`",
    f"- summary: `{root / 'summary.tsv'}`",
    "",
    "## Boundary",
    "",
    "- Training is not run here.",
    "- Segment phase labels are dropped from all exported pairs.",
    "- hickit_like exports are reference/control inputs that mimic hickit adjacent-segment, boundary-coordinate, dup-dist pair handling.",
    "- readgroup exports preserve numeric read_group_id/read_seg columns for later readgroup-aware diagnostics.",
    "",
    "## Export Summary",
    "",
    "| source | export | pairs written | dedup rate | readgroup cols | hickit-like | path |",
    "|---|---|---:|---:|---:|---:|---|",
]
for row in rows:
    lines.append(
        "| {source} | {export} | {n} | {rate} | {rg} | {hk} | `{path}` |".format(
            source=row.get("source_label", "NA"),
            export=row.get("export_name", "NA"),
            n=row.get("n_pairs_written", "NA"),
            rate=row.get("dedup_rate", "NA"),
            rg=row.get("readgroup_cols_written", "NA"),
            hk=row.get("hickit_like_semantics", "NA"),
            path=row.get("out_pairs", "NA"),
        )
    )
lines.extend([
    "",
    "## Next Step",
    "",
    "055 should train/evaluate a minimal matrix using the approved baseline, hickit_like_dup100, readgroup_adjacent pairwise/joint, and readgroup_allcomb pairwise/joint.",
])
(root / "README.md").write_text("\n".join(lines) + "\n")
PY

log_command_line scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT"
scripts/p9016_publish_light_result.sh "$FULL_RUN_ROOT" "$LIGHT_RESULT_ROOT" > "$LOG_ROOT/publish.log" 2>&1

echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$LIGHT_RESULT_ROOT"
