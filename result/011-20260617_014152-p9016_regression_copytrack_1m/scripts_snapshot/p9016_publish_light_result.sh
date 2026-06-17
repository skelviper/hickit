#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
	echo "usage: scripts/p9016_publish_light_result.sh <full_run_root> [dest_result_root]" >&2
	exit 2
fi

FULL_RUN_ROOT=$(cd "$1" && pwd)
DEST=${2:-"$PWD/result/$(basename "$FULL_RUN_ROOT")"}
mkdir -p "$(dirname "$DEST")"
rm -rf "$DEST"

RSYNC_ARGS=(
	-a
	--prune-empty-dirs
	--max-size=10m
	--include='*/'
	--include='README.md'
	--include='LIGHT_COPY_NOTE.md'
	--include='run_manifest.tsv'
	--include='summary.tsv'
	--include='regression_audit.tsv'
	--include='copytrack_summary.tsv'
	--include='commands.log'
	--include='logs/*.log'
	--include='logs/*.txt'
	--include='logs/*.sh'
	--include='outputs/*/*.manifest.tsv'
	--include='outputs/*/*.loop_diag.tsv'
	--include='outputs/*/*.force_class_diag.tsv'
	--include='outputs/*/*.sep_diag.tsv'
	--include='eval/*/*.tsv'
	--include='eval/*/*.md'
	--include='eval/*/*.txt'
	--include='eval/*/*.log'
	--include='eval/*/plots/*.png'
	--include='scripts_snapshot/*'
	--exclude='*.3dg'
	--exclude='*.3dg.gz'
	--exclude='*.pairs'
	--exclude='*.pairs.gz'
	--exclude='*.bam'
	--exclude='*.cram'
	--exclude='*.cool'
	--exclude='*.mcool'
	--exclude='*.npy'
	--exclude='*.npz'
	--exclude='*.bin'
	--exclude='*.o'
	--exclude='*.so'
	--exclude='*.a'
	--exclude='*.tmp'
	--exclude='core.*'
	--exclude='*coord*.tsv'
	--exclude='*posterior*.tsv'
	--exclude='*graph*.tsv'
	--exclude='*'
)

rsync "${RSYNC_ARGS[@]}" "$FULL_RUN_ROOT/" "$DEST/"

RSYNC_PRINT=$(printf 'rsync %q ' "${RSYNC_ARGS[@]}" "$FULL_RUN_ROOT/" "$DEST/")
cat > "$DEST/LIGHT_COPY_NOTE.md" <<EOF
# Lightweight Result Copy

- source full_run_root: \`$FULL_RUN_ROOT\`
- destination: \`$DEST\`
- timestamp: \`$(date -Is)\`
- rsync command: \`$RSYNC_PRINT\`
- exclusion policy: large coordinates, posterior dumps, raw contacts, binary/build artifacts, and files above 10 MB are excluded.
- warning: large files remain only in \`$FULL_RUN_ROOT\`.
EOF

echo "FULL_RESULT_ROOT=$FULL_RUN_ROOT"
echo "LIGHT_RESULT_ROOT=$DEST"
