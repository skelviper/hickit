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
	--include='dscale_semantics_audit.tsv'
	--include='dscale_semantics_headline.txt'
	--include='copytrack_summary.tsv'
	--include='commands.log'
	--include='logs/*.log'
	--include='logs/*.txt'
	--include='logs/*.patch'
	--include='logs/*.sh'
	--include='outputs/*/*.manifest.tsv'
	--include='outputs/*/*.loop_diag.tsv'
	--include='outputs/*/*.force_class_diag.tsv'
	--include='outputs/*/*.sep_diag.tsv'
	--include='eval/*/*.tsv'
	--include='eval/*/*.md'
	--include='eval/*/*.txt'
	--include='eval/*/*.log'
	--include='eval/*/*.json'
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

if [ -d "$FULL_RUN_ROOT/outputs" ]; then
	for out_dir in "$FULL_RUN_ROOT"/outputs/*; do
		[ -e "$out_dir" ] || continue
		config=$(basename "$out_dir")
		mkdir -p "$DEST/outputs/$config"
		for f in \
			p9016_full.manifest.tsv \
			p9016_full.loop_diag.tsv \
			p9016_full.force_class_diag.tsv \
			p9016_full.sep_diag.tsv
		do
			if [ -s "$out_dir/$f" ]; then
				rsync -a --max-size=10m "$out_dir/$f" "$DEST/outputs/$config/"
			fi
		done
	done
fi

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
if [ -s "$FULL_RUN_ROOT/summary.tsv" ]; then
	echo "SUMMARY_TSV=$FULL_RUN_ROOT/summary.tsv"
fi
if [ -s "$FULL_RUN_ROOT/dscale_semantics_audit.tsv" ]; then
	echo "DSCALE_SEMANTICS_AUDIT_TSV=$FULL_RUN_ROOT/dscale_semantics_audit.tsv"
fi
if [ -s "$FULL_RUN_ROOT/dscale_semantics_headline.txt" ]; then
	echo "DSCALE_SEMANTICS_HEADLINE=$(cat "$FULL_RUN_ROOT/dscale_semantics_headline.txt")"
fi
