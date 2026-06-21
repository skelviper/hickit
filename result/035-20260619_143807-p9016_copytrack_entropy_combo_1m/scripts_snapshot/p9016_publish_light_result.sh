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
	--include='gamma_delta.tsv'
	--include='sep_delta.tsv'
	--include='sample_delta.tsv'
	--include='synthetic_contact_summary.tsv'
	--include='trans_delta_summary.tsv'
	--include='trans_top1_delta_summary.tsv'
	--include='entropy_delta_summary.tsv'
	--include='copytrack_entropy_delta_summary.tsv'
	--include='readgroup_delta_summary.tsv'
	--include='readgroup_matched_delta_summary.tsv'
	--include='regression_audit.tsv'
	--include='dscale_semantics_audit.tsv'
	--include='dscale_semantics_headline.txt'
	--include='headline.txt'
	--include='copytrack_summary.tsv'
	--include='commands.log'
	--include='logs/*.log'
	--include='logs/*.txt'
	--include='logs/*.sh'
	--include='logs/*.patch'
	--include='contact_matrix_qc/*.tsv'
	--include='contact_matrix_qc/plots/*.png'
	--include='diagnostics/*/*.tsv'
	--include='diagnostics/*/*.md'
	--include='diagnostics/*/*.txt'
	--include='diagnostics/*/*.log'
	--include='synthetic_pairs/*.metadata.tsv'
	--include='synthetic_pairs/*.json'
	--include='generated_pairs/*.metadata.tsv'
	--include='generated_pairs/*.json'
	--include='outputs/*/*.manifest.tsv'
	--include='outputs/*/*.loop_diag.tsv'
	--include='outputs/*/*.force_class_diag.tsv'
	--include='outputs/*/*.sep_diag.tsv'
	--include='outputs/*/*/*.manifest.tsv'
	--include='outputs/*/*/*.loop_diag.tsv'
	--include='outputs/*/*/*.force_class_diag.tsv'
	--include='outputs/*/*/*.sep_diag.tsv'
	--include='eval/*/*.tsv'
	--include='eval/*/*.md'
	--include='eval/*/*.txt'
	--include='eval/*/*.log'
	--include='eval/*/*.json'
	--include='eval/*/plots/*.png'
	--include='eval/*/*/*.tsv'
	--include='eval/*/*/*.md'
	--include='eval/*/*/*.txt'
	--include='eval/*/*/*.log'
	--include='eval/*/*/*.json'
	--include='eval/*/*/plots/*.png'
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
	while IFS= read -r -d '' out_dir; do
		[ -e "$out_dir" ] || continue
		rel=${out_dir#"$FULL_RUN_ROOT/"}
		mkdir -p "$DEST/$rel"
		for f in \
			p9016_full.manifest.tsv \
			p9016_full.loop_diag.tsv \
			p9016_full.force_class_diag.tsv \
			p9016_full.sep_diag.tsv
		do
			if [ -s "$out_dir/$f" ]; then
				rsync -a --max-size=10m "$out_dir/$f" "$DEST/$rel/"
			fi
		done
	done < <(find -L "$FULL_RUN_ROOT/outputs" -type f -name p9016_full.manifest.tsv -printf '%h\0')
fi

RSYNC_PRINT=$(printf ' %q' "${RSYNC_ARGS[@]}" "$FULL_RUN_ROOT/" "$DEST/")
RSYNC_PRINT="rsync${RSYNC_PRINT}"
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
