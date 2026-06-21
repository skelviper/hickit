#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
	echo "usage: scripts/p9016_common_eval.sh <config_output_dir> <eval_output_dir>" >&2
	exit 2
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

CONFIG_OUTPUT_DIR=$(cd "$1" && pwd)
mkdir -p "$2"
EVAL_OUTPUT_DIR=$(cd "$2" && pwd)

PAIRS=${HK_BLIND_P9016_EVAL_PAIRS:-/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz}
REFERENCE_3DG=${HK_BLIND_P9016_EVAL_TDG:-/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz}
BIN_SIZE=${HK_BLIND_P9016_BIN_SIZE_BP:-1000000}
SAMPLE=${HK_BLIND_SAMPLE:-${HK_BLIND_P9016_SAMPLE:-P9016}}
RECONSTRUCTION="$CONFIG_OUTPUT_DIR/p9016_full.coords.tsv"
POSTERIOR="$CONFIG_OUTPUT_DIR/p9016_full.bpair_posterior.tsv"
MANIFEST="$CONFIG_OUTPUT_DIR/p9016_full.manifest.tsv"

cmd=(
	python "$REPO_ROOT/eval/evaluate_p9016_baseline.py"
	--pairs "$PAIRS"
	--reference-3dg "$REFERENCE_3DG"
	--reconstruction "$RECONSTRUCTION"
	--posterior "$POSTERIOR"
	--train-manifest "$MANIFEST"
	--outdir "$EVAL_OUTPUT_DIR"
	--bin-size "$BIN_SIZE"
	--label "$(basename "$CONFIG_OUTPUT_DIR")"
	--expected-sample "$SAMPLE"
)
if [ "${HK_BLIND_ALLOW_REFERENCE_DERIVED_TRAINING:-0}" = "1" ]; then
	cmd+=(--allow-reference-derived-training)
fi
{
	for i in "${!cmd[@]}"; do
		if [ "$i" -gt 0 ]; then
			printf ' '
		fi
		printf '%q' "${cmd[$i]}"
	done
	printf '\n'
} > "$EVAL_OUTPUT_DIR/eval_command.txt"

if [ ! -s "$PAIRS" ] || [ ! -s "$REFERENCE_3DG" ] || [ ! -s "$RECONSTRUCTION" ] || [ ! -s "$POSTERIOR" ] || [ ! -s "$MANIFEST" ]; then
	{
		echo "status	SKIPPED_EVAL_INPUT_MISSING"
		echo "pairs	$PAIRS"
		echo "reference_3dg	$REFERENCE_3DG"
		echo "sample	$SAMPLE"
		echo "reconstruction	$RECONSTRUCTION"
		echo "posterior	$POSTERIOR"
		echo "manifest	$MANIFEST"
	} > "$EVAL_OUTPUT_DIR/eval_status.tsv"
	echo "SKIPPED_EVAL_INPUT_MISSING" > "$EVAL_OUTPUT_DIR/eval.log"
	exit 0
fi

"${cmd[@]}" > "$EVAL_OUTPUT_DIR/eval.log" 2>&1

if [ ! -s "$EVAL_OUTPUT_DIR/summary.tsv" ]; then
	echo "eval finished but summary.tsv is missing: $EVAL_OUTPUT_DIR" >&2
	exit 1
fi
