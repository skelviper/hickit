#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin

ROOT=${HK_BLIND_P9016_COPYTRACK_SMOKE_ROOT:-/tmp/hk_blind_p9016_copytrack_smoke}
rm -rf "$ROOT"
mkdir -p "$ROOT"

env -i PATH="$PATH" LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}" \
	HK_BLIND_SAMPLE=P9016 \
	HK_BLIND_P9016_PAIRS=testdata/p9016_blind_smoke.pairs \
	HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS=1 \
	HK_BLIND_P9016_OUTPUT_ROOT="$ROOT" \
	HK_BLIND_P9016_CONFIG_NAME=smoke_copytrack_lambda0p01 \
	HK_BLIND_P9016_BIN_SIZE_BP=1000000 \
	HK_BLIND_P9016_MINIMAL_N_ITER=1 \
	HK_BLIND_P9016_MINIMAL_RELAX_STEPS=1 \
	HK_BLIND_P9016_RELAX_BACKEND=cpu \
	HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split \
	HK_BLIND_P9016_INIT_SEED=17 \
	HK_BLIND_P9016_INIT_EPS=0.5 \
	HK_BLIND_P9016_INIT_NOISE_SCALE=0.0 \
	HK_BLIND_P9016_MIN_SEP_UNIT=0.0 \
	HK_BLIND_P9016_LAMBDA_SEP=0.0 \
	HK_BLIND_P9016_LAMBDA_COPYTRACK=0.01 \
	HK_BLIND_P9016_D_SCALE_MODE=posterior_count \
	HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA=1 \
	HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6 \
	./run_blind_p9016_minimal.bin > "$ROOT/train.log" 2>&1

./audit_blind_p9016_full_cpu_output.bin "$ROOT/smoke_copytrack_lambda0p01" > "$ROOT/audit.log" 2>&1

MANIFEST="$ROOT/smoke_copytrack_lambda0p01/p9016_full.manifest.tsv"
LOOP="$ROOT/smoke_copytrack_lambda0p01/p9016_full.loop_diag.tsv"
FORCE="$ROOT/smoke_copytrack_lambda0p01/p9016_full.force_class_diag.tsv"

awk -F'\t' '$1 == "lambda_copytrack" && (($2 + 0) > 0.009999 && ($2 + 0) < 0.010001) { found=1 } END { exit found?0:1 }' "$MANIFEST"
awk -F'\t' '$1 == "copytrack_prior_mode" && $2 == "homolog_vector_continuity" { found=1 } END { exit found?0:1 }' "$MANIFEST"
awk -F'\t' '$1 == "copytrack_uses_phase_labels" && $2 == "0" { found=1 } END { exit found?0:1 }' "$MANIFEST"
awk -F'\t' '$1 == "copytrack_uses_charm_or_reference" && $2 == "0" { found=1 } END { exit found?0:1 }' "$MANIFEST"
head -n 1 "$LOOP" | grep -q 'final_copytrack_force_l1'
head -n 1 "$FORCE" | grep -q 'copytrack_force_l1'

echo "COPYTRACK_SMOKE_OK root=$ROOT"
