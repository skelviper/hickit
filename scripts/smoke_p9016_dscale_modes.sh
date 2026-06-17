#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin

TEST_RES_ROOT=${HK_BLIND_TEST_RES_ROOT:-"$REPO_ROOT/test_res"}
ROOT=${HK_BLIND_P9016_D_SCALE_SMOKE_ROOT:-"$TEST_RES_ROOT/smoke/hk_blind_p9016_dscale_modes_smoke"}
rm -rf "$ROOT"
mkdir -p "$ROOT"

run_one() {
	local config=$1
	local mode=$2
	env -i PATH="$PATH" LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}" \
		HK_BLIND_P9016_PAIRS=testdata/p9016_blind_smoke.pairs \
		HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS=1 \
		HK_BLIND_P9016_OUTPUT_ROOT="$ROOT" \
		HK_BLIND_P9016_CONFIG_NAME="$config" \
		HK_BLIND_P9016_BIN_SIZE_BP=1000000 \
		HK_BLIND_P9016_MINIMAL_N_ITER=1 \
		HK_BLIND_P9016_MINIMAL_RELAX_STEPS=1 \
		HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split \
		HK_BLIND_P9016_INIT_SEED=17 \
		HK_BLIND_P9016_INIT_EPS=0.5 \
		HK_BLIND_P9016_INIT_NOISE_SCALE=0.0 \
		HK_BLIND_P9016_MIN_SEP_UNIT=0.0 \
		HK_BLIND_P9016_LAMBDA_SEP=0.0 \
		HK_BLIND_P9016_D_SCALE_MODE="$mode" \
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6 \
		./run_blind_p9016_minimal.bin > "$ROOT/$config.train.log" 2>&1
	./audit_blind_p9016_full_cpu_output.bin "$ROOT/$config" > "$ROOT/$config.audit.log" 2>&1
}

run_one smoke_raw_count raw_count
run_one smoke_posterior_count posterior_count

RAW_MANIFEST="$ROOT/smoke_raw_count/p9016_full.manifest.tsv"
POST_MANIFEST="$ROOT/smoke_posterior_count/p9016_full.manifest.tsv"

awk -F'\t' '$1=="dscale_effective_count_formula" && $2=="n_raw" { found=1 } END { exit found?0:1 }' "$RAW_MANIFEST"
awk -F'\t' '$1=="dscale_probability_weighted" && $2=="0" { found=1 } END { exit found?0:1 }' "$RAW_MANIFEST"
awk -F'\t' '$1=="d_scale_mode" && $2=="raw_count" { found=1 } END { exit found?0:1 }' "$RAW_MANIFEST"

awk -F'\t' '$1=="dscale_effective_count_formula" && $2=="n_raw*posterior_prob^1" { found=1 } END { exit found?0:1 }' "$POST_MANIFEST"
awk -F'\t' '$1=="dscale_probability_weighted" && $2=="1" { found=1 } END { exit found?0:1 }' "$POST_MANIFEST"
awk -F'\t' '$1=="d_scale_mode" && $2=="posterior_count" { found=1 } END { exit found?0:1 }' "$POST_MANIFEST"
awk -F'\t' '$1=="d_scale_posterior_gamma" && ($2+0)==1 { found=1 } END { exit found?0:1 }' "$POST_MANIFEST"
awk -F'\t' '$1=="legacy_expected_count_alias_used" && $2=="0" { found=1 } END { exit found?0:1 }' "$POST_MANIFEST"

echo "DSCALE_MODE_SMOKE_OK root=$ROOT"
