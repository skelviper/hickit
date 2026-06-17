#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin

ROOT=${HK_BLIND_P9016_D_SCALE_GAMMA_SMOKE_ROOT:-/tmp/hk_blind_p9016_dscale_gamma_smoke}
rm -rf "$ROOT"
mkdir -p "$ROOT"

run_one() {
	local config=$1
	local mode=$2
	local gamma=$3
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
		HK_BLIND_P9016_D_SCALE_POSTERIOR_GAMMA="$gamma" \
		HK_BLIND_P9016_D_SCALE_EPS_COUNT=1e-6 \
		./run_blind_p9016_minimal.bin > "$ROOT/$config.train.log" 2>&1
	./audit_blind_p9016_full_cpu_output.bin "$ROOT/$config" > "$ROOT/$config.audit.log" 2>&1
}

require_manifest_value() {
	local manifest=$1
	local key=$2
	local expected=$3
	awk -F'\t' -v key="$key" -v expected="$expected" \
		'$1 == key && $2 == expected { found=1 } END { exit found?0:1 }' "$manifest"
}

require_manifest_float() {
	local manifest=$1
	local key=$2
	local expected=$3
	awk -F'\t' -v key="$key" -v expected="$expected" \
		'function abs(x){return x<0?-x:x} $1 == key && abs(($2+0)-expected) < 1e-8 { found=1 } END { exit found?0:1 }' "$manifest"
}

run_one smoke_raw_count raw_count 1
run_one smoke_pcgamma0 posterior_count 0
run_one smoke_pcgamma0p5 posterior_count 0.5
run_one smoke_pcgamma1 posterior_count 1

require_manifest_value "$ROOT/smoke_raw_count/p9016_full.manifest.tsv" d_scale_mode raw_count
require_manifest_value "$ROOT/smoke_raw_count/p9016_full.manifest.tsv" dscale_effective_count_formula n_raw
require_manifest_value "$ROOT/smoke_raw_count/p9016_full.manifest.tsv" dscale_probability_weighted 0

require_manifest_value "$ROOT/smoke_pcgamma0/p9016_full.manifest.tsv" d_scale_mode posterior_count
require_manifest_value "$ROOT/smoke_pcgamma0/p9016_full.manifest.tsv" dscale_effective_count_formula 'n_raw*posterior_prob^0'
require_manifest_float "$ROOT/smoke_pcgamma0/p9016_full.manifest.tsv" d_scale_posterior_gamma 0
require_manifest_value "$ROOT/smoke_pcgamma0/p9016_full.manifest.tsv" dscale_probability_weighted 0

require_manifest_value "$ROOT/smoke_pcgamma0p5/p9016_full.manifest.tsv" d_scale_mode posterior_count
require_manifest_value "$ROOT/smoke_pcgamma0p5/p9016_full.manifest.tsv" dscale_effective_count_formula 'n_raw*posterior_prob^0.5'
require_manifest_float "$ROOT/smoke_pcgamma0p5/p9016_full.manifest.tsv" d_scale_posterior_gamma 0.5
require_manifest_value "$ROOT/smoke_pcgamma0p5/p9016_full.manifest.tsv" dscale_probability_weighted 1

require_manifest_value "$ROOT/smoke_pcgamma1/p9016_full.manifest.tsv" d_scale_mode posterior_count
require_manifest_value "$ROOT/smoke_pcgamma1/p9016_full.manifest.tsv" dscale_effective_count_formula 'n_raw*posterior_prob^1'
require_manifest_float "$ROOT/smoke_pcgamma1/p9016_full.manifest.tsv" d_scale_posterior_gamma 1
require_manifest_value "$ROOT/smoke_pcgamma1/p9016_full.manifest.tsv" dscale_probability_weighted 1

echo "DSCALE_GAMMA_SMOKE_OK root=$ROOT"
