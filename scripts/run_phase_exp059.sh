#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate analysis
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/hickit-phase-diagnostics-matplotlib}"
mkdir -p "$MPLCONFIGDIR"
exec python "$repo_dir/scripts/run_phase_exp059.py" "$@"
