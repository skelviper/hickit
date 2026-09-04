#!/usr/bin/env bash
set -euo pipefail

if [[ "${CONDA_DEFAULT_ENV:-}" != "analysis" ]]; then
    echo "run_phase_exp057.sh requires: conda activate analysis" >&2
    exit 2
fi

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/hickit-phase-mpl-${USER:-user}}"
mkdir -p "$MPLCONFIGDIR"
exec python "$repo_dir/scripts/run_phase_exp057.py" "$@"
