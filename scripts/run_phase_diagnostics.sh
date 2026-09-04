#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
output_root="$repo_dir/results/phase_diagnostics"
seed=17
threads=1
mode="quick"
run_audit=0
run_057=0
run_058=0
run_059=0
run_comparison=0
run_tests=0

usage() {
    printf '%s\n' \
        "usage: $0 [--audit] [--exp057] [--exp058] [--exp059] [--model-comparison]" \
        "          [--quick|--full] [--seed N] [--threads N] [--output-dir DIR]"
}

while (($#)); do
    case "$1" in
        --audit) run_audit=1 ;;
        --exp057) run_057=1 ;;
        --exp058) run_058=1 ;;
        --exp059) run_059=1 ;;
        --model-comparison) run_comparison=1 ;;
        --quick)
            mode="quick"
            run_tests=1
            run_audit=1
            run_057=1
            run_058=1
            run_059=1
            run_comparison=1
            ;;
        --full)
            mode="full"
            run_tests=1
            run_audit=1
            run_057=1
            run_058=1
            run_059=1
            run_comparison=1
            ;;
        --seed)
            shift
            seed="${1:?--seed requires an integer}"
            ;;
        --threads)
            shift
            threads="${1:?--threads requires an integer}"
            ;;
        --output-dir)
            shift
            output_root="${1:?--output-dir requires a path}"
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            printf 'unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
    shift
done

if ! [[ "$seed" =~ ^[0-9]+$ && "$threads" =~ ^[1-9][0-9]*$ ]]; then
    printf '%s\n' "--seed must be non-negative and --threads must be positive" >&2
    exit 2
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate analysis
export OMP_NUM_THREADS="$threads"
export MKL_NUM_THREADS="$threads"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/hickit-phase-diagnostics-matplotlib}"
mkdir -p "$MPLCONFIGDIR" "$output_root"

quick_arg=()
if [[ "$mode" == "quick" ]]; then
    quick_arg=(--quick)
fi

if ((run_tests)); then
    test_timestamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    test_commit="$(git -C "$repo_dir" rev-parse HEAD)"
    test_dirty="$(git -C "$repo_dir" status --porcelain | wc -l | tr -d ' ')"
    pytest_log="$output_root/test_python.log"
    ctest_log="$output_root/test_c.log"
    set +e
    pytest -q "$repo_dir/tests/test_phase_diagnostics.py" 2>&1 | tee "$pytest_log"
    pytest_status="${PIPESTATUS[0]}"
    make -C "$repo_dir" --no-print-directory test_blind_estep_scores 2>&1 | tee "$ctest_log"
    ctest_status="${PIPESTATUS[0]}"
    set -e
    {
        printf 'timestamp_utc\tgit_commit\tgit_dirty_file_count\tconda_environment\tcommand\texit_code\tsummary\tlog\n'
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$test_timestamp" "$test_commit" "$test_dirty" "${CONDA_DEFAULT_ENV:-}" \
            "pytest -q tests/test_phase_diagnostics.py" "$pytest_status" \
            "$(tail -n 1 "$pytest_log" | tr '\t' ' ')" "$pytest_log"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$test_timestamp" "$test_commit" "$test_dirty" "${CONDA_DEFAULT_ENV:-}" \
            "make test_blind_estep_scores" "$ctest_status" \
            "$(tail -n 1 "$ctest_log" | tr '\t' ' ')" "$ctest_log"
    } > "$output_root/test_receipt.tsv"
    if ((pytest_status != 0 || ctest_status != 0)); then
        exit 1
    fi
fi
if ((run_audit)); then
    python "$repo_dir/scripts/run_phase_audit.py" --output-dir "$output_root"
fi
if ((run_057)); then
    "$repo_dir/scripts/run_phase_exp057.sh" "${quick_arg[@]}" --seed "$seed" --output-dir "$output_root/057"
fi
if ((run_058)); then
    "$repo_dir/scripts/run_phase_exp058.sh" "${quick_arg[@]}" --seed "$seed" --output-dir "$output_root/058"
fi
if ((run_059)); then
    "$repo_dir/scripts/run_phase_exp059.sh" "${quick_arg[@]}" --seed "$seed" --output-dir "$output_root/059"
fi
if ((run_comparison)); then
    python "$repo_dir/scripts/run_phase_model_comparison.py" "${quick_arg[@]}" \
        --seed "$seed" --results-root "$output_root" --output-dir "$output_root/model_comparison"
fi
