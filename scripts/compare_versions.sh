#!/bin/bash
# Compare optimized vs original hickit performance and correctness

set -e

HICKIT_OPT="/mnt/ssd/zliu/run_charm/hickit_opt/hickit"
HICKIT_ORIG="/mnt/ssd/zliu/run_charm/hickit/hickit"
INPUT_DIR="/mnt/ssd/zliu/run_charm/charm_test/charm_test_13/result/impute/pairs"
WORK_DIR="/mnt/ssd/zliu/run_charm/hickit_opt/comparison"

mkdir -p "$WORK_DIR"

echo "=== Comparing Optimized vs Original Hickit ==="
echo ""

# Test with 6 cells
test_version() {
    local version=$1
    local hickit_bin=$2
    local n_cells=$3

    echo "=== Testing $version with $n_cells cells ==="

    # Start MPS
    export CUDA_MPS_PIPE_DIRECTORY="$WORK_DIR/mps_${version}"
    mkdir -p "$CUDA_MPS_PIPE_DIRECTORY"
    nvidia-cuda-mps-control -d 2>/dev/null || true
    sleep 1

    # Get inputs
    local inputs=($(find "$INPUT_DIR" -maxdepth 1 -name '*.pairs.gz' | sort | head -$n_cells))

    local pids=()
    local start=$(date +%s.%N)

    # Launch all cells
    for input in "${inputs[@]}"; do
        $hickit_bin -s1 -M -i "$input" -Sr1m -c1 -r10m -c2 -b4m -b1m -b200k -D5 -b50k -D5 -b20k -O /dev/null 2>&1 | grep "CPU time" &
        pids+=($!)
    done

    # Wait for completion
    for pid in "${pids[@]}"; do
        wait $pid
    done

    local end=$(date +%s.%N)
    local makespan=$(echo "$end - $start" | bc)

    echo "Makespan: ${makespan}s"
    echo ""

    # Stop MPS
    echo quit | nvidia-cuda-mps-control 2>/dev/null || true
    sleep 2

    echo "$makespan"
}

# Performance comparison
echo "### Performance Comparison ###"
echo ""

opt_6=$(test_version "optimized" "$HICKIT_OPT" 6)
orig_6=$(test_version "original" "$HICKIT_ORIG" 6)

opt_12=$(test_version "optimized" "$HICKIT_OPT" 12)
orig_12=$(test_version "original" "$HICKIT_ORIG" 12)

echo "=== Performance Summary ==="
echo "6 cells:"
echo "  Optimized: ${opt_6}s"
echo "  Original:  ${orig_6}s"
speedup_6=$(echo "scale=2; $orig_6 / $opt_6" | bc)
echo "  Speedup:   ${speedup_6}x"
echo ""
echo "12 cells:"
echo "  Optimized: ${opt_12}s"
echo "  Original:  ${orig_12}s"
speedup_12=$(echo "scale=2; $orig_12 / $opt_12" | bc)
echo "  Speedup:   ${speedup_12}x"
echo ""

# Correctness comparison
echo "### Correctness Comparison ###"
echo ""

test_input=$(find "$INPUT_DIR" -maxdepth 1 -name '*.pairs.gz' | head -1)
test_name=$(basename "$test_input" .impute.pairs.gz)

echo "Testing with: $test_name"

$HICKIT_OPT -s1 -M -i "$test_input" -Sr1m -c1 -r10m -c2 -b4m -b1m -b200k -D5 -b50k -D5 -b20k -O "$WORK_DIR/opt_${test_name}.3dg" 2>&1 | tail -5
$HICKIT_ORIG -s1 -M -i "$test_input" -Sr1m -c1 -r10m -c2 -b4m -b1m -b200k -D5 -b50k -D5 -b20k -O "$WORK_DIR/orig_${test_name}.3dg" 2>&1 | tail -5

echo ""
echo "Comparing output files..."
if diff -q "$WORK_DIR/opt_${test_name}.3dg" "$WORK_DIR/orig_${test_name}.3dg" > /dev/null 2>&1; then
    echo "✓ Output files are identical"
else
    echo "✗ Output files differ"
    echo "Checking coordinate differences..."

    # Compare line counts
    opt_lines=$(wc -l < "$WORK_DIR/opt_${test_name}.3dg")
    orig_lines=$(wc -l < "$WORK_DIR/orig_${test_name}.3dg")
    echo "  Optimized: $opt_lines lines"
    echo "  Original:  $orig_lines lines"

    # Sample first few coordinates
    echo ""
    echo "First 5 lines of optimized:"
    head -5 "$WORK_DIR/opt_${test_name}.3dg"
    echo ""
    echo "First 5 lines of original:"
    head -5 "$WORK_DIR/orig_${test_name}.3dg"
fi

echo ""
echo "=== Comparison Complete ==="
