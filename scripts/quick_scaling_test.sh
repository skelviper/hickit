#!/bin/bash
# Quick test: 6 cells vs 12 cells

set -e

HICKIT="/mnt/ssd/zliu/run_charm/hickit_opt/hickit"
INPUT_DIR="/mnt/ssd/zliu/run_charm/charm_test/charm_test_13/result/impute/pairs"
WORK_DIR="/mnt/ssd/zliu/run_charm/hickit_opt/quick_scaling"

mkdir -p "$WORK_DIR"

run_test() {
    local n_cells=$1
    echo "=== Testing $n_cells cells ==="

    # Start MPS
    export CUDA_MPS_PIPE_DIRECTORY="$WORK_DIR"
    nvidia-cuda-mps-control -d
    sleep 1

    # Get inputs
    local inputs=($(find "$INPUT_DIR" -maxdepth 1 -name '*.pairs.gz' | sort | head -$n_cells))

    # Start GPU monitoring in background
    nvidia-smi dmon -s u -c 100 > "$WORK_DIR/gpu_util_${n_cells}cells.log" 2>&1 &
    local monitor_pid=$!

    local pids=()
    local start=$(date +%s.%N)

    # Launch all cells
    for input in "${inputs[@]}"; do
        $HICKIT -s1 -M -i "$input" -Sr1m -c1 -r10m -c2 -b4m -b1m -b200k -D5 -b50k -D5 -b20k -O /dev/null 2>&1 | grep "CPU time" &
        pids+=($!)
    done

    # Wait for completion
    for pid in "${pids[@]}"; do
        wait $pid
    done

    local end=$(date +%s.%N)
    local makespan=$(echo "$end - $start" | bc)

    # Stop monitoring
    kill $monitor_pid 2>/dev/null || true

    echo "Makespan: ${makespan}s"
    echo "Expected ideal time (if linear): $(echo "scale=2; $makespan / $n_cells" | bc)s per cell"
    echo ""

    # Stop MPS
    echo quit | nvidia-cuda-mps-control
    sleep 2
}

# Test 6 cells
run_test 6

# Test 12 cells
run_test 12

echo "=== Results Summary ==="
echo "Check GPU utilization logs:"
echo "  6 cells:  $WORK_DIR/gpu_util_6cells.log"
echo "  12 cells: $WORK_DIR/gpu_util_12cells.log"
