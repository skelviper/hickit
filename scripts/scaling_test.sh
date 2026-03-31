#!/bin/bash
# Test scaling with different numbers of parallel cells

set -e

HICKIT="/mnt/ssd/zliu/run_charm/hickit_opt/hickit"
INPUT_DIR="/mnt/ssd/zliu/run_charm/charm_test/charm_test_13/result/impute/pairs"
WORK_DIR="/mnt/ssd/zliu/run_charm/hickit_opt/scaling_runs"
MPS_PIPE="$WORK_DIR/mps_pipe"

mkdir -p "$WORK_DIR"

# Get first N input files
get_inputs() {
    local n=$1
    find "$INPUT_DIR" -maxdepth 1 -name '*.pairs.gz' | sort | head -$n
}

# Run benchmark with N cells in parallel
run_benchmark() {
    local n_cells=$1
    echo "=== Testing with $n_cells cells ==="

    # Start MPS
    export CUDA_MPS_PIPE_DIRECTORY="$WORK_DIR"
    nvidia-cuda-mps-control -d

    local inputs=($(get_inputs $n_cells))
    local pids=()
    local start_time=$(date +%s.%N)

    # Launch all cells in parallel
    for input in "${inputs[@]}"; do
        local cell_name=$(basename "$input" .impute.pairs.gz)
        local out_dir="$WORK_DIR/${n_cells}cells/$cell_name"
        mkdir -p "$out_dir"

        (/usr/bin/time -f "%e" $HICKIT -s1 -M -i "$input" \
            -Sr1m -c1 -r10m -c2 -b4m -b1m -b200k -D5 -b50k -D5 -b20k \
            -O "$out_dir/20k.3dg" 2>&1 | tail -1) > "$out_dir/time.txt" &
        pids+=($!)
    done

    # Wait for all to complete
    for pid in "${pids[@]}"; do
        wait $pid
    done

    local end_time=$(date +%s.%N)
    local makespan=$(echo "$end_time - $start_time" | bc)

    # Collect individual times
    local sum_time=0
    for input in "${inputs[@]}"; do
        local cell_name=$(basename "$input" .impute.pairs.gz)
        local time_file="$WORK_DIR/${n_cells}cells/$cell_name/time.txt"
        if [ -f "$time_file" ]; then
            local cell_time=$(cat "$time_file")
            sum_time=$(echo "$sum_time + $cell_time" | bc)
        fi
    done

    local avg_time=$(echo "scale=2; $sum_time / $n_cells" | bc)
    local parallelism=$(echo "scale=2; $sum_time / $makespan" | bc)
    local efficiency=$(echo "scale=1; 100 * $parallelism / $n_cells" | bc)

    echo "Cells: $n_cells"
    echo "Makespan: ${makespan}s"
    echo "Sum of times: ${sum_time}s"
    echo "Avg time: ${avg_time}s"
    echo "Effective parallelism: ${parallelism}x"
    echo "Parallel efficiency: ${efficiency}%"
    echo ""

    # Stop MPS
    echo quit | nvidia-cuda-mps-control
    sleep 2
}

# Test different cell counts
for n in 1 2 4 6 8 12; do
    run_benchmark $n
done

echo "=== Scaling Test Complete ==="
