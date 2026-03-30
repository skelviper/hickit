#!/bin/bash
set -euo pipefail

HICKIT_BIN="${HICKIT_BIN:-/mnt/ssd/zliu/run_charm/hickit_opt/hickit}"
INPUT_DIR="${INPUT_DIR:-/mnt/ssd/zliu/run_charm/charm_test/charm_test_13/result/impute/pairs}"
OUTDIR="${OUTDIR:-/mnt/ssd/zliu/run_charm/hickit_opt/multi_cell_runs}"
JOBS="${JOBS:-4}"
MAX_CELLS="${MAX_CELLS:-$JOBS}"
USE_MPS="${USE_MPS:-0}"
MPS_PIPE_DIR="${MPS_PIPE_DIR:-$OUTDIR/mps_pipe}"
MPS_LOG_DIR="${MPS_LOG_DIR:-$OUTDIR/mps_log}"

if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || [ "$JOBS" -lt 1 ]; then
	echo "JOBS must be a positive integer" >&2
	exit 1
fi
if ! [[ "$MAX_CELLS" =~ ^[0-9]+$ ]] || [ "$MAX_CELLS" -lt 1 ]; then
	echo "MAX_CELLS must be a positive integer" >&2
	exit 1
fi
if [ ! -x "$HICKIT_BIN" ]; then
	echo "HICKIT_BIN is not executable: $HICKIT_BIN" >&2
	exit 1
fi
if [ ! -d "$INPUT_DIR" ]; then
	echo "INPUT_DIR does not exist: $INPUT_DIR" >&2
	exit 1
fi

mapfile -t INPUTS < <(find "$INPUT_DIR" -maxdepth 1 -name '*.pairs.gz' | sort | head -n "$MAX_CELLS")
if [ "${#INPUTS[@]}" -eq 0 ]; then
	echo "No .pairs.gz files found in $INPUT_DIR" >&2
	exit 1
fi

mkdir -p "$OUTDIR"

mps_started=0
cleanup() {
	if [ "$mps_started" -eq 1 ]; then
		echo quit | CUDA_MPS_PIPE_DIRECTORY="$MPS_PIPE_DIR" nvidia-cuda-mps-control >/dev/null 2>&1 || true
	fi
}
trap cleanup EXIT

if [ "$USE_MPS" = "1" ]; then
	if ! command -v nvidia-cuda-mps-control >/dev/null 2>&1; then
		echo "USE_MPS=1 requested, but nvidia-cuda-mps-control is not available" >&2
		exit 1
	fi
	rm -rf "$MPS_PIPE_DIR" "$MPS_LOG_DIR"
	mkdir -p "$MPS_PIPE_DIR" "$MPS_LOG_DIR"
	export CUDA_MPS_PIPE_DIRECTORY="$MPS_PIPE_DIR"
	export CUDA_MPS_LOG_DIRECTORY="$MPS_LOG_DIR"
	nvidia-cuda-mps-control -d
	mps_started=1
	echo "CUDA MPS enabled via $MPS_PIPE_DIR"
fi

run_one_cell() {
	local input="$1"
	local stem
	local cell_dir

	stem="$(basename "$input" .impute.pairs.gz)"
	cell_dir="$OUTDIR/$stem"
	mkdir -p "$cell_dir"

	/usr/bin/time -f "%e" -o "$cell_dir/wall_seconds.txt" \
		"$HICKIT_BIN" -s1 -M \
		-i "$input" -Sr1m -c1 -r10m -c2 \
		-b4m -b1m -O "$cell_dir/1m.3dg" \
		-b200k -O "$cell_dir/200k.3dg" \
		-D5 -b50k -O "$cell_dir/50k.3dg" \
		-D5 -b20k -O "$cell_dir/20k.3dg" \
		>"$cell_dir/stdout.log" 2>"$cell_dir/stderr.log"
}

start_ts="$(date +%s.%N)"
active_jobs=0

printf 'Launching %d cell(s) with concurrency=%d\n' "${#INPUTS[@]}" "$JOBS"
for input in "${INPUTS[@]}"; do
	printf '  %s\n' "$(basename "$input")"
	run_one_cell "$input" &
	active_jobs=$((active_jobs + 1))
	if [ "$active_jobs" -ge "$JOBS" ]; then
		wait -n
		active_jobs=$((active_jobs - 1))
	fi
done
wait
end_ts="$(date +%s.%N)"

python3 - "$OUTDIR" "$start_ts" "$end_ts" <<'PY'
import pathlib
import statistics
import sys

outdir = pathlib.Path(sys.argv[1])
start = float(sys.argv[2])
end = float(sys.argv[3])

rows = []
for wall_path in sorted(outdir.glob("*/wall_seconds.txt")):
    cell = wall_path.parent.name
    wall = float(wall_path.read_text().strip())
    rows.append((cell, wall))

if not rows:
    print("No completed runs found.", file=sys.stderr)
    sys.exit(1)

makespan = end - start
sum_wall = sum(w for _, w in rows)
avg_wall = statistics.mean(w for _, w in rows)

print("")
print("=== Multi-Cell Summary ===")
for cell, wall in rows:
    print(f"{cell}\t{wall:.2f}s")
print(f"cells\t{len(rows)}")
print(f"makespan\t{makespan:.2f}s")
print(f"sum_single_cell_wall\t{sum_wall:.2f}s")
print(f"avg_single_cell_wall\t{avg_wall:.2f}s")
print(f"effective_parallelism\t{sum_wall / makespan:.2f}x")
print(f"cells_per_min\t{len(rows) * 60.0 / makespan:.2f}")
PY
