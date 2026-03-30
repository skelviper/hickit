#!/bin/bash
set -euo pipefail

INPUT="/mnt/ssd/zliu/run_charm/charm_test/charm_test_13/result/impute/pairs/SP0002a01.impute.pairs.gz"
ORIG="/mnt/ssd/zliu/run_charm/hickit/hickit"
OPT="/mnt/ssd/zliu/run_charm/hickit_opt/hickit"
OUTDIR="/mnt/ssd/zliu/run_charm/hickit_opt/bench_results"

mkdir -p "$OUTDIR/orig" "$OUTDIR/opt"

echo "=== Running ORIGINAL (sm_75, linked-list) ==="
/usr/bin/time -f "ORIG wall_clock: %e seconds, max_rss: %M KB" \
    $ORIG -s1 -M \
    -i "$INPUT" -Sr1m -c1 -r10m -c2 \
    -b4m -b1m -O "$OUTDIR/orig/1m.3dg" \
    -b200k -O "$OUTDIR/orig/200k.3dg" \
    -D5 -b50k -O "$OUTDIR/orig/50k.3dg" \
    -D5 -b20k -O "$OUTDIR/orig/20k.3dg" \
    2>&1 | tee "$OUTDIR/orig_log.txt"

echo ""
echo "=== Running OPTIMIZED (sm_89, hash-table) ==="
/usr/bin/time -f "OPT wall_clock: %e seconds, max_rss: %M KB" \
    $OPT -s1 -M \
    -i "$INPUT" -Sr1m -c1 -r10m -c2 \
    -b4m -b1m -O "$OUTDIR/opt/1m.3dg" \
    -b200k -O "$OUTDIR/opt/200k.3dg" \
    -D5 -b50k -O "$OUTDIR/opt/50k.3dg" \
    -D5 -b20k -O "$OUTDIR/opt/20k.3dg" \
    2>&1 | tee "$OUTDIR/opt_log.txt"

echo ""
echo "=== RMSD Comparison ==="
python3 -c "
import sys, math, os

def read_3dg(fn):
    coords = {}
    if not os.path.exists(fn) or os.path.getsize(fn) == 0:
        return coords
    with open(fn) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 5: continue
            key = (parts[0], int(parts[1]))
            coords[key] = (float(parts[2]), float(parts[3]), float(parts[4]))
    return coords

outdir = '$OUTDIR'
for res in ['1m', '200k', '50k', '20k']:
    a = read_3dg(f'{outdir}/orig/{res}.3dg')
    b = read_3dg(f'{outdir}/opt/{res}.3dg')
    if not a or not b:
        print(f'  {res}: SKIP (empty output)')
        continue
    common = set(a.keys()) & set(b.keys())
    if not common:
        print(f'  {res}: No common beads!')
        continue
    ssq = 0.0
    for k in common:
        dx = a[k][0] - b[k][0]
        dy = a[k][1] - b[k][1]
        dz = a[k][2] - b[k][2]
        ssq += dx*dx + dy*dy + dz*dz
    rmsd = math.sqrt(ssq / len(common))
    # Also compute scale (average extent) for context
    all_x = [v[0] for v in a.values()]
    extent = max(all_x) - min(all_x) if all_x else 1.0
    print(f'  {res}: beads={len(common)}, RMSD={rmsd:.6f}, extent~{extent:.2f}, relative_RMSD={rmsd/max(extent,1e-6):.6f}')
"

echo ""
echo "=== Timing Extraction ==="
grep "wall_clock" "$OUTDIR/orig_log.txt" || true
grep "wall_clock" "$OUTDIR/opt_log.txt" || true
