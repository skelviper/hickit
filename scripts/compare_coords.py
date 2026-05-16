#!/usr/bin/env python3
import sys
import math

def parse_3dg_line(line):
    """Parse a coordinate line from .3dg file"""
    parts = line.strip().split('\t')
    if len(parts) < 4 or parts[0].startswith('#'):
        return None
    try:
        chrom = parts[0]
        pos = int(parts[1])
        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
        return (chrom, pos, x, y, z)
    except:
        return None

def euclidean_distance(p1, p2):
    """Calculate Euclidean distance between two 3D points"""
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

def main():
    if len(sys.argv) != 3:
        print("Usage: compare_coords.py <file1.3dg> <file2.3dg>")
        sys.exit(1)

    file1, file2 = sys.argv[1], sys.argv[2]

    coords1, coords2 = [], []

    with open(file1) as f1, open(file2) as f2:
        for line1, line2 in zip(f1, f2):
            c1 = parse_3dg_line(line1)
            c2 = parse_3dg_line(line2)
            if c1 and c2:
                coords1.append(c1)
                coords2.append(c2)

    if len(coords1) != len(coords2):
        print(f"Warning: Different number of coordinates: {len(coords1)} vs {len(coords2)}")

    # Calculate differences
    diffs = []
    max_diff = 0
    max_diff_idx = 0

    for i, (c1, c2) in enumerate(zip(coords1, coords2)):
        if c1[0] != c2[0] or c1[1] != c2[1]:
            print(f"Warning: Mismatch at line {i}: {c1[:2]} vs {c2[:2]}")
            continue

        dist = euclidean_distance(c1[2:], c2[2:])
        diffs.append(dist)

        if dist > max_diff:
            max_diff = dist
            max_diff_idx = i

    # Statistics
    avg_diff = sum(diffs) / len(diffs)
    sorted_diffs = sorted(diffs)
    median_diff = sorted_diffs[len(sorted_diffs)//2]
    p95_diff = sorted_diffs[int(len(sorted_diffs)*0.95)]
    p99_diff = sorted_diffs[int(len(sorted_diffs)*0.99)]

    print(f"=== Coordinate Difference Statistics ===")
    print(f"Total coordinates compared: {len(diffs)}")
    print(f"Average Euclidean distance: {avg_diff:.6f}")
    print(f"Median distance: {median_diff:.6f}")
    print(f"95th percentile: {p95_diff:.6f}")
    print(f"99th percentile: {p99_diff:.6f}")
    print(f"Maximum distance: {max_diff:.6f}")
    print(f"")
    print(f"Maximum difference at index {max_diff_idx}:")
    c1, c2 = coords1[max_diff_idx], coords2[max_diff_idx]
    print(f"  File1: {c1[0]} {c1[1]} ({c1[2]:.6f}, {c1[3]:.6f}, {c1[4]:.6f})")
    print(f"  File2: {c2[0]} {c2[1]} ({c2[2]:.6f}, {c2[3]:.6f}, {c2[4]:.6f})")
    print(f"")

    # Relative error
    avg_coord_magnitude = sum(math.sqrt(c[2]**2 + c[3]**2 + c[4]**2) for c in coords1) / len(coords1)
    relative_error = avg_diff / avg_coord_magnitude * 100
    print(f"Average coordinate magnitude: {avg_coord_magnitude:.6f}")
    print(f"Relative error: {relative_error:.4f}%")

if __name__ == "__main__":
    main()
