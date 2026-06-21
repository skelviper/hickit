# 053 P9016 contacts.seg Source Audit 1Mb

This is a source audit only. No reconstruction was trained here.

## Purpose

The goal is to decide whether upstream `contacts.seg` contains read-level multiway information that is absent from the stripped P9016 `.pairs.gz` file, while explicitly accounting for the fact that the segment file is not deduplicated.

## Why Duplicate Handling Matters

If one read contains repeated or overlapping segment records, expanding all segment pairs can create artificial multiway contacts. This audit therefore reports exact segment duplicates, interval-level duplicates, 1Mb bin duplicates, exact pair duplicates, and collapsed 1Mb bpair duplicates before any training experiment is considered.

## Source Summary

| source | reads | usable >=2 seg | multi-seg frac | multi-chrom frac | multi-seg+multi-chrom frac | exact seg dup frac | 1Mb bpair dup frac | phase labels present |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| sharec | 6085584 | 5903400 | 0.1339 | 0.3291 | 0.0725 | 0.001402 | 0.09341 | 1 |
| archive | 25429920 | 24976813 | 0.3084 | 0.1595 | 0.1093 | 2.375e-05 | 0.2755 | 0 |

## Training Boundary

- This audit reads upstream segment sources but does not train a model.
- Phase-like segment fields, when present, must be dropped for any blind training export. In this audit the sharec source has phase labels, while the archive source did not show phase labels.
- The next training experiment should treat this as an explicit source-boundary change, not as the strict P9016 `.pairs.gz` baseline.

## Proposed Next Step

If these source-level duplicate statistics are acceptable, the next experiment should be a read-level posterior diagnostic: use current model posteriors plus read grouping to test whether read-consistent multi-segment decoding improves trans identity before changing EM/FDG training.
