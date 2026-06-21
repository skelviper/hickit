# 027 P9016 contacts.seg Source Sweep

This experiment tests whether read-level contacts from `contacts.seg.gz` improve full-denominator trans identity when exported as blind-safe raw pairs. Exported pairs preserve readID but write phase labels as `.`; training uses no SNP truth or CHARM/3DG.

- full result root: `/mnt/ssd/zliu/phase3/test_res/027-20260619_030424-p9016_contacts_seg_source_sweep_1m`
- light result root: under `hickit/result` after publish
- headline: `NO_FULL_TRANS_PLUS_0P1`

| config | top1 all | top1 cis | top1 trans | delta trans | same/cross trans | cis Spearman | entropy | pU |
|---|---:|---:|---:|---:|---:|---:|---:|---:|

## Interpretation

- Seg-derived blind-safe raw contacts did not reach full-denominator trans +0.1 in this sweep.
- The read-level source remains useful for future molecule/group constraints, but naive pair export alone is insufficient.
