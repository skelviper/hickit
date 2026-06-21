# 054 P9016 contacts.seg Readgroup/Dedup Audit 1Mb

This run exports blind-safe pairs from upstream contacts.seg sources and audits hickit-like pair-level duplicate semantics. No reconstruction is trained in 054.

## Paths

- full result root: `/mnt/ssd/zliu/phase3/test_res/054-20260620_163922-p9016_contacts_seg_readgroup_dedup_audit_1m`
- summary: `/mnt/ssd/zliu/phase3/test_res/054-20260620_163922-p9016_contacts_seg_readgroup_dedup_audit_1m/summary.tsv`

## Boundary

- Training is not run here.
- Segment phase labels are dropped from all exported pairs.
- hickit_like exports are reference/control inputs that mimic hickit adjacent-segment, boundary-coordinate, dup-dist pair handling.
- readgroup exports preserve numeric read_group_id/read_seg columns for later readgroup-aware diagnostics.

## Export Summary

| source | export | pairs written | dedup rate | readgroup cols | hickit-like | path |
|---|---|---:|---:|---:|---:|---|
| sharec | hickit_like_dup0 | 6658141 | 0 | 0 | 0 | `/mnt/ssd/zliu/phase3/test_res/054-20260620_163922-p9016_contacts_seg_readgroup_dedup_audit_1m/exports/sharec/hickit_like_dup0.pairs.gz` |
| sharec | hickit_like_dup100 | 3043937 | 0.542824791 | 0 | 1 | `/mnt/ssd/zliu/phase3/test_res/054-20260620_163922-p9016_contacts_seg_readgroup_dedup_audit_1m/exports/sharec/hickit_like_dup100.pairs.gz` |
| sharec | readgroup_adjacent_maxseg10_annotated | 6932882 | 0.553691812 | 1 | 0 | `/mnt/ssd/zliu/phase3/test_res/054-20260620_163922-p9016_contacts_seg_readgroup_dedup_audit_1m/exports/sharec/readgroup_adjacent_maxseg10_annotated.pairs.gz` |
| sharec | readgroup_allcomb_maxseg10_annotated | 7954066 | 0.550405918 | 1 | 0 | `/mnt/ssd/zliu/phase3/test_res/054-20260620_163922-p9016_contacts_seg_readgroup_dedup_audit_1m/exports/sharec/readgroup_allcomb_maxseg10_annotated.pairs.gz` |

## Next Step

055 should train/evaluate a minimal matrix using the approved baseline, hickit_like_dup100, readgroup_adjacent pairwise/joint, and readgroup_allcomb pairwise/joint.
