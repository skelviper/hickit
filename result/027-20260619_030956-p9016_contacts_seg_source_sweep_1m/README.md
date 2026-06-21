# 027 P9016 contacts.seg Source Sweep

This experiment tests whether read-level contacts from `contacts.seg.gz` improve full-denominator trans identity when exported as blind-safe raw pairs. Exported pairs preserve readID but write phase labels as `.`; training uses no SNP truth or CHARM/3DG.

- full result root: `/mnt/ssd/zliu/phase3/test_res/027-20260619_030956-p9016_contacts_seg_source_sweep_1m`
- light result root: under `hickit/result` after publish
- headline: `NO_FULL_TRANS_PLUS_0P1`

| config | top1 all | top1 cis | top1 trans | delta trans | same/cross trans | cis Spearman | entropy | pU |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `p9016_seg_all_pcgamma1_ieps0p5_noise0_seed17` | 0.476217 | 0.533829 | 0.353975 | 0.038605 | 0.563828 | 0.797352 | 0.55265522 | 0.398656487 |
| `p9016_approved_baseline_pcgamma1_ieps0p5_noise0_seed17` | 0.465049 | 0.569634 | 0.31537 | 0 | 0.51179 | 0.800612 | 0.607306719 | 0.438079208 |
| `p9016_seg_multiseg_pcgamma1_ieps0p5_noise0_seed17` | 0.448574 | 0.521484 | 0.266229 | -0.049141 | 0.477873 | 0.792498 | 0.596199512 | 0.430067033 |
| `p9016_seg_multichrom_pcgamma1_ieps0p5_noise0_seed17` | 0.254006 | nan | 0.254006 | -0.061364 | 0.502623 | 0.603643 | 0.609811962 | 0.439886332 |
| `p9016_seg_multiseg_multichrom_pcgamma1_ieps0p5_noise0_seed17` | 0.247202 | nan | 0.247202 | -0.068168 | 0.489408 | 0.616479 | 0.6629318 | 0.478204221 |

## Interpretation

- Seg-derived blind-safe raw contacts did not reach full-denominator trans +0.1 in this sweep.
- The read-level source remains useful for future molecule/group constraints, but naive pair export alone is insufficient.
