# 027 P9016 contacts.seg Source Sweep

This experiment tests whether read-level contacts from `contacts.seg.gz` improve full-denominator trans identity when exported as blind-safe raw pairs. Exported pairs preserve readID but write phase labels as `.`; training uses no SNP truth or CHARM/3DG.

- full result root: `/mnt/ssd/zliu/phase3/test_res/027-20260619_034240-p9016_contacts_seg_source_sweep_1m_fast20`
- light result root: under `hickit/result` after publish
- headline: `NO_FULL_TRANS_PLUS_0P1`

| config | top1 all | top1 cis | top1 trans | delta trans | same/cross trans | cis Spearman | entropy | pU |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `p9016_approved_baseline_pcgamma1_ieps0p5_noise0_seed17` | 0.456901 | 0.557409 | 0.313057 | 0 | 0.52948 | 0.805746 | 0.636614501 | 0.45922029 |
| `p9016_seg_all_pcgamma1_ieps0p5_noise0_seed17` | 0.450147 | 0.52912 | 0.282582 | -0.030475 | 0.49392 | 0.789421 | 0.583203077 | 0.420692086 |
| `p9016_seg_multiseg_pcgamma1_ieps0p5_noise0_seed17` | 0.447867 | 0.514744 | 0.280611 | -0.032446 | 0.529227 | 0.787611 | 0.630223632 | 0.454610229 |
| `p9016_seg_multichrom_pcgamma1_ieps0p5_noise0_seed17` | 0.252508 | nan | 0.252508 | -0.060549 | 0.503621 | 0.684903 | 0.743163824 | 0.536079407 |
| `p9016_seg_multiseg_multichrom_pcgamma1_ieps0p5_noise0_seed17` | 0.247202 | nan | 0.247202 | -0.065855 | 0.486649 | 0.653222 | 0.709600985 | 0.511868894 |

## Interpretation

- Seg-derived blind-safe raw contacts did not reach full-denominator trans +0.1 in this sweep.
- The read-level source remains useful for future molecule/group constraints, but naive pair export alone is insufficient.
