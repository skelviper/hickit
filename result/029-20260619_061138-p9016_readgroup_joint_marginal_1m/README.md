# 029 P9016 Readgroup Joint-Marginal Training

This experiment tests a blind readgroup-aware E-step prototype. The runner uses `contacts.seg.gz` read grouping only; phase labels and CHARM/3DG are used only by eval after training.

- full result root: `/mnt/ssd/zliu/phase3/test_res/029-20260619_061138-p9016_readgroup_joint_marginal_1m`
- light result root: under `hickit/result` after publish
- headline: `NO_FULL_TRANS_PLUS_0P1`

| config | readgroup | top1 all | top1 cis | top1 trans | delta trans | same/cross trans | cis Spearman | entropy | pU |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `p9016_seg_all_pairwise_pcgamma1_ieps0p5_noise0_seed17` | off | 0.476217 | 0.533829 | 0.353975 | 0.038605 | 0.563828 | 0.797352 | 0.55265522 | 0.398656487 |
| `p9016_seg_all_readgroup_joint_pcgamma1_ieps0p5_noise0_seed17` | joint_marginal | 0.465985 | 0.528319 | 0.333725 | 0.018355 | 0.532047 | 0.776265 | 0.531492949 | 0.383391142 |
| `p9016_approved_baseline_pcgamma1_ieps0p5_noise0_seed17` | off | 0.465049 | 0.569634 | 0.31537 | 0 | 0.51179 | 0.800612 | 0.607306719 | 0.438079208 |
| `p9016_seg_multiseg_multichrom_readgroup_joint_pcgamma1_ieps0p5_noise0_seed17` | joint_marginal | 0.250306 | nan | 0.250306 | -0.065064 | 0.503516 | 0.279528 | 0.233481124 | 0.16842103 |
| `p9016_seg_multiseg_multichrom_pairwise_pcgamma1_ieps0p5_noise0_seed17` | off | 0.247202 | nan | 0.247202 | -0.068168 | 0.489408 | 0.616479 | 0.6629318 | 0.478204221 |

## Interpretation

- The readgroup-aware prototype did not reach full-denominator trans +0.1 in this run.
- If it improves only slightly, read grouping is not sufficient by itself and the remaining problem is likely global trans copy identity/gauge propagation.
