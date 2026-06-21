# 029 P9016 Readgroup Joint-Marginal Training

This experiment tests a blind readgroup-aware E-step prototype. The runner uses `contacts.seg.gz` read grouping only; phase labels and CHARM/3DG are used only by eval after training.

- full result root: `/tmp/hk_blind_test_res_029_smoke/029-smoke-p9016_readgroup_joint_marginal_1m`
- light result root: under `hickit/result` after publish
- headline: `NO_FULL_TRANS_PLUS_0P1`

| config | readgroup | top1 all | top1 cis | top1 trans | delta trans | same/cross trans | cis Spearman | entropy | pU |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `p9016_seg_all_pairwise_pcgamma1_ieps0p5_noise0_seed17` | off | 0.397044 | 0.410441 | 0.293479 | 0.042224 | 0.502022 | 0.191618 | 1.34952307 | 0.973475099 |
| `p9016_seg_all_readgroup_joint_pcgamma1_ieps0p5_noise0_seed17` | joint_marginal | 0.39588 | 0.409905 | 0.287462 | 0.036207 | 0.51603 | 0.191618 | 1.34952307 | 0.973475099 |
| `p9016_approved_baseline_pcgamma1_ieps0p5_noise0_seed17` | off | 0.356188 | 0.429507 | 0.251255 | 0 | 0.493044 | 0.450647 | 0.818740427 | 0.590596378 |
| `p9016_seg_multiseg_multichrom_pairwise_pcgamma1_ieps0p5_noise0_seed17` | off | 0.249386 | nan | 0.249386 | -0.001869 | 0.490041 | 0.126521 | 1.29423416 | 0.933592618 |
| `p9016_seg_multiseg_multichrom_readgroup_joint_pcgamma1_ieps0p5_noise0_seed17` | joint_marginal | 0.248636 | nan | 0.248636 | -0.002619 | 0.489359 | 0.126521 | 1.29423416 | 0.933592618 |

## Interpretation

- The readgroup-aware prototype did not reach full-denominator trans +0.1 in this run.
- If it improves only slightly, read grouping is not sufficient by itself and the remaining problem is likely global trans copy identity/gauge propagation.
