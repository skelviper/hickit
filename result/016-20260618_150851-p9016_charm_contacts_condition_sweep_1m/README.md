# 016-20260618_150851-p9016_charm_contacts_condition_sweep_1m

Controlled 016 sweep: reuse CHARM20k-derived 1PR/2PR unphased synthetic contacts, then test prior P9016 model-condition knobs that are already implemented and blind-safe with respect to SNP labels. Training remains reference-derived positive control because contacts are generated from CHARM/3DG.

## Paths

- full result root: `/mnt/ssd/zliu/phase3/test_res/016-20260618_150851-p9016_charm_contacts_condition_sweep_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/016-20260618_150851-p9016_charm_contacts_condition_sweep_1m`
- summary.tsv: `/mnt/ssd/zliu/phase3/test_res/016-20260618_150851-p9016_charm_contacts_condition_sweep_1m/summary.tsv`
- condition_matrix.tsv: `/mnt/ssd/zliu/phase3/test_res/016-20260618_150851-p9016_charm_contacts_condition_sweep_1m/condition_matrix.tsv`
- contact matrix QC: `/mnt/ssd/zliu/phase3/test_res/016-20260618_150851-p9016_charm_contacts_condition_sweep_1m/contact_matrix_qc`
- commands.log: `/mnt/ssd/zliu/phase3/test_res/016-20260618_150851-p9016_charm_contacts_condition_sweep_1m/commands.log`

## Boundary

This is not a blind baseline: training contacts are generated from CHARM/3DG. Synthetic pairs remain unphased, `uses_phase_labels=0`, and SNP labels plus P9016.1m.3DG are used only by post-training evaluation.

## Git And Build

- git commit at launch: `1037ee5f79d54e943eef5c37e5aeb2d841698cb8`
- git dirty count at launch: `13`
- WARNING: dirty tree at launch; inspect `logs/git_status.txt` and `logs/git_diff.patch` in full result root.
- run binary sha256: `563303f22b0cbbfa2da75c2af460468b4ce2f0cc3c26136c0857512c1ac5a6f0`
- backend: `gpu`

## Best By Top1 All

| radius | condition | top1 all | top1 cis | top1 trans | same/cross trans | cis Spearman | entropy | pU | min sep | mean sep |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1pr | `pcgamma1_common_msep1p5_lsep1` | 0.616134 | 0.638402 | 0.455087 | 0.574577 | 0.78049 | 0.479795188 | 0.346099079 | 1.14863467 | 5.74935198 |
| 1pr | `pcgamma1_common_ieps1_noise0p05` | 0.616042 | 0.63827 | 0.455277 | 0.56259 | 0.77942 | 0.478756875 | 0.345350087 | 0.804342091 | 5.80347681 |
| 1pr | `expected_common_msep1p5_lsep0p5` | 0.615788 | 0.638008 | 0.455087 | 0.575298 | 0.779509 | 0.479984224 | 0.346235424 | 1.0327388 | 5.75481367 |
| 1pr | `raw_common_ieps1_noise0p05` | 0.602094 | 0.629611 | 0.40308 | 0.531295 | 0.83663 | 0.598924518 | 0.432032704 | 0.0401670635 | 5.63051891 |
| 2pr | `pcgamma1_common_ieps1_noise0p05` | 0.581595 | 0.612049 | 0.380684 | 0.497841 | 0.778746 | 0.538737595 | 0.388617039 | 0.970586777 | 5.33867073 |
| 2pr | `expected_common_msep1p5_lsep0p5` | 0.5814 | 0.611824 | 0.380684 | 0.500978 | 0.778696 | 0.537456691 | 0.387693018 | 1.23422301 | 5.34118509 |
| 2pr | `pcgamma1_common_msep1p5_lsep1` | 0.580689 | 0.610817 | 0.381932 | 0.518554 | 0.776927 | 0.536366284 | 0.386906475 | 1.2473495 | 5.31039095 |
| 2pr | `raw_common_ieps1_noise0p05` | 0.556849 | 0.578884 | 0.411483 | 0.572696 | 0.795088 | 0.703822672 | 0.507700741 | 0.0246444009 | 5.01783752 |

## Best By Top1 Trans

| radius | condition | top1 trans | top1 all | top1 cis | same/cross trans | cis Spearman | entropy | pU | copytrack cos<0 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1pr | `pcgamma1_common_ieps1_noise0p05` | 0.455277 | 0.616042 | 0.63827 | 0.56259 | 0.77942 | 0.478756875 | 0.345350087 | 0.000766577233 |
| 1pr | `expected_common_msep1p5_lsep0p5` | 0.455087 | 0.615788 | 0.638008 | 0.575298 | 0.779509 | 0.479984224 | 0.346235424 | 0.000766577233 |
| 1pr | `pcgamma1_common_msep1p5_lsep1` | 0.455087 | 0.616134 | 0.638402 | 0.574577 | 0.78049 | 0.479795188 | 0.346099079 | 0.000766577233 |
| 1pr | `pcgamma1_msep1p5_lsep1_ieps0p5_noise0` | 0.436158 | 0.555808 | 0.572352 | 0.564676 | 0.581216 | 0.569496095 | 0.4108046 | 0.00766577233 |
| 2pr | `raw_common_ieps1_noise0p05` | 0.411483 | 0.556849 | 0.578884 | 0.572696 | 0.795088 | 0.703822672 | 0.507700741 | 0.00421617478 |
| 2pr | `pcgamma1_msep1p5_lsep1_ieps0p5_noise0` | 0.403623 | 0.488477 | 0.50134 | 0.527695 | 0.499117 | 0.693723202 | 0.500415504 | 0.0256803373 |
| 1pr | `raw_common_ieps1_noise0p05` | 0.40308 | 0.602094 | 0.629611 | 0.531295 | 0.83663 | 0.598924518 | 0.432032704 | 0.00191644308 |
| 2pr | `raw_msep1p5_lsep1_ieps0p5_noise0` | 0.392019 | 0.363949 | 0.359694 | 0.529955 | 0.215679 | 1.11598802 | 0.805015206 | 0.119969337 |

## Best By Cis Spearman

| radius | condition | cis Spearman | top1 all | top1 cis | top1 trans | entropy | pU |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1pr | `raw_common_ieps1_noise0p05` | 0.83663 | 0.602094 | 0.629611 | 0.40308 | 0.598924518 | 0.432032704 |
| 2pr | `raw_common_ieps1_noise0p05` | 0.795088 | 0.556849 | 0.578884 | 0.411483 | 0.703822672 | 0.507700741 |
| 1pr | `pcgamma1_common_msep1p5_lsep1` | 0.78049 | 0.616134 | 0.638402 | 0.455087 | 0.479795188 | 0.346099079 |
| 1pr | `expected_common_msep1p5_lsep0p5` | 0.779509 | 0.615788 | 0.638008 | 0.455087 | 0.479984224 | 0.346235424 |
| 1pr | `pcgamma1_common_ieps1_noise0p05` | 0.77942 | 0.616042 | 0.63827 | 0.455277 | 0.478756875 | 0.345350087 |
| 2pr | `pcgamma1_common_ieps1_noise0p05` | 0.778746 | 0.581595 | 0.612049 | 0.380684 | 0.538737595 | 0.388617039 |
| 2pr | `expected_common_msep1p5_lsep0p5` | 0.778696 | 0.5814 | 0.611824 | 0.380684 | 0.537456691 | 0.387693018 |
| 2pr | `pcgamma1_common_msep1p5_lsep1` | 0.776927 | 0.580689 | 0.610817 | 0.381932 | 0.536366284 | 0.386906475 |

## Interpretation Rules

- A condition is useful only if it improves trans without sacrificing cis Spearman, coverage, entropy/pU, or copytrack continuity.
- Any positive signal here remains a reference-derived positive control, not a blind production baseline.
- If no condition beats 015 robustly, the bottleneck is not fixed by previously implemented dscale/init/sep knobs on clean contacts.
