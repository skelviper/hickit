# 056 P9016 Readchain-Dedup Readgroup Training 1Mb

This experiment tests whole-read/segment-chain deduplication before pair generation. A duplicate chain such as A1-B1-C1 and A2-B2-C2 is collapsed into one representative A-B-C, then adjacent or all-combination contacts are emitted with readgroup metadata.

## Paths

- full result root: `/mnt/ssd/zliu/phase3/test_res/056-20260620_230016-p9016_readchain_dedup_training_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/056-20260620_230016-p9016_readchain_dedup_training_1m`
- summary: `/mnt/ssd/zliu/phase3/test_res/056-20260620_230016-p9016_readchain_dedup_training_1m/summary.tsv`
- source delta summary: `/mnt/ssd/zliu/phase3/test_res/056-20260620_230016-p9016_readchain_dedup_training_1m/source_delta_summary.tsv`
- readgroup matched delta summary: `/mnt/ssd/zliu/phase3/test_res/056-20260620_230016-p9016_readchain_dedup_training_1m/readgroup_matched_delta_summary.tsv`

## Boundary

- `approved_baseline_pairs_only` is the only strict approved `/shared/.../P9016.pairs.gz` baseline row.
- All readchain rows are contacts.seg-derived custom-pairs diagnostics. Phase labels are dropped before training.
- SNP phase and CHARM/3DG are used only by post-training eval through `scripts/p9016_common_eval.sh`.
- copy0/copy1 remain gauge labels.

## Main Results

| config | source | readgroup | n raw | trans raw | top1 all | top1 cis | top1 trans | cis Spearman | entropy | pU | delta trans vs approved |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| approved_baseline_pairs_only | approved_pairs | off | 1703888 | 568434 | 0.467543 | 0.56878 | 0.322656 | 0.800904 | 0.607259214 | 0.438044935 | 0 |
| hickit_like_dup100_pairwise | hickit_like_dup100 | off | 3043937 | 1011325 | 0.426919 | 0.522972 | 0.289449 | 0.775765 | 0.304542035 | 0.219680637 | -0.033207 |
| old_readgroup_adjacent_joint | old_readgroup_adjacent | joint_marginal | 6932882 | 2370651 | 0.426501 | 0.516653 | 0.297478 | 0.731912 | 0.673324943 | 0.485701293 | -0.025178 |
| readchain_adjacent_joint | readchain_adjacent | joint_marginal | 4304745 | 1504017 | 0.419562 | 0.515726 | 0.281934 | 0.769627 | 0.478666574 | 0.345284969 | -0.040722 |
| readchain_adjacent_pairwise | readchain_adjacent | off | 4304745 | 1504017 | 0.414163 | 0.518667 | 0.264598 | 0.786527 | 0.429157823 | 0.309571922 | -0.058058 |
| readchain_allcomb_joint | readchain_allcomb | joint_marginal | 5050272 | 1790830 | 0.38704 | 0.488331 | 0.242073 | 0.688971 | 0.567091167 | 0.409069836 | -0.080583 |
| readchain_allcomb_pairwise | readchain_allcomb | off | 5050272 | 1790830 | 0.407112 | 0.520235 | 0.245213 | 0.732593 | 0.65067488 | 0.469362736 | -0.077443 |

## Matched Readgroup Deltas

| joint config | pairwise config | delta top1 all | delta top1 cis | delta top1 trans | delta cis Spearman |
|---|---|---:|---:|---:|---:|
| readchain_adjacent_joint | readchain_adjacent_pairwise | 0.005399 | -0.002941 | 0.017336 | -0.0169 |
| readchain_allcomb_joint | readchain_allcomb_pairwise | -0.020072 | -0.031904 | -0.00314 | -0.043622 |

## Interpretation

This is a controlled source/interface test. A positive readgroup claim requires the joint row to beat the matched same-source pairwise row without damaging cis accuracy or cis distance Spearman, and ideally to approach or exceed the approved baseline.
