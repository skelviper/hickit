# 018 P9016 Fixed-Posterior Graph Semantics

This is a reference-derived positive-control diagnostic, not a blind baseline. It starts from the best 016 CHARM/3DG-derived 1-pr condition, fixes its binned posterior and coordinates, then tests direct fixed-state-edge M-step variants for one relaxation.

- full result root: `/mnt/ssd/zliu/phase3/test_res/018-20260618_220644-p9016_fixed_posterior_graph_semantics_1m`
- source fixed posterior: `/mnt/ssd/zliu/phase3/test_res/016-20260618_150851-p9016_charm_contacts_condition_sweep_1m/work_outputs/1pr/p9016_charm3dg20k_1pr_pcgamma1_common_msep1p5_lsep1/p9016_charm3dg20k_1pr_pcgamma1_common_msep1p5_lsep1`
- summary: `/mnt/ssd/zliu/phase3/test_res/018-20260618_220644-p9016_fixed_posterior_graph_semantics_1m/summary.tsv`
- supersedes draft: `/mnt/ssd/zliu/phase3/test_res/018-20260618_214516-p9016_fixed_posterior_graph_semantics_1m`

## Boundary

- Training for this diagnostic uses CHARM/3DG-derived fixed posterior and coordinates from the 016 positive control.
- Additional SNP labels and reference reads inside eval scripts are post-relaxation diagnostics only.
- These fixed-state-edge graphs are not native Hickit raw_expected_soft_all replay; the source_016_before_018_relax row is the native-source baseline.
- Because posterior is fixed, posterior-top1 accuracy is expected to be constant across fixed-state configs; geometry metrics are the main signal.

## Configs

| config | graph | native replay | trans dscale | repulsion | min sep | lambda sep |
|---|---:|---:|---:|---:|---:|---:|
| source_016_before_018_relax | source_native_016 | 1 | 1 | 1 | 1.5 | 1 |
| p9016_fixed_state_same_cross_gated | same_cross_gated | 0 | 1 | 1 | 1.5 | 1 |
| p9016_fixed_state_softall | softall | 0 | 1 | 1 | 1.5 | 1 |
| p9016_fixed_state_softall_trans_dscale0p5 | softall | 0 | 0.5 | 1 | 1.5 | 1 |
| p9016_fixed_state_softall_trans_dscale0p75 | softall | 0 | 0.75 | 1 | 1.5 | 1 |
| p9016_fixed_state_top1_only | top1_only | 0 | 1 | 1 | 1.5 | 1 |
| p9016_fixed_state_top2_only | top2_only | 0 | 1 | 1 | 1.5 | 1 |

## Main Results

| config | top1 trans | nearest geometry | truth nearest frac | cis Spearman | trans mean r | trans tail k frac | trans force/k |
|---|---:|---:|---:|---:|---:|---:|---:|
| source_016_before_018_relax | 0.455087 | NA | NA | 0.78049 | NA | NA | NA |
| p9016_fixed_state_same_cross_gated | 0.455087 | 0.456262803 | 0.458880206 | 0.783279 | 2.33212719 | 0.490259882 | 0.806671785 |
| p9016_fixed_state_softall | 0.455087 | 0.456528336 | 0.459183673 | 0.777347 | 2.35265643 | 0.512287046 | 0.832557391 |
| p9016_fixed_state_softall_trans_dscale0p5 | 0.455087 | 0.45679387 | 0.459411274 | 0.775514 | 4.19165357 | 0.83492233 | 1.26293294 |
| p9016_fixed_state_softall_trans_dscale0p75 | 0.455087 | 0.456945604 | 0.459563007 | 0.776217 | 2.95265692 | 0.676027806 | 1.04775608 |
| p9016_fixed_state_top1_only | 0.455087 | 0.456035202 | 0.458690539 | 0.776641 | 2.2740778 | 0.454153475 | 0.760466049 |
| p9016_fixed_state_top2_only | 0.455087 | 0.456186936 | 0.458766406 | 0.779302 | 2.34036013 | 0.502775292 | 0.821609312 |

## Headline

- best nearest-geometry config: `p9016_fixed_state_softall_trans_dscale0p75` with nearest_geometry_state_accuracy=0.456945604
- best cis-Spearman config: `p9016_fixed_state_same_cross_gated` with mean_per_chrom_cis_distance_spearman=0.783279

## Interpretation Rules

- If fixed-state-edge variants improve nearest-geometry accuracy while posterior-top1 accuracy stays fixed, the bottleneck includes M-step graph/force geometry, not only E-step posterior labels.
- If top1/top2/same-cross pruning improves trans geometry but hurts cis Spearman or separation, it is a diagnostic lead, not a baseline candidate.
- If trans dscale multipliers below 1.0 make normalized trans r larger and increase tail-attractive k fraction, that is not evidence for a useful rescue.
- Do not claim model improvement from this run; this is not blind training.
