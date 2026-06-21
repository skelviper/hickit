# 018 P9016 Fixed-Posterior Graph Semantics

This is a reference-derived positive-control diagnostic, not a blind baseline. It starts from the best 016 CHARM/3DG-derived 1-pr condition, fixes its posterior, and changes only the M-step graph/target-distance semantics for one relaxation.

- full result root: `/mnt/ssd/zliu/phase3/test_res/018-20260618_214516-p9016_fixed_posterior_graph_semantics_1m`
- source fixed posterior: `/mnt/ssd/zliu/phase3/test_res/016-20260618_150851-p9016_charm_contacts_condition_sweep_1m/work_outputs/1pr/p9016_charm3dg20k_1pr_pcgamma1_common_msep1p5_lsep1/p9016_charm3dg20k_1pr_pcgamma1_common_msep1p5_lsep1`
- summary: `/mnt/ssd/zliu/phase3/test_res/018-20260618_214516-p9016_fixed_posterior_graph_semantics_1m/summary.tsv`

## Boundary

- Training for this diagnostic uses fixed posterior and coordinates from the 016 CHARM/3DG-derived positive control.
- SNP labels and CHARM reference are still used only in eval/diagnostics after relaxation.
- Because posterior is fixed, posterior-top1 accuracy is expected to be constant across configs; geometry metrics are the main signal.

## Configs

| config | graph | trans dscale | min sep | lambda sep |
|---|---:|---:|---:|---:|
| p9016_fixed_pc_same_cross_gated | same_cross_gated | 1 | 1.5 | 1 |
| p9016_fixed_pc_softall_replay | softall | 1 | 1.5 | 1 |
| p9016_fixed_pc_softall_trans_dscale0p5 | softall | 0.5 | 1.5 | 1 |
| p9016_fixed_pc_softall_trans_dscale0p75 | softall | 0.75 | 1.5 | 1 |
| p9016_fixed_pc_top1_only | top1_only | 1 | 1.5 | 1 |
| p9016_fixed_pc_top2_only | top2_only | 1 | 1.5 | 1 |

## Main Results

| config | top1 trans | nearest geometry | truth nearest frac | cis Spearman | trans mean r | trans tail k frac | trans force/k |
|---|---:|---:|---:|---:|---:|---:|---:|
| p9016_fixed_pc_same_cross_gated | 0.455087 | 0.455845535 | 0.458349139 | 0.784787 | 2.17497217 | 0.422221979 | 0.71993813 |
| p9016_fixed_pc_softall_replay | 0.455087 | 0.457324937 | 0.460056141 | 0.779202 | 2.19755331 | 0.443491461 | 0.747000329 |
| p9016_fixed_pc_softall_trans_dscale0p5 | 0.455087 | 0.457476671 | 0.460169942 | 0.776702 | 3.8095651 | 0.774748731 | 1.18500049 |
| p9016_fixed_pc_softall_trans_dscale0p75 | 0.455087 | 0.457552538 | 0.460283742 | 0.777715 | 2.72260128 | 0.602524111 | 0.958622597 |
| p9016_fixed_pc_top1_only | 0.455087 | 0.456300736 | 0.458994007 | 0.779803 | 2.12642362 | 0.395542598 | 0.681807555 |
| p9016_fixed_pc_top2_only | 0.455087 | 0.457211137 | 0.459828541 | 0.780691 | 2.18519839 | 0.433880832 | 0.735330188 |

## Headline

- best nearest-geometry config: `p9016_fixed_pc_softall_trans_dscale0p75` with nearest_geometry_state_accuracy=0.457552538
- best cis-Spearman config: `p9016_fixed_pc_same_cross_gated` with mean_per_chrom_cis_distance_spearman=0.784787

## Interpretation Rules

- If nearest-geometry accuracy improves while posterior-top1 accuracy stays fixed, the bottleneck includes M-step graph/force geometry, not only E-step posterior labels.
- If top1/top2/same-cross pruning improves trans geometry but hurts cis Spearman or separation, it is a diagnostic lead, not a baseline candidate.
- In this run, trans dscale multipliers below 1.0 made normalized trans r larger and increased tail-attractive k fraction; that is not evidence for a useful rescue.
- Do not claim model improvement from this run; this is not blind training.
