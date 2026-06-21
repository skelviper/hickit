# 017 P9016 Trans Gauge and Force Diagnostics

This is a diagnostic-only run. It does not retrain the model. It reads the existing 016 best positive-control posterior/coordinates and uses phase labels only for eval-only trans diagnostics.

## Paths

- full result root: `/mnt/ssd/zliu/phase3/test_res/017-20260618_210650-p9016_trans_gauge_force_diagnostics`
- source config: `p9016_charm3dg20k_1pr_pcgamma1_common_msep1p5_lsep1`
- diagnostics: `/mnt/ssd/zliu/phase3/test_res/017-20260618_210650-p9016_trans_gauge_force_diagnostics/diagnostics/p9016_charm3dg20k_1pr_pcgamma1_common_msep1p5_lsep1`

## Main Trans Gauge/Geometry Result

`chr_pair_oracle_*` is an eval-only upper-bound diagnostic: within each trans chromosome pair it chooses the relative copy flip using SNP-truth contacts, so it is not a blind selection metric. `nearest_geometry_state_accuracy` uses the nearest copy-pair state in the reconstruction under the whole-chrom SNP gauge; it is not CHARM/3DG nearest-state accuracy and not posterior-top1 accuracy.

| metric | value |
|---|---:|
| n_contacts | 26362 |
| whole_chrom_top1_accuracy | 0.455086867 |
| chr_pair_oracle_top1_accuracy | 0.552310143 |
| chr_pair_oracle_delta_vs_whole | 0.0972232759 |
| nearest_geometry_state_accuracy | 0.455693802 |
| same_cross_accuracy | 0.574577043 |
| pmax90_accuracy | 0.495463623 |
| pmax90_recall | 0.329375616 |
| truth_state_distance_rank_mean | 2.05883469 |
| truth_state_is_nearest_fraction | 0.458273272 |
| posterior_top1_distance_rank_mean | 1.00367954 |
| posterior_top1_is_nearest_fraction | 0.998255064 |
| posterior_top1_distance_gap_mean | 2.90373066 |

## Force Regime Summary

These force-regime rows are reconstructed from final bpair posterior state edges. `sum_abs_force` and `force_per_k` use posterior-edge `|fmag|` for regime comparison; they are not expected to equal the runner's exact `contact_force_l1`, which is reported separately in the source `p9016_full.force_class_diag.tsv` with direction-L1 components on the exact training graph.

| class | sum_k | force_per_k | mean_r | frac_k_zero_force | frac_k_tail_attractive | cancellation_ratio |
|---|---:|---:|---:|---:|---:|---:|
| cis | 44303.0493 | 0.268390321 | 1.41071499 | 0.497475298 | 0.0869401206 | 0.580785932 |
| trans | 4306.00921 | 0.845699361 | 2.39207964 | 0.216065384 | 0.51808262 | 0.37064611 |

## Interpretation Rules

- If chr_pair_oracle_top1_accuracy is much higher than whole_chrom_top1_accuracy, trans is partly a relative chromosome-pair gauge problem.
- If truth_state_is_nearest_fraction remains low, the reconstruction geometry itself does not put SNP-supported trans states nearest.
- Force-regime rows are reconstructed from final bpair posterior state edges, not the unsaved raw softall graph. Use them as a target-distance and cancellation diagnostic, while the source p9016_full.force_class_diag.tsv remains the exact runner graph summary.
- If trans has high force_per_k but high cancellation_ratio or high tail-attractive fraction, the problem is likely long-range force conflict or dscale target semantics rather than absence of trans force.
- This run is not a blind model improvement claim; it is a mechanism diagnostic for choosing the next controlled experiment.
