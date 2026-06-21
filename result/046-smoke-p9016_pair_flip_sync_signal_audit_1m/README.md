# 046 P9016 Pair-Flip Sync Signal Audit

This is a post-training diagnostic. It does not retrain Hickit. It asks whether raw posterior/reconstruction geometry contains a stable blind signal for chromosome-pair copy flips that could explain the remaining trans deficit.

- full result root: `/tmp/hk_blind_test_res_046_smoke/046-smoke-p9016_pair_flip_sync_signal_audit_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/046-smoke-p9016_pair_flip_sync_signal_audit_1m`
- headline: `NO_RAW_PAIR_FLIP_SYNC_PLUS_0P05`
- resolution: 1Mb
- training inputs: inherited source runs, all audited as `input_contact_source=raw_pairs`, `uses_phase_labels=0`, and `uses_charm_or_reference=0`.
- eval inputs: P9016 SNP labels from pairs and P9016 1Mb CHARM/3DG are used only after the blind pair-flip choices are fixed.
- copy policy: copy0/copy1 are gauge labels. The cis whole-chromosome SNP gauge is eval-only alignment; raw pair-flip policies must not use SNP for selection.

## Main Blind Rows

| config_name | policy | base_policy | policy_family | top1_accuracy | delta_top1_vs_cis_selected | same_cross_accuracy | split_agree_pair_fraction | oracle_match_pair_fraction | oracle_match_eval_contact_fraction | target_plus_0p1_met |
|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_all_geom_nearest_gap | geom_nearest_gap | blind_raw_all | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_a_geom_nearest_gap | geom_nearest_gap | blind_raw_split_a | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_b_geom_nearest_gap | geom_nearest_gap | blind_raw_split_b | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_agree_geom_nearest_gap | geom_nearest_gap | blind_raw_split_agree_else_zero | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_all_geom_top_to_nearest | geom_top_to_nearest | blind_raw_all | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_a_geom_top_to_nearest | geom_top_to_nearest | blind_raw_split_a | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_b_geom_top_to_nearest | geom_top_to_nearest | blind_raw_split_b | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_agree_geom_top_to_nearest | geom_top_to_nearest | blind_raw_split_agree_else_zero | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_all_margin_geom_gap | margin_geom_gap | blind_raw_all | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_a_margin_geom_gap | margin_geom_gap | blind_raw_split_a | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_b_margin_geom_gap | margin_geom_gap | blind_raw_split_b | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_agree_margin_geom_gap | margin_geom_gap | blind_raw_split_agree_else_zero | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_all_pmax_geom_gap | pmax_geom_gap | blind_raw_all | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_a_pmax_geom_gap | pmax_geom_gap | blind_raw_split_a | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_b_pmax_geom_gap | pmax_geom_gap | blind_raw_split_b | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_agree_pmax_geom_gap | pmax_geom_gap | blind_raw_split_agree_else_zero | 0.353007772 | 0 | 0.552720226 | 1 | 0.494736842 | 0.532605902 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_all_same_cross_distance_sum | same_cross_distance_sum | blind_raw_all | 0.290421387 | -0.0625863853 | 0.4847302 | 0.994736842 | 0.384210526 | 0.398577551 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_b_same_cross_distance_sum | same_cross_distance_sum | blind_raw_split_b | 0.290421387 | -0.0625863853 | 0.4847302 | 0.994736842 | 0.384210526 | 0.398577551 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_agree_same_cross_distance_sum | same_cross_distance_sum | blind_raw_split_agree_else_zero | 0.290421387 | -0.0625863853 | 0.4847302 | 0.994736842 | 0.384210526 | 0.398577551 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_split_a_same_cross_distance_sum | same_cross_distance_sum | blind_raw_split_a | 0.290233857 | -0.0627739153 | 0.484494051 | 0.994736842 | 0.378947368 | 0.395368705 | 0 |

## Interpretation Rule

If a blind raw policy reaches +0.1 full-denominator trans top1, it is a candidate to turn into a controlled training or post-processing experiment. If only the SNP oracle improves, the data still contain pair-level eval-only space but the current raw posterior/geometry does not expose a blind selector. Split-agreement rows test whether the raw pair-flip signal is stable enough to trust without SNP labels.
