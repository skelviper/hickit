# 048 P9016 Callable-Anchor M-step 1Mb

This controlled blind-training experiment tests whether a sparse warmup callable-anchor M-step from the current blind posterior can reduce trans copy-state failures. The callable subset is chosen only from raw-pair-derived posterior confidence, never from SNP labels or CHARM/3DG reference coordinates.

- full result root: `/mnt/ssd/zliu/phase3/test_res/048-20260620_092323-p9016_callable_anchor_mstep_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/048-20260620_092323-p9016_callable_anchor_mstep_1m`
- headline: `NO_PLUS_0P1`
- training inputs: approved P9016 raw pairs only
- eval inputs: P9016 SNP labels and CHARM/3DG are used only after training
- shared baseline knobs: posterior-count d_scale gamma=1, trans_dscale_multiplier=0.5, lambda_copytrack=0.03, lambda_global_copytrack=0.003, no sep force
- success criterion: +0.1 full-denominator trans top1 over current best blind-safe standard trans, recorded as 0.358585 from prior result audit

## Main Results

| config_name | experiment_family | seed_family | trans_callable_anchor_top_frac | trans_callable_anchor_mix_weight | trans_chr_pair_mstep_lambda | trans_chr_pair_mstep_aggregate_mode | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | delta_trans_vs_same_seed_off | delta_trans_vs_current_best047 | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | copytrack_frac_cos_lt_0 | target_plus_0p1_vs_current_best047_met |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_callable_negent_t0p10_m0p05_w10_seed71_best035 | callable_anchor | seed71_noise0p05 | 0.100000001 | 0.0500000007 | 0 | state4 | 0.483733 | 0.573128 | 0.355793 | 0.000132 | -0.002792 | 0.561437 | 0.747643 | 0.692738771 | 0.499705404 | 0.356060803 | 2.58023596 | 0.00344431688 | 0 |
| p9016_callable_negent_t0p02_m0p05_w10_seed71_best035 | callable_anchor | seed71_noise0p05 | 0.0199999996 | 0.0500000007 | 0 | state4 | 0.483424 | 0.572643 | 0.355737 | 7.6e-05 | -0.002848 | 0.561159 | 0.747986 | 0.692549467 | 0.49956885 | 0.370461911 | 2.57392073 | 0.00344431688 | 0 |
| p9016_callable_negent_t0p10_m0p10_w10_seed71_best035 | callable_anchor | seed71_noise0p05 | 0.100000001 | 0.100000001 | 0 | state4 | 0.483444 | 0.572725 | 0.355668 | 7e-06 | -0.002917 | 0.561152 | 0.747764 | 0.69248271 | 0.499520659 | 0.366522461 | 2.58347106 | 0.00344431688 | 0 |
| p9016_callable_anchor_off_seed71_noise0p05_best035 | control | seed71_noise0p05 | 0 | 0 | 0 | state4 | 0.483299 | 0.572482 | 0.355661 | 0 | -0.002924 | 0.561208 | 0.748009 | 0.692783058 | 0.499737322 | 0.370910287 | 2.57510543 | 0.00344431688 | 0 |
| p9016_callable_negent_t0p01_m0p10_w10_seed71_best035 | callable_anchor | seed71_noise0p05 | 0.00999999978 | 0.100000001 | 0 | state4 | 0.483684 | 0.573147 | 0.355647 | -1.4e-05 | -0.002938 | 0.561312 | 0.748077 | 0.692820668 | 0.499764442 | 0.370443165 | 2.58011484 | 0.00344431688 | 0 |
| p9016_callable_negent_t0p05_m0p05_w10_seed71_best035 | callable_anchor | seed71_noise0p05 | 0.0500000007 | 0.0500000007 | 0 | state4 | 0.483227 | 0.572371 | 0.355647 | -1.4e-05 | -0.002938 | 0.561201 | 0.747949 | 0.692476809 | 0.499516398 | 0.358635813 | 2.58544278 | 0.00344431688 | 0 |
| p9016_callable_negent_t0p01_m0p05_w10_seed71_best035 | callable_anchor | seed71_noise0p05 | 0.00999999978 | 0.0500000007 | 0 | state4 | 0.483336 | 0.572579 | 0.355612 | -4.9e-05 | -0.002973 | 0.561194 | 0.747996 | 0.692750335 | 0.499713749 | 0.370459735 | 2.57420516 | 0.00344431688 | 0 |
| p9016_callable_negent_t0p02_m0p10_w10_seed71_best035 | callable_anchor | seed71_noise0p05 | 0.0199999996 | 0.100000001 | 0 | state4 | 0.483224 | 0.572429 | 0.355557 | -0.000104 | -0.003028 | 0.561242 | 0.748036 | 0.692598164 | 0.499603957 | 0.369866371 | 2.57249475 | 0.00344431688 | 0 |
| p9016_callable_negent_t0p05_m0p10_w10_seed71_best035 | callable_anchor | seed71_noise0p05 | 0.0500000007 | 0.100000001 | 0 | state4 | 0.483173 | 0.572468 | 0.355376 | -0.000285 | -0.003209 | 0.561013 | 0.747819 | 0.692803442 | 0.499752045 | 0.367561787 | 2.58398843 | 0.00344431688 | 0 |
| p9016_samecross_m0p25_p1_w10_seed71_best035 | samecross | seed71_noise0p05 | 0 | 0 | 0.25 | same_cross | 0.484339 | 0.574477 | 0.355335 | -0.000326 | -0.00325 | 0.562909 | 0.763075 | 0.685310841 | 0.494347274 | 1.15381444 | 2.73850608 | 0.00344431688 | 0 |
| p9016_callable_negent_t0p05_m0p05_w10_seed17_best035 | callable_anchor | seed17_noise0 | 0.0500000007 | 0.0500000007 | 0 | state4 | 0.482373 | 0.572754 | 0.353022 | 0.000188 | -0.005563 | 0.552699 | 0.769506 | 0.696921229 | 0.502722383 | 0.464811951 | 2.83231521 | 0.0015308075 | 0 |
| p9016_callable_negent_t0p05_m0p10_w10_seed17_best035 | callable_anchor | seed17_noise0 | 0.0500000007 | 0.100000001 | 0 | state4 | 0.481753 | 0.571706 | 0.353015 | 0.000181 | -0.00557 | 0.552616 | 0.769406 | 0.696930826 | 0.502729356 | 0.458690077 | 2.82583332 | 0.0015308075 | 0 |
| p9016_callable_negent_t0p01_m0p05_w10_seed17_best035 | callable_anchor | seed17_noise0 | 0.00999999978 | 0.0500000007 | 0 | state4 | 0.482265 | 0.572579 | 0.353008 | 0.000174 | -0.005577 | 0.552637 | 0.769446 | 0.696853876 | 0.502673805 | 0.462803274 | 2.83677626 | 0.0015308075 | 0 |
| p9016_callable_negent_t0p01_m0p10_w10_seed17_best035 | callable_anchor | seed17_noise0 | 0.00999999978 | 0.100000001 | 0 | state4 | 0.481979 | 0.572094 | 0.353008 | 0.000174 | -0.005577 | 0.552637 | 0.769483 | 0.696803987 | 0.502637804 | 0.460427642 | 2.83314967 | 0.0015308075 | 0 |
| p9016_callable_negent_t0p02_m0p05_w10_seed17_best035 | callable_anchor | seed17_noise0 | 0.0199999996 | 0.0500000007 | 0 | state4 | 0.48205 | 0.572235 | 0.35298 | 0.000146 | -0.005605 | 0.552623 | 0.769464 | 0.696848571 | 0.50266999 | 0.460000634 | 2.83231759 | 0.0015308075 | 0 |
| p9016_callable_negent_t0p10_m0p10_w10_seed17_best035 | callable_anchor | seed17_noise0 | 0.100000001 | 0.100000001 | 0 | state4 | 0.482116 | 0.57239 | 0.352917 | 8.3e-05 | -0.005668 | 0.552713 | 0.769525 | 0.69686383 | 0.502681017 | 0.460047871 | 2.82815051 | 0.0015308075 | 0 |
| p9016_callable_anchor_off_seed17_noise0_best035 | control | seed17_noise0 | 0 | 0 | 0 | state4 | 0.481959 | 0.572181 | 0.352834 | 0 | -0.005751 | 0.552685 | 0.769483 | 0.69686389 | 0.502681077 | 0.462113857 | 2.83176327 | 0.0015308075 | 0 |
| p9016_callable_negent_t0p02_m0p10_w10_seed17_best035 | callable_anchor | seed17_noise0 | 0.0199999996 | 0.100000001 | 0 | state4 | 0.481907 | 0.572138 | 0.352772 | -6.2e-05 | -0.005813 | 0.552484 | 0.76942 | 0.696946263 | 0.502740502 | 0.463747054 | 2.83394909 | 0.0015308075 | 0 |
| p9016_callable_negent_t0p10_m0p05_w10_seed17_best035 | callable_anchor | seed17_noise0 | 0.100000001 | 0.0500000007 | 0 | state4 | 0.48183 | 0.572138 | 0.352584 | -0.00025 | -0.006001 | 0.552672 | 0.769291 | 0.696723163 | 0.50257951 | 0.466259211 | 2.8360095 | 0.0015308075 | 0 |
| p9016_samecross_m0p50_p1_w10_seed71_best035 | samecross | seed71_noise0p05 | 0 | 0 | 0.5 | same_cross | 0.481116 | 0.573841 | 0.34841 | -0.007251 | -0.010175 | 0.547275 | 0.782815 | 0.674104869 | 0.486263871 | 0.910721719 | 2.99513173 | 0.00229621125 | 0 |
| p9016_samecross_m0p50_p1_w10_seed17_best035 | samecross | seed17_noise0 | 0 | 0 | 0.5 | same_cross | 0.476717 | 0.567848 | 0.346291 | -0.006543 | -0.012294 | 0.540052 | 0.772799 | 0.674608469 | 0.486627162 | 0.847754836 | 2.77142906 | 0.00229621125 | 0 |
| p9016_samecross_m0p25_p1_w10_seed17_best035 | samecross | seed17_noise0 | 0 | 0 | 0.25 | same_cross | 0.477422 | 0.57007 | 0.344826 | -0.008008 | -0.013759 | 0.538906 | 0.766235 | 0.694410324 | 0.500911176 | 0.628052831 | 2.69114518 | 0.00267891313 | 0 |

## Interpretation Boundary

If the headline remains NO_PLUS_0P1, this experiment did not achieve the requested trans gain. A small positive shift should be treated only as a mechanistic lead unless it is consistent across seeds and does not damage cis/Spearman. If the callable-anchor rows improve trans while samecross rows do not, the bottleneck is more likely local confident-state propagation than chromosome-pair same/cross bias. If neither improves, the current raw-only posterior may not contain enough blind-selectable trans identity signal for this M-step family.
