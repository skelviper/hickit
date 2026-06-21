# 040 P9016 Normalized Direction Copytrack 1Mb

This controlled blind-training experiment tests whether adjacent normalized homolog-vector direction continuity improves trans copy identity beyond the best 035 raw-pairs baseline.

- full result root: `/mnt/ssd/zliu/phase3/test_res/040-20260620_010533-p9016_normdir_copytrack_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/040-20260620_010533-p9016_normdir_copytrack_1m`
- headline: `NO_PLUS_0P1`
- training input: raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper
- fixed baseline knobs: scaffold split init, posterior_count gamma1 dscale, trans_dscale_multiplier=0.5, constant rho, temperature 1, no trans chr-pair prior
- changed knob: lambda_normdir_copytrack for normalized(copy1-copy0) adjacent direction continuity; optional msep/lsep guardrail in one row

## Main Results

| config_name | lambda_normdir_copytrack | lambda_copytrack | lambda_global_copytrack | min_sep_unit | lambda_sep | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | delta_trans_vs_best035_replay | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | copytrack_frac_cos_lt_0 | copytrack_frac_projection_sign_switch | normdir_copytrack_force_l1_over_contact_force_l1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_td0p5_normdir0p03_lct0p03_gct0p003 | 0.0299999993 | 0.0299999993 | 0.00300000003 | 0 | 0 | 0.482087 | 0.572235 | 0.35307 | 0.000153 | 0.547532 | 0.760968 | 0.702418208 | 0.506687641 | 0.372794926 | 2.69699645 | 0.0015308075 | 0.00191350938 | 0.000151340294 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035_replay | 0 | 0.0299999993 | 0.00300000003 | 0 | 0 | 0.482062 | 0.572298 | 0.352917 | 0 | 0.552554 | 0.769497 | 0.696893394 | 0.502702355 | 0.462486386 | 2.83474708 | 0.0015308075 | 0.0045924225 | 0 |
| p9016_pcgamma1_td0p5_normdir0p01_lct0p03_gct0p003 | 0.00999999978 | 0.0299999993 | 0.00300000003 | 0 | 0 | 0.481939 | 0.572157 | 0.35282 | -9.7e-05 | 0.547469 | 0.761469 | 0.704168677 | 0.507950306 | 0.388070494 | 2.65233684 | 0.00229621125 | 0.003061615 | 5.1094019e-05 |
| p9016_pcgamma1_td0p5_normdir0p003_lct0p03_gct0p003 | 0.00300000003 | 0.0299999993 | 0.00300000003 | 0 | 0 | 0.478882 | 0.568857 | 0.350111 | -0.002806 | 0.549303 | 0.761056 | 0.703653872 | 0.507578969 | 0.359274387 | 2.62727118 | 0.00229621125 | 0.00420972063 | 1.53949086e-05 |
| p9016_pcgamma1_td0p5_normdir0p1_lct0p03_gct0p003 | 0.100000001 | 0.0299999993 | 0.00300000003 | 0 | 0 | 0.476402 | 0.571051 | 0.340943 | -0.011974 | 0.542198 | 0.749153 | 0.700139105 | 0.505043626 | 0.541474283 | 2.65043998 | 0.0015308075 | 0.00191350938 | 0.000491717077 |
| p9016_pcgamma1_td0p5_normdir0p03_only | 0.0299999993 | 0 | 0 | 0 | 0 | 0.475182 | 0.570119 | 0.339311 | -0.013606 | 0.534391 | 0.759681 | 0.692864478 | 0.499796063 | 0.366554111 | 2.75644827 | 0.00267891313 | 0.00918484501 | 0.000149144508 |
| p9016_pcgamma1_td0p5_normdir0p03_msep1_lsep0p5_lct0p03_gct0p003 | 0.0299999993 | 0.0299999993 | 0.00300000003 | 1 | 0.5 | 0.470477 | 0.562325 | 0.339026 | -0.013891 | 0.537114 | 0.77229 | 0.691283524 | 0.498655647 | 0.804690838 | 2.69480705 | 0.00191350938 | 0.00535782625 | 0.000141902119 |

## Interpretation Boundary

A success requires full-denominator trans top1 accuracy to increase by at least 0.1 over the best035 replay without a major cis/Spearman collapse. If normdir improves copytrack continuity but not trans top1, local homolog-vector direction continuity is not sufficient and the next mechanism should target trans-specific gauge selection rather than local track orientation alone.
