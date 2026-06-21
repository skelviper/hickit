# 035 P9016 Copytrack Entropy Combo 1Mb

This controlled blind-training experiment combines the best 034 trans dscale setting with local/global homolog-vector continuity and entropy-floor rho training.

- full result root: `/mnt/ssd/zliu/phase3/test_res/035-20260619_143807-p9016_copytrack_entropy_combo_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/035-20260619_143807-p9016_copytrack_entropy_combo_1m`
- headline: `NO_PLUS_0P1`
- training input: raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper
- changed knobs: trans_dscale_multiplier, lambda_copytrack, lambda_global_copytrack, optional entropy_with_floor rho, optional homolog separation guardrail
- unchanged model settings: uniform prior, posterior_count gamma1, no chromosome-pair prior, no hard top1 M-step, no phase/CHARM training labels

## Main Results

| config_name | rho_train_mode | rho_train_floor | trans_dscale_multiplier | lambda_copytrack | lambda_global_copytrack | min_sep_unit | lambda_sep | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | delta_trans_vs_baseline | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | copytrack_frac_cos_lt_0 | copytrack_frac_projection_sign_switch |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_td0p5_sep0_lct0p03_gct0p003 | constant | 0 | 0.5 | 0.0299999993 | 0.00300000003 | 0 | 0 | 0.481979 | 0.572094 | 0.353008 | 0.036742 | 0.552637 | 0.769483 | 0.696803987 | 0.502637804 | 0.460427523 | 2.83314991 | 0.0015308075 | 0.0045924225 |
| p9016_pcgamma1_td0p5_entropyfloor0p5_sep0_lct0p01_gct0 | entropy_with_floor | 0.5 | 0.5 | 0.00999999978 | 0 | 0 | 0 | 0.478911 | 0.568003 | 0.351403 | 0.035137 | 0.542295 | 0.763052 | 0.698627293 | 0.503953099 | 0.352698177 | 2.76446962 | 0.0015308075 | 0.00612323 |
| p9016_pcgamma1_td0p5_sep0_lct0p01_gct0 | constant | 0 | 0.5 | 0.00999999978 | 0 | 0 | 0 | 0.479282 | 0.568867 | 0.35107 | 0.034804 | 0.541906 | 0.762935 | 0.698452532 | 0.503827035 | 0.352800637 | 2.75856352 | 0.0015308075 | 0.00612323 |
| p9016_pcgamma1_td0p5_sep0_lct0_gct0p01 | constant | 0 | 0.5 | 0 | 0.00999999978 | 0 | 0 | 0.479911 | 0.570012 | 0.350959 | 0.034693 | 0.551588 | 0.758948 | 0.705884993 | 0.509188354 | 0.348994702 | 2.64297223 | 0.00344431688 | 0.0110983544 |
| p9016_pcgamma1_td0p5_sep0_lct0p01_gct0p003 | constant | 0 | 0.5 | 0.00999999978 | 0.00300000003 | 0 | 0 | 0.476788 | 0.565392 | 0.34998 | 0.033714 | 0.54872 | 0.757005 | 0.703584075 | 0.507528603 | 0.319390893 | 2.69119954 | 0.00191350938 | 0.00382701875 |
| p9016_pcgamma1_td0p5_entropyfloor0p5_sep0_lct0p01_gct0p003 | entropy_with_floor | 0.5 | 0.5 | 0.00999999978 | 0.00300000003 | 0 | 0 | 0.476608 | 0.56512 | 0.349931 | 0.033665 | 0.548553 | 0.757059 | 0.703506231 | 0.507472515 | 0.319245249 | 2.69301796 | 0.00191350938 | 0.00382701875 |
| p9016_pcgamma1_td0p5_sep0_lct0p03_gct0 | constant | 0 | 0.5 | 0.0299999993 | 0 | 0 | 0 | 0.477702 | 0.568702 | 0.347465 | 0.031199 | 0.539267 | 0.762043 | 0.702229738 | 0.506551683 | 0.36555478 | 2.67662477 | 0.0015308075 | 0.00382701875 |
| p9016_pcgamma1_td0p5_entropyfloor0p5_sep0_lct0_gct0 | entropy_with_floor | 0.5 | 0.5 | 0 | 0 | 0 | 0 | 0.475437 | 0.566867 | 0.344583 | 0.028317 | 0.538447 | 0.759551 | 0.700904727 | 0.505595922 | 0.302319199 | 2.71229839 | 0.0015308075 | 0.00420972063 |
| p9016_pcgamma1_td0p5_best034_sep0_lct0_gct0 | constant | 0 | 0.5 | 0 | 0 | 0 | 0 | 0.474368 | 0.565261 | 0.344284 | 0.028018 | 0.538697 | 0.757868 | 0.701750696 | 0.506206155 | 0.30865249 | 2.70723104 | 0.00191350938 | 0.00574052813 |
| p9016_pcgamma1_td0p5_sep0_lct0_gct0p003 | constant | 0 | 0.5 | 0 | 0.00300000003 | 0 | 0 | 0.474585 | 0.566756 | 0.342673 | 0.026407 | 0.539572 | 0.760717 | 0.703910053 | 0.507763743 | 0.367391467 | 2.74015927 | 0.00114810563 | 0.00612323 |
| p9016_pcgamma1_td0p5_entropyfloor0p5_sep0_lct0_gct0p003 | entropy_with_floor | 0.5 | 0.5 | 0 | 0.00300000003 | 0 | 0 | 0.474665 | 0.567018 | 0.342492 | 0.026226 | 0.539149 | 0.762797 | 0.703963101 | 0.50780201 | 0.3746261 | 2.71693945 | 0.00191350938 | 0.00535782625 |
| p9016_pcgamma1_td0p5_entropyfloor0p5_msep1_lsep0p5_lct0p01_gct0p003 | entropy_with_floor | 0.5 | 0.5 | 0.00999999978 | 0.00300000003 | 1 | 0.5 | 0.471851 | 0.567823 | 0.334498 | 0.018232 | 0.529112 | 0.766508 | 0.694540024 | 0.501004696 | 0.888296068 | 2.70557547 | 0.00229621125 | 0.00497512438 |
| p9016_pcgamma1_td0p5_msep1_lsep0p5_lct0p01_gct0p003 | constant | 0 | 0.5 | 0.00999999978 | 0.00300000003 | 1 | 0.5 | 0.472508 | 0.568969 | 0.334456 | 0.01819 | 0.529314 | 0.768668 | 0.694562435 | 0.501020908 | 0.899978101 | 2.71190476 | 0.00229621125 | 0.00497512438 |
| p9016_pcgamma1_td0p5_msep1_lsep1_lct0p01_gct0p003 | constant | 0 | 0.5 | 0.00999999978 | 0.00300000003 | 1 | 1 | 0.469789 | 0.567183 | 0.3304 | 0.014134 | 0.535974 | 0.766694 | 0.672850847 | 0.485359281 | 0.887782514 | 2.7631278 | 0.00114810563 | 0.003061615 |
| p9016_pcgamma1_td1_baseline_sep0_lct0_gct0 | constant | 0 | 1 | 0 | 0 | 0 | 0 | 0.464469 | 0.568022 | 0.316266 | 0 | 0.51352 | 0.801884 | 0.611757576 | 0.441289783 | 1.54554296 | 3.63843393 | 0.00191350938 | 0.003061615 |

## Interpretation Boundary

A real success requires full-denominator trans top1 accuracy to increase by at least 0.1 over the posterior-count baseline without major cis/Spearman collapse. If copytrack discontinuity improves without trans accuracy improvement, the remaining problem is likely cross-chromosome gauge selection rather than local track continuity alone.
# 035 Raw-Pairs Trans Improvement Audit
- baseline config: `p9016_pcgamma1_td1_baseline_sep0_lct0_gct0`
- baseline trans top1: `0.316266`
- best config: `p9016_pcgamma1_td0p5_sep0_lct0p03_gct0p003`
- best trans top1: `0.353008`
- delta: `0.036742`
- target +0.1 met: `no`

## Interpretation

The strongest raw-pairs blind condition so far combines posterior-count gamma1, trans dscale 0.5, local copytrack 0.03, and global copytrack 0.003. It improves trans modestly but does not reach the required +0.1. This argues that force/rho/copytrack guardrails are not enough; the remaining bottleneck is likely cross-chromosome gauge or copy-state identity selection rather than local geometry alone.
