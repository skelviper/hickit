# 032 P9016 Trans Contact Scaling 1Mb

This controlled blind-training experiment tests whether trans-specific M-step edge strength or target-distance scaling can improve four-state trans accuracy.

- full result root: `/mnt/ssd/zliu/phase3/test_res/032-20260619_122253-p9016_trans_contact_scaling_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/032-20260619_122253-p9016_trans_contact_scaling_1m`
- headline: `NO_PLUS_0P1`
- training input: raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper
- unchanged model settings: uniform prior, constant rho_train, temperature 1, posterior_count gamma 1, no chromosome-pair prior, no copytrack force

## Main Results

| config_name | trans_k_multiplier | trans_dscale_multiplier | min_sep_unit | lambda_sep | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | final_mean_sep |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_tk1_td0p5_sep0_ieps0p5_noise0_seed17 | 1 | 0.5 | 0 | 0 | 0.476034 | 0.567358 | 0.345333 | 0.539267 | 0.760296 | 0.699071288 | 0.504273295 | 0.315692484 | 2.7480166 | 4.37885857 |
| p9016_pcgamma1_tk1_td1p25_sep0_ieps0p5_noise0_seed17 | 1 | 1.25 | 0 | 0 | 0.469886 | 0.561738 | 0.338429 | 0.541587 | 0.796829 | 0.570283711 | 0.411372751 | 1.08888483 | 3.72182846 | 5.29072094 |
| p9016_pcgamma1_tk3_td1_sep0_ieps0p5_noise0_seed17 | 3 | 1 | 0 | 0 | 0.455273 | 0.543932 | 0.328386 | 0.529967 | 0.767364 | 0.697331846 | 0.503018618 | 0.977159083 | 3.18713737 | 4.462327 |
| p9016_pcgamma1_tk2_td0p5_sep0_ieps0p5_noise0_seed17 | 2 | 0.5 | 0 | 0 | 0.456827 | 0.548377 | 0.325802 | 0.526299 | 0.682825 | 0.772300422 | 0.557097018 | 0.512166083 | 2.25653887 | 3.89710522 |
| p9016_pcgamma1_tk1_td1_msep1_lsep1_ieps0p5_noise0_seed17 | 1 | 1 | 1 | 1 | 0.467915 | 0.567741 | 0.325045 | 0.527487 | 0.794501 | 0.596827805 | 0.430520266 | 1.35799825 | 3.58317113 | 4.93934774 |
| p9016_pcgamma1_tk1p5_td0p75_sep0_ieps0p5_noise0_seed17 | 1.5 | 0.75 | 0 | 0 | 0.465863 | 0.564902 | 0.324121 | 0.53146 | 0.771768 | 0.694347382 | 0.500865757 | 0.813229501 | 2.86855555 | 4.36475754 |
| p9016_pcgamma1_tk1_td0p75_sep0_ieps0p5_noise0_seed17 | 1 | 0.75 | 0 | 0 | 0.462841 | 0.561141 | 0.322156 | 0.514645 | 0.783085 | 0.651093245 | 0.469664484 | 1.10920179 | 2.96020865 | 4.62464523 |
| p9016_pcgamma1_tk2_td0p75_sep0_ieps0p5_noise0_seed17 | 2 | 0.75 | 0 | 0 | 0.448396 | 0.537943 | 0.320239 | 0.53105 | 0.756175 | 0.714820385 | 0.515633881 | 1.07665956 | 2.94707084 | 4.35529661 |
| p9016_pcgamma1_tk1_td0p9_sep0_ieps0p5_noise0_seed17 | 1 | 0.899999976 | 0 | 0 | 0.466638 | 0.569143 | 0.319933 | 0.53328 | 0.792774 | 0.61406517 | 0.442954391 | 1.31809092 | 3.48826861 | 4.92637205 |
| p9016_pcgamma1_tk2_td1_sep0_ieps0p5_noise0_seed17 | 2 | 1 | 0 | 0 | 0.456399 | 0.551983 | 0.3196 | 0.517048 | 0.776025 | 0.66403687 | 0.479001343 | 1.41076636 | 3.30917549 | 4.57871914 |
| p9016_pcgamma1_tk1p5_td1_sep0_ieps0p5_noise0_seed17 | 1.5 | 1 | 0 | 0 | 0.456893 | 0.55325 | 0.318988 | 0.521319 | 0.781341 | 0.637193799 | 0.459638149 | 1.86375415 | 3.26747608 | 4.752388 |
| p9016_pcgamma1_tk1_td1_sep0_ieps0p5_noise0_seed17 | 1 | 1 | 0 | 0 | 0.464792 | 0.56925 | 0.315293 | 0.511755 | 0.800822 | 0.607318997 | 0.438088059 | 1.51055276 | 3.72313046 | 5.01565552 |
| p9016_pcgamma1_tk3_td0p75_sep0_ieps0p5_noise0_seed17 | 3 | 0.75 | 0 | 0 | 0.44082 | 0.531023 | 0.311723 | 0.51661 | 0.719835 | 0.744914889 | 0.537342489 | 1.18458843 | 2.6116693 | 4.06597757 |
| p9016_pcgamma1_tk2_td0p75_msep1_lsep1_ieps0p5_noise0_seed17 | 2 | 0.75 | 1 | 1 | 0.450565 | 0.547912 | 0.311244 | 0.508039 | 0.758491 | 0.718498468 | 0.518287063 | 1.00033987 | 2.82052922 | 4.23078489 |
| p9016_pcgamma1_tk1p25_td1_sep0_ieps0p5_noise0_seed17 | 1.25 | 1 | 0 | 0 | 0.453062 | 0.553861 | 0.308799 | 0.523855 | 0.784836 | 0.626615703 | 0.452007681 | 2.13693619 | 3.5486033 | 4.89378881 |

## Interpretation Boundary

A real success requires full-denominator trans top1 accuracy to increase by at least 0.1 over the no-scaling baseline without a major cis accuracy or cis distance-Spearman collapse. These knobs change only trans M-step contact edges; they do not solve copy-gauge synchronization directly.
