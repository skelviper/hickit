# 022 P9016 Trans Pair Sharpen/Warmup 1Mb

This controlled blind-training experiment tests whether a posterior-derived trans chromosome-pair prior needs warmup and power sharpening to escape the symmetric trans fixed point.

- full result root: `/mnt/ssd/zliu/phase3/test_res/022-20260619_004818-p9016_trans_pair_sharpen_warmup_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/022-20260619_004818-p9016_trans_pair_sharpen_warmup_1m`
- headline: `NO_PLUS_0P1`
- training input: raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper

## Main Results

| config_name | d_scale_posterior_gamma | min_sep_unit | lambda_sep | lambda_copytrack | trans_chr_pair_prior_lambda | trans_chr_pair_prior_power | trans_chr_pair_prior_warmup_iter | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0 | 1 | 1 | 0.5 | 0.0299999993 | 0 | 1 | 0 | 0.476197 | 0.569731 | 0.342332 | 0.539343 | 0.797915 | 0.595835328 | 0.429804355 |
| p9016_pcgamma1_prior0p5_power1_warmup0_sep_off | 1 | 0 | 0 | 0 | 0.5 | 1 | 0 | 0.470837 | 0.561578 | 0.340971 | 0.534558 | 0.813097 | 0.394255072 | 0.28439492 |
| p9016_pcgamma1_prior0p35_power4_warmup10_sep_off | 1 | 0 | 0 | 0 | 0.349999994 | 4 | 10 | 0.461129 | 0.553221 | 0.32933 | 0.53055 | 0.821153 | 0.434677601 | 0.313553602 |
| p9016_pcgamma1_prior0p5_power2_warmup10_sep_off | 1 | 0 | 0 | 0 | 0.5 | 2 | 10 | 0.458478 | 0.552066 | 0.324538 | 0.526626 | 0.813624 | 0.36807999 | 0.265513599 |
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p5_power4_warmup10 | 1 | 1 | 0.5 | 0.0299999993 | 0.5 | 4 | 10 | 0.46091 | 0.560855 | 0.31787 | 0.529307 | 0.809986 | 0.376300067 | 0.271443099 |
| p9016_pcgamma1_prior0p65_power4_warmup10_sep_off | 1 | 0 | 0 | 0 | 0.649999976 | 4 | 10 | 0.450962 | 0.54564 | 0.31546 | 0.516207 | 0.818329 | 0.307407916 | 0.221747935 |
| p9016_pcgamma1_baseline_sep_off_copytrack0_prior0 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0.465049 | 0.569634 | 0.31537 | 0.51179 | 0.800612 | 0.607306719 | 0.438079208 |
| p9016_pcgamma1_prior0p5_power4_warmup20_sep_off | 1 | 0 | 0 | 0 | 0.5 | 4 | 20 | 0.452636 | 0.553531 | 0.308237 | 0.515818 | 0.813577 | 0.392216802 | 0.282924622 |
| p9016_pcgamma0p5_prior0p5_power4_warmup10_sep_off | 0.5 | 0 | 0 | 0 | 0.5 | 4 | 10 | 0.452899 | 0.555269 | 0.306389 | 0.520917 | 0.824543 | 0.355875105 | 0.256709635 |
| p9016_pcgamma1_prior0p5_power8_warmup10_sep_off | 1 | 0 | 0 | 0 | 0.5 | 8 | 10 | 0.445414 | 0.546975 | 0.300062 | 0.512193 | 0.815216 | 0.372925282 | 0.269008726 |
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p5_power2_warmup10 | 1 | 1 | 0.5 | 0.0299999993 | 0.5 | 2 | 10 | 0.448325 | 0.552051 | 0.299874 | 0.506643 | 0.807221 | 0.373840153 | 0.269668669 |
| p9016_pcgamma1_prior0p5_power4_warmup10_sep_off | 1 | 0 | 0 | 0 | 0.5 | 4 | 10 | 0.443454 | 0.546116 | 0.296527 | 0.511318 | 0.814129 | 0.372400433 | 0.268630147 |

## Interpretation Boundary

The sharpened chromosome-pair prior is blind because it is estimated from previous posterior probabilities and raw contact counts only. It does not use phase labels or CHARM/3DG during training. If the +0.1 target is not met, the eval-only pair oracle is not self-bootstrappable from the current posterior family.
