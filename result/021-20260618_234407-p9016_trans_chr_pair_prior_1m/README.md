# 021 P9016 Trans Chromosome-Pair Prior 1Mb

This controlled blind-training experiment tests whether a posterior-derived trans chromosome-pair prior can recover the pair-specific gauge space identified by the eval-only oracle diagnostic.

- full result root: `/mnt/ssd/zliu/phase3/test_res/021-20260618_234407-p9016_trans_chr_pair_prior_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/021-20260618_234407-p9016_trans_chr_pair_prior_1m`
- headline: `NO_PLUS_0P1`
- training input: raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper

## Main Results

| config_name | d_scale_posterior_gamma | min_sep_unit | lambda_sep | lambda_copytrack | trans_chr_pair_prior_lambda | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0 | 1 | 1 | 0.5 | 0.0299999993 | 0 | 0.476197 | 0.569731 | 0.342332 | 0.539343 | 0.797915 | 0.595835328 | 0.429804355 |
| p9016_pcgamma1_prior0p5_sep_off | 1 | 0 | 0 | 0 | 0.5 | 0.470837 | 0.561578 | 0.340971 | 0.534558 | 0.813097 | 0.394255072 | 0.28439492 |
| p9016_pcgamma1_prior0p25_sep_off | 1 | 0 | 0 | 0 | 0.25 | 0.462152 | 0.557428 | 0.325795 | 0.531835 | 0.803608 | 0.492683023 | 0.355395675 |
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p5 | 1 | 1 | 0.5 | 0.0299999993 | 0.5 | 0.462029 | 0.557496 | 0.325399 | 0.519597 | 0.814805 | 0.402104855 | 0.290057331 |
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p25 | 1 | 1 | 0.5 | 0.0299999993 | 0.25 | 0.462561 | 0.564742 | 0.316321 | 0.523264 | 0.799927 | 0.488389373 | 0.352298468 |
| p9016_pcgamma1_baseline_sep_off_copytrack0_prior0 | 1 | 0 | 0 | 0 | 0 | 0.465049 | 0.569634 | 0.31537 | 0.51179 | 0.800612 | 0.607306719 | 0.438079208 |
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p75 | 1 | 1 | 0.5 | 0.0299999993 | 0.75 | 0.443863 | 0.549299 | 0.292963 | 0.509602 | 0.815416 | 0.307506353 | 0.221818954 |
| p9016_pcgamma1_prior0p75_sep_off | 1 | 0 | 0 | 0 | 0.75 | 0.438452 | 0.542287 | 0.289845 | 0.499385 | 0.806051 | 0.310993999 | 0.224334747 |

## Interpretation Boundary

The trans chromosome-pair prior is blind because it is estimated from previous posterior probabilities and raw contact counts only. It does not use phase labels or CHARM/3DG during training. If the +0.1 target is not met, the eval-only pair oracle is not self-bootstrappable from the current posterior family.
