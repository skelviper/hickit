# 021 P9016 Trans Chromosome-Pair Prior 1Mb

This controlled blind-training experiment tests whether a posterior-derived trans chromosome-pair prior can recover the pair-specific gauge space identified by the eval-only oracle diagnostic.

- full result root: `/tmp/hk_blind_test_res_021_smoke/021-smoke-p9016_trans_chr_pair_prior_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/021-smoke-p9016_trans_chr_pair_prior_1m`
- headline: `NO_PLUS_0P1`
- training input: raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper

## Main Results

| config_name | d_scale_posterior_gamma | min_sep_unit | lambda_sep | lambda_copytrack | trans_chr_pair_prior_lambda | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma0p5_msep1_lsep0p5_copytrack0p01_prior0p5 | 0.5 | 1 | 0.5 | 0.00999999978 | 0.5 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |
| p9016_pcgamma0p5_prior0p5_sep_off | 0.5 | 0 | 0 | 0 | 0.5 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |
| p9016_pcgamma1_baseline_sep_off_copytrack0_prior0 | 1 | 0 | 0 | 0 | 0 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0 | 1 | 1 | 0.5 | 0.0299999993 | 0 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p25 | 1 | 1 | 0.5 | 0.0299999993 | 0.25 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p5 | 1 | 1 | 0.5 | 0.0299999993 | 0.5 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p75 | 1 | 1 | 0.5 | 0.0299999993 | 0.75 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior1 | 1 | 1 | 0.5 | 0.0299999993 | 1 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |
| p9016_pcgamma1_prior0p25_sep_off | 1 | 0 | 0 | 0 | 0.25 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |
| p9016_pcgamma1_prior0p5_sep_off | 1 | 0 | 0 | 0 | 0.5 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |
| p9016_pcgamma1_prior0p75_sep_off | 1 | 0 | 0 | 0 | 0.75 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |
| p9016_pcgamma1_prior1_sep_off | 1 | 0 | 0 | 0 | 1 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 |

## Interpretation Boundary

The trans chromosome-pair prior is blind because it is estimated from previous posterior probabilities and raw contact counts only. It does not use phase labels or CHARM/3DG during training. If the +0.1 target is not met, the eval-only pair oracle is not self-bootstrappable from the current posterior family.
