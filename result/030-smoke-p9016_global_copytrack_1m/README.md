# 030 P9016 Global Copytrack 1Mb

This controlled blind-training experiment tests whether local and chromosome-mean homolog-vector gauge alignment can improve P9016 trans four-state accuracy.

- full result root: `/tmp/hk_blind_test_res_030_smoke/030-smoke-p9016_global_copytrack_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/030-smoke-p9016_global_copytrack_1m`
- summary: `/tmp/hk_blind_test_res_030_smoke/030-smoke-p9016_global_copytrack_1m/summary.tsv`
- trans delta summary: `/tmp/hk_blind_test_res_030_smoke/030-smoke-p9016_global_copytrack_1m/trans_delta_summary.tsv`
- headline: `NO_PLUS_0P1`
- training input: approved raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper
- backend: CPU for this experiment, because GPU extra-force semantics use operator splitting and are not used for conclusions here

## Main Results

| config_name | lambda_copytrack | lambda_global_copytrack | min_sep_unit | lambda_sep | model_top1_accuracy_genome_trans | delta_trans_vs_baseline | model_top1_accuracy_genome_cis | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | final_mean_sep | copytrack_frac_cos_lt_0 | global_copytrack_force_l1_over_contact_force_l1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_sep0_lct0_gct0_ieps0p5_noise0_seed17 | 0 | 0 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0_gct0p001_ieps0p5_noise0_seed17 | 0 | 0.00100000005 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0_gct0p003_ieps0p5_noise0_seed17 | 0 | 0.00300000003 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0_gct0p01_ieps0p5_noise0_seed17 | 0 | 0.00999999978 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0_gct0p03_ieps0p5_noise0_seed17 | 0 | 0.0299999993 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0_gct0p1_ieps0p5_noise0_seed17 | 0 | 0.100000001 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0p001_gct0_ieps0p5_noise0_seed17 | 0.00100000005 | 0 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0p003_gct0_ieps0p5_noise0_seed17 | 0.00300000003 | 0 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0p003_gct0p003_ieps0p5_noise0_seed17 | 0.00300000003 | 0.00300000003 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0p003_gct0p01_ieps0p5_noise0_seed17 | 0.00300000003 | 0.00999999978 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0p01_gct0_ieps0p5_noise0_seed17 | 0.00999999978 | 0 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0p01_gct0p003_ieps0p5_noise0_seed17 | 0.00999999978 | 0.00300000003 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0p01_gct0p01_ieps0p5_noise0_seed17 | 0.00999999978 | 0.00999999978 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0p01_gct0p03_ieps0p5_noise0_seed17 | 0.00999999978 | 0.0299999993 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0p03_gct0_ieps0p5_noise0_seed17 | 0.0299999993 | 0 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0p03_gct0p01_ieps0p5_noise0_seed17 | 0.0299999993 | 0.00999999978 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |
| p9016_pcgamma1_sep0_lct0p1_gct0_ieps0p5_noise0_seed17 | 0.100000001 | 0 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 | 0 | NA |

## Interpretation Boundary

Global copytrack is an internal gauge-alignment force on chromosome-mean copy1-minus-copy0 vectors. It does not assign maternal/paternal identity and does not use phase labels or CHARM/3DG during training. A real success requires trans top1 to increase by at least 0.1 without a major cis/Spearman collapse.
