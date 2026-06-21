# 031 P9016 Confidence-Boosted Trans Pairs 1Mb

This controlled blind-training experiment tests whether high-confidence trans contacts from a baseline posterior can act as a self-training scaffold.

- full result root: `/tmp/hk_blind_test_res_031_smoke2/031-20260619_105804-p9016_confidence_boost_pairs_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/031-20260619_105804-p9016_confidence_boost_pairs_1m`
- summary: `/tmp/hk_blind_test_res_031_smoke2/031-20260619_105804-p9016_confidence_boost_pairs_1m/summary.tsv`
- trans delta summary: `/tmp/hk_blind_test_res_031_smoke2/031-20260619_105804-p9016_confidence_boost_pairs_1m/trans_delta_summary.tsv`
- headline: `NO_PLUS_0P1`
- training input: P9016 pairs only; boosted pairs duplicate original raw contacts selected by blind baseline posterior confidence
- training boundary: SNP phase labels and CHARM/3DG are not used for training; generated training pairs write phase columns as `.`
- eval-only inputs: original P9016 pairs with SNP labels and CHARM/3DG are used only by the standard eval wrapper after training

## Main Results

| config_name | score_mode | trans_top_frac | boost_copies | min_sep_unit | lambda_sep | model_top1_accuracy_genome_trans | delta_trans_vs_baseline | model_top1_accuracy_genome_cis | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | final_mean_sep |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_approved_baseline_ieps0p5_noise0_seed17 | NA | NA | NA | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 |
| p9016_pcgamma1_approved_msep1_lsep1_ieps0p5_noise0_seed17 | NA | NA | NA | 1 | 1 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 |
| p9016_pcgamma1_confboost_pmax_margin_top0p01_b2_sep0_ieps0p5_noise0_seed17 | pmax_margin | 0.01 | 2 | 0 | 0 | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 |

## Interpretation Boundary

A real success requires full-denominator trans top1 accuracy to increase by at least 0.1 over the approved baseline without collapsing cis accuracy or cis distance Spearman. This experiment tests blind self-training geometry only; it does not use oracle posterior labels, CHARM/3DG, or SNP truth during training.
