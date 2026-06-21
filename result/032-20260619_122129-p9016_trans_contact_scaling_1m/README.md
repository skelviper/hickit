# 032 P9016 Trans Contact Scaling 1Mb

This controlled blind-training experiment tests whether trans-specific M-step edge strength or target-distance scaling can improve four-state trans accuracy.

- full result root: `/tmp/hk_blind_test_res_032_smoke/032-20260619_122129-p9016_trans_contact_scaling_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/032-20260619_122129-p9016_trans_contact_scaling_1m`
- headline: `NO_PLUS_0P1`
- training input: raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper
- unchanged model settings: uniform prior, constant rho_train, temperature 1, posterior_count gamma 1, no chromosome-pair prior, no copytrack force

## Main Results

| config_name | trans_k_multiplier | trans_dscale_multiplier | min_sep_unit | lambda_sep | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | final_mean_sep |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_transscale_smoke_tk1_td1 | 1 | 1 | 0 | 0 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 |
| p9016_transscale_smoke_tk2_td0p75 | 2 | 0.75 | 0 | 0 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 1.00231004 | 1.00237644 |

## Interpretation Boundary

A real success requires full-denominator trans top1 accuracy to increase by at least 0.1 over the no-scaling baseline without a major cis accuracy or cis distance-Spearman collapse. These knobs change only trans M-step contact edges; they do not solve copy-gauge synchronization directly.
