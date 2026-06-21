# 043 P9016 Raw-Pairs Remaining Ablation 1Mb

This controlled blind-training experiment tests the last small set of raw-pairs-only mechanisms that remained after prior P9016 trans experiments: coarse-to-fine 4Mb-to-1Mb initialization path, weak trans chromosome-pair E-step prior, and a mild pmax+margin trans M-step gate.

- full result root: `/mnt/ssd/zliu/phase3/test_res/043-20260620_031947-p9016_raw_pairs_remaining_ablation_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/043-20260620_031947-p9016_raw_pairs_remaining_ablation_1m`
- headline: `NO_PLUS_0P1`
- training input: approved P9016 raw pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG reference are used only after training through the standard eval wrapper
- copy0/copy1 are gauge labels; this run does not use maternal/paternal labels in training
- success criterion: full-denominator trans top1 accuracy must improve by at least 0.1 over the td1 baseline
- 4Mb intermediate chain stages are retained in `summary.tsv`; the main table below reports only final 1Mb evaluated stages.
- post-run reporting note: `POST_RUN_REPORTING_NOTE.md` records the review-driven reporting cleanup; training and eval outputs were not modified.

## Why These Configs

- `chain4m1m`: tests whether a coarse-to-fine basin changes trans copy-gauge identity when local copytrack alone does not.
- `chrpair_estep_l0p10_w20`: tests a weak posterior-derived chromosome-pair E-step prior without changing M-step graph topology.
- `pmaxmargin_p045_m002`: fills the mild pmax+margin gate hole left after stronger confidence gates failed.

## Main Results

| config_name | trans_chr_pair_prior_lambda | trans_chr_pair_prior_warmup_iter | trans_gate_mode | trans_gate_min_pmax | trans_gate_min_margin | trans_dscale_multiplier | lambda_copytrack | lambda_global_copytrack | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | delta_trans_vs_td1_baseline | delta_trans_vs_best035 | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | copytrack_frac_cos_lt_0 | copytrack_frac_projection_sign_switch |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_rawremain_best035_td0p5_lct0p03_gct0p003 | 0 | 0 | off | 1 | 1 | 0.5 | 0.0299999993 | 0.00300000003 | 0.482056 | 0.572337 | 0.352848 | 0.037541 | 0 | 0.552651 | 0.769549 | 0.696871698 | 0.502686679 | 0.458690643 | 2.83252883 | 0.0015308075 | 0.0045924225 |
| p9016_rawremain_chrpair_estep_l0p10_w20_best035 | 0.100000001 | 20 | off | 1 | 1 | 0.5 | 0.0299999993 | 0.00300000003 | 0.484876 | 0.577762 | 0.351938 | 0.036631 | -0.00091 | 0.54822 | 0.782694 | 0.666954219 | 0.481105775 | 0.846141398 | 2.88482285 | 0.0015308075 | 0.00420972063 |
| p9016_rawremain_pmaxmargin_p045_m002_best035 | 0 | 0 | pmax_margin | 0.449999988 | 0.0199999996 | 0.5 | 0.0299999993 | 0.00300000003 | 0.478762 | 0.580703 | 0.332866 | 0.017559 | -0.019982 | 0.535718 | 0.802833 | 0.657553792 | 0.474324793 | 0.812747002 | 3.00020385 | 0.00229621125 | 0.00267891313 |
| p9016_rawremain_td1_baseline | 0 | 0 | off | 1 | 1 | 1 | 0 | 0 | 0.464812 | 0.569275 | 0.315307 | 0 | -0.037541 | 0.511769 | 0.80079 | 0.60726285 | 0.438047558 | 1.52262414 | 3.72325063 | 0.00191350938 | 0.003061615 |
| p9016_rawremain_chain4m1m_chrpair_l0p10_w20_best035 | 0.100000001 | 20 | off | 1 | 1 | 0.5 | 0.0299999993 | 0.00300000003 | 0.427596 | 0.518983 | 0.296804 | -0.018503 | -0.056044 | 0.514089 | 0.652946 | 0.651883364 | 0.470234424 | 0.482694596 | 2.39523578 | 0.00382701875 | 0.0110983544 |
| p9016_rawremain_chain4m1m_best035 | 0 | 0 | off | 1 | 1 | 0.5 | 0.0299999993 | 0.00300000003 | 0.415271 | 0.505486 | 0.286157 | -0.02915 | -0.066691 | 0.505789 | 0.609572 | 0.658562362 | 0.475052327 | 0.715660155 | 2.41946864 | 0.0168388825 | 0.0221967088 |

## Interpretation Boundary

This is not a new baseline unless the full-denominator trans top1 criterion is met without degrading cis accuracy or cis distance Spearman. In this run, best035 remains the best raw-pairs-only condition; weak chromosome-pair E-step prior, mild pmax+margin gate, and 4Mb-to-1Mb chain do not recover the missing +0.1.
