# 031 P9016 Confidence-Boosted Trans Pairs 1Mb

This controlled blind-training experiment tests whether high-confidence trans contacts from a baseline posterior can act as a self-training scaffold.

- full result root: `/mnt/ssd/zliu/phase3/test_res/031-20260619_105814-p9016_confidence_boost_pairs_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/031-20260619_105814-p9016_confidence_boost_pairs_1m`
- summary: `/mnt/ssd/zliu/phase3/test_res/031-20260619_105814-p9016_confidence_boost_pairs_1m/summary.tsv`
- trans delta summary: `/mnt/ssd/zliu/phase3/test_res/031-20260619_105814-p9016_confidence_boost_pairs_1m/trans_delta_summary.tsv`
- headline: `NO_PLUS_0P1`
- training input: P9016 pairs only; boosted pairs duplicate original raw contacts selected by blind baseline posterior confidence
- training boundary: SNP phase labels and CHARM/3DG are not used for training; generated training pairs write phase columns as `.`
- eval-only inputs: original P9016 pairs with SNP labels and CHARM/3DG are used only by the standard eval wrapper after training

## Main Results

| config_name | score_mode | trans_top_frac | boost_copies | min_sep_unit | lambda_sep | model_top1_accuracy_genome_trans | delta_trans_vs_baseline | model_top1_accuracy_genome_cis | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | final_mean_sep |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_confboost_pmax_margin_top0p005_b2_sep0_ieps0p5_noise0_seed17 | pmax_margin | 0.005 | 2 | 0 | 0 | 0.329685 | 0.012093 | 0.570371 | 0.808165 | 0.612581909 | 0.441884458 | 1.58095813 | 3.54065919 | 4.95014334 |
| p9016_pcgamma1_confboost_pmax_top0p02_b2_sep0_ieps0p5_noise0_seed17 | pmax | 0.02 | 2 | 0 | 0 | 0.32817 | 0.010578 | 0.551949 | 0.739462 | 0.609425187 | 0.439607322 | 1.30304587 | 2.9890933 | 4.69488001 |
| p9016_pcgamma1_confboost_pmax_margin_entropy_top0p05_b2_sep0_ieps0p5_noise0_seed17 | pmax_margin_entropy | 0.05 | 2 | 0 | 0 | 0.327254 | 0.009662 | 0.555599 | 0.742854 | 0.604432225 | 0.436005682 | 0.446749866 | 2.94931841 | 4.80511236 |
| p9016_pcgamma1_approved_msep1_lsep1_ieps0p5_noise0_seed17 | NA | NA | NA | 1 | 1 | 0.323489 | 0.005897 | 0.568061 | 0.79463 | 0.59696579 | 0.430619806 | 1.33821261 | 3.58990312 | 4.93882084 |
| p9016_pcgamma1_confboost_pmax_margin_entropy_top0p02_b2_sep0_ieps0p5_noise0_seed17 | pmax_margin_entropy | 0.02 | 2 | 0 | 0 | 0.322961 | 0.005369 | 0.557283 | 0.747853 | 0.61669004 | 0.444847822 | 0.911422431 | 2.95037198 | 4.61535692 |
| p9016_pcgamma1_confboost_pmax_margin_top0p02_b2_sep0_ieps0p5_noise0_seed17 | pmax_margin | 0.02 | 2 | 0 | 0 | 0.321496 | 0.003904 | 0.562403 | 0.781876 | 0.618372977 | 0.44606182 | 1.13790894 | 3.11447215 | 4.70014858 |
| p9016_pcgamma1_confboost_pmax_top0p05_b5_sep0_ieps0p5_noise0_seed17 | pmax | 0.05 | 5 | 0 | 0 | 0.320232 | 0.00264 | 0.54564 | 0.755212 | 0.652778745 | 0.4708803 | 0.926455319 | 2.77900171 | 4.6450119 |
| p9016_pcgamma1_approved_baseline_ieps0p5_noise0_seed17 | NA | NA | NA | 0 | 0 | 0.317592 | 0 | 0.569799 | 0.801886 | 0.611727953 | 0.441268414 | 1.51510656 | 3.66769004 | 5.01181269 |
| p9016_pcgamma1_confboost_pmax_margin_top0p02_b2_msep1_lsep1_ieps0p5_noise0_seed17 | pmax_margin | 0.02 | 2 | 1 | 1 | 0.317328 | -0.000264 | 0.559141 | 0.767741 | 0.627714634 | 0.452800423 | 0.965634406 | 2.90483856 | 4.60769844 |
| p9016_pcgamma1_confboost_pmax_top0p02_b5_sep0_ieps0p5_noise0_seed17 | pmax | 0.02 | 5 | 0 | 0 | 0.316044 | -0.001548 | 0.549052 | 0.763133 | 0.621080041 | 0.448014557 | 0.922621965 | 3.07649088 | 4.98615265 |
| p9016_pcgamma1_confboost_pmax_margin_top0p01_b2_sep0_ieps0p5_noise0_seed17 | pmax_margin | 0.01 | 2 | 0 | 0 | 0.316037 | -0.001555 | 0.556642 | 0.751412 | 0.629301131 | 0.453944832 | 1.17500842 | 2.98041725 | 4.75986052 |
| p9016_pcgamma1_confboost_pmax_margin_top0p01_b5_sep0_ieps0p5_noise0_seed17 | pmax_margin | 0.01 | 5 | 0 | 0 | 0.313543 | -0.004049 | 0.557006 | 0.753231 | 0.628001213 | 0.453007132 | 0.899563432 | 2.72828317 | 4.55162144 |
| p9016_pcgamma1_confboost_pmax_top0p05_b2_sep0_ieps0p5_noise0_seed17 | pmax | 0.05 | 2 | 0 | 0 | 0.311008 | -0.006584 | 0.548115 | 0.755258 | 0.628883719 | 0.453643739 | 1.26932907 | 2.93119526 | 4.60191107 |
| p9016_pcgamma1_confboost_pmax_margin_top0p02_b5_sep0_ieps0p5_noise0_seed17 | pmax_margin | 0.02 | 5 | 0 | 0 | 0.308702 | -0.00889 | 0.546625 | 0.756678 | 0.631068289 | 0.455219537 | 1.4714272 | 3.18092632 | 4.74794674 |
| p9016_pcgamma1_confboost_pmax_margin_top0p05_b2_msep1_lsep1_ieps0p5_noise0_seed17 | pmax_margin | 0.05 | 2 | 1 | 1 | 0.306848 | -0.010744 | 0.537264 | 0.73027 | 0.636308193 | 0.458999336 | 1.3922317 | 3.13709331 | 4.66782475 |
| p9016_pcgamma1_confboost_pmax_margin_entropy_top0p05_b5_sep0_ieps0p5_noise0_seed17 | pmax_margin_entropy | 0.05 | 5 | 0 | 0 | 0.299235 | -0.018357 | 0.538817 | 0.735932 | 0.642053246 | 0.463143498 | 1.43721569 | 3.05684185 | 4.59194803 |
| p9016_pcgamma1_confboost_pmax_margin_top0p02_b5_msep1_lsep1_ieps0p5_noise0_seed17 | pmax_margin | 0.02 | 5 | 1 | 1 | 0.298548 | -0.019044 | 0.544922 | 0.762145 | 0.633439779 | 0.45693022 | 0.950161338 | 3.03735232 | 4.79556894 |
| p9016_pcgamma1_confboost_pmax_margin_entropy_top0p02_b5_sep0_ieps0p5_noise0_seed17 | pmax_margin_entropy | 0.02 | 5 | 0 | 0 | 0.298131 | -0.019461 | 0.539579 | 0.760632 | 0.635539055 | 0.458444506 | 0.319049299 | 2.7222743 | 4.78499794 |
| p9016_pcgamma1_confboost_pmax_margin_top0p1_b5_sep0_ieps0p5_noise0_seed17 | pmax_margin | 0.1 | 5 | 0 | 0 | 0.295707 | -0.021885 | 0.528761 | 0.72688 | 0.697956562 | 0.503469229 | 0.92659384 | 2.62419438 | 4.20302868 |
| p9016_pcgamma1_confboost_pmax_margin_top0p1_b2_msep1_lsep1_ieps0p5_noise0_seed17 | pmax_margin | 0.1 | 2 | 1 | 1 | 0.294082 | -0.02351 | 0.53394 | 0.722743 | 0.646851718 | 0.466604888 | 1.10235345 | 3.1242311 | 4.59238863 |
| p9016_pcgamma1_confboost_pmax_margin_top0p05_b5_sep0_ieps0p5_noise0_seed17 | pmax_margin | 0.05 | 5 | 0 | 0 | 0.294054 | -0.023538 | 0.531479 | 0.73203 | 0.651181877 | 0.46972844 | 1.24166203 | 2.83278298 | 4.50050592 |
| p9016_pcgamma1_confboost_pmax_margin_top0p05_b2_sep0_ieps0p5_noise0_seed17 | pmax_margin | 0.05 | 2 | 0 | 0 | 0.291956 | -0.025636 | 0.543854 | 0.739751 | 0.638324976 | 0.460454136 | 1.30524433 | 3.0872221 | 4.64400434 |
| p9016_pcgamma1_confboost_pmax_margin_top0p05_b5_msep1_lsep1_ieps0p5_noise0_seed17 | pmax_margin | 0.05 | 5 | 1 | 1 | 0.291477 | -0.026115 | 0.540797 | 0.743593 | 0.641856551 | 0.463001609 | 1.3190912 | 2.84231377 | 4.56488895 |
| p9016_pcgamma1_confboost_pmax_margin_top0p1_b2_sep0_ieps0p5_noise0_seed17 | pmax_margin | 0.1 | 2 | 0 | 0 | 0.29063 | -0.026962 | 0.529737 | 0.704715 | 0.653384328 | 0.471317172 | 1.13342643 | 2.95843148 | 4.47316217 |
| p9016_pcgamma1_confboost_pmax_margin_top0p2_b2_sep0_ieps0p5_noise0_seed17 | pmax_margin | 0.2 | 2 | 0 | 0 | 0.289164 | -0.028428 | 0.513081 | 0.759403 | 0.645794034 | 0.465841949 | 1.23188651 | 3.17833638 | 4.51755333 |

## Interpretation Boundary

A real success requires full-denominator trans top1 accuracy to increase by at least 0.1 over the approved baseline without collapsing cis accuracy or cis distance Spearman. This experiment tests blind self-training geometry only; it does not use oracle posterior labels, CHARM/3DG, or SNP truth during training.
