# 036 P9016 Trans Confidence Rho 1Mb

This controlled blind-training experiment tests whether low-confidence trans bpair forces are harming the M-step. The new rho modes keep cis bpair training weight at full strength and downweight trans bpair weight by posterior confidence only during blind training.

- full result root: `/mnt/ssd/zliu/phase3/test_res/036-20260619_212753-p9016_trans_confidence_rho_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/036-20260619_212753-p9016_trans_confidence_rho_1m`
- headline: `NO_PLUS_0P1`
- training input: raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper
- changed knobs: rho_train_mode in {trans_entropy, trans_entropy_with_floor}, rho floor, trans dscale, copytrack/global-copytrack combinations
- unchanged model settings: uniform prior, posterior_count gamma1, no chromosome-pair prior, no hard top1 M-step, no phase/CHARM training labels

## Main Results

| config_name | rho_train_mode | rho_train_floor | trans_dscale_multiplier | lambda_copytrack | lambda_global_copytrack | min_sep_unit | lambda_sep | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | delta_trans_vs_baseline | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | copytrack_frac_cos_lt_0 | copytrack_frac_projection_sign_switch |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | constant | 0 | 0.5 | 0.0299999993 | 0.00300000003 | 0 | 0 | 0.482013 | 0.57223 | 0.352897 | 0.037611 | 0.552567 | 0.769441 | 0.696904838 | 0.502710581 | 0.459362566 | 2.82940531 | 0.0015308075 | 0.0045924225 |
| p9016_pcgamma1_td0p5_transfloor0p25_lct0p03_gct0p003 | trans_entropy_with_floor | 0.25 | 0.5 | 0.0299999993 | 0.00300000003 | 0 | 0 | 0.492584 | 0.591564 | 0.350924 | 0.035638 | 0.548977 | 0.766898 | 0.626250327 | 0.451744109 | 0.60813117 | 2.9773519 | 0.00382701875 | 0.00688863375 |
| p9016_pcgamma1_td0p5_transfloor0p25_sep0 | trans_entropy_with_floor | 0.25 | 0.5 | 0 | 0 | 0 | 0 | 0.488221 | 0.585396 | 0.349146 | 0.03386 | 0.547685 | 0.765617 | 0.624772608 | 0.45067817 | 0.529953122 | 3.10519791 | 0.003061615 | 0.0076540375 |
| p9016_pcgamma1_td0p5_transfloor0p5_sep0 | trans_entropy_with_floor | 0.5 | 0.5 | 0 | 0 | 0 | 0 | 0.478562 | 0.571051 | 0.346194 | 0.030908 | 0.552929 | 0.800016 | 0.626305461 | 0.451783866 | 1.43201578 | 3.21148229 | 0.00076540375 | 0.00191350938 |
| p9016_pcgamma1_td0p5_baseline_sep0_lct0_gct0 | constant | 0 | 0.5 | 0 | 0 | 0 | 0 | 0.476034 | 0.567358 | 0.345333 | 0.030047 | 0.539267 | 0.760296 | 0.699071348 | 0.504273355 | 0.315695494 | 2.74802518 | 0.00191350938 | 0.00420972063 |
| p9016_pcgamma1_td1_transfloor0p5_sep0 | trans_entropy_with_floor | 0.5 | 1 | 0 | 0 | 0 | 0 | 0.474805 | 0.573734 | 0.33322 | 0.017934 | 0.537634 | 0.792767 | 0.581895292 | 0.419748753 | 0.94514966 | 3.67826438 | 0.00114810563 | 0.0015308075 |
| p9016_pcgamma1_td0p5_transfloor0p5_lct0p03_gct0p003 | trans_entropy_with_floor | 0.5 | 0.5 | 0.0299999993 | 0.00300000003 | 0 | 0 | 0.471294 | 0.568105 | 0.332741 | 0.017455 | 0.539711 | 0.785175 | 0.629778743 | 0.454289347 | 1.10118473 | 3.01590252 | 0.0015308075 | 0.003061615 |
| p9016_pcgamma1_td0p5_transfloor0p75_lct0p03_gct0p003 | trans_entropy_with_floor | 0.75 | 0.5 | 0.0299999993 | 0.00300000003 | 0 | 0 | 0.478256 | 0.58047 | 0.33197 | 0.016684 | 0.536718 | 0.793719 | 0.650631666 | 0.469331563 | 0.68178916 | 2.76012778 | 0.00191350938 | 0.00344431688 |
| p9016_pcgamma1_td0p5_transfloor0p5_msep1_lsep0p5_lct0p01_gct0p003 | trans_entropy_with_floor | 0.5 | 0.5 | 0.00999999978 | 0.00300000003 | 1 | 0.5 | 0.475779 | 0.577238 | 0.330574 | 0.015288 | 0.535905 | 0.798418 | 0.630545914 | 0.454842716 | 1.03724945 | 3.00380659 | 0.0015308075 | 0.003061615 |
| p9016_pcgamma1_td0p5_transfloor0p1_sep0 | trans_entropy_with_floor | 0.100000001 | 0.5 | 0 | 0 | 0 | 0 | 0.480051 | 0.585653 | 0.328914 | 0.013628 | 0.531877 | 0.765002 | 0.634980679 | 0.458041728 | 0.876264095 | 3.22376037 | 0.00688863375 | 0.0103329506 |
| p9016_pcgamma1_td0p5_transentropy_sep0 | trans_entropy | 0 | 0.5 | 0 | 0 | 0 | 0 | 0.469252 | 0.567896 | 0.328073 | 0.012787 | 0.526334 | 0.74236 | 0.647221923 | 0.466871947 | 0.532030702 | 2.75811791 | 0.0045924225 | 0.00918484501 |
| p9016_pcgamma1_td1_transfloor0p25_sep0 | trans_entropy_with_floor | 0.25 | 1 | 0 | 0 | 0 | 0 | 0.472254 | 0.57785 | 0.321128 | 0.005842 | 0.528834 | 0.789651 | 0.584746361 | 0.421805352 | 1.09124362 | 3.65861225 | 0.003061615 | 0.00420972063 |
| p9016_pcgamma1_td1_transfloor0p75_sep0 | trans_entropy_with_floor | 0.75 | 1 | 0 | 0 | 0 | 0 | 0.465352 | 0.569386 | 0.31646 | 0.001174 | 0.518395 | 0.792164 | 0.589335799 | 0.425115913 | 1.6642921 | 3.58475852 | 0.00114810563 | 0.00191350938 |
| p9016_pcgamma1_td1_baseline_sep0_lct0_gct0 | constant | 0 | 1 | 0 | 0 | 0 | 0 | 0.464675 | 0.569056 | 0.315286 | 0 | 0.511735 | 0.800914 | 0.607293606 | 0.438069761 | 1.52243781 | 3.72219133 | 0.0015308075 | 0.003061615 |
| p9016_pcgamma1_td0p5_transfloor0p75_sep0 | trans_entropy_with_floor | 0.75 | 0.5 | 0 | 0 | 0 | 0 | 0.468323 | 0.575836 | 0.314453 | -0.000833 | 0.516506 | 0.79455 | 0.64328146 | 0.464029461 | 0.655248463 | 3.05205321 | 0.00191350938 | 0.0045924225 |
| p9016_pcgamma1_td1_transentropy_sep0 | trans_entropy | 0 | 1 | 0 | 0 | 0 | 0 | 0.464586 | 0.57042 | 0.313119 | -0.002167 | 0.520312 | 0.769132 | 0.587744892 | 0.423968285 | 0.918976307 | 3.54183364 | 0.00420972063 | 0.00535782625 |

## Interpretation Boundary

A real success requires full-denominator trans top1 accuracy to increase by at least 0.1 over the posterior-count baseline without major cis/Spearman collapse. If high-confidence trans downweighting improves only pmax-filtered calls but not all-trans top1, the raw-pairs full-denominator trans identity problem remains unresolved.
