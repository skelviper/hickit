# 044 P9016 Multiseed Basin Audit 1Mb

This controlled blind-training experiment tests whether the remaining trans gap is merely a seed/basin issue. It does not add a new model mechanism. Training uses only the approved P9016 raw pairs file; phase labels and CHARM/3DG are eval-only.

- full result root: `/mnt/ssd/zliu/phase3/test_res/044-20260620_040939-p9016_multiseed_basin_audit_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/044-20260620_040939-p9016_multiseed_basin_audit_1m`
- headline: `NO_PLUS_0P1`
- success criterion: full-denominator trans top1 accuracy must improve by at least 0.1 over the td1 seed17/noise0 baseline
- tested families: td1 anchor, best035 anchor/multiseed, and two chr-pair M-step semantics whose 038 eval-only pair oracle had +0.1 ceiling
- copy0/copy1 are gauge labels; no maternal/paternal labels are used in training

## Main Results

| config_name | basin_family | seed | init_noise_short | trans_chr_pair_mstep_lambda | trans_chr_pair_mstep_power | trans_chr_pair_mstep_warmup_iter | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | delta_trans_vs_td1_seed17 | delta_trans_vs_best035_seed17 | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | copytrack_frac_cos_lt_0 | target_plus_0p1_vs_td1_met |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_multiseed_best035_seed71_noise0p05 | best035 | 71 | 0.05 | 0 | 1 | 0 | 0.483482 | 0.572808 | 0.35564 | 0.038131 | 0.002757 | 0.561305 | 0.748036 | 0.692817748 | 0.499762356 | 0.370847613 | 2.58046293 | 0.00344431688 | 0 |
| p9016_multiseed_best035_seed17_noise0 | best035 | 17 | 0 | 0 | 1 | 0 | 0.481905 | 0.572055 | 0.352883 | 0.035374 | 0 | 0.55263 | 0.769516 | 0.696735561 | 0.502588511 | 0.458484918 | 2.83278489 | 0.0015308075 | 0 |
| p9016_multiseed_chrpair_m0p5_p4_w10_seed71_noise0p05 | chrpair_m0p5_p4_w10 | 71 | 0.05 | 0.5 | 4 | 10 | 0.469189 | 0.556463 | 0.344284 | 0.026775 | -0.008599 | 0.554616 | 0.796257 | 0.641890824 | 0.463026345 | 1.97402585 | 3.0470264 | 0.0015308075 | 0 |
| p9016_multiseed_best035_seed31_noise0p05 | best035 | 31 | 0.05 | 0 | 1 | 0 | 0.470286 | 0.559539 | 0.342548 | 0.025039 | -0.010335 | 0.54299 | 0.725411 | 0.697113812 | 0.502861321 | 0.237134248 | 2.35090852 | 0.00535782625 | 0 |
| p9016_multiseed_best035_seed17_noise0p05 | best035 | 17 | 0.05 | 0 | 1 | 0 | 0.476177 | 0.570818 | 0.340728 | 0.023219 | -0.012155 | 0.53969 | 0.752456 | 0.710712731 | 0.512670875 | 0.375790328 | 2.69662166 | 0.00191350938 | 0 |
| p9016_multiseed_chrpair_m0p5_p4_w10_seed31_noise0p05 | chrpair_m0p5_p4_w10 | 31 | 0.05 | 0.5 | 4 | 10 | 0.480142 | 0.578078 | 0.339978 | 0.022469 | -0.012905 | 0.544393 | 0.809394 | 0.59278506 | 0.427604049 | 1.59609532 | 3.03872943 | 0 | 0 |
| p9016_multiseed_chrpair_m1_p1_w0_seed47_noise0p05 | chrpair_m1_p1_w0 | 47 | 0.05 | 1 | 1 | 0 | 0.473843 | 0.567537 | 0.339749 | 0.02224 | -0.013134 | 0.540878 | 0.801155 | 0.592511952 | 0.427407026 | 0.90223515 | 3.12432122 | 0.00114810563 | 0 |
| p9016_multiseed_chrpair_m1_p1_w0_seed31_noise0p05 | chrpair_m1_p1_w0 | 31 | 0.05 | 1 | 1 | 0 | 0.483133 | 0.584712 | 0.337755 | 0.020246 | -0.015128 | 0.54092 | 0.819437 | 0.642468393 | 0.463442981 | 0.83201611 | 2.85311103 | 0.00076540375 | 0 |
| p9016_multiseed_chrpair_m1_p1_w0_seed71_noise0p05 | chrpair_m1_p1_w0 | 71 | 0.05 | 1 | 1 | 0 | 0.477819 | 0.575938 | 0.337394 | 0.019885 | -0.015489 | 0.541309 | 0.826292 | 0.607097626 | 0.437928379 | 1.9653976 | 3.34487486 | 0 | 0 |
| p9016_multiseed_chrpair_m0p5_p4_w10_seed17_noise0p05 | chrpair_m0p5_p4_w10 | 17 | 0.05 | 0.5 | 4 | 10 | 0.476319 | 0.575224 | 0.334769 | 0.01726 | -0.018114 | 0.529154 | 0.797472 | 0.628831267 | 0.45360589 | 1.93848681 | 3.01703 | 0 | 0 |
| p9016_multiseed_chrpair_m1_p1_w0_seed83_noise0p05 | chrpair_m1_p1_w0 | 83 | 0.05 | 1 | 1 | 0 | 0.480068 | 0.582732 | 0.333137 | 0.015628 | -0.019746 | 0.550692 | 0.809694 | 0.645166755 | 0.465389431 | 1.52334809 | 2.88464952 | 0 | 0 |
| p9016_multiseed_chrpair_m0p5_p4_w10_seed83_noise0p05 | chrpair_m0p5_p4_w10 | 83 | 0.05 | 0.5 | 4 | 10 | 0.464761 | 0.561112 | 0.326865 | 0.009356 | -0.026018 | 0.540788 | 0.792351 | 0.661380768 | 0.477085382 | 1.54491055 | 2.95907044 | 0.00076540375 | 0 |
| p9016_multiseed_chrpair_m1_p1_w0_seed17_noise0p05 | chrpair_m1_p1_w0 | 17 | 0.05 | 1 | 1 | 0 | 0.47404 | 0.577796 | 0.325545 | 0.008036 | -0.027338 | 0.525695 | 0.823713 | 0.626831949 | 0.452163666 | 1.25513124 | 2.9374547 | 0 | 0 |
| p9016_multiseed_chrpair_m0p5_p4_w10_seed59_noise0p05 | chrpair_m0p5_p4_w10 | 59 | 0.05 | 0.5 | 4 | 10 | 0.47228 | 0.575433 | 0.324649 | 0.00714 | -0.028234 | 0.533794 | 0.785623 | 0.667563856 | 0.481545508 | 1.36636722 | 2.96427655 | 0.00076540375 | 0 |
| p9016_multiseed_chrpair_m0p5_p4_w10_seed23_noise0p05 | chrpair_m0p5_p4_w10 | 23 | 0.05 | 0.5 | 4 | 10 | 0.459361 | 0.556768 | 0.319954 | 0.002445 | -0.032929 | 0.538919 | 0.775214 | 0.602353454 | 0.434506178 | 0.591273129 | 3.07136846 | 0.000382701875 | 0 |
| p9016_multiseed_chrpair_m1_p1_w0_seed23_noise0p05 | chrpair_m1_p1_w0 | 23 | 0.05 | 1 | 1 | 0 | 0.469249 | 0.575234 | 0.317565 | 5.6e-05 | -0.035318 | 0.53462 | 0.797971 | 0.607577801 | 0.438274741 | 1.10074508 | 3.24164462 | 0 | 0 |
| p9016_multiseed_td1_seed17_noise0 | td1 | 17 | 0 | 0 | 1 | 0 | 0.466095 | 0.569915 | 0.317509 | 0 | -0.035374 | 0.514631 | 0.80181 | 0.6116696 | 0.441226333 | 1.51912522 | 3.66526151 | 0.0015308075 | 0 |
| p9016_multiseed_best035_seed83_noise0p05 | best035 | 83 | 0.05 | 0 | 1 | 0 | 0.453039 | 0.548635 | 0.316224 | -0.001285 | -0.036659 | 0.522993 | 0.708704 | 0.694479823 | 0.500961304 | 1.01912582 | 2.47752666 | 0.00229621125 | 0 |
| p9016_multiseed_best035_seed47_noise0p05 | best035 | 47 | 0.05 | 0 | 1 | 0 | 0.457258 | 0.55621 | 0.315641 | -0.001868 | -0.037242 | 0.514624 | 0.709321 | 0.685047984 | 0.494157672 | 0.422151715 | 2.41585326 | 0.00880214313 | 0 |
| p9016_multiseed_best035_seed97_noise0p05 | best035 | 97 | 0.05 | 0 | 1 | 0 | 0.436543 | 0.523981 | 0.311404 | -0.006105 | -0.041479 | 0.523903 | 0.630451 | 0.6872527 | 0.495748043 | 0.749864101 | 2.06341004 | 0.0045924225 | 0 |
| p9016_multiseed_best035_seed59_noise0p05 | best035 | 59 | 0.05 | 0 | 1 | 0 | 0.448796 | 0.544806 | 0.31139 | -0.006119 | -0.041493 | 0.528939 | 0.704147 | 0.690476596 | 0.498073578 | 0.253713846 | 2.51502442 | 0.00612323 | 0 |
| p9016_multiseed_chrpair_m0p5_p4_w10_seed47_noise0p05 | chrpair_m0p5_p4_w10 | 47 | 0.05 | 0.5 | 4 | 10 | 0.450748 | 0.550823 | 0.307521 | -0.009988 | -0.045362 | 0.509984 | 0.76474 | 0.603851438 | 0.435586751 | 1.25704026 | 2.91185498 | 0.00267891313 | 0 |
| p9016_multiseed_best035_seed23_noise0p05 | best035 | 23 | 0.05 | 0 | 1 | 0 | 0.449842 | 0.552153 | 0.303417 | -0.014092 | -0.049466 | 0.526167 | 0.710717 | 0.668001294 | 0.481861055 | 0.295423239 | 2.63665986 | 0.00076540375 | 0 |
| p9016_multiseed_chrpair_m0p5_p4_w10_seed97_noise0p05 | chrpair_m0p5_p4_w10 | 97 | 0.05 | 0.5 | 4 | 10 | 0.438817 | 0.535095 | 0.301027 | -0.016482 | -0.051856 | 0.501594 | 0.756854 | 0.619956613 | 0.447204143 | 0.875445366 | 2.81355071 | 0.00344431688 | 0 |
| p9016_multiseed_chrpair_m1_p1_w0_seed97_noise0p05 | chrpair_m1_p1_w0 | 97 | 0.05 | 1 | 1 | 0 | 0.451279 | 0.557283 | 0.299569 | -0.01794 | -0.053314 | 0.506859 | 0.803181 | 0.547770858 | 0.395133168 | 1.09148681 | 2.92598677 | 0.00076540375 | 0 |
| p9016_multiseed_chrpair_m1_p1_w0_seed59_noise0p05 | chrpair_m1_p1_w0 | 59 | 0.05 | 1 | 1 | 0 | 0.452233 | 0.561932 | 0.295235 | -0.022274 | -0.057648 | 0.512075 | 0.768932 | 0.600555122 | 0.433208942 | 1.05137479 | 2.95582581 | 0.00076540375 | 0 |

## Interpretation Boundary

If no seed reaches the +0.1 full-denominator trans target, this argues against a simple basin-selection explanation under the current raw-pairs-only 1Mb training boundary. A positive result here would still need a blind selection rule, because eval labels are used only after training.

## Blind Selection Audits

Two post-training selection audits were added after the main run:

- `blind_selection_audit.tsv` ranks configs by blind-visible run diagnostics from `summary.tsv`, including loop diagnostics, force ratios, separation, and copytrack geometry.
- `raw_posterior_summary.tsv` and `raw_posterior_selection_audit.tsv` summarize each final `p9016_full.bpair_posterior.tsv` using only raw posterior/contact fields, then ask whether those final-posterior summaries can select a high-trans basin.

These are not heldout-validation results. They are post-training audits of whether training-visible or final raw-posterior quantities could have selected the better basin before looking at SNP/CHARM evaluation. The selection metrics do not use SNP truth, phase labels, CHARM/3DG, or reference structures. Columns such as `selected_trans`, `eval_best_trans`, `selected_rank_by_trans`, and `spearman_with_eval_trans` are eval-side annotations used only to score what a blind/raw-posterior selector would have picked.

Both audits are negative:

| audit | selector rows | selected eval-best rows | rows meeting +0.1 target | best selected trans |
|---|---:|---:|---:|---:|
| `blind_selection_audit.tsv` | 52 | 0 | 0 | 0.352883 |
| `raw_posterior_selection_audit.tsv` | 108 | 0 | 0 | 0.352883 |

The eval-best basin is `p9016_multiseed_best035_seed71_noise0p05` with trans top1 `0.355640`, but neither audit selects it. The best blind/raw-posterior selected trans value remains `0.352883`, below the +0.1 target `0.417509`. This strengthens the conclusion that the observed multiseed variation is not presently blind-selectable under the P9016 raw-pairs-only 1Mb boundary.
