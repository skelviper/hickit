# 034 P9016 Entropy Rho Temperature 1Mb

This controlled blind-training experiment tests whether constant full-weight training of uncertain contacts reinforces collapsed/symmetric trans copy assignments.

- full result root: `/mnt/ssd/zliu/phase3/test_res/034-20260619_140108-p9016_entropy_rho_temperature_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/034-20260619_140108-p9016_entropy_rho_temperature_1m`
- headline: `NO_PLUS_0P1`
- training input: raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper
- changed knobs: rho_train_mode/rho_train_floor and temperature schedule; posterior_count gamma1 remains the baseline d_scale semantics
- unchanged model settings: uniform prior, no chromosome-pair prior, no copytrack force, no hard top1 M-step, no phase/CHARM training labels

## Main Results

| config_name | rho_train_mode | rho_train_floor | temperature_start | temperature_end | trans_dscale_multiplier | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_mean_rho_train_bpair | final_min_sep | sep_p05 | final_mean_sep | copytrack_frac_cos_lt_0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_rho_entropy_floor0p5_temp1_td0p5_sep0 | entropy_with_floor | 0.5 | 1 | 1 | 0.5 | 0.475859 | 0.567333 | 0.344944 | 0.538746 | 0.75923 | 0.700257838 | 0.505129218 | NA | 0.325333238 | 2.73246861 | 4.37077999 | 0.0015308075 |
| p9016_pcgamma1_rho_constant_temp1_td0p5_sep0 | constant | 0 | 1 | 1 | 0.5 | 0.475042 | 0.565979 | 0.344895 | 0.539371 | 0.758404 | 0.699464262 | 0.504556835 | NA | 0.296572924 | 2.72949338 | 4.37399006 | 0.0015308075 |
| p9016_pcgamma1_rho_entropy_floor0p25_temp1_td0p5_sep0 | entropy_with_floor | 0.25 | 1 | 1 | 0.5 | 0.474314 | 0.565149 | 0.344312 | 0.538683 | 0.758114 | 0.701603591 | 0.506099999 | NA | 0.310703486 | 2.7110374 | 4.36681652 | 0.00191350938 |
| p9016_pcgamma1_rho_entropy_floor0p5_temp2to1_td0p5_sep0 | entropy_with_floor | 0.5 | 2 | 1 | 0.5 | 0.450508 | 0.524797 | 0.344187 | 0.54283 | 0.55646 | 0.762722671 | 0.550188124 | NA | 0.237204537 | 1.6234889 | 3.86999917 | 0.0133945656 |
| p9016_pcgamma1_rho_entropy_floor0p25_temp2to1_td0p5_sep0 | entropy_with_floor | 0.25 | 2 | 1 | 0.5 | 0.449871 | 0.526257 | 0.340547 | 0.537558 | 0.566328 | 0.76397711 | 0.551092982 | NA | 0.272409618 | 1.60313642 | 3.86955619 | 0.01224646 |
| p9016_pcgamma1_rho_entropy_floor0p5_temp2to1_td1_sep0 | entropy_with_floor | 0.5 | 2 | 1 | 1 | 0.438015 | 0.518536 | 0.322774 | 0.520361 | 0.590908 | 0.661246836 | 0.476988763 | NA | 0.919389367 | 2.11179495 | 4.37761974 | 0.0130118638 |
| p9016_pcgamma1_rho_entropy_temp2to1_td1_sep0 | entropy | 0 | 2 | 1 | 1 | 0.438069 | 0.518667 | 0.322718 | 0.520278 | 0.590851 | 0.661292911 | 0.477022022 | NA | 0.916508615 | 2.10796475 | 4.37731791 | 0.0130118638 |
| p9016_pcgamma1_rho_entropy_cis_floor_trans0p5_temp1_td1_sep0 | entropy_cis_floor_trans | 0.5 | 1 | 1 | 1 | 0.468003 | 0.569541 | 0.322683 | 0.518729 | 0.800553 | 0.607251942 | 0.43803969 | NA | 1.49391007 | 3.72215891 | 5.0159502 | 0.0015308075 |
| p9016_pcgamma1_rho_entropy_floor0p5_temp1_td1_sep0 | entropy_with_floor | 0.5 | 1 | 1 | 1 | 0.467543 | 0.56878 | 0.322656 | 0.519076 | 0.800904 | 0.607259214 | 0.438044935 | NA | 1.50978088 | 3.71906877 | 5.01709318 | 0.00191350938 |
| p9016_pcgamma1_rho_entropy_cis_floor_trans0p25_temp1_td1_sep0 | entropy_cis_floor_trans | 0.25 | 1 | 1 | 1 | 0.467843 | 0.569328 | 0.3226 | 0.518708 | 0.800545 | 0.607269883 | 0.438052654 | NA | 1.4862839 | 3.72238922 | 5.01571178 | 0.0015308075 |
| p9016_pcgamma1_rho_entropy_floor0p25_temp2to1_td1_sep0 | entropy_with_floor | 0.25 | 2 | 1 | 1 | 0.437106 | 0.517473 | 0.322086 | 0.519694 | 0.59082 | 0.661805749 | 0.477391928 | NA | 0.913538158 | 2.09659505 | 4.38045073 | 0.0141599694 |
| p9016_pcgamma1_rho_constant_temp1_td1_sep0_ieps0p5_noise0_seed17 | constant | 0 | 1 | 1 | 1 | 0.464421 | 0.567178 | 0.317356 | 0.514353 | 0.80206 | 0.611631513 | 0.441198856 | NA | 1.53237367 | 3.65533304 | 5.0175581 | 0.00191350938 |
| p9016_pcgamma1_rho_entropy_floor0p25_temp1_td1_sep0 | entropy_with_floor | 0.25 | 1 | 1 | 1 | 0.464663 | 0.569037 | 0.315286 | 0.5117 | 0.800898 | 0.607277393 | 0.438058019 | NA | 1.52097571 | 3.72318769 | 5.01568031 | 0.00191350938 |
| p9016_pcgamma1_rho_entropy_temp1_td1_sep0 | entropy | 0 | 1 | 1 | 1 | 0.464729 | 0.569148 | 0.315286 | 0.511741 | 0.800824 | 0.60732168 | 0.438089967 | NA | 1.51003408 | 3.7230227 | 5.01561832 | 0.0015308075 |

## Interpretation Boundary

A real success requires full-denominator trans top1 accuracy to increase by at least 0.1 over the constant-rho baseline without a major cis accuracy or cis distance-Spearman collapse. If entropy-aware weighting improves rho/uncertainty but not trans identity, the next controlled mechanism is homolog-vector continuity rather than stronger trans contact force.
