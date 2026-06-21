# 030 P9016 Global Copytrack 1Mb

This controlled blind-training experiment tests whether local and chromosome-mean homolog-vector gauge alignment can improve P9016 trans four-state accuracy.

- full result root: `/mnt/ssd/zliu/phase3/test_res/030-20260619_081645-p9016_global_copytrack_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/030-20260619_081645-p9016_global_copytrack_1m`
- summary: `/mnt/ssd/zliu/phase3/test_res/030-20260619_081645-p9016_global_copytrack_1m/summary.tsv`
- trans delta summary: `/mnt/ssd/zliu/phase3/test_res/030-20260619_081645-p9016_global_copytrack_1m/trans_delta_summary.tsv`
- headline: `NO_PLUS_0P1`
- training input: approved raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper
- backend: CPU for this experiment, because GPU extra-force semantics use operator splitting and are not used for conclusions here

## Main Results

| config_name | lambda_copytrack | lambda_global_copytrack | min_sep_unit | lambda_sep | model_top1_accuracy_genome_trans | delta_trans_vs_baseline | model_top1_accuracy_genome_cis | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | final_mean_sep | copytrack_frac_cos_lt_0 | global_copytrack_force_l1_over_contact_force_l1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_msep1_lsep1_lct0_gct0p01_ieps0p5_noise0_seed17 | 0 | 0.00999999978 | 1 | 1 | 0.33738 | 0.02201 | 0.571827 | 0.79352 | 0.600699723 | 0.43331328 | 1.57954323 | 3.38730145 | 4.83461523 | 0.00076540375 | 0.00788051522 |
| p9016_pcgamma1_sep0_lct0p01_gct0p003_ieps0p5_noise0_seed17 | 0.00999999978 | 0.00300000003 | 0 | 0 | 0.332838 | 0.017468 | 0.56761 | 0.793358 | 0.610727489 | 0.440546781 | 1.51625073 | 3.56462717 | 4.97009754 | 0.0015308075 | 0.00250406203 |
| p9016_pcgamma1_sep0_lct0p1_gct0_ieps0p5_noise0_seed17 | 0.100000001 | 0 | 0 | 0 | 0.332303 | 0.016933 | 0.564659 | 0.790981 | 0.612564385 | 0.441871822 | 1.65871501 | 3.51053071 | 4.97388268 | 0.0015308075 | 0 |
| p9016_pcgamma1_sep0_lct0p01_gct0_ieps0p5_noise0_seed17 | 0.00999999978 | 0 | 0 | 0 | 0.331838 | 0.016468 | 0.56366 | 0.79502 | 0.610551178 | 0.440419555 | 1.39938819 | 3.56279898 | 5.02834511 | 0.00191350938 | 0 |
| p9016_pcgamma1_sep0_lct0p003_gct0p01_ieps0p5_noise0_seed17 | 0.00300000003 | 0.00999999978 | 0 | 0 | 0.331581 | 0.016211 | 0.570473 | 0.794388 | 0.617691517 | 0.44557023 | 1.31515574 | 3.48987579 | 4.88200665 | 0.0015308075 | 0.00794899385 |
| p9016_pcgamma1_sep0_lct0_gct0p01_ieps0p5_noise0_seed17 | 0 | 0.00999999978 | 0 | 0 | 0.331393 | 0.016023 | 0.57057 | 0.794273 | 0.618445218 | 0.446113944 | 1.37579799 | 3.49611449 | 4.87934637 | 0.0015308075 | 0.00795340053 |
| p9016_pcgamma1_sep0_lct0p003_gct0p003_ieps0p5_noise0_seed17 | 0.00300000003 | 0.00300000003 | 0 | 0 | 0.331275 | 0.015905 | 0.569051 | 0.799291 | 0.609913051 | 0.439959258 | 1.61267376 | 3.65783191 | 4.94982815 | 0.0015308075 | 0.0025034852 |
| p9016_pcgamma1_sep0_lct0p01_gct0p01_ieps0p5_noise0_seed17 | 0.00999999978 | 0.00999999978 | 0 | 0 | 0.330122 | 0.014752 | 0.571337 | 0.794973 | 0.618248343 | 0.445971906 | 1.55240119 | 3.48404813 | 4.87223482 | 0.0015308075 | 0.00796116052 |
| p9016_pcgamma1_sep0_lct0p001_gct0_ieps0p5_noise0_seed17 | 0.00100000005 | 0 | 0 | 0 | 0.330032 | 0.014662 | 0.566678 | 0.79793 | 0.606311202 | 0.437361062 | 1.61915672 | 3.68939471 | 5.02152395 | 0.0015308075 | 0 |
| p9016_pcgamma1_sep0_lct0p03_gct0_ieps0p5_noise0_seed17 | 0.0299999993 | 0 | 0 | 0 | 0.330032 | 0.014662 | 0.565033 | 0.795811 | 0.610972822 | 0.440723717 | 1.38533425 | 3.59001327 | 5.00980568 | 0.0015308075 | 0 |
| p9016_pcgamma1_msep1_lsep0p5_lct0p01_gct0p01_ieps0p5_noise0_seed17 | 0.00999999978 | 0.00999999978 | 1 | 0.5 | 0.328004 | 0.012634 | 0.574724 | 0.793963 | 0.605196118 | 0.436556727 | 1.31942332 | 3.48908162 | 4.84167242 | 0.0015308075 | 0.00796258192 |
| p9016_pcgamma1_msep1_lsep0p5_lct0_gct0p01_ieps0p5_noise0_seed17 | 0 | 0.00999999978 | 1 | 0.5 | 0.327434 | 0.012064 | 0.571628 | 0.791643 | 0.605228722 | 0.436580271 | 1.53278816 | 3.51447892 | 4.85355091 | 0.0015308075 | 0.00797217484 |
| p9016_pcgamma1_sep0_lct0p003_gct0_ieps0p5_noise0_seed17 | 0.00300000003 | 0 | 0 | 0 | 0.327122 | 0.011752 | 0.567862 | 0.798299 | 0.607007027 | 0.437863022 | 1.52308083 | 3.68156242 | 5.01408863 | 0.0015308075 | 0 |
| p9016_pcgamma1_msep1_lsep1_lct0p01_gct0p01_ieps0p5_noise0_seed17 | 0.00999999978 | 0.00999999978 | 1 | 1 | 0.324781 | 0.009411 | 0.575147 | 0.792007 | 0.601511538 | 0.433898866 | 1.60109282 | 3.39084244 | 4.83074951 | 0.00076540375 | 0.00787159719 |
| p9016_pcgamma1_sep0_lct0_gct0p003_ieps0p5_noise0_seed17 | 0 | 0.00300000003 | 0 | 0 | 0.323795 | 0.008425 | 0.568687 | 0.797256 | 0.610846996 | 0.440632969 | 1.42230642 | 3.57992506 | 4.97040129 | 0.00191350938 | 0.00250506236 |
| p9016_pcgamma1_sep0_lct0p03_gct0p01_ieps0p5_noise0_seed17 | 0.0299999993 | 0.00999999978 | 0 | 0 | 0.323461 | 0.008091 | 0.571225 | 0.791962 | 0.620948434 | 0.447919607 | 1.39030683 | 3.43103409 | 4.87455988 | 0.00191350938 | 0.00805960709 |
| p9016_pcgamma1_msep1_lsep0p5_lct0p03_gct0p01_ieps0p5_noise0_seed17 | 0.0299999993 | 0.00999999978 | 1 | 0.5 | 0.321551 | 0.006181 | 0.572196 | 0.791861 | 0.607792139 | 0.438429356 | 1.42068315 | 3.32314348 | 4.80211735 | 0.00114810563 | 0.00770996681 |
| p9016_pcgamma1_msep1_lsep1_lct0p03_gct0p01_ieps0p5_noise0_seed17 | 0.0299999993 | 0.00999999978 | 1 | 1 | 0.320489 | 0.005119 | 0.571007 | 0.788995 | 0.605444014 | 0.436735511 | 1.71426082 | 3.23360848 | 4.796422 | 0.00076540375 | 0.00772642471 |
| p9016_pcgamma1_sep0_lct0_gct0p03_ieps0p5_noise0_seed17 | 0 | 0.0299999993 | 0 | 0 | 0.320114 | 0.004744 | 0.559748 | 0.781911 | 0.636047959 | 0.458811611 | 1.25102425 | 3.33570457 | 4.71781111 | 0.0015308075 | 0.0223580739 |
| p9016_pcgamma1_sep0_lct0p01_gct0p03_ieps0p5_noise0_seed17 | 0.00999999978 | 0.0299999993 | 0 | 0 | 0.319093 | 0.003723 | 0.563626 | 0.779025 | 0.637080669 | 0.45955655 | 1.14614439 | 3.26086378 | 4.7031579 | 0.00344431688 | 0.0222840917 |
| p9016_pcgamma1_msep1p5_lsep1_lct0p01_gct0p01_ieps0p5_noise0_seed17 | 0.00999999978 | 0.00999999978 | 1.5 | 1 | 0.315668 | 0.000298 | 0.571832 | 0.797739 | 0.604443133 | 0.436013579 | 1.45236945 | 3.34523559 | 4.8393507 | 0.00076540375 | 0.00792100587 |
| p9016_pcgamma1_sep0_lct0_gct0_ieps0p5_noise0_seed17 | 0 | 0 | 0 | 0 | 0.31537 | 0 | 0.569634 | 0.800612 | 0.607306719 | 0.438079208 | 1.53509235 | 3.72130752 | 5.014678 | 0.00191350938 | 0 |
| p9016_pcgamma1_sep0_lct0_gct0p001_ieps0p5_noise0_seed17 | 0 | 0.00100000005 | 0 | 0 | 0.31389 | -0.00148 | 0.567324 | 0.797198 | 0.611304224 | 0.440962762 | 1.56193113 | 3.62275386 | 4.99303818 | 0.00191350938 | 0.000843288735 |
| p9016_pcgamma1_sep0_lct0_gct0p1_ieps0p5_noise0_seed17 | 0 | 0.100000001 | 0 | 0 | 0.305799 | -0.009571 | 0.537376 | 0.6899 | 0.634370327 | 0.457601458 | 1.21122575 | 3.06141543 | 4.64760113 | 0.00535782625 | 0.0442287684 |

## Interpretation Boundary

Global copytrack is an internal gauge-alignment force on chromosome-mean copy1-minus-copy0 vectors. It does not assign maternal/paternal identity and does not use phase labels or CHARM/3DG during training. A real success requires trans top1 to increase by at least 0.1 without a major cis/Spearman collapse.
