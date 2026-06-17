# 012-20260617_111813-p9016_dscale_semantics_regression_1m

This controlled run tests whether historical 003/004 1 Mb scaffold softall behavior corresponds to posterior-count d_scale semantics.

## Paths

- full result root: `/mnt/ssd/zliu/phase3/test_res/012-20260617_111813-p9016_dscale_semantics_regression_1m`
- lightweight hickit result: `/mnt/ssd/zliu/phase3/hickit/result/012-20260617_111813-p9016_dscale_semantics_regression_1m`
- summary TSV: `/mnt/ssd/zliu/phase3/test_res/012-20260617_111813-p9016_dscale_semantics_regression_1m/summary.tsv`
- d_scale semantics audit TSV: `/mnt/ssd/zliu/phase3/test_res/012-20260617_111813-p9016_dscale_semantics_regression_1m/dscale_semantics_audit.tsv`
- d_scale semantics headline: `HISTORICAL_D_SCALE_SEMANTICS_CONFIRMED`

## Tree State

- hickit git dirty entries at start: `6`
- WARNING: this run started from a dirty tree; inspect `logs/git_status.txt` and `logs/git_diff.patch`.

## Training Boundary

Training reads only raw P9016 pairs. Phase labels and CHARM/3DG are used only by post-training eval.

## Configs

| config | backend | init_eps | init_noise | d_scale_mode | formula | min_sep_unit | lambda_sep |
| --- | --- | ---: | ---: | --- | --- | ---: | ---: |
| `current_raw_cpu_exact_ieps0p5_noise0_seed17` | cpu | 0.5 | 0 | raw_count | n_raw | 0 | 0 |
| `current_raw_gpu_exact_ieps0p5_noise0_seed17` | gpu | 0.5 | 0 | raw_count | n_raw | 0 | 0 |
| `posterior_count_common_msep1p5_lsep1_ieps1_noise0p05_seed17` | gpu | 1 | 0.0500000007 | posterior_count | n_raw*posterior_prob | 1.5 | 1 |
| `posterior_count_common_sep_off_ieps1_noise0p05_seed17` | gpu | 1 | 0.0500000007 | posterior_count | n_raw*posterior_prob | 0 | 0 |
| `posterior_count_cpu_exact_ieps0p5_noise0_seed17` | cpu | 0.5 | 0 | posterior_count | n_raw*posterior_prob | 0 | 0 |
| `posterior_count_gpu_exact_ieps0p5_noise0_seed17` | gpu | 0.5 | 0 | posterior_count | n_raw*posterior_prob | 0 | 0 |

## Main Results

| config | backend | top1 all | top1 cis | top1 trans | same/cross cis | cis Spearman | entropy | pU | min sep | mean sep | cos<0 | projection switch |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `current_raw_cpu_exact_ieps0p5_noise0_seed17` | cpu | 0.409357 | 0.494038 | 0.288164 | 0.871589 | 0.572039 | 0.80723244 | 0.58229512 | 0.0909356773 | 3.42471027 | 0.0103329506 | 0.0344431688 |
| `current_raw_gpu_exact_ieps0p5_noise0_seed17` | gpu | 0.416077 | 0.49531 | 0.30268 | 0.874345 | 0.573947 | 0.805009067 | 0.580691278 | 0.129438356 | 3.43346024 | 0.0103329506 | 0.0264064294 |
| `posterior_count_common_msep1p5_lsep1_ieps1_noise0p05_seed17` | gpu | 0.436823 | 0.54979 | 0.275148 | 0.974958 | 0.816021 | 0.618878543 | 0.446426481 | 2.94410539 | 5.5891428 | 0 | 0 |
| `posterior_count_common_sep_off_ieps1_noise0p05_seed17` | gpu | 0.437218 | 0.550382 | 0.275259 | 0.974958 | 0.815783 | 0.618496358 | 0.446150839 | 2.95507336 | 5.59178591 | 0 | 0 |
| `posterior_count_cpu_exact_ieps0p5_noise0_seed17` | cpu | 0.465049 | 0.569634 | 0.31537 | 0.969562 | 0.800612 | 0.607306719 | 0.438079208 | 1.53509235 | 5.014678 | 0.00191350938 | 0.003061615 |
| `posterior_count_gpu_exact_ieps0p5_noise0_seed17` | gpu | 0.464469 | 0.568022 | 0.316266 | 0.969882 | 0.801884 | 0.611757576 | 0.441289783 | 1.54554307 | 5.01241636 | 0.00191350938 | 0.003061615 |

## d_scale Semantic Audit

`HISTORICAL_D_SCALE_SEMANTICS_CONFIRMED`

Interpretation rules: if posterior_count exact matches 003/004 while raw_count exact fails, the 011 regression is explained by d_scale semantic drift. If both fail, d_scale semantics alone is insufficient and inputs/init/backend/compiler/dirty-tree differences should be checked. This run is a regression audit, not a model improvement claim.
