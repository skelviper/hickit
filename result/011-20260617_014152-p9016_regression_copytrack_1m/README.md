# 011-20260617_014152-p9016_regression_copytrack_1m

This is a controlled P9016 blind-diploid infrastructure and regression audit. It does not introduce a new model.

## Paths

- full result root: `/mnt/ssd/zliu/phase3/test_res/011-20260617_014152-p9016_regression_copytrack_1m`
- lightweight hickit result: `/mnt/ssd/zliu/phase3/hickit/result/011-20260617_014152-p9016_regression_copytrack_1m`
- summary TSV: `/mnt/ssd/zliu/phase3/test_res/011-20260617_014152-p9016_regression_copytrack_1m/summary.tsv`
- regression audit TSV: `/mnt/ssd/zliu/phase3/test_res/011-20260617_014152-p9016_regression_copytrack_1m/regression_audit.tsv`
- copytrack summary TSV: `/mnt/ssd/zliu/phase3/test_res/011-20260617_014152-p9016_regression_copytrack_1m/copytrack_summary.tsv`

## Training Boundary

Training reads only the P9016 raw pairs file. Phase labels and CHARM/3DG are passed only to post-training eval. Copy labels are gauge labels, not parental labels.

## Configs

| config | backend | init_eps | init_noise_scale | min_sep_unit | lambda_sep | d_scale_mode |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| legacy_1m_cpu_exact_ieps0p5_noise0_seed17 | cpu | 0.5 | 0.0 | 0.0 | 0.0 | raw_count |
| legacy_1m_gpu_exact_ieps0p5_noise0_seed17 | gpu | 0.5 | 0.0 | 0.0 | 0.0 | raw_count |
| common_raw_sep_off_ieps1_noise0p05_seed17 | gpu | 1.0 | 0.05 | 0.0 | 0.0 | raw_count |
| common_raw_msep1p5_lsep1_ieps1_noise0p05_seed17 | gpu | 1.0 | 0.05 | 1.5 | 1.0 | raw_count |

## Regression Audit

REGRESSION_FAIL: at least one exact legacy baseline metric or manifest field differs from the historical 003/004-style baseline beyond tolerance. Debug this before interpreting later rows as model changes.

## Headline Metrics

| config | backend | top1 all | top1 cis | top1 trans | same/cross cis | mean cis Spearman | entropy | pU | min sep | mean sep | frac cos<0 | projection switch frac |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `legacy_1m_cpu_exact_ieps0p5_noise0_seed17` | cpu | 0.409357 | 0.494038 | 0.288164 | 0.871589 | 0.572039 | 0.80723244 | 0.58229512 | 0.0909356773 | 3.42471027 | 0.0103329506 | 0.0344431688 |
| `legacy_1m_gpu_exact_ieps0p5_noise0_seed17` | gpu | 0.409786 | 0.495188 | 0.28756 | 0.872317 | 0.571805 | 0.806435645 | 0.581720352 | 0.0815229639 | 3.42777348 | 0.0110983544 | 0.0359739763 |
| `common_raw_sep_off_ieps1_noise0p05_seed17` | gpu | 0.449325 | 0.549596 | 0.30582 | 0.972294 | 0.804125 | 0.656980932 | 0.473911583 | 0.160607338 | 5.22866869 | 0 | 0 |
| `common_raw_msep1p5_lsep1_ieps1_noise0p05_seed17` | gpu | 0.455113 | 0.54416 | 0.32767 | 0.973954 | 0.825303 | 0.656348825 | 0.473455578 | 1.0009135 | 5.24730158 | 0 | 0 |

## Interpretation

This run is a diagnostic and regression audit, not a final model. Stronger separation is only a guardrail: it should increase min/p05 separation, but it is not expected to solve 00-vs-11 identity alone. High `copytrack_frac_cos_lt_0` or high projection sign-switch rate supports local homolog-vector flipping as a bottleneck. If four-state top1 accuracy is low while same/cross cis is high, and copy-track discontinuity is high, the next controlled experiment should test a small homolog-vector continuity prior. If copy-track discontinuity is low, the next diagnostic should target oracle posterior or E-step prior behavior.
