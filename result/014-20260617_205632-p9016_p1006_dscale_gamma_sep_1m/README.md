# 014-20260617_205632-p9016_p1006_dscale_gamma_sep_1m

Posterior-count d_scale gamma sweep across P9016 and P1006, plus modest homolog-separation guardrail tests at 1 Mb.

## Paths

- full result root: `/mnt/ssd/zliu/phase3/test_res/014-20260617_205632-p9016_p1006_dscale_gamma_sep_1m`
- lightweight hickit result: `/mnt/ssd/zliu/phase3/hickit/result/014-20260617_205632-p9016_p1006_dscale_gamma_sep_1m`
- summary.tsv: `/mnt/ssd/zliu/phase3/test_res/014-20260617_205632-p9016_p1006_dscale_gamma_sep_1m/summary.tsv`
- gamma_delta.tsv: `/mnt/ssd/zliu/phase3/test_res/014-20260617_205632-p9016_p1006_dscale_gamma_sep_1m/gamma_delta.tsv`
- sep_delta.tsv: `/mnt/ssd/zliu/phase3/test_res/014-20260617_205632-p9016_p1006_dscale_gamma_sep_1m/sep_delta.tsv`
- sample_delta.tsv: `/mnt/ssd/zliu/phase3/test_res/014-20260617_205632-p9016_p1006_dscale_gamma_sep_1m/sample_delta.tsv`
- commands.log: `/mnt/ssd/zliu/phase3/test_res/014-20260617_205632-p9016_p1006_dscale_gamma_sep_1m/commands.log`

## Training Boundary

Training uses only raw sample pairs. Phase labels and CHARM/3DG are eval-only through `scripts/p9016_common_eval.sh`. Manifests should report `uses_phase_labels=0` and `uses_charm_or_reference=0`.

## Git And Build

- git commit: `0ec1ebad10b1d7eea082841b6bfbd93135d9de49`
- git dirty count at start: `10`
- WARNING: dirty tree at launch; inspect `logs/git_status.txt` and `logs/git_diff.patch` in full test_res.
- run binary sha256: `b02d8094a49f72567b685870a6496be4c78f49cb012ba8fb9e73ca5b3887c95f`
- backend: `gpu`

## Configs

| sample | config | gamma | min_sep_unit | lambda_sep | backend |
| --- | --- | ---: | ---: | ---: | --- |
| P1006 | `p1006_pcgamma0_sep_off_ieps0p5_noise0_seed17` | 0 | 0 | 0 | gpu |
| P1006 | `p1006_pcgamma0p25_sep_off_ieps0p5_noise0_seed17` | 0.25 | 0 | 0 | gpu |
| P1006 | `p1006_pcgamma0p5_msep1_lsep0p5_ieps0p5_noise0_seed17` | 0.5 | 1 | 0.5 | gpu |
| P1006 | `p1006_pcgamma0p5_msep1_lsep1_ieps0p5_noise0_seed17` | 0.5 | 1 | 1 | gpu |
| P1006 | `p1006_pcgamma0p5_sep_off_ieps0p5_noise0_seed17` | 0.5 | 0 | 0 | gpu |
| P1006 | `p1006_pcgamma0p75_sep_off_ieps0p5_noise0_seed17` | 0.75 | 0 | 0 | gpu |
| P1006 | `p1006_pcgamma1_msep1_lsep0p5_ieps0p5_noise0_seed17` | 1 | 1 | 0.5 | gpu |
| P1006 | `p1006_pcgamma1_msep1_lsep1_ieps0p5_noise0_seed17` | 1 | 1 | 1 | gpu |
| P1006 | `p1006_pcgamma1_msep1p5_lsep1_ieps0p5_noise0_seed17` | 1 | 1.5 | 1 | gpu |
| P1006 | `p1006_pcgamma1_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | gpu |
| P9016 | `p9016_pcgamma0_sep_off_ieps0p5_noise0_seed17` | 0 | 0 | 0 | gpu |
| P9016 | `p9016_pcgamma0p25_sep_off_ieps0p5_noise0_seed17` | 0.25 | 0 | 0 | gpu |
| P9016 | `p9016_pcgamma0p5_msep1_lsep0p5_ieps0p5_noise0_seed17` | 0.5 | 1 | 0.5 | gpu |
| P9016 | `p9016_pcgamma0p5_msep1_lsep1_ieps0p5_noise0_seed17` | 0.5 | 1 | 1 | gpu |
| P9016 | `p9016_pcgamma0p5_sep_off_ieps0p5_noise0_seed17` | 0.5 | 0 | 0 | gpu |
| P9016 | `p9016_pcgamma0p75_sep_off_ieps0p5_noise0_seed17` | 0.75 | 0 | 0 | gpu |
| P9016 | `p9016_pcgamma1_msep1_lsep0p5_ieps0p5_noise0_seed17` | 1 | 1 | 0.5 | gpu |
| P9016 | `p9016_pcgamma1_msep1_lsep1_ieps0p5_noise0_seed17` | 1 | 1 | 1 | gpu |
| P9016 | `p9016_pcgamma1_msep1p5_lsep1_ieps0p5_noise0_seed17` | 1 | 1.5 | 1 | gpu |
| P9016 | `p9016_pcgamma1_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | gpu |

## Main Results

| sample | config | gamma | min_sep | lambda_sep | top1 all | top1 cis | top1 trans | same/cross cis | cis Spearman | entropy | pU | min sep | sep_p05 | mean sep | cos<0 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| P1006 | `p1006_pcgamma0_sep_off_ieps0p5_noise0_seed17` | 0 | 0 | 0 | 0.437179 | 0.498543 | 0.365062 | 0.875782 | 0.561463 | 0.750889957 | 0.54165262 | 0.0725905448 | 0.750393569 | 3.75432897 | 0.0164561806 |
| P1006 | `p1006_pcgamma0p25_sep_off_ieps0p5_noise0_seed17` | 0.25 | 0 | 0 | 0.438857 | 0.503976 | 0.362327 | 0.898671 | 0.586205 | 0.723253071 | 0.521716833 | 0.251400173 | 1.36867046 | 4.02481794 | 0.0141599694 |
| P1006 | `p1006_pcgamma0p5_msep1_lsep0p5_ieps0p5_noise0_seed17` | 0.5 | 1 | 0.5 | 0.460673 | 0.519754 | 0.391239 | 0.910969 | 0.614501 | 0.683582127 | 0.493100286 | 0.560481191 | 1.97284436 | 4.35219145 | 0.0168388825 |
| P1006 | `p1006_pcgamma0p5_msep1_lsep1_ieps0p5_noise0_seed17` | 0.5 | 1 | 1 | 0.459843 | 0.518706 | 0.390666 | 0.916211 | 0.60475 | 0.675743222 | 0.487445682 | 0.627615035 | 2.0245564 | 4.3936553 | 0.0156907769 |
| P1006 | `p1006_pcgamma0p5_sep_off_ieps0p5_noise0_seed17` | 0.5 | 0 | 0 | 0.449939 | 0.512146 | 0.376832 | 0.900157 | 0.582439 | 0.69119221 | 0.498589754 | 0.195361316 | 1.90029037 | 4.25061798 | 0.0164561806 |
| P1006 | `p1006_pcgamma0p75_sep_off_ieps0p5_noise0_seed17` | 0.75 | 0 | 0 | 0.448798 | 0.519818 | 0.365334 | 0.908706 | 0.597272 | 0.66863203 | 0.482316047 | 0.68348968 | 2.11805868 | 4.21392727 | 0.0179869881 |
| P1006 | `p1006_pcgamma1_msep1_lsep0p5_ieps0p5_noise0_seed17` | 1 | 1 | 0.5 | 0.445512 | 0.518317 | 0.35995 | 0.926191 | 0.633733 | 0.64821738 | 0.467589974 | 0.941831589 | 2.57705212 | 4.50541019 | 0.0126291619 |
| P1006 | `p1006_pcgamma1_msep1_lsep1_ieps0p5_noise0_seed17` | 1 | 1 | 1 | 0.446884 | 0.521905 | 0.358718 | 0.930069 | 0.643198 | 0.647177219 | 0.466839701 | 0.965927005 | 2.59566569 | 4.50228691 | 0.0107156525 |
| P1006 | `p1006_pcgamma1_msep1p5_lsep1_ieps0p5_noise0_seed17` | 1 | 1.5 | 1 | 0.450686 | 0.527318 | 0.360626 | 0.93779 | 0.658392 | 0.635517478 | 0.458428979 | 1.13204575 | 2.71032453 | 4.52591038 | 0.00574052813 |
| P1006 | `p1006_pcgamma1_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | 0.45229 | 0.516639 | 0.376664 | 0.921324 | 0.619827 | 0.652609468 | 0.4707582 | 0.654938996 | 2.4321084 | 4.47344398 | 0.0145426713 |
| P9016 | `p9016_pcgamma0_sep_off_ieps0p5_noise0_seed17` | 0 | 0 | 0 | 0.410629 | 0.496261 | 0.288074 | 0.87484 | 0.57727 | 0.804043949 | 0.579995096 | 0.0775339454 | 0.487044543 | 3.43938279 | 0.00956754688 |
| P9016 | `p9016_pcgamma0p25_sep_off_ieps0p5_noise0_seed17` | 0.25 | 0 | 0 | 0.42309 | 0.504749 | 0.306223 | 0.886405 | 0.607277 | 0.725981593 | 0.523685038 | 0.106479555 | 1.00029325 | 3.8365283 | 0.00841944126 |
| P9016 | `p9016_pcgamma0p5_msep1_lsep0p5_ieps0p5_noise0_seed17` | 0.5 | 1 | 0.5 | 0.446251 | 0.541675 | 0.309681 | 0.942467 | 0.751263 | 0.661850154 | 0.477423966 | 0.46203199 | 2.11252952 | 4.45830822 | 0.00650593188 |
| P9016 | `p9016_pcgamma0p5_msep1_lsep1_ieps0p5_noise0_seed17` | 0.5 | 1 | 1 | 0.45209 | 0.547737 | 0.315203 | 0.955163 | 0.77686 | 0.638747096 | 0.460758626 | 0.544398665 | 2.08213711 | 4.61278582 | 0.003061615 |
| P9016 | `p9016_pcgamma0p5_sep_off_ieps0p5_noise0_seed17` | 0.5 | 0 | 0 | 0.446245 | 0.536662 | 0.316842 | 0.934062 | 0.727488 | 0.644453108 | 0.464874685 | 0.239052653 | 2.15219688 | 4.55422401 | 0.003061615 |
| P9016 | `p9016_pcgamma0p75_sep_off_ieps0p5_noise0_seed17` | 0.75 | 0 | 0 | 0.457318 | 0.563034 | 0.306021 | 0.968989 | 0.797652 | 0.617134631 | 0.445168525 | 1.12510908 | 3.51288509 | 4.96309948 | 0.00076540375 |
| P9016 | `p9016_pcgamma1_msep1_lsep0p5_ieps0p5_noise0_seed17` | 1 | 1 | 0.5 | 0.472694 | 0.56745 | 0.337082 | 0.970722 | 0.798876 | 0.597418725 | 0.430946529 | 1.38624048 | 3.67749 | 4.96726274 | 0.0015308075 |
| P9016 | `p9016_pcgamma1_msep1_lsep1_ieps0p5_noise0_seed17` | 1 | 1 | 1 | 0.467915 | 0.567741 | 0.325045 | 0.969877 | 0.794501 | 0.596827805 | 0.430520266 | 1.35799825 | 3.58317113 | 4.93934774 | 0.00191350938 |
| P9016 | `p9016_pcgamma1_msep1p5_lsep1_ieps0p5_noise0_seed17` | 1 | 1.5 | 1 | 0.46956 | 0.572123 | 0.322774 | 0.969722 | 0.796762 | 0.597607255 | 0.431082517 | 1.31952798 | 3.53292322 | 4.92633867 | 0.00191350938 |
| P9016 | `p9016_pcgamma1_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | 0.464812 | 0.569275 | 0.315307 | 0.969591 | 0.80079 | 0.607262909 | 0.438047558 | 1.52262402 | 3.72325063 | 5.01577425 | 0.00191350938 |

## Interpretation

- Best P9016 gamma by top1_all: `1`.
- Best P1006 gamma by top1_all: `1`.
- Best P9016 gamma by cis Spearman: `1`.
- Best P1006 gamma by cis Spearman: `1`.
- Treat gamma=1 as the historical posterior-count baseline. If another gamma is not better by more than about 0.005 top1_all, prefer gamma=1 for continuity.
- Separation guardrails should be interpreted as geometry controls: accept them only if sep_p05/min separation improves without a meaningful top1 or Spearman cost.
- This run does not introduce priors, entropy-aware rho, annealing, sharpening, random init, copytrack smoothness, oracle training, 200 kb, or multiresolution changes.
