# 019-20260618_224538-p9016_copytrack_continuity_1m

Controlled blind P9016 copy-track continuity experiment at 1 Mb. The experiment tests whether a homolog-vector continuity force can reduce local copy-track flipping and improve trans phase accuracy.

## Paths

- full result root: `/mnt/ssd/zliu/phase3/test_res/019-20260618_224538-p9016_copytrack_continuity_1m`
- lightweight hickit result: `/mnt/ssd/zliu/phase3/hickit/result/019-20260618_224538-p9016_copytrack_continuity_1m`
- summary.tsv: `/mnt/ssd/zliu/phase3/test_res/019-20260618_224538-p9016_copytrack_continuity_1m/summary.tsv`
- trans_delta.tsv: `/mnt/ssd/zliu/phase3/test_res/019-20260618_224538-p9016_copytrack_continuity_1m/trans_delta.tsv`
- commands.log: `/mnt/ssd/zliu/phase3/test_res/019-20260618_224538-p9016_copytrack_continuity_1m/commands.log`

## Training Boundary

Training uses only P9016 pairs. SNP/phase labels and CHARM/3DG are used only in eval through `scripts/p9016_common_eval.sh`. The copy-track force is blind geometry: adjacent same-chromosome homolog vectors are encouraged to be continuous.

## Git And Build

- git commit: `1037ee5f79d54e943eef5c37e5aeb2d841698cb8`
- git dirty count at start: `31`
- WARNING: dirty tree at launch; inspect `logs/git_status.txt` and `logs/git_diff.patch`.
- run binary sha256: `9b23e86098d665ac796a44f5bce19502fec2e8e4c9f13f52577991831c0e8056`
- backend: `gpu`

## Main Results

| config | gamma | min_sep | lambda_sep | lambda_copytrack | top1 all | top1 cis | top1 trans | delta trans | cis Spearman | entropy | pU | sep_p05 | cos<0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `p9016_pcgamma0p5_copytrack0p003_sep_off_ieps0p5_noise0_seed17` | 0.5 | 0 | 0 | 0.00300000003 | 0.445545 | 0.53706 | 0.314571 | -0.002778 | 0.728268 | 0.644095719 | 0.464616835 | 2.12205195 | 0.00267891313 |
| `p9016_pcgamma0p5_copytrack0p01_sep_off_ieps0p5_noise0_seed17` | 0.5 | 0 | 0 | 0.00999999978 | 0.454664 | 0.542903 | 0.328379 | 0.01103 | 0.723658 | 0.654084742 | 0.471822411 | 2.28488183 | 0.0045924225 |
| `p9016_pcgamma0p5_copytrack0p03_sep_off_ieps0p5_noise0_seed17` | 0.5 | 0 | 0 | 0.0299999993 | 0.445034 | 0.536381 | 0.3143 | -0.003049 | 0.721601 | 0.647433579 | 0.467024595 | 2.07748103 | 0.00267891313 |
| `p9016_pcgamma0p5_msep1_lsep0p5_copytrack0p01_ieps0p5_noise0_seed17` | 0.5 | 1 | 0.5 | 0.00999999978 | 0.447717 | 0.543568 | 0.310536 | -0.006813 | 0.751398 | 0.663760662 | 0.478802085 | 2.22871518 | 0.00688863375 |
| `p9016_pcgamma0p5_msep1_lsep1_copytrack0p01_ieps0p5_noise0_seed17` | 0.5 | 1 | 1 | 0.00999999978 | 0.451222 | 0.54565 | 0.316078 | -0.001271 | 0.7772 | 0.636720955 | 0.459297091 | 2.39688277 | 0.00267891313 |
| `p9016_pcgamma1_copytrack0_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | 0 | 0.464258 | 0.566906 | 0.317349 | 0 | 0.802046 | 0.611699879 | 0.441248208 | 3.65642905 | 0.00191350938 |
| `p9016_pcgamma1_copytrack0p001_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | 0.00100000005 | 0.463115 | 0.568416 | 0.312411 | -0.004938 | 0.798059 | 0.60670501 | 0.437645137 | 3.70150208 | 0.0015308075 |
| `p9016_pcgamma1_copytrack0p003_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | 0.00300000003 | 0.464355 | 0.565538 | 0.319544 | 0.002195 | 0.796921 | 0.60938853 | 0.439580888 | 3.66589355 | 0.00191350938 |
| `p9016_pcgamma1_copytrack0p01_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | 0.00999999978 | 0.465726 | 0.564101 | 0.324934 | 0.007585 | 0.796511 | 0.61022681 | 0.440185577 | 3.62008357 | 0.00191350938 |
| `p9016_pcgamma1_copytrack0p03_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | 0.0299999993 | 0.4691 | 0.566314 | 0.329969 | 0.01262 | 0.79623 | 0.610959709 | 0.4407143 | 3.59635711 | 0.0015308075 |
| `p9016_pcgamma1_copytrack0p1_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | 0.100000001 | 0.469063 | 0.56877 | 0.326365 | 0.009016 | 0.789224 | 0.612320185 | 0.441695631 | 3.50605226 | 0.0015308075 |
| `p9016_pcgamma1_msep1_lsep0p5_copytrack0p01_ieps0p5_noise0_seed17` | 1 | 1 | 0.5 | 0.00999999978 | 0.472437 | 0.566363 | 0.338012 | 0.020663 | 0.800769 | 0.596784353 | 0.430488914 | 3.67958283 | 0.0015308075 |
| `p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_ieps0p5_noise0_seed17` | 1 | 1 | 0.5 | 0.0299999993 | 0.476474 | 0.571172 | 0.340943 | 0.023594 | 0.800932 | 0.595618963 | 0.42964825 | 3.68792653 | 0.0015308075 |
| `p9016_pcgamma1_msep1_lsep1_copytrack0p01_ieps0p5_noise0_seed17` | 1 | 1 | 1 | 0.00999999978 | 0.469317 | 0.571138 | 0.323593 | 0.006244 | 0.796419 | 0.59716326 | 0.430762231 | 3.61683774 | 0.00191350938 |
| `p9016_pcgamma1_msep1_lsep1_copytrack0p03_ieps0p5_noise0_seed17` | 1 | 1 | 1 | 0.0299999993 | 0.475628 | 0.570803 | 0.339415 | 0.022066 | 0.803514 | 0.598020494 | 0.4313806 | 3.68552732 | 0.00076540375 |
| `p9016_pcgamma1_msep1p5_lsep1_copytrack0p01_ieps0p5_noise0_seed17` | 1 | 1.5 | 1 | 0.00999999978 | 0.4715 | 0.57355 | 0.325448 | 0.008099 | 0.79733 | 0.596729457 | 0.430449337 | 3.53024435 | 0.0015308075 |

## Interpretation

- Best observed trans config: `p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_ieps0p5_noise0_seed17` with trans top1 `0.340943` and delta `0.023594` versus baseline.
- The +0.1 trans target was not met in this sweep; continue with chromosome-level gauge synchronization or E-step prior diagnostics.
- Accept copy-track continuity only if trans improves without a large cis/Spearman/entropy penalty. A drop in copytrack cos<0 without trans improvement means internal chromosome copy tracks stabilized but cross-chromosome gauge remains unresolved.
