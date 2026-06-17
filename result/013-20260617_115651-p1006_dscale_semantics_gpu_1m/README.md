# 013-20260617_115651-p1006_dscale_semantics_gpu_1m

P1006 1 Mb direct GPU scaffold-split diagnostic, analogous to the 003-style P9016 GPU baseline but using P1006 pairs. Training remains blind and uses only P1006 pairs; P1006 phase labels and CHARM/3DG are eval-only.

- full result root: `/mnt/ssd/zliu/phase3/test_res/013-20260617_115651-p1006_dscale_semantics_gpu_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/013-20260617_115651-p1006_dscale_semantics_gpu_1m`

| config | sample | backend | d_scale_mode | formula | top1 all | top1 cis | top1 trans | same/cross cis | cis Spearman | entropy | pU | min sep | mean sep |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `p1006_current_raw_gpu_exact_ieps0p5_noise0_seed17` | P1006 | gpu | raw_count | n_raw | 0.440298 | 0.498652 | 0.371719 | 0.878799 | 0.566779 | 0.7545439 | 0.544288397 | 0.0769124925 | 3.74474788 |
| `p1006_posterior_count_gpu_exact_ieps0p5_noise0_seed17` | P1006 | gpu | posterior_count | n_raw*posterior_prob | 0.451744 | 0.515881 | 0.376369 | 0.921309 | 0.619957 | 0.65221709 | 0.470475167 | 0.624530017 | 4.47352648 |
