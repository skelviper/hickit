# 019-20260618_224526-p9016_copytrack_continuity_1m

Controlled blind P9016 copy-track continuity experiment at 1 Mb. The experiment tests whether a homolog-vector continuity force can reduce local copy-track flipping and improve trans phase accuracy.

## Paths

- full result root: `/tmp/hk_blind_test_res_019_smoke/019-20260618_224526-p9016_copytrack_continuity_1m`
- lightweight hickit result: `/mnt/ssd/zliu/phase3/hickit/result/019-20260618_224526-p9016_copytrack_continuity_1m`
- summary.tsv: `/tmp/hk_blind_test_res_019_smoke/019-20260618_224526-p9016_copytrack_continuity_1m/summary.tsv`
- trans_delta.tsv: `/tmp/hk_blind_test_res_019_smoke/019-20260618_224526-p9016_copytrack_continuity_1m/trans_delta.tsv`
- commands.log: `/tmp/hk_blind_test_res_019_smoke/019-20260618_224526-p9016_copytrack_continuity_1m/commands.log`

## Training Boundary

Training uses only P9016 pairs. SNP/phase labels and CHARM/3DG are used only in eval through `scripts/p9016_common_eval.sh`. The copy-track force is blind geometry: adjacent same-chromosome homolog vectors are encouraged to be continuous.

## Git And Build

- git commit: `1037ee5f79d54e943eef5c37e5aeb2d841698cb8`
- git dirty count at start: `30`
- WARNING: dirty tree at launch; inspect `logs/git_status.txt` and `logs/git_diff.patch`.
- run binary sha256: `7c4fc219ea7c722bd132164b7f4e7a48fe4cdf45cae31e682ab56c8c586f4ded`
- backend: `cpu`

## Main Results

| config | gamma | min_sep | lambda_sep | lambda_copytrack | top1 all | top1 cis | top1 trans | delta trans | cis Spearman | entropy | pU | sep_p05 | cos<0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `p9016_pcgamma1_copytrack0_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | 0 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 0 |
| `p9016_pcgamma1_copytrack0p01_sep_off_ieps0p5_noise0_seed17` | 1 | 0 | 0 | 0.00999999978 | NA | NA | NA | NA | NA | 1.32029605 | 0.95239228 | 1.00231004 | 0 |

## Interpretation

- The +0.1 trans target was not met in this sweep; continue with chromosome-level gauge synchronization or E-step prior diagnostics.
- Accept copy-track continuity only if trans improves without a large cis/Spearman/entropy penalty. A drop in copytrack cos<0 without trans improvement means internal chromosome copy tracks stabilized but cross-chromosome gauge remains unresolved.
