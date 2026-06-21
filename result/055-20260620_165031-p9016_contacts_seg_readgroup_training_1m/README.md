# 055 P9016 contacts.seg Readgroup Training 1Mb

This experiment tests contacts.seg-derived hickit-like/readgroup sources for blind P9016 1Mb reconstruction. It is a source/interface diagnostic, not a replacement for the strict approved P9016.pairs.gz baseline.

## Paths

- full result root: `/mnt/ssd/zliu/phase3/test_res/055-20260620_165031-p9016_contacts_seg_readgroup_training_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/055-20260620_165031-p9016_contacts_seg_readgroup_training_1m`
- summary: `/mnt/ssd/zliu/phase3/test_res/055-20260620_165031-p9016_contacts_seg_readgroup_training_1m/summary.tsv`
- source delta summary: `/mnt/ssd/zliu/phase3/test_res/055-20260620_165031-p9016_contacts_seg_readgroup_training_1m/source_delta_summary.tsv`
- readgroup matched delta summary: `/mnt/ssd/zliu/phase3/test_res/055-20260620_165031-p9016_contacts_seg_readgroup_training_1m/readgroup_matched_delta_summary.tsv`

## Boundary

- `approved_baseline_pairs_only` is the only strict approved `/shared/.../P9016.pairs.gz` baseline row.
- All `hickit_like_*` and `readgroup_*` rows use 054 contacts.seg-derived custom pairs with phase labels dropped. They are blind-safe upstream-contact diagnostics, but they are not approved-pairs-only baseline rows.
- SNP phase and CHARM/3DG are used only by post-training eval through `scripts/p9016_common_eval.sh`.
- copy0/copy1 remain gauge labels.
- `readgroup_adjacent_pairwise` originally stalled under GPU backend and was rerun with CPU backend using the same run-local binary and model settings; compare its numeric result with that backend caveat.

## Main Results

| config | source | readgroup | backend | top1 all | top1 cis | top1 trans | same/cross cis | cis Spearman | entropy | pU | delta trans vs approved |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| approved_baseline_pairs_only | approved_pairs | off | gpu | 0.466095 | 0.569915 | 0.317509 | 0.969974 | 0.80181 | 0.61166954 | 0.441226333 | 0 |
| hickit_like_dup100_pairwise | hickit_like_dup100 | off | gpu | 0.424893 | 0.524292 | 0.282635 | 0.97028 | 0.7801 | 0.378401071 | 0.272958666 | -0.034874 |
| readgroup_adjacent_joint | readgroup_adjacent | joint_marginal | gpu | 0.422456 | 0.522302 | 0.279559 | 0.97195 | 0.799923 | 0.772842586 | 0.557488084 | -0.03795 |
| readgroup_adjacent_pairwise | readgroup_adjacent | off | cpu | 0.409486 | 0.520235 | 0.250985 | 0.967558 | 0.789474 | 0.793431222 | 0.572339654 | -0.066524 |
| readgroup_allcomb_joint | readgroup_allcomb | joint_marginal | gpu | 0.396436 | 0.487477 | 0.26614 | 0.929 | 0.723715 | 0.829507172 | 0.598362923 | -0.051369 |
| readgroup_allcomb_pairwise | readgroup_allcomb | off | gpu | 0.399584 | 0.517138 | 0.231343 | 0.969489 | 0.779874 | 0.650420904 | 0.469179511 | -0.086166 |

## Matched Readgroup Deltas

| joint config | pairwise config | joint backend | pairwise backend | delta top1 all | delta top1 cis | delta top1 trans | delta cis Spearman |
|---|---|---|---|---:|---:|---:|---:|
| readgroup_adjacent_joint | readgroup_adjacent_pairwise | gpu | cpu | 0.01297 | 0.002067 | 0.028574 | 0.010449 |
| readgroup_allcomb_joint | readgroup_allcomb_pairwise | gpu | gpu | -0.003148 | -0.029661 | 0.034797 | -0.056159 |

## Interpretation

The contacts.seg-derived rows do not support a trans improvement claim over the approved baseline. The readgroup joint rows must be judged against matched same-source pairwise controls: joint marginal improves trans relative to the contacts.seg pairwise controls, but both joint rows remain below the approved baseline. The allcomb joint row also damages cis accuracy and cis distance Spearman, and the adjacent matched comparison has a backend caveat because the pairwise control had to be rerun on CPU after the GPU path stalled.

The hickit-like source and readgroup source change the upstream contact distribution. Any difference versus approved baseline is a source effect unless isolated by the matched pairwise/joint comparison.
