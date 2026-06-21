# p9016_pcgamma1_confboost_pmax_margin_top0p05_b5_msep1_lsep1_ieps0p5_noise0_seed17

- label: `p9016_pcgamma1_confboost_pmax_margin_top0p05_b5_msep1_lsep1_ieps0p5_noise0_seed17`
- created_at: `2026-06-19T12:07:22`
- training input: raw P9016 pairs: `/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- reconstruction: `/mnt/ssd/zliu/phase3/test_res/031-20260619_105814-p9016_confidence_boost_pairs_1m/outputs/p9016_pcgamma1_confboost_pmax_margin_top0p05_b5_msep1_lsep1_ieps0p5_noise0_seed17/p9016_full.coords.tsv`
- CHARM/3DG eval reference: `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz`
- train manifest: `/mnt/ssd/zliu/phase3/test_res/031-20260619_105814-p9016_confidence_boost_pairs_1m/outputs/p9016_pcgamma1_confboost_pmax_margin_top0p05_b5_msep1_lsep1_ieps0p5_noise0_seed17/p9016_full.manifest.tsv`
- boundary: training used raw P9016 contact information only; phase labels and CHARM/3DG were read only by this post-training evaluator.
- copy gauge: structure plots and distance metrics use per-chromosome cis distance-matrix Spearman correlation to select a geometry gauge. Contact identity metrics report the reconstruction under a whole-chromosome SNP cis-top1 oracle gauge, which is eval-only and exists because copy0/copy1 are gauge labels that can be swapped independently per chromosome.
- contact denominator: contact accuracy uses eval-only raw contacts with both `phase0` and `phase1`, excluding same-bin contacts, and requiring a matching posterior bpair.
- alignment: 3D scatter plots and Procrustes RMSD use rigid alignment only: translation and rotation are fitted, reconstruction scale is not fitted to CHARM/3DG. The reported similarity scale is diagnostic only and is not applied.

## Allele Separation

| source | mean copy0/copy1 separation |
| --- | ---: |
| CHARM/3DG | 4.65535 |
| reconstruction | 4.56489 |

## Contact Accuracy

| scope | top1 accuracy |
| --- | ---: |
| all | 0.438246 |
| cis | 0.540797 |
| trans | 0.291477 |

| scope | pmax >= 0.9 accuracy | pmax >= 0.9 recall |
| --- | ---: | ---: |
| all | 0.465415 | 0.169084 |
| cis | 0.58048 | 0.191714 |
| trans | 0.332944 | 0.136695 |

## Cis Distance Matrix Correlation

Mean per-chromosome Spearman correlation, split by CHARM/3DG copy and reconstruction copy.

| CHARM/3DG copy \\ reconstruction copy | copy0 | copy1 |
| --- | ---: | ---: |
| copy0 | 0.424983 | 0.406421 |
| copy1 | 0.417139 | 0.455562 |

## Contact Distance Distribution

Uses two contact sets. The SNP-labeled row uses eval-only binned contacts with SNP phase labels on both ends, excluding same-bin contacts, and requiring a matching posterior bpair. For cis SNP-labeled contacts, the SNP phase selects the copy pair after the whole-chromosome SNP gauge is applied to reconstruction. For trans SNP-labeled contacts, SNP 0/1 is not a shared cross-chromosome copy gauge, so the plotted distance is the closest copy-pair distance for that chromosome pair in each structure. The posterior-top1 row uses all non-same-bin posterior binned contacts, weighted by raw contact count. For posterior-top1 contacts, cis uses the model top1 copy pair with the same eval-only contact gauge for CHARM/3DG; trans uses the model top1 copy-pair distance in reconstruction and the closest CHARM/3DG copy-pair distance as a gauge-safe reference.

| contact set | scope | contacts | CHARM/3DG mean | reconstruction mean | CHARM/3DG median | reconstruction median |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| SNP-labeled | cis | 207084 | 1.05109 | 0.94862 | 0.890171 | 0.601489 |
| SNP-labeled | trans | 143977 | 2.26114 | 1.93124 | 2.18208 | 1.92785 |
| posterior top1 | cis | 675367 | 0.943408 | 0.527395 | 0.674384 | 0.471368 |
| posterior top1 | trans | 948013 | 1.80283 | 1.52615 | 1.44928 | 1.34379 |

## Output Tables

- `cis_distance_correlation_matrix.tsv`: per-chromosome 2x2 cis distance-matrix Pearson/Spearman table for CHARM/3DG copy0/copy1 against reconstruction copy0/copy1.
- `cis_distance_correlations.tsv`: per-chromosome cis distance-matrix Pearson/Spearman for the two whole-chromosome swap choices, with the selected geometry swap marked.
- `contact_accuracy.tsv`: four-state top1 accuracy, same/cross accuracy, pmax >= 0.9 accuracy, called fraction, and recall for all/cis/trans contacts plus per-chromosome cis contacts.
- `contact_distance_distribution.tsv`: cis/trans 3D distance summaries for SNP-labeled contacts and posterior-top1 contacts in CHARM/3DG and reconstruction.
- `copy_separation.tsv`: per-chromosome and genome mean/median distance between copy0 and copy1 of the same bin.

## Plots

- `plots/chr1_distance_maps.png`: rows are CHARM/3DG and reconstruction; columns are copy0 and copy1; larger distances are blue.
- `plots/all_chrom_3d_scatter.png`: each point is one bin; chromosome+copy states are colored separately; reconstruction is rigid-Procrustes-aligned into the CHARM/3DG coordinate frame with no scale fitting, and both panels share one 3D coordinate range.
- `plots/chr1_copy_3d_scatter.png`: chr1 copy0/copy1 bins highlighted on top of low-alpha non-chr1 background points, using the same global rigid transform as the all-chromosome scatter plot.
- `plots/contact_distance_histograms.png`: 2x2 cis/trans histograms for SNP-labeled and posterior-top1 contact distances; CHARM/3DG and reconstruction use different colors.
