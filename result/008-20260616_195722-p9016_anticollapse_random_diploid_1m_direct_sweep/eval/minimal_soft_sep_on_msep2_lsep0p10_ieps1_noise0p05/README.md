# minimal_soft_sep_on_msep2_lsep0p10_ieps1_noise0p05

- label: `P9016 anticollapse random_diploid 1Mb direct minimal_soft_sep_on_msep2_lsep0p10_ieps1_noise0p05`
- created_at: `2026-06-16T20:46:48`
- training input: `/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- reconstruction: `/mnt/ssd/zliu/phase3/test_res/008-20260616_195722-p9016_anticollapse_random_diploid_1m_direct_sweep/outputs/minimal_soft_sep_on_msep2_lsep0p10_ieps1_noise0p05/p9016_full.coords.tsv`
- CHARM/3DG eval reference: `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz`
- train manifest: `/mnt/ssd/zliu/phase3/test_res/008-20260616_195722-p9016_anticollapse_random_diploid_1m_direct_sweep/outputs/minimal_soft_sep_on_msep2_lsep0p10_ieps1_noise0p05/p9016_full.manifest.tsv`
- boundary: training used raw P9016 contact information only; phase labels and CHARM/3DG were read only by this post-training evaluator.
- copy gauge: structure plots and distance metrics use per-chromosome cis distance-matrix Spearman correlation to select a geometry gauge. Contact identity metrics report the reconstruction under a whole-chromosome SNP cis-top1 oracle gauge, which is eval-only and exists because copy0/copy1 are gauge labels that can be swapped independently per chromosome.
- contact denominator: contact accuracy uses eval-only raw contacts with both `phase0` and `phase1`, excluding same-bin contacts, and requiring a matching posterior bpair.
- alignment: 3D scatter plots and Procrustes RMSD use rigid alignment only: translation and rotation are fitted, reconstruction scale is not fitted to CHARM/3DG. The reported similarity scale is diagnostic only and is not applied.

## Allele Separation

| source | mean copy0/copy1 separation |
| --- | ---: |
| CHARM/3DG | 4.65535 |
| reconstruction | 3.43074 |

## Contact Accuracy

| scope | top1 accuracy |
| --- | ---: |
| all | 0.305996 |
| cis | 0.338732 |
| trans | 0.259146 |

| scope | pmax >= 0.9 accuracy | pmax >= 0.9 recall |
| --- | ---: | ---: |
| all | 0.316902 | 0.120743 |
| cis | 0.340209 | 0.146149 |
| trans | 0.270899 | 0.0843815 |

## Cis Distance Matrix Correlation

Mean per-chromosome Spearman correlation, split by CHARM/3DG copy and reconstruction copy.

| CHARM/3DG copy \\ reconstruction copy | copy0 | copy1 |
| --- | ---: | ---: |
| copy0 | 0.23174 | 0.208367 |
| copy1 | 0.192229 | 0.237372 |

## Contact Distance Distribution

Uses two contact sets. The SNP-labeled row uses eval-only binned contacts with SNP phase labels on both ends, excluding same-bin contacts, and requiring a matching posterior bpair. For cis SNP-labeled contacts, the SNP phase selects the copy pair after the whole-chromosome SNP gauge is applied to reconstruction. For trans SNP-labeled contacts, SNP 0/1 is not a shared cross-chromosome copy gauge, so the plotted distance is the closest copy-pair distance for that chromosome pair in each structure. The posterior-top1 row uses all non-same-bin posterior binned contacts, weighted by raw contact count. For posterior-top1 contacts, cis uses the model top1 copy pair with the same eval-only contact gauge for CHARM/3DG; trans uses the model top1 copy-pair distance in reconstruction and the closest CHARM/3DG copy-pair distance as a gauge-safe reference.

| contact set | scope | contacts | CHARM/3DG mean | reconstruction mean | CHARM/3DG median | reconstruction median |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| SNP-labeled | cis | 207084 | 1.05109 | 1.74995 | 0.890171 | 1.18366 |
| SNP-labeled | trans | 143977 | 2.26114 | 2.33458 | 2.18208 | 2.18068 |
| posterior top1 | cis | 677648 | 2.24959 | 0.562146 | 1.12327 | 0.505439 |
| posterior top1 | trans | 545353 | 2.34973 | 2.41093 | 2.28405 | 2.27279 |

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
