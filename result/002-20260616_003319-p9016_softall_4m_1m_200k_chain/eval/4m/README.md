# 4m

- label: `P9016 softall chain 4m`
- created_at: `2026-06-16T10:59:56`
- training input: `/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- reconstruction: `/mnt/ssd/zliu/phase3/test_res/002-20260616_003319-p9016_softall_4m_1m_200k_chain/outputs/4m/minimal_soft_sep_off/p9016_full.coords.tsv`
- CHARM/3DG eval reference: `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz`
- train manifest: `/mnt/ssd/zliu/phase3/test_res/002-20260616_003319-p9016_softall_4m_1m_200k_chain/outputs/4m/minimal_soft_sep_off/p9016_full.manifest.tsv`
- boundary: training used raw P9016 contact information only; phase labels and CHARM/3DG were read only by this post-training evaluator.
- copy gauge: evaluation first selects one copy swap per chromosome by per-chromosome cis distance-matrix Spearman correlation. Contact top1 and pmax metrics then use this fixed geometry-selected gauge when comparing four-state probabilities to SNP phase truth.
- contact denominator: contact accuracy uses eval-only raw contacts with both `phase0` and `phase1`, excluding same-bin contacts, and requiring a matching posterior bpair.
- CHARM/3DG probability baseline: `charm3dg_uniform_prior_fdg` recomputes four-state probabilities from CHARM/3DG distances using hickit FDG contact energy, posterior `base_d_scale/base_k`, CHARM/3DG reference copy labels, and uniform four-state prior because posterior log-priors are not exported.

- alignment: 3D scatter plots and Procrustes RMSD use rigid alignment only: translation and rotation are fitted, reconstruction scale is not fitted to CHARM/3DG. The reported similarity scale is diagnostic only and is not applied.

## Quantitative Results

| metric | value |
| --- | ---: |
| bin_size_bp | 4000000 |
| pairs_total | 1703888 |
| pairs_cis | 1135454 |
| pairs_trans | 568434 |
| pairs_with_any_phase_eval_only | 1302203 |
| pairs_with_both_phases_eval_only | 496021 |
| reference_points_aggregated | 1304 |
| reconstruction_points | 1318 |
| shared_points_per_chrom_best | 1295 |
| chr1_shared_points_per_chrom_best | 98 |
| mean_per_chrom_cis_distance_spearman | 0.683784 |
| chr1_cis_distance_spearman | 0.804614 |
| model_top1_accuracy_genome_all | 0.352892 |
| model_pmax90_accuracy_genome_all | 0.362528 |
| model_pmax90_recall_genome_all | 0.194183 |
| model_top1_accuracy_genome_cis | 0.455591 |
| model_pmax90_accuracy_genome_cis | 0.441315 |
| model_pmax90_recall_genome_cis | 0.274932 |
| model_top1_accuracy_genome_trans | 0.23812 |
| model_pmax90_accuracy_genome_trans | 0.237298 |
| model_pmax90_recall_genome_trans | 0.103941 |
| charm3dg_top1_accuracy_genome_all | 0.214292 |
| charm3dg_pmax90_accuracy_genome_all | 0.142124 |
| charm3dg_pmax90_recall_genome_all | 0.0783571 |
| charm3dg_top1_accuracy_genome_cis | 0.255442 |
| charm3dg_pmax90_accuracy_genome_cis | 0.159808 |
| charm3dg_pmax90_recall_genome_cis | 0.0917688 |
| charm3dg_top1_accuracy_genome_trans | 0.168305 |
| charm3dg_pmax90_accuracy_genome_trans | 0.120537 |
| charm3dg_pmax90_recall_genome_trans | 0.0633688 |
| reconstruction_mean_copy01_separation | 3.11714 |
| charm3dg_mean_copy01_separation | 4.6402 |
| reconstruction_genome_volume | 240.892 |
| charm3dg_genome_volume | 301.031 |

## Output Tables

- `cis_distance_correlations.tsv`: per-chromosome cis distance-matrix Pearson/Spearman for both copy swaps, with the selected per-chrom swap marked.
- `contact_accuracy.tsv`: four-state top1 accuracy, pmax >= 0.9 accuracy, called fraction, and recall for all/cis/trans contacts plus per-chromosome cis contacts. The `copy_swap_policy` column records whether rows use the reconstruction's fixed cis-distance-selected gauge or CHARM/3DG reference copy labels.
- `copy_separation.tsv`: per-chromosome and genome mean/median distance between copy0 and copy1 of the same bin.
- `per_chrom_volume.tsv`: per-chromosome and genome convex-hull volumes for CHARM/3DG and reconstruction.

## Plots

- `plots/chr1_distance_maps.png`: rows are CHARM/3DG and reconstruction; columns are copy0 and copy1; larger distances are blue.
- `plots/all_chrom_3d_scatter.png`: each point is one bin; chromosome+copy states are colored separately; reconstruction is rigid-Procrustes-aligned into the CHARM/3DG coordinate frame with no scale fitting, and both panels share one 3D coordinate range.
- `plots/chr1_copy_3d_scatter.png`: chr1 copy0/copy1 bins highlighted on top of low-alpha non-chr1 background points, using the same global rigid transform as the all-chromosome scatter plot.
