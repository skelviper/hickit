# 200k

- label: `P9016 softall chain 200k`
- created_at: `2026-06-16T11:01:08`
- training input: `/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- reconstruction: `/mnt/ssd/zliu/phase3/test_res/002-20260616_003319-p9016_softall_4m_1m_200k_chain/outputs/200k/minimal_soft_sep_off/p9016_full.coords.tsv`
- CHARM/3DG eval reference: `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.200k.3dg.gz`
- train manifest: `/mnt/ssd/zliu/phase3/test_res/002-20260616_003319-p9016_softall_4m_1m_200k_chain/outputs/200k/minimal_soft_sep_off/p9016_full.manifest.tsv`
- boundary: training used raw P9016 contact information only; phase labels and CHARM/3DG were read only by this post-training evaluator.
- copy gauge: evaluation first selects one copy swap per chromosome by per-chromosome cis distance-matrix Spearman correlation. Contact top1 and pmax metrics then use this fixed geometry-selected gauge when comparing four-state probabilities to SNP phase truth.
- contact denominator: contact accuracy uses eval-only raw contacts with both `phase0` and `phase1`, excluding same-bin contacts, and requiring a matching posterior bpair.
- CHARM/3DG probability baseline: `charm3dg_uniform_prior_fdg` recomputes four-state probabilities from CHARM/3DG distances using hickit FDG contact energy, posterior `base_d_scale/base_k`, CHARM/3DG reference copy labels, and uniform four-state prior because posterior log-priors are not exported.

- alignment: 3D scatter plots and Procrustes RMSD use rigid alignment only: translation and rotation are fitted, reconstruction scale is not fitted to CHARM/3DG. The reported similarity scale is diagnostic only and is not applied.

## Quantitative Results

| metric | value |
| --- | ---: |
| bin_size_bp | 200000 |
| pairs_total | 1703888 |
| pairs_cis | 1135454 |
| pairs_trans | 568434 |
| pairs_with_any_phase_eval_only | 1302203 |
| pairs_with_both_phases_eval_only | 496021 |
| reference_points_aggregated | 24757 |
| reconstruction_points | 26338 |
| shared_points_per_chrom_best | 24757 |
| chr1_shared_points_per_chrom_best | 1898 |
| mean_per_chrom_cis_distance_spearman | 0.694382 |
| chr1_cis_distance_spearman | 0.762259 |
| model_top1_accuracy_genome_all | 0.384278 |
| model_pmax90_accuracy_genome_all | 0.319269 |
| model_pmax90_recall_genome_all | 0.117764 |
| model_top1_accuracy_genome_cis | 0.468564 |
| model_pmax90_accuracy_genome_cis | 0.438143 |
| model_pmax90_recall_genome_cis | 0.11574 |
| model_top1_accuracy_genome_trans | 0.229915 |
| model_pmax90_accuracy_genome_trans | 0.216682 |
| model_pmax90_recall_genome_trans | 0.12147 |
| charm3dg_top1_accuracy_genome_all | 0.191512 |
| charm3dg_pmax90_accuracy_genome_all | 0.0698055 |
| charm3dg_pmax90_recall_genome_all | 0.0274965 |
| charm3dg_top1_accuracy_genome_cis | 0.206266 |
| charm3dg_pmax90_accuracy_genome_cis | 0.00856731 |
| charm3dg_pmax90_recall_genome_cis | 0.00246224 |
| charm3dg_top1_accuracy_genome_trans | 0.164491 |
| charm3dg_pmax90_accuracy_genome_trans | 0.124534 |
| charm3dg_pmax90_recall_genome_trans | 0.0733444 |
| reconstruction_mean_copy01_separation | 9.49272 |
| charm3dg_mean_copy01_separation | 9.55054 |
| reconstruction_genome_volume | 5451.07 |
| charm3dg_genome_volume | 3317.86 |

## Output Tables

- `cis_distance_correlations.tsv`: per-chromosome cis distance-matrix Pearson/Spearman for both copy swaps, with the selected per-chrom swap marked.
- `contact_accuracy.tsv`: four-state top1 accuracy, pmax >= 0.9 accuracy, called fraction, and recall for all/cis/trans contacts plus per-chromosome cis contacts. The `copy_swap_policy` column records whether rows use the reconstruction's fixed cis-distance-selected gauge or CHARM/3DG reference copy labels.
- `copy_separation.tsv`: per-chromosome and genome mean/median distance between copy0 and copy1 of the same bin.
- `per_chrom_volume.tsv`: per-chromosome and genome convex-hull volumes for CHARM/3DG and reconstruction.

## Plots

- `plots/chr1_distance_maps.png`: rows are CHARM/3DG and reconstruction; columns are copy0 and copy1; larger distances are blue.
- `plots/all_chrom_3d_scatter.png`: each point is one bin; chromosome+copy states are colored separately; reconstruction is rigid-Procrustes-aligned into the CHARM/3DG coordinate frame with no scale fitting, and both panels share one 3D coordinate range.
- `plots/chr1_copy_3d_scatter.png`: chr1 copy0/copy1 bins highlighted on top of low-alpha non-chr1 background points, using the same global rigid transform as the all-chromosome scatter plot.
