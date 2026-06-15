# 001-20260615_214533-p9016_softall_4m_baseline

- label: `P9016 softall baseline 4000000bp`
- created_at: `2026-06-15T23:17:26`
- training input: `/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- reconstruction: `test_res/001-20260615_214533-p9016_softall_4m_baseline/outputs/minimal_soft_sep_off/p9016_full.coords.tsv`
- CHARM/3DG eval reference: `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz`
- train manifest: `test_res/001-20260615_214533-p9016_softall_4m_baseline/outputs/minimal_soft_sep_off/p9016_full.manifest.tsv`
- boundary: training used raw P9016 contact information only; phase labels and CHARM/3DG were read only by this post-training evaluator.
- copy gauge: metrics were computed for global swaps; headline metrics and standard plots use the best global swap. A per-chromosome local best-swap row is included only as an eval-only gauge-sensitivity diagnostic.

- alignment: 3D scatter plots and Procrustes RMSD use rigid alignment only: translation and rotation are fitted, reconstruction scale is not fitted to CHARM/3DG. The reported similarity scale is diagnostic only and is not applied.

## Quantitative Results

| metric | value |
| --- | ---: |
| bin_size_bp | 4000000 |
| pairs_total | 1703888 |
| pairs_cis | 1135454 |
| pairs_trans | 568434 |
| pairs_with_any_phase_eval_only | 1302203 |
| reference_points_aggregated | 1304 |
| reconstruction_points | 1318 |
| shared_points_best_global_swap | 1295 |
| chr1_shared_points_best_global_swap | 98 |
| best_global_copy_swap | 1 |
| chr1_copy_swap_for_standard_plots | 1 |
| headline_copy_swap_policy | global_swap_1 |
| best_genome_distance_spearman | 0.167402 |
| best_genome_distance_rmse_rg_norm | 0.714589 |
| best_genome_rigid_procrustes_rmsd_rg_norm | 1.01952 |
| best_genome_similarity_scale_to_reference_diagnostic | 0.543722 |
| best_chr1_distance_spearman | 0.703538 |
| best_chr1_distance_rmse_rg_norm | 0.529354 |
| best_chr1_rigid_procrustes_rmsd_rg_norm | 0.607649 |
| best_chr1_similarity_scale_to_reference_diagnostic | 0.845744 |

## Swap-Aware Metrics

| scope | copy_swap_policy | n_points | distance_spearman | distance_rmse_rg_norm | rigid_procrustes_rmsd_rg_norm | similarity_scale_to_reference_diagnostic |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| genome | global_swap_0 | 1295 | 0.103151 | 0.730969 | 1.09196 | 0.386927 |
| chr1 | global_swap_0 | 98 | 0.696762 | 0.542847 | 0.645452 | 0.819779 |
| genome | global_swap_1 | 1295 | 0.167402 | 0.714589 | 1.01952 | 0.543722 |
| chr1 | global_swap_1 | 98 | 0.703538 | 0.529354 | 0.607649 | 0.845744 |
| genome | per_chrom_best | 1295 | 0.158531 | 0.716315 | 1.0811 | 0.445728 |
| chr1 | per_chrom_best | 98 | 0.703538 | 0.529354 | 0.607649 | 0.845744 |

## Plots

- `plots/chr1_distance_maps.png`: rows are CHARM/3DG and reconstruction; columns are copy0 and copy1; larger distances are blue.
- `plots/all_chrom_3d_scatter.png`: each point is one bin; chromosome+copy states are colored separately; reconstruction is rigid-Procrustes-aligned into the CHARM/3DG coordinate frame with no scale fitting, and both panels share one 3D coordinate range.
- `plots/chr1_copy_3d_scatter.png`: chr1 copy0/copy1 bins highlighted on top of low-alpha non-chr1 background points, using the same global rigid transform as the all-chromosome scatter plot.
