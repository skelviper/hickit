# 047 P9016 Init-Coordinate Anchor 1Mb

This controlled blind-training experiment tests whether preserving the initial scaffold-split coordinate gauge can reduce local copy-track flipping and improve full-denominator trans top1 accuracy. It reuses the existing parent-anchor force with an identity bmap and a snapshot of the direct initial diploid coordinates.

- full result root: `/mnt/ssd/zliu/phase3/test_res/047-20260620_082138-p9016_init_coord_anchor_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/047-20260620_082138-p9016_init_coord_anchor_1m`
- headline: `NO_PLUS_0P1`
- training inputs: approved P9016 raw pairs only
- eval inputs: P9016 SNP labels and CHARM/3DG are used only after training
- anchor source: blind direct initial diploid coordinates, not phase labels or reference coordinates
- success criterion: +0.1 full-denominator trans top1 over the best035 seed17 anchor-off replay

## Main Results

| config_name | seed_family | init_coord_anchor_k | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | delta_trans_vs_same_seed_anchor_off | delta_trans_vs_best035_seed17_anchor_off | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | copytrack_frac_cos_lt_0 | anchor_force_l1_over_contact_force_l1 | target_plus_0p1_vs_best035_seed17_met |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_initanchor_k0p001_seed71_noise0p05_best035 | seed71_noise0p05 | 0.00100000005 | 0.485804 | 0.574695 | 0.358585 | 0.002924 | 0.005584 | 0.563729 | 0.744353 | 0.693435669 | 0.50020808 | 0.379367381 | 2.56006312 | 0.00344431688 | 0.000155974728 | 0 |
| p9016_initanchor_k0p003_seed71_noise0p05_best035 | seed71_noise0p05 | 0.00300000003 | 0.484162 | 0.572327 | 0.357981 | 0.00232 | 0.00498 | 0.561902 | 0.745886 | 0.69200927 | 0.499179155 | 0.438126266 | 2.59001517 | 0.00382701875 | 0.000468624667 | 0 |
| p9016_initanchor_k0p01_seed71_noise0p05_best035 | seed71_noise0p05 | 0.00999999978 | 0.483233 | 0.572191 | 0.355918 | 0.000257 | 0.002917 | 0.561874 | 0.749965 | 0.694105268 | 0.500691116 | 0.349551201 | 2.59364271 | 0.003061615 | 0.00154862853 | 0 |
| p9016_initanchor_off_seed71_noise0p05_best035 | seed71_noise0p05 | 0 | 0.483053 | 0.572065 | 0.355661 | 0 | 0.00266 | 0.56111 | 0.748189 | 0.69240886 | 0.499467403 | 0.371240407 | 2.58296108 | 0.00344431688 | 0 | 0 |
| p9016_initanchor_k0p0003_seed71_noise0p05_best035 | seed71_noise0p05 | 0.000300000014 | 0.484073 | 0.574477 | 0.354689 | -0.000972 | 0.001688 | 0.560478 | 0.750541 | 0.692574203 | 0.499586672 | 0.365777731 | 2.5639205 | 0.00344431688 | 4.65595716e-05 | 0 |
| p9016_initanchor_off_seed17_noise0_best035 | seed17_noise0 | 0 | 0.482105 | 0.572313 | 0.353001 | 0 | 0 | 0.552755 | 0.769509 | 0.69679147 | 0.502628803 | 0.460410267 | 2.83495998 | 0.0015308075 | 0 | 0 |
| p9016_initanchor_k0p0003_seed17_noise0_best035 | seed17_noise0 | 0.000300000014 | 0.480442 | 0.570017 | 0.352244 | -0.000757 | -0.000757 | 0.551977 | 0.767643 | 0.697908223 | 0.50343436 | 0.408795953 | 2.79330826 | 0.00229621125 | 4.58167388e-05 | 0 |
| p9016_initanchor_k0p03_seed17_noise0_best035 | seed17_noise0 | 0.0299999993 | 0.478836 | 0.568843 | 0.350021 | -0.00298 | -0.00298 | 0.542406 | 0.769066 | 0.704008698 | 0.507834911 | 0.478495926 | 2.71031857 | 0.0015308075 | 0.00453577259 | 0 |
| p9016_initanchor_k0p01_seed17_noise0_best035 | seed17_noise0 | 0.00999999978 | 0.478616 | 0.568891 | 0.349417 | -0.003584 | -0.003584 | 0.5411 | 0.759412 | 0.701791406 | 0.50623548 | 0.334618449 | 2.70151734 | 0.0015308075 | 0.0015183714 | 0 |
| p9016_initanchor_k0p001_seed17_noise0_best035 | seed17_noise0 | 0.00100000005 | 0.478016 | 0.5681 | 0.34909 | -0.003911 | -0.003911 | 0.548546 | 0.762925 | 0.699313104 | 0.504447758 | 0.379511237 | 2.7558136 | 0.0015308075 | 0.000152788572 | 0 |
| p9016_initanchor_k0p03_seed71_noise0p05_best035 | seed71_noise0p05 | 0.0299999993 | 0.476671 | 0.566455 | 0.348174 | -0.007487 | -0.004827 | 0.547664 | 0.777702 | 0.694042861 | 0.500646055 | 0.522811532 | 2.90213394 | 0.003061615 | 0.0046025395 | 0 |
| p9016_initanchor_k0p003_seed17_noise0_best035 | seed17_noise0 | 0.00300000003 | 0.478336 | 0.571488 | 0.34502 | -0.007981 | -0.007981 | 0.540517 | 0.767616 | 0.699067175 | 0.504270375 | 0.408593953 | 2.79019547 | 0.0015308075 | 0.000457137128 | 0 |

## Interpretation Boundary

A positive result would mean the current trans gap is sensitive to training-time gauge drift from the initial scaffold split. It would not by itself establish a final model, because anchor_k selection still needs a blind rule. A negative result means simply preserving the initial coordinate gauge is not enough to recover the missing trans identity signal.
