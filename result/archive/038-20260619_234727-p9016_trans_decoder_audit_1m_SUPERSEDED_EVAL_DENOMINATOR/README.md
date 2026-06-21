# 038 P9016 Trans Decoder Audit

This is a post-training diagnostic. It does not retrain Hickit. It asks whether current reconstruction/posterior quantities contain a blind-selectable trans copy-state decoder that can close the +0.1 trans top1 gap.

- full result root: `/mnt/ssd/zliu/phase3/test_res/038-20260619_234727-p9016_trans_decoder_audit_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/038-20260619_234727-p9016_trans_decoder_audit_1m`
- headline: `NO_BLIND_DECODER_PLUS_0P1`
- training boundary: no training is run here; SNP labels and CHARM/3DG are eval-only.
- blind decoder policies use only reconstruction geometry and posterior confidence for selection.
- oracle decoder policy uses SNP truth for selection and is an eval-only ceiling.

## Main Table

| config_name | policy | flip_source | top1_accuracy | same_cross_accuracy | delta_top1_vs_cis_selected | target_plus_0p1_met |
|---|---|---|---|---|---|---|
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | pair_independent_4state_oracle | snp_truth_oracle | 0.418212631 | 0.607215041 | 0.0652048591 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.353007772 | 0.552720226 | 0 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_geom_nearest_to_posterior_top | blind_posterior_reconstruction_geometry | 0.353007772 | 0.552720226 | 0 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_geom_nearest_gap_weighted | blind_posterior_reconstruction_geometry | 0.353007772 | 0.552720226 | 0 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_margin_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.353007772 | 0.552720226 | 0 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_pmax_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.353007772 | 0.552720226 | 0 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_same_cross_distance_sum | blind_posterior_reconstruction_geometry | 0.290233857 | 0.484494051 | -0.0627739153 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_chrom_center_same_closer | blind_reconstruction_geometry | 0.28940039 | 0.483868951 | -0.0636073817 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_posterior_mass_same_cross | blind_posterior_reconstruction_geometry | 0.286434639 | 0.484591289 | -0.0665731332 | 0 |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_chrom_homolog_vector_dot | blind_reconstruction_geometry | 0.282204797 | 0.477194274 | -0.0708029755 | 0 |
| p9016_pcgamma1_td0p5_mstep0_sep0 | pair_independent_4state_oracle | snp_truth_oracle | 0.417309709 | 0.603804774 | 0.0729352605 | 0 |
| p9016_pcgamma1_td0p5_mstep0_sep0 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.344374449 | 0.538391549 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_geom_nearest_to_posterior_top | blind_posterior_reconstruction_geometry | 0.344374449 | 0.538391549 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_geom_nearest_gap_weighted | blind_posterior_reconstruction_geometry | 0.344374449 | 0.538391549 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_margin_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.344374449 | 0.538391549 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_pmax_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.344374449 | 0.538391549 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_same_cross_distance_sum | blind_posterior_reconstruction_geometry | 0.29353994 | 0.487765407 | -0.0508345083 | 0 |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_chrom_center_same_closer | blind_reconstruction_geometry | 0.292637018 | 0.486716628 | -0.0517374303 | 0 |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_posterior_mass_same_cross | blind_posterior_reconstruction_geometry | 0.288337721 | 0.486779138 | -0.0560367281 | 0 |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_chrom_homolog_vector_dot | blind_reconstruction_geometry | 0.282392327 | 0.47757628 | -0.0619821221 | 0 |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | pair_independent_4state_oracle | snp_truth_oracle | 0.446501872 | 0.628739313 | 0.109211888 | 1 |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.337289984 | 0.53584253 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_geom_nearest_to_posterior_top | blind_posterior_reconstruction_geometry | 0.337289984 | 0.53584253 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_geom_nearest_gap_weighted | blind_posterior_reconstruction_geometry | 0.337289984 | 0.53584253 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_margin_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.337289984 | 0.53584253 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_pmax_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.337289984 | 0.53584253 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_chrom_homolog_vector_dot | blind_reconstruction_geometry | 0.290060218 | 0.511171923 | -0.0472297659 | 0 |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_chrom_center_same_closer | blind_reconstruction_geometry | 0.285330296 | 0.502205213 | -0.051959688 | 0 |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_same_cross_distance_sum | blind_posterior_reconstruction_geometry | 0.2850108 | 0.493155157 | -0.0522791835 | 0 |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_posterior_mass_same_cross | blind_posterior_reconstruction_geometry | 0.281357439 | 0.482132563 | -0.0559325448 | 0 |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | pair_independent_4state_oracle | snp_truth_oracle | 0.444571008 | 0.630496538 | 0.108878501 | 1 |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.335692506 | 0.537169131 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_geom_nearest_to_posterior_top | blind_posterior_reconstruction_geometry | 0.335692506 | 0.537169131 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_geom_nearest_gap_weighted | blind_posterior_reconstruction_geometry | 0.335692506 | 0.537169131 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_margin_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.335692506 | 0.537169131 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_pmax_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.335692506 | 0.537169131 | 0 | 0 |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_chrom_homolog_vector_dot | blind_reconstruction_geometry | 0.298450447 | 0.508879891 | -0.0372420595 | 0 |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_chrom_center_same_closer | blind_reconstruction_geometry | 0.292227231 | 0.500246567 | -0.0434652757 | 0 |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_same_cross_distance_sum | blind_posterior_reconstruction_geometry | 0.283788383 | 0.494599832 | -0.0519041236 | 0 |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_posterior_mass_same_cross | blind_posterior_reconstruction_geometry | 0.263305945 | 0.460934733 | -0.0723865617 | 0 |

## Interpretation

If a blind decoder reaches +0.1, the next step is to turn that decoder into a controlled training or post-processing candidate. If only the SNP oracle improves, the current posterior/geometry does not contain enough blind information to select the needed trans copy-state gauge.
