# 038 P9016 Trans Decoder Audit

This is a post-training diagnostic. It does not retrain Hickit. It asks whether current reconstruction/posterior quantities contain a blind-selectable trans copy-state decoder that can close the +0.1 trans top1 gap.

- full result root: `/mnt/ssd/zliu/phase3/test_res/038-20260620_000849-p9016_trans_decoder_audit_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/038-20260620_000849-p9016_trans_decoder_audit_1m`
- headline: `NO_RAW_BLIND_DECODER_PLUS_0P1`
- training boundary: no training is run here; source runs are required to have `input_contact_source=raw_pairs`, `uses_phase_labels=0`, and `uses_charm_or_reference=0`.
- base gauge: `cis_selected_whole_chrom` uses the standard eval-only whole-chromosome SNP cis oracle gauge, matching the existing evaluator headline.
- blind decoder boundary: blind policies select chromosome-pair swaps using raw posterior trans bpair weights (`n_raw`) and reconstruction geometry/posterior fields only.
- scoring boundary: SNP labels and CHARM/3DG are used only after decoder selection for scoring and shared denominator matching.
- oracle decoder policy uses SNP truth for pair selection and is an eval-only ceiling.

## Main Table

| config_name | policy | flip_source | top1_accuracy | same_cross_accuracy | delta_top1_vs_cis_selected | target_plus_0p1_met | uses_snp_for_pair_decoder_selection | uses_snp_labeled_denominator_for_selection_weight | decoder_selection_denominator |
|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | pair_independent_4state_oracle | snp_truth_oracle | 0.418212631 | 0.607215041 | 0.0652048591 | 0 | 1 | 1 | snp_truth_eval_contacts |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.353007772 | 0.552720226 | 0 | 0 | 0 | 0 | none_cis_snp_gauge_only |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_geom_nearest_to_posterior_top | blind_posterior_reconstruction_geometry | 0.353007772 | 0.552720226 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_geom_nearest_gap_weighted | blind_posterior_reconstruction_geometry | 0.353007772 | 0.552720226 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_margin_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.353007772 | 0.552720226 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_pmax_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.353007772 | 0.552720226 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_same_cross_distance_sum | blind_posterior_reconstruction_geometry | 0.290421387 | 0.4847302 | -0.0625863853 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_chrom_center_same_closer | blind_reconstruction_geometry | 0.28940039 | 0.483868951 | -0.0636073817 | 0 | 0 | 0 | raw_reconstruction_chrom_centers_plus_raw_trans_pair_set |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_posterior_mass_same_cross | blind_posterior_reconstruction_geometry | 0.286434639 | 0.484591289 | -0.0665731332 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | blind_chrom_homolog_vector_dot | blind_reconstruction_geometry | 0.282204797 | 0.477194274 | -0.0708029755 | 0 | 0 | 0 | raw_reconstruction_chrom_centers_plus_raw_trans_pair_set |
| p9016_pcgamma1_td0p5_mstep0_sep0 | pair_independent_4state_oracle | snp_truth_oracle | 0.417309709 | 0.603804774 | 0.0729352605 | 0 | 1 | 1 | snp_truth_eval_contacts |
| p9016_pcgamma1_td0p5_mstep0_sep0 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.344374449 | 0.538391549 | 0 | 0 | 0 | 0 | none_cis_snp_gauge_only |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_geom_nearest_to_posterior_top | blind_posterior_reconstruction_geometry | 0.344374449 | 0.538391549 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_geom_nearest_gap_weighted | blind_posterior_reconstruction_geometry | 0.344374449 | 0.538391549 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_margin_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.344374449 | 0.538391549 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_pmax_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.344374449 | 0.538391549 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_same_cross_distance_sum | blind_posterior_reconstruction_geometry | 0.293532995 | 0.487716788 | -0.0508414538 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_chrom_center_same_closer | blind_reconstruction_geometry | 0.292637018 | 0.486716628 | -0.0517374303 | 0 | 0 | 0 | raw_reconstruction_chrom_centers_plus_raw_trans_pair_set |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_posterior_mass_same_cross | blind_posterior_reconstruction_geometry | 0.289073949 | 0.486681901 | -0.0553004994 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | blind_chrom_homolog_vector_dot | blind_reconstruction_geometry | 0.282392327 | 0.47757628 | -0.0619821221 | 0 | 0 | 0 | raw_reconstruction_chrom_centers_plus_raw_trans_pair_set |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | pair_independent_4state_oracle | snp_truth_oracle | 0.446501872 | 0.628739313 | 0.109211888 | 1 | 1 | 1 | snp_truth_eval_contacts |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.337289984 | 0.53584253 | 0 | 0 | 0 | 0 | none_cis_snp_gauge_only |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_geom_nearest_to_posterior_top | blind_posterior_reconstruction_geometry | 0.337289984 | 0.53584253 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_geom_nearest_gap_weighted | blind_posterior_reconstruction_geometry | 0.337289984 | 0.53584253 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_margin_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.337289984 | 0.53584253 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_pmax_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.337289984 | 0.53584253 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_chrom_homolog_vector_dot | blind_reconstruction_geometry | 0.290060218 | 0.511171923 | -0.0472297659 | 0 | 0 | 0 | raw_reconstruction_chrom_centers_plus_raw_trans_pair_set |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_chrom_center_same_closer | blind_reconstruction_geometry | 0.285330296 | 0.502205213 | -0.051959688 | 0 | 0 | 0 | raw_reconstruction_chrom_centers_plus_raw_trans_pair_set |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_same_cross_distance_sum | blind_posterior_reconstruction_geometry | 0.284906617 | 0.494245609 | -0.0523833668 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | blind_posterior_mass_same_cross | blind_posterior_reconstruction_geometry | 0.281357439 | 0.482132563 | -0.0559325448 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | pair_independent_4state_oracle | snp_truth_oracle | 0.444571008 | 0.630496538 | 0.108878501 | 1 | 1 | 1 | snp_truth_eval_contacts |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.335692506 | 0.537169131 | 0 | 0 | 0 | 0 | none_cis_snp_gauge_only |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_geom_nearest_to_posterior_top | blind_posterior_reconstruction_geometry | 0.335692506 | 0.537169131 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_geom_nearest_gap_weighted | blind_posterior_reconstruction_geometry | 0.335692506 | 0.537169131 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_margin_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.335692506 | 0.537169131 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_pmax_geom_gap_weighted | blind_posterior_reconstruction_geometry | 0.335692506 | 0.537169131 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_chrom_homolog_vector_dot | blind_reconstruction_geometry | 0.298450447 | 0.508879891 | -0.0372420595 | 0 | 0 | 0 | raw_reconstruction_chrom_centers_plus_raw_trans_pair_set |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_chrom_center_same_closer | blind_reconstruction_geometry | 0.292227231 | 0.500246567 | -0.0434652757 | 0 | 0 | 0 | raw_reconstruction_chrom_centers_plus_raw_trans_pair_set |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_same_cross_distance_sum | blind_posterior_reconstruction_geometry | 0.282697931 | 0.490703376 | -0.0529945755 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | blind_posterior_mass_same_cross | blind_posterior_reconstruction_geometry | 0.264021337 | 0.461122263 | -0.0716711697 | 0 | 0 | 0 | raw_posterior_trans_bpair |

## Interpretation

If a blind decoder reaches +0.1, the next step is to turn that decoder into a controlled training or post-processing candidate. If only the SNP oracle improves, the current posterior/geometry does not contain enough blind information to select the needed trans copy-state gauge.
