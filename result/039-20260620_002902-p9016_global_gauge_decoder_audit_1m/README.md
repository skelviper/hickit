# 039 P9016 Global Gauge Decoder Audit

This is a post-training diagnostic. It does not retrain Hickit. It tests whether raw posterior/reconstruction signals can choose one chromosome copy-swap bit per chromosome and improve full-denominator trans top1.

- full result root: `/mnt/ssd/zliu/phase3/test_res/039-20260620_002902-p9016_global_gauge_decoder_audit_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/039-20260620_002902-p9016_global_gauge_decoder_audit_1m`
- headline: `NO_RAW_GLOBAL_GAUGE_PLUS_0P1`
- training boundary: no training is run here; source runs are required to have `input_contact_source=raw_pairs`, `uses_phase_labels=0`, and `uses_charm_or_reference=0`.
- base gauge: `cis_selected_whole_chrom` uses the standard eval-only whole-chromosome SNP cis oracle gauge, matching the existing evaluator headline.
- blind decoder boundary: `raw_global_*` policies select chromosome swaps using raw posterior trans bpair weights (`n_raw`) and reconstruction geometry/posterior fields only.
- scoring boundary: SNP labels and CHARM/3DG are used only after decoder selection for scoring and shared denominator matching.
- oracle decoder policy uses SNP truth for global chromosome gauge selection and is an eval-only ceiling.

## Main Table

| config_name | policy | flip_source | top1_accuracy | same_cross_accuracy | delta_top1_vs_cis_selected | target_plus_0p1_met | uses_snp_for_global_gauge_selection | uses_snp_labeled_denominator_for_selection_weight | decoder_selection_denominator |
|---|---|---|---|---|---|---|---|---|---|
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | global_chrom_snp_oracle | snp_truth_global_chrom_oracle | 0.355980469 | 0.551789522 | 0.00297269703 | 0 | 1 | 1 | snp_truth_eval_contacts |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.353007772 | 0.552720226 | 0 | 0 | 0 | 0 | none_cis_snp_gauge_only |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_global_geom_nearest_to_posterior_top | blind_raw_posterior_reconstruction_global_gauge | 0.353007772 | 0.552720226 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_global_geom_nearest_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.353007772 | 0.552720226 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_global_margin_geom_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.353007772 | 0.552720226 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_global_pmax_geom_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.353007772 | 0.552720226 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_global_same_cross_distance_sum | blind_raw_posterior_reconstruction_global_gauge | 0.272425457 | 0.491647972 | -0.0805823152 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_global_posterior_same_cross_mass | blind_raw_posterior_reconstruction_global_gauge | 0.26581329 | 0.489355939 | -0.0871944825 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_global_chrom_center_same_closer | blind_raw_posterior_reconstruction_global_gauge | 0.248428568 | 0.483459164 | -0.104579204 | 0 | 0 | 0 | raw_reconstruction_chrom_centers |
| p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035 | raw_global_homolog_vector_dot | blind_raw_posterior_reconstruction_global_gauge | 0.248428568 | 0.483459164 | -0.104579204 | 0 | 0 | 0 | raw_reconstruction_chrom_centers |
| p9016_pcgamma1_td0p5_mstep0_sep0 | global_chrom_snp_oracle | snp_truth_global_chrom_oracle | 0.354098224 | 0.547851393 | 0.00972377533 | 0 | 1 | 1 | snp_truth_eval_contacts |
| p9016_pcgamma1_td0p5_mstep0_sep0 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.344374449 | 0.538391549 | 0 | 0 | 0 | 0 | none_cis_snp_gauge_only |
| p9016_pcgamma1_td0p5_mstep0_sep0 | raw_global_geom_nearest_to_posterior_top | blind_raw_posterior_reconstruction_global_gauge | 0.344374449 | 0.538391549 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | raw_global_geom_nearest_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.344374449 | 0.538391549 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | raw_global_margin_geom_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.344374449 | 0.538391549 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | raw_global_pmax_geom_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.344374449 | 0.538391549 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | raw_global_same_cross_distance_sum | blind_raw_posterior_reconstruction_global_gauge | 0.277558221 | 0.500420206 | -0.0668162276 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | raw_global_posterior_same_cross_mass | blind_raw_posterior_reconstruction_global_gauge | 0.27007786 | 0.498044827 | -0.074296589 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0_sep0 | raw_global_chrom_center_same_closer | blind_raw_posterior_reconstruction_global_gauge | 0.254964335 | 0.492703696 | -0.0894101141 | 0 | 0 | 0 | raw_reconstruction_chrom_centers |
| p9016_pcgamma1_td0p5_mstep0_sep0 | raw_global_homolog_vector_dot | blind_raw_posterior_reconstruction_global_gauge | 0.254964335 | 0.492703696 | -0.0894101141 | 0 | 0 | 0 | raw_reconstruction_chrom_centers |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | global_chrom_snp_oracle | snp_truth_global_chrom_oracle | 0.344812019 | 0.54574689 | 0.00752203477 | 0 | 1 | 1 | snp_truth_eval_contacts |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.337289984 | 0.53584253 | 0 | 0 | 0 | 0 | none_cis_snp_gauge_only |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | raw_global_geom_nearest_to_posterior_top | blind_raw_posterior_reconstruction_global_gauge | 0.337289984 | 0.53584253 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | raw_global_geom_nearest_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.337289984 | 0.53584253 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | raw_global_margin_geom_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.337289984 | 0.53584253 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | raw_global_pmax_geom_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.337289984 | 0.53584253 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | raw_global_same_cross_distance_sum | blind_raw_posterior_reconstruction_global_gauge | 0.240837078 | 0.48441765 | -0.0964529057 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | raw_global_chrom_center_same_closer | blind_raw_posterior_reconstruction_global_gauge | 0.240837078 | 0.48441765 | -0.0964529057 | 0 | 0 | 0 | raw_reconstruction_chrom_centers |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | raw_global_homolog_vector_dot | blind_raw_posterior_reconstruction_global_gauge | 0.240837078 | 0.48441765 | -0.0964529057 | 0 | 0 | 0 | raw_reconstruction_chrom_centers |
| p9016_pcgamma1_td0p5_mstep0p5_p4_w10_sep0 | raw_global_posterior_same_cross_mass | blind_raw_posterior_reconstruction_global_gauge | 0.230203435 | 0.476305243 | -0.107086549 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | global_chrom_snp_oracle | snp_truth_global_chrom_oracle | 0.347187398 | 0.564694361 | 0.0114948915 | 0 | 1 | 1 | snp_truth_eval_contacts |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | cis_selected_whole_chrom | cis_snp_truth_whole_chrom | 0.335692506 | 0.537169131 | 0 | 0 | 0 | 0 | none_cis_snp_gauge_only |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | raw_global_geom_nearest_to_posterior_top | blind_raw_posterior_reconstruction_global_gauge | 0.335692506 | 0.537169131 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | raw_global_geom_nearest_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.335692506 | 0.537169131 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | raw_global_margin_geom_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.335692506 | 0.537169131 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | raw_global_pmax_geom_gap_weighted | blind_raw_posterior_reconstruction_global_gauge | 0.335692506 | 0.537169131 | 0 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | raw_global_same_cross_distance_sum | blind_raw_posterior_reconstruction_global_gauge | 0.263972718 | 0.512338776 | -0.0717197886 | 0 | 0 | 0 | raw_posterior_trans_bpair |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | raw_global_chrom_center_same_closer | blind_raw_posterior_reconstruction_global_gauge | 0.263972718 | 0.512338776 | -0.0717197886 | 0 | 0 | 0 | raw_reconstruction_chrom_centers |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | raw_global_homolog_vector_dot | blind_raw_posterior_reconstruction_global_gauge | 0.263972718 | 0.512338776 | -0.0717197886 | 0 | 0 | 0 | raw_reconstruction_chrom_centers |
| p9016_pcgamma1_td0p5_mstep1_p1_w0_sep0 | raw_global_posterior_same_cross_mass | blind_raw_posterior_reconstruction_global_gauge | 0.221785424 | 0.464768678 | -0.113907082 | 0 | 0 | 0 | raw_posterior_trans_bpair |

## Interpretation

Current results show that even the eval-only `global_chrom_snp_oracle` is far below +0.1 over the cis-selected baseline, so a single whole-chromosome copy-swap gauge cannot explain the trans deficit in these source configs. This does not rule out finer chromosome-pair, local, or contact-level gauge errors, but it argues against whole-chrom global gauge as the next training mechanism.
