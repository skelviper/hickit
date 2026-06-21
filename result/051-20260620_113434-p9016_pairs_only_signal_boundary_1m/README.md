# 051 P9016 Pairs-Only Signal Boundary 1Mb

This report consolidates 024/028/029/038/039/046/049/050 evidence to decide whether another strict P9016.pairs.gz-only knob sweep is likely to deliver +0.1 full-denominator trans top1.

- full result root: `/mnt/ssd/zliu/phase3/test_res/051-20260620_113434-p9016_pairs_only_signal_boundary_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/051-20260620_113434-p9016_pairs_only_signal_boundary_1m`
- current best strict pairs-only trans: `0.358585`
- +0.1 target: `0.458585`
- pairs-only read linkage available: `0`
- pairs-only strand information available: `0`
- recommended next step: `stop_pairs_only_knob_sweeps; request boundary change to upstream read-level source or report callable/abstention only`

## Evidence Table

| evidence_name | category | strict_pairs_only_candidate | trans_top1 | target_plus_0p1 | meets_plus_0p1 | evidence_path | note |
|---|---|---|---|---|---|---|---|
| current_best_047_init_anchor | strict_pairs_only_model | 1 | 0.358585 | 0.458585 | 0 | result/047-20260620_082138-p9016_init_coord_anchor_1m/summary.tsv | Current best strict pairs-only standard trans result used as active target reference. |
| 049_trans_dscale_gamma_best | strict_pairs_only_model | 1 | 0.355661 | 0.458585 | 0 | result/049-20260620_103248-p9016_trans_dscale_gamma_1m/summary.tsv | Trans-specific dscale gamma failed to exceed current best. |
| 050_pairs_only_identifiability_baseline | strict_pairs_only_audit | 1 | 0.355660974 | 0.458585 | 0 | result/050-20260620_113123-p9016_pairs_identifiability_audit_1m/summary.tsv | Approved P9016 pairs have no read linkage or strand variation; high callable accuracy has low recall. |
| 050_best_blind_callable_subset_le50pct | strict_pairs_only_callable_subset | 0 | 0.613148789 | 0.458585 | 1 | result/050-20260620_113123-p9016_pairs_identifiability_audit_1m/diagnostics/pairs_identifiability/coverage_accuracy_curve.tsv | High top1 on selected subset only; recall is not full denominator model performance. |
| 038_pair_independent_snp_oracle | eval_only_oracle | 0 | 0.446501872 | 0.458585 | 0 | result/038-20260620_000849-p9016_trans_decoder_audit_1m/summary.tsv | SNP oracle decoder improves only about 0.065 over its cis-selected base and still misses +0.1 over current best. |
| 039_global_chrom_snp_oracle | eval_only_oracle | 0 | 0.355980469 | 0.458585 | 0 | result/039-20260620_002902-p9016_global_gauge_decoder_audit_1m/summary.tsv | Whole chromosome gauge oracle gives only small improvement. |
| 029_contacts_seg_pairwise | source_boundary_violation | 0 | 0.353975 | 0.41537 | 0 | result/029-20260619_061138-p9016_readgroup_joint_marginal_1m/summary.tsv | Uses contacts.seg-derived pairs, not approved P9016 pairs; improves but not +0.1. |
| 029_contacts_seg_readgroup_joint | source_boundary_violation | 0 | 0.333725 | 0.41537 | 0 | result/029-20260619_061138-p9016_readgroup_joint_marginal_1m/summary.tsv | Uses contacts.seg read grouping and is worse than seg pairwise. |

## Interpretation

The strict candidate space has not produced a +0.1 full-denominator trans gain. The approved P9016 pairs file has no usable readID/read-group linkage and no strand variation; phase columns are eval-only. High-accuracy trans calls exist only as low-coverage callable subsets. Contacts.seg-derived experiments are useful source-boundary controls but are not valid strict pairs-only candidates. Continuing EM knob sweeps inside the same pairs-only signal surface is therefore unlikely to meet the requested +0.1 goal.
