# 052 P9016 Manifest-Strict Boundary Audit 1Mb

This diagnostic audits every numbered `test_res` summary row with a manifest-aware strict P9016 pairs-only boundary.

## Problem

A summary-only leaderboard can misclassify positive controls as valid blind candidates when the summary table does not carry all boundary fields. This matters because several reference-derived controls have high trans accuracy and can look close to the +0.1 goal unless their manifests are checked.

Small conceptual example: a row can report `model_top1_accuracy_genome_trans=0.455`, but its manifest may say `input_contact_source=charm3dg_derived_pairs` and `uses_charm_for_training=1`. That row is useful as a positive control, but it is not a strict P9016.pairs.gz-only training result.

Validation: this audit joins each summary row to `outputs/**/p9016_full.manifest.tsv`, then requires raw P9016 pairs, sample P9016, no training phase/reference flags, no fixed posterior, and an approved P9016 pairs path.

## Headline

- baseline trans used for +0.1 target: `0.358585`
- target trans: `0.458585`
- audited rows: `418`
- strict pairs-only rows: `308`
- non-strict/control rows: `110`
- best strict config: `p9016_initanchor_k0p001_seed71_noise0p05_best035`
- best strict run: `047-20260620_082138-p9016_init_coord_anchor_1m`
- best strict trans: `0.358585`
- best strict delta vs baseline: `0`
- reaches +0.1 target: `0`
- best excluded control config: `p9016_charm3dg20k_1pr_pcgamma1_common_ieps1_noise0p05`
- best excluded control trans: `0.455277`
- exclusion reason: `input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1`

## Top Strict Rows

| trans | config | run | cis | same/cross trans | cis Spearman |
|---:|---|---|---:|---:|---:|
| 0.358585 | `p9016_initanchor_k0p001_seed71_noise0p05_best035` | `047-20260620_082138-p9016_init_coord_anchor_1m` | 0.574695 | 0.563729 | 0.744353 |
| 0.357981 | `p9016_initanchor_k0p003_seed71_noise0p05_best035` | `047-20260620_082138-p9016_init_coord_anchor_1m` | 0.572327 | 0.561902 | 0.745886 |
| 0.355918 | `p9016_initanchor_k0p01_seed71_noise0p05_best035` | `047-20260620_082138-p9016_init_coord_anchor_1m` | 0.572191 | 0.561874 | 0.749965 |
| 0.355793 | `p9016_callable_negent_t0p10_m0p05_w10_seed71_best035` | `048-20260620_092323-p9016_callable_anchor_mstep_1m` | 0.573128 | 0.561437 | 0.747643 |
| 0.355737 | `p9016_callable_negent_t0p02_m0p05_w10_seed71_best035` | `048-20260620_092323-p9016_callable_anchor_mstep_1m` | 0.572643 | 0.561159 | 0.747986 |
| 0.355668 | `p9016_callable_negent_t0p10_m0p10_w10_seed71_best035` | `048-20260620_092323-p9016_callable_anchor_mstep_1m` | 0.572725 | 0.561152 | 0.747764 |
| 0.355661 | `p9016_initanchor_off_seed71_noise0p05_best035` | `047-20260620_082138-p9016_init_coord_anchor_1m` | 0.572065 | 0.56111 | 0.748189 |
| 0.355661 | `p9016_callable_anchor_off_seed71_noise0p05_best035` | `048-20260620_092323-p9016_callable_anchor_mstep_1m` | 0.572482 | 0.561208 | 0.748009 |
| 0.355661 | `p9016_transgamma_g1_td0p5_seed71_best035_replay` | `049-20260620_103248-p9016_trans_dscale_gamma_1m` | 0.572482 | 0.561208 | 0.748009 |
| 0.355647 | `p9016_callable_negent_t0p01_m0p10_w10_seed71_best035` | `048-20260620_092323-p9016_callable_anchor_mstep_1m` | 0.573147 | 0.561312 | 0.748077 |
| 0.355647 | `p9016_callable_negent_t0p05_m0p05_w10_seed71_best035` | `048-20260620_092323-p9016_callable_anchor_mstep_1m` | 0.572371 | 0.561201 | 0.747949 |
| 0.35564 | `p9016_multiseed_best035_seed71_noise0p05` | `044-20260620_040939-p9016_multiseed_basin_audit_1m` | 0.572808 | 0.561305 | 0.748036 |

## Top Excluded Controls

| trans | config | run | reason |
|---:|---|---|---|
| 0.455277 | `p9016_charm3dg20k_1pr_pcgamma1_common_ieps1_noise0p05` | `016-20260618_150851-p9016_charm_contacts_condition_sweep_1m` | input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1 |
| 0.455087 | `p9016_charm3dg20k_1pr_expected_common_msep1p5_lsep0p5` | `016-20260618_150851-p9016_charm_contacts_condition_sweep_1m` | input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1 |
| 0.455087 | `p9016_charm3dg20k_1pr_pcgamma1_common_msep1p5_lsep1` | `016-20260618_150851-p9016_charm_contacts_condition_sweep_1m` | input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1 |
| 0.455087 | `p9016_fixed_pc_same_cross_gated` | `018-20260618_214516-p9016_fixed_posterior_graph_semantics_1m` | input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1;fixed_posterior=1 |
| 0.455087 | `p9016_fixed_pc_softall_replay` | `018-20260618_214516-p9016_fixed_posterior_graph_semantics_1m` | input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1;fixed_posterior=1 |
| 0.455087 | `p9016_fixed_pc_softall_trans_dscale0p5` | `018-20260618_214516-p9016_fixed_posterior_graph_semantics_1m` | input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1;fixed_posterior=1 |
| 0.455087 | `p9016_fixed_pc_softall_trans_dscale0p75` | `018-20260618_214516-p9016_fixed_posterior_graph_semantics_1m` | input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1;fixed_posterior=1 |
| 0.455087 | `p9016_fixed_pc_top1_only` | `018-20260618_214516-p9016_fixed_posterior_graph_semantics_1m` | input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1;fixed_posterior=1 |
| 0.455087 | `p9016_fixed_pc_top2_only` | `018-20260618_214516-p9016_fixed_posterior_graph_semantics_1m` | input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1;fixed_posterior=1 |
| 0.455087 | `source_016_before_018_relax` | `018-20260618_220644-p9016_fixed_posterior_graph_semantics_1m` | manifest_missing;sample=missing;input_contact_source=missing |
| 0.455087 | `p9016_fixed_state_same_cross_gated` | `018-20260618_220644-p9016_fixed_posterior_graph_semantics_1m` | input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1;fixed_posterior=1 |
| 0.455087 | `p9016_fixed_state_softall` | `018-20260618_220644-p9016_fixed_posterior_graph_semantics_1m` | input_contact_source=charm3dg_derived_pairs;uses_charm_or_reference=1;uses_charm_for_training=1;fixed_posterior=1 |

## Interpretation

The high 018 rows are positive controls, not blind candidates: their manifests trace back to CHARM/3DG-derived contacts and fixed posteriors. Under the strict P9016.pairs.gz-only boundary, the best verified full-denominator trans result remains below the requested +0.1 target.

The useful scientific lesson is that reference-derived clean contacts can nearly reach the target, but the approved pairs-only training surface has not exposed enough blind signal to reproduce that behavior. The next validation step is either a boundary change to an upstream read-level source, or a different objective that reports callable/abstention accuracy and recall.
