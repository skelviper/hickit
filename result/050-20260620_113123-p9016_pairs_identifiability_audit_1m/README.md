# 050 P9016 Pairs Identifiability Audit 1Mb

This is a post-hoc diagnostic, not a new model. It audits whether fields visible in the approved P9016 pairs file and blind-model confidence scores contain enough signal to support a +0.1 full-denominator trans top1 improvement.

- full result root: `/mnt/ssd/zliu/phase3/test_res/050-20260620_113123-p9016_pairs_identifiability_audit_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/050-20260620_113123-p9016_pairs_identifiability_audit_1m`
- source config output: `/mnt/ssd/zliu/phase3/test_res/049-20260620_103248-p9016_trans_dscale_gamma_1m/outputs/p9016_transgamma_g1_td0p5_seed71_best035_replay/p9016_transgamma_g1_td0p5_seed71_best035_replay`
- baseline trans top1: `0.355660974`
- pairs-visible read linkage available: `0`
- pairs-visible strand information available: `0`
- best blind score <=50% coverage: `neg_entropy`
- best blind called fraction <=50% coverage: `0.0100363252`
- best blind top1 accuracy <=50% coverage: `0.613148789`
- best blind recall <=50% coverage: `0.00615376067`

## Top Blind Coverage Rows

| score | called_fraction | top1_accuracy | top1_recall | same_cross_accuracy | delta_top1_vs_baseline |
|---|---|---|---|---|---|
| neg_entropy | 0.0100363252 | 0.613148789 | 0.00615376067 | 0.72733564 | 0.257487815 |
| pmax | 0.0200309772 | 0.613037448 | 0.0122797391 | 0.731969487 | 0.257376474 |
| neg_entropy | 0.0200518138 | 0.610322134 | 0.0122380658 | 0.729130585 | 0.25466116 |
| neg_entropy | 0.0500566063 | 0.610240044 | 0.0305465456 | 0.71236298 | 0.254579071 |
| margin | 0.0200448683 | 0.609494109 | 0.0122172291 | 0.728343728 | 0.253833136 |
| pmax | 0.0500010418 | 0.60883456 | 0.0304423623 | 0.711904431 | 0.253173587 |
| margin | 0.0100363252 | 0.607612457 | 0.00609819624 | 0.716955017 | 0.251951483 |
| margin | 0.0501191162 | 0.606014412 | 0.0303729068 | 0.709534368 | 0.250353439 |
| margin | 0.100009029 | 0.601013959 | 0.0601068226 | 0.707549135 | 0.245352986 |
| pmax | 0.100015975 | 0.600625 | 0.0600720948 | 0.707430556 | 0.244964026 |

## Interpretation

The approved P9016 pairs file has phase columns, but those are eval-only and cannot enter blind training. The pairs-visible readID/strand fields provide no useful read-level linkage for training. Therefore, any remaining blind selection signal must come from raw count or model confidence. If high accuracy appears only at small coverage and recall stays low, that supports a callable/abstention interpretation rather than a full-denominator +0.1 model improvement.
