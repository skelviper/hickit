# 024 P9016 Pairs Source Signal Audit

This is a diagnostic run, not a training experiment. It asks whether the current P9016 `.pairs.gz` file retains blind read/molecule-group information that could anchor trans copy identity beyond single binned pairwise contacts.

## Paths

- full result root: `/mnt/ssd/zliu/phase3/test_res/024-20260619_024741-p9016_pairs_source_signal_audit_1m`
- diagnostics: `/mnt/ssd/zliu/phase3/test_res/024-20260619_024741-p9016_pairs_source_signal_audit_1m/diagnostics/pairs_source_signal`
- light result root: written under `hickit/result` by `scripts/p9016_publish_light_result.sh`

## Boundary

- Training uses no SNP, phase, CHARM/3DG, or reference information here because no training is run.
- Phase labels are used only for the eval-only split-half truth-stability audit.
- CHARM/3DG is not used.

## Headline

- `NO_READID_GROUP_ANCHOR`
- readID status: `UNUSABLE_ALL_READID_MISSING`
- total rows: `1703888`
- non-missing readID rows: `0`
- repeated readID groups: `0`
- trans singleton contact fraction at 1 Mb: `0.702155161`

## Split-Half Truth Stability

| scope | truth_n bucket | contact weight | 4-state agreement | same/cross agreement |
|---|---:|---:|---:|---:|
| trans | 11-20 | 5715 | 0.996150481 | 0.996150481 |
| trans | 2 | 9864 | 0.822992701 | 0.882400649 |
| trans | 21-50 | 1010 | 1 | 1 |
| trans | 3 | 3486 | 0.871772806 | 0.909638554 |
| trans | 4 | 2472 | 0.959546926 | 0.966019417 |
| trans | 5 | 1980 | 0.972222222 | 0.982323232 |
| trans | 6-10 | 7603 | 0.991319216 | 0.993949757 |

## Interpretation

- The current P9016 pairs file does not preserve usable readID or molecule-group identifiers. A blind multi-contact/read-group trans anchor cannot be built from this file alone.
- Most SNP-labeled trans binned pairs are singleton at 1 Mb, so full-contact four-state trans accuracy is heavily affected by low-count bpair strata.
- The next trainable route needs a new blind information source or a changed objective such as high-confidence call/no-call; more sep/dscale tuning alone is unlikely to give a stable +0.1 all-contact trans gain.
