# 025 P9016 Callable Trans Confidence Audit

This is a diagnostic run, not a new training experiment. It asks whether blind confidence scores can select a trans subset with substantially higher four-state top1 accuracy than the all-contact trans denominator.

## Boundary

- Call scores use only posterior probabilities, reconstruction geometry, and raw contact count.
- Call scores do not use SNP labels, phase labels, CHARM/3DG, or reference structure.
- SNP labels are used only after scoring to compute accuracy under the standard whole-chromosome cis SNP gauge.
- CHARM/3DG is used only as the same shared-denominator filter used by the standard evaluator.

## Headline

- `CALLABLE_SUBSET_PLUS_0P1`
- best config: `p9016_pcgamma1_prior0p5_power1_warmup0_sep_off`
- all-trans top1: `0.340971127`
- best called top1 <=50% coverage: `0.595626519`
- best delta: `0.254655391`
- best called fraction: `0.0200101405`

## Top Rows

| config | all trans top1 | best score | called frac | called top1 | delta | top1 at 10% | top1 at 20% | top1 at 30% |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| `p9016_pcgamma1_prior0p5_power1_warmup0_sep_off` | 0.340971127 | pmax | 0.0200101405 | 0.595626519 | 0.254655391 | 0.561180556 | 0.506667593 | 0.445617447 |
| `p9016_pcgamma1_prior0p35_power4_warmup10_sep_off` | 0.329330379 | distance_gap | 0.0100849441 | 0.556473829 | 0.22714345 | 0.488405998 | 0.45475066 | 0.413390749 |
| `p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0` | 0.342332456 | margin | 0.0200518138 | 0.56356079 | 0.221228334 | 0.53992501 | 0.515696625 | 0.455827198 |
| `p9016_pcgamma0p5_prior0p5_power4_warmup10_sep_off` | 0.306389215 | margin | 0.0101196719 | 0.51132464 | 0.204935425 | 0.46673149 | 0.423704681 | 0.381595092 |
| `p9016_pcgamma1_baseline_sep_off_copytrack0_prior0` | 0.315369816 | neg_entropy | 0.0101543997 | 0.510259918 | 0.194890102 | 0.449652778 | 0.419849979 | 0.387290196 |
| `p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p5_power4_warmup10` | 0.317870215 | pmax | 0.0500079874 | 0.475833333 | 0.157963118 | 0.473819444 | 0.442179469 | 0.401643709 |
| `p9016_pcgamma1_prior0p5_power4_warmup20_sep_off` | 0.308236732 | distance_gap | 0.0100015975 | 0.465277778 | 0.157041046 | 0.424375 | 0.395888318 | 0.367603482 |
| `p9016_pcgamma1_prior0p5_power2_warmup10_sep_off` | 0.324537947 | n_raw | 0.100036811 | 0.462820246 | 0.138282299 | 0.462820246 | 0.450236144 | 0.40475529 |
| `p9016_pcgamma1_prior0p5_power8_warmup10_sep_off` | 0.300061815 | n_raw | 0.100036811 | 0.403318753 | 0.103256938 | 0.403318753 | 0.389081817 | 0.353914896 |
| `p9016_pcgamma1_prior0p65_power4_warmup10_sep_off` | 0.315460108 | neg_entropy | 0.0100363252 | 0.415224913 | 0.0997648053 | 0.41387211 | 0.407626059 | 0.374010279 |
| `p9016_pcgamma1_prior0p5_power4_warmup10_sep_off` | 0.296526529 | n_raw | 0.100036811 | 0.383253489 | 0.0867269603 | 0.383253489 | 0.373003195 | 0.342362365 |
| `p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0p5_power2_warmup10` | 0.299874285 | n_raw | 0.100036811 | 0.386169548 | 0.0862952625 | 0.386169548 | 0.375677177 | 0.348049543 |

## Interpretation

- A blind call/no-call objective can produce a trans subset above all-contact trans top1 by more than 0.1.
- This does not solve the original full-denominator trans identity problem; it identifies a reliable subset.
- The next trainable/reporting step should expose a calibrated trans callability score and predeclare coverage levels.
- Because 024 found no usable readID/molecule-group anchor in the current pairs file, full-denominator +0.1 likely needs a new raw data source or external/reference information.
