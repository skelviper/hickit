# 041 P9016 Callable Trans Abstention Audit

This is a diagnostic run, not a new training experiment. It asks whether raw-blind confidence scores can select a trans subset whose four-state top1 accuracy exceeds the all-contact trans denominator by at least 0.1.

## Boundary

- Call scores use only posterior probabilities, reconstruction geometry, and raw contact count.
- Call scores do not use SNP labels, phase labels, CHARM/3DG, or reference structure.
- SNP labels are used only after scoring to compute accuracy under the standard whole-chromosome cis SNP gauge.
- CHARM/3DG is used only as the same shared-denominator filter used by the standard evaluator.

## Headline

- `CALLABLE_SUBSET_PLUS_0P1`
- best config: `p9016_pcgamma1_td0p5_normdir0p1_lct0p03_gct0p003`
- all-trans top1: `0.340943345`
- best called top1 <=50% coverage: `0.748611111`
- best delta: `0.407667766`
- best called fraction: `0.0100015975`

## Predeclared Policy

- primary score: `neg_entropy`
- fixed coverage levels: `0.01, 0.02, 0.05, 0.1, 0.2, 0.3`

| coverage | top1 | delta vs all-trans | same/cross | called contacts | called bin pairs | chrom pairs | top chrom-pair fraction |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | 0.740484429 | 0.387566949 | 0.776470588 | 1445 | 125 | 32 | 0.142560554 |
| 0.02 | 0.692839848 | 0.339922368 | 0.755447942 | 2891 | 277 | 54 | 0.0927014874 |
| 0.05 | 0.663286568 | 0.310369088 | 0.736491179 | 7199 | 942 | 86 | 0.0670926518 |
| 0.1 | 0.620320856 | 0.267403376 | 0.701923745 | 14399 | 2679 | 114 | 0.0467393569 |
| 0.2 | 0.544103348 | 0.191185868 | 0.657070426 | 28796 | 11014 | 182 | 0.0417766356 |
| 0.3 | 0.484187619 | 0.131270139 | 0.623165254 | 43194 | 22757 | 187 | 0.0340556559 |


## Final Scope Verdict

- Full-denominator trans top1 is not improved by +0.1 in this run; this is not a full trans rescue.
- Predeclared blind abstention policy does improve the called trans subset by +0.1: primary score `neg_entropy`, fixed 10% coverage, best035 replay trans top1 `0.620320856` versus all-trans `0.35291748`.
- The called 10% subset covers 114 chromosome pairs and the largest chromosome pair contributes 4.67% of called contacts, so the signal is not from a single chromosome pair.
- SNP labels are used only for scoring after call selection; CHARM/3DG is used only for the same shared-denominator filter as the standard evaluator.

## Top Rows

| config | all trans top1 | best score | called frac | called top1 | delta | top1 at 10% | top1 at 20% | top1 at 30% |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| `p9016_pcgamma1_td0p5_normdir0p1_lct0p03_gct0p003` | 0.340943345 | distance_gap | 0.0100015975 | 0.748611111 | 0.407667766 | 0.566219876 | 0.516113349 | 0.455006366 |
| `p9016_pcgamma1_td0p5_lct0p03_gct0p003_best035_replay` | 0.35291748 | neg_entropy | 0.0100363252 | 0.740484429 | 0.387566949 | 0.620320856 | 0.546256425 | 0.484187619 |
| `p9016_pcgamma1_td0p5_normdir0p03_only` | 0.33931114 | pmax | 0.0100154886 | 0.710818308 | 0.371507168 | 0.580347222 | 0.509376302 | 0.454345511 |
| `p9016_pcgamma1_td0p5_normdir0p01_lct0p03_gct0p003` | 0.352820242 | distance_gap | 0.0100154886 | 0.684466019 | 0.331645777 | 0.612121633 | 0.550111127 | 0.477635783 |
| `p9016_pcgamma1_td0p5_normdir0p003_lct0p03_gct0p003` | 0.350111476 | neg_entropy | 0.0200170861 | 0.67869535 | 0.328583874 | 0.602402945 | 0.537505209 | 0.47080613 |
| `p9016_pcgamma1_td0p5_normdir0p03_lct0p03_gct0p003` | 0.353070282 | neg_entropy | 0.0200101405 | 0.65706352 | 0.303993238 | 0.608556744 | 0.549972218 | 0.477658934 |
| `p9016_pcgamma1_td0p5_normdir0p03_msep1_lsep0p5_lct0p03_gct0p003` | 0.339026372 | neg_entropy | 0.010008543 | 0.598195697 | 0.259169325 | 0.553495799 | 0.505417419 | 0.448279854 |

## Interpretation

- A blind call/no-call objective can produce a trans subset above all-contact trans top1 by more than 0.1.
- This does not solve the original full-denominator trans identity problem; it identifies a reliable subset.
- The next trainable/reporting step should expose a calibrated trans callability score and predeclare coverage levels.
- Because 024 found no usable readID/molecule-group anchor in the current pairs file, full-denominator +0.1 likely needs a new raw data source or external/reference information.
