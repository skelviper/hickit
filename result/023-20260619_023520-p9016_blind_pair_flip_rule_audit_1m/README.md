# 023 P9016 Blind Pair Flip Rule Audit 1Mb

This is an eval-stage diagnostic for the current trans failure. It asks whether the +0.1 trans gain seen in pair-specific SNP oracle space can be approximated by blind rules that choose a relative chromosome-pair flip from reconstruction geometry and posterior confidence only.

## Boundary

- Training is not run in this experiment.
- The two input structures come from 022 outputs.
- SNP labels are used only to score selected rules after the flip rule is chosen.
- CHARM/3DG is used only to keep the same standard eval denominator.
- `pair_independent_trans_oracle` is eval-only and uses SNP truth to select flips.
- Policies starting with `blind_` do not use SNP/CHARM to select flips.

## Main Result

The blind rules do not improve over cis-selected whole-chrom gauge. For both tested configs, all geometry/posterior blind flip policies give the same trans top1 as cis-selected baseline, while the eval-only SNP oracle is higher.

| source config | cis-selected trans | eval-only pair oracle | best blind rule |
|---|---:|---:|---:|
| p9016_pcgamma1_baseline_sep_off_copytrack0_prior0 | 0.315370 | 0.391097 | 0.315370 |
| p9016_pcgamma1_msep1_lsep0p5_copytrack0p03_prior0 | 0.342332 | 0.395403 | 0.342332 |

## Interpretation

This means the available reconstruction/posterior geometry does not contain a blind-selectable chromosome-pair flip rule that recovers the eval-only pair oracle. Combined with 021/022, this argues against continuing to strengthen posterior-derived pair priors: the prior would amplify the current posterior/geometry choice, which is not aligned with measured SNP truth.

The practical next step is not another average/sharpened pair prior. To get a real +0.1 under blind training, we would need a new source of blind information, such as true higher-order/multiway contacts, split-read group structure, or a formally different objective that can validate pair-level copy orientation without SNP truth.
