# Experiment 057: ORACLE decomposition

> Formal alias 119 is the corrected result freeze. Historical alias 116 is
> retained only as superseded provenance because its synthetic near-truth
> contact metrics did not separate raw from geometry-gauge-aligned values.

## Status and scope

Logical experiment 057 is complete at 1 Mb for P9016 and for the explicitly
named synthetic near-truth control. Its corrected formal result-folder alias is
experiment 119 because prior sequence numbers were already occupied. The
machine-readable run contract is in
`results/phase_diagnostics/057/run_manifest.tsv`.

This experiment does not produce a new blind P9016 model. It asks four narrower
questions:

1. Given fixed CHARM/3DG coordinates, can the current E-step identify the
   measured contact state?
2. Given measured contact states, can native Hickit recover unordered homolog
   geometry?
3. Given an initialization near synthetic truth, does the joint loop remain in
   or return to that basin?
4. How much reported accuracy is explained by a prior or a randomized control?

## Truth boundary

P9016 phase and CHARM/3DG coordinates enter only rows marked
`ORACLE_COORDINATES` or `ORACLE_CONTACT_PHASE`. Synthetic truth enters
initialization only in rows marked `ORACLE_NEAR_TRUTH_INITIALIZATION`. The
uniform, configured cis-parity, unlabeled chromosome-pair, maximum-entropy, and
random-posterior baselines do not use truth in inference. No SNP label is used
to choose a blind seed, tune a score, construct a prior, select a chromosome
flip, or filter a training contact.

## 057A: fixed ORACLE coordinates and current E-step

For observation m and state s, the tested implementation computes

```text
z_ms = [-E(d_ms / (unit * d_scale_m)) + log prior_ms] / temperature
       + log_rate_weight_ms
p_ms = softmax_s(z_ms)
```

All six score modes implemented by the corrected C API were evaluated on the
same 350,034 contacts: `fdg_flat`, `dist2`, `logdist2`,
`monotonic_powerlaw`, `monotone_log`, and `monotone_logistic`. The Python
diagnostic reproduces their formulas with a fixed base coefficient because the
historical reference-coordinate neighbor coefficient was not serialized.

### Accuracy and proper scores

| score | scope | N | exact top1 | same/cross | exact log loss | same/cross log loss | exact ECE |
|---|---|---:|---:|---:|---:|---:|---:|
| FDG-flat | all | 350,034 | 0.198392 | 0.847369 | 3.374511 | 0.464089 | 0.533804 |
| FDG-flat | cis | 206,057 | 0.221094 | 0.977035 | 1.981686 | 0.180895 | 0.479043 |
| FDG-flat | trans | 143,977 | 0.165901 | 0.661793 | 5.367894 | 0.869390 | 0.612177 |
| DIST2 | all | 350,034 | 0.188702 | 0.847732 | 12.233863 | 1.874182 | 0.726571 |
| LOGDIST2 | all | 350,034 | 0.188702 | 0.847438 | 1.907065 | 0.294315 | 0.498156 |
| monotone power law | all | 350,034 | 0.188702 | 0.847398 | 1.690982 | 0.295519 | 0.451472 |
| monotone log | all | 350,034 | 0.188702 | 0.847455 | 1.879695 | 0.294083 | 0.493202 |
| monotone logistic | all | 350,034 | 0.188702 | 0.847789 | 8.434970 | 1.070444 | 0.716191 |

Uniform four-state prediction has exact top1 0.419999 and exact log loss
1.386294 because exact top1 tie-breaking always chooses state 00 and the truth
states are imbalanced. This is precisely why top1 alone is not an adequate
metric. The monotone score has worse exact log loss and accuracy than the
uniform distribution despite seeing reference geometry.

### Geometry copying

The truth state is the geometrically closest of the four candidates for only
18.8702% of all contacts and 16.5860% of trans contacts. Monotone-log posterior
top1 is the geometrically closest state for 100% of contacts. FDG-flat posterior
top1 is closest for 95.9024% overall and 99.8416% for trans. Therefore the
posterior is primarily a decoder of the current geometry, not an independent
recovery of measured exact phase.

FDG-flat is additionally unsuitable as a contact observation likelihood. It
penalizes distances below 0.5 layout units, is exactly flat from 0.5 to 1.5,
and then increases. Numerical curves at multiple temperatures are in
`estep_score_curves.tsv` and `estep_score_curves.pdf`. DIST2, LOGDIST2,
monotone power law, monotone log, and monotone logistic are finite at zero and
non-increasing as contact log scores; FDG-flat fails monotonicity.

### Calibration and selective calling

Exact-state callable coverage is zero at target precision 0.55, 0.60, 0.70,
0.80, and 0.90 for every main fixed-coordinate scorer. This result uses complete
confidence ties and cannot split equal-confidence calls to manufacture
coverage.

Monotone-log same/cross coverage at target precision 0.70, 0.80, and 0.90 is
1.000, 1.000, and 0.862 overall, but cis contacts dominate. For trans the
corresponding values are 0.803, 0.473, and 0.281. The trans same/cross ECE is
0.0637 for monotone-log versus 0.1875 for FDG-flat. These are evaluation curves,
not a blind-calibrated deployment threshold, and they do not imply globally
labeled exact state.

The reliability table is `calibration.tsv`. Accuracy/precision, coverage,
correct-recall, and minimum confidence on a fixed threshold grid are in
`selective_calling.tsv`.

## 057B: ORACLE contact phase and current M-step

The current source was run in three explicit arms for seeds 17, 29, and 43.
Every input contact state was fixed to measured phase before the M-step; these
are `ORACLE_CONTACT_PHASE` experiments, not blind models. The native arms use
the current CPU Hickit optimizer. The third arm starts from each corresponding
cis-only result, freezes its chromosome shapes, and optimizes rigid chromosome
placement from trans contacts in the isolated Python diagnostic.

| ORACLE-phase arm | unordered distance Spearman, mean +/- SD | Procrustes RMSD, mean +/- SD | homolog-separation Spearman | centroid-placement Spearman | runtime |
|---|---:|---:|---:|---:|---:|
| native Hickit cis only | 0.8024 +/- 0.0066 | 1.4082 +/- 0.0262 | 0.1686 | -0.0668 | 21.80 s |
| native Hickit cis + trans | 0.9661 +/- 0.0155 | 0.7533 +/- 0.3735 | 0.8792 | 0.9724 | 48.12 s |
| fixed cis shapes, rigid trans placement | 0.8024 +/- 0.0066 | 1.4082 +/- 0.0262 | 0.1686 | 0.1497 | 9.24 s placement only |

Joint cis+trans Hickit strongly reconstructs unordered geometry when phase is
known. Cis-only reconstruction recovers chromosome-internal distance order but
not homolog separation or interchromosomal placement. The rigid-body trans
stage does not recover placement from those frozen cis-only shapes: its mean
centroid-distance Spearman is only 0.150. Its end-to-end runtime is the cis-only
runtime plus 9.24 seconds, about 31.03 seconds on average. Thus the current
M-step is capable under ORACLE phase, but the proposed rigid staged placement
is not validated. Historical experiment 083 remains a hash-closed lineage
comparator and is not substituted for these three current arms.

## 057C: ORACLE near-truth initialization

The joint synthetic loop was initialized at truth plus Gaussian noise measured
in median adjacent-bin distance units. Each level used seeds 17, 29, and 43.

| noise | exact top1 mean, geometry-aligned | same/cross mean, geometry-aligned | unordered distance Spearman mean | global RMSD mean | converged fraction |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.3186 | 0.5510 | 0.9415 | 1.2740 | 0.667 |
| 0.05 | 0.3224 | 0.5504 | 0.9487 | 1.2357 | 0.667 |
| 0.10 | 0.3226 | 0.5549 | 0.9510 | 1.2033 | 0.667 |
| 0.25 | 0.3323 | 0.5697 | 0.9525 | 1.1573 | 0.667 |
| 0.50 | 0.3313 | 0.5850 | 0.9314 | 1.0795 | 0.667 |
| 1.00 | 0.3084 | 0.5643 | 0.8274 | 1.3632 | 0.667 |
| 2.00 | 0.2721 | 0.5428 | 0.5909 | 1.9106 | 0.000 |

Starting at truth does not yield high phase accuracy: evaluation-only
geometry-aligned exact top1 is 0.319 and same/cross is 0.551 at zero noise.
Structure stays near truth through moderate noise but degrades at noise 1-2,
and convergence remains 2/3 through noise 1 before falling to zero.
Random initialization is therefore a secondary optimization problem, not the
sole explanation: exact decoding is weak even in the truth basin.

## 057D: prior-only baselines

The uniform, unlabeled chromosome-pair, and unlabeled state-frequency baselines
are identical because no labeled training target can break the four-state
symmetry. Their apparent all-contact exact accuracy of 0.4200 and same/cross
accuracy of 0.7713 are majority/tie-breaking effects. The configured 0.9 cis
same prior improves log loss for the highly imbalanced cis set but uses no
geometry. A full model must beat these proper scores and stratified accuracies,
not merely the raw top1 number.

## 057E: randomized controls

Random geometry and random p4 produce approximately random exact and parity
accuracy (about 0.25 and 0.50). Stratum-preserving endpoint shuffling retains
all-contact exact 0.3622 and same/cross 0.8252 because it retains chromosome,
genomic-distance, and class-imbalance structure; it does not recover the true
contact state. Permuting chromosome labels for trans reduces trans same/cross
to 0.5012. These controls demonstrate that large cis parity numbers can survive
destruction of contact-specific information.

## Failure decomposition

| component | evidence | conclusion |
|---|---|---|
| E-step statistical rule | ORACLE coordinates give exact 0.189-0.198 and worse exact log loss than uniform | dominant exact-state failure; current geometry is not sufficient |
| E-step implementation | monotone modes are stable; FDG-flat is nonmonotone and flat | FDG-flat is misspecified as an observation likelihood |
| M-step | current ORACLE-phase joint cis+trans unordered Spearman 0.966 +/- 0.016 | geometry engine works well when contact state is known |
| staged trans placement | fixed-cis centroid-distance Spearman 0.150 +/- 0.067 | rigid trans placement is not sufficient in this implementation |
| optimization/basin | near-truth runs drift and convergence is seed/noise dependent | substantial secondary failure |
| chromosome gauge | trans label permutation removes parity signal; exact labels require declared alignment | necessary ambiguity, but not the only trans failure |
| score/data ambiguity | truth is nearest for only 16.6% of trans contacts | the tested nearest-distance decoder cannot identify exact trans state; this is not an algorithm-independent information bound |

## Diagnostic-supported correction

Concrete problem: the E-step imports collision repulsion and a flat shell from
the layout optimizer. This matters because a contact at distance 0.1 can score
worse than one at distance 0.5 even though the observation model has no capture
mechanism that predicts this inversion. For example, two state candidates at
0.2 and 1.0 receive different physical shell terms but both can lie in a
layout's acceptable region; a contact likelihood should rank 0.2 no lower.
The implemented correction exposes a monotone likelihood independently of the
M-step potential. Validation consists of numerical monotonicity/zero-distance
tests, ORACLE-coordinate log loss and ECE, and blocked held-out pseudo-NLL.

Concrete problem: posterior enters both the edge coefficient and the
state-specific target distance in the historical M-step. This matters because
a small early preference can be amplified more than once. For example, changing
p from 0.5 to 0.95 changes both `k*p` and a target proportional to
`p^(-1/3)`. The `posterior_once` mode keeps the target based on raw count and
uses posterior only in the coefficient. Validation is the per-observation mass
table and symmetry/unit tests, followed by the fixed comparison matrix.

## Required conclusions

- Is the E-step statistically informative under ORACLE coordinates? It is
  informative for selected same/cross subsets, especially cis, but not for
  exact four-state prediction. Exact proper scores are worse than uniform.
- Can the M-step recover geometry under ORACLE phase? Yes. Unordered homolog
  distance Spearman is 0.966 +/- 0.016 for current joint cis+trans Hickit.
  Cis-only is 0.802, while the fixed-shape trans-placement prototype fails to
  recover chromosome placement.
- Is random initialization the main problem? No. It is a secondary source of
  drift and non-convergence, but exact inference is weak even at truth.
- Does the full model outperform prior-only baselines? Not on exact
  ORACLE-coordinate P9016 metrics; the apparent prior accuracy is class
  imbalance, while the full score has worse exact log loss.
- Which failure is dominant? Exact-state likelihood/data mismatch is dominant
  in this decomposition, with joint optimization instability secondary. The
  M-step is not dominant when supplied ORACLE phase. Generator-conditional
  information limits are quantified separately in experiment 059.
