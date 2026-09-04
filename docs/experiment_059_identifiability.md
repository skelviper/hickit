# Experiment 059: synthetic identifiability phase diagram

> **Corrected diagnostic freeze.** Formal alias 121 contains the post-audit
> results. Historical alias 118 is retained only as superseded provenance because
> it used fixed total anchor weight, unpaired anchor data, depth-dependent prior
> strength, raw-gauge headline metrics, and local pair rather than joint molecule
> marginals.

## Status and scope

Logical experiment 059 is complete as a corrected synthetic diagnostic. Its
formal result-folder alias is 121. The full run contains 98 conditions,
including eight targeted transition-refinement conditions, with three random
replicates per condition. The four scientific slices are not a full factorial sweep.
Every condition, seed, generator family, likelihood-match status, and truth
boundary is recorded in `results/phase_diagnostics/059/conditions.tsv`.
Values below are three-replicate means. The corresponding sample standard
deviations are retained in every `*_sd` column of `replicate_summary.tsv`.

The nominal experiment resolution is 1,000,000 bp (1 Mb), with three
chromosomes and ten bins per chromosome. This small system is intended to
separate implementation, optimization, gauge, and observation limits; it is
not a quantitative surrogate for a full P9016 cell.

## Generator families

The well-specified family samples endpoint states with weight

```text
L_s(d_s) = (d_s + 0.1)^(-2)
```

and inference uses the matching normalized monotone-log endpoint likelihood.
The model-mismatched family uses saturated logistic capture, negative-binomial
count variation where configured, and random-ligation background. It does not
reuse the inference score.

Both families generate polymer midpoint trajectories and smooth homolog
differences with the symmetric representation

```text
X_i^0 = M_i + Delta_i
X_i^1 = M_i - Delta_i.
```

Multi-segment molecules are sampled from spatial neighborhoods. Pair expansion
is retained as a comparator. The molecule-level mode enumerates all `2^m`
shared copy assignments for the tested molecule lengths and divides pair terms
by `choose(m,2)`, giving every source molecule unit total mass. Reported pair
posteriors are marginals of the same shared-assignment posterior used for
training, not independent pairwise softmax values.

## Truth boundary

Synthetic truth generates observations and evaluates results. It does not
enter zero-anchor initialization, loss, seed selection, or Z2 synchronization.
Every nonzero anchor condition is explicitly
`ORACLE_EXTERNAL_PHASE_ANCHORS`; every ORACLE-coordinate row is separately
`ORACLE_COORDINATES`. No truth-derived chromosome spin is used in a blind row.

## Controlled slices

The API supports homolog separation, contacts per cell, cis/trans fraction,
background fraction, resolution, chromosome count, chromosome length,
readchain length, molecule count, count overdispersion, and anchor fraction.
The executed grid is:

1. contacts 300, 900, 3,000, and 12,000 by separation 0.1, 0.75, and 2.0;
2. targeted refinements at contacts 1,500 and 6,000 and separations 1.25 and
   1.6;
3. trans fractions 0.1, 0.5, and 0.9 by background 0, 0.2, and 0.5;
4. contacts 500, 2,500, and 10,000 by anchor 0, 0.001, 0.005, 0.01, and 0.05;
5. readchain lengths 2, 4, and 6 by molecule counts 0, 20, and 80.

Resolution and chromosome dimensions remain configurable but were not swept.

## Contacts by homolog separation

At 12,000 contacts, the zero-anchor blind solver does not improve reliably as
homolog separation increases:

All contact metrics in this table use per-chromosome geometry-selected copy
swaps for evaluation only; the raw arbitrary-gauge values remain in the TSV.

| family | separation | exact, geometry-aligned | same/cross, geometry-aligned | cis exact | cis same/cross | trans exact | trans same/cross | unordered distance Spearman | converged fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| matched | 0.10 | 0.2607 | 0.5006 | 0.2689 | 0.5011 | 0.2450 | 0.4995 | 0.8368 | 1.000 |
| matched | 0.75 | 0.2952 | 0.5196 | 0.3014 | 0.5237 | 0.2834 | 0.5119 | 0.8702 | 1.000 |
| matched | 2.00 | 0.2615 | 0.5303 | 0.2746 | 0.5477 | 0.2371 | 0.4979 | 0.7637 | 1.000 |
| mismatched | 0.10 | 0.2572 | 0.5059 | 0.2599 | 0.5040 | 0.2524 | 0.5093 | 0.8842 | 1.000 |
| mismatched | 0.75 | 0.2701 | 0.4946 | 0.2725 | 0.4937 | 0.2658 | 0.4962 | 0.8675 | 1.000 |
| mismatched | 2.00 | 0.2588 | 0.4786 | 0.2647 | 0.4687 | 0.2479 | 0.4969 | 0.7694 | 1.000 |

The same runs can recover a high unordered distance-matrix correlation while
remaining at random contact-state accuracy. This is direct evidence that
unordered structure, parity, and exact state are different inference targets.
Because the matched high-depth regime also fails, there is an implementation
or optimization limitation, not model mismatch alone.

### ORACLE coordinates and generator-conditional decision limit

The corrected run separates a monotone-log inference scorer at truth geometry
from a generator-conditional posterior. The latter knows the synthetic
generator and therefore estimates the conditional Bayes decision limit for
one stochastic contact. It does not estimate a P9016 limit.

| family | separation | generator-conditional exact | generator-conditional same/cross | trans exact | trans same/cross |
|---|---:|---:|---:|---:|---:|
| matched | 0.10 | 0.2818 | 0.4997 | 0.2635 | 0.4891 |
| matched | 0.75 | 0.3627 | 0.5532 | 0.3054 | 0.5202 |
| matched | 2.00 | 0.4119 | 0.6422 | 0.3476 | 0.5305 |
| mismatched | 0.10 | 0.2689 | 0.5131 | 0.2548 | 0.5062 |
| mismatched | 0.75 | 0.3360 | 0.5686 | 0.3052 | 0.5211 |
| mismatched | 2.00 | 0.4169 | 0.7157 | 0.2593 | 0.5035 |

Matched generator-conditional and monotone-log rows agree by construction.
Across the mismatched grid, using the actual generator rather than monotone-log
improves mean exact and same/cross log loss by 0.0403 and 0.0275. In the
12,000-contact separation slice shown above, even truth coordinates and the
correct generator leave exact state below 0.42 and trans same/cross near
0.50-0.53. Thus the blind solver gap is real, but substantial per-contact
ambiguity is also present in these generator conditions.

## Trans fraction and background

Across the tested trans/background slice, geometry-aligned blind trans
same/cross remains near chance: 0.484-0.522 for the matched family and
0.430-0.508 under mismatch. Unordered structure Spearman ranges 0.779-0.937
when matched and 0.866-0.933 under mismatch, without a monotonic background
response. A higher background fraction can change calibration and fitted signal
fraction, but it does not create phase information. The optional fifth
component is useful as an outlier sink; it is not evidence for one of the four
homolog states.

Background metrics are reported separately from conditional four-state
metrics. This prevents a model from appearing well phased merely by assigning
large probability to random ligation.

## ORACLE anchor fraction

At 10,000 contacts:

| family | anchor | actual anchor N | raw exact | raw same/cross | geometry-aligned trans same/cross | gauge recovery |
|---|---:|---:|---:|---:|---:|---:|
| matched | 0 | 0 | 0.2248 | 0.5122 | 0.5074 | 0.8889 |
| matched | 0.001 | 10 | 0.2619 | 0.5388 | 0.5009 | 0.7778 |
| matched | 0.005 | 50 | 0.2856 | 0.5389 | 0.5067 | 0.7778 |
| matched | 0.010 | 100 | 0.3063 | 0.5403 | 0.5037 | 0.8889 |
| matched | 0.050 | 500 | 0.3786 | 0.5807 | 0.5076 | 0.7778 |
| mismatched | 0 | 0 | 0.2262 | 0.4808 | 0.4946 | 0.7778 |
| mismatched | 0.001 | 10 | 0.2366 | 0.4766 | 0.4995 | 0.8889 |
| mismatched | 0.005 | 50 | 0.2554 | 0.4829 | 0.4950 | 0.7778 |
| mismatched | 0.010 | 100 | 0.2654 | 0.4822 | 0.4975 | 0.7778 |
| mismatched | 0.050 | 500 | 0.2957 | 0.4747 | 0.4968 | 0.8889 |

All nonzero-anchor rows are ORACLE. Within each family/depth/replicate, the
geometry and observations are identical and anchor masks are nested prefixes;
anchored observations are excluded from evaluation. Anchor evidence increases
with actual anchor mass rather than receiving a fixed total weight. Five
percent improves raw exact accuracy, especially when matched, but does not
produce trans parity above chance. Gauge recovery is nonmonotonic. No small-
anchor phase transition is supported.

Nominal and realized fractions differ at low depth because anchor counts are
integer and every nonzero nominal condition receives at least one anchor. At
500 contacts the five actual counts are `0,1,2,5,25`, corresponding to realized
fractions `0,0.002,0.004,0.01,0.05`; at 2,500 contacts they are
`0,2,12,25,125`, or `0,0.0008,0.0048,0.01,0.05`. The corresponding realized
fractions at 10,000 contacts equal their nominal fractions.

## Readchain likelihood

For lengths 4 and 6 and nonzero molecule counts, molecule-level minus normalized
pair-expansion averages are:

| family | geometry-aligned exact delta | geometry-aligned same/cross delta | unordered structure Spearman delta | runtime ratio | molecule replicate failure fraction |
|---|---:|---:|---:|---:|---:|
| matched | +0.0160 | -0.0046 | +0.0127 | 5.20 | 0.667 |
| mismatched | +0.0007 | -0.0065 | +0.0089 | 4.42 | 0.583 |

The joint molecule posterior has no stable parity gain and only a small
structure gain in this prototype, costs 4.4-5.2 times more, and has a high
replicate-level optimizer failure fraction: 8/12 matched and 7/12 mismatched
molecule fits fail across the four length-4/6 conditions. The paired normalized
pair-expansion baseline fails in 7/12 replicates in each family, so small deltas
must not be overinterpreted. For the reported conditional four-state and
structure metrics, length-two molecule and pair modes are identical by
construction; background/signal-aware coverage and runtime are not identical.
Real P9016 readchain benefit cannot be measured because its pairs file has no
usable molecule identifier.

One finite-count exception is retained rather than imputed: the three
`model_mismatched, readchain_length=6, molecule_count=80` replicates contain 9,
85, and 0 trans signal contacts. Its 12 trans summary metrics therefore report
`n_finite_replicates=2`; all replicate counts remain explicit in the TSV.

## Z2 synchronization

The conservative comparison has no calibrated blind relative-gauge edges.
The Z2 solver therefore reports zero edges, one unresolved component per
chromosome set, unresolved fraction 1.0, and leaves coordinates/posteriors
unchanged. This is the correct no-evidence behavior. It does not demonstrate
that synchronization helps; it shows that spins must not be invented from
truth labels or uncalibrated pair-state preferences.

## Interpretation rules

1. Well-specified high-depth contact recovery fails despite convergence. This
   identifies an implementation/optimization gap.
2. Generator-conditional ORACLE-coordinate exact accuracy remains only
   0.282-0.417 in the high-depth slice, and trans parity remains 0.489-0.531.
   These synthetic generators therefore have substantial per-contact
   information limits even when geometry and the capture model are known.
3. Small ORACLE anchor fractions do not cause a trans phase transition. Gauge
   ambiguity is necessary but does not explain the entire trans failure.
4. Experiment 057 shows that current joint Hickit can recover geometry under
   ORACLE phase. The fully informed M-step is not the dominant failure.
5. Exact state is much less recoverable than unordered structure. Same/cross
   becomes informative under favorable ORACLE geometry, but not in the blind
   trans solver.
6. The defensible output is not globally labeled `00/01/10/11` phasing.

## Improvement directions and validation

Concrete problem: the blind solver recovers a distance matrix while its state
posterior stays random. This matters because optimizing geometry can look
successful without learning the latent contact assignments. For example, both
homolog curves may occupy the correct chromosome territory while their Delta
directions and individual contact states rotate among equivalent or local
basins. Validation must require matched high-depth recovery of exact, parity,
and structure simultaneously, not structure alone.

Concrete problem: sparse trans contacts can deform a useful cis result. This
matters because a few noisy links should move a chromosome territory before
bending every homolog trajectory. The staged prototype freezes or weakly
deforms cis shapes, then fits rigid placement. Experiment 057B shows the current
rigid stage is insufficient, so validation requires improved ORACLE fixed-shape
centroid correlation before any blind launch.

Concrete problem: an uncalibrated chromosome-pair preference can demand an
impossible spin cycle, such as `r12=+1`, `r23=+1`, `r13=-1`. This matters because
forcing a solution hides missing gauge information. Validation must use blind
q-values, cycle/frustration diagnostics, positive per-spin margins, and null
tests that leave unsupported components unresolved.

Concrete problem: pair expansion gives a six-segment molecule 15 pair terms.
This matters because molecule length becomes an unintended force weight. The
unit-mass shared-assignment model fixes that accounting. Validation combines
mass assertions with paired synthetic comparisons and, only when provenance-
controlled molecule IDs exist, a real-data comparison.

## Scientific conclusion

Experiment 059 rejects Claim A. It supports an implementation-level ability to
recover some unordered chromosome shape, but not stable blind diploid contact
phase. Claim B is established strongly only under ORACLE phase, not for blind
P9016. Claim C remains a possible calibrated selective output under favorable
geometry, but current blind trans evidence does not validate it. The strongest
current operational claim is Claim D: the current implementation and evidence
do not demonstrate stable blind per-cell diploid reconstruction. This synthetic
diagnostic is not an algorithm-independent proof that P9016 is insufficient.
