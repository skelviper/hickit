# SNP-free diploid phase model implementation audit

> Formal aliases 119-121 are the corrected result freezes. Historical aliases
> 116-118 remain only as reviewed provenance; independent review found blocking
> synthetic mass, anchor, gauge, and molecule-posterior issues in alias 118.

## Scope and evidence boundary

This audit describes branch `phase` at commit
`cfa590306e953cdc3a3ea35860b803a5669b895c`, plus the explicitly listed
diagnostic corrections through source commit
`fb1eebf522578ed4c31ba49d4d7ffbc8e3bba33e` on branch
`phase-diagnostics-057-059`. The training surface is the raw contact portion of
the P9016 `.pairs.gz` file at 1 Mb.
`phase0`, `phase1`, SNP labels, and CHARM/3DG coordinates are excluded from
blind training. They are read only by evaluation code and by experiments whose
configuration is marked `ORACLE`.

The historical experiment 056 `approved_baseline_pairs_only` is retained as a
lineage comparator, not as evidence that the intended objective is correct. It
was produced from a dirty tree and its terminal audit preceded a later manifest
mutation. Later force-correctness work also invalidated its old mechanical force
scale. The exact baseline contract is in
`results/phase_diagnostics/baseline_manifest.tsv`.

The three inference targets are reported separately throughout:

1. two unordered homolog trajectories per chromosome;
2. same-copy versus cross-copy contact parity;
3. exact `00/01/10/11` state after a declared chromosome-gauge alignment.

`copy0` and `copy1` are gauge labels. Neither is called maternal or paternal in
blind output.

## Execution graph

The active path starts in `run_blind_p9016_minimal.c` and enters the scheduled
loop in `blind.c`:

```text
P9016 raw contact fields
  -> parse chr/position only
  -> bin at configured resolution
  -> canonicalize haploid endpoint order
  -> aggregate duplicate haploid bin pairs and retain raw multiplicity n_raw
  -> initialize two coordinate copies per bin
  -> repeat EM-like schedule
       -> optional blind chromosome-pair prior update
       -> E-step: four distances -> four scores -> normalized p4
       -> optional readgroup joint posterior update
       -> uncertainty, entropy, rho, and optional training gate
       -> M-step: expand each binned pair into four weighted diploid edges
       -> add polymer and independently switchable regularizers
       -> FDG relaxation with endpoint backtracking
       -> temporal per-chromosome gauge stabilization and p4 relabeling
  -> write coordinates, binned p4, optional raw p4, diagnostics, manifest
  -> post-training evaluator may read SNP phase and CHARM/3DG
```

The training structure `hk_blind_pair` contains chromosome indices and
positions, not phase. The standard evaluator reads phase-bearing pairs only
after the training artifacts exist. The prototype diagnostic framework keeps
the same boundary in `phase_diagnostics/p9016.py`: the module and every caller
are explicitly named ORACLE.

## Implemented E-step

For binned observation (m=(i,j)) and state
(s=(a,b)\in\{00,01,10,11\}), the code computes

\[
d_{ms}=\lVert X_i^a-X_j^b\rVert,
\qquad
r_{ms}=\frac{d_{ms}}{u\,a_m},
\]

where (u) is the layout unit and (a_m) is the contact-specific distance
scale. In `hk_blind_posterior_from_energy_with_log_norm`, the unnormalized
state score is

\[
z_{ms}=\frac{-E(r_{ms})+\log \pi_{ms}}{T}+\log \lambda_{ms},
\qquad
p_{ms}=\frac{\exp z_{ms}}{\sum_t\exp z_{mt}}.
\]

The subtraction of the maximum score makes the softmax numerically stable.
`log_rate_weight` is intentionally outside the temperature term. The state
order in `hk_blind_bpair_posterior_from_coords_score_mode_with_log_norm` is
exactly `00=i0-j0`, `01=i0-j1`, `10=i1-j0`, `11=i1-j1`.

The implemented energy functions are:

\[
E_{\rm dist2}(r)=k r^2,
\]

\[
E_{\rm logdist2}(r)=k\log(1+r^2/10^{-6}),
\]

\[
E_{\rm monotonic\_powerlaw}(r)=k\log(1+r^2),
\]

\[
E_{\rm monotone\_log}(r)=2k\log(r+0.05),
\]

and

\[
E_{\rm monotone\_logistic}(r)
=k\operatorname{softplus}\left(\frac{r-1}{0.25}\right).
\]

Thus the corresponding monotone-log observation score is
(-2k\log(r+0.05)), and the logistic observation score is
\(\log\sigma((1-r)/0.25)\) up to multiplication by (k).

The legacy `fdg_flat` E-step instead reuses the layout shell potential:

\[
E_{\rm FDG}(r)/k=
\begin{cases}
(0.5-r)^2, & r<0.5,\\
0, & 0.5\le r\le1.5,\\
(r-1.5)^2, & 1.5<r\le2,\\
1.5(r-2)+0.125/(r-1.5), & r>2.
\end{cases}
\]

It is not a monotone contact likelihood: it penalizes very close points and is
exactly flat over (0.5\le r\le1.5). The historical production default used
`fdg_flat`; the corrected runner defaults to `monotone_log`, exposes all modes
through `--estep-score` and `HK_BLIND_P9016_ESTEP_SCORE_MODE`, and retains
`legacy_fdg_flat` solely to reproduce omitted-environment historical semantics. Numerical curves and
monotonicity flags are in `results/phase_diagnostics/estep_score_curves.tsv` and
the companion PDF.

## Implemented M-step and optimizer

For aggregate binned pair (m), `hk_blind_bpair_expand_weighted_edges_mode_gamma_ex`
creates four layout edges. Ignoring optional gates for one line, state (s)
gets

\[
k_{ms}=k_m\,\rho_m\,p_{ms}.
\]

In the historical `posterior_count` mode,

\[
a_{ms}=\left(\max(\epsilon,n_m p_{ms}^{\gamma})\right)^{-1/3}.
\]

With the historical gamma=1, posterior therefore changes both the edge
coefficient and its target-distance scale. In the corrected `raw_count` /
`posterior_once` mode,

\[
a_m=\left(\max(\epsilon,n_m)\right)^{-1/3},
\]

which is independent of state posterior. The posterior then enters exactly
once through (k_{ms}). Because the FDG derivative also depends on
(d/(u a)), the legacy effective force is not merely linear in (p); away
from flat regions its scale contains approximately (p^{4/3}) when gamma=1.

The M-step remains a force-directed shell potential plus polymer and optional
regularizers. It is an optimizer/layout objective, not an asserted
`log P(contact | distance)`. The corrected E-step likelihood configuration is
therefore independent from the M-step potential configuration.

The CPU physical-force path applies the chain-rule factor
`1/(unit*d_scale)`. Endpoint backtracking rejects a relaxation proposal when
the monitored objective increases. This protects against a single gross
overshoot, but not against a wrong objective, self-reinforcing EM updates, or a
poor local basin.

## Where p4 enters

`p4` is created as a normalized E-step posterior. It is subsequently used as:

- an expected state fraction in duplicate-bin-pair aggregation;
- a force multiplier in every four-state M-step edge;
- historically, a modifier of count-derived target distance;
- an input to entropy, `pmax`, margin, and `rho_output`;
- an optional training inclusion/downweighting signal;
- an input to chromosome-pair priors and later E-steps in configurations that
  enable those features;
- a confidence/calibration object in serialization and evaluation.

`pU=H(p4)/log(4)` is normalized entropy. It is not an outlier probability and
must not be interpreted as a background mixture state. The full semantic trace
and one-observation formulas are in `docs/p4_semantics.md`.

## Raw count and force mass

Raw multiplicity `n_raw` enters at least three distinct computations:

1. the base distance scale, conventionally (n_{raw}^{-1/3});
2. neighbor-density heuristics used to construct and clip `base_k`;
3. `posterior_count`, where (n_{raw}p_s^\gamma) is treated as a state-specific
   effective count.

Expanding one binned pair into four states does not by itself multiply mass by
four because sum_s p_s=1. Entropy gating changes total coefficient mass to rho,
which is intentional but must be reported. Neighbor-dependent `base_k`,
state-specific target distances, aggregation, and resolution changes mean that
coefficient conservation is not equivalent to force or energy conservation.
The new audit therefore reports all three separately in
`results/phase_diagnostics/mass_conservation.tsv`.

The historical C readgroup path enumerates shared segment-copy assignments but
then emits every derived pair at unit mass. A molecule with (m) segments can
therefore contribute O(m^2) pair mass. The isolated molecule likelihood now gives
one unit total mass to the molecule and divides pairwise terms by
`choose(m,2)` during exact `2^m` enumeration. P9016 `.pairs.gz` has `readID='.'`,
so that correction cannot be evaluated on P9016 pairs without a separate,
provenance-controlled molecule source.

## Endpoint order and canonicalization

Canonicalization stores whether the input endpoints were swapped. Conversion
between raw and canonical p4 is the involution

```text
00 -> 00
01 -> 10
10 -> 01
11 -> 11
```

No 01/10 error was found in the inspected production path. The new tests cover
both input orders and the full
`parse -> canonicalize -> aggregate -> p4 -> serialize -> reload -> evaluate`
round trip. They fail on any state-order permutation.

One concrete bug was found and corrected: when readgroup joint inference stored
per-raw `raw_p4`, the raw posterior columns were serialized alongside aggregate
binned-pair entropy, margin, `pmax`, and `rho_output`. The writer now recomputes
all confidence columns from the same per-raw posterior written on that line.

## Truth-reading inventory

Blind training reads only chromosome and position fields. Truth enters at these
declared boundaries:

- `eval/evaluate_p9016_baseline.py` reads `phase0/phase1` and CHARM/3DG after
  training for contact and structure evaluation;
- `phase_diagnostics/p9016.py` reads them only for logical experiment 057A,
  which is marked `ORACLE_COORDINATES`;
- historical experiment 083 supplies the `ORACLE_CONTACT_PHASE` M-step control;
- synthetic truth is generated inside experiment 059 and is used for evaluation
  or explicitly configured anchor contacts only;
- near-truth coordinates in 057C are explicitly `ORACLE_NEAR_TRUTH_INITIALIZATION`.

No truth is used to select a blind seed, construct a blind chromosome-pair
prior, choose chromosome flips in a blind result, or filter ordinary blind
contacts. Historical post hoc SNP-based chromosome swaps are labeled eval-only,
not blind gauge recovery.

## Copy-swap symmetry audit

With no external anchor, flipping both coordinate copies for chromosome (c)
and relabeling every incident state leaves a valid SNP-free objective unchanged.
For an edge between chromosomes (c,d), the transformed state is
((a\oplus f_c,b\oplus f_d)). Same/cross parity changes only when exactly one
endpoint chromosome flips; exact labels always relabel. The objective must be
compared after that corresponding relabeling.

The following components preserve whole-chromosome Z2 symmetry:

- contact likelihood and four-state expected edge expansion;
- midpoint backbone and bending terms;
- local copy-track continuity under a whole-chromosome flip;
- centroid or pointwise homolog separation magnitudes;
- symmetric Delta smoothness and shrinkage.

The following components can inject or retain arbitrary copy identity:

- global copy-track alignment, which aligns chromosome Delta directions and
  reduces independent chromosome gauges to one genome-wide flip;
- a fixed displacement direction used when homolog coordinates overlap;
- temporal gauge stabilization, which is truth-free but makes the reported
  label depend on initialization and optimization history;
- unrestricted chromosome-pair four-state priors, which need not factor as
  `r_cd=s_c*s_d` and can be cycle-inconsistent;
- any truth-derived post hoc chromosome swap, which is evaluation only.

Hard minimum separation is Z2-invariant, but it can create unsupported Delta
when data prefer homolog collapse. For that reason it is disabled in the
conservative comparison even though it does not itself select copy0.

The deterministic tests flip one, several, and all chromosomes and compare the
E-step, M-step, total objective, parity after the mathematically required
relabeling, and exact labels. Explicit Z2 synchronization is implemented with
connected components, margins, cycle inconsistency, and unresolved components.

## Held-out and seed-selection audit

The historical external split groups exact binned pairs and available molecule
keys, but does not block nearby genomic loci. The legacy C held-out score also
uses an FDG-flat wrapper rather than the configured E-step mode, so it is not a
likelihood for the active observation model.

The new split module supports:

- `raw_observation`;
- `molecule_or_readchain`;
- `exact_haploid_bin_pair`;
- `genomic_block` at configurable 1, 5, or 10 Mb;
- `chromosome_pair`.

Union-find grouping guarantees that a source molecule cannot cross the split.
Leakage is audited separately for duplicate rows, exact bin pairs, molecules,
and neighboring genomic blocks. Held-out contact scores are called
`pseudo-NLL` when the non-contact exposure denominator is absent.

## Metric definitions

The centralized definitions are in `phase_diagnostics/metrics.py`:

- exact top1: `argmax_s p_s == truth_s` under the stated gauge;
- exact log loss: `-mean(log p_truth)`;
- exact Brier: mean four-class squared error;
- same/cross posterior: `(p00+p11, p01+p10)` before taking argmax;
- same/cross log loss and Brier: computed from that two-class posterior;
- entropy: `-sum_s p_s log p_s`;
- ECE: weighted absolute confidence-accuracy gap over fixed bins;
- selective coverage at precision (q): longest confidence-ranked prefix whose
  cumulative observed precision is at least (q);
- unordered chromosome RMSD: minimum over the two whole-chromosome copy swaps
  after Procrustes alignment;
- distance correlation: Pearson/Spearman over within-structure pair distances;
- gauge recovery: reported only on synthetic or eval-only truth;
- seed agreement: exact labels after optimal chromosome flips, parity agreement,
  distance-matrix correlation, and posterior Jensen-Shannon divergence.

The old standard evaluator had computed same/cross by taking four-state argmax
first and then its parity. That is not equivalent to comparing summed parity
posteriors. It is corrected here. Historical same/cross values produced by the
old evaluator are marked legacy and are not silently mixed with centralized
metrics.

## Historical experiments 001-056

The reproducible mapping is
`results/phase_diagnostics/experiment_lineage_001_056.tsv`. It records every
physical run, the selected canonical run where IDs were reused, training source,
truth boundary, entrypoint, and provenance grade. Important exceptions are:

- 013 is a P1006 control, not P9016 evidence;
- 015, 016, and 018 are ORACLE/reference-derived training controls;
- 020 lacks the required README and manifest and is marked
  `MISSING_REQUIRED_PROVENANCE`;
- 018, 021, 024, 027-029, 031, 041, and 054 have multiple physical runs;
- 053-056 include contacts.seg/readgroup/readchain source experiments;
- only the explicitly named arm inside 056 is the historical approved
  pairs-only comparator.

Existing logical experiment numbers 057-059 were already occupied in the
historical repository. The new requested experiments therefore write their
scientific outputs under `results/phase_diagnostics/057`, `058`, and `059`, but
their corrected formal `test_res` aliases are 119, 120, and 121. Historical
aliases 116-118 were not overwritten and are superseded. The corrected aliases
were assigned only after fail-closed packaging validation.

## Failure-mode decomposition before new experiments

### E-step

Confirmed risks are the non-monotone FDG-flat score, geometry-nearest-state
copying, uncalibrated distance scale, and priors/count-density terms mixed into
what is treated as a posterior. Experiment 057A asks whether state inference is
informative even with oracle coordinates.

### M-step

Confirmed risks are duplicated posterior amplification, layout-shell semantics
standing in for a likelihood, and trans forces deforming sparse cis shapes.
Experiment 057B asks whether the native optimizer reconstructs structure when
phase is fixed to truth.

### Optimization

Backtracking prevents some numerical explosions but does not identify the
correct basin. Random initialization, EM self-reinforcement, and seed-to-seed
variation remain plausible. Experiments 057C and 058D test those mechanisms.

### Gauge

Independent chromosome Z2 labels are unidentifiable without anchors. Unfactored
pair priors and global copy-track can impose arbitrary labels, while weak trans
evidence can leave synchronization unresolved. Experiments 058D and 059's
anchor phase diagram separate gauge-equivalent variation from different
structures.

### Information limit

P9016 contacts contain no recoverable molecule IDs, exact state has four labels
while parity has two, and sparse trans contacts must estimate relative
chromosome gauge and placement simultaneously. Experiment 059 tests depth,
homolog separation, trans fraction, background, molecule length, and external
anchor fraction without expanding into an indiscriminate hyperparameter sweep.

## Minimal correction hypotheses and validation

Each correction addresses a concrete failure rather than being a generic model
extension:

- **Monotone E-step.** Problem: FDG collision repulsion makes an extremely close
  contact less likely. This matters because E-step confidence then reflects a
  mechanical shell, not capture probability. Example: distances 0.1 and 0.8
  should not favor 0.8 solely because 0.1 lies inside a repulsive shell.
  Validation: curve monotonicity tests and 057A oracle-coordinate calibration.
- **Posterior-once weighting.** Problem: a posterior 0.95 both strengthens an
  edge and shortens its target, creating superlinear positive feedback. This
  matters in EM because an early random preference can become geometry.
  Example: changing p from 0.5 to 0.95 should change expected mass by 1.9-fold,
  not also alter the observation's physical count. Validation: mass table,
  symmetry tests, held-out pseudo-NLL, and synthetic seed stability.
- **Background component.** Problem: every observation is forced into one of
  four homolog states even when it is random ligation. This matters most for
  sparse trans data. Example: one distant outlier should increase background
  probability instead of pulling a chromosome pair together. Validation:
  mismatched synthetic background phase diagram and calibration.
- **Midpoint/Delta parameterization.** Problem: arbitrary copy tracks and hard
  separation can manufacture identity. This matters because weak data should
  permit homolog differences to shrink. Example: `X0=M-Delta`, `X1=M+Delta`
  allows `Delta=0` without choosing a privileged copy. Validation: objective
  symmetry, high-depth synthetic recovery, and homolog-separation correlation.
- **Staged cis then trans.** Problem: sparse trans forces can destroy an already
  useful chromosome shape. Example: optimize a chromosome's internal distance
  matrix first, then allow trans contacts mainly to move the chromosome as a
  body. Validation: compare held-out score and within-chromosome distance
  preservation against joint optimization.
- **Z2 synchronization.** Problem: pairwise relative-copy suggestions can be
  mutually inconsistent. Example: `r12=+1`, `r23=+1`, `r13=-1` is a frustrated
  cycle and cannot define chromosome spins. Validation: cycle/frustration,
  component confidence, and anchor-fraction phase diagram without truth-based
  spin selection.
- **Molecule likelihood.** Problem: all-pairs expansion gives a long molecule
  quadratic force mass. Example: a six-segment molecule must not count as 15
  independent molecules. Validation: exact `2^m` likelihood on synthetic
  molecules, mass assertions, runtime, and comparison to normalized pairwise
  expansion.
