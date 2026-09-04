# SNP-free diploid 3D reconstruction diagnostics, experiments 057-059

> **Corrected final freeze.** Formal aliases 119, 120, and 121 contain the
> post-audit 057, 058, and 059 results at 1 Mb. Historical aliases 116-118 are
> retained only as superseded provenance. The corrected run is pinned to source
> commit `fb1eebf`; logical results and formal packages are committed as
> `9fdd348` and `38c1aa4`.

## 1. Executive conclusion

The current evidence does not support full SNP-free exact `00/01/10/11`
contact phasing or globally labeled diploid reconstruction for P9016. The
strongest defensible operational claim is **Claim D**: the current
implementation and evidence have not demonstrated stable blind per-cell
diploid reconstruction. This is not yet an algorithm-independent proof that
P9016 itself lacks sufficient information.

The dominant exact-state failure is upstream of the fully informed geometry
optimizer. With P9016 ORACLE coordinates, the truth state is the nearest of the
four candidate distances for only 18.87% of observations and 16.59% of trans
observations. A monotone distance E-step chooses that nearest state for 100% of
contacts, so it mostly reproduces its current geometry rather than recovering
measured exact phase. Its exact log loss, 1.880, is worse than the uniform
four-state value 1.386.

The current joint M-step is capable when contact phase is known. Across seeds
17, 29, and 43, native cis+trans Hickit obtains mean per-chromosome unordered
distance Spearman 0.966 +/- 0.016. Cis-only reconstruction obtains 0.802 +/-
0.007, but a rigid trans placement stage applied to those fixed cis shapes has
chromosome-centroid distance Spearman only 0.150 +/- 0.067. Thus phase-known
joint optimization works, while the current staged rigid-body prototype does
not.

Optimization remains a second major problem. In the matched synthetic family,
blind exact and parity accuracy remain near 0.25 and 0.50 even at 12,000
contacts, while unordered distance correlation is 0.764-0.870. This is a
failure in a well-specified high-depth regime. At the same time, the corrected
generator-conditional ORACLE-coordinate exact accuracy is only 0.282-0.417,
and trans parity is 0.489-0.531. The result therefore decomposes into both a
blind inference/optimization gap and substantial per-contact ambiguity in the
tested synthetic generators.

The historical FDG-flat E-step is scientifically misspecified as a contact
likelihood: it penalizes very small distances and is exactly flat over a broad
distance interval. The corrected implementation separates E-step observation
score from M-step layout potential, exposes monotone scores, conserves contact
mass with posterior-once weighting, and disables arbitrary gauge-forcing
regularizers in the conservative baseline. These fixes correct the objective;
they do not create information absent from the contacts.

## 2. Current model and mathematical formulation

For haploid bins `i` and `j`, the model creates coordinates
`X_i^0, X_i^1, X_j^0, X_j^1`. A contact has latent state
`s=(a,b)` in `00,01,10,11` and candidate distance

```text
d_ms = ||X_i^a - X_j^b||
r_ms = d_ms / (unit * d_scale_m).
```

The implemented E-step is

```text
z_ms = [-E_estep(r_ms) + log prior_ms] / temperature
       + log_rate_weight_ms
p_ms = softmax_s(z_ms).
```

This is posterior-like only when `-E_estep` is a justified observation log
likelihood. The corrected monotone-log option is

```text
log L(contact | d) = -alpha * log(d + d0),
```

with the C score represented as `alpha=2*k`, `d0=0.05` in normalized layout
units. The monotone-logistic option is

```text
log L(contact | d) = -k * softplus((d-r0)/scale),
```

with `r0=1`, `scale=0.25`. Both are finite at zero and non-increasing with
distance.

The corrected M-step treats the four states as expected edges:

```text
k_ms = raw_mass_m * gate_m * p_ms * base_k_m
d_scale_m = max(epsilon, n_raw_m)^(-1/3).
```

Posterior probability enters the edge coefficient once and does not alter the
state target distance. The M-step still uses an FDG layout potential plus
polymer terms. It is not claimed to equal the observation likelihood.

The isolated symmetric prototype uses

```text
X_i^0 = M_i + Delta_i
X_i^1 = M_i - Delta_i,
```

with separate midpoint, Delta smoothness, and Delta shrinkage terms. This is a
coordinate parameterization, not a biological copy convention.

## 3. Repository implementation audit

The active real-data entry point is `run_blind_p9016_minimal.c`. `blind.c`
contains contact aggregation, E-step scheduling, p4 updates, four-state graph
expansion, regularizers, FDG relaxation, raw serialization, and temporal gauge
stabilization. The standard post-training evaluator is
`eval/evaluate_p9016_baseline.py`. The complete execution graph, formulas,
truth-reading inventory, and experiment 001-056 map are in
`docs/phase_model_audit.md`.

The P9016 input contains 1,703,888 contacts: 1,135,454 cis and 568,434 trans.
The phase-plus-coordinate evaluation subset contains 350,034 contacts: 206,057
cis and 143,977 trans. Resolution is 1 Mb.

The approved pairs-only experiment 056 is retained only as historical lineage.
An exact current replay is not accepted as new evidence because its original
tree was dirty, its terminal audit preceded a later manifest change, its force
scale was subsequently invalidated, and the required GPU runner is unavailable
in this checkout. `baseline_manifest.tsv` records the exact input hashes,
compiler, parameters, and a reproducible command template.

Nine concrete implementation, diagnostic, or reporting defects were found and
corrected:

1. raw readgroup p4 rows were serialized with confidence columns from their
   aggregate binned pair rather than the same raw posterior;
2. same/cross evaluation used parity of the four-state argmax rather than
   comparing `p00+p11` with `p01+p10`;
3. selective coverage could split an equal-confidence tie and manufacture
   coverage;
4. the initial high-depth synthetic negative sampler could exhaust a finite
   unobserved-pair set rather than sampling exposures with replacement;
5. the comparison-table builder initially allowed missing or mislabeled 057B input
   to be silently omitted or relabeled ORACLE. It now fails closed on truth
   boundary, seed set, replicate count, chromosome count, and cross-table
   mismatches;
6. synthetic regularizers originally retained fixed mean-loss strength as depth
   increased, so the contact-depth axis did not increase observation evidence
   relative to the structural prior;
7. all nonzero anchor fractions originally received the same mean anchor-loss
   weight, used independently generated data, and included anchors in evaluation;
   anchors now replace the ordinary observation term under unit mass, use paired
   nested masks, and are excluded from evaluation;
8. molecule training used a shared assignment but output local pairwise softmax
   posteriors; reported pair p4 now marginalizes the same joint `2^m` assignment
   posterior;
9. provenance and packaging used `.strip()` on `git status --porcelain`, which
   destroyed the first status column and falsely classified a result file as
   dirty source. Both readers now remove only terminal newlines and have
   regression tests.

No `01`/`10` endpoint-canonicalization bug was found. Both endpoint orders and
all four pure states pass parse, aggregation, serialization, reload, and
evaluation tests.

## 4. E-step score behavior

The legacy FDG-flat layout energy is

```text
E(r)/k = (0.5-r)^2                    r < 0.5
         0                            0.5 <= r <= 1.5
         (r-1.5)^2                    1.5 < r <= 2
         1.5*(r-2)+0.125/(r-1.5)     r > 2.
```

Because E-step score is `-E`, a distance approaching zero becomes less likely,
and all distances from 0.5 to 1.5 tie. This is a collision/layout shell, not a
monotonic contact-capture likelihood. The production runner previously reached
this mode by default even though DIST2 and LOGDIST2 surfaces existed.

The runner now exposes `fdg_flat`, `dist2`, `logdist2`,
`monotonic_powerlaw`, `monotone_log`, and `monotone_logistic` through
`--estep-score` and `HK_BLIND_P9016_ESTEP_SCORE_MODE`. The score-curve table
covers distance zero, far distances, and multiple temperatures. Every monotone
mode is finite at zero and passes the non-increasing test; only FDG-flat fails.

On P9016 ORACLE coordinates:

| score | exact top1 | exact log loss | exact ECE | same/cross | trans same/cross |
|---|---:|---:|---:|---:|---:|
| uniform prior | 0.4200 | 1.3863 | 0.1700 | 0.7713 | 0.4798 |
| FDG-flat | 0.1984 | 3.3745 | 0.5338 | 0.8474 | 0.6618 |
| DIST2 | 0.1887 | 12.2339 | 0.7266 | 0.8477 | 0.6622 |
| LOGDIST2 | 0.1887 | 1.9071 | 0.4982 | 0.8474 | 0.6617 |
| monotone power law | 0.1887 | 1.6910 | 0.4515 | 0.8474 | 0.6617 |
| monotone log | 0.1887 | 1.8797 | 0.4932 | 0.8475 | 0.6617 |
| monotone logistic | 0.1887 | 8.4350 | 0.7162 | 0.8478 | 0.6623 |

Uniform top1 is inflated by tie-breaking toward prevalent state 00. Proper
scores show that no tested geometric E-step beats uniform exact log loss.

## 5. p4 and force-mass audit

At E-step output p4 is a normalized posterior-like vector in order
`00,01,10,11`. Downstream the historical pipeline also treats it as expected
count fraction, force multiplier, target-distance modifier, confidence object,
training gate input, and next-iteration prior input.

For one observation of mass `w` and gate `g`, posterior-once weighting obeys

```text
c_ms = w * g * p_ms
sum_s c_ms = w * g.
```

The diagnostic confirms total M-step coefficient one for pmax 0.25, 0.50,
0.75, 0.95, and 0.999. Historical posterior-count weighting also changes
`d_scale` as approximately `p^(-1/3)`, so early posterior differences are
amplified twice.

For an m-segment molecule, normalized pair expansion assigns each of
`choose(m,2)` pairs mass `1/choose(m,2)`. The exact shared-assignment prototype
uses the same unit molecule mass. Historical unnormalized expansion has total
mass 3, 6, and 28 for m=3, 4, and 8. The mass table deliberately marks these
legacy rows as failures; m=2 passes because it creates one pair.

## 6. Symmetry and gauge analysis

Without external anchors, every chromosome has an independent Z2 copy gauge.
Flipping both coordinate copies of chromosome c and relabeling all incident
states leaves the SNP-free objective unchanged. Exact state always relabels;
trans same/cross also changes when exactly one endpoint chromosome flips.

Deterministic tests apply flips to one chromosome, several chromosomes, and all
chromosomes. The E-step, M-step, and total objective remain invariant to
floating-point tolerance after the corresponding state relabeling.

The following mechanisms inject or retain arbitrary copy identity:

- global copy-track alignment across chromosomes;
- a fixed displacement direction when homolog points overlap;
- temporal gauge stabilization, which preserves initialization convention;
- unrestricted chromosome-pair four-state priors that do not factor through
  `r_cd=s_c*s_d`;
- truth-derived flips, which are permitted only in evaluation.

Hard minimum homolog separation is gauge-invariant but can invent unsupported
Delta. Global copy-track and hard separation are disabled in the conservative
baseline.

The Z2 solver maximizes

```text
sum_cd w_cd * q_cd * s_c * s_d
```

over connected components and reports cycle inconsistency, frustrated edges,
spin margins, and unresolved components. In the current blind comparison there
are no calibrated relative-gauge edges. The solver correctly leaves all
chromosomes unresolved and is a no-op.

## 7. Experiment 057 results

With fixed P9016 ORACLE coordinates, monotone posterior top1 equals the nearest
geometric state for every contact, but the truth state is nearest for only
18.87% overall and 16.59% trans. Exact target-precision coverage is zero at
0.55, 0.60, 0.70, 0.80, and 0.90 for every principal scorer.

Same/cross is more informative under ORACLE geometry. Monotone-log trans parity
has ECE 0.0637 and evaluation-only coverage 0.803, 0.473, and 0.281 at precision
0.70, 0.80, and 0.90. Cis prevalence inflates all-contact parity, so these
strata must not be merged into an exact-phasing claim.

Current explicit ORACLE-phase M-step results are:

| arm | unordered distance Spearman | homolog-separation Spearman | centroid-placement Spearman |
|---|---:|---:|---:|
| native cis only | 0.802 +/- 0.007 | 0.169 | -0.067 |
| native cis + trans | 0.966 +/- 0.016 | 0.879 | 0.972 |
| fixed cis, rigid trans placement | 0.802 +/- 0.007 | 0.169 | 0.150 |

The near-truth synthetic loop also separates initialization from information.
At zero coordinate noise it produces evaluation-only geometry-aligned exact
0.319, parity 0.551, and unordered structure Spearman 0.942. At noise 2.0 these
become 0.272, 0.543, and 0.591, with zero converged replicates. Random
initialization hurts, but exact contact state is weak even at truth
initialization.

Prior and negative controls show why raw accuracy is misleading. Uniform exact
top1 is 0.420 because state 00 is prevalent. Stratum-preserving endpoint
shuffling retains exact 0.362 and parity 0.825 while destroying contact-specific
state. Randomized geometry and p4 return approximately 0.25 exact and 0.50
parity.

## 8. Experiment 058 results

The P9016 split audit scans all 1,703,888 rows without reading phase columns.
Raw-observation splitting leaks 69,132 exact bin-pair groups and 66,281
neighboring block groups. Exact-pair blocking removes exact leaks but leaves
57,909 neighboring-block leaks. Chromosome-pair and genomic-block splits at 1,
5, and 10 Mb pass the implemented leakage invariants.

P9016 `readID` is `.` throughout the file, so molecule atomicity cannot be
audited. Repeated P9016 held-out fitting and seed selection were not run.
Experiment 058 is therefore **PARTIAL** for real data.

On the seven-seed model-mismatched synthetic diagnostic, blind pseudo-NLL
Spearman ranges -0.750 to 0.607 and the pmax-surprisal score ranges -0.786 to
0.847 across targets/splits. Mean correlations are 0.064 and -0.163. Either
score selects a truth top-quartile seed in only 28.6% of cases. No tested score
is a reliable blind seed selector; 43 of 49 fits converged and failed fits are
retained.

Across 21 seed pairs, raw same/cross agreement is 0.632; agreement after
chromosome flips is 0.644, with cis 0.701 and trans 0.521. Exact-state
agreement rises from 0.277 to only 0.372 after gauge alignment. Mean
chromosome-structure agreement is 0.787, but only 0.222% of contacts are
consistently callable as the same exact state.

## 9. Experiment 059 results

The full synthetic run contains 98 conditions, two generator families, and
three replicates. It adds focused transition points without a broad Cartesian
sweep. Reported values are replicate means; sample standard deviations are in
the machine-readable `replicate_summary.tsv`.

At 12,000 contacts, blind exact/parity remain near random. In the matched
family, geometry-aligned exact ranges 0.261-0.295 and parity 0.501-0.530 across
homolog separation 0.1-2.0. Under mismatch, exact ranges 0.257-0.270 and parity
0.479-0.506. Unordered structure remains 0.764-0.884. High structure quality
therefore does not imply recovered contact phase.

At the same depth and separation 2.0, the generator-conditional ORACLE reaches
exact/parity 0.412/0.642 for the matched family and 0.417/0.716 under mismatch.
However, trans parity is only 0.531 and 0.503. The corrected analysis therefore
shows both an inference gap and stochastic contact-state overlap; unlike the
superseded analysis, it does not call a mismatched monotone-log row a Bayes
ceiling.

At 10,000 contacts, paired nested 0-5% ORACLE anchors increase matched raw
exact/parity from 0.225/0.512 to 0.379/0.581 and mismatched raw exact/parity from
0.226/0.481 to 0.296/0.475. Anchors are excluded from evaluation. At 5%,
geometry-aligned trans parity is still only 0.508 matched and 0.497 mismatched,
and gauge recovery is nonmonotonic. There is no supported small-anchor phase
transition. At 500 and 2,500 contacts, integer rounding produces actual anchor
counts `0,1,2,5,25` and `0,2,12,25,125`; machine tables report both nominal and
realized fractions.

The joint molecule posterior changes geometry-aligned exact accuracy by +0.0160
matched and +0.0007 mismatched relative to normalized pair expansion. Parity
changes by -0.0046 and -0.0065; unordered structure changes by only +0.0127 and
+0.0089. It costs 5.20-fold and 4.42-fold runtime. The molecule optimizer fails
in 8/12 matched and 7/12 mismatched replicates across the four selected
length-4/6 conditions; the paired normalized pair baseline fails in 7/12
replicates in each family. Readchains do not add stable contact-phase
information after molecule-level normalization in this prototype. In the
mismatched length-6, 80-molecule condition, one replicate has zero trans signal
contacts, so 12 trans summaries explicitly use two rather than three finite
replicates.

## 10. Minimal model correction results

| correction | implementation | diagnostic outcome | status |
|---|---|---|---|
| exposed E-step mode | C CLI and environment | removes production hardcoding | ready |
| monotone log/logistic | C and Python | correct monotonicity; exact phase still weak | ready, scale not calibrated |
| separate E/M functions | C configuration | collision shell no longer default likelihood | ready |
| posterior-once M-step | C raw-count target mode | exact coefficient mass conservation | ready |
| background component | isolated Python solver | useful outlier sink; no robust phase gain | experimental |
| conservative regularizers | switches/prototype | preserves independent chromosome gauges | recommended default |
| midpoint/Delta | isolated Python solver | symmetric and permits Delta shrinkage | experimental C port pending |
| staged cis then trans | isolated Python solver | blind structure competitive; ORACLE rigid placement fails | not production-ready |
| explicit Z2 | Python component solver | no calibrated blind edges, correctly unresolved | diagnostic only |
| molecule likelihood | exact `2^m` small-m enumeration | small structure gain, no parity gain, 4.4-5.2x runtime | experimental |

In the fixed synthetic comparison, staged reconstruction has unordered
structure Spearman 0.862, geometry-aligned exact 0.258, parity 0.514, and blind
held-out pseudo-NLL 6.0193. Conservative monotone plus background has structure
0.863, exact 0.236, parity 0.497, and pseudo-NLL 6.0195. This is not a meaningful
winner. The molecule row reaches structure 0.896 and parity 0.495, but uses a
different composite objective and readchain data scope, so its held-out value is
not comparable to the pairwise rows. These are scoped synthetic comparisons,
not truth-selected production winners.

Selected rows from the 28-row machine-readable comparison follow. Synthetic
contact metrics use per-chromosome geometry-selected swaps for evaluation only;
the raw arbitrary-gauge values are separate columns in
`comparison_matrix.tsv`. Callable definitions and structure-score definitions
are not silently equated across evidence scopes.

| configuration | truth_used_in_inference | cis_exact | cis_same_cross | trans_exact | trans_same_cross | callable_coverage | calibration_error | structure_score | heldout_score | runtime | notes |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| historical 056 approved pairs-only | none | 0.5688 | NC | 0.3227 | NC | 0.3405 | NA | 0.8009 | NA | NA | Historical legacy force; callable and parity definitions not comparable |
| corrected historical 078 K1 | none | 0.5265 | 0.9738 | 0.2799 | 0.5026 | 0.0285 | 0.1939 | 0.8125 | NC | NA | Best correctness-audited P9016 evidence; reused final test |
| FDG-flat synthetic proxy | none | 0.1333 | 0.2889 | 0.2595 | 0.4915 | 0.0000 | 0.5068 | 0.6892 | 6.6005 | 0.643 | Layout shell used as endpoint score |
| monotone posterior-once | none | 0.2444 | 0.5556 | 0.2369 | 0.4988 | 0.0000 | 0.0616 | 0.8463 | 6.1049 | 0.582 | Correct mass, but phase remains near chance |
| monotone + background, conservative | none | 0.2222 | 0.4889 | 0.2369 | 0.4976 | 0.0000 | 0.0645 | 0.8627 | 6.0195 | 0.662 | Copy-track and hard separation disabled |
| staged cis then trans | none | 0.2444 | 0.5556 | 0.2582 | 0.5128 | 0.0000 | 0.0464 | 0.8616 | 6.0193 | 0.726 | No reliable advantage over conservative joint fit |
| staged + Z2 | none | 0.2444 | 0.5556 | 0.2582 | 0.5128 | 0.0000 | 0.0464 | 0.8616 | 6.0193 | 0.726 | Zero calibrated edges; all three chromosomes unresolved; exact no-op |
| molecule-level normalized composite | none | 0.2222 | 0.5556 | 0.2570 | 0.4933 | 0.0000 | 0.0792 | 0.8956 | 2.8680 NC | 5.680 | Different readchain data and composite objective; held-out score not comparable |
| ORACLE coordinates, monotone-log proxy | ORACLE_COORDINATES | 0.2000 | 0.4667 | 0.2527 | 0.5128 | 0.0000 | 0.0324 | NA | NA | 0 | Fixed synthetic comparison condition only |
| ORACLE contact phase, joint Hickit | ORACLE_CONTACT_PHASE | NA | NA | NA | NA | NA | NA | 0.9661 | NA | 48.123 | Structure-only 057B result, three seeds |

`NC` means not comparable under the stated metric contract, not missing zero.

## 11. Blind model-selection reliability

No new P9016 configuration is declared best from truth. The corrected
historical 078 K1 arm remains the strongest correctness-audited blind real-data
evidence because its phase labels were opened only after the blind scorer
closed. It is not a newly selected winner, and its final-test denominator is
reused rather than a model-selection set.

K1 has cis/trans exact accuracy 0.5265/0.2799, cis/trans parity
0.9738/0.5026, mean per-chromosome cis distance Spearman 0.8125, and exact ECE
0.1939 on 31,973 evaluated contacts. Cis parity matches the 0.975 cis-same
prevalence; trans parity is chance. At pmax >=0.9 it calls 2.846% of contacts at
60.66% accuracy, with 1.726% correct-call recall.

The recommended objective is chosen from symmetry, mass, and likelihood
validity, not by selecting the best truth row. A production seed selector still
requires a preregistered blocked score that correlates consistently with exact,
parity, and structure targets on independent cells.

## 12. What is identifiable without SNPs

The likelihood can represent two unordered homolog trajectories per
chromosome. Under ORACLE contact phase, current joint Hickit reconstructs this
target strongly. Synthetic blind runs can also recover high per-chromosome
distance correlation, although they do not recover exact state.

Same/cross can be more informative than exact state when geometry and homolog
separation are favorable. P9016 ORACLE-coordinate trans parity is 0.662 while
exact state is 0.166. This is a fixed-reference-geometry scorer diagnostic, not
a formal P9016 upper bound; current blind P9016 and synthetic trans parity remain
near chance.

Calibrated selective parity is therefore a scientifically plausible future
output. It is not yet validated as a robust blind P9016 product.

## 13. What is not identifiable without external anchors

Absolute copy0/copy1 identity for each chromosome is not identifiable under the
SNP-free likelihood. A globally aligned exact `00/01/10/11` state requires one
gauge choice per connected chromosome component. Those choices cannot be
called biological copy identity without external anchors.

Sparse trans observations do not automatically determine relative chromosome
gauges. The blind comparison supplies no calibrated Z2 edges, and small ORACLE
anchor fractions do not create a reproducible trans phase transition.

An individual stochastic contact also need not originate from the closest of
four reference distances. In the controlled synthetic generators, even the
generator-conditional ORACLE-coordinate exact decision accuracy is at most
0.417 in the tested high-depth slice, and trans parity remains near chance.
Gauge alignment cannot repair this observation ambiguity. The P9016
fixed-coordinate result is diagnostic of score/data disagreement, not an
algorithm-independent information bound.

## 14. Recommended production objective

Use a monotone, saturating contact likelihood in the E-step, independently
configured from the M-step layout potential. Start with

```text
log L(contact | d) = -alpha * log(d + d0)
```

plus an optional constrained background component. Use gauge-symmetric priors
only. In the M-step use raw contact mass times posterior exactly once; keep
target distance independent of posterior. Disable global copy-track alignment
and hard minimum separation unless a preregistered diagnostic demonstrates a
held-out benefit without breaking symmetry or inventing Delta.

Retain midpoint/Delta as an experimental symmetric parameterization. Do not
promote rigid cis-then-trans placement until it succeeds in the ORACLE
fixed-shape control. Run Z2 only when calibrated blind relative edges exist;
otherwise report unresolved components. Use molecule likelihood only when real
molecule IDs and a source-mass contract are available.

Conceptual example: if one contact has posterior `(0.70,0.10,0.10,0.10)`, its
total M-step coefficient remains one observation, not four, and state 00 gets
70% of that mass. A closer candidate is never penalized by collision repulsion
inside the E-step. Validation requires score monotonicity, mass invariance,
ORACLE decomposition, blocked held-out stability, and negative controls before
any real-data promotion.

## 15. Recommended output representation

The primary structure output should be two **unordered homolog trajectories per
chromosome**, with explicit chromosome-component gauge labels and uncertainty.
Do not call either trajectory maternal or paternal.

Contact output should contain the normalized four-state distribution only as a
gauge-dependent latent representation, plus separately reported same/cross
probability, background probability, entropy, margin, and calibrated callable
flag. Every table must state the gauge policy, cis/trans stratum, evaluated N,
coverage, and whether a threshold was calibrated blind or only evaluated with
truth.

Globally aligned exact-state metrics belong in ORACLE/evaluation tables, not in
the biological identity of a blind reconstruction.

## 16. Remaining risks and next experiments

1. **Real blocked validation is missing.** Problem: 058 did not fit repeated
   P9016 seeds under block splits. Why it matters: synthetic score correlation
   may not transfer. Example: repeated 1 Mb contacts can leak across a raw
   split. Validation: preregister one 5 or 10 Mb split, train a fixed small seed
   ensemble, and keep phase sealed until scores close.
2. **The observation likelihood lacks an exposure model.** Problem: contact-only
   pseudo-NLL omits noncontacts. Why it matters: a score can reward collapsing
   all distances. Example: every observed edge becomes likely with no penalty
   for predicting unobserved edges. Validation: add sampled, chromosome- and
   genomic-distance-matched exposures and test calibration on held-out blocks.
3. **Rigid staged placement fails its ORACLE control.** Problem: frozen cis
   shapes plus trans contacts give centroid Spearman 0.150. Why it matters:
   staging cannot protect cis shape if it cannot place territories. Example:
   sparse trans links permit many rotations/translations. Validation: test
   analytic rigid gradients, multiple starts, and only then weak deformation,
   without using truth to select the run.
4. **Blind Z2 evidence is absent.** Problem: no calibrated q-values connect
   chromosomes. Why it matters: forcing spins would create copy identity.
   Example: a disconnected chromosome has exactly two equivalent labels.
   Validation: calibrate pair evidence on synthetic nulls, require positive
   margins and cycle consistency, and leave weak components unresolved.
5. **Real molecule IDs are unavailable.** Problem: P9016 pairs use `readID='.'`.
   Why it matters: the molecule model cannot be evaluated or mass-normalized on
   this input. Validation: recover a provenance-controlled readchain source and
   compare pair versus molecule likelihood at fixed source-molecule mass.

### Reproduction and result freeze

All Python commands run in Conda environment `analysis`. The exact full
diagnostic command was:

```bash
conda activate analysis
MPLCONFIGDIR=/tmp/phase-matplotlib-fb1eebf \
  ./scripts/run_phase_diagnostics.sh --full --seed 17 --threads 4 \
  --output-dir results/phase_diagnostics
```

The non-overwriting formal packaging command was:

```bash
conda activate analysis
python scripts/package_phase_formal_results.py --execute \
  --sequence-start 119 --timestamp 20260904_140536 --seed 17
```

The resulting folders are
`test_res/119-20260904_140536-phase_exp057_oracle_decomposition_1m`,
`test_res/120-20260904_140536-phase_exp058_heldout_invariants_partial_1m`,
and
`test_res/121-20260904_140536-phase_exp059_identifiability_synthetic_1m`.
Each contains a README, plots, machine-readable tables, a package receipt, and
an artifact hash manifest. Repeating formal packaging requires a new sequence
and timestamp because overwrite is intentionally rejected.

## Explicit answers

- **Does the current posterior mainly reproduce the nearest current geometry?**
  Yes. Monotone top1 equals the nearest state for 100% of P9016 ORACLE-coordinate
  contacts, while truth is nearest for only 18.87%.
- **Is the current E-step score a valid monotonic contact likelihood?** No for
  historical FDG-flat. The new monotone modes satisfy the shape constraint but
  do not solve exact-state identifiability.
- **Can the current M-step reconstruct structure when contact phase is known?**
  Yes for joint cis+trans Hickit: unordered distance Spearman 0.966 +/- 0.016.
  Cis-only is weaker, and fixed-shape rigid trans placement fails.
- **Is poor performance mainly local minima or lack of information?** Both.
  Matched high-depth failure demonstrates optimization/implementation limits;
  generator-conditional truth-geometry performance demonstrates stochastic
  per-contact ambiguity independently of score misspecification.
- **Is same/cross substantially more identifiable than exact state?** Under
  favorable or ORACLE geometry, yes; in current blind trans inference, no.
- **Does chromosome-level Z2 alignment explain trans failure?** No. It is a
  necessary ambiguity, but no calibrated blind edges exist and anchors do not
  rescue trans parity.
- **Does a small number of phased anchors cause a phase transition?** No
  reproducible transition was observed through 5% ORACLE anchors.
- **Can a blind-visible held-out metric reliably choose a better seed?** Not in
  the synthetic audit. Real P9016 repeated held-out fitting remains `NOT_RUN`.
- **Do readchains add information after molecule normalization?** They improve
  prototype structure by only 0.009-0.013 Spearman and do not improve parity,
  at 4.4-5.2x runtime with 15/24 molecule-fit replicates failing across the two
  families. P9016 cannot be tested from the current pairs file.
- **Which current regularizers inject arbitrary copy identity?** Global
  copy-track alignment, fixed displacement direction, temporal gauge
  convention, and unrestricted pair-state priors. Hard separation does not
  choose a label but can invent homolog separation.
- **What fraction can be called at useful precision?** P9016 K1 pmax >=0.9
  covers 2.846% at only 60.66% precision. Even ORACLE-coordinate exact target
  coverage is zero at precision 0.55-0.90. ORACLE-coordinate trans parity has
  descriptive coverage 0.803/0.473/0.281 at precision 0.70/0.80/0.90.
- **What should the project claim?** Operational Claim D: the current method and
  evidence do not demonstrate stable blind diploid reconstruction. Claim B is
  demonstrated as an ORACLE-phase capability, Claim C remains a future
  selective target, and Claim A is unsupported. This does not prove an
  algorithm-independent impossibility for P9016.
