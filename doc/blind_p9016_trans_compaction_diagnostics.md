# P9016 Blind Diploid Cis/Trans Diagnostics

This note describes diagnostics and small ablations for the P9016 1Mb blind /
unphased Hi-C driven diploid structural decomposition pipeline with a latent
copy-pair posterior under gauge symmetry.  The latent copy labels `copy0` and
`copy1` have gauge symmetry and must not be interpreted as fixed maternal or
paternal labels, and this pipeline should not be described as genetic haplotype
phasing.  CHARM `mat/pat` labels are used only in post-hoc evaluation scripts,
never in training, initialization, prior construction, seed selection, stopping,
or grid execution.

## Interpretation guardrail: trans ranks are relative diagnostics

The concrete problem is that trans bead-pair Spearman and centroid-distance
Spearman ranks can look better than another run while still being far from a
successful trans reconstruction.  These metrics should be read as relative,
least-bad comparisons inside a diagnostic grid, not as evidence that the blind
unphased Hi-C driven decomposition has recovered the true inter-chromosome
geometry.

This matters because the current 1Mb P9016 setting uses only the pairs file and
has no phase labels or external diploid structural target during training.
Sparse, ambiguous trans contacts can still leave many global layouts nearly
equivalent under the objective, so a higher rank correlation can simply mean
one failed layout is less bad than another failed layout.

A small conceptual example is three chromosome territories A, B, and C.  A
model may correctly learn that A is usually closer to B than to C, giving a
positive centroid-rank signal, while placing all three centroids at the wrong
absolute scale and orientation and assigning no stable biological meaning to
`copy0` versus `copy1`.  That is useful diagnostic signal, but it is not
successful trans reconstruction.

Validation should report trans/centroid Spearman as relative grid-ranking
metrics alongside absolute residuals, null baselines, and compaction measures.
A run should not be called successful on trans structure unless it beats simple
nulls and also satisfies predeclared absolute acceptability checks for
distance scale, residual distribution, and cis preservation.

## Problem: relative vs absolute acceptability

The concrete problem is that rank-based improvements can hide unacceptable
absolute geometry.  A run can improve trans Spearman or centroid Spearman
relative to a previous configuration while still producing chromosome
territories that are too compact, too far apart, or inconsistent with observed
trans contact intensities.

This matters because model selection based only on relative rank can reward a
configuration that is merely the least bad in a small grid.  In the current
pipeline, that would risk over-interpreting blind, unphased structural
decomposition as if it had achieved reliable diploid inter-chromosome placement.

A small conceptual example is two candidate reconstructions.  Run A has
trans Spearman `0.02`; run B has trans Spearman `0.08`.  Run B is better as a
relative diagnostic, but if both have median trans distance residuals several
times larger than within-chromosome distances and both fail contact-decay
residual checks, neither should be considered acceptable.

Validation should separate two decision layers.  First, use relative metrics to
rank ablations within a controlled grid.  Second, apply absolute acceptability
criteria such as calibrated distance-slope ranges, bounded residual quantiles,
reasonable radius-of-gyration ratios, stable homolog separation, and no major
loss of all-chromosome cis metrics.  The document and manifests should preserve
that distinction explicitly.

## Problem: baseline and null diagnostics

The concrete problem is that raw trans/centroid correlations have no
interpretation without baselines.  A positive Spearman value may be meaningful,
trivial, or noise depending on how it compares with chromosome-wise random
layouts, shuffled contacts, centroid-only heuristics, or existing cis-dominated
reconstructions.

This matters because blind evaluation against CHARM after training can easily
look more persuasive than it is.  If a diagnostic does not beat simple nulls,
then it should not be used to claim trans reconstruction, even if it is larger
than zero or larger than one previous failed setting.

A small conceptual example is a shuffled-centroid null that keeps each
chromosome's internal structure but randomly permutes chromosome centroids
across runs.  If the candidate model's centroid Spearman is inside that null
distribution, the model has not shown evidence of specific trans placement
beyond generic chromosome-territory structure.

Validation should include truth-free and eval-only baselines.  Truth-free
baselines include observed trans contact count versus model inverse-distance
correlation, posterior entropy summaries, and contact residuals under shuffled
or permuted trans edges.  Eval-only baselines can include CHARM centroid
permutation, random rigid placement of per-chromosome reconstructions,
cis-preserving trans shuffle, and naive whole-genome scaling controls.  Report
candidate-minus-null effect sizes and null percentiles rather than only raw
correlations.

## Problem: trans residual diagnostics

The concrete problem is that trans rank metrics do not reveal where the model
fails.  A single Spearman number can hide chromosome-pair-specific outliers,
distance-scale bias, copy-pair posterior ambiguity, or a small set of high-count
trans contacts with very poor fitted distances.

This matters because the current optimization balances cis contacts, sparse
trans contacts, backbone, homolog separation, and repulsion.  Without residual
diagnostics, a change that improves centroid rank could still be worsening the
actual trans contact fit or solving one chromosome pair by damaging many others.

A small conceptual example is a model that improves the rank ordering of
chromosome-pair centroids but places chr1-chr5 much too close and chr2-chr8 much
too far away.  The global centroid Spearman may rise, while per-pair residuals
show that the trans structure is not broadly reliable.

Validation should summarize trans residuals by chromosome pair, genomic
distance-free contact-count bins, posterior entropy bins, and copy-pair
posterior mode.  Useful outputs include observed-count versus inverse-distance
calibration curves, signed and absolute residual quantiles, high-count trans
edge residual lists, per-chromosome-pair residual heatmaps, and comparisons to
cis-preserving trans-shuffle nulls.  Improvements should be accepted only when
residual diagnostics improve without creating self-compaction or damaging cis
structure.

## Problem: chr1 cis correlation can hide global failure

The concrete problem is that a high chr1 within-copy distance-map correlation
can be achieved even when other chromosomes or inter-chromosome layout are poor.
This matters because the current model can learn chromosome-internal structure
while leaving relative chromosome placement underconstrained.

A small conceptual example is two chromosomes whose internal bead order and
distances are exactly right, but whose centroids are translated to the wrong
locations.  Per-chromosome cis correlations remain high, yet trans distances and
whole-genome Procrustes RMSD are poor.

Validation uses eval-only CHARM comparison after reconstruction finishes:
per-chromosome cis metrics, centroid distance matrix correlation, sampled trans
bead-pair distance correlation, whole-genome distance summaries, and global vs
per-chromosome Procrustes RMSD.  Trans and centroid correlations are relative
diagnostics only; they are not by themselves evidence of successful trans
reconstruction.  These metrics are written by `evaluate_blind_p9016_grid.py`
and must not feed back into training.

## Problem: per-chromosome self-compaction

The concrete problem is that chromosomes may shrink into compact territories
while preserving distance rank order.  This matters because Pearson/Spearman
correlation can stay high when all distances are scaled down.

A small conceptual example is a model chromosome equal to `0.5 * reference`
around its centroid.  Its distance-map correlation is nearly one, but radius of
gyration ratio and distance-regression slope are about 0.5.

Validation uses radius of gyration ratios, within-chromosome distance quantile
ratios, and regressions of `d_model = slope * d_ref + intercept`.  A consistent
`median_rg_ratio < 1` or `median_cis_distance_slope < 1` supports shrinkage even
when cis correlation is decent.

## Problem: entropy rho can underweight trans contacts

The concrete problem is that trans copy-state posterior can be high-entropy, so
`rho_train_mode=entropy` may assign very low effective contact stiffness to
trans binned contacts.  This matters because sparse trans contacts then cannot
stabilize global chromosome placement against cis contacts, backbone, homolog
separation, and repulsion.

A small conceptual example is a trans bpair with posterior
`p00=p01=p10=p11=0.25`.  Under entropy weighting its confidence is near zero and
its total effective wedge k is near zero; under constant rho it still contributes
force, while remaining copy-state ambiguous.

Validation uses posterior cis/trans summaries: mean and median `pU`,
`rho_train_bpair`, effective k, cis/trans effective-k ratio, observed contact
count versus model inverse-distance Spearman correlations, and trans residual
summaries.  These are truth-free diagnostics based on model coords and unphased
contacts.  Higher trans contact Spearman should be interpreted as a relative
fit signal unless it also beats null baselines and residual checks.

## Small Ablation: rho_train cis/trans handling

The concrete problem addressed is entropy suppression of uncertain trans
contacts.  The small ablations keep defaults unchanged and expose explicit,
manifest-recorded modes such as `constant`, `entropy`, `entropy_with_floor`,
`entropy_cis_constant_trans`, and `entropy_cis_floor_trans`.

The intuitive example is to treat copy assignment confidence separately from
contact evidence: a trans contact can be ambiguous about which copy pair is
responsible while still providing information that two chromosome regions should
not be arbitrarily far apart.

Validation should compare a small diagnostic grid around existing good main
settings.  A direction is promising only if trans/centroid relative ranks,
baseline-adjusted effects, and contact residuals improve without catastrophic
loss of all-chromosome cis metrics or worsening compaction.

## Small Ablation: trans contact multiplier

The concrete problem addressed is weak trans contact force relative to cis,
backbone, and repulsion.  The ablation adds explicit
`contact_k_multiplier_cis` and `contact_k_multiplier_trans` fields; defaults are
both `1.0`.

The intuitive example is a graph with one cis edge and one trans edge.  Setting
the trans multiplier to `2.0` doubles only the trans contact stiffness while
leaving cis, backbone, repulsion, and homolog separation unchanged.

Validation uses force-class diagnostics, trans contact residuals, centroid
metrics, null-adjusted relative ranks, Rg ratios, and cis correlations.  A
useful multiplier should improve trans fit relative to baselines without simply
collapsing chromosomes or destroying cis structure.

## Small Ablation: repulsion multiplier

The concrete problem addressed is force balance: repulsion may separate compact
chromosome territories while weak trans contacts cannot pull them into the
right global arrangement.  The planned mini-grid varies the existing repulsion
multiplier only and keeps current edge-blocking behavior unchanged.

The intuitive example is lowering repulsion enough that trans contacts can
matter, but not so far that homologs or chromosomes collapse into one mass.

Validation uses repulsion-to-trans force ratio, Rg ratio, centroid/trans
relative metrics, homolog separation, null comparisons, and contact residuals.
No new repulsion blocking threshold is enabled by default.

## Small Ablation: expected-count d_scale

The concrete problem addressed is that `raw_count` d_scale gives low-posterior
states the same short target distance as high-posterior states.  The
`expected_count` ablation uses `n_raw * p_ab` for each state distance scale,
with an explicit epsilon floor.

The intuitive example is a bpair with one state at `p=0.9` and another at
`p=0.01`.  The high-posterior state keeps a shorter target distance, while the
low-posterior state gets a larger, weaker geometric target.

Validation compares contact residuals, Rg, homolog separation, trans/centroid
relative metrics, null-adjusted effects, and cis correlation.  This is an
ablation, not a guaranteed improvement.

## Disabled Design: aggregate trans constraints

The concrete problem is that individual trans contacts are sparse and
copy-state ambiguous.  A larger future model change could add truth-free
aggregate constraints, such as weak chromosome-pair centroid attraction or
10Mb-domain aggregate attraction estimated from unphased observed contacts.

The intuitive example is many weak inter-chromosomal contacts between two
chromosomes: each single bpair is ambiguous, but the aggregate density may
contain stable information about relative chromosome placement.

Validation would require improved centroid distance correlation and observed
trans contact versus inverse-distance correlation relative to null baselines,
with improved residual diagnostics, while Rg and cis metrics stay reasonable.
This design must remain disabled until diagnostics and small ablations show it
is needed.

## Disabled Design: anti-collapse constraints

The concrete problem is chromosome self-collapse driven by dense cis contacts.
Future disabled designs include Rg regularization to an unphased scaffold scale,
polymer-scale regularization from unphased contact decay, or class-specific
contact/repulsion balance.

The intuitive example is a chromosome whose internal ranks are right but whose
long-range distances are compressed.  A truth-free scaffold-scale regularizer
could discourage collapse without using CHARM/reference Rg.

Validation would require Rg ratio and long-range cis slope improvements without
hurting short-range contact residuals, trans residuals, or null-adjusted
trans/centroid metrics.  No CHARM, reference, phase, or `mat/pat` values may
enter these training terms.
