# p4 semantics and force-mass audit

## Definition

For one canonical haploid bin pair (m=(i,j)), `p4` is ordered as

```text
index 0: P(i copy0, j copy0 | current coordinates and configured score)
index 1: P(i copy0, j copy1 | current coordinates and configured score)
index 2: P(i copy1, j copy0 | current coordinates and configured score)
index 3: P(i copy1, j copy1 | current coordinates and configured score)
```

At E-step output it is a normalized four-state posterior under the implemented
score, prior, temperature, and optional rate weights:

\[
p_s=\operatorname{softmax}_s\left[
\frac{-E(d_s/(u a))+\log\pi_s}{T}+\log\lambda_s
\right].
\]

Calling this quantity a posterior does not prove that (E) is a calibrated
contact likelihood. In particular, the legacy FDG-flat score is a layout
potential with close-distance repulsion and a flat interval.

## Semantic trace

| Stage | Interpretation of p4 | Mathematical use | Caveat |
|---|---|---|---|
| E-step | posterior-like probability | four-state softmax | only as valid as score and prior |
| duplicate aggregation | expected count fraction / mixture weight | raw p4 values are endpoint-canonicalized and averaged/accumulated | raw observations sharing a bin pair receive the aggregate geometry |
| uncertainty | confidence distribution | entropy, pmax, top-two margin | `pU` is normalized entropy, not background probability |
| entropy gating | training multiplier | rho may be multiplied by `1-H/log(4)` | changes total training mass |
| M-step edge expansion | force multiplier | `k_s=base_k*rho*p_s` | posterior enters every edge coefficient |
| historical target conversion | target-distance modifier | `a_s=(n_raw*p_s^gamma)^(-1/3)` | duplicates posterior influence |
| corrected target conversion | no target role | `a=(n_raw)^(-1/3)` | posterior enters coefficient once |
| chromosome-pair update | next-iteration prior source | aggregate p4 by chromosome pair | unrestricted four-state aggregate may violate Z2 cycle consistency |
| readgroup update | shared-assignment mixture | combine pair p4 values across a molecule | historical path reuses posterior rather than a direct observation likelihood |
| serialization | probability and confidence score | p00-p11, parity sums, entropy, margin, pmax | raw rows must compute all fields from raw p4 |
| evaluation | calibrated prediction | exact and summed-parity metrics | exact label needs a declared gauge policy |

## One raw pair in the M-step

Let a raw observation have statistical mass (w_m), posterior (p_{ms}), and
optional confidence gate (g_m\in[0,1]). Before graph aggregation, the corrected
posterior-once coefficient is

\[
c_{ms}=w_m g_m p_{ms},
\qquad
\sum_s c_{ms}=w_m g_m.
\]

The state target scale is independent of posterior:

\[
a_m=\max(\epsilon,n_m)^{-1/3}.
\]

The contribution of one observation to the layout objective is therefore

\[
L_m^{\rm layout}(X)
=\sum_s w_m g_m p_{ms}
  V_{\rm FDG}\left(\frac{d_{ms}}{u a_m}\right).
\]

This formula conserves coefficient mass. It does not assert that
(-V_{FDG}\) is a contact log likelihood.

In the historical posterior-count path,

\[
L_{m,legacy}^{\rm layout}(X)
=\sum_s w_m g_m p_{ms}
  V_{\rm FDG}\left(
  \frac{d_{ms}}{u\max(\epsilon,n_m p_{ms}^{\gamma})^{-1/3}}
  \right).
\]

Here p4 controls both the multiplier and the argument of the force law. With
gamma=1, the target scale changes as (p_s^{-1/3}); in non-flat regimes an early
posterior preference is therefore amplified more strongly than its expected
count fraction alone.

## Aggregated binned pairs

If raw observations (r\in R_m) map to one binned pair, endpoint order is first
canonicalized. A coherent expected-count construction is

\[
C_{ms}=\sum_{r\in R_m} w_r p_{rs},
\qquad
\sum_s C_{ms}=\sum_{r\in R_m}w_r.
\]

If one shared binned posterior (p_{ms}) is used instead, then
(C_{ms}=n_m p_{ms}). Either representation conserves raw mass if it is used
once. It ceases to be a pure mass statement after neighbor-based `base_k`,
state-specific target-distance conversion, hard filtering, or resampling.
Those stages must therefore report both coefficient mass and the actual
effective layout coefficient.

## Gating and sharpening

Posterior sharpening must renormalize:

\[
\tilde p_s=\frac{p_s^{1/T_s}}{\sum_t p_t^{1/T_s}}.
\]

It preserves sum_s p_s=1 but changes uncertainty and thus can indirectly change
training mass if entropy gating is active. A hard gate that drops an
observation changes its mass from (w_m) to zero. A soft gate changes it to
(w_m g_m). Both are permitted only when the lost mass is explicit in the
audit; neither can be described as mass-conserving with respect to the original
observation.

## Molecule-level mass

For molecule (q) with (m) genomic segments and shared copy assignments
(z_1,\ldots,z_m\), the normalized prototype uses

\[
\log W_q(z)=\frac{1}{\binom m2}
\sum_{a<b}\log P(y_{ab}\mid z_a,z_b,X),
\]

followed by exact normalization over the (2^m) assignments for small (m).
The entire molecule has mass one. The `1/choose(m,2)` factor prevents a
six-segment molecule from receiving 15 times the force of a two-segment
molecule merely because it has 15 derived pairs.

The historical C pair-expansion baseline assigns mass one to each derived pair:

\[
M_q^{legacy}=\binom m2.
\]

This is intentionally retained only as a comparator. The failing rows in
`results/phase_diagnostics/mass_conservation.tsv` document this behavior rather
than hiding it.

## Confidence serialization correction

When `set->raw_p4` exists, `hk_blind_write_raw_contact_posterior_tsv` now derives
the following from that same raw posterior row:

```text
pU, psame_raw, pcross_raw, entropy, margin, pmax, rho_output
```

Previously only p00-p11 were raw while the confidence fields came from the
aggregate binned pair. This was an implementation bug in readgroup-enabled
output. It did not affect the historical approved pairs-only arm, where
readgroups were disabled.

## Invariants and expected assertions

The automated suite asserts:

- p4 is finite, non-negative, and sums to one for four-state conditional output;
- endpoint reversal swaps only 01 and 10;
- posterior sharpening preserves unit state mass;
- posterior-once expansion has total coefficient `raw_mass * gate`;
- target distance is posterior-independent in posterior-once mode;
- a molecule has total mass one regardless of segment count;
- raw p4 and serialized confidence fields survive a round trip;
- chromosome flips preserve the appropriately relabeled objective.

The audit table reports the intentionally non-conserving legacy molecule rows as
`pass=0`. A table that contained only passing corrected rows would conceal the
failure being diagnosed.
