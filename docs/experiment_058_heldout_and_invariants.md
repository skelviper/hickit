# Experiment 058: held-out, invariant, and model-selection audit

> Formal alias 120 is the corrected result freeze. Historical alias 117 is
> retained only as superseded provenance; the corrected synthetic metrics use
> explicitly labeled evaluation-only geometry gauge alignment while raw
> arbitrary-gauge values remain separate.

## Status and scope

Logical experiment 058 is **PARTIAL**. The corrected formal result-folder alias is 120 and
the resolution is 1 Mb. The P9016 read-only split/leakage audit is complete. A
fixed model-mismatched synthetic data set is used for repeated held-out scoring
and truth evaluation because repeated P9016 baseline fits, real-data seed
stability, and blind seed selection were not run. No new P9016 blind model was
selected or trained from evaluation results.

Outputs are under `results/phase_diagnostics/058/`:

- `split_manifest.tsv`: grouping and leakage counts;
- `heldout_metrics.tsv`: blind scores and evaluation-only outcomes;
- `model_selection.tsv`: rank correlations and seed-selection regret;
- `seed_stability.tsv`: pairwise seed agreement;
- `figures/`: split and stability summaries;
- `provenance.json`: command, seed, versions, environment whitelist, and git
  state.

## Truth and gauge boundary

P9016 split keys use chromosome, position, optional read/molecule ID, and the
configured block size. They do not read phase columns. Synthetic truth metrics
are computed only after a candidate's held-out pseudo-NLL and contact energy
exist. The selection rule never reads truth.

Exact-state seed agreement is shown both in raw gauge and after an optimal
chromosome flip for agreement analysis. The latter does not label either copy
as maternal or paternal and is not used to select a model. Same/cross is
reported separately because a one-chromosome flip changes trans parity even
though it leaves the SNP-free objective equivalent after state relabeling.

## 058A: current split leakage

The P9016 scan covers all 1,703,888 contacts. Its `readID` field is `.`; there
are zero nonempty molecule identifiers. A molecule-aware split therefore
cannot recover molecule grouping from this file.

| split | train | held-out | duplicate leak | exact pair leak | molecule leak | neighbor-block leak | pass |
|---|---:|---:|---:|---:|---:|---:|---|
| raw observation | 1,363,110 | 340,778 | 0 | 69,132 | 0 | 66,281 | no |
| molecule/readchain | 1,363,110 | 340,778 | 0 | 68,932 | 0 | 65,835 | no |
| exact haploid bin pair | 1,363,110 | 340,778 | 0 | 0 | 0 | 57,909 | no |
| chromosome pair | 1,362,548 | 341,340 | 0 | 0 | 0 | 0 | yes |
| genomic block 1 Mb | 1,363,106 | 340,782 | 0 | 0 | 0 | 0 | yes |
| genomic block 5 Mb | 1,363,108 | 340,780 | 0 | 0 | 0 | 0 | yes |
| genomic block 10 Mb | 1,363,108 | 340,780 | 0 | 0 | 0 | 0 | yes |

The zero duplicate-leak count does not rescue raw splitting: many distinct raw
rows map to the same 1 Mb pair. Exact-pair blocking removes this leak but still
places neighboring blocks on both sides. The union-find blocked splitter keeps
every available source molecule atomic and makes all requested constraints
transitive.

## 058B: blocked split definitions

`raw_observation` hashes each row independently, except that source molecule
atomicity is always imposed when an ID exists. `molecule_or_readchain` groups by
source molecule. `exact_haploid_bin_pair` groups canonical binned endpoints.
`genomic_block` groups observations connected through endpoint blocks of 1, 5,
or 10 Mb. `chromosome_pair` groups the complete canonical chromosome pair.

Every split is deterministic for an explicit seed. The same molecule cannot be
split even when another grouping mode is selected. This invariant is tested by
a case in which molecule and bin-pair relations form a transitive chain.

## 058C: blind-visible scores versus truth

Seven fixed seeds were evaluated under each split. Two blind-visible scores
were frozen before truth evaluation:

```text
held-out pseudo-NLL = -mean log sum_s pi_s * L(contact | distance_s)
held-out pmax surprisal = -mean log max_s p_s
```

They are called pseudo-NLL/energy because the real P9016 pipeline does not have
a complete non-contact exposure denominator. Across seven split modes and
three truth targets:

| blind score | Spearman range | mean Spearman | Kendall range | mean top-quartile selection probability |
|---|---:|---:|---:|---:|
| pseudo-NLL | -0.7500 to 0.6071 | 0.0639 | -0.5238 to 0.4880 | 0.2857 |
| pmax surprisal | -0.7857 to 0.8469 | -0.1630 | -0.6190 to 0.6831 | 0.2857 |

An isolated favorable correlation is not reliable. For example, the 5 Mb
blocked pmax-surprisal score correlates strongly with parity (Spearman 0.8469)
but negatively with exact state (-0.2143) and only weakly with unordered
structure (0.1786). At 10 Mb, pseudo-NLL correlates -0.7500 with structure.
The target and split change the sign, so neither score defines a production
seed selector. Of the 49 new synthetic fits, 43 meet the optimizer convergence
criterion; failures are retained rather than silently filtered.

Top-k regret is provided in the machine table using the truth outcome only
after selection. No truth metric is combined with the blind score, and no
configuration is rerun based on regret.

## 058D: seed stability

Twenty-one pairs among seven fixed seeds yield:

| metric | mean | minimum | maximum |
|---|---:|---:|---:|
| same/cross agreement, raw gauge | 0.631587 | 0.562333 | 0.726333 |
| same/cross agreement, chromosome flips aligned | 0.644365 | 0.578000 | 0.713333 |
| aligned same/cross agreement, cis | 0.701244 | 0.594055 | 0.809942 |
| aligned same/cross agreement, trans | 0.521248 | 0.367089 | 0.620253 |
| exact-state agreement, raw gauge | 0.276683 | 0.101667 | 0.439000 |
| exact-state agreement after chromosome flips | 0.371587 | 0.235667 | 0.564667 |
| gauge-equivalent agreement gain | 0.094905 | 0.000000 | 0.300667 |
| distance-matrix Spearman | 0.413693 | 0.099413 | 0.646383 |
| mean chromosome structure Spearman | 0.787161 | 0.656823 | 0.858883 |
| posterior Jensen-Shannon divergence after alignment | 0.060225 | 0.036519 | 0.082315 |
| consistently callable same exact state | 0.002222 | 0.000000 | 0.015333 |

Chromosome-wise shapes can be similar while global placements and contact
states disagree. Gauge alignment explains only part of exact-state variation:
mean exact agreement remains 0.372 after optimal chromosome flips.

## Objective and symmetry invariants

The deterministic test suite flips one chromosome, several chromosomes, and
all chromosomes. Coordinates and incident p4 states are relabeled together.
The E-step, M-step, and total objective remain unchanged within floating-point
tolerance. Global copy-track alignment violates independent chromosome gauge;
hard minimum separation is gauge-invariant but can create unsupported homolog
separation. An unrestricted chromosome-pair four-state prior can violate the
factorization `r_cd = s_c * s_d` and become cycle-inconsistent.

The Z2 solver therefore reports components, cycle inconsistency, frustrated
edge fraction, per-spin margin, and unresolved chromosomes rather than forcing
a label on disconnected evidence.

## Intuitive example

Consider two raw contacts from chr1 bin 10 to chr2 bin 20 and one from chr1 bin
11 to chr2 bin 20. A raw split can place identical bin pair 10-20 in train and
test. Exact-pair blocking fixes that, but 11-20 remains a nearly identical local
constraint. A 5 Mb block split keeps all three together. If the held-out score
then ranks seed A above seed B while exact state ranks B above A, the score is
not a valid blind selector for exact phasing.

## Improvement direction and validation

Concrete problem: the existing held-out surface is not matched to the active
observation model and can leak loci. This matters because a lower score may
measure memorization of repeated genomic pairs. The correction is to use
genomic-block grouping and compute the configured E-step likelihood on held-out
contacts, with an explicit exposure model before calling it a full NLL.
Validation requires zero leakage in the manifest and consistent positive rank
correlation across independent seed ensembles, not one split.

Concrete problem: a single seed can appear best under one truth metric and
worst under another. This matters because production has no truth metric. For
example, 5 Mb pmax surprisal correlates positively with parity but negatively
with exact state in this run. Validation should preregister the blind score,
repeat on additional synthetic replicates and cells, and require low top-k
regret for exact, parity, and unordered structure simultaneously.

## Conclusion

The current row-level and exact-pair held-out designs are insufficiently
blocked. Genomic-block and chromosome-pair splits satisfy the audited leakage
invariants. On the synthetic diagnostic, the tested blind-visible scores do not
reliably select a better truth basin. That negative model-selection result must
not be promoted to a completed P9016 held-out study: repeated P9016 fits remain
unresolved. Production should report a preregistered seed ensemble and
stability rather than claim that one held-out-selected seed is optimal.
