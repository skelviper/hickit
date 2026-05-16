# Blind P9016 Multiresolution Prototype Design

This document describes the first blind/unphased diploid structural
decomposition prototype for:

```text
4Mb coarse reconstruction -> 1Mb lifted initialization -> anchored 1Mb refinement
```

The design follows standard Hickit where possible and records blind-specific
divergences explicitly.

## Scientific Guardrails

- This is blind/unphased diploid structural decomposition, not genetic haplotype
  phasing.
- `copy0` and `copy1` labels are gauge-symmetric internal labels. They must not
  be interpreted as maternal/paternal.
- P9016 pairs may contain `phase0`/`phase1`, but phase fields are forbidden for
  training, initialization, prior construction, seed/config selection, stopping,
  and model selection.
- CHARM 3DG is post-hoc evaluation only. It must not be used for training,
  initialization, seed/config selection, stopping, priors, or model selection.

## Problem Being Tested

Concrete problem: current high-cis 1Mb blind runs fit cis structure better than
early runs but show severe chromosome scale compression and still fail absolute
trans/centroid acceptability flags.

Why it matters: if the 1Mb objective is locally cis-dense and scale-compressive,
starting directly at 1Mb may destroy global chromosome layout before trans or
centroid information can organize the structure.

Conceptual example: a 4Mb bead summarizes several 1Mb beads. At 4Mb, fewer
within-chromosome edges compete with trans/global constraints. If a 4Mb coarse
scaffold preserves broader chromosome scale, 1Mb children can be initialized near
their coarse parents and refined locally instead of discovering both global and
local layout from scratch.

Experimental validation: compare direct 1Mb baselines to 4Mb->1Mb refinements
using both truth-free heldout contact metrics and post-hoc CHARM evaluation. A
valid improvement must preserve or improve cis quality without worsening scale,
homolog separation, or absolute trans/centroid flags.

## Standard Hickit Components Reused

The prototype should reuse these standard Hickit concepts:

- Bmap construction through `hk_bmap_gen()`.
- Explicit resolution/bin size through the `size` argument passed to
  `hk_bmap_gen()`.
- Merge semantics through the final `bmap_skip_merge_flag` argument.
- Standard unphased FDG haploid scaffold through `hk_fdg()` when
  `standard_hickit_unphased_fdg_4mb` is requested.
- Standard coarse-to-fine intuition from repeated `-b` runs: previous-resolution
  coordinates seed the next resolution.
- Standard distance normalization form:
  `r = (distance / unit) / d_scale`.

## Blind-Specific Components

The prototype remains blind-specific in these places:

- Raw P9016 pairs are converted to `hk_blind_pair`, which stores only
  chromosome, position, and strand fields.
- Latent copy-pair posterior is computed over four copy states.
- Weighted M-step contact edges are expanded from each bpair posterior.
- Copy labels remain gauge-only.
- Prior, rho, d_scale, homolog separation, and anchor terms operate without
  phase labels or reference coordinates.
- Current baseline unit remains `1.0` unless a run explicitly records a different
  scale-normalization choice.

## Part 1: Resolution-Aware Blind Bmap

Concrete problem: current P9016 blind runners hard-code 1Mb bmaps, preventing
4Mb coarse reconstruction.

Why it matters: without a parameterized bmap, there is no controlled way to
compare direct 1Mb training to 4Mb coarse training or to lift 4Mb coordinates
into 1Mb initialization.

Conceptual example: `chr1:0-4Mb` should become one coarse bead at 4Mb and four
fine beads at 1Mb, except at chromosome ends where the child count may be fewer.

Validation:

- 1Mb mode reproduces existing bmap/bpair counts.
- 4Mb mode yields fewer haploid and diploid beads.
- Same-bin filtering is valid at both resolutions.
- Manifest records `bin_size_bp`, resolution label, merge mode, and skip-merge
  argument.
- Bmap summary records bead intervals and child counts where applicable.

## Part 2: 4Mb-To-1Mb Coarse/Fine Map

Concrete problem: the lift needs a true parent-child map. Assuming exactly four
children per coarse bead is incorrect at chromosome ends and under merged-bead
settings.

Why it matters: wrong mapping would silently mis-anchor 1Mb beads, corrupting
the biological interpretation of scale and residual diagnostics.

Conceptual example: if a chromosome has 10.5Mb, the final 4Mb coarse bead may
cover fewer than four 1Mb children. The map must derive children by interval
containment/overlap, not by integer division alone.

Validation:

- Every fine bead has exactly one coarse parent.
- Child order follows genomic coordinate order.
- Terminal beads with fewer children are handled.
- No phase labels are present in the mapping.

## Part 3: Resolution-Aware d_scale

Concrete problem: raw 4Mb contact counts aggregate multiple 1Mb child-pair
opportunities. Applying `n_raw^(-1/3)` directly may make coarse high-exposure
contacts artificially short.

Why it matters: recent diagnostics implicate scale/contact-balance issues. A
4Mb run must avoid introducing a new coarse-resolution count-density artifact.

Conceptual example: if a 4Mb-by-4Mb bpair represents 16 possible 1Mb child
pairs, `n_raw=160` should be comparable to density `10` per child-pair, not to a
single 1Mb edge with count 160.

Validation:

- `raw_count` exactly reproduces old `n_raw^(-1/3)` behavior.
- `density_normalized_raw_count` uses `n_density = n_raw / exposure`.
- Terminal coarse beads use their actual child counts.
- No NaN/inf occurs for low or zero effective counts.
- Diagnostics report d_scale distributions by cis/trans, genomic-distance
  stratum, raw-count bin, and posterior-confidence bin.

## Part 4: 4Mb Initialization

Concrete problem: a coarse stage needs initialization modes that are either
standard Hickit-derived or explicitly labeled random/null controls.

Why it matters: calling an ad hoc scaffold “Hickit” would make mechanism
interpretation invalid.

Conceptual example: `standard_hickit_unphased_fdg_4mb` should call standard
`hk_fdg()` on unphased 4Mb contacts, then use the existing blind split function
to make copy0/copy1 internal labels.

Validation:

- Same seed is reproducible.
- Different seed changes random modes.
- Standard FDG path does not read `hk_pair.phase`.
- `toy_split` is used only when explicitly requested and is labeled as null/smoke.

## Part 5: 4Mb Coarse Scan

Concrete problem: seed/config selection cannot use CHARM because that would leak
reference geometry into a blind training workflow.

Why it matters: otherwise a multiresolution pipeline could look good only
because post-hoc reference information selected the coarse seed.

Conceptual example: choose between two random 4Mb seeds by heldout contact
residual and finite/non-collapsed sanity, not by CHARM trans correlation.

Validation:

- Deterministic train/heldout split by bpair-key hash.
- Selection report includes only truth-free metrics.
- CHARM metrics, if computed, appear only in post-hoc eval summary.
- Regression review verifies no CHARM metric is read by seed-selection code.

## Part 6: Lift From 4Mb To 1Mb

Concrete problem: a lifted 1Mb initialization should inherit coarse global
layout while giving children non-identical local positions.

Why it matters: identical child coordinates produce immediate local degeneracy;
arbitrary random offsets can erase the coarse layout being tested.

Conceptual example: for a coarse bead with four 1Mb children, place the children
along the local coarse chromosome axis around the parent coordinate, with offset
scale tied to neighboring coarse distances or recorded unit scale.

Validation:

- Every fine bead is mapped.
- Child offsets are deterministic for a fixed parent config and seed.
- Chromosome ends work.
- Manifest records parent config, map path, offset mode, offset scale, and the
  gauge-only status of inherited copy labels.

## Part 7: Anchored 1Mb Refinement

Concrete problem: even if 4Mb layout is reasonable, the 1Mb refinement objective
may collapse scale during release.

Why it matters: the experiment needs to distinguish “4Mb cannot make a useful
global scaffold” from “4Mb scaffold is useful but 1Mb refinement destroys it.”

Conceptual example: for each 4Mb parent/copy, a weak spring keeps the centroid
of its 1Mb children near the parent 4Mb coordinate. Local 1Mb geometry remains
free within the parent.

Validation:

- Anchor disabled reproduces old non-lift behavior.
- On a toy example, anchor force pulls child centroid toward the parent.
- Manifest records anchor schedule and strength.
- Diagnostics include anchor energy/force and final deviation from parent
  layout.
- Reject any run that improves trans Spearman only by destroying cis or scale.

## Required Output Families

- `mr4mb_scan_summary.tsv`
- `mr4mb_truth_free_selection.tsv`
- `mr4mb_eval_summary.tsv`
- `mr4mb_seed_selection_report.md`
- `coarse_to_fine_map.tsv`
- `lifted_1mb_init.coords.tsv.gz`
- `lift_sanity.tsv`
- `mr4mb_to_1mb_refinement_report.md`

## Interpretation Rules

Use these terms:

- blind/unphased diploid structural decomposition
- copy labels remain gauge-symmetric
- 4Mb coarse scaffold
- truth-free seed selection
- post-hoc CHARM evaluation

Avoid these terms unless absolute acceptability flags pass and the statement is
carefully qualified:

- maternal/paternal recovered
- genetic haplotype phasing
- trans-good
- successful trans reconstruction
