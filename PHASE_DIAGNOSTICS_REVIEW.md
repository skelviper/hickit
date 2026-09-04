# Phase diagnostics 057-059: compact review package

This branch is a deliberately compact package for scientific and code review.
It is based on the public `origin/phase` commit
`2df02c60261995b5297132852513917d7ceb0b29` and contains only the diagnostic
implementation, tests, technical reports, summary-level result tables, and key
figures needed to assess experiments 057-059.

The package is not a new blind P9016 production result and is not a claim that
the C changes can be merged directly into the public `phase` branch. The full
local development history contains 163 unpublished commits between the public
base and the source tree used for the diagnostics. Publishing that history
would add unrelated code and approximately 86 MiB of Git objects, so it is
intentionally excluded here.

## Review boundaries

- P9016 training uses raw contact fields only. Phase labels and CHARM/3DG
  coordinates are used only in explicitly marked ORACLE analyses and final
  evaluation.
- `copy0` and `copy1` are chromosome-level gauge labels. They are not maternal
  and paternal identities.
- Exact `00/01/10/11`, same/cross parity, unordered chromosome structure, and
  chromosome-gauge recovery are reported as different targets.
- Experiment 058 is PARTIAL for real P9016: the leakage audit is complete, but
  repeated blocked P9016 fits and blind seed-selection validation were not run.
- Experiment 059 is synthetic. Its identifiability boundaries must not be
  reported as empirical P9016 limits.
- The strongest supported current claim is Claim D: stable blind per-cell
  diploid reconstruction has not been demonstrated with the current data and
  implementation.

## Provenance

| role | commit | meaning |
|---|---|---|
| public review base | `2df02c60261995b5297132852513917d7ceb0b29` | live `origin/phase` when this review branch was created |
| local C-patch base | `cfa590306e953cdc3a3ea35860b803a5669b895c` | unpublished source parent used to isolate the C diagnostic changes |
| corrected run source | `fb1eebf522578ed4c31ba49d4d7ffbc8e3bba33e` | source recorded by the corrected 057-059 run |
| corrected logical results | `9fdd34855b5363f94d81d36a37d7a385a5499f53` | regenerated result commit in the full local branch |
| formal-package commit | `38c1aa48e9a0a7a8c53ca6316b76ba11e52a0384` | local 119-121 package commit, not copied here |
| full local diagnostic head | `e94ca15787a5c55f5a1e70489eb6bb240ec57163` | source of the files curated into this branch |
| compact code commit | `bccb04c285af15dba216198b1a2f79dddfafad17` | Python modules, scripts, tests, and review-only C patch |

The C delta is preserved as
`patches/phase_diagnostics_c_core.patch`. It was generated from
`cfa590306e953cdc3a3ea35860b803a5669b895c` to the full diagnostic head and is
included for review of the 250-line core change. It is not asserted to apply to
the public base, because the required unpublished parent code is intentionally
not bundled. It reverse-applies cleanly to the full local diagnostic head.

## Included evidence

- `phase_diagnostics/`: centralized likelihood, mass, metrics, blocked-split,
  synthetic, structure, serialization, and Z2 modules.
- `scripts/run_phase_exp057.py`, `run_phase_exp058.py`, and
  `run_phase_exp059.py`: experiment drivers with explicit seeds and provenance.
- `scripts/run_phase_diagnostics.sh`: the full-tree entry point. Its pure
  Python portions are reviewable here; C targets require the full source tree.
- `tests/test_phase_diagnostics.py`: 47 deterministic Python tests covering
  endpoint ordering, p4 semantics, mass, symmetry, likelihood behavior,
  blocked splits, serialization, synthetic recovery, negative controls,
  packaging, and reporting invariants.
- `docs/phase_diagnostics_final_report.md`: the scientific synthesis.
- `docs/phase_model_audit.md` and `docs/p4_semantics.md`: execution graph,
  formulas, truth-read audit, symmetry risks, double-counting risks, and p4
  force semantics.
- `results/phase_diagnostics/057/`: complete compact 057 tables and figures.
- `results/phase_diagnostics/058/`: complete compact 058 tables and figures.
- `results/phase_diagnostics/059/replicate_summary.tsv`: all condition-level
  means, standard deviations, convergence rates, and finite-replicate counts.
- `results/phase_diagnostics/model_comparison/comparison_matrix.tsv`: the
  28-row comparison with evidence scope and truth boundary on every row.
- E-step curves, mass-conservation audit, experiment 001-056 lineage map, and
  test receipts under `results/phase_diagnostics/`.

## Intentionally omitted

- `results/phase_diagnostics/059/metrics_long.tsv`: a 35,382,711-byte
  replicate-by-metric table. Its scientific content is aggregated in
  `replicate_summary.tsv`, while `conditions.tsv` and `provenance.json` retain
  condition and run provenance. The omitted file SHA-256 is
  `a3188e748e4a87e4486b8b678808c399b543a4a323ef126d4647625b9f7ebe93`.
- `test_res/116-*` through `test_res/121-*`: superseded packages and duplicate
  formal aliases. The reports retain their original provenance references.
- Raw `.pairs.gz` inputs, CHARM/3DG reference inputs, reconstructed `.3dg`
  intermediates, logs beyond the compact test receipts, binaries, caches, and
  temporary work directories.
- The 163 unpublished commits preceding the diagnostic work.

## Key scientific findings

1. P9016 ORACLE-coordinate truth is the nearest geometric state for only
   18.87% of evaluated contacts and 16.59% of trans contacts. Monotone E-step
   top1 exactly follows the nearest geometry, so exact phase remains poor.
2. The historical FDG-flat score is not a monotone contact likelihood: it
   penalizes distances near zero and is flat over an interval.
3. With ORACLE contact phase, native cis+trans Hickit reaches mean
   chromosome-wise unordered distance Spearman `0.966 +/- 0.016`; the current
   joint M-step therefore works when phase is supplied.
4. In the well-specified 12,000-contact synthetic regime, blind exact and
   parity recovery remain near `0.25` and `0.50`, despite unordered structure
   Spearman `0.764-0.870`. This exposes an inference/optimization failure in
   addition to contact-level ambiguity.
5. Generator-conditional ORACLE-coordinate exact accuracy is only
   `0.282-0.417` in the tested high-depth slice; trans parity is
   `0.489-0.531`. The synthetic data also have a real observation-information
   limit.
6. No tested blind-visible held-out score reliably selects better synthetic
   seeds. Real P9016 held-out fitting remains unperformed.
7. The blind Z2 graph has no calibrated relative-gauge edges, so all components
   correctly remain unresolved and synchronization is a no-op.
8. Normalized molecule-level inference gives no stable parity gain, has only a
   small structure gain, costs 4.4-5.2 times more, and cannot be evaluated on
   P9016 because usable molecule IDs are absent.

## Verification

Run the self-contained tests from the repository root:

```bash
conda run -n analysis python -m pytest -q tests/test_phase_diagnostics.py
```

The compact branch passes 47 tests. The C test evidence in
`results/phase_diagnostics/test_c.log` was produced in the full source tree;
it is retained as provenance, not presented as a C build performed against the
public review base.

For review, start with this file, then read
`docs/phase_diagnostics_final_report.md`,
`results/phase_diagnostics/model_comparison/comparison_matrix.tsv`, and the
three experiment reports. Use the machine-readable tables for quantitative
claims rather than extracting values from prose or figures.
