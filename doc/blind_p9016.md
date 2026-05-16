# Blind P9016 1Mb Diploid Structural Decomposition

This note documents the current blind 1Mb P9016 runner behavior as implemented
in the C code and tests. Treat the C implementation as the source of truth.

## Interpretation

The blind pipeline is an unphased Hi-C driven diploid structural decomposition.
It estimates a latent copy-pair posterior (`p00`, `p01`, `p10`, `p11`) for each
binned contact. These copy labels are not guaranteed maternal/paternal genetic
haplotypes.

Copy labels have chromosome-level gauge symmetry: swapping copy 0 and copy 1
for a chromosome can describe the same structural solution. Reported copy labels
must therefore be interpreted as latent structural copies unless an explicitly
isolated eval-only mode is added.

The blind input adapter intentionally drops `phase0`, `phase1`, and
`hk_pair.phase[]`. New initialization, prior, training, and grid-search modes
must not use phase labels for training, model selection, thresholds, priors,
or initialization.

## Effective Modes

Initialization modes:

- `unphased_scaffold_split`: the preferred main mode. It builds a consensus
  haploid scaffold from unphased binned contacts with Hickit's FDG machinery,
  then splits the two copies with a chromosome-level deterministic direction.
- `random_diploid`: a deterministic-by-seed random diploid coordinate baseline.
  This is an ablation/null baseline and may be unstable.
- `random_haploid_split`: a deterministic-by-seed random haploid scaffold,
  split into two copies.
- `toy_split`: the old deterministic toy scaffold plus split. This remains
  available only as an explicit smoke/null mode, not as a biological default.

Prior modes:

- `uniform`: the old four-state uniform prior.
- `cis_inter_ratio`: a truth-free same-chromosome same-copy prior. It estimates
  an inter-chromosomal/background contact density over all possible inter-chrom
  bin pairs, estimates smoothed cis density by genomic distance over all
  possible same-chromosome bin pairs, and sets same-copy versus cross-copy mass
  from their ratio. The current smoother is a three-bin possible-weighted
  moving average, and `alpha_cross` is clamped to `[eps_prior, 0.5]`. Trans
  pairs remain uniform.

Training uncertainty modes:

- `rho_train_mode=constant`: old behavior. The schedule value is used for every
  binned pair.
- `rho_train_mode=entropy`: posterior confidence participates in training.
  The effective per-bpair weight is `schedule_rho * (1 - entropy / log(4))`,
  where the existing `rho_output` confidence factor is clamped to `[0, 1]`
  before multiplication by the schedule value.

Distance-scale modes:

- `d_scale_mode=raw_count`: old behavior. All four expanded state edges use
  `n_raw^(-1/3)`.
- `d_scale_mode=expected_count`: ablation. State `ab` uses
  `max(n_raw * p_ab, eps_count)^(-1/3)`, while stiffness remains
  `base_k * rho_train_bpair * p_ab`.

Same-bin handling:

- `same_bin_filter_enabled=1` filters canonical binned contacts with
  `bid0 == bid1` from posterior training and contact-force training.
- Same-bin bpair and raw posterior rows are retained in the schema as
  uniform-unknown rows (`p00..p11 = 0.25`, `pU = 1`, `rho_output = 0`).
- Manifest fields record `n_raw_same_bin_excluded`,
  `n_bpair_same_bin_excluded`, and `raw_posterior_same_bin_policy`.

Posterior refresh:

- After the final relax and final gauge stabilization, the scheduled loop runs
  a no-relax E-step on the final stabilized coordinates.
- Manifest fields record `posterior_refreshed_after_final_relax`,
  `posterior_refresh_temperature`, `posterior_refresh_prior_mode`, mean KL,
  top-state switch fraction, and mean `pU` before/after refresh.

Repulsion:

- Repulsion blocking is intentionally unchanged in this task. The manifest
  records `repulsion_blocking_mode=current_edge_blocking`.

Base stiffness:

- Binned blind pairs currently have effective `base_k=1.0`.
- Runner schedule structs may still carry a legacy/global value for historical
  compatibility, but manifests record it only as `legacy_base_k_unused`.
- Effective manifest fields are `base_k_mode`, `base_k_effective`,
  `base_k_min`, `base_k_mean`, `base_k_max`, and `base_k_n_nonfinite`.
- Initialization manifests use effective fields: `init_eps_effective`,
  `init_noise_scale_effective`, and `init_split_params_used`. For
  `random_diploid`, split params are recorded as unused/effectively zero.

## Why These Changes Matter

Initialization problem: a deterministic toy scaffold can look like a meaningful
biological starting point even though it is only a synthetic coordinate pattern.
This matters because EM-like reconstruction can preserve initialization bias.
Conceptual example: if every chromosome starts as two copies split around a toy
curve, a successful run may partly reflect that curve rather than the contacts.
Validation: compare `unphased_scaffold_split`, `random_diploid`, and `toy_split`
in the matrix runner while checking finite coordinates, manifest fields, and
output audit status.

Uncertainty-weighting problem: high-entropy contacts previously contributed
full total training force spread across four states. This matters because
ambiguous contacts can pull several incompatible copy-pair states at once.
Conceptual example: a uniform posterior `[0.25,0.25,0.25,0.25]` should convey
little confidence about which copy pair is close. Validation: tests assert that
uniform posteriors receive near-zero effective `rho_train_bpair` under entropy
mode, while confident posteriors retain the scheduled weight.

Final-posterior problem: the last written posterior could correspond to
pre-final-relax coordinates. This matters because the final coordinates are the
scientific object being interpreted. Conceptual example: if final relaxation
moves copies apart, stale posteriors may still report the pre-relax state.
Validation: loop diagnostics record the no-relax refresh, KL change, top-state
switch rate, and before/after mean `pU`; tests assert refresh fields are written
and that stored posterior values match recomputation from final coordinates.

Same-bin problem: a binned self-contact can otherwise create artificial
cross-copy homolog attraction. This matters most at 1Mb because many raw
contacts can collapse into the same bead. Conceptual example: contact `(i,i)`
should not become a force pulling `(i,copy0)` toward `(i,copy1)`. Validation:
tests assert same-bin bpairs produce no training wedges and auditors require
uniform-unknown posterior rows for same-bin outputs.

Cis-prior problem: a uniform prior is weak for short-range same-chromosome
contacts, but truth labels are forbidden. This matters because local cis
contacts are often much denser than trans/background contacts. Conceptual
example: if short-range cis density is 100x background, the prior can favor
same-copy states without asserting genetic truth. Validation: synthetic tests
check strong same-copy prior when cis density exceeds inter density, uniform
limit when densities match, trans-pair uniformity, and finite/clamped prior
metadata in the manifest/auditor path.

Expected-count distance-scale problem: low-posterior states previously inherited
the short target distance implied by the full raw count. This matters because an
unlikely state could still exert a short-distance target, even with low
stiffness. Conceptual example: with `n_raw=8` and `p_ab=0.125`, the expected
count target is longer than the target for `p_ab=0.875`. Validation: tests check
larger `d_scale` for low posterior states and finite behavior near zero.

Grid-search problem: new modes need reproducible coverage without selecting
configs from truth labels. This matters because ablations should compare model
assumptions, not hidden phase information. Conceptual example: a matrix can vary
`init_mode`, `prior_mode`, `rho_train_mode`, and `d_scale_mode` while recording
all effective parameters. Validation: `test_blind_grid_smoke` runs two configs,
writes distinct manifests and a summary TSV, and runs the output auditor.

## Output/Audit Expectations

Every full or matrix run should write:

- `p9016_full.bpair_posterior.tsv`
- `p9016_full.coords.tsv`
- `p9016_full.loop_diag.tsv`
- `p9016_full.manifest.tsv`
- optional `p9016_full.raw_posterior.tsv`

The auditor checks:

- no `phase0`, `phase1`, `truth`, or `oracle` strings in outputs
- manifest effective parameter fields
- prior estimator metadata, including smoothing and clamp fields
- bpair posterior probability ranges and sums
- same-bin uniform-unknown posterior rows
- final loop diagnostics, including posterior refresh fields
- coordinate row completeness and finite coordinates

The final scientific interpretation should remain conservative: this pipeline
reports latent copy-pair structural posteriors from unphased Hi-C contacts, not
validated genetic haplotype phasing.
