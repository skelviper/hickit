# 010-20260617_001012-p9016_dscale_sep_sweep

This experiment compares the active P9016 `softall` baseline under two
distance-scale interpretations:

- `raw_count`: contact `d_scale` follows the observed raw/effective contact
  count, preserving the Hickit-style count-to-distance assumption.
- `expected_count`: contact `d_scale` follows the posterior-weighted expected
  count for each latent copy-state edge.

The sweep also checks whether homolog-copy separation changes the behavior
under the two `d_scale` modes.

- Training input: `/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- Resolution: 1 Mb
- Initialization: `unphased_scaffold_split`, seed 17. The common-init configs
  use `init_eps=1.0` and `init_noise_scale=0.05`; the legacy raw-count control
  uses `init_eps=0.5` and `init_noise_scale=0`.
- Iterations: 100 EM iterations, 100 relax steps per iteration
- Backend: GPU
- Baseline settings: `softall`, `raw_expected_soft_all` M-step graph, uniform
  prior, constant rho_train, temperature 1.0, FDG-flat E-step scoring.
- Blind-training boundary: no SNP labels, phase labels, CHARM/3DG reference,
  multi-cell statistics, bulk reference, or extra omics were used during
  training. All configs have `uses_phase_labels=0`, `uses_charm_or_reference=0`,
  and `input_contact_source=raw_pairs` in the training manifests.
- Eval-only inputs: `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz`
  plus phase columns from the P9016 pairs file. These are used only by
  `hickit/eval/evaluate_p9016_baseline.py` after training.
- Copy-swap policy: structure plots and distance metrics use per-chromosome
  cis distance-matrix Spearman correlation to select a geometry gauge. Contact
  identity headline metrics use eval-only whole-chromosome SNP cis-top1 oracle
  gauge. `copy0/copy1` are gauge labels, not maternal/paternal labels.

All 7 training configs completed with `status=OK`, and all 7 eval reports were
generated.

## Quantitative Summary

Contact identity metrics use the eval-only whole-chromosome SNP cis-top1 gauge
for each chromosome. The distance metric uses the per-chromosome cis
distance-matrix Spearman geometry gauge. The p90 columns report accuracy and
recall among contacts whose maximum posterior state probability is at least
0.9.

| config | d_scale | min sep | lambda sep | top1 all | top1 cis | top1 trans | same/cross cis | p90 acc all | p90 rec all | mean cis Spearman | mean copy sep |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| expected_sep_off | expected_count | 0 | 0 | 0.441331 | 0.549397 | 0.286671 | 0.974958 | 0.452380 | 0.142400 | 0.816514 | 5.589254 |
| expected_lsep0.5 | expected_count | 1.5 | 0.5 | 0.436392 | 0.549440 | 0.274599 | 0.974958 | 0.441071 | 0.140015 | 0.815823 | 5.589593 |
| expected_lsep1 | expected_count | 1.5 | 1 | 0.436575 | 0.549518 | 0.274933 | 0.974958 | 0.439806 | 0.139635 | 0.815763 | 5.590745 |
| legacy_raw_sep_off | raw_count | 0 | 0 | 0.417357 | 0.498886 | 0.300673 | 0.869711 | 0.440056 | 0.070279 | 0.577551 | 3.481423 |
| raw_sep_off_common_init | raw_count | 0 | 0 | 0.449325 | 0.549596 | 0.305820 | 0.972294 | 0.413743 | 0.081035 | 0.804125 | 5.228669 |
| raw_lsep0.5 | raw_count | 1.5 | 0.5 | 0.445168 | 0.550387 | 0.294582 | 0.973231 | 0.397747 | 0.077967 | 0.812928 | 5.260802 |
| raw_lsep1 | raw_count | 1.5 | 1 | 0.455033 | 0.543990 | 0.327719 | 0.973973 | 0.449996 | 0.086860 | 0.825276 | 5.247381 |

CHARM/3DG uniform-prior FDG probability baseline, computed from the reference
structure only: top1 all/cis/trans = 0.198369 / 0.221055 / 0.165901; p90
accuracy all/cis/trans = 0.050782 / 0.011865 / 0.088384; p90 recall
all/cis/trans = 0.014481 / 0.002824 / 0.031165.

## Training-Side Diagnostics

| config | final entropy | final pU | final min sep | final mean sep | sep force L1 | contact force L1 | sep/contact force |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| expected_sep_off | 0.618535 | 0.446178 | 2.957557 | 5.589254 | 0.000000 | 87420.679700 | 0.000000 |
| expected_lsep0.5 | 0.618755 | 0.446338 | 2.947668 | 5.589593 | 0.000000 | 87298.328100 | 0.000000 |
| expected_lsep1 | 0.618786 | 0.446360 | 2.950811 | 5.590745 | 0.000000 | 87306.851600 | 0.000000 |
| legacy_raw_sep_off | 0.806347 | 0.581656 | 0.068249 | 3.481423 | 0.000000 | 108387.781000 | 0.000000 |
| raw_sep_off_common_init | 0.656981 | 0.473912 | 0.160607 | 5.228669 | 0.000000 | 106458.703000 | 0.000000 |
| raw_lsep0.5 | 0.656192 | 0.473343 | 0.371016 | 5.260802 | 80.988762 | 106396.781000 | 0.000761 |
| raw_lsep1 | 0.656317 | 0.473432 | 1.001274 | 5.247381 | 15.147427 | 106398.594000 | 0.000142 |

## Interpretation

The common-init `raw_count` configs are stronger than the `expected_count`
configs on four-state contact identity, especially trans identity. The best
overall row is `raw_lsep1`: top1 all/cis/trans = 0.455033 / 0.543990 /
0.327719, with mean cis distance Spearman 0.825276. `expected_count` gives
higher p90 all-contact recall than the common-init `raw_count` rows, but it
does not improve top1 identity and is weaker on trans.

The legacy raw-count control is not an apples-to-apples baseline for the
common-init comparison because it used `init_eps=0.5` and no initialization
noise. It is kept here as a continuity check against earlier scaffold-init
runs.

This experiment supports keeping `raw_count` as the conservative Hickit-faithful
baseline. `expected_count` is a useful ablation, but these results do not
support replacing the baseline with it.

## Files

- `dscale_sep_sweep_summary.tsv`: combined manifest, training-diagnostic, and
  eval metric table for all 7 configs.
- `eval/<config>/README.md`: per-config eval report with standard quantitative
  metrics.
- `eval/<config>/summary.tsv`: per-config machine-readable eval summary.
- `eval/<config>/contact_accuracy.tsv`: top1 and p90 contact accuracy/recall
  for all/cis/trans plus per-chromosome cis rows.
- `eval/<config>/cis_distance_correlation_matrix.tsv`: per-chromosome 2x2
  cis distance-matrix correlation table.
- `eval/<config>/contact_distance_distribution.tsv`: cis/trans 3D distance
  summaries for SNP-labeled contacts and posterior-top1 contacts.
- `eval/<config>/copy_separation.tsv`: copy0/copy1 separation by chromosome.
- `eval/<config>/per_chrom_volume.tsv`: per-chromosome volumes for CHARM/3DG
  and reconstruction.
- `plots/`: root-level symlinks to each config's chr1 distance map, chr1
  highlighted 3D scatter, all-chromosome 3D scatter, and contact-distance
  histogram. The same plots also live under each `eval/<config>/plots/`
  directory.
- `logs/`: training, audit, and eval logs.
