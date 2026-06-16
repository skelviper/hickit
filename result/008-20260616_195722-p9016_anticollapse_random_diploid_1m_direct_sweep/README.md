# 008-20260616_195722-p9016_anticollapse_random_diploid_1m_direct_sweep

This experiment repeats the 007 homolog-copy anti-collapse sweep from a
`random_diploid` initialization, using direct 1 Mb training.

- Training input: `/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- Resolution: 1 Mb
- Initialization: `random_diploid`, seed 17, scale 10
- Iterations: 100 EM iterations, 100 relax steps per iteration
- Backend: GPU
- Baseline settings preserved: uniform prior, constant rho_train, raw-count d_scale,
  temperature 1.0, FDG-flat E-step scoring, raw_expected_soft_all graph mode.
- Blind-training boundary: no SNP labels, phase labels, CHARM/3DG reference,
  multi-cell statistics, bulk reference, or extra omics were used during
  training.
- Eval-only inputs: `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz`
  plus phase columns from the P9016 pairs file. These are used only by
  `hickit/eval/evaluate_p9016_baseline.py` after training.
- Copy-swap policy: structure plots and distance metrics use per-chromosome
  cis distance-matrix Spearman correlation to select a geometry gauge. Contact
  identity headline metrics use eval-only whole-chromosome SNP cis-top1 oracle
  gauge. `copy0/copy1` are gauge labels, not maternal/paternal labels.

## Quantitative Summary

### Training-Side Anti-Collapse Diagnostics

| config | min_sep_unit | lambda_sep | final_min_sep | final_mean_sep | final_mean_entropy | final_mean_pU | final_sep_force_l1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline_off | 0 | 0 | 0.471240 | 3.476672 | 0.729346 | 0.526112 | 0.000000 |
| m1p5_l0p05 | 1.5 | 0.0500000007 | 0.688637 | 3.493239 | 0.729021 | 0.525878 | 8.108814 |
| m2_l0p05 | 2 | 0.0500000007 | 0.490184 | 3.496930 | 0.731348 | 0.527556 | 49.914402 |
| m2_l0p10 | 2 | 0.100000001 | 0.711376 | 3.430741 | 0.735878 | 0.530823 | 105.252190 |

Unlike the scaffold-init 007 run, the separation force remains nonzero in the
sep-on final diagnostics. However, it does not produce a clear contact-identity
improvement.

### Eval-Only Contact and Structure Metrics

Contact identity metrics use the eval-only whole-chromosome SNP cis-top1 gauge
for each chromosome. The distance metric uses the per-chromosome cis
distance-matrix Spearman geometry gauge. The p90 columns report accuracy and
recall among contacts whose maximum posterior state probability is at least
0.9.

| config | allele sep | top1 all | top1 cis | top1 trans | p90 acc all | p90 acc cis | p90 acc trans | p90 rec all | p90 rec cis | p90 rec trans | mean cis Spearman |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline_off | 3.476670 | 0.304748 | 0.337300 | 0.258159 | 0.316341 | 0.339563 | 0.268988 | 0.123391 | 0.150963 | 0.083930 | 0.146686 |
| m1p5_l0p05 | 3.493240 | 0.304670 | 0.334437 | 0.262070 | 0.307047 | 0.322229 | 0.276329 | 0.118354 | 0.141204 | 0.085653 | 0.143038 |
| m2_l0p05 | 3.496930 | 0.307253 | 0.337751 | 0.263605 | 0.313575 | 0.329434 | 0.281988 | 0.119303 | 0.141742 | 0.087188 | 0.148472 |
| m2_l0p10 | 3.430740 | 0.305996 | 0.338732 | 0.259146 | 0.316902 | 0.340209 | 0.270899 | 0.120743 | 0.146149 | 0.084381 | 0.159075 |

CHARM/3DG uniform-prior FDG probability baseline, computed from the reference
structure only: top1 all/cis/trans = 0.198369 / 0.221055 / 0.165901; p90
accuracy all/cis/trans = 0.050782 / 0.011865 / 0.088384; p90 recall
all/cis/trans = 0.014481 / 0.002824 / 0.031165.

### Mean Cis Distance Correlation Matrix

Mean per-chromosome Spearman correlation between CHARM/3DG cis distance
matrices and reconstruction cis distance matrices.

| config | CHARM c0 vs recon c0 | CHARM c0 vs recon c1 | CHARM c1 vs recon c0 | CHARM c1 vs recon c1 |
| --- | ---: | ---: | ---: | ---: |
| baseline_off | 0.216732 | 0.198047 | 0.194158 | 0.230897 |
| m1p5_l0p05 | 0.214423 | 0.199035 | 0.184019 | 0.230175 |
| m2_l0p05 | 0.218746 | 0.188697 | 0.182216 | 0.222038 |
| m2_l0p10 | 0.231740 | 0.208367 | 0.192229 | 0.237372 |

### Contact Distance Histograms

This eval-only plot uses raw contacts whose two ends both have SNP phase labels,
excludes same-bin contacts, and requires the phased endpoint coordinates to
exist in both CHARM/3DG and reconstruction. Distances are split into cis/trans.
The reconstruction copy labels use the eval-only whole-chromosome SNP contact
gauge.

| config | cis contacts | CHARM cis mean/median | recon cis mean/median | trans contacts | CHARM trans mean/median | recon trans mean/median |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline_off | 207084 | 1.051090 / 0.890171 | 1.786640 / 1.226610 | 143977 | 2.261140 / 2.182080 | 2.323990 / 2.175000 |
| m1p5_l0p05 | 207084 | 1.051090 / 0.890171 | 1.798010 / 1.228230 | 143977 | 2.261140 / 2.182080 | 2.326400 / 2.183280 |
| m2_l0p05 | 207084 | 1.051090 / 0.890171 | 1.795200 / 1.205340 | 143977 | 2.261140 / 2.182080 | 2.329530 / 2.180760 |
| m2_l0p10 | 207084 | 1.051090 / 0.890171 | 1.749950 / 1.183660 | 143977 | 2.261140 / 2.182080 | 2.334580 / 2.180680 |

## Files

- `outputs/`: runner outputs for each config.
- `eval/<config>/README.md`: per-config eval report with headline metrics.
- `eval/<config>/summary.tsv`: per-config machine-readable eval summary.
- `eval/<config>/contact_accuracy.tsv`: top1 and p90 contact accuracy/recall
  for all/cis/trans plus per-chromosome cis rows.
- `eval/<config>/cis_distance_correlation_matrix.tsv`: per-chromosome 2x2
  cis distance-matrix correlation table.
- `eval/<config>/contact_distance_distribution.tsv`: cis/trans 3D distance
  summaries for SNP-labeled contacts and posterior-top1 contacts.
- `eval/<config>/copy_separation.tsv`: copy0/copy1 separation by chromosome.
- `logs/sweep.log`: full training and audit log.
- `logs/eval_<config>.log`: eval logs for each config.
- `anticollapse_sweep_summary.tsv`: parsed manifest/loop/force diagnostics.
- `plots/`: root-level copies of each config's chr1 distance map, chr1
  highlighted 3D scatter, all-chromosome 3D scatter, and contact-distance
  histogram. The same plots also live under each `eval/<config>/plots/`
  directory.

## Interpretation

The random-diploid 1 Mb direct baseline is much weaker than the scaffold-init
1 Mb runs. Its top1 all/cis/trans is 0.305 / 0.337 / 0.258, and mean cis
distance Spearman is 0.147.

The sep-on configs do not give a robust improvement over the random-diploid
baseline. The best top1 all is `m2_l0p05` at 0.307253, a very small absolute
increase over 0.304748. The best mean cis Spearman is `m2_l0p10` at 0.159075,
but its top1 trans is similar to baseline and its p90 recall is lower than the
baseline on cis. The separation force changes final geometry diagnostics, but
the contact identity metrics remain near the random-init baseline.

This experiment supports the narrow conclusion that anti-collapse terms can
alter random-init geometry under 1 Mb direct training. It does not support a
model-improvement claim for posterior assignment or contact identity.
