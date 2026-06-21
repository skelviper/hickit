# 049 P9016 Trans Dscale Gamma 1Mb

This controlled blind-training experiment tests whether trans-specific posterior-count dscale gamma can fix the softall trans geometry bottleneck without changing posterior calculation, priors, rho, temperature, or candidate contact endpoints. The knob intentionally changes the trans count-to-distance semantics used by the M-step graph.

- full result root: `/mnt/ssd/zliu/phase3/test_res/049-20260620_103248-p9016_trans_dscale_gamma_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/049-20260620_103248-p9016_trans_dscale_gamma_1m`
- headline: `NO_PLUS_0P1`
- training inputs: approved P9016 raw pairs only
- eval inputs: P9016 SNP labels and CHARM/3DG are used only after training
- fixed knobs: posterior_count global gamma=1, lambda_copytrack=0.03, lambda_global_copytrack=0.003, no sep force, no callable anchor, no chr-pair M-step
- controlled knobs: trans_d_scale_posterior_gamma and trans_dscale_multiplier only
- success criterion: +0.1 full-denominator trans top1 over current best blind-safe standard trans, recorded as 0.358585
- strict matrix audit: `/mnt/ssd/zliu/phase3/test_res/049-20260620_103248-p9016_trans_dscale_gamma_1m/strict_matrix_audit.tsv` (`PASS`)

## Main Results

| config_name | seed_family | trans_d_scale_posterior_gamma | trans_dscale_multiplier | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | delta_trans_vs_same_seed_g1_td0p5 | delta_trans_vs_current_best047 | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | force_diag_total_contact_force_l1 | copytrack_frac_cos_lt_0 | target_plus_0p1_vs_current_best047_met |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_transgamma_g1_td0p5_seed71_best035_replay | seed71_noise0p05 | 1 | 0.5 | 0.483299 | 0.572482 | 0.355661 | 0 | -0.002924 | 0.561208 | 0.748009 | 0.692783058 | 0.499737322 | 133210.594 | 0.00344431688 | 0 |
| p9016_transgamma_g0p75_td0p75_seed71_best035 | seed71_noise0p05 | 0.75 | 0.75 | 0.479519 | 0.567794 | 0.353181 | -0.00248 | -0.005404 | 0.556554 | 0.789754 | 0.657840312 | 0.474531442 | 111173.25 | 0.003061615 | 0 |
| p9016_transgamma_g1_td0p5_seed17_best035_replay | seed17_noise0 | 1 | 0.5 | 0.481133 | 0.570929 | 0.352619 | 0 | -0.005966 | 0.55256 | 0.76939 | 0.697600842 | 0.503212631 | 132769.453 | 0.0015308075 | 0 |
| p9016_transgamma_g0p5_td0p75_seed71_best035 | seed71_noise0p05 | 0.5 | 0.75 | 0.474745 | 0.56116 | 0.35107 | -0.004591 | -0.007515 | 0.564396 | 0.780464 | 0.661346436 | 0.477060616 | 112149.461 | 0.003061615 | 0 |
| p9016_transgamma_g0p25_td0p75_seed71_best035 | seed71_noise0p05 | 0.25 | 0.75 | 0.480753 | 0.574433 | 0.34668 | -0.008981 | -0.011905 | 0.561423 | 0.765422 | 0.665140331 | 0.479797333 | 113820.234 | 0.00267891313 | 0 |
| p9016_transgamma_g0_td1_seed71_best035 | seed71_noise0p05 | 0 | 1 | 0.478094 | 0.570439 | 0.34593 | -0.009731 | -0.012655 | 0.547803 | 0.779605 | 0.628072739 | 0.45305872 | 93922.9688 | 0.003061615 | 0 |
| p9016_transgamma_g1_td0p75_seed71_best035 | seed71_noise0p05 | 1 | 0.75 | 0.469323 | 0.557268 | 0.343458 | -0.012203 | -0.015127 | 0.543635 | 0.78117 | 0.648466766 | 0.467769891 | 110257.477 | 0.00344431688 | 0 |
| p9016_transgamma_g0_td0p75_seed71_best035 | seed71_noise0p05 | 0 | 0.75 | 0.475577 | 0.568396 | 0.342735 | -0.012926 | -0.01585 | 0.540538 | 0.778882 | 0.662590384 | 0.477957934 | 115804.867 | 0.00229621125 | 0 |
| p9016_transgamma_g0p25_td1_seed71_best035 | seed71_noise0p05 | 0.25 | 1 | 0.465775 | 0.559025 | 0.332317 | -0.023344 | -0.026268 | 0.537912 | 0.771652 | 0.619978607 | 0.447219998 | 92534.7969 | 0.0015308075 | 0 |
| p9016_transgamma_g1_td1_seed17_best035 | seed17_noise0 | 1 | 1 | 0.470751 | 0.568794 | 0.330435 | -0.022184 | -0.02815 | 0.538565 | 0.794001 | 0.612831295 | 0.442064345 | 87555.625 | 0.0015308075 | 0 |
| p9016_transgamma_g0p5_td0p75_seed17_best035 | seed17_noise0 | 0.5 | 0.75 | 0.466489 | 0.563621 | 0.327476 | -0.025143 | -0.031109 | 0.532557 | 0.783198 | 0.661525667 | 0.477189898 | 113237.367 | 0.00114810563 | 0 |
| p9016_transgamma_g0p75_td0p75_seed17_best035 | seed17_noise0 | 0.75 | 0.75 | 0.467043 | 0.565208 | 0.326552 | -0.026067 | -0.032033 | 0.528987 | 0.779458 | 0.658691585 | 0.475145549 | 110810.016 | 0.00191350938 | 0 |
| p9016_transgamma_g1_td0p75_seed17_best035 | seed17_noise0 | 1 | 0.75 | 0.462118 | 0.558384 | 0.324343 | -0.028276 | -0.034242 | 0.52398 | 0.779609 | 0.653337836 | 0.471283615 | 110132.5 | 0.0015308075 | 0 |
| p9016_transgamma_g0p25_td0p75_seed17_best035 | seed17_noise0 | 0.25 | 0.75 | 0.462235 | 0.559059 | 0.323663 | -0.028956 | -0.034922 | 0.531627 | 0.765965 | 0.672164202 | 0.484863997 | 113892.18 | 0.00191350938 | 0 |
| p9016_transgamma_g0_td0p75_seed17_best035 | seed17_noise0 | 0 | 0.75 | 0.463321 | 0.561757 | 0.32244 | -0.030179 | -0.036145 | 0.529341 | 0.763758 | 0.67823112 | 0.489240348 | 115525.977 | 0.00114810563 | 0 |
| p9016_transgamma_g0_td1_seed17_best035 | seed17_noise0 | 0 | 1 | 0.470171 | 0.575748 | 0.319072 | -0.033547 | -0.039513 | 0.519798 | 0.80583 | 0.62447536 | 0.450463742 | 90809.4922 | 0.00114810563 | 0 |
| p9016_transgamma_g0p75_td1_seed17_best035 | seed17_noise0 | 0.75 | 1 | 0.461055 | 0.562587 | 0.315745 | -0.036874 | -0.04284 | 0.533217 | 0.795562 | 0.621849716 | 0.448569775 | 88398.5938 | 0.00229621125 | 0 |
| p9016_transgamma_g0p5_td1_seed71_best035 | seed71_noise0p05 | 0.5 | 1 | 0.453493 | 0.551275 | 0.31355 | -0.042111 | -0.045035 | 0.538454 | 0.77753 | 0.633038104 | 0.456640482 | 89927.5547 | 0.0015308075 | 0 |
| p9016_transgamma_g0p75_td1_seed71_best035 | seed71_noise0p05 | 0.75 | 1 | 0.453522 | 0.551561 | 0.31321 | -0.042451 | -0.045375 | 0.521215 | 0.770784 | 0.599611759 | 0.432528436 | 91995.7656 | 0.0015308075 | 0 |
| p9016_transgamma_g0p5_td1_seed17_best035 | seed17_noise0 | 0.5 | 1 | 0.452865 | 0.55176 | 0.311328 | -0.041291 | -0.047257 | 0.523146 | 0.795267 | 0.617253423 | 0.445254207 | 89930.0391 | 0.000382701875 | 0 |
| p9016_transgamma_g0p25_td1_seed17_best035 | seed17_noise0 | 0.25 | 1 | 0.451485 | 0.556424 | 0.301298 | -0.051321 | -0.057287 | 0.504713 | 0.800419 | 0.625678539 | 0.451331675 | 91594.8359 | 0.00191350938 | 0 |
| p9016_transgamma_g1_td1_seed71_best035 | seed71_noise0p05 | 1 | 1 | 0.434212 | 0.543961 | 0.277141 | -0.07852 | -0.081444 | 0.491669 | 0.791948 | 0.594022334 | 0.42849651 | 90070.2344 | 0.00114810563 | 0 |

## Interpretation Boundary

A positive result would mean the current trans gap is sensitive to count-to-distance semantics in softall, not merely posterior assignment. A negative result means the largest known dscale lever is insufficient; then the remaining problem is more likely raw-pair identifiability or missing data linkage rather than another posterior self-training rule.

## Post-Run Audit Note

The initial run completed with the runner snapshot stored in `scripts_snapshot/`. After completion, `scripts/audit_p9016_049_strict.py` was run as a post-hoc result-integrity audit. It checks the exact 22-config matrix, finite headline metrics, per-config manifest key/value invariants, and one replay baseline per seed family. It passed.
