# 045 P9016 Raw-Heldout Selection Audit 1Mb

This controlled blind-training diagnostic tests whether a deterministic heldout split of raw P9016 contacts provides a reference-free selector for better trans basins. Training uses only the training split of raw pairs; phase labels and CHARM/3DG are eval-only.

- full result root: `/mnt/ssd/zliu/phase3/test_res/045-20260620_055343-p9016_raw_heldout_selection_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/045-20260620_055343-p9016_raw_heldout_selection_1m`
- headline: `HELDOUT_SELECTOR_EVAL_BEST`
- training pairs: `/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- eval pairs: `/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- eval CHARM/3DG: `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz`
- bin size bp: `1000000`
- heldout fraction / seed: `0.1` / `17`
- success criterion: a blind-visible heldout selector chooses a config with full-denominator trans top1 at least +0.1 over `p9016_heldout_td1_seed17_noise0`
- heldout split: deterministic hash of raw contact fields; no SNP, CHARM/3DG, or eval labels enter the split or training
- training boundary: blind training uses only raw P9016 pairs; SNP phase labels and CHARM/3DG are used only after training by the standard eval wrapper
- gauge policy: `copy0`/`copy1` are arbitrary homolog labels, not maternal/paternal labels; eval metrics use the established swap-aware/per-chromosome best orientation where applicable
- plots: per-config eval plots, including contact-distance and structure diagnostics, are under `eval/<config>/plots/` in the full result root
- caveat: heldout diagnostics are binned raw-contact fit diagnostics, not SNP truth; eval columns are post hoc annotations

## Main Results

| config_name | basin_family | seed | heldout_mean_expected_energy | heldout_mean_entropy | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | delta_trans_vs_heldout_td1_seed17 | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | sep_p05 | copytrack_frac_cos_lt_0 | target_plus_0p1_vs_heldout_td1_met |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_heldout_chrpair_m1_p1_w0_seed23_noise0p05 | chrpair_m1_p1_w0 | 23 | 0.45075330268595293 | 0.54548101704750362 | 0.489025 | 0.586062 | 0.34251 | 0.054642 | 0.805144 | 0.473873705 | 0.341827631 | 3.41710949 | 0.000382701875 | 0 |
| p9016_heldout_best035_seed31_noise0p05 | best035 | 31 | 0.59526042929892342 | 0.66927190809234183 | 0.473447 | 0.561044 | 0.341186 | 0.053318 | 0.693584 | 0.678898871 | 0.489722013 | 2.2382021 | 0.0045924225 | 0 |
| p9016_heldout_best035_seed97_noise0p05 | best035 | 97 | 0.61959068100484782 | 0.67449606380039384 | 0.477687 | 0.572064 | 0.335188 | 0.04732 | 0.718532 | 0.694950283 | 0.501300633 | 2.39567804 | 0.00650593188 | 0 |
| p9016_heldout_chrpair_m0p5_p4_w10_seed47_noise0p05 | chrpair_m0p5_p4_w10 | 47 | 0.5390414457424364 | 0.64324284056303105 | 0.470621 | 0.560946 | 0.334239 | 0.046371 | 0.799091 | 0.643758833 | 0.464373827 | 2.85796452 | 0 | 0 |
| p9016_heldout_chrpair_m1_p1_w0_seed59_noise0p05 | chrpair_m1_p1_w0 | 59 | 0.51267231415680636 | 0.62563836807537243 | 0.481959 | 0.581485 | 0.331685 | 0.043817 | 0.805928 | 0.618219733 | 0.445951253 | 3.23213959 | 0.000382701875 | 0 |
| p9016_heldout_chrpair_m0p5_p4_w10_seed59_noise0p05 | chrpair_m0p5_p4_w10 | 59 | 0.51890010158479383 | 0.63464712161856696 | 0.477552 | 0.575207 | 0.330103 | 0.042235 | 0.810218 | 0.640726507 | 0.462186456 | 3.17356396 | 0.00114810563 | 0 |
| p9016_heldout_chrpair_m0p5_p4_w10_seed31_noise0p05 | chrpair_m0p5_p4_w10 | 31 | 0.50425303332815574 | 0.61766831841010428 | 0.486164 | 0.589605 | 0.329978 | 0.04211 | 0.757152 | 0.610025346 | 0.44004029 | 3.21430326 | 0.00076540375 | 0 |
| p9016_heldout_chrpair_m1_p1_w0_seed97_noise0p05 | chrpair_m1_p1_w0 | 97 | 0.51102187324140147 | 0.6262064082353217 | 0.481648 | 0.582684 | 0.329094 | 0.041226 | 0.819105 | 0.602039695 | 0.434279829 | 2.99255848 | 0 | 0 |
| p9016_heldout_chrpair_m0p5_p4_w10_seed97_noise0p05 | chrpair_m0p5_p4_w10 | 97 | 0.53617345591630927 | 0.64508872354029345 | 0.488644 | 0.595717 | 0.326975 | 0.039107 | 0.80719 | 0.643081725 | 0.463885427 | 2.98208499 | 0.00267891313 | 0 |
| p9016_heldout_chrpair_m0p5_p4_w10_seed23_noise0p05 | chrpair_m0p5_p4_w10 | 23 | 0.5097347072188273 | 0.62494118390202236 | 0.470084 | 0.56757 | 0.322891 | 0.035023 | 0.795929 | 0.613536239 | 0.442572862 | 2.91190004 | 0.00076540375 | 0 |
| p9016_heldout_best035_seed17_noise0 | best035 | 17 | 0.58999189845420186 | 0.66325748308069032 | 0.457128 | 0.546646 | 0.321963 | 0.034095 | 0.693877 | 0.671137571 | 0.484123409 | 2.33239985 | 0.0045924225 | 0 |
| p9016_heldout_chrpair_m0p5_p4_w10_seed83_noise0p05 | chrpair_m0p5_p4_w10 | 83 | 0.53353872625536014 | 0.63783460928764113 | 0.467428 | 0.564489 | 0.320874 | 0.033006 | 0.737694 | 0.628151298 | 0.453115344 | 2.76607251 | 0.00344431688 | 0 |
| p9016_heldout_chrpair_m1_p1_w0_seed31_noise0p05 | chrpair_m1_p1_w0 | 31 | 0.53570369858608535 | 0.63833216178583074 | 0.476781 | 0.581461 | 0.318725 | 0.030857 | 0.772307 | 0.632941186 | 0.456570566 | 2.74594688 | 0.00267891313 | 0 |
| p9016_heldout_best035_seed47_noise0p05 | best035 | 47 | 0.58076698583858744 | 0.65042982818394557 | 0.453204 | 0.543059 | 0.317533 | 0.029665 | 0.708303 | 0.65909791 | 0.475438625 | 2.60854506 | 0.00382701875 | 0 |
| p9016_heldout_chrpair_m1_p1_w0_seed83_noise0p05 | chrpair_m1_p1_w0 | 83 | 0.51033453292041553 | 0.63216148904706171 | 0.468926 | 0.569573 | 0.316959 | 0.029091 | 0.778843 | 0.614757001 | 0.443453461 | 2.94251776 | 0 | 0 |
| p9016_heldout_best035_seed59_noise0p05 | best035 | 59 | 0.60714163939309262 | 0.6733167905937083 | 0.461499 | 0.558075 | 0.315679 | 0.027811 | 0.713955 | 0.685909152 | 0.494778872 | 2.26173711 | 0.00535782625 | 0 |
| p9016_heldout_best035_seed23_noise0p05 | best035 | 23 | 0.57748061224455427 | 0.65421077934126703 | 0.455052 | 0.547913 | 0.31484 | 0.026972 | 0.738765 | 0.661873639 | 0.477440923 | 2.61536121 | 0.00382701875 | 0 |
| p9016_heldout_chrpair_m1_p1_w0_seed17_noise0p05 | chrpair_m1_p1_w0 | 17 | 0.48350053546931049 | 0.60082599836750628 | 0.468428 | 0.570709 | 0.313993 | 0.026125 | 0.808542 | 0.585657537 | 0.422462612 | 3.32729292 | 0.000382701875 | 0 |
| p9016_heldout_chrpair_m0p5_p4_w10_seed17_noise0p05 | chrpair_m0p5_p4_w10 | 17 | 0.49374818259914743 | 0.60716433860653651 | 0.472286 | 0.577786 | 0.312993 | 0.025125 | 0.774067 | 0.597333968 | 0.430885375 | 3.24738002 | 0 | 0 |
| p9016_heldout_best035_seed83_noise0p05 | best035 | 83 | 0.61696741834192181 | 0.69102311482575574 | 0.441723 | 0.528203 | 0.311145 | 0.023277 | 0.601046 | 0.694328547 | 0.500852168 | 1.8441211 | 0.00918484501 | 0 |
| p9016_heldout_chrpair_m0p5_p4_w10_seed71_noise0p05 | chrpair_m0p5_p4_w10 | 71 | 0.5120494605999798 | 0.61288674889432004 | 0.453609 | 0.548703 | 0.310027 | 0.022159 | 0.734778 | 0.591317594 | 0.426545471 | 2.99650431 | 0.00382701875 | 0 |
| p9016_heldout_best035_seed17_noise0p05 | best035 | 17 | 0.56682710780915713 | 0.65788766041965463 | 0.450322 | 0.544229 | 0.308533 | 0.020665 | 0.689667 | 0.660164595 | 0.476208091 | 2.34140182 | 0.00420972063 | 0 |
| p9016_heldout_best035_seed71_noise0p05 | best035 | 71 | 0.60597141395149701 | 0.67156176098419584 | 0.436175 | 0.524962 | 0.302116 | 0.014248 | 0.633453 | 0.669095397 | 0.48265031 | 2.11841083 | 0.00880214313 | 0 |
| p9016_heldout_chrpair_m1_p1_w0_seed47_noise0p05 | chrpair_m1_p1_w0 | 47 | 0.51874210773125906 | 0.62794584680272991 | 0.456913 | 0.562228 | 0.297899 | 0.010031 | 0.810883 | 0.629774392 | 0.454286188 | 3.21292996 | 0.000382701875 | 0 |
| p9016_heldout_td1_seed17_noise0 | td1 | 17 | 0.59585976924760442 | 0.59760553816909245 | 0.441775 | 0.543707 | 0.287868 | 0 | 0.770306 | 0.588657796 | 0.424626857 | 3.4680357 | 0.0045924225 | 0 |
| p9016_heldout_chrpair_m1_p1_w0_seed71_noise0p05 | chrpair_m1_p1_w0 | 71 | 0.49918795283121131 | 0.6150928650832832 | 0.455582 | 0.571425 | 0.280671 | -0.007197 | 0.791879 | 0.598560393 | 0.431770056 | 3.15025091 | 0.00114810563 | 0 |

## Heldout Selector Audit

| metric | selection_direction | spearman_with_eval_trans | selected_config | selected_trans | selected_delta_vs_baseline | selected_is_eval_best | target_met |
|---|---|---|---|---|---|---|---|
| heldout_mean_expected_energy | low | -0.0577777778 | p9016_heldout_chrpair_m1_p1_w0_seed23_noise0p05 | 0.34251 | 0.054642 | 1 | 0 |
| heldout_mean_expected_energy | high | -0.0577777778 | p9016_heldout_best035_seed97_noise0p05 | 0.335188 | 0.04732 | 0 | 0 |
| heldout_mean_min_energy | low | -0.109059829 | p9016_heldout_chrpair_m1_p1_w0_seed23_noise0p05 | 0.34251 | 0.054642 | 1 | 0 |
| heldout_mean_min_energy | high | -0.109059829 | p9016_heldout_td1_seed17_noise0 | 0.287868 | 0 | 0 | 0 |
| heldout_mean_entropy | low | 0.0974358974 | p9016_heldout_chrpair_m1_p1_w0_seed23_noise0p05 | 0.34251 | 0.054642 | 1 | 0 |
| heldout_mean_entropy | high | 0.0974358974 | p9016_heldout_best035_seed83_noise0p05 | 0.311145 | 0.023277 | 0 | 0 |
| heldout_mean_pU | low | 0.0974358974 | p9016_heldout_chrpair_m1_p1_w0_seed23_noise0p05 | 0.34251 | 0.054642 | 1 | 0 |
| heldout_mean_pU | high | 0.0974358974 | p9016_heldout_best035_seed83_noise0p05 | 0.311145 | 0.023277 | 0 | 0 |
| heldout_mean_best_normalized_distance | low | -0.131623932 | p9016_heldout_chrpair_m1_p1_w0_seed23_noise0p05 | 0.34251 | 0.054642 | 1 | 0 |
| heldout_mean_best_normalized_distance | high | -0.131623932 | p9016_heldout_td1_seed17_noise0 | 0.287868 | 0 | 0 | 0 |
| heldout_short_distance_frac | low | 0.298461538 | p9016_heldout_td1_seed17_noise0 | 0.287868 | 0 | 0 | 0 |
| heldout_short_distance_frac | high | 0.298461538 | p9016_heldout_chrpair_m1_p1_w0_seed23_noise0p05 | 0.34251 | 0.054642 | 1 | 0 |

## Interpretation Boundary

The headline `HELDOUT_SELECTOR_EVAL_BEST` means that at least one blind-visible heldout/raw-fit selector picked the eval-best config within this candidate set. It does not mean the predeclared +0.1 trans improvement target was met; check `target_met` and `target_plus_0p1_vs_heldout_td1_met` for that stricter criterion.

If heldout raw-contact fit cannot select the eval-best or +0.1 trans basin, then raw-contact heldout likelihood is not sufficient as a blind selector under the current P9016 1Mb softall family. That would support the earlier diagnosis that the remaining trans identity problem is not solved by blind seed/basin selection alone.
