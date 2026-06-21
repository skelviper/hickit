# 042 P9016 Trans Gate M-step 1Mb

This controlled blind-training experiment tests whether low-confidence trans bpair contacts should be removed from the softall M-step graph. The gate is computed only from the current blind bpair posterior and is applied after the E-step and before softall graph construction.

- full result root: `/mnt/ssd/zliu/phase3/test_res/042-20260620_022516-p9016_trans_gate_mstep_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/042-20260620_022516-p9016_trans_gate_mstep_1m`
- headline: `NO_PLUS_0P1`
- training input: raw P9016 pairs only
- eval-only inputs: SNP phase labels and CHARM/3DG are used only after training through the standard eval wrapper
- unchanged settings: scaffold split init, posterior_count gamma1 dscale, uniform prior, constant rho, temperature 1, no trans chr-pair prior, no top1 M-step hardening
- changed knobs: trans gate mode/threshold, trans dscale controls, local/global copytrack controls, plus one separation stress-control row

## Main Results

| config_name | trans_gate_mode | trans_gate_min_neg_entropy | final_n_softall_gate_skip_raw | n_raw_trans | frac_trans_raw_gated | trans_dscale_multiplier | lambda_copytrack | lambda_global_copytrack | model_top1_accuracy_genome_all | model_top1_accuracy_genome_cis | model_top1_accuracy_genome_trans | delta_trans_vs_td1_baseline | delta_trans_vs_best035 | model_same_cross_accuracy_genome_trans | mean_per_chrom_cis_distance_spearman | final_mean_entropy | final_mean_pU | final_min_sep | sep_p05 | copytrack_frac_cos_lt_0 | copytrack_frac_projection_sign_switch |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p9016_transgate_off_best035_td0p5_lct0p03_gct0p003 | off | -1.38629436 | 0 | 568434 | 0 | 0.5 | 0.0299999993 | 0.00300000003 | 0.481885 | 0.572114 | 0.352751 | 0.03536 | 0 | 0.552484 | 0.769473 | 0.696720541 | 0.502577662 | 0.459196299 | 2.83363295 | 0.0015308075 | 0.0045924225 |
| p9016_transgate_off_td0p5 | off | -1.38629436 | 0 | 568434 | 0 | 0.5 | 0 | 0 | 0.474777 | 0.565892 | 0.344374 | 0.026983 | -0.008377 | 0.538392 | 0.75945 | 0.699788868 | 0.504790962 | 0.324222177 | 2.73347139 | 0.00191350938 | 0.00574052813 |
| p9016_transgate_ne0p20_td0p5_no_copytrack | neg_entropy | -0.157075733 | 388247 | 568434 | 0.683011572 | 0.5 | 0 | 0 | 0.48339 | 0.588856 | 0.332449 | 0.015058 | -0.020302 | 0.532224 | 0.736288 | 0.539320409 | 0.38903746 | 0.519316018 | 3.34258127 | 0.00420972063 | 0.0103329506 |
| p9016_transgate_pmax0p5_best035 | pmax | -1.38629436 | 54873 | 568434 | 0.0965336345 | 0.5 | 0.0299999993 | 0.00300000003 | 0.477371 | 0.578811 | 0.332192 | 0.014801 | -0.020559 | 0.530335 | 0.784239 | 0.670489788 | 0.483656168 | 1.00741613 | 2.83499002 | 0.00076540375 | 0.00267891313 |
| p9016_transgate_ne0p20_best035 | neg_entropy | -0.157075733 | 385559 | 568434 | 0.678282791 | 0.5 | 0.0299999993 | 0.00300000003 | 0.485184 | 0.592258 | 0.331942 | 0.014551 | -0.020809 | 0.540184 | 0.74861 | 0.534083724 | 0.385259986 | 0.589244723 | 3.27710485 | 0.00229621125 | 0.00612323 |
| p9016_transgate_ne0p20_best035_msep1_lsep0p5 | neg_entropy | -0.157075733 | 388919 | 568434 | 0.684193767 | 0.5 | 0.0299999993 | 0.00300000003 | 0.48849 | 0.6008 | 0.327754 | 0.010363 | -0.024997 | 0.538774 | 0.761912 | 0.552224159 | 0.39834553 | 0.876912475 | 2.70791554 | 0.00420972063 | 0.00497512438 |
| p9016_transgate_ne0p20_td1_no_copytrack | neg_entropy | -0.157075733 | 368992 | 568434 | 0.649137807 | 1 | 0 | 0 | 0.476039 | 0.580111 | 0.327094 | 0.009703 | -0.025657 | 0.530758 | 0.770559 | 0.513487339 | 0.370402843 | 0.740771234 | 3.59059691 | 0.00344431688 | 0.00880214313 |
| p9016_transgate_ne0p30_best035 | neg_entropy | -0.326260984 | 338171 | 568434 | 0.594916912 | 0.5 | 0.0299999993 | 0.00300000003 | 0.481342 | 0.589157 | 0.327038 | 0.009647 | -0.025713 | 0.528862 | 0.790367 | 0.575293064 | 0.414986223 | 0.68467474 | 3.2904644 | 0.00267891313 | 0.00574052813 |
| p9016_transgate_pmax0p4_best035 | pmax | -1.38629436 | 9787 | 568434 | 0.0172174782 | 0.5 | 0.0299999993 | 0.00300000003 | 0.469986 | 0.570929 | 0.325517 | 0.008126 | -0.027234 | 0.532147 | 0.798421 | 0.674645007 | 0.486653507 | 0.725424588 | 2.96730685 | 0.00114810563 | 0.00344431688 |
| p9016_transgate_ne0p30_td0p5_no_copytrack | neg_entropy | -0.326260984 | 331654 | 568434 | 0.583452081 | 0.5 | 0 | 0 | 0.482445 | 0.592443 | 0.325017 | 0.007626 | -0.027734 | 0.522146 | 0.809395 | 0.572727382 | 0.413135469 | 0.622862041 | 3.68679976 | 0.00420972063 | 0.0045924225 |
| p9016_transgate_ne0p30_td1_no_copytrack | neg_entropy | -0.326260984 | 308117 | 568434 | 0.542045339 | 1 | 0 | 0 | 0.46796 | 0.571589 | 0.319648 | 0.002257 | -0.033103 | 0.528425 | 0.773537 | 0.529510617 | 0.381961167 | 0.917859614 | 3.64256287 | 0.00344431688 | 0.00382701875 |
| p9016_transgate_off_td1_baseline | off | -1.38629436 | 0 | 568434 | 0 | 1 | 0 | 0 | 0.464661 | 0.567561 | 0.317391 | 0 | -0.03536 | 0.514353 | 0.802138 | 0.611671507 | 0.441227734 | 1.53334069 | 3.65623307 | 0.00191350938 | 0.003061615 |
| p9016_transgate_ne0p10_best035 | neg_entropy | -0.0116029605 | 520035 | 568434 | 0.914855551 | 0.5 | 0.0299999993 | 0.00300000003 | 0.431435 | 0.547261 | 0.265667 | -0.051724 | -0.087084 | 0.497614 | 0.68524 | 0.667218566 | 0.48129645 | 0.254818916 | 1.55842245 | 0.00688863375 | 0.0107156525 |
| p9016_transgate_ne0p05_best035 | neg_entropy | -0.000426874933 | 564760 | 568434 | 0.993536629 | 0.5 | 0.0299999993 | 0.00300000003 | 0.416288 | 0.533784 | 0.24813 | -0.069261 | -0.104621 | 0.496774 | 0.632907 | 0.836944163 | 0.603727639 | 0.506086946 | 1.29345071 | 0.00612323 | 0.0110983544 |

## Interpretation Boundary

A success requires full-denominator trans top1 accuracy to increase by at least 0.1 over the ungated posterior-count td1 baseline. Callable-subset gains are not sufficient for this headline. If the gate skips many trans raw contacts without full-denominator trans improvement, the bottleneck is not simply low-confidence trans force pollution.
