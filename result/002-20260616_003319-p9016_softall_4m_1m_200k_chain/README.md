# 002-20260616_003319-p9016_softall_4m_1m_200k_chain

- baseline: `softall` / `raw_expected_soft_all`
- resolution chain: `4000000,1000000,200000`
- training input: `/mnt/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- output root: `test_res/002-20260616_003319-p9016_softall_4m_1m_200k_chain/outputs`
- relax backend: `gpu`
- n_iter: `100`
- relax_steps: `100`
- relax_step: `0.012`
- child_offset_step: `0.05`
- parent_anchor_k: `0.05`

训练只读取 P9016 raw pairs。phase labels 和 CHARM/3DG reference 只在每个 stage 的 post-training eval 中读取。4m stage 初始化一次；后续 1m 和 200k stage 从上一层坐标 coarse-to-fine lift 后继续优化，不重新初始化。

Contact identity 使用 eval-only whole-chromosome SNP cis-top1 gauge；distance metrics 使用 per-chromosome cis distance-matrix Spearman geometry gauge。

## Stage Summary

| stage | bin_size_bp | allele_sep_recon | allele_sep_charm | top1_all | top1_cis | top1_trans | p90_acc_all | p90_acc_cis | p90_acc_trans | p90_recall_all | p90_recall_cis | p90_recall_trans | dist_corr c0-r0 | c0-r1 | c1-r0 | c1-r1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `4m` | 4000000 | 3.11714 | 4.6402 | 0.384297 | 0.492105 | 0.263817 | 0.404194 | 0.487839 | 0.271243 | 0.216501 | 0.303916 | 0.11881 | 0.280978 | 0.331958 | 0.324848 | 0.311408 |
| `1m` | 1000000 | 4.93077 | 4.65535 | 0.420891 | 0.520802 | 0.277899 | 0.415719 | 0.526726 | 0.296891 | 0.158367 | 0.176228 | 0.132806 | 0.405983 | 0.43083 | 0.452764 | 0.434127 |
| `200k` | 200000 | 9.49272 | 9.55054 | 0.423974 | 0.501143 | 0.282646 | 0.398305 | 0.513653 | 0.29876 | 0.146916 | 0.135687 | 0.167482 | 0.248615 | 0.25783 | 0.266481 | 0.256162 |

## Output Layout

- `4m`: training output in `outputs/4m/minimal_soft_sep_off/`; eval output in `eval/4m/`.
- `1m`: training output in `outputs/1m/minimal_soft_sep_off/`; eval output in `eval/1m/`.
- `200k`: training output in `outputs/200k/minimal_soft_sep_off/`; eval output in `eval/200k/`.

Each eval folder contains the compact README, `summary.tsv`, `contact_accuracy.tsv`, `copy_separation.tsv`, `cis_distance_correlation_matrix.tsv`, and retina-DPI plots.
