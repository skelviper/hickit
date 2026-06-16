# 006-20260616_161250-p9016_softall_random_diploid_4m_1m_chain

- baseline: `softall` / `raw_expected_soft_all`
- resolution chain: `4000000,1000000`
- training input: `/mnt/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- output root: `/mnt/ssd/zliu/phase3/test_res/006-20260616_161250-p9016_softall_random_diploid_4m_1m_chain/outputs`
- relax backend: `gpu`
- n_iter: `100`
- relax_steps: `100`
- relax_step: `0.012`
- init_mode: `random_diploid`
- init_scale: `10`
- init_seed: `17`
- child_offset_step: `0.05`
- parent_anchor_k: `0.05`

训练只读取 P9016 raw pairs。phase labels 和 CHARM/3DG reference 只在每个 stage 的 post-training eval 中读取。第一个 stage 使用 `init_mode` 初始化一次；后续 stage 从上一层坐标 coarse-to-fine lift 后继续优化，不重新初始化。

结构图和 distance metrics 使用每条染色体 cis distance-matrix Spearman correlation 选定 reconstruction 到 CHARM/3DG 的 geometry gauge。contact identity headline metrics 使用 eval-only whole-chrom SNP cis-top1 oracle gauge；geometry-gauge contact metrics 保留为 `model_geometry_*` 诊断项。phase labels 和 CHARM/3DG reference 都只在 post-training eval 中读取。

## Stage Summary

| stage | bin_size_bp | allele_sep_recon | allele_sep_charm | top1_all | top1_cis | top1_trans | p90_acc_all | p90_acc_cis | p90_acc_trans | p90_recall_all | p90_recall_cis | p90_recall_trans | dist_corr_c0_r0 | dist_corr_c0_r1 | dist_corr_c1_r0 | dist_corr_c1_r1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `4m` | 4000000 | 2.21536 | 4.6402 | 0.270969 | 0.28812 | 0.251802 | 0.277651 | 0.295183 | 0.249743 | 0.161141 | 0.199364 | 0.118426 | 0.148377 | 0.144385 | 0.169385 | 0.145326 |
| `1m` | 1000000 | 3.55894 | 4.65535 | 0.311627 | 0.348544 | 0.258791 | 0.299162 | 0.323504 | 0.266358 | 0.126773 | 0.133681 | 0.116887 | 0.210974 | 0.206734 | 0.237897 | 0.212121 |

## Output Layout

- `4m`: training output in `outputs/4m/minimal_soft_sep_off/`; eval output in `eval/4m/`.
- `1m`: training output in `outputs/1m/minimal_soft_sep_off/`; eval output in `eval/1m/`.

Each eval folder contains `README.md`, `summary.tsv`, quantitative tables, and retina-DPI plots.
