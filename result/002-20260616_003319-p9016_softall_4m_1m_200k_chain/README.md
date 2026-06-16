# 002-20260616_003319-p9016_softall_4m_1m_200k_chain

- baseline: `softall` / `raw_expected_soft_all`
- resolution chain: `4000000,1000000,200000`
- training input: `/mnt/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz`
- output root: `/mnt/ssd/zliu/phase3/test_res/002-20260616_003319-p9016_softall_4m_1m_200k_chain/outputs`
- relax backend: `gpu`
- n_iter: `100`
- relax_steps: `100`
- relax_step: `0.012`
- child_offset_step: `0.05`
- parent_anchor_k: `0.05`

训练只读取 P9016 raw pairs。phase labels 和 CHARM/3DG reference 只在每个 stage 的 post-training eval 中读取。4m stage 初始化一次；后续 1m 和 200k stage 从上一层坐标 coarse-to-fine lift 后继续优化，不重新初始化。

## Stage Summary

| stage | bin_size_bp | reference_3dg | shared_points | mean_cis_spearman | model_top1_all | model_pmax90_acc_all | model_pmax90_recall_all | charm_top1_all |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `4m` | 4000000 | `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz` | 1295 | 0.683784 | 0.384297 | 0.404194 | 0.216501 | 0.214292 |
| `1m` | 1000000 | `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz` | 4947 | 0.71501 | 0.426059 | 0.425573 | 0.162121 | 0.198369 |
| `200k` | 200000 | `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.200k.3dg.gz` | 24757 | 0.694382 | 0.423974 | 0.398305 | 0.146916 | 0.191512 |

## Output Layout

- `4m`: training output in `outputs/4m/minimal_soft_sep_off/`; eval output in `eval/4m/`.
- `1m`: training output in `outputs/1m/minimal_soft_sep_off/`; eval output in `eval/1m/`.
- `200k`: training output in `outputs/200k/minimal_soft_sep_off/`; eval output in `eval/200k/`.

Each eval folder contains `README.md`, `summary.tsv`, quantitative tables, and retina-DPI plots.
