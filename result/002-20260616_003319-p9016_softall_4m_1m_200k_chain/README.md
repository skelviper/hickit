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

Contact accuracy 先用每条染色体 cis distance-matrix Spearman correlation 选定 reconstruction 到 CHARM/3DG 的 copy gauge，然后在这个固定 gauge 下比较 posterior top1/pmax 和 SNP phase truth；SNP phase 不参与 copy gauge 选择。

`model_top1_*` 是完整四状态 00/01/10/11 accuracy；`model_samecross_*` 只评价 same-copy(00/11) vs cross-copy(01/10)。cis 中 `model_samecross_cis` 很高，说明低于 ~0.49 的四状态 cis top1 主要来自 00 vs 11 内部标签选择，而不是 same/cross territory 判断失败。

## Stage Summary

| stage | bin_size_bp | reference_3dg | shared_points | mean_cis_spearman | truth_same_cis | truth_majority_cis | truth_samecross_rand4_cis | model_top1_all | model_top1_cis | model_top1_trans | model_samecross_all | model_samecross_cis | model_samecross_trans | model_pmax90_acc_all | model_pmax90_acc_cis | model_pmax90_acc_trans | model_pmax90_recall_all | model_pmax90_recall_cis | model_pmax90_recall_trans | charm_top1_all | charm_top1_cis | charm_top1_trans |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `4m` | 4000000 | `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz` | 1295 | 0.683784 | 0.96983 | 0.531728 | 0.484915 | 0.352892 | 0.455591 | 0.23812 | 0.723493 | 0.913945 | 0.510654 | 0.362528 | 0.441315 | 0.237298 | 0.194183 | 0.274932 | 0.103941 | 0.214292 | 0.255442 | 0.168305 |
| `1m` | 1000000 | `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz` | 4947 | 0.71501 | 0.974824 | 0.534348 | 0.487412 | 0.374549 | 0.465692 | 0.244108 | 0.759175 | 0.935654 | 0.506602 | 0.345532 | 0.439129 | 0.245342 | 0.13163 | 0.146921 | 0.109747 | 0.198369 | 0.221055 | 0.165901 |
| `200k` | 200000 | `/shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.200k.3dg.gz` | 24757 | 0.694382 | 0.978969 | 0.533319 | 0.489484 | 0.384278 | 0.468564 | 0.229915 | 0.802625 | 0.957772 | 0.518486 | 0.319269 | 0.438143 | 0.216682 | 0.117764 | 0.11574 | 0.12147 | 0.191512 | 0.206266 | 0.164491 |

## Output Layout

- `4m`: training output in `outputs/4m/minimal_soft_sep_off/`; eval output in `eval/4m/`.
- `1m`: training output in `outputs/1m/minimal_soft_sep_off/`; eval output in `eval/1m/`.
- `200k`: training output in `outputs/200k/minimal_soft_sep_off/`; eval output in `eval/200k/`.

Each eval folder contains `README.md`, `summary.tsv`, quantitative tables, and retina-DPI plots.
