# FDG 三维建模（CPU/GPU）概览

## 1. 从 pairs 到三维结构
- **分箱与聚合**：`hk_bmap_gen`（`bin.c`）按照 `-b` 设定的宽度把基因组切成 beads，并合并相邻但没有任何接触的 bin，以减小节点数。pairs 坐标被映射到 bead id，累加为 `hk_bpair`（`n` 接触计数，`max_nei` 邻域密度指标）。
- **初始坐标**：默认在立方体 `[-target_radius, target_radius]^3` 内均匀随机采样（`hk_fdg_init`）。若提供粗分辨率结构（`-I`），`hk_fdg_copy_x` 会按基因组位置插值并加小扰动，将旧结构放大/缩放到新的 bead 网格。
- **距离缩放**：定义
$$
  r = \|x_i - x_j\| , \qquad \hat r = \frac{r}{\text{unit} \, d_\text{scale}}
$$

  其中 `unit = target_radius / n_beads^{1/3}`。`d_scale` 为每条边的期望尺度：
  - backbone：$d_\text{scale} = \big((\ell_i+\ell_{i+1})/(2 \cdot \text{mid\_dist})\big)^{1/3}$，`mid_dist` 为 bead 长度中位数。
  - 接触：$d_\text{scale} = n_\text{eff}^{-1/3}$，其中 $n_\text{eff}$ 为接触的有效计数（可被 GC 校正，见 §3）。
  接触的弹性系数 `k` 会随局部稠密度衰减：

  $$
  k =
  \begin{cases}
  1, & m \ge \mathrm{median}(\text{max\_nei}) \\
  \big(m/\mathrm{median}(\text{max\_nei})\big)^{1/3}, & \text{otherwise}
  \end{cases}
  $$

- **吸引力势函数**：backbone 与接触共用分段势能
  
  $$
  U(\hat r) =
  \begin{cases}
    k(d_1-\hat r)^2, & \hat r < d_1 \\
    0, & d_1 \le \hat r \le d_2 \\
    k(\hat r-d_2)^2, & d_2 < \hat r \le d_3 \\
    k\big[c_1(\hat r-d_3) + c_2/(\hat r-d_2)\big], & \hat r > d_3
  \end{cases},
  \quad c_1 = 3(d_3-d_2),\ c_2 = (d_3-d_2)^3
  $$
  力 $\mathbf{F} = -\partial U/\partial \hat r\,\hat{\mathbf{r}}$ 对称施加并在整合前由 `max_f` 截顶。默认参数：$(d_{b1}, d_{b2}) = (0.1, 1.1)$，$(d_{c1}, d_{c2}, d_{c3}) = (0.5, 1.5, 2.0)$ 分别用于 backbone 与接触。

- **排斥项**：对未在吸引列表中的 bead 对，若 $\hat r < d_r$ 则
  
  $$
  U_\text{rep}(\hat r) = k_\text{rep}(d_r-\hat r)^2,\qquad
  k_\text{rep} = k_\text{rel\_rep}\,\cdot\,\text{rel\_rep\_k}
  $$
  
  $$
  \text{rel\_rep\_k} = \frac{1}{1 + e^{-\alpha(t-\tau)}} ,\qquad t = \frac{\text{iter}+1}{n_\text{iter}}
  $$
  
  其中 $\alpha=10,\ \tau=1/3$；随迭代逐步打开排斥强度。所有 backbone/接触边会进入 blocklist，避免重复施加排斥。

- **积分与多级建模**：每轮更新
  
  $$
  x_{t+1} = x_t + \text{coef\_moment}\,(x_t - x_{t-1}) + \text{step}\,\mathbf{F}_t,
  \quad \text{step} = \text{opt.step}\,\cdot\,\text{unit}\,(n_\text{beads}/1500)^{1/3}
  $$
  
  并记录 RMS force 最小的结构作为最佳结果。可多次使用 `-b` 由粗到细重复建模，后一级以 `-I` 读取前一级结构作为初始值。

## 2. CPU 与 GPU 实现差异
- **CPU (`fdg.c:hk_fdg1_cpu`)**：每轮重新遍历 backbone/接触边计算吸引力；排斥使用“按 X 排序 + Y 轴 AVL” 的滑动窗口找近邻（近似 $O(n\log n)$）。blocklist 由哈希保存，确保已存在吸引力的边不参与排斥。
- **GPU (`fdg_gpu.cu:hk_fdg1_gpu`)**：首次迭代把 backbone/接触边打包上传，后续复用；`hk_fdg_gpu_compute` 调度多核 kernel：`fdg_force_kernel` 计算吸引力，`fdg_repulsion_kernel` 建立统一网格并在 3×3×3 相邻 cell 中并行查找排斥对（附带 blocklist），`fdg_update_kernel` 负责动量积分与 RMS force 统计。若 GPU 不可用或出错会自动退回 CPU；`--fdg-backend=cpu|gpu|auto` 控制选择。
- **差异要点**：GPU 缓存吸引力边、并行网格邻域搜索和原子累加统计，适合高节点场景；CPU 每轮重算吸引力、序贯邻域搜索，代码路径简单但在大规模下更耗时。

## 3. GC/CpG correction（接触权重预处理）
- **使用**：`-g/--gc-corr` 打开校正，`-G/--cpg <bedGraph>` 提供 CpG 曲线。
- **流程 (`bin.c:hk_bmap_apply_gc_correction`)**：
  1. 读取 CpG 密度并映射到 beads。
  2. 计算每个 bead 的接触覆盖度 $\text{cov}_i$，对 $x=\text{CpG}$ 拟合 $\text{cov} \approx a + b x + c x^2$。
  3. 得到 $\text{bias}[i]=\max(a + b x_i + c x_i^2, 10^{-6})$，再缩放使均值为 1。
  4. 每条接触的期望偏置 $\text{gc\_exp} = \text{bias}[i]\,\text{bias}[j]$，有效计数 $\text{gc\_norm\_n} = n / \text{gc\_exp}$；设置 `gc_corrected=1`，可用 `HK_GC_BIAS_OUT` 导出 bias。
- **与 FDG 交互**：`fdg_contact_count` 在建模时优先用 $\text{gc\_norm\_n}$ 作为 $n_\text{eff}$，否则回退到原始计数 `n`，从而在距离尺度 $d_\text{scale} = n_\text{eff}^{-1/3}$ 上抑制覆盖度偏高的接触、放大偏低的接触。
