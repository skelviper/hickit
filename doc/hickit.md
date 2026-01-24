# Force-Directed Model and GC Normalization

This note summarizes how the current code implements the force model and the optional GC/CpG-based normalization (`-g` / `--gc-corr`, `-G` / `--cpg`).

## Distance normalization
- Physical distance: `dist = ||x_i - x_j||`.
- Model distance: `dist_norm = dist / unit`, where `unit` is the current bead scaling.
- Each pair carries a distance scale `d_scale`:
  - Backbone: `d_scale = ((bead_len_avg) / mid_dist)^{1/3}` (per adjacent bead pair).
  - Contact: `d_scale = n_eff^{-1/3}`, where `n_eff` is the effective contact count (GC-normalized if available; see below).
- The force law uses `r = dist_norm / d_scale`.

## Piecewise potential and force
The potential for attractive/backbone and contact pairs follows the same shape parameters `(d1,d2,d3,k)`:

```
U(r) =
  k(d1 - r)^2                     , r < d1
  0                               , d1 ≤ r ≤ d2
  k(r - d2)^2                     , d2 < r ≤ d3
  k[c1 (r - d3) + c2 / (r - d2)]  , r > d3
```

with `c1 = 3(d3 - d2)` and `c2 = (d3 - d2)^3`.

In code:
- Backbone uses `(d_b1, d_b2, - , k=1)` and only the inner quadratic + outer quadratic parts.
- Contacts use `(d_c1, d_c2, d_c3, k)` with the full four-piece form.
- Default contact thresholds: `d_c1=0.5`, `d_c2=1.5`, `d_c3=2.0` (`fdg.c:161`).
- Force is the **negative gradient** of `U`:
  - `r < d1`: `F = +2k(d1 - r) * r_hat`
  - `d1 ≤ r ≤ d2`: `F = 0`
  - `d2 < r ≤ d3`: `F = -2k(r - d2) * r_hat`
  - `r > d3`: `F = -k[c1 - c2/(r - d2)^2] * r_hat`
- CPU implementation: `update_force` (`fdg.c:248-287`); GPU mirror: `fdg_force_kernel` (`fdg_gpu.cu:352-376`).
- Forces are applied symmetrically (`+F` to `i`, `-F` to `j`), with optional clipping to `max_f` during integration.

## Force coefficients
- Contact stiffness `k` is scaled by local neighborhood: `k = 1` when `max_nei >= max_nei_threshold`; otherwise `k = (max_nei / threshold)^{1/3}` (`fdg.c:323-327`).
- Repulsion uses a separate piecewise quadratic cutoff at `d_r` (not detailed here; see `fdg.c` and `fdg_gpu.cu`).

## Repulsion (collision avoidance)
- Potential: `U_rep(r) = k_rep (d_r - r)^2` for `r < d_r`, else 0; force `F = +2 k_rep (d_r - r) * r_hat`.
- Parameters: default `d_r = 2.0`, base stiffness `k_rel_rep = 0.05` (`fdg.c:157`), multiplied by an annealed factor `rel_rep_k` per iteration: `rel_rep_k = 1 / (1 + exp(-alpha * (t - turning)))` with `alpha=10`, `turning=1/3`, `t = (iter+1)/n_iter` (`fdg.c:600-618`). Effective `k_rep = k_rel_rep * rel_rep_k`.
- Cutoff radius for neighbor search: `rep_radius = d_r * unit`. Only pairs closer than this are considered.
- Pair exclusion: repulsion is skipped for backbone/contact pairs (CPU checks a hash of attractive edges; GPU uses a blocklist).
- CPU neighbor search: sweep-sort on X plus AVL tree on Y to find candidates within `rep_radius`; force/energy computed via `update_force(..., FORCE_REPEL, d_scale=1)` (`fdg.c:332-377`).
- GPU neighbor search: build a uniform grid with cell size `unit * rep_radius`; kernel scans 3x3x3 neighboring cells, skips blocked pairs, and applies the same quadratic repulsion (`fdg_gpu.cu:430-509`). Forces accumulated symmetrically; stats track active repulsive pairs.

## GC/CpG normalization (`-g` and `-G/--cpg`)
1. Load CpG content per bead (`hk_bmap_load_cpg` via `hk_bmap_apply_gc_correction` in `bin.c`).
2. Fit a quadratic model of coverage vs. CpG fraction: `cov ≈ a + b*x + c*x^2`, using beads with CpG data (`bin.c:326-349`).
3. Predict per-bead bias `bias[i] = a + b*x + c*x^2`, floor at `1e-6`, then renormalize so mean bias over used beads is 1 (`bin.c:351-365`).
4. For each contact between beads `i,j`:
   - Expected multiplicative bias: `gc_exp = bias[i] * bias[j]` (floored at `1e-6`).
   - GC-normalized count: `gc_norm_n = raw_count / gc_exp` (`bin.c:367-373`).
5. Mark the bead map as GC-corrected (`gc_corrected=1`) and store per-bead bias.
6. During force computation (`fdg_contact_count` in `fdg.c:131-136`), if `gc_corrected` is set and `gc_norm_n > 0`, the contact weight `n_eff` uses `gc_norm_n`; otherwise it falls back to the raw count.
7. The distance scale for the contact is `d_scale = n_eff^{-1/3}`; this is the only place the normalized count enters the force law.

## CLI switches
- `-g` / `--gc-corr`: enable GC/CpG-based normalization for subsequent `-b` binning/modeling.
- `-G <file>` / `--cpg=<file>`: provide the CpG track (bedGraph-like) required by `-g`.
- If `-g` is set without a CpG file when binning, the program aborts with an error (`main.c:208-220`). 
