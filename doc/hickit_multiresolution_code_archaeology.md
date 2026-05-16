# Hickit Multiresolution Code Archaeology

This report records the actual code paths inspected before adding a blind
multiresolution prototype. It distinguishes standard Hickit behavior from the
current blind/unphased diploid structural decomposition path.

Scientific guardrails:

- The blind path is a Hi-C driven diploid structural decomposition with latent
  copy-pair posterior, not genetic haplotype phasing.
- `copy0`/`copy1` are gauge-symmetric internal labels and must not be read as
  maternal/paternal.
- P9016 pair phase labels and CHARM/3DG reference coordinates are forbidden for
  training, initialization, prior construction, seed/config selection, stopping,
  and model selection.
- CHARM 3DG may be used only for post-hoc evaluation.

## Files Inspected

- Standard Hickit:
  - [main.c](/mnt/ssd/zliu/phase3/hickit/main.c)
  - [fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c)
  - [bin.c](/mnt/ssd/zliu/phase3/hickit/bin.c)
  - [hickit.h](/mnt/ssd/zliu/phase3/hickit/hickit.h)
  - [hkpriv.h](/mnt/ssd/zliu/phase3/hickit/hkpriv.h)
- Blind path:
  - [blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c)
  - [run_blind_p9016_full_cpu.c](/mnt/ssd/zliu/phase3/hickit/run_blind_p9016_full_cpu.c)
  - [run_blind_p9016_full_cpu_matrix.c](/mnt/ssd/zliu/phase3/hickit/run_blind_p9016_full_cpu_matrix.c)
  - [candidate_rerun_diagnostics.py](/mnt/ssd/zliu/phase3/hickit/candidate_rerun_diagnostics.py)
  - `test_blind_*.c` and `test_blind_eval.py`

## 1. Resolution And Bmap

### Standard Hickit Behavior

The command-line `-b NUM` path in [main.c](/mnt/ssd/zliu/phase3/hickit/main.c:240)
parses a bin size, counts neighbor support if needed, then calls:

```c
b = hk_bmap_gen(m->d, m->n_pairs, m->pairs, bin_size, bmap_skip_merge_flag);
hk_fdg(&fdg_opt, b, d3, &rng);
```

The last argument of `hk_bmap_gen()` is the bead merge skip flag. In
[bin.c](/mnt/ssd/zliu/phase3/hickit/bin.c:147), `bmap_skip_merge_flag == 0`
runs `hk_bmap_merge_beads()` twice; `1` skips merge. The standard CLI default
is `0`, while `-M` sets it to `1` in [main.c](/mnt/ssd/zliu/phase3/hickit/main.c:279).

Uniform beads are generated in
[hk_bmap_gen_beads_uniform](/mnt/ssd/zliu/phase3/hickit/bin.c:53). A bead is
nominally `size` bp, but chromosome ends are folded into the final bead when the
remaining segment is shorter than `1.5 * size`:

```c
l = size * 1.5 < len - st ? size : len - st;
```

So standard Hickit uses fixed-size beads plus terminal variable-size beads, and
optionally merges beads with no non-self, non-adjacent contact support. Bead
sizes are stored directly as `hk_bead.st`/`hk_bead.en`
([hickit.h](/mnt/ssd/zliu/phase3/hickit/hickit.h:97) and bmap fields around
[hickit.h](/mnt/ssd/zliu/phase3/hickit/hickit.h:141)).

Bead size affects standard backbone scaling. In
[fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:267), `hk_fdg_bead_size()` computes
the median bead size. Standard FDG then sets a per-backbone-edge scale
`d_opt = pow(d / mid_dist, 1/3)` where `d` is the mean bp size of adjacent beads
([fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:503)).

### Current Blind Behavior

The full P9016 blind runners currently hard-code 1Mb and skip merge:

- [run_blind_p9016_full_cpu.c](/mnt/ssd/zliu/phase3/hickit/run_blind_p9016_full_cpu.c:613):
  `hk_bmap_gen(..., HK_BLIND_P9016_FULL_RESOLUTION, 1)`
- [run_blind_p9016_full_cpu_matrix.c](/mnt/ssd/zliu/phase3/hickit/run_blind_p9016_full_cpu_matrix.c:1586):
  `hk_bmap_gen(..., HK_BLIND_P9016_MATRIX_RESOLUTION, 1)`

`HK_BLIND_P9016_MATRIX_RESOLUTION` is currently `1000000`
([run_blind_p9016_full_cpu_matrix.c](/mnt/ssd/zliu/phase3/hickit/run_blind_p9016_full_cpu_matrix.c:15)).

The existing blind manifest writes a `resolution` field but not a general
`bin_size_bp`, `resolution_label`, or merge-mode field
([run_blind_p9016_full_cpu_matrix.c](/mnt/ssd/zliu/phase3/hickit/run_blind_p9016_full_cpu_matrix.c:1048)).

## 2. Contact Scaling

### Standard Hickit Behavior

Standard Hickit converts binned contact count to contact target scale as:

```c
n_eff = fdg_contact_count(m, p, contact_scale);
d_scale = pow(n_eff, -1.0 / 3.0);
```

This is in both CPU/GPU paths; the GPU buffer creation is shown at
[fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:515). `hk_fdg_contact_energy_dist()`
then normalizes distance as:

```c
r = (distance / unit) / d_scale;
```

([fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:225)).

The count is not always just raw `n`. Standard Hickit can use GC-normalized
counts when `m->gc_norm_pairs` is present, and can apply a global
`contact_target` scaling through `fdg_contact_scale()`. The FDG contact
stiffness also uses neighbor support:

```c
k = p->max_nei >= max_nei ? 1.0f : powf((double)p->max_nei / max_nei, 1.0/3.0);
```

where `max_nei` is the median of binned pair `max_nei`
([fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:617) and
[fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:519)).

Different resolutions are handled by rebuilding a bmap at the requested `-b`
resolution. Standard count scaling remains count-driven, while backbone scaling
uses median bead size and adjacent bead sizes.

### Current Blind Behavior

The blind bpair builder in
[blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c:3856) constructs canonical
haploid bead-pairs from raw pairs and sets:

```c
base_d_scale = powf(n_raw, -1.0f / 3.0f);
base_k = 1.0f;
```

([blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c:49) and
[blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c:3905)).

The blind path records cis/trans bpair and raw counts in `hk_blind_bpair_set`
([hickit.h](/mnt/ssd/zliu/phase3/hickit/hickit.h:156)) and classifies bpairs by
chromosome during construction
([blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c:3888)).

Blind weighted-edge expansion creates four copy-state edges with:

```c
k = base_k * contact_k_multiplier * rho_eff * p4[state];
d_scale = base_d_scale;
```

or, for `expected_count`, state-specific:

```c
n_eff = n_raw * p4[state];
d_scale = max(n_eff, eps_count)^(-1/3);
```

([blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c:1235)).

Important divergence: the blind path currently does not reuse standard Hickit
neighbor-aware `max_nei` stiffness. The code explicitly says this is deferred
until blind bpairs carry `max_nei`/median context
([blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c:3906)).

## 3. FDG Initialization

### Standard Hickit Behavior

Standard `hk_fdg()` initializes coordinates randomly if no source bmap is
provided:

```c
m->x = hk_fdg_init(rng, m->n_beads, opt->target_radius);
unit = target_radius / n_beads^(1/3);
```

([fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:635)).

The seed comes from the CLI `-s` option; the main RNG is seeded through
`kr_srand_r(&rng, seed)` ([main.c](/mnt/ssd/zliu/phase3/hickit/main.c:138) and
[main.c](/mnt/ssd/zliu/phase3/hickit/main.c:276)).

If a source bmap exists, standard `hk_fdg()` uses `hk_fdg_copy_x()` to initialize
the new bmap from prior coordinates:

```c
src_dist = hk_fdg_copy_x(m, src, rng);
unit = src_dist / pow(m->n_beads / src->n_beads, 1/3);
```

([fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:631)).
`hk_fdg_copy_x()` maps destination beads to source beads by chromosome/start
coordinate, interpolates to the next same-chromosome source bead when possible,
and adds random jitter scaled by the inherited unit
([fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:278)).

Standard Hickit therefore already has a multistage coarse-to-fine pattern: each
successive `-b` operation uses the previous `d3` coordinates as the source for
the next resolution. The usage string in `main.c` includes a repeated
`-b4m -b1m -b200k ...` workflow
([main.c](/mnt/ssd/zliu/phase3/hickit/main.c:240)).

### Current Blind Behavior

Blind initialization supports:

- `toy_split`
- `unphased_scaffold_split`
- `random_diploid`
- `random_haploid_split`

The names and enum are in [hickit.h](/mnt/ssd/zliu/phase3/hickit/hickit.h:56)
and [blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c:55).

`unphased_scaffold_split` reuses standard Hickit FDG through
`hk_blind_init_haploid_scaffold_from_bmap_fdg()`, which calls `hk_fdg()` on the
unphased bmap and copies `bmap->x` into a haploid scaffold
([blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c:689)). The matrix runner uses
this path with CPU backend and `scaffold_fdg_n_iter`
([run_blind_p9016_full_cpu_matrix.c](/mnt/ssd/zliu/phase3/hickit/run_blind_p9016_full_cpu_matrix.c:597)).

The haploid scaffold is then duplicated into a gauge-symmetric diploid split by
`hk_blind_init_diploid_coords_from_haploid()`
([blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c:588)). The split is geometric
and seed-driven; it does not use phase labels.

## 4. Unit And Scale

### Standard Hickit Behavior

In standard FDG, `unit` is a coordinate normalization scale, not directly the
bp bin size:

- Random initialization: `unit = target_radius / n_beads^(1/3)`
  ([fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:637)).
- Source/lift initialization: `unit` is inherited from source average backbone
  distance and adjusted by bead count ratio
  ([fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:631)).

Contact energy normalizes distance as:

```c
r = (distance / unit) / d_scale;
```

([fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:225)).

Backbone `d_scale` accounts for actual bead size relative to median bead size
([fdg.c](/mnt/ssd/zliu/phase3/hickit/fdg.c:503)).

### Current Blind Behavior

The P9016 blind matrix runner currently sets:

```c
HK_BLIND_P9016_MATRIX_UNIT = 1.0f
HK_BLIND_P9016_MATRIX_D_SCALE = 1.0f
```

([run_blind_p9016_full_cpu_matrix.c](/mnt/ssd/zliu/phase3/hickit/run_blind_p9016_full_cpu_matrix.c:16)).

The blind posterior update uses each bpair's `base_d_scale` and `base_k`
([blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c:544)). Weighted contacts then
use each wedge's `d_scale`/`k` in the M-step force functions
([blind.c](/mnt/ssd/zliu/phase3/hickit/blind.c:1299)).

For a first blind multiresolution prototype, keeping `unit=1.0` across 4Mb and
1Mb preserves comparability with existing blind results. This is a deliberate
blind-specific divergence from standard Hickit and must be recorded in manifests.
If a standard FDG scaffold is used, the scaffold coordinates originate in
standard Hickit units; the prototype must record that source and any subsequent
normalization.

## 5. Blind-Specific Divergence Confirmed

The following current blind assumptions are confirmed in actual code:

- Fixed 1Mb bmap in the existing P9016 blind runners:
  `hk_bmap_gen(..., 1000000, 1)` with merge skipped.
- Blind raw pair conversion strips phase labels. `hk_blind_pair` stores only
  `chr`, `pos`, and `strand` ([hickit.h](/mnt/ssd/zliu/phase3/hickit/hickit.h:126)).
  `hk_blind_pair_from_pair()` copies only those fields
  ([hickit.h](/mnt/ssd/zliu/phase3/hickit/hickit.h:427)).
- Blind bpair `base_d_scale = n_raw^(-1/3)`.
- Blind bpair `base_k = 1.0`, with standard neighbor-aware stiffness not yet
  carried over.
- Blind M-step expands four copy-state weighted edges using posterior `p4` and
  training weight `rho_train`.
- `expected_count` d_scale currently affects weighted-edge expansion; posterior
  scoring still uses each bpair's `base_d_scale`.
- The standard phased split path `-S` in [main.c](/mnt/ssd/zliu/phase3/hickit/main.c:232)
  calls `hk_pair_split_phase()` and is not valid for blind training.

## Standard Hickit Behavior Vs Current Blind Behavior

| Topic | Standard Hickit | Current Blind P9016 |
|---|---|---|
| Resolution | CLI `-b NUM`, arbitrary bin size | Hard-coded 1Mb in current runners |
| Merge mode | Default merges contactless beads twice; `-M` skips | Always passes skip-merge `1` |
| Bead sizes | Fixed-size plus terminal variable-size, optional merged beads | Same bmap machinery, but fixed 1Mb skip-merge in runners |
| Contact d_scale | `n_eff^(-1/3)` where `n_eff` may include GC/contact target scaling | `n_raw^(-1/3)` per bpair; optional state `expected_count` in M-step |
| Contact stiffness | Neighbor-aware `max_nei` stiffness | Uniform `base_k=1.0`; cis/trans multipliers optional |
| Coarse-to-fine | Repeated `-b`, source coords copied/interpolated by `hk_fdg_copy_x()` | No current 4Mb-to-1Mb blind lift |
| Initialization | Random global coordinates or copied source bmap | Random diploid, random haploid split, toy split, or standard FDG haploid scaffold then blind split |
| Unit | Derived from coordinate scale and bead count/source | Fixed `unit=1.0` in current P9016 blind runners |
| Phase labels | Standard diploid `-S` can use phase | Forbidden; blind pair conversion strips phase |
| Reference/CHARM | Not part of standard training path | Eval-only; must not enter training/selection |

## Implementation Implications

Concrete problem: current blind reconstruction cannot test whether a 4Mb coarse
layout protects 1Mb refinement from self-collapse because the runner assumes a
single 1Mb bmap and has no coarse-to-fine lift or anchor.

Why it matters: recent diagnostics show severe scale compression and failed
absolute trans/centroid criteria. A multiresolution prototype tests whether the
global layout can be established at a coarser, less cis-dense scale before 1Mb
local refinement.

Small conceptual example: a 4Mb parent bead may contain four 1Mb child beads. If
the 4Mb parent copies are already arranged in a plausible global layout, the 1Mb
children can start near their parent instead of all chromosomes optimizing from a
pure 1Mb cis-heavy objective.

Validation plan:

- 1Mb bmap path must reproduce old 1Mb bmap and bpair counts.
- 4Mb bmap must have fewer beads, valid canonical bpair mapping, and no phase
  fields in blind inputs.
- 4Mb-to-1Mb map must be derived from bmap intervals, not from an assumption of
  exactly four children.
- Truth-free 4Mb seed selection must use heldout/contact diagnostics only.
- CHARM metrics may be computed post-hoc but must not select seeds or configs.
- Regression tests must scan training code paths for forbidden phase/reference
  usage.
