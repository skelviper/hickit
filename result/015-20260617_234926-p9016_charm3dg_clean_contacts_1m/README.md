# 015-20260617_234926-p9016_charm3dg_clean_contacts_1m

Positive-control experiment: contacts were generated from P9016 CHARM 20kb 3DG geometry at radius thresholds, downsampled to the observed contact level, then used as unphased input for the existing blind softall runner.

## Paths

- full result root: `/mnt/ssd/zliu/phase3/test_res/015-20260617_234926-p9016_charm3dg_clean_contacts_1m`
- light result root: `/mnt/ssd/zliu/phase3/hickit/result/015-20260617_234926-p9016_charm3dg_clean_contacts_1m`
- summary.tsv: `/mnt/ssd/zliu/phase3/test_res/015-20260617_234926-p9016_charm3dg_clean_contacts_1m/summary.tsv`
- synthetic pairs: `/mnt/ssd/zliu/phase3/test_res/015-20260617_234926-p9016_charm3dg_clean_contacts_1m/synthetic_pairs`
- commands.log: `/mnt/ssd/zliu/phase3/test_res/015-20260617_234926-p9016_charm3dg_clean_contacts_1m/commands.log`

## Boundary

This is not a blind baseline. Training contacts are reference-derived from CHARM/3DG, but no SNP phase labels are written into the synthetic pairs. Manifests must report `input_contact_source=charm3dg_derived_pairs`, `uses_phase_labels=0`, and `uses_charm_for_training=1`.

## Git And Build

- git commit at launch: `1037ee5f79d54e943eef5c37e5aeb2d841698cb8`
- git dirty count at launch: `8`
- WARNING: dirty tree at launch; inspect `logs/git_status.txt` and `logs/git_diff.patch` in the full result root.
- run binary sha256: `563303f22b0cbbfa2da75c2af460468b4ce2f0cc3c26136c0857512c1ac5a6f0`
- backend: `gpu`

## Results

| config | radius_pr | contacts | trans contacts | top1 all | top1 cis | top1 trans | same/cross trans | cis Spearman | entropy | pU | min sep | mean sep | cos<0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `p9016_charm3dg20k_1pr_clean_contacts_pcgamma1_ieps0p5_noise0_seed17` | 1.0 | 1265114 | 108119 | 0.453938 | 0.468187 | 0.350884 | 0.540892 | 0.40785 | 0.657832682 | 0.474525958 | 0.214037076 | 2.92738295 | 0.0394787275 |
| `p9016_charm3dg20k_2pr_clean_contacts_pcgamma1_ieps0p5_noise0_seed17` | 2.0 | 1265114 | 162836 | 0.425772 | 0.435357 | 0.362535 | 0.540008 | 0.377934 | 0.76606375 | 0.552598178 | 0.145999104 | 2.69791913 | 0.06094289 |

## Interpretation Rules

- If trans improves strongly here, the raw-pairs failure is consistent with noise/sparsity or misleading contacts being a major bottleneck.
- If trans still fails with CHARM-derived clean contacts, the bottleneck is likely in gauge/copy-track identifiability or EM/FDG dynamics rather than contact noise alone.
- Do not treat this as a production model improvement because training used CHARM/3DG-derived contacts.
