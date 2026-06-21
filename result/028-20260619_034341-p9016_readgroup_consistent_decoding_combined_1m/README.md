# 028 P9016 Read-Group Consistent Decoding

This diagnostic uses contacts.seg read grouping after training to decode a read-consistent copy assignment. It is not training; construction uses no SNP phase labels and no CHARM/3DG. SNP labels are used only for evaluation.

- full result root: `/mnt/ssd/zliu/phase3/test_res/028-20260619_034341-p9016_readgroup_consistent_decoding_combined_1m`
- headline: `NO_READGROUP_DECODING_PLUS_0P1`

| source config | policy | trans top1 | delta | same/cross | decoded fraction |
|---|---|---:|---:|---:|---:|
| baseline | baseline | 0.315369816 | 0 | 0.511790078 | 0 |
| baseline | hybrid_lambda0.25 | 0.315619856 | 0.00025004 | 0.511776186 | 0.983101468 |
| baseline | hybrid_lambda0.5 | 0.316314411 | 0.000944595 | 0.511859533 | 0.983101468 |
| baseline | hybrid_lambda0.75 | 0.31630052 | 0.000930704 | 0.511498364 | 0.983101468 |
| baseline | hybrid_lambda1 | 0.316633907 | 0.001264091 | 0.511685894 | 0.983101468 |
| baseline | readgroup_only | 0.317385407 | 0.002015591 | 0.511819646 | 1 |
| sep_copytrack | baseline | 0.342332456 | 0 | 0.53934309 | 0 |
| sep_copytrack | hybrid_lambda0.25 | 0.342415803 | 8.3347e-05 | 0.539176396 | 0.983101468 |
| sep_copytrack | hybrid_lambda0.5 | 0.34313814 | 0.000805684 | 0.539259743 | 0.983101468 |
| sep_copytrack | hybrid_lambda0.75 | 0.343436799 | 0.001104343 | 0.539190287 | 0.983101468 |
| sep_copytrack | hybrid_lambda1 | 0.343756294 | 0.001423838 | 0.539259743 | 0.983101468 |
| sep_copytrack | readgroup_only | 0.345298988 | 0.002966532 | 0.539895721 | 1 |

## Interpretation

- Read-group consistent decoding covers most trans eval contacts but improves full trans top1 only marginally.
- This argues against a simple missing read-group consistency decoder being sufficient for +0.1 full-denominator trans.
