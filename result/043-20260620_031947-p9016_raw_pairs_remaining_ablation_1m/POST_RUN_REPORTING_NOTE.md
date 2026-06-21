# Post-Run Reporting Note

- Timestamp: 2026-06-20T03:58:00+08:00
- Scope: reporting/provenance cleanup only.
- Training outputs were not modified.
- Eval outputs were not modified.

After independent review, the reporting files were regenerated from the existing
training/eval outputs so that `raw_pairs_remaining_delta_summary.tsv`,
`headline.txt`, and the README use only final 1 Mb rows with eval metrics. This
excludes intermediate 4 Mb chain-refinement rows from the delta table and
headline calculation.

The 043 launcher in `scripts_snapshot/` was also synchronized to the current
repository script so that the archived script matches the final reporting logic.
This change does not affect any model run, audit log, manifest, coordinate
output, posterior output, or eval output.
