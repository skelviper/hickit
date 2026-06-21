# Superseded Result

This 018 light result is superseded by:

`/mnt/ssd/zliu/phase3/hickit/result/018-20260618_220644-p9016_fixed_posterior_graph_semantics_1m`

Reason:

- The original run used a direct fixed binned-posterior state-edge graph but the manifest/README wording made it look like a native Hickit `raw_expected_soft_all` replay.
- The corrected run records `runner_family=p9016_fixed_posterior_diagnostic`, `mstep_graph_mode=fixed_bpair_state_edges`, `native_softall_replay=0`, and includes a `source_016_before_018_relax` row.
- The corrected run also uses `repulsion_multiplier=1`, matching the 016 source setting.

Do not use this directory for scientific interpretation except as an archived draft.
