#!/usr/bin/env python3
"""Strict post-hoc audit for experiment 049.

This script validates the completed 049 matrix without touching training
artifacts. It is intentionally independent from the runner so it can be run
after a long job that used an older inline audit.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


EXPECTED_CONFIGS = [
    "p9016_transgamma_g0_td0p75_seed17_best035",
    "p9016_transgamma_g0p25_td0p75_seed17_best035",
    "p9016_transgamma_g0p5_td0p75_seed17_best035",
    "p9016_transgamma_g0p75_td0p75_seed17_best035",
    "p9016_transgamma_g1_td0p75_seed17_best035",
    "p9016_transgamma_g0_td1_seed17_best035",
    "p9016_transgamma_g0p25_td1_seed17_best035",
    "p9016_transgamma_g0p5_td1_seed17_best035",
    "p9016_transgamma_g0p75_td1_seed17_best035",
    "p9016_transgamma_g1_td1_seed17_best035",
    "p9016_transgamma_g1_td0p5_seed17_best035_replay",
    "p9016_transgamma_g0_td0p75_seed71_best035",
    "p9016_transgamma_g0p25_td0p75_seed71_best035",
    "p9016_transgamma_g0p5_td0p75_seed71_best035",
    "p9016_transgamma_g0p75_td0p75_seed71_best035",
    "p9016_transgamma_g1_td0p75_seed71_best035",
    "p9016_transgamma_g0_td1_seed71_best035",
    "p9016_transgamma_g0p25_td1_seed71_best035",
    "p9016_transgamma_g0p5_td1_seed71_best035",
    "p9016_transgamma_g0p75_td1_seed71_best035",
    "p9016_transgamma_g1_td1_seed71_best035",
    "p9016_transgamma_g1_td0p5_seed71_best035_replay",
]


REQUIRED_FINITE_SUMMARY_COLUMNS = [
    "model_top1_accuracy_genome_all",
    "model_top1_accuracy_genome_cis",
    "model_top1_accuracy_genome_trans",
    "model_same_cross_accuracy_genome_trans",
    "mean_per_chrom_cis_distance_spearman",
    "final_mean_entropy",
    "final_mean_pU",
    "final_min_sep",
    "final_mean_sep",
    "trans_d_scale_posterior_gamma",
    "trans_dscale_multiplier",
]


REQUIRED_MANIFEST_VALUES = {
    "sample": "P9016",
    "input_contact_source": "raw_pairs",
    "approved_p9016_raw_pairs_realpath": "1",
    "uses_phase_labels": "0",
    "uses_charm_or_reference": "0",
    "copy_labels_are_gauge_only": "1",
    "dscale_mode": "posterior_count",
    "d_scale_mode": "posterior_count",
    "d_scale_posterior_gamma": "1",
    "trans_k_multiplier": "1",
    "rho_train_mode": "constant",
    "rho_train": "1",
    "rho_train_floor": "0",
    "temperature_start": "1",
    "temperature_end": "1",
    "trans_top1_mstep_mode": "off",
    "trans_gate_mode": "off",
    "trans_chr_pair_prior_mode": "off",
    "trans_chr_pair_mstep_mode": "off",
    "trans_callable_anchor_mode": "off",
    "copytrack_uses_phase_labels": "0",
    "copytrack_uses_charm_or_reference": "0",
    "global_copytrack_uses_phase_labels": "0",
    "global_copytrack_uses_charm_or_reference": "0",
    "trans_contact_scaling_uses_phase_labels": "0",
    "trans_contact_scaling_uses_charm_or_reference": "0",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_key_value_tsv(path: Path) -> tuple[dict[str, str], list[str]]:
    data: dict[str, str] = {}
    duplicates: list[str] = []
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            key = row[0]
            value = row[1] if len(row) > 1 else ""
            if key in data:
                duplicates.append(key)
            data[key] = value
    return data, duplicates


def finite_float(value: str) -> bool:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(parsed)


def close_float_string(actual: str, expected: str, tol: float = 1e-8) -> bool:
    try:
        return abs(float(actual) - float(expected)) <= tol
    except (TypeError, ValueError):
        return actual == expected


def expected_from_config(config: str) -> tuple[str, str, str, str]:
    seed = "71" if "_seed71_" in config else "17"
    noise = "0.05" if seed == "71" else "0"
    if "_td0p5_" in config:
        td = "0.5"
    elif "_td0p75_" in config:
        td = "0.75"
    else:
        td = "1"
    gamma_token = config.split("_transgamma_", 1)[1].split("_td", 1)[0]
    gamma = {
        "g0": "0",
        "g0p25": "0.25",
        "g0p5": "0.5",
        "g0p75": "0.75",
        "g1": "1",
    }[gamma_token]
    return seed, noise, gamma, td


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: scripts/audit_p9016_049_strict.py <full_run_root>", file=sys.stderr)
        return 2
    root = Path(argv[1]).resolve()
    summary_path = root / "summary.tsv"
    if not summary_path.exists():
        raise SystemExit(f"missing summary.tsv: {summary_path}")

    rows = read_tsv(summary_path)
    by_config = {row.get("config_name", ""): row for row in rows}
    audit_rows: list[dict[str, str]] = []

    def add(check: str, config: str, status: str, detail: str) -> None:
        audit_rows.append({"check": check, "config_name": config, "status": status, "detail": detail})

    missing = sorted(set(EXPECTED_CONFIGS) - set(by_config))
    extra = sorted(set(by_config) - set(EXPECTED_CONFIGS))
    add("expected_config_set", "ALL", "FAIL" if missing or extra else "PASS",
        f"missing={','.join(missing) or 'NA'} extra={','.join(extra) or 'NA'}")
    add("expected_config_count", "ALL", "PASS" if len(rows) == len(EXPECTED_CONFIGS) else "FAIL",
        f"observed={len(rows)} expected={len(EXPECTED_CONFIGS)}")

    for config in EXPECTED_CONFIGS:
        row = by_config.get(config)
        if row is None:
            continue
        if row.get("status") != "OK":
            add("summary_status_ok", config, "FAIL", row.get("status", ""))
        else:
            add("summary_status_ok", config, "PASS", "OK")
        for col in REQUIRED_FINITE_SUMMARY_COLUMNS:
            value = row.get(col, "")
            add(f"summary_finite:{col}", config, "PASS" if finite_float(value) else "FAIL", value or "MISSING")

        seed, noise, gamma, td = expected_from_config(config)
        expected_summary = {
            "init_seed": seed,
            "init_noise_scale": noise,
            "d_scale_posterior_gamma": "1",
            "trans_d_scale_posterior_gamma": gamma,
            "trans_dscale_multiplier": td,
            "d_scale_mode": "posterior_count",
            "rho_train_mode": "constant",
            "temperature_start": "1",
            "temperature_end": "1",
            "trans_top1_mstep_mode": "off",
            "trans_gate_mode": "off",
            "input_contact_source": "raw_pairs",
            "uses_phase_labels": "0",
            "uses_charm_or_reference": "0",
        }
        for key, expected in expected_summary.items():
            actual = row.get(key, "")
            ok = close_float_string(actual, expected)
            add(f"summary_value:{key}", config, "PASS" if ok else "FAIL",
                f"actual={actual or 'MISSING'} expected={expected}")

        output_dir = Path(row.get("output_dir", ""))
        manifest_path = output_dir / "p9016_full.manifest.tsv"
        if not manifest_path.exists():
            add("manifest_exists", config, "FAIL", str(manifest_path))
            continue
        add("manifest_exists", config, "PASS", str(manifest_path))
        manifest, duplicates = read_key_value_tsv(manifest_path)
        add("manifest_no_duplicate_keys", config, "PASS" if not duplicates else "FAIL",
            ",".join(sorted(set(duplicates))) or "NA")
        required_manifest = dict(REQUIRED_MANIFEST_VALUES)
        required_manifest.update({
            "config_name": config,
            "init_seed": seed,
            "init_noise_scale": noise,
            "trans_d_scale_posterior_gamma": gamma,
            "trans_dscale_multiplier": td,
            "trans_dscale_effective_count_formula": f"n_raw*posterior_prob^{gamma}",
            "trans_dscale_probability_weighted": "0" if gamma == "0" else "1",
        })
        expected_scaling_mode = "off" if gamma == "1" and td == "1" else "trans_edge_scaling"
        required_manifest["trans_contact_scaling_mode"] = expected_scaling_mode
        for key, expected in required_manifest.items():
            actual = manifest.get(key, "")
            ok = close_float_string(actual, expected)
            add(f"manifest_value:{key}", config, "PASS" if ok else "FAIL",
                f"actual={actual or 'MISSING'} expected={expected}")

    for seed_name, seed in (("seed17", "17"), ("seed71", "71")):
        replay = [row for row in rows if row.get("init_seed") == seed and
                  close_float_string(row.get("trans_d_scale_posterior_gamma", ""), "1") and
                  close_float_string(row.get("trans_dscale_multiplier", ""), "0.5")]
        add("baseline_replay_exactly_one", seed_name,
            "PASS" if len(replay) == 1 else "FAIL", f"count={len(replay)}")

    out_path = root / "strict_matrix_audit.tsv"
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t",
                                fieldnames=["check", "config_name", "status", "detail"],
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit_rows)
    failures = [row for row in audit_rows if row["status"] != "PASS"]
    headline = "PASS" if not failures else "FAIL"
    (root / "strict_matrix_audit_headline.txt").write_text(headline + "\n")
    print(f"STRICT_MATRIX_AUDIT_TSV={out_path}")
    print(f"STRICT_MATRIX_AUDIT_HEADLINE={headline}")
    if failures:
        print(f"STRICT_MATRIX_AUDIT_FAILURES={len(failures)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
