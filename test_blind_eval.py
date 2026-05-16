#!/usr/bin/env python3
"""Unit tests for eval-only blind diploid diagnostics."""

from __future__ import annotations

import gzip
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from blind_eval import (
    PERMUTATIONS,
    deterministic_pair_indices,
    evaluate_chromosome_permutation,
    load_charm_3dg,
    load_model_coords,
    mapped_point_table,
    procrustes_rmsd,
)
from evaluate_blind_p9016_grid import (
    annotate_acceptability,
    baseline_diagnostics,
    best_mapping_for_chrom,
    procrustes_rows,
    posterior_contact_diagnostics,
    read_force_class_diag,
    rg_and_scale_rows,
    rho_train_for_bpair,
    sampled_trans_null_diagnostics,
    trans_chrpair_mean_rows,
    trans_metrics,
    warning_rows,
    write_rankings,
)
from candidate_rerun_diagnostics import (
    DSCALE_ABLATION_ORDER,
    PRIOR_RHO_2X2_ORDER,
    SCAFFOLD_STAGE_ORDER,
    expected_config_order,
    normalize_candidate_summary,
    scaffold_stage_summary_rows,
    validate_candidate_file_schemas,
    write_dscale_diag_output,
    write_force_cancellation_output,
    write_force_normalized_output,
    write_residual_stratified_output,
)


STARTS = [0, 1_000_000, 2_000_000, 3_000_000]
REF_COORDS = {
    ("chr1", "mat"): np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.1], [2.0, 0.4, 0.0], [3.0, 0.0, 0.2]]),
    ("chr1", "pat"): np.array([[0.0, 2.0, 0.0], [1.4, 2.2, 0.2], [1.7, 3.0, 0.8], [3.3, 2.1, 0.1]]),
    ("chr2", "mat"): np.array([[10.0, 0.0, 0.0], [11.0, 0.2, 0.0], [12.0, 0.6, 0.2], [13.0, 0.1, 0.1]]),
    ("chr2", "pat"): np.array([[10.0, 2.0, 0.0], [11.2, 2.4, 0.1], [12.5, 2.0, 0.9], [13.0, 2.8, 0.0]]),
}


def write_charm(path: Path) -> None:
    with gzip.open(path, "wt") if str(path).endswith(".gz") else open(path, "wt") as handle:
        for (chrom, hap), coords in REF_COORDS.items():
            for start, xyz in zip(STARTS, coords):
                handle.write(f"{chrom}({hap})\t{start}\t{xyz[0]}\t{xyz[1]}\t{xyz[2]}\n")


def model_rows(mapping: dict[tuple[str, int], tuple[str, str]], scale: float = 1.0, chr2_extra_shift: float = 0.0) -> list[dict[str, object]]:
    rows = []
    bid = 0
    for chrom in ["chr1", "chr2"]:
        for start_idx, start in enumerate(STARTS):
            for copy in [0, 1]:
                ref_chrom, hap = mapping[(chrom, copy)]
                coords = REF_COORDS[(ref_chrom, hap)].copy()
                center = coords.mean(axis=0)
                xyz = center + scale * (coords[start_idx] - center)
                if chrom == "chr2":
                    xyz = xyz + np.array([chr2_extra_shift, 0.0, 0.0])
                rows.append({
                    "chr": chrom,
                    "start": start,
                    "end": start + 1_000_000,
                    "bid": bid,
                    "copy": copy,
                    "diploid_bid": bid * 2 + copy,
                    "x": xyz[0],
                    "y": xyz[1],
                    "z": xyz[2],
                })
            bid += 1
    return rows


def write_model(path: Path, rows: list[dict[str, object]]) -> None:
    frame = pd.DataFrame(rows)
    if str(path).endswith(".gz"):
        frame.to_csv(path, sep="\t", index=False, compression="gzip")
    else:
        frame.to_csv(path, sep="\t", index=False)


class BlindEvalTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.ref_path = self.root / "tiny.3dg.gz"
        write_charm(self.ref_path)

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def test_load_model_coords_tsv_and_gz(self) -> None:
        rows = model_rows({("chr1", 0): ("chr1", "mat"), ("chr1", 1): ("chr1", "pat"), ("chr2", 0): ("chr2", "mat"), ("chr2", 1): ("chr2", "pat")})
        tsv = self.root / "coords.tsv"
        gz = self.root / "coords.tsv.gz"
        write_model(tsv, rows)
        write_model(gz, rows)
        self.assertEqual(len(load_model_coords(tsv)), len(rows))
        self.assertEqual(len(load_model_coords(gz)), len(rows))

    def test_load_charm_3dg(self) -> None:
        ref = load_charm_3dg(self.ref_path)
        self.assertEqual(set(ref["hap"]), {"mat", "pat"})
        self.assertEqual(set(ref["chr"]), {"chr1", "chr2"})
        six_col = self.root / "tiny_feature.3dg"
        with open(six_col, "wt") as handle:
            handle.write("chr1(mat)\t0\t0\t0\t0\tfeature\n")
            handle.write("chr1(mat)\t1000000\t1\t0\t0\tfeature\n")
            handle.write("chr1(pat)\t0\t0\t1\t0\tfeature\n")
            handle.write("chr1(pat)\t1000000\t1\t1\t0\tfeature\n")
        self.assertEqual(len(load_charm_3dg(six_col)), 4)

    def test_duplicate_and_nan_rejected(self) -> None:
        mapping = {("chr1", 0): ("chr1", "mat"), ("chr1", 1): ("chr1", "pat"), ("chr2", 0): ("chr2", "mat"), ("chr2", 1): ("chr2", "pat")}
        rows = model_rows(mapping)
        dup_path = self.root / "dup.tsv"
        write_model(dup_path, rows + [dict(rows[0])])
        with self.assertRaises(ValueError):
            load_model_coords(dup_path)
        rows[0]["x"] = float("nan")
        nan_path = self.root / "nan.tsv"
        write_model(nan_path, rows)
        with self.assertRaises(ValueError):
            load_model_coords(nan_path)

    def test_missing_copy_is_strict_invalid(self) -> None:
        mapping = {("chr1", 0): ("chr1", "mat"), ("chr1", 1): ("chr1", "pat"), ("chr2", 0): ("chr2", "mat"), ("chr2", 1): ("chr2", "pat")}
        rows = [row for row in model_rows(mapping) if not (row["chr"] == "chr1" and row["copy"] == 1)]
        path = self.root / "missing.tsv"
        write_model(path, rows)
        model = load_model_coords(path)
        ref = load_charm_3dg(self.ref_path)
        result = evaluate_chromosome_permutation(model, ref, "chr1", PERMUTATIONS[0])
        self.assertFalse(result["valid"])

    def test_copy_swap_mapping_selected(self) -> None:
        mapping = {("chr1", 0): ("chr1", "pat"), ("chr1", 1): ("chr1", "mat"), ("chr2", 0): ("chr2", "pat"), ("chr2", 1): ("chr2", "mat")}
        path = self.root / "swapped.tsv"
        write_model(path, model_rows(mapping))
        best = best_mapping_for_chrom(load_model_coords(path), load_charm_3dg(self.ref_path), "chr1")
        self.assertEqual(best["mapping"], "copy0_pat_copy1_mat")
        self.assertGreater(best["min_copy_pearson"], 0.99)

    def test_compaction_catches_scaled_chromosome(self) -> None:
        mapping = {("chr1", 0): ("chr1", "mat"), ("chr1", 1): ("chr1", "pat"), ("chr2", 0): ("chr2", "mat"), ("chr2", 1): ("chr2", "pat")}
        path = self.root / "scaled.tsv"
        write_model(path, model_rows(mapping, scale=0.5))
        model = load_model_coords(path)
        ref = load_charm_3dg(self.ref_path)
        best = [best_mapping_for_chrom(model, ref, "chr1")]
        rg_rows, _, scale_rows = rg_and_scale_rows("scaled", best)
        self.assertTrue(all(abs(row["rg_ratio"] - 0.5) < 1e-6 for row in rg_rows))
        self.assertTrue(all(abs(row["slope"] - 0.5) < 1e-6 for row in scale_rows))

    def test_per_chr_procrustes_detects_bad_global_placement(self) -> None:
        mapping = {("chr1", 0): ("chr1", "mat"), ("chr1", 1): ("chr1", "pat"), ("chr2", 0): ("chr2", "mat"), ("chr2", 1): ("chr2", "pat")}
        path = self.root / "translated.tsv"
        write_model(path, model_rows(mapping, chr2_extra_shift=20.0))
        model = load_model_coords(path)
        ref = load_charm_3dg(self.ref_path)
        per_chr_mapping = {"chr1": "copy0_mat_copy1_pat", "chr2": "copy0_mat_copy1_pat"}
        points = mapped_point_table(model, ref, per_chr_mapping)
        _, summary = procrustes_rows("translated", points)
        self.assertLess(summary["procrustes_per_chr_rmsd_mean"], 1e-6)
        self.assertGreater(summary["procrustes_global_rmsd"], 0.1)

    def test_procrustes_does_not_allow_reflection(self) -> None:
        ref = np.array([[0.0, 0.0, 0.0], [1.0, 0.1, 0.0], [0.2, 1.0, 0.3], [0.1, 0.2, 1.0]])
        reflected = ref * np.array([-1.0, 1.0, 1.0])
        self.assertGreater(procrustes_rmsd(reflected, ref), 0.1)

    def test_trans_sampling_deterministic(self) -> None:
        a1, b1 = deterministic_pair_indices(10, 12, 20, 17)
        a2, b2 = deterministic_pair_indices(10, 12, 20, 17)
        a3, b3 = deterministic_pair_indices(10, 12, 20, 18)
        self.assertTrue(np.array_equal(a1, a2))
        self.assertTrue(np.array_equal(b1, b2))
        self.assertFalse(np.array_equal(a1, a3) and np.array_equal(b1, b3))

    def test_force_class_diag_reader(self) -> None:
        path = self.root / "p9016_full.force_class_diag.tsv"
        path.write_text(
            "class\tn_wedges\tsum_wedge_k\tcontact_energy\tcontact_force_l1\t"
            "backbone_force_l1\trepulsion_force_l1\thomolog_sep_force_l1\tn_nonfinite\n"
            "cis\t2\t3\t5\t7\t0\t0\t0\t0\n"
            "trans\t1\t2\t4\t2\t0\t0\t0\t0\n"
            "total\t3\t5\t9\t9\t11\t13\t17\t0\n"
        )
        summary = read_force_class_diag(path)
        self.assertEqual(summary["force_class_diag_available"], 1)
        self.assertEqual(summary["final_contact_energy_cis"], 5.0)
        self.assertEqual(summary["final_contact_energy_trans"], 4.0)
        self.assertEqual(summary["force_l1_cis_trans_ratio"], 3.5)
        self.assertEqual(summary["final_backbone_force_l1"], 11.0)
        self.assertEqual(summary["final_repulsion_force_l1"], 13.0)
        self.assertEqual(summary["final_homolog_sep_force_l1"], 17.0)
        self.assertEqual(summary["repulsion_to_trans_force_ratio"], 6.5)

    def test_bpair_posterior_schema_and_force_schema(self) -> None:
        posterior = self.root / "p9016_full.bpair_posterior.tsv"
        posterior.write_text(
            "chr1\tstart1\tend1\tchr2\tstart2\tend2\tbid1\tbid2\tn_raw\t"
            "base_d_scale\tbase_k\tp00\tp01\tp10\tp11\tpU\tpsame_raw\tpcross_raw\t"
            "entropy\tmargin\tpmax\trho_output\tcontact_class\n"
            "chr1\t0\t1000000\tchr1\t1000000\t2000000\t0\t1\t8\t0.5\t1\t"
            "0.7\t0.1\t0.1\t0.1\t0.2\t0.8\t0.2\t0.5\t0.6\t0.7\t0.8\tcis\n"
        )
        force = self.root / "p9016_full.force_class_diag.tsv"
        force.write_text(
            "class\tn_wedges\tsum_wedge_k\tcontact_energy\tcontact_force_l1\t"
            "backbone_force_l1\trepulsion_force_l1\thomolog_sep_force_l1\tn_nonfinite\n"
            "cis\t1\t1\t2\t3\t0\t0\t0\t0\n"
            "trans\t1\t1\t4\t5\t0\t0\t0\t0\n"
            "total\t2\t2\t6\t8\t11\t13\t17\t0\n"
        )
        frame = pd.DataFrame([{
            "config_name": "tiny",
            "output_bpair_posterior": str(posterior),
            "output_force_class_diag": str(force),
        }])
        schema = validate_candidate_file_schemas(frame)
        self.assertEqual(int(schema.iloc[0]["posterior_schema_ok"]), 1)
        self.assertEqual(int(schema.iloc[0]["force_schema_ok"]), 1)

    def test_candidate_force_normalized_schema(self) -> None:
        run_dir = self.root / "tiny"
        eval_dir = self.root / "eval" / "tiny"
        run_dir.mkdir()
        eval_dir.mkdir(parents=True)
        force = run_dir / "p9016_full.force_class_diag.tsv"
        force.write_text(
            "class\tn_wedges\tsum_wedge_k\tcontact_energy\tcontact_force_l1\t"
            "backbone_force_l1\trepulsion_force_l1\thomolog_sep_force_l1\tn_nonfinite\n"
            "cis\t4\t3\t5\t8\t0\t0\t0\t0\n"
            "trans\t2\t1\t4\t2\t0\t0\t0\t0\n"
            "total\t6\t4\t9\t10\t11\t13\t17\t0\n"
        )
        (eval_dir / "eval.posterior_cis_trans_summary.tsv").write_text(
            "config_name\tclass\tavailable\tn_bpair\tn_raw_sum\tsum_effective_k\n"
            "tiny\tcis\t1\t2\t16\t3\n"
            "tiny\ttrans\t1\t1\t4\t1\n"
        )
        summary = pd.DataFrame([{
            "config_name": "tiny",
            "output_force_class_diag": str(force),
        }])
        eval_summary = pd.DataFrame([{
            "config_name": "tiny",
            "eval_dir": str(eval_dir),
        }])
        out = write_force_normalized_output(self.root, summary, eval_summary, prefix="prior_rho_2x2")
        expected = {
            "n_bpair_cis",
            "n_raw_cis",
            "n_wedge_cis",
            "sum_effective_k_cis",
            "mean_effective_k_per_bpair_cis",
            "mean_effective_k_per_raw_cis",
            "mean_effective_k_per_wedge_cis",
            "force_l1_cis",
            "force_l1_per_bpair_cis",
            "force_l1_per_raw_cis",
            "force_l1_per_wedge_cis",
            "n_bpair_trans",
            "n_raw_trans",
            "n_wedge_trans",
            "sum_effective_k_trans",
            "mean_effective_k_per_bpair_trans",
            "mean_effective_k_per_raw_trans",
            "mean_effective_k_per_wedge_trans",
            "force_l1_trans",
            "force_l1_per_bpair_trans",
            "force_l1_per_raw_trans",
            "force_l1_per_wedge_trans",
        }
        self.assertTrue(expected.issubset(set(out.columns)))
        self.assertTrue((self.root / "prior_rho_2x2_force_normalized.tsv").exists())
        self.assertAlmostEqual(float(out.iloc[0]["force_l1_per_bpair_cis"]), 4.0)
        self.assertAlmostEqual(float(out.iloc[0]["mean_effective_k_per_wedge_trans"]), 0.5)

    def test_candidate_residual_stratified_schema(self) -> None:
        mapping = {("chr1", 0): ("chr1", "mat"), ("chr1", 1): ("chr1", "pat"), ("chr2", 0): ("chr2", "mat"), ("chr2", 1): ("chr2", "pat")}
        model_path = self.root / "coords.tsv"
        write_model(model_path, model_rows(mapping))
        posterior = self.root / "p9016_full.bpair_posterior.tsv"
        posterior.write_text(
            "chr1\tstart1\tend1\tchr2\tstart2\tend2\tbid1\tbid2\tn_raw\t"
            "base_d_scale\tbase_k\tp00\tp01\tp10\tp11\tpU\tpsame_raw\tpcross_raw\t"
            "entropy\tmargin\tpmax\trho_output\tcontact_class\n"
            "chr1\t0\t1000000\tchr1\t1000000\t2000000\t0\t1\t8\t0.5\t1\t"
            "0.7\t0.1\t0.1\t0.1\t0.2\t0.8\t0.2\t0.5\t0.6\t0.7\t0.8\tcis\n"
            "chr1\t0\t1000000\tchr1\t3000000\t4000000\t0\t3\t4\t0.8\t1\t"
            "0.25\t0.25\t0.25\t0.25\t0.9\t0.5\t0.5\t1.38\t0.0\t0.25\t0.1\tcis\n"
            "chr1\t0\t1000000\tchr2\t0\t1000000\t0\t4\t1\t1\t1\t"
            "0.25\t0.25\t0.25\t0.25\t0.99\t0.5\t0.5\t1.38\t0.0\t0.25\t0.00001\ttrans\n"
        )
        summary = pd.DataFrame([{
            "config_name": "tiny",
            "coords_gz": str(model_path),
            "output_bpair_posterior": str(posterior),
            "rho_train_mode": "entropy",
            "rho_train_end": 1.0,
            "rho_train_floor": 0.5,
            "unit": 1.0,
            "d_scale_mode": "raw_count",
            "contact_k_multiplier_cis": 1.0,
            "contact_k_multiplier_trans": 1.0,
        }])
        out = write_residual_stratified_output(self.root, summary)
        expected = {
            "config_name",
            "stratification",
            "stratum",
            "available",
            "n_bpair",
            "n_raw_sum",
            "mean_pU",
            "median_pmax",
            "median_expected_residual",
            "sum_effective_k",
            "contact_count_vs_inv_expected_distance_spearman",
        }
        self.assertTrue(expected.issubset(set(out.columns)))
        self.assertTrue((self.root / "candidate_rerun_residual_stratified.tsv").exists())
        self.assertIn("cis_short_range", set(out["stratum"]))
        self.assertIn("trans", set(out["stratum"]))
        self.assertIn("n_raw[1,2)", set(out["stratum"]))
        self.assertIn("pU[0.99,1)", set(out["stratum"]))

    def test_candidate_force_cancellation_schema(self) -> None:
        mapping = {("chr1", 0): ("chr1", "mat"), ("chr1", 1): ("chr1", "pat"), ("chr2", 0): ("chr2", "mat"), ("chr2", 1): ("chr2", "pat")}
        model_path = self.root / "coords.tsv"
        write_model(model_path, model_rows(mapping))
        posterior = self.root / "p9016_full.bpair_posterior.tsv"
        posterior.write_text(
            "chr1\tstart1\tend1\tchr2\tstart2\tend2\tbid1\tbid2\tn_raw\t"
            "base_d_scale\tbase_k\tp00\tp01\tp10\tp11\tpU\tpsame_raw\tpcross_raw\t"
            "entropy\tmargin\tpmax\trho_output\tcontact_class\n"
            "chr1\t0\t1000000\tchr1\t1000000\t2000000\t0\t1\t8\t0.5\t1\t"
            "0.7\t0.1\t0.1\t0.1\t0.2\t0.8\t0.2\t0.5\t0.6\t0.7\t0.8\tcis\n"
            "chr1\t0\t1000000\tchr2\t0\t1000000\t0\t4\t1\t1\t1\t"
            "0.25\t0.25\t0.25\t0.25\t0.99\t0.5\t0.5\t1.38\t0.0\t0.25\t0.00001\ttrans\n"
        )
        summary = pd.DataFrame([{
            "config_name": "tiny",
            "coords_gz": str(model_path),
            "output_bpair_posterior": str(posterior),
            "rho_train_mode": "constant",
            "rho_train_end": 1.0,
            "rho_train_floor": 0.5,
            "unit": 1.0,
            "d_scale_mode": "raw_count",
            "contact_k_multiplier_cis": 1.0,
            "contact_k_multiplier_trans": 1.0,
        }])
        out = write_force_cancellation_output(self.root, summary)
        expected = {
            "config_name",
            "available",
            "force_l1_cis",
            "force_l1_trans",
            "net_force_norm_cis",
            "net_force_norm_trans",
            "force_cancellation_ratio_cis",
            "force_cancellation_ratio_trans",
            "force_cancellation_ratio_all",
        }
        self.assertTrue(expected.issubset(set(out.columns)))
        self.assertTrue((self.root / "candidate_rerun_force_cancellation.tsv").exists())
        self.assertEqual(int(out.iloc[0]["available"]), 1)
        self.assertGreaterEqual(float(out.iloc[0]["force_cancellation_ratio_trans"]), 0.0)

    def test_posterior_diagnostics_split_cis_trans(self) -> None:
        mapping = {("chr1", 0): ("chr1", "mat"), ("chr1", 1): ("chr1", "pat"), ("chr2", 0): ("chr2", "mat"), ("chr2", 1): ("chr2", "pat")}
        model_path = self.root / "coords.tsv"
        write_model(model_path, model_rows(mapping))
        model = load_model_coords(model_path)
        posterior = self.root / "posterior.tsv"
        posterior.write_text(
            "chr1\tstart1\tend1\tchr2\tstart2\tend2\tbid1\tbid2\tn_raw\t"
            "base_d_scale\tbase_k\tp00\tp01\tp10\tp11\tpU\tpsame_raw\tpcross_raw\t"
            "entropy\tmargin\tpmax\trho_output\tcontact_class\n"
            "chr1\t0\t1000000\tchr1\t1000000\t2000000\t0\t1\t8\t0.5\t1\t"
            "0.7\t0.1\t0.1\t0.1\t0.2\t0.8\t0.2\t0.5\t0.6\t0.7\t0.8\tcis\n"
            "chr1\t0\t1000000\tchr2\t0\t1000000\t0\t4\t1\t1\t1\t"
            "0.25\t0.25\t0.25\t0.25\t0.99\t0.5\t0.5\t1.38\t0.0\t0.25\t0.00001\ttrans\n"
        )
        summary = pd.Series({"rho_train_mode": "entropy", "rho_train_end": 1.0, "rho_train_floor": 0.5})
        out = posterior_contact_diagnostics("tiny", model, posterior, summary, self.root / "eval")
        self.assertEqual(out["posterior_diag_available"], 1)
        self.assertAlmostEqual(out["mean_pU_cis"], 0.2)
        self.assertAlmostEqual(out["mean_pU_trans"], 0.99)
        self.assertAlmostEqual(out["mean_pmax_cis"], 0.7)
        self.assertAlmostEqual(out["mean_rho_train_cis"], 0.8)
        self.assertAlmostEqual(out["mean_rho_train_trans"], 0.00001)
        self.assertAlmostEqual(out["fraction_trans_contacts_with_rho_train_near_zero"], 1.0)
        self.assertGreater(out["effective_k_cis_trans_ratio"], 1000.0)

    def test_constant_rho_uses_schedule_independent_of_entropy(self) -> None:
        self.assertAlmostEqual(rho_train_for_bpair("constant", 0.7, 0.0, "trans", 0.5), 0.7)
        self.assertAlmostEqual(rho_train_for_bpair("constant", 0.7, 1.0, "cis", 0.5), 0.7)
        self.assertAlmostEqual(rho_train_for_bpair("entropy", 0.7, 0.2, "trans", 0.5), 0.14)

    def test_prior_rho_2x2_expected_config_order(self) -> None:
        order = expected_config_order("prior_rho_2x2")
        self.assertEqual(list(order), PRIOR_RHO_2X2_ORDER)
        self.assertEqual(len(order), 8)
        self.assertIn("prior_rho2x2_n100_rs100_uniform_constant", order)

    def test_scaffold_and_dscale_expected_config_order(self) -> None:
        self.assertEqual(list(expected_config_order("scaffold_stage")), SCAFFOLD_STAGE_ORDER)
        self.assertEqual(len(expected_config_order("scaffold_stage")), 4)
        self.assertEqual(list(expected_config_order("dscale_ablation")), DSCALE_ABLATION_ORDER)
        self.assertEqual(len(expected_config_order("dscale_ablation")), 6)
        self.assertIn("dscale_n80_rs100_expected", expected_config_order("dscale_ablation"))

    def test_scaffold_stage_coord_schema_is_loadable(self) -> None:
        run_dir = self.root / "main_rep1_n100_rs100"
        run_dir.mkdir()
        mapping = {("chr1", 0): ("chr1", "mat"), ("chr1", 1): ("chr1", "pat"), ("chr2", 0): ("chr2", "mat"), ("chr2", 1): ("chr2", "pat")}
        for name in [
            "init_haploid_scaffold.coords.tsv.gz",
            "init_diploid_split.coords.tsv.gz",
            "after_iter1.coords.tsv.gz",
            "final.coords.tsv.gz",
        ]:
            write_model(run_dir / name, model_rows(mapping))
        summary = pd.DataFrame([{
            "config_name": "main_rep1_n100_rs100",
            "output_dir": str(run_dir),
            "coords_gz": str(run_dir / "p9016_full.coords.tsv.gz"),
            "output_bpair_posterior": str(run_dir / "p9016_full.bpair_posterior.tsv"),
            "output_force_class_diag": str(run_dir / "p9016_full.force_class_diag.tsv"),
        }])
        stages = scaffold_stage_summary_rows(summary)
        self.assertEqual(len(stages), 4)
        self.assertEqual(set(stages["stage"]), {"init_haploid_scaffold", "init_diploid_split", "after_iter1", "final"})
        for path in stages["coords_gz"]:
            loaded = load_model_coords(Path(path))
            self.assertEqual(set(loaded.columns), {"chr", "start", "end", "bid", "copy", "diploid_bid", "x", "y", "z"})
            self.assertEqual(len(loaded), 16)

    def test_dscale_diag_raw_and_expected_count_behavior(self) -> None:
        posterior = self.root / "p9016_full.bpair_posterior.tsv"
        posterior.write_text(
            "chr1\tstart1\tend1\tchr2\tstart2\tend2\tbid1\tbid2\tn_raw\t"
            "base_d_scale\tbase_k\tp00\tp01\tp10\tp11\tpU\tpsame_raw\tpcross_raw\t"
            "entropy\tmargin\tpmax\trho_output\tcontact_class\n"
            "chr1\t0\t1000000\tchr1\t1000000\t2000000\t0\t1\t8\t0.5\t1\t"
            "0.999999\t0.000001\t0.0\t0.0\t0.1\t1.0\t0.0\t0.0\t0.999998\t0.999999\t1.0\tcis\n"
        )
        expected_summary = pd.DataFrame([{
            "config_name": "tiny_expected",
            "output_bpair_posterior": str(posterior),
            "d_scale_mode": "expected_count",
            "d_scale_eps_count": 1e-6,
        }])
        expected = write_dscale_diag_output(self.root, expected_summary, prefix="dscale_ablation")
        low = expected.loc[
            (expected["config_name"] == "tiny_expected")
            & (expected["stratification"] == "state_p_bin")
            & (expected["stratum"] == "state_p[0,1e-06)")
        ].iloc[0]
        high = expected.loc[
            (expected["config_name"] == "tiny_expected")
            & (expected["stratification"] == "state_p_bin")
            & (expected["stratum"] == "state_p[0.25,1)")
        ].iloc[0]
        self.assertGreater(float(low["median_state_d_scale"]), float(high["median_state_d_scale"]))
        self.assertTrue(np.isfinite(expected["median_state_d_scale"].dropna().astype(float)).all())

        raw_summary = pd.DataFrame([{
            "config_name": "tiny_raw",
            "output_bpair_posterior": str(posterior),
            "d_scale_mode": "raw_count",
            "d_scale_eps_count": 1e-6,
        }])
        raw = write_dscale_diag_output(self.root, raw_summary, prefix="raw_dscale")
        all_row = raw.loc[(raw["config_name"] == "tiny_raw") & (raw["stratification"] == "all")].iloc[0]
        self.assertAlmostEqual(float(all_row["median_state_d_scale"]), 0.5)
        self.assertTrue((self.root / "raw_dscale_dscale_diag.tsv").exists())

    def test_candidate_summary_normalization_uses_manifest_outputs(self) -> None:
        run_dir = self.root / "main_rep1_n5_rs5"
        run_dir.mkdir()
        (self.root / "matrix_summary.tsv").write_text(
            "config_id\toutput_dir\tinit_mode\tprior_mode\trho_train_mode\td_scale_mode\tn_iter\trelax_steps\taudit_status\n"
            f"main_rep1_n5_rs5\t{run_dir}\tunphased_scaffold_split\tcis_inter_ratio\tentropy\traw_count\t5\t5\tOK\n"
        )
        (run_dir / "p9016_full.manifest.tsv").write_text(
            "key\tvalue\n"
            "sample\tP9016\n"
            "output_dir\t%s\n"
            "output_bpair_posterior\t%s\n"
            "output_coords\t%s\n"
            "output_loop_diag\t%s\n"
            "output_force_class_diag\t%s\n"
            "repulsion_multiplier\t1\n"
            "k_rel_rep_effective\t0.05\n"
            "relax_step\t0.0005\n"
            "temperature_start\t2\n"
            "temperature_end\t1\n"
            "rho_train_start\t1\n"
            "rho_train_end\t1\n"
            "unit\t1\n"
            "d_scale\t1\n"
            "legacy_base_k_unused\t2\n"
            "min_sep_unit\t0.25\n"
            "lambda_sep\t0.05\n"
            "init_eps_effective\t0.5\n"
            "init_noise_scale_effective\t0\n"
            "init_seed\t17\n"
            "init_scale\t0\n"
            "init_mode\tunphased_scaffold_split\n"
            "prior_mode\tcis_inter_ratio\n"
            "rho_train_mode\tentropy\n"
            "d_scale_mode\traw_count\n"
            "scaffold_source\tunphased_fdg_1mb\n"
            "same_bin_filter_enabled\t1\n"
            "n_raw_same_bin_excluded\t0\n"
            "n_bpair_same_bin_excluded\t0\n"
            "posterior_refreshed_after_final_relax\t1\n"
            "status\tOK\n"
            % (
                run_dir,
                run_dir / "p9016_full.bpair_posterior.tsv",
                run_dir / "p9016_full.coords.tsv",
                run_dir / "p9016_full.loop_diag.tsv",
                run_dir / "p9016_full.force_class_diag.tsv",
            )
        )
        normalized = normalize_candidate_summary(self.root)
        self.assertEqual(normalized.iloc[0]["config_name"], "main_rep1_n5_rs5")
        self.assertEqual(normalized.iloc[0]["coords_gz"], str(run_dir / "p9016_full.coords.tsv"))
        self.assertEqual(normalized.iloc[0]["init_seed_effective"], "17")
        self.assertEqual(normalized.iloc[0]["scaffold_source"], "unphased_fdg_1mb")
        with self.assertRaisesRegex(ValueError, "candidate set mismatch"):
            normalize_candidate_summary(self.root, require_exact_candidates=True)

    def test_candidate_rerun_training_source_has_no_reference_reader(self) -> None:
        source = (Path(__file__).resolve().parent / "run_blind_p9016_full_cpu_matrix.c").read_text()
        forbidden = [
            "CHARM",
            "P9016.1m.3dg",
            "/shared/",
            "load_charm",
            "3dg.gz",
            "phase0",
            "phase1",
            "phase[",
            "hk_pair.phase",
            "truth",
            "oracle",
        ]
        for token in forbidden:
            self.assertNotIn(token, source)

    def test_trans_relative_metric_is_not_absolute_pass(self) -> None:
        frame = pd.DataFrame([
            {
                "config_name": "least_bad_trans",
                "trans_relative_metric": 0.35,
                "trans_distance_slope": 0.20,
                "trans_median_distance_ratio": 0.56,
                "centroid_relative_metric": 0.43,
                "posterior_diag_available": 0,
                "all_chr_cis_min_copy_pearson_mean": 0.17,
                "median_rg_ratio": 0.67,
                "median_cis_distance_slope": 0.15,
                "trans_null_permutation_spearman_p95": 0.02,
            },
            {
                "config_name": "worse_trans",
                "trans_relative_metric": 0.10,
                "trans_distance_slope": 0.10,
                "trans_median_distance_ratio": 0.40,
                "centroid_relative_metric": 0.10,
                "posterior_diag_available": 0,
                "all_chr_cis_min_copy_pearson_mean": 0.05,
                "median_rg_ratio": 0.50,
                "median_cis_distance_slope": 0.05,
                "trans_null_permutation_spearman_p95": 0.02,
            },
        ])
        annotated = annotate_acceptability(frame)
        top = annotated.loc[annotated["config_name"] == "least_bad_trans"].iloc[0]
        self.assertEqual(top["trans_metric_rank_best"], 1)
        self.assertEqual(top["trans_metric_absolute_pass"], 0)
        self.assertIn("trans_slope_not_severely_compressed", top["trans_metric_absolute_fail_reasons"])
        warnings = warning_rows(annotated)
        self.assertTrue(any(row["warning_type"] == "relative_top_trans_not_absolute_pass" for row in warnings))

    def test_trans_baselines_and_null_are_eval_only_tables(self) -> None:
        mapping = {("chr1", 0): ("chr1", "mat"), ("chr1", 1): ("chr1", "pat"), ("chr2", 0): ("chr2", "mat"), ("chr2", 1): ("chr2", "pat")}
        path = self.root / "coords.tsv"
        write_model(path, model_rows(mapping))
        model = load_model_coords(path)
        ref = load_charm_3dg(self.ref_path)
        points = mapped_point_table(model, ref, {"chr1": "copy0_mat_copy1_pat", "chr2": "copy0_mat_copy1_pat"})
        trans_rows, trans_summary, payload = trans_metrics("tiny", points, max_pairs=8, sample_seed=17)
        chrpair_rows, chrpair_summary = trans_chrpair_mean_rows("tiny", trans_rows)
        baseline_rows, baseline_summary = baseline_diagnostics("tiny", points, max_pairs=4, sample_seed=17)
        null_rows, null_summary = sampled_trans_null_diagnostics("tiny", payload, chrpair_rows, sample_seed=17)
        self.assertIn("trans_relative_metric", trans_summary)
        self.assertIn("trans_centroid_residual_spearman", trans_summary)
        self.assertIn("trans_chrpair_mean_spearman", chrpair_summary)
        self.assertEqual(
            {row["baseline_name"] for row in baseline_rows},
            {"random_diploid", "shuffled_chromosome_labels", "centroid_only", "compact_ball_per_chromosome", "per_chromosome_translated_scaffold"},
        )
        self.assertIn("baseline_random_diploid_trans_relative_metric", baseline_summary)
        self.assertTrue(any(row["null_type"] == "sampled_trans_ref_distance_permutation" for row in null_rows))
        self.assertIn("trans_null_permutation_spearman_p95", null_summary)

    def test_relative_rank_table_names_replace_legacy_names(self) -> None:
        (self.root / "rank_by_trans.tsv").write_text("stale\n")
        (self.root / "rank_by_centroid.tsv").write_text("stale\n")
        frame = pd.DataFrame([
            {
                "config_name": "a",
                "mode_group": "main",
                "n_iter": 1,
                "relax_steps": 1,
                "multiplier": 1.0,
                "chr1_cis_min_copy_pearson": 0.2,
                "chr1_copy_identity_margin_min_pearson": 0.0,
                "all_chr_cis_min_copy_pearson_mean": 0.3,
                "all_chr_cis_min_copy_spearman_mean": 0.3,
                "all_chr_cis_pooled_spearman": 0.3,
                "trans_relative_metric": 0.1,
                "trans_metric_rank": 2,
                "trans_metric_rank_best": 0,
                "trans_metric_absolute_pass": 0,
                "trans_metric_absolute_fail_reasons": "not_absolute",
                "trans_distance_spearman": 0.1,
                "trans_distance_slope": 0.2,
                "trans_median_distance_ratio": 0.5,
                "centroid_relative_metric": 0.4,
                "centroid_metric_rank": 1,
                "centroid_metric_rank_best": 1,
                "centroid_metric_absolute_pass": 0,
                "centroid_metric_absolute_fail_reasons": "not_absolute",
                "centroid_distance_spearman": 0.4,
                "centroid_distance_slope": 0.2,
                "centroid_median_distance_ratio": 0.5,
                "median_rg_ratio": 0.8,
                "median_cis_distance_slope": 0.2,
                "procrustes_global_over_per_chr_ratio": 2.0,
                "eval_dir": "eval/a",
            },
            {
                "config_name": "b",
                "mode_group": "main",
                "n_iter": 2,
                "relax_steps": 2,
                "multiplier": 1.0,
                "chr1_cis_min_copy_pearson": 0.1,
                "chr1_copy_identity_margin_min_pearson": 0.0,
                "all_chr_cis_min_copy_pearson_mean": 0.2,
                "all_chr_cis_min_copy_spearman_mean": 0.2,
                "all_chr_cis_pooled_spearman": 0.2,
                "trans_relative_metric": 0.5,
                "trans_metric_rank": 1,
                "trans_metric_rank_best": 1,
                "trans_metric_absolute_pass": 0,
                "trans_metric_absolute_fail_reasons": "not_absolute",
                "trans_distance_spearman": 0.5,
                "trans_distance_slope": 0.2,
                "trans_median_distance_ratio": 0.5,
                "centroid_relative_metric": 0.2,
                "centroid_metric_rank": 2,
                "centroid_metric_rank_best": 0,
                "centroid_metric_absolute_pass": 0,
                "centroid_metric_absolute_fail_reasons": "not_absolute",
                "centroid_distance_spearman": 0.2,
                "centroid_distance_slope": 0.2,
                "centroid_median_distance_ratio": 0.5,
                "median_rg_ratio": 0.9,
                "median_cis_distance_slope": 0.2,
                "procrustes_global_over_per_chr_ratio": 2.0,
                "eval_dir": "eval/b",
            },
        ])
        write_rankings(self.root, frame)
        self.assertFalse((self.root / "rank_by_trans.tsv").exists())
        self.assertFalse((self.root / "rank_by_centroid.tsv").exists())
        trans_rank = pd.read_csv(self.root / "rank_by_trans_relative_metric.tsv", sep="\t")
        centroid_rank = pd.read_csv(self.root / "rank_by_centroid_relative_metric.tsv", sep="\t")
        self.assertEqual(trans_rank.iloc[0]["config_name"], "b")
        self.assertEqual(centroid_rank.iloc[0]["config_name"], "a")
        self.assertEqual(trans_rank.iloc[0]["rank_metric_kind"], "relative_spearman")
        self.assertEqual(centroid_rank.iloc[0]["rank_metric_kind"], "relative_spearman")


if __name__ == "__main__":
    unittest.main()
