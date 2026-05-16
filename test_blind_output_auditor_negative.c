#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "hickit.h"

#define HK_NEG_N_RAW 2
#define HK_NEG_N_BPAIR 2
#define HK_NEG_N_BEADS 8
#define HK_NEG_OUTPUT_MAX 32768

enum corruption_case {
	CORRUPT_NONE = 0,
	CORRUPT_BPAIR_P4_SUM,
	CORRUPT_BPAIR_FIVE_STATE,
	CORRUPT_COORD_DIPLOID_BID,
	CORRUPT_COORD_MISSING_COPY,
	CORRUPT_RAW_SWAPPED,
	CORRUPT_RAW_P4_SUM,
	CORRUPT_FORBIDDEN_STRING,
	CORRUPT_MANIFEST_NUMERIC,
	CORRUPT_BPAIR_BASE_D_SCALE,
	CORRUPT_BPAIR_BASE_K
};

struct audit_run {
	int status;
	char output[HK_NEG_OUTPUT_MAX];
};

struct negative_case {
	enum corruption_case kind;
	const char *name;
	int write_raw;
	const char *expect1;
	const char *expect2;
};

static int check_true(const char *label, int pred)
{
	if (!pred) {
		fprintf(stderr, "%s: predicate failed\n", label);
		return 1;
	}
	return 0;
}

static int file_exists(const char *path)
{
	FILE *fp = fopen(path, "rb");
	if (fp == 0) return 0;
	fclose(fp);
	return 1;
}

static void path_join(char *dst, size_t dst_size, const char *dir, const char *name)
{
	size_t n = strlen(dir);
	if (n > 0 && dir[n - 1] == '/')
		snprintf(dst, dst_size, "%s%s", dir, name);
	else
		snprintf(dst, dst_size, "%s/%s", dir, name);
}

static int make_case_dir(char *dir, size_t dir_size, const char *case_name)
{
	int i;
	for (i = 0; i < 100; ++i) {
		snprintf(dir, dir_size, "/tmp/hk_blind_auditor_negative_%ld_%s_%d",
				 (long)getpid(), case_name, i);
		if (mkdir(dir, 0700) == 0)
			return 0;
		if (errno != EEXIST)
			return -1;
	}
	return -1;
}

static void remove_case_dir(const char *dir)
{
	char path[1024];
	if (dir == 0 || dir[0] == 0) return;
	path_join(path, sizeof(path), dir, "p9016_full.manifest.tsv");
	remove(path);
	path_join(path, sizeof(path), dir, "p9016_full.bpair_posterior.tsv");
	remove(path);
	path_join(path, sizeof(path), dir, "p9016_full.coords.tsv");
	remove(path);
	path_join(path, sizeof(path), dir, "p9016_full.loop_diag.tsv");
	remove(path);
	path_join(path, sizeof(path), dir, "p9016_full.raw_posterior.tsv");
	remove(path);
	rmdir(dir);
}

static int write_manifest(const char *path, const char *dir, int write_raw,
						  enum corruption_case kind)
{
	FILE *fp = fopen(path, "w");
	char bpair_path[1024], coords_path[1024], diag_path[1024], raw_path[1024];
	if (fp == 0) return -1;
	path_join(bpair_path, sizeof(bpair_path), dir, "p9016_full.bpair_posterior.tsv");
	path_join(coords_path, sizeof(coords_path), dir, "p9016_full.coords.tsv");
	path_join(diag_path, sizeof(diag_path), dir, "p9016_full.loop_diag.tsv");
	path_join(raw_path, sizeof(raw_path), dir, "p9016_full.raw_posterior.tsv");
	if (fprintf(fp,
				"key\tvalue\n"
				"sample\tP9016\n"
				"runner_family\ttest_fixture\n"
				"runner_version\t2026-04-30\n"
				"default_profile\tp9016_auditor_negative_fixture_v1\n"
				"input_path\t../pairs/P9016.pairs.gz\n"
				"output_dir\t%s\n",
				dir) < 0) goto fail;
	if (kind == CORRUPT_MANIFEST_NUMERIC) {
		if (fprintf(fp, "n_raw\tnot_a_number\n") < 0) goto fail;
	} else {
		if (fprintf(fp, "n_raw\t%d\n", HK_NEG_N_RAW) < 0) goto fail;
	}
	if (fprintf(fp,
				"n_bpair\t%d\n"
				"n_beads\t%d\n"
				"resolution\t1000000\n"
				"n_iter\t3\n"
				"unit\t1\n"
				"d_scale\t1\n"
				"base_k_mode\tuniform\n"
				"base_k_effective\t1\n"
				"base_k_min\t1\n"
				"base_k_mean\t1\n"
				"base_k_max\t1\n"
				"base_k_n_nonfinite\t0\n"
				"legacy_base_k_unused\t2\n"
				"init_mode\tunphased_scaffold_split\n"
				"init_scale\t1\n"
				"prior_mode\tuniform\n"
				"prior_eps\t1e-06\n"
				"prior_inter_density\t0\n"
				"prior_observed_inter\t0\n"
				"prior_possible_inter\t0\n"
				"prior_n_distance_bins\t0\n"
				"prior_n_alpha\t0\n"
				"prior_alpha_min\t0\n"
				"prior_alpha_median\t0\n"
				"prior_alpha_max\t0\n"
				"prior_smoothing_method\tnone\n"
				"prior_alpha_clamp_min\t1e-06\n"
				"prior_alpha_clamp_max\t0.5\n"
				"rho_train_mode\tconstant\n"
				"d_scale_mode\traw_count\n"
				"d_scale_eps_count\t1e-06\n"
				"same_bin_filter_enabled\t1\n"
				"n_raw_same_bin_excluded\t0\n"
				"n_bpair_same_bin_excluded\t0\n"
				"raw_posterior_same_bin_policy\tuniform_unknown_rows\n"
				"min_sep_unit\t0.25\n"
				"lambda_sep\t0.05\n"
				"relax_step\t0.001\n"
				"relax_steps\t5\n"
				"temperature_start\t2\n"
				"temperature_end\t1\n"
				"rho_train_start\t1\n"
				"rho_train_end\t1\n"
				"init_eps_effective\t0.5\n"
				"init_noise_scale_effective\t0\n"
				"init_split_params_used\t1\n"
				"init_seed\t17\n"
				"enable_repulsion\t1\n"
				"repulsion_mode\t2\n"
				"repulsion_blocking_mode\tcurrent_edge_blocking\n"
				"write_raw_posterior\t%d\n"
				"output_bpair_posterior\t%s\n"
				"output_coords\t%s\n"
				"output_loop_diag\t%s\n",
				HK_NEG_N_BPAIR, HK_NEG_N_BEADS, write_raw,
				bpair_path, coords_path, diag_path) < 0) goto fail;
	if (write_raw) {
		if (fprintf(fp, "output_raw_posterior\t%s\n", raw_path) < 0) goto fail;
	}
	if (fprintf(fp,
				"final_mean_entropy\t1.0\n"
				"final_mean_pU\t0.2\n"
				"final_sum_wedge_k\t3.5\n"
				"final_repulsion_energy\t0.75\n"
				"posterior_refreshed_after_final_relax\t1\n"
				"posterior_refresh_temperature\t1\n"
				"posterior_refresh_prior_mode\tuniform\n"
				"posterior_refresh_mean_kl\t0\n"
				"posterior_refresh_top_state_switch_frac\t0\n"
				"posterior_refresh_mean_pU_before\t0.2\n"
				"posterior_refresh_mean_pU_after\t0.2\n"
				"n_bad_iter\t0\n"
				"n_relax_nonfinite_iter\t0\n"
				"n_coord_nonfinite\t0\n"
				"status\tOK\n") < 0) goto fail;
	if (fclose(fp) != 0) return -1;
	return 0;

fail:
	fclose(fp);
	return -1;
}

static void bpair_probs(enum corruption_case kind, int row,
						double *p00, double *p01, double *p10, double *p11,
						double *pU, double *rho)
{
	*p00 = 0.25; *p01 = 0.25; *p10 = 0.25; *p11 = 0.25;
	*pU = 0.20; *rho = 0.80;
	if (row == 0 && kind == CORRUPT_BPAIR_P4_SUM) {
		*p00 = 0.50; *p01 = 0.50; *p10 = 0.50; *p11 = 0.50;
	}
	if (row == 0 && kind == CORRUPT_BPAIR_FIVE_STATE) {
		*pU = 0.50;
	}
}

static int write_bpair_file(const char *path, enum corruption_case kind)
{
	FILE *fp = fopen(path, "w");
	int row;
	if (fp == 0) return -1;
	if (fprintf(fp,
				"chr1\tstart1\tend1\tchr2\tstart2\tend2\tbid1\tbid2\tn_raw\t"
				"base_d_scale\tbase_k\tp00\tp01\tp10\tp11\tpU\tpsame_raw\tpcross_raw\t"
				"entropy\tmargin\tpmax\trho_output") < 0)
		goto fail;
	if (kind == CORRUPT_FORBIDDEN_STRING) {
		if (fprintf(fp, "\tphase0") < 0) goto fail;
	}
	if (fprintf(fp, "\n") < 0) goto fail;
	for (row = 0; row < HK_NEG_N_BPAIR; ++row) {
		double p00, p01, p10, p11, pU, rho;
		int n_raw = row == 0? 1 : 8;
		double base_d_scale = pow((double)n_raw, -1.0 / 3.0);
		double base_k = 1.0;
		int bid1 = row == 0? 0 : 2;
		int bid2 = row == 0? 1 : 7;
		bpair_probs(kind, row, &p00, &p01, &p10, &p11, &pU, &rho);
		if (kind == CORRUPT_BPAIR_BASE_D_SCALE && row == 1)
			base_d_scale = 1.0;
		if (kind == CORRUPT_BPAIR_BASE_K && row == 0)
			base_k = 2.0;
		if (fprintf(fp,
					"chr1\t%d\t%d\tchr1\t%d\t%d\t%d\t%d\t%d\t"
					"%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t1.38629436\t0\t0.25\t%.9g\n",
					bid1 * 1000000, bid1 * 1000000 + 1000000,
					bid2 * 1000000, bid2 * 1000000 + 1000000,
					bid1, bid2, n_raw, base_d_scale, base_k, p00, p01, p10, p11, pU,
					p00 + p11, p01 + p10, rho) < 0)
			goto fail;
	}
	if (fclose(fp) != 0) return -1;
	return 0;

fail:
	fclose(fp);
	return -1;
}

static int write_coords_file(const char *path, enum corruption_case kind)
{
	FILE *fp = fopen(path, "w");
	int bid, copy;
	if (fp == 0) return -1;
	if (fprintf(fp, "chr\tstart\tend\tbid\tcopy\tdiploid_bid\tx\ty\tz\n") < 0)
		goto fail;
	for (bid = 0; bid < HK_NEG_N_BEADS; ++bid) {
		for (copy = 0; copy < HK_DIPLOID_N_COPY; ++copy) {
			int diploid_bid;
			if (kind == CORRUPT_COORD_MISSING_COPY && bid == 7 && copy == 1)
				continue;
			diploid_bid = hk_diploid_bid(bid, copy);
			if (kind == CORRUPT_COORD_DIPLOID_BID && bid == 7 && copy == 1)
				diploid_bid = 999;
			if (fprintf(fp, "chr1\t%d\t%d\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\n",
						bid * 1000000, bid * 1000000 + 1000000,
						bid, copy, diploid_bid,
						0.1 * (double)bid, 0.2 * (double)copy,
						0.01 * (double)(bid + copy)) < 0)
				goto fail;
		}
	}
	if (fclose(fp) != 0) return -1;
	return 0;

fail:
	fclose(fp);
	return -1;
}

static int write_loop_diag_file(const char *path)
{
	FILE *fp = fopen(path, "w");
	if (fp == 0) return -1;
	if (fprintf(fp,
				"n_iter\tn_completed\tinitial_mean_entropy\tfinal_mean_entropy\t"
				"initial_mean_pU\tfinal_mean_pU\tinitial_temperature\tfinal_temperature\t"
				"initial_rho_train\tfinal_rho_train\ttotal_chr_flipped\tn_bad_iter\t"
				"n_relax_nonfinite_iter\tn_coord_nonfinite\tfinal_mean_sep\tfinal_min_sep\t"
				"final_max_sep\tfinal_sum_wedge_k\tfinal_mean_rho_train_bpair\t"
				"final_min_rho_train_bpair\tfinal_max_rho_train_bpair\t"
				"final_n_skipped_same_bin_bpairs\tfinal_repulsion_energy\t"
				"final_repulsion_force_l1\tn_repulsion_nonfinite_step\trepulsion_mode\t"
				"posterior_refreshed_after_final_relax\tposterior_refresh_temperature\t"
				"posterior_refresh_mean_kl\tposterior_refresh_top_state_switch_frac\t"
				"posterior_refresh_mean_pU_before\tposterior_refresh_mean_pU_after\n"
				"3\t3\t1.1\t1.0\t0.25\t0.2\t2\t1\t1\t1\t0\t0\t0\t0\t1.0\t0.5\t1.5\t3.5\t1\t1\t1\t0\t0.75\t1.25\t0\t2\t1\t1\t0\t0\t0.2\t0.2\n") < 0) {
		fclose(fp);
		return -1;
	}
	if (fclose(fp) != 0) return -1;
	return 0;
}

static void raw_probs(enum corruption_case kind, int row,
					  double *p00, double *p01, double *p10, double *p11,
					  double *pU, double *rho)
{
	*p00 = 0.25; *p01 = 0.25; *p10 = 0.25; *p11 = 0.25;
	*pU = 0.20; *rho = 0.80;
	if (row == 0 && kind == CORRUPT_RAW_P4_SUM) {
		*p00 = 0.50; *p01 = 0.50; *p10 = 0.50; *p11 = 0.50;
	}
}

static int write_raw_file(const char *path, enum corruption_case kind)
{
	FILE *fp = fopen(path, "w");
	int row;
	if (fp == 0) return -1;
	if (fprintf(fp,
				"raw_id\tchr1\tpos1\tchr2\tpos2\tbid1_raw\tbid2_raw\tbpair_id\t"
				"bid1_canonical\tbid2_canonical\tswapped\tp00\tp01\tp10\tp11\tpU\t"
				"psame_raw\tpcross_raw\tentropy\tmargin\tpmax\trho_output\n") < 0)
		goto fail;
	for (row = 0; row < HK_NEG_N_RAW; ++row) {
		double p00, p01, p10, p11, pU, rho;
		int bpair_id = row;
		int bid1_can = row == 0? 0 : 2;
		int bid2_can = row == 0? 1 : 7;
		int bid1_raw = row == 0? bid1_can : bid2_can;
		int bid2_raw = row == 0? bid2_can : bid1_can;
		int swapped = row == 0? 0 : 1;
		if (kind == CORRUPT_RAW_SWAPPED && row == 0)
			swapped = 3;
		raw_probs(kind, row, &p00, &p01, &p10, &p11, &pU, &rho);
		if (fprintf(fp,
					"%d\tchr1\t%d\tchr1\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t"
					"%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t1.38629436\t0\t0.25\t%.9g\n",
					row, bid1_raw * 1000000 + 123, bid2_raw * 1000000 + 456,
					bid1_raw, bid2_raw, bpair_id, bid1_can, bid2_can, swapped,
					p00, p01, p10, p11, pU, p00 + p11, p01 + p10, rho) < 0)
			goto fail;
	}
	if (fclose(fp) != 0) return -1;
	return 0;

fail:
	fclose(fp);
	return -1;
}

static int write_fixture_dir(const char *dir, int write_raw, enum corruption_case kind)
{
	char manifest_path[1024], bpair_path[1024], coords_path[1024], diag_path[1024], raw_path[1024];
	path_join(manifest_path, sizeof(manifest_path), dir, "p9016_full.manifest.tsv");
	path_join(bpair_path, sizeof(bpair_path), dir, "p9016_full.bpair_posterior.tsv");
	path_join(coords_path, sizeof(coords_path), dir, "p9016_full.coords.tsv");
	path_join(diag_path, sizeof(diag_path), dir, "p9016_full.loop_diag.tsv");
	path_join(raw_path, sizeof(raw_path), dir, "p9016_full.raw_posterior.tsv");
	if (write_manifest(manifest_path, dir, write_raw, kind) != 0) return -1;
	if (write_bpair_file(bpair_path, kind) != 0) return -1;
	if (write_coords_file(coords_path, kind) != 0) return -1;
	if (write_loop_diag_file(diag_path) != 0) return -1;
	if (write_raw && write_raw_file(raw_path, kind) != 0) return -1;
	return 0;
}

static int run_auditor_capture(const char *dir, struct audit_run *run)
{
	char cmd[1400];
	FILE *fp;
	size_t used = 0;
	int c;

	memset(run, 0, sizeof(*run));
	snprintf(cmd, sizeof(cmd), "./audit_blind_p9016_full_cpu_output.bin \"%s\" 2>&1", dir);
	fp = popen(cmd, "r");
	if (fp == 0) {
		snprintf(run->output, sizeof(run->output), "popen failed");
		run->status = -1;
		return -1;
	}
	while ((c = fgetc(fp)) != EOF) {
		if (used + 1 < sizeof(run->output))
			run->output[used++] = (char)c;
	}
	run->output[used] = 0;
	run->status = pclose(fp);
	return 0;
}

static int contains_expected(const char *output, const char *a, const char *b)
{
	if (a != 0 && strstr(output, a) != 0) return 1;
	if (b != 0 && strstr(output, b) != 0) return 1;
	return 0;
}

static int run_valid_control(void)
{
	char dir[512] = {0};
	struct audit_run run;
	int failed = 0;

	failed |= check_true("auditor negative valid temp dir",
						 make_case_dir(dir, sizeof(dir), "valid") == 0);
	if (failed) goto cleanup;
	failed |= check_true("auditor negative valid fixture",
						 write_fixture_dir(dir, 1, CORRUPT_NONE) == 0);
	if (failed) goto cleanup;
	failed |= check_true("auditor negative valid capture", run_auditor_capture(dir, &run) == 0);
	failed |= check_true("auditor negative valid status", run.status == 0);
	failed |= check_true("auditor negative valid OK", strstr(run.output, "audit status: OK") != 0);
	if (failed)
		fprintf(stderr, "valid auditor output:\n%s\n", run.output);

cleanup:
	remove_case_dir(dir);
	return failed != 0;
}

static int run_negative_case(const struct negative_case *tc)
{
	char dir[512] = {0};
	struct audit_run run;
	int failed = 0;

	failed |= check_true("auditor negative temp dir",
						 make_case_dir(dir, sizeof(dir), tc->name) == 0);
	if (failed) goto cleanup;
	failed |= check_true("auditor negative fixture",
						 write_fixture_dir(dir, tc->write_raw, tc->kind) == 0);
	if (failed) goto cleanup;
	failed |= check_true("auditor negative capture", run_auditor_capture(dir, &run) == 0);
	failed |= check_true("auditor negative nonzero status", run.status != 0);
	failed |= check_true("auditor negative status fail",
						 strstr(run.output, "audit status: FAIL") != 0);
	failed |= check_true("auditor negative keyword",
						 contains_expected(run.output, tc->expect1, tc->expect2));
	if (failed) {
		fprintf(stderr, "negative case %s auditor output:\n%s\n", tc->name, run.output);
	} else {
		fprintf(stderr, "auditor negative case %s: observed expected failure keyword\n", tc->name);
	}

cleanup:
	remove_case_dir(dir);
	return failed != 0;
}

int main(void)
{
	static const struct negative_case cases[] = {
		{ CORRUPT_BPAIR_P4_SUM, "bpair_p4_sum", 0, "p4 sum", "posterior_bad_rows" },
		{ CORRUPT_BPAIR_FIVE_STATE, "bpair_five_state", 0, "five-state", "posterior_bad_rows" },
		{ CORRUPT_COORD_DIPLOID_BID, "coord_diploid_bid", 0, "diploid_bid", "coord_bad_rows" },
		{ CORRUPT_COORD_MISSING_COPY, "coord_missing_copy", 0, "copy mask", "row count" },
		{ CORRUPT_RAW_SWAPPED, "raw_swapped", 1, "swapped", "raw_bad_rows" },
		{ CORRUPT_RAW_P4_SUM, "raw_p4_sum", 1, "p4 sum", "raw_bad_rows" },
		{ CORRUPT_FORBIDDEN_STRING, "forbidden_string", 0, "forbidden string", "forbidden_string_hits" },
		{ CORRUPT_MANIFEST_NUMERIC, "manifest_numeric", 0, "bad numeric value", "manifest" },
		{ CORRUPT_BPAIR_BASE_D_SCALE, "bpair_base_d_scale", 0, "base_d_scale", "posterior_bad_rows" },
		{ CORRUPT_BPAIR_BASE_K, "bpair_base_k", 0, "base_k", "posterior_bad_rows" }
	};
	int failed = 0;
	size_t i;

	if (!file_exists("./audit_blind_p9016_full_cpu_output.bin")) {
		fprintf(stderr, "auditor binary missing; build audit_blind_p9016_full_cpu_output.bin first\n");
		return 1;
	}
	failed |= run_valid_control();
	for (i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
		failed |= run_negative_case(&cases[i]);
	return failed != 0;
}
