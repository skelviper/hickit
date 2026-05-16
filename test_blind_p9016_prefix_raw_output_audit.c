#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "hickit.h"

#define HK_BLIND_P9016_RAW_AUDIT_PATH "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_RAW_AUDIT_PREFIX_N 10000
#define HK_BLIND_P9016_RAW_AUDIT_RESOLUTION 1000000
#define HK_BLIND_P9016_RAW_AUDIT_N_ITER 3
#define HK_BLIND_P9016_RAW_AUDIT_UNIT 1.0f
#define HK_BLIND_P9016_RAW_AUDIT_D_SCALE 1.0f
#define HK_BLIND_P9016_RAW_AUDIT_LEGACY_BASE_K_UNUSED 2.0f
#define HK_BLIND_P9016_RAW_AUDIT_TEMPERATURE_START 2.0f
#define HK_BLIND_P9016_RAW_AUDIT_TEMPERATURE_END 1.0f
#define HK_BLIND_P9016_RAW_AUDIT_RHO_TRAIN_START 1.0f
#define HK_BLIND_P9016_RAW_AUDIT_RHO_TRAIN_END 1.0f
#define HK_BLIND_P9016_RAW_AUDIT_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_RAW_AUDIT_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_RAW_AUDIT_STEP 0.001f
#define HK_BLIND_P9016_RAW_AUDIT_RELAX_STEPS 5
#define HK_BLIND_P9016_RAW_AUDIT_INIT_EPS 0.5f
#define HK_BLIND_P9016_RAW_AUDIT_INIT_NOISE_SCALE 0.0f
#define HK_BLIND_P9016_RAW_AUDIT_INIT_SEED 17ULL

static int check_i32(const char *label, int32_t got, int32_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %d, expected %d\n", label, got, expected);
		return 1;
	}
	return 0;
}

static int check_close(const char *label, double got, double expected)
{
	double tol = 1e-5;
	double scale = fabs(expected) > 1.0? fabs(expected) : 1.0;
	if (fabs(got - expected) > tol * scale) {
		fprintf(stderr, "%s: got %.12g, expected %.12g\n", label, got, expected);
		return 1;
	}
	return 0;
}

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

static long file_size_or_negative(const char *path)
{
	struct stat st;
	if (stat(path, &st) != 0) return -1;
	return (long)st.st_size;
}

static int count_data_rows(const char *path, int64_t *out_rows)
{
	FILE *fp = fopen(path, "r");
	char buf[8192];
	int64_t rows = -1; /* Header line makes the first newline non-data. */
	int saw_any = 0;
	if (fp == 0) return -1;
	while (fgets(buf, sizeof(buf), fp)) {
		saw_any = 1;
		++rows;
	}
	if (ferror(fp)) {
		fclose(fp);
		return -1;
	}
	fclose(fp);
	if (!saw_any) rows = 0;
	if (rows < 0) rows = 0;
	*out_rows = rows;
	return 0;
}

static int make_temp_dir(char *dir, size_t dir_size)
{
	int i;
	for (i = 0; i < 100; ++i) {
		snprintf(dir, dir_size, "/tmp/hk_blind_p9016_raw_audit_%ld_%d", (long)getpid(), i);
		if (mkdir(dir, 0700) == 0)
			return 0;
		if (errno != EEXIST)
			return -1;
	}
	return -1;
}

static void remove_outputs(const char *out_dir)
{
	char path[768];
	if (out_dir == 0 || out_dir[0] == 0) return;
	snprintf(path, sizeof(path), "%s/p9016_full.bpair_posterior.tsv", out_dir);
	remove(path);
	snprintf(path, sizeof(path), "%s/p9016_full.raw_posterior.tsv", out_dir);
	remove(path);
	snprintf(path, sizeof(path), "%s/p9016_full.coords.tsv", out_dir);
	remove(path);
	snprintf(path, sizeof(path), "%s/p9016_full.loop_diag.tsv", out_dir);
	remove(path);
	snprintf(path, sizeof(path), "%s/p9016_full.manifest.tsv", out_dir);
	remove(path);
	rmdir(out_dir);
}

static void set_haploid_scaffold(const struct hk_bmap *bmap, fvec3_t *haploid)
{
	int32_t i;
	for (i = 0; i < bmap->n_beads; ++i) {
		haploid[i][0] = 0.10f * (float)(i % 97);
		haploid[i][1] = 0.07f * (float)((i / 97) % 97);
		haploid[i][2] = 0.03f * (float)(i % 17);
	}
}

static void copy_coords(fvec3_t *dst, const fvec3_t *src, int32_t n)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			dst[i][a] = src[i][a];
}

static int coords_any_changed(const fvec3_t *before, const fvec3_t *after, int32_t n, float tol)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			if (fabsf(after[i][a] - before[i][a]) > tol)
				return 1;
	return 0;
}

static int check_coords_finite(const fvec3_t *coords, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("raw audit coord finite", isfinite(coords[i][a]));
	return failed;
}

static void set_schedule_conf(struct hk_blind_iter_schedule_conf *conf)
{
	conf->n_iter = HK_BLIND_P9016_RAW_AUDIT_N_ITER;
	conf->base_conf.unit = HK_BLIND_P9016_RAW_AUDIT_UNIT;
	conf->base_conf.d_scale = HK_BLIND_P9016_RAW_AUDIT_D_SCALE;
	conf->base_conf.base_k = HK_BLIND_P9016_RAW_AUDIT_LEGACY_BASE_K_UNUSED;
	conf->base_conf.temperature = HK_BLIND_P9016_RAW_AUDIT_TEMPERATURE_START;
	conf->base_conf.rho_train = HK_BLIND_P9016_RAW_AUDIT_RHO_TRAIN_START;
	conf->base_conf.min_sep_unit = HK_BLIND_P9016_RAW_AUDIT_MIN_SEP_UNIT;
	conf->base_conf.lambda_sep = HK_BLIND_P9016_RAW_AUDIT_LAMBDA_SEP;
	conf->base_conf.relax_step = HK_BLIND_P9016_RAW_AUDIT_STEP;
	conf->base_conf.relax_steps = HK_BLIND_P9016_RAW_AUDIT_RELAX_STEPS;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->base_conf.d_scale_eps_count = 1e-6f;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = HK_BLIND_P9016_RAW_AUDIT_TEMPERATURE_START;
	conf->temperature_end = HK_BLIND_P9016_RAW_AUDIT_TEMPERATURE_END;
	conf->rho_train_start = HK_BLIND_P9016_RAW_AUDIT_RHO_TRAIN_START;
	conf->rho_train_end = HK_BLIND_P9016_RAW_AUDIT_RHO_TRAIN_END;
}

static int check_binned_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
							int32_t max_raw)
{
	int failed = 0;
	int32_t i;

	failed |= check_true("raw audit set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_true("raw audit raw positive", set->n_raw > 0);
	failed |= check_true("raw audit raw bounded", set->n_raw <= max_raw);
	failed |= check_true("raw audit bpair positive", set->n_bpairs > 0);
	for (i = 0; i < set->n_raw; ++i) {
		failed |= check_true("raw audit raw2binned id lower", set->raw2binned[i].bpair_id >= 0);
		failed |= check_true("raw audit raw2binned id upper", set->raw2binned[i].bpair_id < set->n_bpairs);
		failed |= check_true("raw audit raw2binned swapped",
							 set->raw2binned[i].swapped == 0 || set->raw2binned[i].swapped == 1);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		failed |= check_true("raw audit key sorted", bp->key.bid[0] <= bp->key.bid[1]);
		failed |= check_true("raw audit key bid0 valid", bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		failed |= check_true("raw audit key bid1 valid", bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
	}
	return failed;
}

static int check_loop_health(const struct hk_blind_iter_loop_diag *diag,
							 const struct hk_blind_single_iter_diag *per_iter)
{
	int failed = 0;
	int32_t i;

	failed |= check_i32("raw audit loop n_iter", diag->n_iter, HK_BLIND_P9016_RAW_AUDIT_N_ITER);
	failed |= check_i32("raw audit loop completed", diag->n_completed, HK_BLIND_P9016_RAW_AUDIT_N_ITER);
	failed |= check_i32("raw audit loop bad iter", diag->n_bad_iter, 0);
	failed |= check_i32("raw audit loop relax bad", diag->n_relax_nonfinite_iter, 0);
	failed |= check_i32("raw audit loop coord bad", diag->n_coord_nonfinite, 0);
	failed |= check_close("raw audit initial temperature", diag->initial_temperature,
						  HK_BLIND_P9016_RAW_AUDIT_TEMPERATURE_START);
	failed |= check_close("raw audit final temperature", diag->final_temperature,
						  HK_BLIND_P9016_RAW_AUDIT_TEMPERATURE_END);
	failed |= check_close("raw audit initial rho", diag->initial_rho_train,
						  HK_BLIND_P9016_RAW_AUDIT_RHO_TRAIN_START);
	failed |= check_close("raw audit final rho", diag->final_rho_train,
						  HK_BLIND_P9016_RAW_AUDIT_RHO_TRAIN_END);
	failed |= check_true("raw audit final entropy finite", isfinite(diag->final_mean_entropy));
	failed |= check_true("raw audit final entropy range", diag->final_mean_entropy >= -1e-6f &&
						 diag->final_mean_entropy <= logf((float)HK_BLIND_N_STATE) + 1e-6f);
	failed |= check_true("raw audit final pU finite", isfinite(diag->final_mean_pU));
	failed |= check_true("raw audit final pU range",
						 diag->final_mean_pU >= -1e-6f && diag->final_mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("raw audit final sep finite", isfinite(diag->final_mean_sep));
	failed |= check_true("raw audit final sum k finite", isfinite(diag->final_sum_wedge_k));
	failed |= check_true("raw audit final sum k nonnegative", diag->final_sum_wedge_k >= 0.0);
	for (i = 0; i < HK_BLIND_P9016_RAW_AUDIT_N_ITER; ++i) {
		const struct hk_blind_iter_diag *pre = &per_iter[i].pre_relax_diag;
		failed |= check_i32("raw audit posterior nonfinite", pre->n_posterior_nonfinite, 0);
		failed |= check_i32("raw audit posterior bad sum", pre->n_posterior_bad_sum, 0);
		failed |= check_i32("raw audit posterior out of range", pre->n_posterior_out_of_range, 0);
		failed |= check_i32("raw audit uncertainty nonfinite", pre->n_uncertainty_nonfinite, 0);
		failed |= check_i32("raw audit uncertainty out of range", pre->n_uncertainty_out_of_range, 0);
		failed |= check_i32("raw audit five-state bad sum", pre->n_five_state_bad_sum, 0);
		failed |= check_i32("raw audit wedge nonfinite", pre->n_wedge_nonfinite, 0);
		failed |= check_i32("raw audit wedge bad k", pre->n_wedge_bad_k, 0);
		failed |= check_i32("raw audit wedge bad d_scale", pre->n_wedge_bad_d_scale, 0);
		failed |= check_i32("raw audit sep force nonfinite", pre->sep_force_nonfinite, 0);
		failed |= check_i32("raw audit relax completed", per_iter[i].relax_diag.n_completed,
							HK_BLIND_P9016_RAW_AUDIT_RELAX_STEPS);
		failed |= check_i32("raw audit relax nonfinite step", per_iter[i].relax_diag.n_nonfinite_step, 0);
		failed |= check_i32("raw audit relax coord bad", per_iter[i].relax_diag.n_coord_nonfinite, 0);
		failed |= check_i32("raw audit backbone bad", per_iter[i].relax_diag.n_backbone_nonfinite_step, 0);
		failed |= check_i32("raw audit repulsion bad", per_iter[i].relax_diag.n_repulsion_nonfinite_step, 0);
		failed |= check_true("raw audit relax total finite",
							 isfinite(per_iter[i].relax_diag.final_total_energy));
		failed |= check_true("raw audit repulsion finite",
							 isfinite(per_iter[i].relax_diag.final_repulsion_energy));
	}
	return failed;
}

static int write_outputs(const char *posterior_path, const char *raw_path,
						 const char *coords_path, const char *diag_path,
						 const struct hk_blind_pair *raw, int32_t n_raw,
						 const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
						 const fvec3_t *coords, const struct hk_blind_iter_loop_diag *loop_diag)
{
	FILE *fp;
	int ret;

	fp = fopen(posterior_path, "w");
	if (fp == 0) return -1;
	ret = hk_blind_write_bpair_posterior_tsv(fp, bmap, set);
	if (fclose(fp) != 0) ret = -1;
	if (ret != 0) return -1;

	fp = fopen(raw_path, "w");
	if (fp == 0) return -1;
	ret = hk_blind_write_raw_contact_posterior_tsv(fp, raw, n_raw, bmap, set);
	if (fclose(fp) != 0) ret = -1;
	if (ret != 0) return -1;

	fp = fopen(coords_path, "w");
	if (fp == 0) return -1;
	ret = hk_blind_write_diploid_coords_tsv(fp, bmap, coords);
	if (fclose(fp) != 0) ret = -1;
	if (ret != 0) return -1;

	fp = fopen(diag_path, "w");
	if (fp == 0) return -1;
	ret = hk_blind_write_iter_loop_diag_tsv(fp, loop_diag);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int write_manifest(const char *manifest_path, const char *out_dir,
						  const char *posterior_path, const char *raw_path,
						  const char *coords_path, const char *diag_path,
						  const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
						  const struct hk_blind_iter_loop_diag *loop_diag)
{
	FILE *fp = fopen(manifest_path, "w");
	struct hk_blind_base_k_stats base_k_stats;
	const char *status;
	if (fp == 0) return -1;
	hk_blind_bpair_set_base_k_stats(set, &base_k_stats);
	status = (loop_diag->n_bad_iter == 0 &&
			  loop_diag->n_relax_nonfinite_iter == 0 &&
			  loop_diag->n_coord_nonfinite == 0)? "OK" : "WARN";
	if (fprintf(fp,
				"key\tvalue\n"
				"sample\tP9016\n"
				"runner_family\ttest_fixture\n"
				"runner_version\t2026-04-30\n"
				"default_profile\tp9016_raw_audit_fixture_v1\n"
				"input_path\t%s\n"
				"output_dir\t%s\n"
				"n_raw\t%d\n"
				"n_bpair\t%d\n"
				"n_beads\t%d\n"
				"resolution\t%d\n"
				"n_iter\t%d\n"
				"unit\t%.9g\n"
				"d_scale\t%.9g\n"
				"base_k_mode\t%s\n"
				"base_k_effective\t%.9g\n"
				"base_k_min\t%.9g\n"
				"base_k_mean\t%.9g\n"
				"base_k_max\t%.9g\n"
				"base_k_n_nonfinite\t%d\n"
				"legacy_base_k_unused\t%.9g\n"
				"init_mode\ttoy_split\n"
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
				"same_bin_filter_enabled\t%d\n"
				"n_raw_same_bin_excluded\t%lld\n"
				"n_bpair_same_bin_excluded\t%d\n"
				"raw_posterior_same_bin_policy\tuniform_unknown_rows\n"
				"min_sep_unit\t%.9g\n"
				"lambda_sep\t%.9g\n"
				"relax_step\t%.9g\n"
				"relax_steps\t%d\n"
				"temperature_start\t%.9g\n"
				"temperature_end\t%.9g\n"
				"rho_train_start\t%.9g\n"
				"rho_train_end\t%.9g\n"
				"init_eps_effective\t%.9g\n"
				"init_noise_scale_effective\t%.9g\n"
				"init_split_params_used\t1\n"
				"init_seed\t%llu\n"
				"enable_repulsion\t1\n"
				"repulsion_mode\t%d\n"
				"repulsion_blocking_mode\tcurrent_edge_blocking\n"
				"write_raw_posterior\t1\n"
				"output_bpair_posterior\t%s\n"
				"output_coords\t%s\n"
				"output_loop_diag\t%s\n"
				"output_raw_posterior\t%s\n"
				"final_mean_entropy\t%.9g\n"
				"final_mean_pU\t%.9g\n"
				"final_sum_wedge_k\t%.17g\n"
				"posterior_refreshed_after_final_relax\t%d\n"
				"posterior_refresh_temperature\t%.9g\n"
				"posterior_refresh_prior_mode\tuniform\n"
				"posterior_refresh_mean_kl\t%.17g\n"
				"posterior_refresh_top_state_switch_frac\t%.9g\n"
				"posterior_refresh_mean_pU_before\t%.9g\n"
				"posterior_refresh_mean_pU_after\t%.9g\n"
				"n_bad_iter\t%d\n"
				"n_relax_nonfinite_iter\t%d\n"
				"n_coord_nonfinite\t%d\n"
				"status\t%s\n",
				HK_BLIND_P9016_RAW_AUDIT_PATH, out_dir, set->n_raw, set->n_bpairs, bmap->n_beads,
				HK_BLIND_P9016_RAW_AUDIT_RESOLUTION, loop_diag->n_completed,
				HK_BLIND_P9016_RAW_AUDIT_UNIT, HK_BLIND_P9016_RAW_AUDIT_D_SCALE,
				base_k_stats.n > 0 && base_k_stats.min == base_k_stats.max? "uniform" : "variable",
				base_k_stats.mean, base_k_stats.min, base_k_stats.mean, base_k_stats.max,
				base_k_stats.n_nonfinite, HK_BLIND_P9016_RAW_AUDIT_LEGACY_BASE_K_UNUSED,
				set->same_bin_filter_enabled, (long long)set->n_raw_same_bin_excluded,
				set->n_bpair_same_bin_excluded,
				HK_BLIND_P9016_RAW_AUDIT_MIN_SEP_UNIT,
				HK_BLIND_P9016_RAW_AUDIT_LAMBDA_SEP, HK_BLIND_P9016_RAW_AUDIT_STEP,
				HK_BLIND_P9016_RAW_AUDIT_RELAX_STEPS,
				HK_BLIND_P9016_RAW_AUDIT_TEMPERATURE_START,
				HK_BLIND_P9016_RAW_AUDIT_TEMPERATURE_END,
				HK_BLIND_P9016_RAW_AUDIT_RHO_TRAIN_START,
				HK_BLIND_P9016_RAW_AUDIT_RHO_TRAIN_END,
				HK_BLIND_P9016_RAW_AUDIT_INIT_EPS,
				HK_BLIND_P9016_RAW_AUDIT_INIT_NOISE_SCALE,
				(unsigned long long)HK_BLIND_P9016_RAW_AUDIT_INIT_SEED,
				HK_BLIND_REPULSION_CELL,
				posterior_path, coords_path, diag_path, raw_path,
				loop_diag->final_mean_entropy, loop_diag->final_mean_pU,
				loop_diag->final_sum_wedge_k,
				loop_diag->posterior_refreshed_after_final_relax,
				loop_diag->posterior_refresh_temperature,
				loop_diag->posterior_refresh_mean_kl,
				loop_diag->posterior_refresh_top_state_switch_frac,
				loop_diag->posterior_refresh_mean_pU_before,
				loop_diag->posterior_refresh_mean_pU_after,
				loop_diag->n_bad_iter,
				loop_diag->n_relax_nonfinite_iter, loop_diag->n_coord_nonfinite,
				status) < 0) {
		fclose(fp);
		return -1;
	}
	if (fclose(fp) != 0) return -1;
	return 0;
}

static int validate_output_files(const char *posterior_path, const char *raw_path,
								 const char *coords_path, const char *diag_path,
								 const char *manifest_path, const struct hk_bmap *bmap,
								 const struct hk_blind_bpair_set *set)
{
	int failed = 0;
	int64_t rows = 0;

	failed |= check_true("raw audit bpair output exists", file_exists(posterior_path));
	failed |= check_true("raw audit raw output exists", file_exists(raw_path));
	failed |= check_true("raw audit coords output exists", file_exists(coords_path));
	failed |= check_true("raw audit diag output exists", file_exists(diag_path));
	failed |= check_true("raw audit manifest exists", file_exists(manifest_path));
	failed |= check_true("raw audit bpair output nonempty", file_size_or_negative(posterior_path) > 0);
	failed |= check_true("raw audit raw output nonempty", file_size_or_negative(raw_path) > 0);
	failed |= check_true("raw audit coords output nonempty", file_size_or_negative(coords_path) > 0);
	failed |= check_true("raw audit diag output nonempty", file_size_or_negative(diag_path) > 0);
	failed |= check_true("raw audit manifest nonempty", file_size_or_negative(manifest_path) > 0);
	if (failed) return failed;

	failed |= check_i32("raw audit bpair row counter", count_data_rows(posterior_path, &rows), 0);
	failed |= check_true("raw audit bpair rows fit", rows == set->n_bpairs);
	failed |= check_i32("raw audit raw row counter", count_data_rows(raw_path, &rows), 0);
	failed |= check_true("raw audit raw rows fit", rows == set->n_raw);
	failed |= check_i32("raw audit coord row counter", count_data_rows(coords_path, &rows), 0);
	failed |= check_true("raw audit coord rows fit", rows == 2LL * bmap->n_beads);
	failed |= check_i32("raw audit loop diag row counter", count_data_rows(diag_path, &rows), 0);
	failed |= check_true("raw audit diag one row", rows == 1);
	return failed;
}

static int run_auditor(const char *out_dir)
{
	char cmd[1024];
	int ret;
	snprintf(cmd, sizeof(cmd), "./audit_blind_p9016_full_cpu_output.bin \"%s\"", out_dir);
	ret = system(cmd);
	if (ret != 0) {
		fprintf(stderr, "raw audit auditor command failed with status %d\n", ret);
		return 1;
	}
	return 0;
}

int main(void)
{
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	struct hk_blind_bpair_set *set = 0;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	struct hk_blind_single_iter_diag *per_iter = 0;
	fvec3_t *haploid = 0, *diploid = 0, *diploid_before = 0;
	char out_dir[512] = {0};
	char posterior_path[768], raw_path[768], coords_path[768], diag_path[768], manifest_path[768];
	int32_t n_prefix, n_diploid;
	int32_t i;
	float max_force_l1 = 0.0f;
	int moved = 0;
	int failed = 0;

	if (!file_exists(HK_BLIND_P9016_RAW_AUDIT_PATH)) {
		fprintf(stderr, "SKIP: %s not found\n", HK_BLIND_P9016_RAW_AUDIT_PATH);
		return 0;
	}

	hk_verbose = 0;
	m = hk_map_read(HK_BLIND_P9016_RAW_AUDIT_PATH);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", HK_BLIND_P9016_RAW_AUDIT_PATH);
		return 1;
	}
	n_prefix = m->n_pairs < HK_BLIND_P9016_RAW_AUDIT_PREFIX_N?
		m->n_pairs : HK_BLIND_P9016_RAW_AUDIT_PREFIX_N;
	failed |= check_true("raw audit prefix positive", n_prefix > 0);
	failed |= check_true("raw audit prefix bounded", n_prefix <= HK_BLIND_P9016_RAW_AUDIT_PREFIX_N);
	if (failed) goto cleanup;

	failed |= check_i32("raw audit temp dir", make_temp_dir(out_dir, sizeof(out_dir)), 0);
	if (failed) goto cleanup;
	snprintf(posterior_path, sizeof(posterior_path), "%s/p9016_full.bpair_posterior.tsv", out_dir);
	snprintf(raw_path, sizeof(raw_path), "%s/p9016_full.raw_posterior.tsv", out_dir);
	snprintf(coords_path, sizeof(coords_path), "%s/p9016_full.coords.tsv", out_dir);
	snprintf(diag_path, sizeof(diag_path), "%s/p9016_full.loop_diag.tsv", out_dir);
	snprintf(manifest_path, sizeof(manifest_path), "%s/p9016_full.manifest.tsv", out_dir);

	bmap = hk_bmap_gen(m->d, n_prefix, m->pairs, HK_BLIND_P9016_RAW_AUDIT_RESOLUTION, 1);
	raw = (struct hk_blind_pair*)calloc((size_t)n_prefix, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate P9016 raw-audit bmap/raw data\n");
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n_prefix; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);

	set = hk_blind_bpair_set_build(bmap, n_prefix, raw);
	failed |= check_binned_set(bmap, set, HK_BLIND_P9016_RAW_AUDIT_PREFIX_N);
	failed |= check_true("raw audit bead count positive", bmap && bmap->n_beads > 0);
	if (failed) goto cleanup;

	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	haploid = (fvec3_t*)calloc((size_t)bmap->n_beads, sizeof(*haploid));
	diploid = (fvec3_t*)calloc((size_t)n_diploid, sizeof(*diploid));
	diploid_before = (fvec3_t*)calloc((size_t)n_diploid, sizeof(*diploid_before));
	per_iter = (struct hk_blind_single_iter_diag*)calloc(HK_BLIND_P9016_RAW_AUDIT_N_ITER, sizeof(*per_iter));
	if (haploid == 0 || diploid == 0 || diploid_before == 0 || per_iter == 0) {
		fprintf(stderr, "failed to allocate P9016 raw-audit coordinate/diagnostic data\n");
		failed = 1;
		goto cleanup;
	}

	set_haploid_scaffold(bmap, haploid);
	failed |= check_i32("raw audit init diploid",
						hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
																  diploid, HK_BLIND_P9016_RAW_AUDIT_INIT_EPS,
																  HK_BLIND_P9016_RAW_AUDIT_INIT_NOISE_SCALE,
																  HK_BLIND_P9016_RAW_AUDIT_INIT_SEED), 0);
	failed |= check_coords_finite(diploid, n_diploid);
	if (failed) goto cleanup;
	copy_coords(diploid_before, diploid, n_diploid);

	hk_fdg_conf_init(&fdg_conf);
	set_schedule_conf(&schedule_conf);
	failed |= check_i32("raw audit scheduled ret",
						hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &fdg_conf, diploid, 0,
															 &schedule_conf, per_iter, &loop_diag), 0);
	failed |= check_loop_health(&loop_diag, per_iter);
	for (i = 0; i < HK_BLIND_P9016_RAW_AUDIT_N_ITER; ++i)
		if (per_iter[i].relax_diag.max_force_l1 > max_force_l1)
			max_force_l1 = per_iter[i].relax_diag.max_force_l1;
	failed |= check_coords_finite(diploid, n_diploid);
	moved = coords_any_changed(diploid_before, diploid, n_diploid, 1e-8f);
	if (max_force_l1 > 1e-12f)
		failed |= check_true("raw audit coordinate movement", moved);
	else
		fprintf(stderr, "NOTE: P9016 raw-output audit force_l1 is zero; no movement expected\n");
	if (failed) goto cleanup;

	failed |= check_i32("raw audit write outputs",
						write_outputs(posterior_path, raw_path, coords_path, diag_path,
									  raw, n_prefix, bmap, set, diploid, &loop_diag), 0);
	failed |= check_i32("raw audit write manifest",
						write_manifest(manifest_path, out_dir, posterior_path, raw_path,
									   coords_path, diag_path, bmap, set, &loop_diag), 0);
	failed |= validate_output_files(posterior_path, raw_path, coords_path, diag_path,
									manifest_path, bmap, set);
	if (failed) goto cleanup;

	failed |= run_auditor(out_dir);
	if (failed) goto cleanup;

	fprintf(stderr,
			"P9016 raw-output audit smoke: n_raw=%d n_bpair=%d n_beads=%d "
			"raw_rows=%d final_entropy=%.8g final_pU=%.8g "
			"final_sum_wedge_k=%.8g moved=%d auditor=OK\n",
			set->n_raw, set->n_bpairs, bmap->n_beads, set->n_raw,
			loop_diag.final_mean_entropy, loop_diag.final_mean_pU,
			loop_diag.final_sum_wedge_k, moved);

cleanup:
	remove_outputs(out_dir);
	free(haploid);
	free(diploid);
	free(diploid_before);
	free(per_iter);
	hk_blind_bpair_set_destroy(set);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
