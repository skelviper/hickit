#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include "hickit.h"

#define HK_BLIND_P9016_CAND_PATH "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_CAND_RESOLUTION 1000000
#define HK_BLIND_P9016_CAND_N_ITER 3
#define HK_BLIND_P9016_CAND_RELAX_STEPS 50
#define HK_BLIND_P9016_CAND_RELAX_STEP 0.0005f
#define HK_BLIND_P9016_CAND_TEMPERATURE_START 2.0f
#define HK_BLIND_P9016_CAND_TEMPERATURE_END 1.0f
#define HK_BLIND_P9016_CAND_RHO_TRAIN_START 1.0f
#define HK_BLIND_P9016_CAND_RHO_TRAIN_END 1.0f
#define HK_BLIND_P9016_CAND_UNIT 1.0f
#define HK_BLIND_P9016_CAND_D_SCALE 1.0f
#define HK_BLIND_P9016_CAND_BASE_K 2.0f
#define HK_BLIND_P9016_CAND_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_CAND_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_CAND_INIT_EPS 0.5f
#define HK_BLIND_P9016_CAND_INIT_NOISE_SCALE 0.0f
#define HK_BLIND_P9016_CAND_INIT_SEED 17ULL
#define HK_BLIND_P9016_CAND_PREFIX_SCAN 8192

struct candidate_config {
	const char *name;
	float multiplier;
};

struct bpair_metrics {
	float mean_pU, median_pU, frac_pU_lt_0_75;
	float mean_pmax, median_pmax;
	float mean_margin, median_margin;
	float mean_psame, mean_pcross;
};

struct sep_metrics {
	float mean, median, min, max;
};

struct candidate_result {
	const struct candidate_config *config;
	char output_dir[768];
	int32_t n_raw;
	int32_t n_bpair;
	int32_t n_beads;
	int same_bin_filter_enabled;
	int64_t n_raw_same_bin_excluded;
	int32_t n_bpair_same_bin_excluded;
	float k_rel_rep;
	struct hk_blind_iter_loop_diag loop_diag;
	struct bpair_metrics bpair;
	struct sep_metrics sep;
	float final_contact_energy;
	float final_backbone_energy;
	float final_repulsion_energy;
	float final_sep_energy;
	float final_total_energy;
	float final_force_l1;
	int audit_ok;
	int status_ok;
	double elapsed_cpu_sec;
};

static const struct candidate_config candidate_configs[] = {
	{ "rep10_rs50", 10.0f },
	{ "rep12_rs50", 12.0f }
};

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

static char *read_file_prefix(const char *path, size_t max_len)
{
	FILE *fp = fopen(path, "rb");
	char *buf;
	size_t n_read;

	if (fp == 0) return 0;
	buf = (char*)malloc(max_len + 1);
	if (buf == 0) {
		fclose(fp);
		return 0;
	}
	n_read = fread(buf, 1, max_len, fp);
	buf[n_read] = 0;
	fclose(fp);
	return buf;
}

static int make_output_root(char *dir, size_t dir_size)
{
	long pid = (long)getpid();
	long stamp = (long)time(0);
	int i;

	for (i = 0; i < 100; ++i) {
		snprintf(dir, dir_size, "/tmp/hk_blind_p9016_candidate_outputs_%ld_%ld_%d",
				 pid, stamp, i);
		if (mkdir(dir, 0700) == 0)
			return 0;
		if (errno != EEXIST)
			return -1;
	}
	return -1;
}

static void path_join(char *dst, size_t dst_size, const char *dir, const char *name)
{
	size_t n = strlen(dir);
	if (n > 0 && dir[n - 1] == '/')
		snprintf(dst, dst_size, "%s%s", dir, name);
	else
		snprintf(dst, dst_size, "%s/%s", dir, name);
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

static float dist3(const fvec3_t a, const fvec3_t b)
{
	float dx = a[0] - b[0];
	float dy = a[1] - b[1];
	float dz = a[2] - b[2];
	return sqrtf(dx * dx + dy * dy + dz * dz);
}

static int float_cmp(const void *a_, const void *b_)
{
	const float a = *(const float*)a_;
	const float b = *(const float*)b_;
	return (a > b) - (a < b);
}

static float median_float(float *a, int32_t n)
{
	assert(a);
	assert(n > 0);
	qsort(a, (size_t)n, sizeof(*a), float_cmp);
	if ((n & 1) != 0)
		return a[n >> 1];
	return 0.5f * (a[(n >> 1) - 1] + a[n >> 1]);
}

static int check_coords_finite(const fvec3_t *coords, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("candidate coord finite", isfinite(coords[i][a]));
	return failed;
}

static void set_schedule_conf(struct hk_blind_iter_schedule_conf *conf)
{
	conf->n_iter = HK_BLIND_P9016_CAND_N_ITER;
	conf->base_conf.unit = HK_BLIND_P9016_CAND_UNIT;
	conf->base_conf.d_scale = HK_BLIND_P9016_CAND_D_SCALE;
	conf->base_conf.base_k = HK_BLIND_P9016_CAND_BASE_K;
	conf->base_conf.temperature = HK_BLIND_P9016_CAND_TEMPERATURE_START;
	conf->base_conf.rho_train = HK_BLIND_P9016_CAND_RHO_TRAIN_START;
	conf->base_conf.min_sep_unit = HK_BLIND_P9016_CAND_MIN_SEP_UNIT;
	conf->base_conf.lambda_sep = HK_BLIND_P9016_CAND_LAMBDA_SEP;
	conf->base_conf.relax_step = HK_BLIND_P9016_CAND_RELAX_STEP;
	conf->base_conf.relax_steps = HK_BLIND_P9016_CAND_RELAX_STEPS;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->base_conf.d_scale_eps_count = 1e-6f;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = HK_BLIND_P9016_CAND_TEMPERATURE_START;
	conf->temperature_end = HK_BLIND_P9016_CAND_TEMPERATURE_END;
	conf->rho_train_start = HK_BLIND_P9016_CAND_RHO_TRAIN_START;
	conf->rho_train_end = HK_BLIND_P9016_CAND_RHO_TRAIN_END;
}

static int check_binned_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
							int32_t n_raw)
{
	int failed = 0;
	int32_t i;

	failed |= check_true("candidate set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_i32("candidate raw count", set->n_raw, n_raw);
	failed |= check_true("candidate bpair positive", set->n_bpairs > 0);
	for (i = 0; i < set->n_raw; ++i) {
		failed |= check_true("candidate raw2binned id lower", set->raw2binned[i].bpair_id >= 0);
		failed |= check_true("candidate raw2binned id upper", set->raw2binned[i].bpair_id < set->n_bpairs);
		failed |= check_true("candidate raw2binned swapped",
							 set->raw2binned[i].swapped == 0 || set->raw2binned[i].swapped == 1);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		failed |= check_true("candidate key sorted", bp->key.bid[0] <= bp->key.bid[1]);
		failed |= check_true("candidate key bid0 valid", bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		failed |= check_true("candidate key bid1 valid", bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
	}
	return failed;
}

static int check_loop_health(const struct hk_blind_iter_loop_diag *diag,
							 const struct hk_blind_single_iter_diag *per_iter)
{
	int failed = 0;
	int32_t i;

	failed |= check_i32("candidate loop n_iter", diag->n_iter, HK_BLIND_P9016_CAND_N_ITER);
	failed |= check_i32("candidate loop completed", diag->n_completed, HK_BLIND_P9016_CAND_N_ITER);
	failed |= check_i32("candidate loop bad iter", diag->n_bad_iter, 0);
	failed |= check_i32("candidate loop relax bad", diag->n_relax_nonfinite_iter, 0);
	failed |= check_i32("candidate loop coord bad", diag->n_coord_nonfinite, 0);
	failed |= check_close("candidate initial temperature", diag->initial_temperature,
						  HK_BLIND_P9016_CAND_TEMPERATURE_START);
	failed |= check_close("candidate final temperature", diag->final_temperature,
						  HK_BLIND_P9016_CAND_TEMPERATURE_END);
	failed |= check_close("candidate initial rho", diag->initial_rho_train,
						  HK_BLIND_P9016_CAND_RHO_TRAIN_START);
	failed |= check_close("candidate final rho", diag->final_rho_train,
						  HK_BLIND_P9016_CAND_RHO_TRAIN_END);
	failed |= check_true("candidate final entropy finite", isfinite(diag->final_mean_entropy));
	failed |= check_true("candidate final pU finite", isfinite(diag->final_mean_pU));
	failed |= check_true("candidate final pU range",
						 diag->final_mean_pU >= -1e-6f && diag->final_mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("candidate final sep finite", isfinite(diag->final_mean_sep));
	failed |= check_true("candidate final sum k finite", isfinite(diag->final_sum_wedge_k));
	failed |= check_true("candidate final sum k nonnegative", diag->final_sum_wedge_k >= 0.0);
	for (i = 0; i < HK_BLIND_P9016_CAND_N_ITER; ++i) {
		const struct hk_blind_iter_diag *pre = &per_iter[i].pre_relax_diag;
		failed |= check_i32("candidate posterior nonfinite", pre->n_posterior_nonfinite, 0);
		failed |= check_i32("candidate posterior bad sum", pre->n_posterior_bad_sum, 0);
		failed |= check_i32("candidate posterior out of range", pre->n_posterior_out_of_range, 0);
		failed |= check_i32("candidate uncertainty nonfinite", pre->n_uncertainty_nonfinite, 0);
		failed |= check_i32("candidate uncertainty out of range", pre->n_uncertainty_out_of_range, 0);
		failed |= check_i32("candidate five-state bad sum", pre->n_five_state_bad_sum, 0);
		failed |= check_i32("candidate wedge nonfinite", pre->n_wedge_nonfinite, 0);
		failed |= check_i32("candidate wedge bad k", pre->n_wedge_bad_k, 0);
		failed |= check_i32("candidate wedge bad d_scale", pre->n_wedge_bad_d_scale, 0);
		failed |= check_i32("candidate sep force nonfinite", pre->sep_force_nonfinite, 0);
		failed |= check_i32("candidate relax completed", per_iter[i].relax_diag.n_completed,
							HK_BLIND_P9016_CAND_RELAX_STEPS);
		failed |= check_i32("candidate relax nonfinite step", per_iter[i].relax_diag.n_nonfinite_step, 0);
		failed |= check_i32("candidate relax coord bad", per_iter[i].relax_diag.n_coord_nonfinite, 0);
		failed |= check_i32("candidate backbone bad", per_iter[i].relax_diag.n_backbone_nonfinite_step, 0);
		failed |= check_i32("candidate repulsion bad", per_iter[i].relax_diag.n_repulsion_nonfinite_step, 0);
		failed |= check_true("candidate relax total finite", isfinite(per_iter[i].relax_diag.final_total_energy));
		failed |= check_true("candidate repulsion finite", isfinite(per_iter[i].relax_diag.final_repulsion_energy));
	}
	return failed;
}

static int validate_final_bpair_set(const struct hk_blind_bpair_set *set)
{
	struct hk_blind_iter_diag diag;
	int failed = 0;

	hk_blind_iter_diag_init(&diag);
	hk_blind_iter_diag_validate_bpair_set(set, &diag);
	failed |= check_i32("candidate final posterior nonfinite", diag.n_posterior_nonfinite, 0);
	failed |= check_i32("candidate final posterior bad sum", diag.n_posterior_bad_sum, 0);
	failed |= check_i32("candidate final posterior out of range", diag.n_posterior_out_of_range, 0);
	failed |= check_i32("candidate final uncertainty nonfinite", diag.n_uncertainty_nonfinite, 0);
	failed |= check_i32("candidate final uncertainty out of range", diag.n_uncertainty_out_of_range, 0);
	failed |= check_i32("candidate final five-state bad sum", diag.n_five_state_bad_sum, 0);
	return failed;
}

static int compute_bpair_metrics(const struct hk_blind_bpair_set *set, struct bpair_metrics *out)
{
	float *pU = 0, *pmax = 0, *margin = 0;
	double sum_pU = 0.0, sum_pmax = 0.0, sum_margin = 0.0;
	double sum_psame = 0.0, sum_pcross = 0.0;
	int32_t n, i;
	int failed = 0;

	assert(set);
	assert(out);
	memset(out, 0, sizeof(*out));
	n = set->n_bpairs;
	if (n <= 0)
		return 1;
	pU = (float*)malloc((size_t)n * sizeof(*pU));
	pmax = (float*)malloc((size_t)n * sizeof(*pmax));
	margin = (float*)malloc((size_t)n * sizeof(*margin));
	if (pU == 0 || pmax == 0 || margin == 0) {
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		float psame = bp->p4[HK_BLIND_STATE_00] + bp->p4[HK_BLIND_STATE_11];
		float pcross = bp->p4[HK_BLIND_STATE_01] + bp->p4[HK_BLIND_STATE_10];
		pU[i] = bp->pU;
		pmax[i] = bp->pmax;
		margin[i] = bp->margin;
		failed |= check_true("candidate metric pU finite", isfinite(pU[i]));
		failed |= check_true("candidate metric pmax finite", isfinite(pmax[i]));
		failed |= check_true("candidate metric margin finite", isfinite(margin[i]));
		failed |= check_true("candidate metric psame finite", isfinite(psame));
		failed |= check_true("candidate metric pcross finite", isfinite(pcross));
		if (pU[i] < 0.75f) out->frac_pU_lt_0_75 += 1.0f;
		sum_pU += pU[i];
		sum_pmax += pmax[i];
		sum_margin += margin[i];
		sum_psame += psame;
		sum_pcross += pcross;
	}
	if (failed) goto cleanup;
	out->mean_pU = (float)(sum_pU / n);
	out->mean_pmax = (float)(sum_pmax / n);
	out->mean_margin = (float)(sum_margin / n);
	out->mean_psame = (float)(sum_psame / n);
	out->mean_pcross = (float)(sum_pcross / n);
	out->median_pU = median_float(pU, n);
	out->median_pmax = median_float(pmax, n);
	out->median_margin = median_float(margin, n);
	out->frac_pU_lt_0_75 /= (float)n;
	failed |= check_true("candidate mean pU range", out->mean_pU >= -1e-6f && out->mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("candidate mean pmax range",
						 out->mean_pmax >= 0.25f - 1e-6f && out->mean_pmax <= 1.0f + 1e-6f);
	failed |= check_true("candidate mean margin range",
						 out->mean_margin >= -1e-6f && out->mean_margin <= 1.0f + 1e-6f);
	failed |= check_true("candidate mean psame range",
						 out->mean_psame >= -1e-6f && out->mean_psame <= 1.0f + 1e-6f);
	failed |= check_true("candidate mean pcross range",
						 out->mean_pcross >= -1e-6f && out->mean_pcross <= 1.0f + 1e-6f);

cleanup:
	free(pU);
	free(pmax);
	free(margin);
	return failed;
}

static int compute_sep_metrics(const struct hk_bmap *bmap, const fvec3_t *coords,
							   struct sep_metrics *out)
{
	float *sep = 0;
	double sum = 0.0;
	int32_t n, i;
	int failed = 0;

	assert(bmap);
	assert(out);
	memset(out, 0, sizeof(*out));
	n = bmap->n_beads;
	if (n <= 0)
		return 1;
	sep = (float*)malloc((size_t)n * sizeof(*sep));
	if (sep == 0)
		return 1;
	out->min = INFINITY;
	for (i = 0; i < n; ++i) {
		int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		sep[i] = dist3(coords[b0], coords[b1]);
		failed |= check_true("candidate sep finite", isfinite(sep[i]));
		if (sep[i] < out->min) out->min = sep[i];
		if (sep[i] > out->max) out->max = sep[i];
		sum += sep[i];
	}
	if (failed) goto cleanup;
	out->mean = (float)(sum / n);
	out->median = median_float(sep, n);
	failed |= check_true("candidate sep mean finite", isfinite(out->mean));
	failed |= check_true("candidate sep median finite", isfinite(out->median));
	failed |= check_true("candidate sep min finite", isfinite(out->min));
	failed |= check_true("candidate sep max finite", isfinite(out->max));

cleanup:
	free(sep);
	return failed;
}

static int validate_result_metrics(const struct candidate_result *result)
{
	int failed = 0;
	failed |= check_true("candidate final contact energy finite", isfinite(result->final_contact_energy));
	failed |= check_true("candidate final backbone energy finite", isfinite(result->final_backbone_energy));
	failed |= check_true("candidate final repulsion energy finite", isfinite(result->final_repulsion_energy));
	failed |= check_true("candidate final sep energy finite", isfinite(result->final_sep_energy));
	failed |= check_true("candidate final total energy finite", isfinite(result->final_total_energy));
	failed |= check_true("candidate final force finite", isfinite(result->final_force_l1));
	failed |= check_true("candidate mean pU finite", isfinite(result->bpair.mean_pU));
	failed |= check_true("candidate median pU finite", isfinite(result->bpair.median_pU));
	failed |= check_true("candidate mean pmax finite", isfinite(result->bpair.mean_pmax));
	failed |= check_true("candidate median pmax finite", isfinite(result->bpair.median_pmax));
	failed |= check_true("candidate mean margin finite", isfinite(result->bpair.mean_margin));
	failed |= check_true("candidate median margin finite", isfinite(result->bpair.median_margin));
	failed |= check_true("candidate sep mean finite result", isfinite(result->sep.mean));
	failed |= check_true("candidate sep median finite result", isfinite(result->sep.median));
	return failed;
}

static int write_required_outputs(const char *posterior_path, const char *coords_path,
								  const char *diag_path, const struct hk_bmap *bmap,
								  const struct hk_blind_bpair_set *set, const fvec3_t *coords,
								  const struct hk_blind_iter_loop_diag *loop_diag)
{
	FILE *fp;
	int ret;

	fp = fopen(posterior_path, "w");
	if (fp == 0) return -1;
	ret = hk_blind_write_bpair_posterior_tsv(fp, bmap, set);
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

static int write_optional_raw_output(const char *raw_path, const struct hk_blind_pair *raw,
									 int32_t n_raw, const struct hk_bmap *bmap,
									 const struct hk_blind_bpair_set *set)
{
	FILE *fp;
	int ret;

	fp = fopen(raw_path, "w");
	if (fp == 0) return -1;
	ret = hk_blind_write_raw_contact_posterior_tsv(fp, raw, n_raw, bmap, set);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int write_manifest(const char *manifest_path, const char *out_dir,
						  const char *posterior_path, const char *coords_path,
						  const char *diag_path, const char *raw_path,
						  int write_raw, const struct candidate_result *result)
{
	const struct candidate_config *config = result->config;
	const struct hk_blind_iter_loop_diag *loop_diag = &result->loop_diag;
	FILE *fp = fopen(manifest_path, "w");
	const char *status;
	if (fp == 0) return -1;
	status = (loop_diag->n_bad_iter == 0 &&
			  loop_diag->n_relax_nonfinite_iter == 0 &&
			  loop_diag->n_coord_nonfinite == 0)? "OK" : "WARN";
	if (fprintf(fp,
				"key\tvalue\n"
				"sample\tP9016\n"
				"runner_family\tcandidate_outputs\n"
				"runner_version\t2026-04-30\n"
				"default_profile\tp9016_candidate_outputs_v1\n"
				"input_path\t%s\n"
				"output_dir\t%s\n"
				"candidate_name\t%s\n"
				"k_rel_rep_multiplier\t%.9g\n"
				"k_rel_rep\t%.9g\n"
				"n_raw\t%d\n"
				"n_bpair\t%d\n"
				"n_beads\t%d\n"
				"resolution\t%d\n"
				"n_iter\t%d\n"
				"unit\t%.9g\n"
				"d_scale\t%.9g\n"
				"base_k_mode\tuniform\n"
				"base_k_effective\t1\n"
				"base_k_min\t1\n"
				"base_k_mean\t1\n"
				"base_k_max\t1\n"
				"base_k_n_nonfinite\t0\n"
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
				"rho_train_mode\t%s\n"
				"d_scale_mode\t%s\n"
				"d_scale_eps_count\t1e-06\n"
				"same_bin_filter_enabled\t%d\n"
				"n_raw_same_bin_excluded\t%lld\n"
				"n_bpair_same_bin_excluded\t%d\n"
				"raw_posterior_same_bin_policy\tuniform_unknown_rows\n"
				"min_sep_unit\t%.9g\n"
				"lambda_sep\t%.9g\n"
				"relax_step\t%.9g\n"
				"relax_steps\t%d\n"
				"enable_repulsion\t%d\n"
				"repulsion_mode\t%d\n"
				"temperature_start\t%.9g\n"
				"temperature_end\t%.9g\n"
				"rho_train_start\t%.9g\n"
				"rho_train_end\t%.9g\n"
				"init_eps_effective\t%.9g\n"
				"init_noise_scale_effective\t%.9g\n"
				"init_split_params_used\t1\n"
				"init_seed\t%llu\n"
				"repulsion_blocking_mode\tcurrent_edge_blocking\n"
				"write_raw_posterior\t%d\n"
				"output_bpair_posterior\t%s\n"
				"output_coords\t%s\n"
				"output_loop_diag\t%s\n",
				HK_BLIND_P9016_CAND_PATH, out_dir, config->name,
				config->multiplier, result->k_rel_rep,
				result->n_raw, result->n_bpair, result->n_beads,
				HK_BLIND_P9016_CAND_RESOLUTION, HK_BLIND_P9016_CAND_N_ITER,
				HK_BLIND_P9016_CAND_UNIT, HK_BLIND_P9016_CAND_D_SCALE,
				HK_BLIND_P9016_CAND_BASE_K,
				hk_blind_rho_train_mode_name(HK_BLIND_RHO_TRAIN_CONSTANT),
				hk_blind_d_scale_mode_name(HK_BLIND_D_SCALE_RAW_COUNT),
				result->same_bin_filter_enabled,
				(long long)result->n_raw_same_bin_excluded,
				result->n_bpair_same_bin_excluded,
				HK_BLIND_P9016_CAND_MIN_SEP_UNIT,
				HK_BLIND_P9016_CAND_LAMBDA_SEP, HK_BLIND_P9016_CAND_RELAX_STEP,
				HK_BLIND_P9016_CAND_RELAX_STEPS, 1, HK_BLIND_REPULSION_CELL,
				HK_BLIND_P9016_CAND_TEMPERATURE_START, HK_BLIND_P9016_CAND_TEMPERATURE_END,
				HK_BLIND_P9016_CAND_RHO_TRAIN_START, HK_BLIND_P9016_CAND_RHO_TRAIN_END,
				HK_BLIND_P9016_CAND_INIT_EPS, HK_BLIND_P9016_CAND_INIT_NOISE_SCALE,
				(unsigned long long)HK_BLIND_P9016_CAND_INIT_SEED, write_raw,
				posterior_path, coords_path, diag_path) < 0) {
		fclose(fp);
		return -1;
	}
	if (write_raw && fprintf(fp, "output_raw_posterior\t%s\n", raw_path) < 0) {
		fclose(fp);
		return -1;
	}
	if (fprintf(fp,
				"final_mean_entropy\t%.9g\n"
				"final_mean_pU\t%.9g\n"
				"final_sum_wedge_k\t%.17g\n"
				"final_repulsion_energy\t%.9g\n"
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
				loop_diag->final_mean_entropy, loop_diag->final_mean_pU,
				loop_diag->final_sum_wedge_k, loop_diag->final_repulsion_energy,
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

static int validate_file_prefix(const char *label, const char *path, const char *header_token)
{
	char *prefix;
	int failed = 0;

	failed |= check_true(label, file_exists(path));
	failed |= check_true("candidate output size positive", file_size_or_negative(path) > 0);
	prefix = read_file_prefix(path, HK_BLIND_P9016_CAND_PREFIX_SCAN);
	failed |= check_true("candidate output prefix readable", prefix != 0);
	if (prefix == 0) return 1;
	failed |= check_true("candidate output header token", strstr(prefix, header_token) != 0);
	free(prefix);
	return failed;
}

static int validate_manifest_prefix(const char *path, int write_raw, const char *candidate_name)
{
	char *prefix = read_file_prefix(path, HK_BLIND_P9016_CAND_PREFIX_SCAN);
	int failed = 0;
	char expected_candidate[128];

	snprintf(expected_candidate, sizeof(expected_candidate), "candidate_name\t%s", candidate_name);
	failed |= check_true("candidate manifest exists", file_exists(path));
	failed |= check_true("candidate manifest size positive", file_size_or_negative(path) > 0);
	failed |= check_true("candidate manifest readable", prefix != 0);
	if (prefix == 0) return 1;
	failed |= check_true("candidate manifest sample", strstr(prefix, "sample\tP9016") != 0);
	failed |= check_true("candidate manifest candidate", strstr(prefix, expected_candidate) != 0);
	failed |= check_true("candidate manifest n_raw", strstr(prefix, "n_raw") != 0);
	failed |= check_true("candidate manifest n_bpair", strstr(prefix, "n_bpair") != 0);
	failed |= check_true("candidate manifest n_beads", strstr(prefix, "n_beads") != 0);
	failed |= check_true("candidate manifest status OK", strstr(prefix, "status\tOK") != 0);
	failed |= check_true("candidate manifest raw flag",
						 write_raw? strstr(prefix, "write_raw_posterior\t1") != 0
								  : strstr(prefix, "write_raw_posterior\t0") != 0);
	free(prefix);
	return failed;
}

static int run_auditor(const char *out_dir)
{
	char cmd[1024];
	int ret;
	snprintf(cmd, sizeof(cmd), "./audit_blind_p9016_full_cpu_output.bin \"%s\"", out_dir);
	ret = system(cmd);
	if (ret != 0) {
		fprintf(stderr, "candidate audit failed for %s with status %d\n", out_dir, ret);
		return 1;
	}
	return 0;
}

static int write_summary_header(FILE *fp)
{
	return fprintf(fp,
				   "candidate_name\toutput_dir\tmultiplier\tk_rel_rep\tn_raw\tn_bpair\tn_beads\t"
				   "n_iter\trelax_steps\trelax_step\tfinal_mean_entropy\tfinal_mean_pU\t"
				   "mean_pU\tmedian_pU\tfrac_pU_lt_0_75\tmean_pmax\tmedian_pmax\t"
				   "mean_margin\tmedian_margin\tmean_psame\tmean_pcross\tsep_mean\tsep_median\t"
				   "sep_min\tsep_max\tfinal_contact_energy\tfinal_backbone_energy\t"
				   "final_repulsion_energy\tfinal_sep_energy\tfinal_total_energy\t"
				   "repulsion_over_contact\tcontact_energy_per_wedge_k\tn_bad_iter\t"
				   "n_relax_nonfinite_iter\tn_coord_nonfinite\taudit_status\tstatus\t"
				   "elapsed_cpu_sec\n") < 0? -1 : 0;
}

static int append_summary_row(FILE *fp, const struct candidate_result *result)
{
	const struct candidate_config *config = result->config;
	const struct hk_blind_iter_loop_diag *d = &result->loop_diag;
	const double repulsion_over_contact = result->final_contact_energy > 0.0f?
		(double)result->final_repulsion_energy / (double)result->final_contact_energy : 0.0;
	const double contact_energy_per_wedge_k = d->final_sum_wedge_k > 0.0?
		(double)result->final_contact_energy / d->final_sum_wedge_k : 0.0;
	const char *audit_status = result->audit_ok? "OK" : "FAIL";
	const char *status = result->status_ok? "OK" : "FAIL";

	return fprintf(fp,
				   "%s\t%s\t%.9g\t%.9g\t%d\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%d\t%d\t%d\t%s\t%s\t%.3f\n",
				   config->name, result->output_dir, config->multiplier, result->k_rel_rep,
				   result->n_raw, result->n_bpair, result->n_beads,
				   HK_BLIND_P9016_CAND_N_ITER, HK_BLIND_P9016_CAND_RELAX_STEPS,
				   HK_BLIND_P9016_CAND_RELAX_STEP, d->final_mean_entropy, d->final_mean_pU,
				   result->bpair.mean_pU, result->bpair.median_pU,
				   result->bpair.frac_pU_lt_0_75, result->bpair.mean_pmax,
				   result->bpair.median_pmax, result->bpair.mean_margin,
				   result->bpair.median_margin, result->bpair.mean_psame,
				   result->bpair.mean_pcross, result->sep.mean, result->sep.median,
				   result->sep.min, result->sep.max, result->final_contact_energy,
				   result->final_backbone_energy, result->final_repulsion_energy,
				   result->final_sep_energy, result->final_total_energy,
				   repulsion_over_contact, contact_energy_per_wedge_k,
				   d->n_bad_iter, d->n_relax_nonfinite_iter, d->n_coord_nonfinite,
				   audit_status, status, result->elapsed_cpu_sec) < 0? -1 : 0;
}

static int validate_summary_file(const char *summary_path, int n_expected)
{
	FILE *fp = fopen(summary_path, "r");
	char line[8192];
	int n_rows = -1;
	int n_fail = 0;
	int failed = 0;

	if (fp == 0) return 1;
	while (fgets(line, sizeof(line), fp)) {
		if (n_rows < 0) {
			failed |= check_true("candidate summary header name", strstr(line, "candidate_name") != 0);
			failed |= check_true("candidate summary header audit", strstr(line, "audit_status") != 0);
			failed |= check_true("candidate summary header ratio", strstr(line, "repulsion_over_contact") != 0);
		} else if (strstr(line, "\tFAIL\t") != 0 || strstr(line, "\tFAIL\n") != 0) {
			++n_fail;
		}
		++n_rows;
	}
	if (ferror(fp)) failed = 1;
	fclose(fp);
	failed |= check_i32("candidate summary row count", n_rows, n_expected);
	failed |= check_i32("candidate summary fail count", n_fail, 0);
	return failed;
}

static int run_one_candidate(const struct candidate_config *config, const struct hk_bmap *bmap,
							 const struct hk_blind_pair *raw, int32_t n_raw,
							 const fvec3_t *haploid, const char *root_dir,
							 int write_raw, FILE *summary_fp,
							 struct candidate_result *result)
{
	struct hk_blind_bpair_set *set = 0;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_single_iter_diag *per_iter = 0;
	fvec3_t *diploid = 0;
	char posterior_path[1024], coords_path[1024], diag_path[1024], raw_path[1024], manifest_path[1024];
	int32_t n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	int failed = 0;
	clock_t t0 = clock(), t1;

	memset(result, 0, sizeof(*result));
	result->config = config;
	path_join(result->output_dir, sizeof(result->output_dir), root_dir, config->name);
	if (mkdir(result->output_dir, 0700) != 0) {
		fprintf(stderr, "failed to create candidate directory %s\n", result->output_dir);
		failed = 1;
		goto cleanup;
	}
	path_join(posterior_path, sizeof(posterior_path), result->output_dir, "p9016_full.bpair_posterior.tsv");
	path_join(coords_path, sizeof(coords_path), result->output_dir, "p9016_full.coords.tsv");
	path_join(diag_path, sizeof(diag_path), result->output_dir, "p9016_full.loop_diag.tsv");
	path_join(raw_path, sizeof(raw_path), result->output_dir, "p9016_full.raw_posterior.tsv");
	path_join(manifest_path, sizeof(manifest_path), result->output_dir, "p9016_full.manifest.tsv");

	fprintf(stderr, "P9016 candidate %s: multiplier=%.4g relax_steps=%d\n",
			config->name, config->multiplier, HK_BLIND_P9016_CAND_RELAX_STEPS);
	set = hk_blind_bpair_set_build(bmap, n_raw, raw);
	failed |= check_binned_set(bmap, set, n_raw);
	if (set) {
		result->n_raw = set->n_raw;
		result->n_bpair = set->n_bpairs;
		result->same_bin_filter_enabled = set->same_bin_filter_enabled;
		result->n_raw_same_bin_excluded = set->n_raw_same_bin_excluded;
		result->n_bpair_same_bin_excluded = set->n_bpair_same_bin_excluded;
	}
	result->n_beads = bmap->n_beads;
	if (failed) goto cleanup;

	diploid = (fvec3_t*)calloc((size_t)n_diploid, sizeof(*diploid));
	per_iter = (struct hk_blind_single_iter_diag*)calloc((size_t)HK_BLIND_P9016_CAND_N_ITER, sizeof(*per_iter));
	if (diploid == 0 || per_iter == 0) {
		fprintf(stderr, "failed to allocate candidate coordinate/diagnostic data for %s\n",
				config->name);
		failed = 1;
		goto cleanup;
	}
	failed |= check_i32("candidate init diploid",
						hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
																  diploid, HK_BLIND_P9016_CAND_INIT_EPS,
																  HK_BLIND_P9016_CAND_INIT_NOISE_SCALE,
																  HK_BLIND_P9016_CAND_INIT_SEED), 0);
	failed |= check_coords_finite(diploid, n_diploid);
	if (failed) goto cleanup;

	hk_fdg_conf_init(&fdg_conf);
	fdg_conf.k_rel_rep *= config->multiplier;
	result->k_rel_rep = fdg_conf.k_rel_rep;
	set_schedule_conf(&schedule_conf);
	failed |= check_i32("candidate scheduled ret",
						hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &fdg_conf, diploid, 0,
															 &schedule_conf, per_iter,
															 &result->loop_diag), 0);
	failed |= check_loop_health(&result->loop_diag, per_iter);
	failed |= check_coords_finite(diploid, n_diploid);
	failed |= validate_final_bpair_set(set);
	if (failed) goto cleanup;

	result->final_contact_energy = per_iter[HK_BLIND_P9016_CAND_N_ITER - 1].relax_diag.final_contact_energy;
	result->final_backbone_energy = per_iter[HK_BLIND_P9016_CAND_N_ITER - 1].relax_diag.final_backbone_energy;
	result->final_repulsion_energy = per_iter[HK_BLIND_P9016_CAND_N_ITER - 1].relax_diag.final_repulsion_energy;
	result->final_sep_energy = per_iter[HK_BLIND_P9016_CAND_N_ITER - 1].relax_diag.final_sep_energy;
	result->final_total_energy = per_iter[HK_BLIND_P9016_CAND_N_ITER - 1].relax_diag.final_total_energy;
	result->final_force_l1 = per_iter[HK_BLIND_P9016_CAND_N_ITER - 1].relax_diag.final_force_l1;
	failed |= compute_bpair_metrics(set, &result->bpair);
	failed |= compute_sep_metrics(bmap, diploid, &result->sep);
	failed |= validate_result_metrics(result);
	if (failed) goto cleanup;

	fprintf(stderr, "P9016 candidate %s: writing outputs to %s\n",
			config->name, result->output_dir);
	failed |= check_i32("candidate write required outputs",
						write_required_outputs(posterior_path, coords_path, diag_path,
											   bmap, set, diploid, &result->loop_diag), 0);
	if (write_raw) {
		fprintf(stderr, "P9016 candidate %s: HK_BLIND_WRITE_RAW=1, writing raw posterior\n",
				config->name);
		failed |= check_i32("candidate write raw output",
							write_optional_raw_output(raw_path, raw, n_raw, bmap, set), 0);
	}
	failed |= check_i32("candidate write manifest",
						write_manifest(manifest_path, result->output_dir, posterior_path, coords_path,
									   diag_path, raw_path, write_raw, result), 0);
	if (failed) goto cleanup;

	failed |= validate_file_prefix("candidate bpair posterior", posterior_path, "p00");
	failed |= validate_file_prefix("candidate coords", coords_path, "diploid_bid");
	failed |= validate_file_prefix("candidate loop diag", diag_path, "final_mean_entropy");
	failed |= validate_manifest_prefix(manifest_path, write_raw, config->name);
	if (write_raw)
		failed |= validate_file_prefix("candidate raw posterior", raw_path, "raw_id");
	if (failed) goto cleanup;

	result->audit_ok = run_auditor(result->output_dir) == 0;
	failed |= check_true("candidate audit OK", result->audit_ok);

cleanup:
	t1 = clock();
	result->elapsed_cpu_sec = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
	result->status_ok = failed == 0 &&
		result->loop_diag.n_bad_iter == 0 &&
		result->loop_diag.n_relax_nonfinite_iter == 0 &&
		result->loop_diag.n_coord_nonfinite == 0 &&
		result->audit_ok;
	if (append_summary_row(summary_fp, result) != 0)
		failed = 1;
	fflush(summary_fp);
	fprintf(stderr,
			"P9016 candidate summary: name=%s mean_pU=%.6g mean_pmax=%.6g "
			"mean_margin=%.6g sep_mean=%.6g contact_energy_per_wedge_k=%.6g "
			"repulsion_over_contact=%.6g output_dir=%s audit=%s status=%s\n",
			config->name, result->bpair.mean_pU, result->bpair.mean_pmax,
			result->bpair.mean_margin, result->sep.mean,
			result->loop_diag.final_sum_wedge_k > 0.0?
			(double)result->final_contact_energy / result->loop_diag.final_sum_wedge_k : 0.0,
			result->final_contact_energy > 0.0f?
			(double)result->final_repulsion_energy / (double)result->final_contact_energy : 0.0,
			result->output_dir, result->audit_ok? "OK" : "FAIL",
			result->status_ok? "OK" : "FAIL");
	if (failed)
		fprintf(stderr, "P9016 candidate %s failed; partial output remains in %s\n",
				config->name, result->output_dir);
	free(diploid);
	free(per_iter);
	hk_blind_bpair_set_destroy(set);
	return failed != 0;
}

int main(void)
{
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	fvec3_t *haploid = 0;
	struct candidate_result results[sizeof(candidate_configs) / sizeof(candidate_configs[0])];
	char root_dir[512] = {0};
	char summary_path[768];
	FILE *summary_fp = 0;
	int write_raw = 0;
	int32_t n_raw;
	int32_t i;
	size_t c;
	int failed = 0;
	clock_t t0 = clock(), t1;
	double elapsed_sec;

	if (!file_exists(HK_BLIND_P9016_CAND_PATH)) {
		fprintf(stderr, "ERROR: %s not found; candidate output run did not start\n",
				HK_BLIND_P9016_CAND_PATH);
		return 1;
	}
	if (!file_exists("./audit_blind_p9016_full_cpu_output.bin")) {
		fprintf(stderr, "ERROR: audit_blind_p9016_full_cpu_output.bin not found; build target dependency first\n");
		return 1;
	}
	write_raw = (getenv("HK_BLIND_WRITE_RAW") != 0 &&
				 strcmp(getenv("HK_BLIND_WRITE_RAW"), "1") == 0);

	hk_verbose = 0;
	fprintf(stderr, "P9016 candidate outputs: reading %s\n", HK_BLIND_P9016_CAND_PATH);
	m = hk_map_read(HK_BLIND_P9016_CAND_PATH);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", HK_BLIND_P9016_CAND_PATH);
		return 1;
	}
	n_raw = m->n_pairs;
	if (n_raw <= 0) {
		fprintf(stderr, "candidate input raw positive: predicate failed\n");
		failed = 1;
		goto cleanup;
	}

	fprintf(stderr, "P9016 candidate outputs: building 1Mb bmap from %d raw pairs\n", n_raw);
	bmap = hk_bmap_gen(m->d, n_raw, m->pairs, HK_BLIND_P9016_CAND_RESOLUTION, 1);
	raw = (struct hk_blind_pair*)calloc((size_t)n_raw, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate candidate bmap/raw data\n");
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n_raw; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);
	failed |= check_true("candidate bead count positive", bmap->n_beads > 0);
	if (failed) goto cleanup;

	haploid = (fvec3_t*)calloc((size_t)bmap->n_beads, sizeof(*haploid));
	if (haploid == 0) {
		fprintf(stderr, "failed to allocate candidate haploid scaffold\n");
		failed = 1;
		goto cleanup;
	}
	set_haploid_scaffold(bmap, haploid);

	failed |= check_i32("candidate output root", make_output_root(root_dir, sizeof(root_dir)), 0);
	if (failed) goto cleanup;
	path_join(summary_path, sizeof(summary_path), root_dir, "candidate_summary.tsv");
	summary_fp = fopen(summary_path, "w");
	if (summary_fp == 0) {
		fprintf(stderr, "failed to open candidate summary %s\n", summary_path);
		failed = 1;
		goto cleanup;
	}
	failed |= check_i32("candidate summary header", write_summary_header(summary_fp), 0);
	if (failed) goto cleanup;

	fprintf(stderr, "P9016 candidate outputs: output root %s write_raw_posterior=%d\n",
			root_dir, write_raw);
	memset(results, 0, sizeof(results));
	for (c = 0; c < sizeof(candidate_configs) / sizeof(candidate_configs[0]); ++c) {
		failed |= run_one_candidate(&candidate_configs[c], bmap, raw, n_raw, haploid,
									root_dir, write_raw, summary_fp, &results[c]);
	}
	if (fclose(summary_fp) != 0) {
		summary_fp = 0;
		failed = 1;
		goto cleanup;
	}
	summary_fp = 0;
	failed |= validate_summary_file(summary_path,
									(int)(sizeof(candidate_configs) / sizeof(candidate_configs[0])));
	if (failed) goto cleanup;

	t1 = clock();
	elapsed_sec = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
	fprintf(stderr, "P9016 candidate outputs: output_root=%s elapsed_cpu_sec=%.3f status=OK\n",
			root_dir, elapsed_sec);
	for (c = 0; c < sizeof(candidate_configs) / sizeof(candidate_configs[0]); ++c) {
		fprintf(stderr,
				"P9016 candidate final: name=%s mean_pU=%.6g mean_pmax=%.6g "
				"mean_margin=%.6g sep_mean=%.6g audit=%s output_dir=%s\n",
				results[c].config->name, results[c].bpair.mean_pU,
				results[c].bpair.mean_pmax, results[c].bpair.mean_margin,
				results[c].sep.mean, results[c].audit_ok? "OK" : "FAIL",
				results[c].output_dir);
	}

cleanup:
	if (summary_fp) fclose(summary_fp);
	if (failed)
		fprintf(stderr, "P9016 candidate outputs failed%s%s\n",
				root_dir[0]? "; partial output remains in " : "",
				root_dir[0]? root_dir : "");
	free(haploid);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
