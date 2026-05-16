#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include "hickit.h"

#define HK_BLIND_P9016_SWEEP_PATH "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_SWEEP_RESOLUTION 1000000
#define HK_BLIND_P9016_SWEEP_UNIT 1.0f
#define HK_BLIND_P9016_SWEEP_D_SCALE 1.0f
#define HK_BLIND_P9016_SWEEP_BASE_K 2.0f
#define HK_BLIND_P9016_SWEEP_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_SWEEP_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_SWEEP_INIT_EPS 0.5f
#define HK_BLIND_P9016_SWEEP_INIT_NOISE_SCALE 0.0f
#define HK_BLIND_P9016_SWEEP_INIT_SEED 17ULL

struct sweep_config {
	const char *name;
	int32_t n_iter;
	int32_t relax_steps;
	float relax_step;
	float temperature_start;
	float temperature_end;
	float rho_train_start;
	float rho_train_end;
};

struct bpair_dist_metrics {
	float mean_pU, median_pU;
	float frac_pU_lt_0_90, frac_pU_lt_0_75, frac_pU_lt_0_50;
	float mean_pmax, median_pmax;
	float frac_pmax_gt_0_50, frac_pmax_gt_0_75;
	float mean_margin, median_margin;
	float frac_margin_gt_0_10, frac_margin_gt_0_25;
	float mean_psame, mean_pcross;
	float frac_psame_gt_0_75, frac_pcross_gt_0_75;
};

struct sep_dist_metrics {
	float mean, median, min, max;
	float frac_lt_0_25, frac_lt_0_50, frac_lt_1_00;
};

struct sweep_result {
	const struct sweep_config *config;
	int32_t n_raw;
	int32_t n_bpair;
	int32_t n_beads;
	struct hk_blind_iter_loop_diag loop_diag;
	struct bpair_dist_metrics bpair_metrics;
	struct sep_dist_metrics sep_metrics;
	float final_contact_energy;
	float final_backbone_energy;
	float final_repulsion_energy;
	float final_sep_energy;
	float final_total_energy;
	float final_force_l1;
	int status_ok;
	double elapsed_cpu_sec;
};

static const struct sweep_config sweep_configs[] = {
	{ "n3_rs20_step0p0005_T2to1", 3, 20, 0.0005f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n5_rs20_step0p0005_T2to1", 5, 20, 0.0005f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n3_rs50_step0p0005_T2to1", 3, 50, 0.0005f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n5_rs50_step0p0005_T2to1", 5, 50, 0.0005f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n3_rs100_step0p0005_T2to1", 3, 100, 0.0005f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n5_rs100_step0p0005_T2to1", 5, 100, 0.0005f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n3_rs200_step0p0005_T2to1", 3, 200, 0.0005f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n5_rs200_step0p0005_T2to1", 5, 200, 0.0005f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n3_rs20_step0p001_T2to1", 3, 20, 0.001f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n5_rs20_step0p001_T2to1", 5, 20, 0.001f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n3_rs50_step0p001_T2to1", 3, 50, 0.001f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n5_rs50_step0p001_T2to1", 5, 50, 0.001f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n3_rs100_step0p001_T2to1", 3, 100, 0.001f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n5_rs100_step0p001_T2to1", 5, 100, 0.001f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n3_rs200_step0p001_T2to1", 3, 200, 0.001f, 2.0f, 1.0f, 1.0f, 1.0f },
	{ "n5_rs200_step0p001_T2to1", 5, 200, 0.001f, 2.0f, 1.0f, 1.0f, 1.0f }
};

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

static int env_flag_is_1(const char *name)
{
	const char *v = getenv(name);
	return v != 0 && strcmp(v, "1") == 0;
}

static int parse_sweep_jobs(size_t n_configs)
{
	const char *v = getenv("HK_BLIND_SWEEP_JOBS");
	char *end = 0;
	long jobs;

	if (v == 0 || *v == 0)
		return 1;
	errno = 0;
	jobs = strtol(v, &end, 10);
	if (errno != 0 || end == v || *end != 0 || jobs <= 0) {
		fprintf(stderr, "ERROR: invalid HK_BLIND_SWEEP_JOBS=%s\n", v);
		return -1;
	}
	if ((size_t)jobs > n_configs)
		jobs = (long)n_configs;
	return (int)jobs;
}

static int check_true(const char *label, int pred)
{
	if (!pred) {
		fprintf(stderr, "%s: predicate failed\n", label);
		return 1;
	}
	return 0;
}

static int check_i32(const char *label, int32_t got, int32_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %d, expected %d\n", label, got, expected);
		return 1;
	}
	return 0;
}

static int make_output_root(char *dir, size_t dir_size)
{
	long pid = (long)getpid();
	long stamp = (long)time(0);
	int i;

	for (i = 0; i < 100; ++i) {
		snprintf(dir, dir_size, "/tmp/hk_blind_p9016_full_cpu_sweep_%ld_%ld_%d",
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
			failed |= check_true("sweep coord finite", isfinite(coords[i][a]));
	return failed;
}

static void set_schedule_conf(struct hk_blind_iter_schedule_conf *conf,
							  const struct sweep_config *config)
{
	conf->n_iter = config->n_iter;
	conf->base_conf.unit = HK_BLIND_P9016_SWEEP_UNIT;
	conf->base_conf.d_scale = HK_BLIND_P9016_SWEEP_D_SCALE;
	conf->base_conf.base_k = HK_BLIND_P9016_SWEEP_BASE_K;
	conf->base_conf.temperature = config->temperature_start;
	conf->base_conf.rho_train = config->rho_train_start;
	conf->base_conf.min_sep_unit = HK_BLIND_P9016_SWEEP_MIN_SEP_UNIT;
	conf->base_conf.lambda_sep = HK_BLIND_P9016_SWEEP_LAMBDA_SEP;
	conf->base_conf.relax_step = config->relax_step;
	conf->base_conf.relax_steps = config->relax_steps;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->base_conf.d_scale_eps_count = 1e-6f;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = config->temperature_start;
	conf->temperature_end = config->temperature_end;
	conf->rho_train_start = config->rho_train_start;
	conf->rho_train_end = config->rho_train_end;
}

static int check_binned_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
							int32_t n_raw)
{
	int failed = 0;
	int32_t i;

	failed |= check_true("sweep set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_i32("sweep raw count", set->n_raw, n_raw);
	failed |= check_true("sweep bpair positive", set->n_bpairs > 0);
	for (i = 0; i < set->n_raw; ++i) {
		failed |= check_true("sweep raw2binned id lower", set->raw2binned[i].bpair_id >= 0);
		failed |= check_true("sweep raw2binned id upper", set->raw2binned[i].bpair_id < set->n_bpairs);
		failed |= check_true("sweep raw2binned swapped",
							 set->raw2binned[i].swapped == 0 || set->raw2binned[i].swapped == 1);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		failed |= check_true("sweep key sorted", bp->key.bid[0] <= bp->key.bid[1]);
		failed |= check_true("sweep key bid0 valid", bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		failed |= check_true("sweep key bid1 valid", bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
	}
	return failed;
}

static int validate_final_bpair_set(const struct hk_blind_bpair_set *set)
{
	struct hk_blind_iter_diag diag;
	int failed = 0;

	hk_blind_iter_diag_init(&diag);
	hk_blind_iter_diag_validate_bpair_set(set, &diag);
	failed |= check_i32("sweep final posterior nonfinite", diag.n_posterior_nonfinite, 0);
	failed |= check_i32("sweep final posterior bad sum", diag.n_posterior_bad_sum, 0);
	failed |= check_i32("sweep final posterior out of range", diag.n_posterior_out_of_range, 0);
	failed |= check_i32("sweep final uncertainty nonfinite", diag.n_uncertainty_nonfinite, 0);
	failed |= check_i32("sweep final uncertainty out of range", diag.n_uncertainty_out_of_range, 0);
	failed |= check_i32("sweep final five-state bad sum", diag.n_five_state_bad_sum, 0);
	return failed;
}

static int validate_loop_health(const struct sweep_config *config,
								const struct hk_blind_iter_loop_diag *diag,
								const struct hk_blind_single_iter_diag *per_iter)
{
	const struct hk_blind_iter_diag *last_pre;
	int failed = 0;

	failed |= check_i32("sweep completed", diag->n_completed, config->n_iter);
	failed |= check_i32("sweep bad iter", diag->n_bad_iter, 0);
	failed |= check_i32("sweep relax nonfinite iter", diag->n_relax_nonfinite_iter, 0);
	failed |= check_i32("sweep coord nonfinite", diag->n_coord_nonfinite, 0);
	failed |= check_true("sweep entropy finite", isfinite(diag->final_mean_entropy));
	failed |= check_true("sweep pU finite", isfinite(diag->final_mean_pU));
	failed |= check_true("sweep pU range", diag->final_mean_pU >= -1e-6f && diag->final_mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("sweep sum wedge finite", isfinite(diag->final_sum_wedge_k));
	failed |= check_true("sweep sum wedge nonnegative", diag->final_sum_wedge_k >= 0.0);
	if (config->n_iter <= 0) return failed;
	last_pre = &per_iter[config->n_iter - 1].pre_relax_diag;
	failed |= check_i32("sweep final posterior nonfinite diag", last_pre->n_posterior_nonfinite, 0);
	failed |= check_i32("sweep final posterior bad sum diag", last_pre->n_posterior_bad_sum, 0);
	failed |= check_i32("sweep final posterior out range diag", last_pre->n_posterior_out_of_range, 0);
	failed |= check_i32("sweep final uncertainty nonfinite diag", last_pre->n_uncertainty_nonfinite, 0);
	failed |= check_i32("sweep final uncertainty out range diag", last_pre->n_uncertainty_out_of_range, 0);
	failed |= check_i32("sweep final five-state bad sum diag", last_pre->n_five_state_bad_sum, 0);
	failed |= check_i32("sweep final repulsion nonfinite step",
						per_iter[config->n_iter - 1].relax_diag.n_repulsion_nonfinite_step, 0);
	return failed;
}

static int compute_bpair_metrics(const struct hk_blind_bpair_set *set,
								 struct bpair_dist_metrics *out)
{
	float *pU = 0, *pmax = 0, *margin = 0, *psame = 0, *pcross = 0;
	double sum_pU = 0.0, sum_pmax = 0.0, sum_margin = 0.0, sum_psame = 0.0, sum_pcross = 0.0;
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
	psame = (float*)malloc((size_t)n * sizeof(*psame));
	pcross = (float*)malloc((size_t)n * sizeof(*pcross));
	if (pU == 0 || pmax == 0 || margin == 0 || psame == 0 || pcross == 0) {
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		float same = bp->p4[HK_BLIND_STATE_00] + bp->p4[HK_BLIND_STATE_11];
		float cross = bp->p4[HK_BLIND_STATE_01] + bp->p4[HK_BLIND_STATE_10];
		pU[i] = bp->pU;
		pmax[i] = bp->pmax;
		margin[i] = bp->margin;
		psame[i] = same;
		pcross[i] = cross;
		failed |= check_true("sweep metric pU finite", isfinite(pU[i]));
		failed |= check_true("sweep metric pmax finite", isfinite(pmax[i]));
		failed |= check_true("sweep metric margin finite", isfinite(margin[i]));
		failed |= check_true("sweep metric psame finite", isfinite(same));
		failed |= check_true("sweep metric pcross finite", isfinite(cross));
		if (failed) continue;
		sum_pU += pU[i];
		sum_pmax += pmax[i];
		sum_margin += margin[i];
		sum_psame += same;
		sum_pcross += cross;
		if (pU[i] < 0.90f) out->frac_pU_lt_0_90 += 1.0f;
		if (pU[i] < 0.75f) out->frac_pU_lt_0_75 += 1.0f;
		if (pU[i] < 0.50f) out->frac_pU_lt_0_50 += 1.0f;
		if (pmax[i] > 0.50f) out->frac_pmax_gt_0_50 += 1.0f;
		if (pmax[i] > 0.75f) out->frac_pmax_gt_0_75 += 1.0f;
		if (margin[i] > 0.10f) out->frac_margin_gt_0_10 += 1.0f;
		if (margin[i] > 0.25f) out->frac_margin_gt_0_25 += 1.0f;
		if (same > 0.75f) out->frac_psame_gt_0_75 += 1.0f;
		if (cross > 0.75f) out->frac_pcross_gt_0_75 += 1.0f;
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
	out->frac_pU_lt_0_90 /= (float)n;
	out->frac_pU_lt_0_75 /= (float)n;
	out->frac_pU_lt_0_50 /= (float)n;
	out->frac_pmax_gt_0_50 /= (float)n;
	out->frac_pmax_gt_0_75 /= (float)n;
	out->frac_margin_gt_0_10 /= (float)n;
	out->frac_margin_gt_0_25 /= (float)n;
	out->frac_psame_gt_0_75 /= (float)n;
	out->frac_pcross_gt_0_75 /= (float)n;
	failed |= check_true("sweep mean pU range", out->mean_pU >= -1e-6f && out->mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("sweep median pU range", out->median_pU >= -1e-6f && out->median_pU <= 1.0f + 1e-6f);
	failed |= check_true("sweep mean pmax range", out->mean_pmax >= 0.25f - 1e-6f && out->mean_pmax <= 1.0f + 1e-6f);
	failed |= check_true("sweep mean margin range", out->mean_margin >= -1e-6f && out->mean_margin <= 1.0f + 1e-6f);
	failed |= check_true("sweep mean psame range", out->mean_psame >= -1e-6f && out->mean_psame <= 1.0f + 1e-6f);
	failed |= check_true("sweep mean pcross range", out->mean_pcross >= -1e-6f && out->mean_pcross <= 1.0f + 1e-6f);

cleanup:
	free(pU);
	free(pmax);
	free(margin);
	free(psame);
	free(pcross);
	return failed;
}

static int compute_sep_metrics(const struct hk_bmap *bmap, const fvec3_t *coords,
							   struct sep_dist_metrics *out)
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
	out->max = 0.0f;
	for (i = 0; i < n; ++i) {
		int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		sep[i] = dist3(coords[b0], coords[b1]);
		failed |= check_true("sweep sep finite", isfinite(sep[i]));
		if (!isfinite(sep[i]))
			continue;
		sum += sep[i];
		if (sep[i] < out->min) out->min = sep[i];
		if (sep[i] > out->max) out->max = sep[i];
		if (sep[i] < 0.25f) out->frac_lt_0_25 += 1.0f;
		if (sep[i] < 0.50f) out->frac_lt_0_50 += 1.0f;
		if (sep[i] < 1.00f) out->frac_lt_1_00 += 1.0f;
	}
	if (failed) goto cleanup;
	out->mean = (float)(sum / n);
	out->median = median_float(sep, n);
	out->frac_lt_0_25 /= (float)n;
	out->frac_lt_0_50 /= (float)n;
	out->frac_lt_1_00 /= (float)n;
	failed |= check_true("sweep sep mean finite", isfinite(out->mean));
	failed |= check_true("sweep sep median finite", isfinite(out->median));
	failed |= check_true("sweep sep min finite", isfinite(out->min));
	failed |= check_true("sweep sep max finite", isfinite(out->max));

cleanup:
	free(sep);
	return failed;
}

static int validate_result_metrics(const struct sweep_result *result)
{
	const struct bpair_dist_metrics *b = &result->bpair_metrics;
	const struct sep_dist_metrics *s = &result->sep_metrics;
	int failed = 0;

	failed |= check_true("sweep final contact energy finite", isfinite(result->final_contact_energy));
	failed |= check_true("sweep final backbone energy finite", isfinite(result->final_backbone_energy));
	failed |= check_true("sweep final repulsion energy finite", isfinite(result->final_repulsion_energy));
	failed |= check_true("sweep final sep energy finite", isfinite(result->final_sep_energy));
	failed |= check_true("sweep final total energy finite", isfinite(result->final_total_energy));
	failed |= check_true("sweep final force finite", isfinite(result->final_force_l1));
	failed |= check_true("sweep mean pU finite", isfinite(b->mean_pU));
	failed |= check_true("sweep median pU finite", isfinite(b->median_pU));
	failed |= check_true("sweep mean pmax finite", isfinite(b->mean_pmax));
	failed |= check_true("sweep mean margin finite", isfinite(b->mean_margin));
	failed |= check_true("sweep mean psame finite", isfinite(b->mean_psame));
	failed |= check_true("sweep mean pcross finite", isfinite(b->mean_pcross));
	failed |= check_true("sweep sep mean finite result", isfinite(s->mean));
	failed |= check_true("sweep sep median finite result", isfinite(s->median));
	failed |= check_true("sweep sep min finite result", isfinite(s->min));
	failed |= check_true("sweep sep max finite result", isfinite(s->max));
	return failed;
}

static int write_bpair_output(const char *path, const struct hk_bmap *bmap,
							  const struct hk_blind_bpair_set *set)
{
	FILE *fp = fopen(path, "w");
	int ret;
	if (fp == 0) return -1;
	ret = hk_blind_write_bpair_posterior_tsv(fp, bmap, set);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int write_coords_output(const char *path, const struct hk_bmap *bmap, const fvec3_t *coords)
{
	FILE *fp = fopen(path, "w");
	int ret;
	if (fp == 0) return -1;
	ret = hk_blind_write_diploid_coords_tsv(fp, bmap, coords);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int write_loop_diag_output(const char *path, const struct hk_blind_iter_loop_diag *diag)
{
	FILE *fp = fopen(path, "w");
	int ret;
	if (fp == 0) return -1;
	ret = hk_blind_write_iter_loop_diag_tsv(fp, diag);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int write_raw_output(const char *path, const struct hk_blind_pair *raw, int32_t n_raw,
							const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set)
{
	FILE *fp = fopen(path, "w");
	int ret;
	if (fp == 0) return -1;
	ret = hk_blind_write_raw_contact_posterior_tsv(fp, raw, n_raw, bmap, set);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int write_manifest(const char *manifest_path, const char *out_dir,
						  const char *posterior_path, const char *coords_path,
						  const char *diag_path, const char *raw_path,
						  int write_raw, const struct sweep_config *config,
						  const struct hk_bmap *bmap,
						  const struct hk_blind_bpair_set *set,
						  const struct hk_blind_iter_loop_diag *loop_diag,
						  const char *status)
{
	FILE *fp = fopen(manifest_path, "w");
	struct hk_blind_base_k_stats base_k_stats;
	if (fp == 0) return -1;
	hk_blind_bpair_set_base_k_stats(set, &base_k_stats);
	if (fprintf(fp,
				"key\tvalue\n"
				"sample\tP9016\n"
				"runner_family\tfull_cpu_sweep\n"
				"runner_version\t2026-04-30\n"
				"default_profile\tp9016_full_cpu_sweep_v1\n"
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
				HK_BLIND_P9016_SWEEP_PATH, out_dir, set->n_raw, set->n_bpairs, bmap->n_beads,
				HK_BLIND_P9016_SWEEP_RESOLUTION, config->n_iter,
				HK_BLIND_P9016_SWEEP_UNIT, HK_BLIND_P9016_SWEEP_D_SCALE,
				base_k_stats.n > 0 && base_k_stats.min == base_k_stats.max? "uniform" : "variable",
				base_k_stats.mean, base_k_stats.min, base_k_stats.mean, base_k_stats.max,
				base_k_stats.n_nonfinite, HK_BLIND_P9016_SWEEP_BASE_K,
				hk_blind_rho_train_mode_name(HK_BLIND_RHO_TRAIN_CONSTANT),
				hk_blind_d_scale_mode_name(HK_BLIND_D_SCALE_RAW_COUNT),
				set->same_bin_filter_enabled, (long long)set->n_raw_same_bin_excluded,
				set->n_bpair_same_bin_excluded, HK_BLIND_P9016_SWEEP_MIN_SEP_UNIT,
				HK_BLIND_P9016_SWEEP_LAMBDA_SEP, config->relax_step,
				config->relax_steps, 1, HK_BLIND_REPULSION_CELL,
				config->temperature_start, config->temperature_end,
				config->rho_train_start, config->rho_train_end,
				HK_BLIND_P9016_SWEEP_INIT_EPS,
				HK_BLIND_P9016_SWEEP_INIT_NOISE_SCALE,
				(unsigned long long)HK_BLIND_P9016_SWEEP_INIT_SEED, write_raw,
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

static int run_auditor(const char *out_dir)
{
	char cmd[1024];
	int ret;
	snprintf(cmd, sizeof(cmd), "./audit_blind_p9016_full_cpu_output.bin \"%s\"", out_dir);
	ret = system(cmd);
	if (ret != 0) {
		fprintf(stderr, "sweep audit failed for %s with status %d\n", out_dir, ret);
		return 1;
	}
	return 0;
}

static int write_full_outputs(const char *root_dir, const struct sweep_config *config,
							  const struct hk_bmap *bmap,
							  const struct hk_blind_pair *raw, int32_t n_raw,
							  const struct hk_blind_bpair_set *set,
							  const fvec3_t *coords,
							  const struct hk_blind_iter_loop_diag *loop_diag,
							  int write_raw, const char *status)
{
	char out_dir[768], posterior_path[1024], coords_path[1024], diag_path[1024], raw_path[1024], manifest_path[1024];
	int failed = 0;

	path_join(out_dir, sizeof(out_dir), root_dir, config->name);
	if (mkdir(out_dir, 0700) != 0) {
		fprintf(stderr, "failed to create sweep full-output directory %s\n", out_dir);
		return 1;
	}
	path_join(posterior_path, sizeof(posterior_path), out_dir, "p9016_full.bpair_posterior.tsv");
	path_join(coords_path, sizeof(coords_path), out_dir, "p9016_full.coords.tsv");
	path_join(diag_path, sizeof(diag_path), out_dir, "p9016_full.loop_diag.tsv");
	path_join(raw_path, sizeof(raw_path), out_dir, "p9016_full.raw_posterior.tsv");
	path_join(manifest_path, sizeof(manifest_path), out_dir, "p9016_full.manifest.tsv");
	failed |= check_i32("sweep full bpair output", write_bpair_output(posterior_path, bmap, set), 0);
	failed |= check_i32("sweep full coords output", write_coords_output(coords_path, bmap, coords), 0);
	failed |= check_i32("sweep full loop diag output", write_loop_diag_output(diag_path, loop_diag), 0);
	if (write_raw)
		failed |= check_i32("sweep full raw output", write_raw_output(raw_path, raw, n_raw, bmap, set), 0);
	failed |= check_i32("sweep full manifest", write_manifest(manifest_path, out_dir, posterior_path,
															  coords_path, diag_path, raw_path,
															  write_raw, config, bmap, set, loop_diag, status), 0);
	failed |= check_true("sweep full posterior exists", file_size_or_negative(posterior_path) > 0);
	failed |= check_true("sweep full coords exists", file_size_or_negative(coords_path) > 0);
	failed |= check_true("sweep full loop diag exists", file_size_or_negative(diag_path) > 0);
	failed |= check_true("sweep full manifest exists", file_size_or_negative(manifest_path) > 0);
	if (write_raw)
		failed |= check_true("sweep full raw exists", file_size_or_negative(raw_path) > 0);
	if (!failed && strcmp(status, "OK") == 0)
		failed |= run_auditor(out_dir);
	return failed;
}

static int write_summary_header(FILE *fp)
{
	return fprintf(fp,
				   "config_name\tn_raw\tn_bpair\tn_beads\tn_iter\trelax_steps\trelax_step\t"
				   "temperature_start\ttemperature_end\trho_train_start\trho_train_end\t"
				   "unit\td_scale\tlegacy_base_k_unused\tmin_sep_unit\tlambda_sep\tinit_eps\tinit_noise_scale\tinit_seed\t"
				   "write_full_outputs\twrite_raw_posterior\trepulsion_mode\tn_completed\tn_bad_iter\t"
				   "n_relax_nonfinite_iter\tn_coord_nonfinite\ttotal_chr_flipped\t"
				   "final_mean_entropy\tfinal_mean_pU\tfinal_sum_wedge_k\t"
				   "mean_pU\tmedian_pU\tfrac_pU_lt_0_90\tfrac_pU_lt_0_75\tfrac_pU_lt_0_50\t"
				   "mean_pmax\tmedian_pmax\tfrac_pmax_gt_0_50\tfrac_pmax_gt_0_75\t"
				   "mean_margin\tmedian_margin\tfrac_margin_gt_0_10\tfrac_margin_gt_0_25\t"
				   "mean_psame\tmean_pcross\tfrac_psame_gt_0_75\tfrac_pcross_gt_0_75\t"
				   "sep_mean\tsep_median\tsep_min\tsep_max\tfrac_sep_lt_0_25\t"
				   "frac_sep_lt_0_50\tfrac_sep_lt_1_00\tfinal_contact_energy\t"
				   "final_backbone_energy\tfinal_repulsion_energy\tfinal_sep_energy\tfinal_total_energy\t"
				   "final_force_l1\tstatus\telapsed_cpu_sec\n") < 0? -1 : 0;
}

static int append_summary_row(FILE *fp, const struct sweep_result *result,
							  int write_full, int write_raw)
{
	const struct sweep_config *c = result->config;
	const struct hk_blind_iter_loop_diag *d = &result->loop_diag;
	const struct bpair_dist_metrics *b = &result->bpair_metrics;
	const struct sep_dist_metrics *s = &result->sep_metrics;
	const char *status = result->status_ok? "OK" : "FAIL";

	return fprintf(fp,
				   "%s\t%d\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%llu\t%d\t%d\t%d\t"
				   "%d\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.17g\t"
					   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
					   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%s\t%.3f\n",
				   c->name, result->n_raw, result->n_bpair, result->n_beads,
				   c->n_iter, c->relax_steps, c->relax_step,
				   c->temperature_start, c->temperature_end,
				   c->rho_train_start, c->rho_train_end,
				   HK_BLIND_P9016_SWEEP_UNIT, HK_BLIND_P9016_SWEEP_D_SCALE,
				   HK_BLIND_P9016_SWEEP_BASE_K, HK_BLIND_P9016_SWEEP_MIN_SEP_UNIT,
				   HK_BLIND_P9016_SWEEP_LAMBDA_SEP, HK_BLIND_P9016_SWEEP_INIT_EPS,
				   HK_BLIND_P9016_SWEEP_INIT_NOISE_SCALE,
				   (unsigned long long)HK_BLIND_P9016_SWEEP_INIT_SEED, write_full, write_raw,
				   d->repulsion_mode,
				   d->n_completed, d->n_bad_iter, d->n_relax_nonfinite_iter,
				   d->n_coord_nonfinite, d->total_chr_flipped,
				   d->final_mean_entropy, d->final_mean_pU, d->final_sum_wedge_k,
				   b->mean_pU, b->median_pU, b->frac_pU_lt_0_90,
				   b->frac_pU_lt_0_75, b->frac_pU_lt_0_50,
				   b->mean_pmax, b->median_pmax, b->frac_pmax_gt_0_50,
				   b->frac_pmax_gt_0_75, b->mean_margin, b->median_margin,
				   b->frac_margin_gt_0_10, b->frac_margin_gt_0_25,
				   b->mean_psame, b->mean_pcross, b->frac_psame_gt_0_75,
				   b->frac_pcross_gt_0_75, s->mean, s->median, s->min, s->max,
				   s->frac_lt_0_25, s->frac_lt_0_50, s->frac_lt_1_00,
				   result->final_contact_energy, result->final_backbone_energy,
				   result->final_repulsion_energy, result->final_sep_energy, result->final_total_energy,
				   result->final_force_l1, status, result->elapsed_cpu_sec) < 0? -1 : 0;
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
			failed |= check_true("sweep summary header config", strstr(line, "config_name") != 0);
			failed |= check_true("sweep summary header mean pU", strstr(line, "mean_pU") != 0);
			failed |= check_true("sweep summary header sep median", strstr(line, "sep_median") != 0);
			failed |= check_true("sweep summary header status", strstr(line, "status") != 0);
		} else if (strstr(line, "\tFAIL\t") != 0) {
			++n_fail;
		}
		++n_rows;
	}
	if (ferror(fp)) failed = 1;
	fclose(fp);
	failed |= check_i32("sweep summary row count", n_rows, n_expected);
	failed |= check_i32("sweep summary fail count", n_fail, 0);
	return failed;
}

static int run_one_config(const struct sweep_config *config, const struct hk_bmap *bmap,
						  const struct hk_blind_pair *raw, int32_t n_raw,
						  const fvec3_t *haploid, const char *root_dir,
						  int write_full, int write_raw, FILE *summary_fp,
						  struct sweep_result *result)
{
	struct hk_blind_bpair_set *set = 0;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_single_iter_diag *per_iter = 0;
	fvec3_t *diploid = 0;
	int32_t n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	int failed = 0;
	int ret;
	clock_t t0 = clock(), t1;
	const char *status;

	memset(result, 0, sizeof(*result));
	result->config = config;
	result->n_raw = n_raw;
	result->n_beads = bmap->n_beads;

	fprintf(stderr, "P9016 sweep config %s: n_iter=%d relax_steps=%d relax_step=%.4g\n",
			config->name, config->n_iter, config->relax_steps, config->relax_step);
	set = hk_blind_bpair_set_build(bmap, n_raw, raw);
	failed |= check_binned_set(bmap, set, n_raw);
	if (set) result->n_bpair = set->n_bpairs;
	if (failed) goto finish;

	diploid = (fvec3_t*)calloc((size_t)n_diploid, sizeof(*diploid));
	per_iter = (struct hk_blind_single_iter_diag*)calloc((size_t)config->n_iter, sizeof(*per_iter));
	if (diploid == 0 || per_iter == 0) {
		fprintf(stderr, "failed to allocate sweep coordinate/diagnostic data for %s\n", config->name);
		failed = 1;
		goto finish;
	}
	failed |= check_i32("sweep init diploid",
						hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
																  diploid, HK_BLIND_P9016_SWEEP_INIT_EPS,
																  HK_BLIND_P9016_SWEEP_INIT_NOISE_SCALE,
																  HK_BLIND_P9016_SWEEP_INIT_SEED), 0);
	failed |= check_coords_finite(diploid, n_diploid);
	if (failed) goto finish;

	hk_fdg_conf_init(&fdg_conf);
	set_schedule_conf(&schedule_conf, config);
	ret = hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &fdg_conf, diploid, 0,
											   &schedule_conf, per_iter, &result->loop_diag);
	failed |= check_i32("sweep scheduled ret", ret, 0);
	if (ret == 0) {
		failed |= validate_loop_health(config, &result->loop_diag, per_iter);
		failed |= check_coords_finite(diploid, n_diploid);
		failed |= validate_final_bpair_set(set);
		if (!failed) {
			const struct hk_blind_relax_diag *rd = &per_iter[config->n_iter - 1].relax_diag;
			result->final_contact_energy = rd->final_contact_energy;
			result->final_backbone_energy = rd->final_backbone_energy;
			result->final_repulsion_energy = rd->final_repulsion_energy;
			result->final_sep_energy = rd->final_sep_energy;
			result->final_total_energy = rd->final_total_energy;
			result->final_force_l1 = rd->final_force_l1;
			failed |= compute_bpair_metrics(set, &result->bpair_metrics);
			failed |= compute_sep_metrics(bmap, diploid, &result->sep_metrics);
			failed |= validate_result_metrics(result);
		}
	}

finish:
	t1 = clock();
	result->elapsed_cpu_sec = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
	if (write_full && set && diploid)
		failed |= write_full_outputs(root_dir, config, bmap, raw, n_raw, set, diploid,
									 &result->loop_diag, write_raw, failed == 0? "OK" : "FAIL");
	result->status_ok = failed == 0;
	status = result->status_ok? "OK" : "FAIL";
	if (append_summary_row(summary_fp, result, write_full, write_raw) != 0)
		failed = 1;
	result->status_ok = failed == 0;
	status = result->status_ok? "OK" : "FAIL";
	fflush(summary_fp);
	fprintf(stderr,
			"%s n_iter=%d relax_steps=%d relax_step=%.4g mean_pU=%.6g "
			"frac_pU_lt_0.75=%.6g mean_pmax=%.6g mean_margin=%.6g "
			"sep_mean=%.6g status=%s\n",
			config->name, config->n_iter, config->relax_steps, config->relax_step,
			result->bpair_metrics.mean_pU, result->bpair_metrics.frac_pU_lt_0_75,
			result->bpair_metrics.mean_pmax, result->bpair_metrics.mean_margin,
			result->sep_metrics.mean, status);

	free(diploid);
	free(per_iter);
	hk_blind_bpair_set_destroy(set);
	return failed != 0;
}

static int run_one_config_to_row_file(size_t config_idx, const struct hk_bmap *bmap,
									  const struct hk_blind_pair *raw, int32_t n_raw,
									  const fvec3_t *haploid, const char *root_dir,
									  int write_full, int write_raw, const char *row_path)
{
	FILE *fp = fopen(row_path, "w");
	struct sweep_result result;
	int ret;

	if (fp == 0) {
		fprintf(stderr, "failed to open sweep summary row %s\n", row_path);
		return 1;
	}
	ret = run_one_config(&sweep_configs[config_idx], bmap, raw, n_raw, haploid,
						 root_dir, write_full, write_raw, fp, &result);
	if (fclose(fp) != 0)
		ret = 1;
	return ret != 0;
}

static int append_summary_row_file(FILE *summary_fp, const char *row_path)
{
	FILE *fp = fopen(row_path, "r");
	char line[8192];
	int failed = 0;

	if (fp == 0) {
		fprintf(stderr, "failed to open completed summary row %s\n", row_path);
		return 1;
	}
	while (fgets(line, sizeof(line), fp)) {
		if (fputs(line, summary_fp) == EOF) {
			failed = 1;
			break;
		}
	}
	if (ferror(fp)) failed = 1;
	fclose(fp);
	return failed;
}

static int run_configs_parallel(int jobs, const struct hk_bmap *bmap,
								const struct hk_blind_pair *raw, int32_t n_raw,
								const fvec3_t *haploid, const char *root_dir,
								int write_full, int write_raw, FILE *summary_fp)
{
	size_t n_configs = sizeof(sweep_configs) / sizeof(sweep_configs[0]);
	char (*row_paths)[1024] = 0;
	pid_t *pids = 0;
	int *done = 0;
	size_t next = 0, running = 0;
	int failed = 0;
	size_t i;

	assert(jobs > 0);
	row_paths = (char (*)[1024])calloc(n_configs, sizeof(*row_paths));
	pids = (pid_t*)calloc(n_configs, sizeof(*pids));
	done = (int*)calloc(n_configs, sizeof(*done));
	if (row_paths == 0 || pids == 0 || done == 0) {
		fprintf(stderr, "failed to allocate sweep parallel bookkeeping\n");
		free(row_paths);
		free(pids);
		free(done);
		return 1;
	}
	for (i = 0; i < n_configs; ++i)
		snprintf(row_paths[i], sizeof(row_paths[i]), "%s/%s.summary_row.tsv", root_dir, sweep_configs[i].name);

	while (next < n_configs || running > 0) {
		while (next < n_configs && running < (size_t)jobs) {
			size_t idx = next++;
			pid_t pid = fork();
			if (pid < 0) {
				fprintf(stderr, "failed to fork sweep config %s\n", sweep_configs[idx].name);
				failed = 1;
				continue;
			}
			if (pid == 0) {
				int ret = run_one_config_to_row_file(idx, bmap, raw, n_raw, haploid,
													 root_dir, write_full, write_raw, row_paths[idx]);
				_exit(ret? 1 : 0);
			}
			pids[idx] = pid;
			++running;
		}
		if (running > 0) {
			int status = 0;
			pid_t pid = wait(&status);
			if (pid < 0) {
				fprintf(stderr, "wait failed during sweep parallel run\n");
				failed = 1;
				break;
			}
			--running;
			for (i = 0; i < n_configs; ++i) {
				if (pids[i] == pid) {
					done[i] = 1;
					if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
						fprintf(stderr, "sweep config %s failed in child\n", sweep_configs[i].name);
						failed = 1;
					}
					break;
				}
			}
		}
	}
	while (running > 0) {
		int status = 0;
		pid_t pid = wait(&status);
		if (pid < 0)
			break;
		--running;
		for (i = 0; i < n_configs; ++i) {
			if (pids[i] == pid) {
				done[i] = 1;
				if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
					failed = 1;
				break;
			}
		}
	}
	for (i = 0; i < n_configs; ++i) {
		if (!done[i]) {
			fprintf(stderr, "sweep config %s did not complete\n", sweep_configs[i].name);
			failed = 1;
			continue;
		}
		if (append_summary_row_file(summary_fp, row_paths[i]) != 0) {
			fprintf(stderr, "failed to append summary row for %s\n", sweep_configs[i].name);
			failed = 1;
		}
	}
	free(row_paths);
	free(pids);
	free(done);
	return failed != 0;
}

int main(void)
{
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	fvec3_t *haploid = 0;
	char root_dir[512] = {0};
	char summary_path[768];
	FILE *summary_fp = 0;
	int write_full = env_flag_is_1("HK_BLIND_SWEEP_WRITE_FULL");
	int write_raw = env_flag_is_1("HK_BLIND_WRITE_RAW");
	int jobs = parse_sweep_jobs(sizeof(sweep_configs) / sizeof(sweep_configs[0]));
	int failed = 0;
	int32_t n_raw;
	int32_t i;
	size_t c;
	time_t t0 = time(0), t1;
	double elapsed_sec;

	if (jobs < 0)
		return 1;
	if (!write_full)
		write_raw = 0;
	if (!file_exists(HK_BLIND_P9016_SWEEP_PATH)) {
		fprintf(stderr, "ERROR: %s not found; full P9016 sweep did not start\n",
				HK_BLIND_P9016_SWEEP_PATH);
		return 1;
	}
	if (write_full && !file_exists("./audit_blind_p9016_full_cpu_output.bin")) {
		fprintf(stderr, "ERROR: audit_blind_p9016_full_cpu_output.bin not found; build target dependency first\n");
		return 1;
	}

	hk_verbose = 0;
	fprintf(stderr, "P9016 full CPU sweep: reading %s\n", HK_BLIND_P9016_SWEEP_PATH);
	m = hk_map_read(HK_BLIND_P9016_SWEEP_PATH);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", HK_BLIND_P9016_SWEEP_PATH);
		return 1;
	}
	n_raw = m->n_pairs;
	if (n_raw <= 0) {
		fprintf(stderr, "sweep input raw positive: predicate failed\n");
		failed = 1;
		goto cleanup;
	}

	fprintf(stderr, "P9016 full CPU sweep: building 1Mb bmap from %d raw pairs\n", n_raw);
	bmap = hk_bmap_gen(m->d, n_raw, m->pairs, HK_BLIND_P9016_SWEEP_RESOLUTION, 1);
	raw = (struct hk_blind_pair*)calloc((size_t)n_raw, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate sweep bmap/raw data\n");
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n_raw; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);
	failed |= check_true("sweep bead count positive", bmap->n_beads > 0);
	if (failed) goto cleanup;

	haploid = (fvec3_t*)calloc((size_t)bmap->n_beads, sizeof(*haploid));
	if (haploid == 0) {
		fprintf(stderr, "failed to allocate sweep haploid scaffold\n");
		failed = 1;
		goto cleanup;
	}
	set_haploid_scaffold(bmap, haploid);

	failed |= check_i32("sweep output root", make_output_root(root_dir, sizeof(root_dir)), 0);
	if (failed) goto cleanup;
	path_join(summary_path, sizeof(summary_path), root_dir, "sweep_summary.tsv");
	summary_fp = fopen(summary_path, "w");
	if (summary_fp == 0) {
		fprintf(stderr, "failed to open sweep summary %s\n", summary_path);
		failed = 1;
		goto cleanup;
	}
	failed |= check_i32("sweep summary header", write_summary_header(summary_fp), 0);
	if (failed) goto cleanup;
	fflush(summary_fp);

	fprintf(stderr, "P9016 full CPU sweep: output root %s write_full_outputs=%d write_raw_posterior=%d jobs=%d\n",
			root_dir, write_full, write_raw, jobs);
	if (jobs == 1) {
		for (c = 0; c < sizeof(sweep_configs) / sizeof(sweep_configs[0]); ++c) {
			struct sweep_result result;
			int ret = run_one_config(&sweep_configs[c], bmap, raw, n_raw, haploid,
									 root_dir, write_full, write_raw, summary_fp, &result);
			if (ret != 0)
				failed = 1;
		}
	} else {
		failed |= run_configs_parallel(jobs, bmap, raw, n_raw, haploid,
									   root_dir, write_full, write_raw, summary_fp);
	}
	if (fclose(summary_fp) != 0) {
		summary_fp = 0;
		failed = 1;
		goto cleanup;
	}
	summary_fp = 0;
	failed |= validate_summary_file(summary_path, (int)(sizeof(sweep_configs) / sizeof(sweep_configs[0])));
	if (failed) goto cleanup;

	t1 = time(0);
	elapsed_sec = difftime(t1, t0);
	fprintf(stderr, "P9016 full CPU sweep: output_root=%s elapsed_wall_sec=%.3f status=OK\n",
			root_dir, elapsed_sec);

cleanup:
	if (summary_fp) fclose(summary_fp);
	if (failed)
		fprintf(stderr, "P9016 full CPU sweep failed%s%s\n",
				root_dir[0]? "; partial output remains in " : "",
				root_dir[0]? root_dir : "");
	free(haploid);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
