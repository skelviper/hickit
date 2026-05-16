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

#define HK_BLIND_P9016_REP_SWEEP_PATH "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_REP_SWEEP_RESOLUTION 1000000
#define HK_BLIND_P9016_REP_SWEEP_N_ITER 3
#define HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS 20
#define HK_BLIND_P9016_REP_SWEEP_RELAX_STEP 0.0005f
#define HK_BLIND_P9016_REP_SWEEP_TEMPERATURE_START 2.0f
#define HK_BLIND_P9016_REP_SWEEP_TEMPERATURE_END 1.0f
#define HK_BLIND_P9016_REP_SWEEP_RHO_TRAIN_START 1.0f
#define HK_BLIND_P9016_REP_SWEEP_RHO_TRAIN_END 1.0f
#define HK_BLIND_P9016_REP_SWEEP_UNIT 1.0f
#define HK_BLIND_P9016_REP_SWEEP_D_SCALE 1.0f
#define HK_BLIND_P9016_REP_SWEEP_BASE_K 2.0f
#define HK_BLIND_P9016_REP_SWEEP_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_REP_SWEEP_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_REP_SWEEP_INIT_EPS 0.5f
#define HK_BLIND_P9016_REP_SWEEP_INIT_NOISE_SCALE 0.0f
#define HK_BLIND_P9016_REP_SWEEP_INIT_SEED 17ULL

#if defined(HK_BLIND_REPULSION_TRADEOFF_SWEEP) || defined(HK_BLIND_REPULSION_DEPTH_SWEEP)
#define HK_BLIND_REPULSION_DETAILED_SWEEP 1
#endif

#ifdef HK_BLIND_REPULSION_DEPTH_SWEEP
#define HK_BLIND_P9016_REP_SWEEP_LABEL "P9016 repulsion depth sweep"
#define HK_BLIND_P9016_REP_SWEEP_ROOT_PREFIX "/tmp/hk_blind_p9016_repulsion_depth_sweep"
#define HK_BLIND_P9016_REP_SWEEP_SUMMARY_NAME "repulsion_depth_sweep_summary.tsv"
#elif defined(HK_BLIND_REPULSION_TRADEOFF_SWEEP)
#define HK_BLIND_P9016_REP_SWEEP_LABEL "P9016 repulsion tradeoff sweep"
#define HK_BLIND_P9016_REP_SWEEP_ROOT_PREFIX "/tmp/hk_blind_p9016_repulsion_tradeoff_sweep"
#define HK_BLIND_P9016_REP_SWEEP_SUMMARY_NAME "repulsion_tradeoff_sweep_summary.tsv"
#else
#define HK_BLIND_P9016_REP_SWEEP_LABEL "P9016 repulsion strength sweep"
#define HK_BLIND_P9016_REP_SWEEP_ROOT_PREFIX "/tmp/hk_blind_p9016_repulsion_strength_sweep"
#define HK_BLIND_P9016_REP_SWEEP_SUMMARY_NAME "repulsion_strength_sweep_summary.tsv"
#endif

struct sweep_config {
	float multiplier;
	int32_t relax_steps;
};

struct bpair_dist_metrics {
	float mean_pU, median_pU;
	float pU_p10, pU_p90;
	float frac_pU_lt_0_90, frac_pU_lt_0_75, frac_pU_lt_0_50;
	float mean_pmax, median_pmax;
	float pmax_p10, pmax_p90;
	float frac_pmax_gt_0_50, frac_pmax_gt_0_75;
	float mean_margin, median_margin;
	float margin_p10, margin_p90;
	float frac_margin_gt_0_10, frac_margin_gt_0_25;
	float mean_psame, mean_pcross;
	float frac_psame_gt_0_75, frac_pcross_gt_0_75;
};

struct sep_dist_metrics {
	float mean, median, min, max;
	float p10, p90;
	float frac_lt_0_25, frac_lt_0_50, frac_lt_1_00;
};

struct sweep_result {
	const struct sweep_config *config;
	int32_t n_raw;
	int32_t n_bpair;
	int32_t n_beads;
	float k_rel_rep;
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
#ifdef HK_BLIND_REPULSION_DEPTH_SWEEP
	{ 8.0f, 20 },
	{ 8.0f, 50 },
	{ 8.0f, 100 },
	{ 10.0f, 20 },
	{ 10.0f, 50 },
	{ 10.0f, 100 },
	{ 12.0f, 20 },
	{ 12.0f, 50 },
	{ 12.0f, 100 }
#elif defined(HK_BLIND_REPULSION_TRADEOFF_SWEEP)
	{ 1.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 2.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 3.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 4.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 6.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 8.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 10.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 12.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS }
#else
	{ 0.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 0.5f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 1.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 2.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 4.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS },
	{ 8.0f, HK_BLIND_P9016_REP_SWEEP_RELAX_STEPS }
#endif
};

static int file_exists(const char *path)
{
	FILE *fp = fopen(path, "rb");
	if (fp == 0) return 0;
	fclose(fp);
	return 1;
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

static int check_close_float(const char *label, float got, float expected)
{
	const float tol = 1e-6f;
	if (fabsf(got - expected) > tol) {
		fprintf(stderr, "%s: got %.9g, expected %.9g\n", label, got, expected);
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
		snprintf(dir, dir_size, "%s_%ld_%ld_%d",
				 HK_BLIND_P9016_REP_SWEEP_ROOT_PREFIX, pid, stamp, i);
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

static float nearest_rank_quantile_sorted_float(const float *a, int32_t n, float q)
{
	double rank;
	int32_t idx;

	assert(a);
	assert(n > 0);
	assert(q >= 0.0f && q <= 1.0f);
	// Nearest-rank quantile: sort values, use ceil(q * n) - 1, then clamp.
	rank = ceil((double)q * (double)n);
	idx = (int32_t)rank - 1;
	if (idx < 0) idx = 0;
	if (idx >= n) idx = n - 1;
	return a[idx];
}

static int test_quantile_helper(void)
{
	float a[5] = { 3.0f, 1.0f, 4.0f, 0.0f, 2.0f };
	qsort(a, 5, sizeof(a[0]), float_cmp);
	return check_close_float("quantile p10", nearest_rank_quantile_sorted_float(a, 5, 0.10f), 0.0f) |
		   check_close_float("quantile p90", nearest_rank_quantile_sorted_float(a, 5, 0.90f), 4.0f);
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
	conf->n_iter = HK_BLIND_P9016_REP_SWEEP_N_ITER;
	conf->base_conf.unit = HK_BLIND_P9016_REP_SWEEP_UNIT;
	conf->base_conf.d_scale = HK_BLIND_P9016_REP_SWEEP_D_SCALE;
	conf->base_conf.base_k = HK_BLIND_P9016_REP_SWEEP_BASE_K;
	conf->base_conf.temperature = HK_BLIND_P9016_REP_SWEEP_TEMPERATURE_START;
	conf->base_conf.rho_train = HK_BLIND_P9016_REP_SWEEP_RHO_TRAIN_START;
	conf->base_conf.min_sep_unit = HK_BLIND_P9016_REP_SWEEP_MIN_SEP_UNIT;
	conf->base_conf.lambda_sep = HK_BLIND_P9016_REP_SWEEP_LAMBDA_SEP;
	conf->base_conf.relax_step = HK_BLIND_P9016_REP_SWEEP_RELAX_STEP;
	conf->base_conf.relax_steps = config->relax_steps;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->base_conf.d_scale_eps_count = 1e-6f;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = HK_BLIND_P9016_REP_SWEEP_TEMPERATURE_START;
	conf->temperature_end = HK_BLIND_P9016_REP_SWEEP_TEMPERATURE_END;
	conf->rho_train_start = HK_BLIND_P9016_REP_SWEEP_RHO_TRAIN_START;
	conf->rho_train_end = HK_BLIND_P9016_REP_SWEEP_RHO_TRAIN_END;
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

	(void)config;
	failed |= check_i32("rep strength completed", diag->n_completed, HK_BLIND_P9016_REP_SWEEP_N_ITER);
	failed |= check_i32("sweep bad iter", diag->n_bad_iter, 0);
	failed |= check_i32("sweep relax nonfinite iter", diag->n_relax_nonfinite_iter, 0);
	failed |= check_i32("sweep coord nonfinite", diag->n_coord_nonfinite, 0);
	failed |= check_true("sweep entropy finite", isfinite(diag->final_mean_entropy));
	failed |= check_true("sweep pU finite", isfinite(diag->final_mean_pU));
	failed |= check_true("sweep pU range", diag->final_mean_pU >= -1e-6f && diag->final_mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("sweep sum wedge finite", isfinite(diag->final_sum_wedge_k));
	failed |= check_true("sweep sum wedge nonnegative", diag->final_sum_wedge_k >= 0.0);
	last_pre = &per_iter[HK_BLIND_P9016_REP_SWEEP_N_ITER - 1].pre_relax_diag;
	failed |= check_i32("sweep final posterior nonfinite diag", last_pre->n_posterior_nonfinite, 0);
	failed |= check_i32("sweep final posterior bad sum diag", last_pre->n_posterior_bad_sum, 0);
	failed |= check_i32("sweep final posterior out range diag", last_pre->n_posterior_out_of_range, 0);
	failed |= check_i32("sweep final uncertainty nonfinite diag", last_pre->n_uncertainty_nonfinite, 0);
	failed |= check_i32("sweep final uncertainty out range diag", last_pre->n_uncertainty_out_of_range, 0);
	failed |= check_i32("sweep final five-state bad sum diag", last_pre->n_five_state_bad_sum, 0);
	failed |= check_i32("sweep final repulsion nonfinite step",
						per_iter[HK_BLIND_P9016_REP_SWEEP_N_ITER - 1].relax_diag.n_repulsion_nonfinite_step, 0);
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
	out->pU_p10 = nearest_rank_quantile_sorted_float(pU, n, 0.10f);
	out->pU_p90 = nearest_rank_quantile_sorted_float(pU, n, 0.90f);
	out->median_pmax = median_float(pmax, n);
	out->pmax_p10 = nearest_rank_quantile_sorted_float(pmax, n, 0.10f);
	out->pmax_p90 = nearest_rank_quantile_sorted_float(pmax, n, 0.90f);
	out->median_margin = median_float(margin, n);
	out->margin_p10 = nearest_rank_quantile_sorted_float(margin, n, 0.10f);
	out->margin_p90 = nearest_rank_quantile_sorted_float(margin, n, 0.90f);
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
	out->p10 = nearest_rank_quantile_sorted_float(sep, n, 0.10f);
	out->p90 = nearest_rank_quantile_sorted_float(sep, n, 0.90f);
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
	failed |= check_true("sweep pU p10 finite", isfinite(b->pU_p10));
	failed |= check_true("sweep pU p90 finite", isfinite(b->pU_p90));
	failed |= check_true("sweep mean pmax finite", isfinite(b->mean_pmax));
	failed |= check_true("sweep pmax p10 finite", isfinite(b->pmax_p10));
	failed |= check_true("sweep pmax p90 finite", isfinite(b->pmax_p90));
	failed |= check_true("sweep mean margin finite", isfinite(b->mean_margin));
	failed |= check_true("sweep margin p10 finite", isfinite(b->margin_p10));
	failed |= check_true("sweep margin p90 finite", isfinite(b->margin_p90));
	failed |= check_true("sweep mean psame finite", isfinite(b->mean_psame));
	failed |= check_true("sweep mean pcross finite", isfinite(b->mean_pcross));
	failed |= check_true("sweep sep mean finite result", isfinite(s->mean));
	failed |= check_true("sweep sep median finite result", isfinite(s->median));
	failed |= check_true("sweep sep p10 finite result", isfinite(s->p10));
	failed |= check_true("sweep sep p90 finite result", isfinite(s->p90));
	failed |= check_true("sweep sep min finite result", isfinite(s->min));
	failed |= check_true("sweep sep max finite result", isfinite(s->max));
	return failed;
}

static int write_summary_header(FILE *fp)
{
#ifdef HK_BLIND_REPULSION_DETAILED_SWEEP
	return fprintf(fp,
				   "multiplier\tk_rel_rep\tn_raw\tn_bpair\tn_beads\tn_iter\trelax_steps\t"
				   "relax_step\ttemperature_start\ttemperature_end\tfinal_mean_entropy\t"
				   "final_mean_pU\tmean_pU\tmedian_pU\tpU_p10\tpU_p90\t"
				   "frac_pU_lt_0_90\tfrac_pU_lt_0_75\tfrac_pU_lt_0_50\t"
				   "mean_pmax\tmedian_pmax\tpmax_p10\tpmax_p90\tmean_margin\t"
				   "median_margin\tmargin_p10\tmargin_p90\tmean_psame\tmean_pcross\t"
				   "sep_mean\tsep_median\tsep_p10\tsep_p90\tsep_min\tsep_max\t"
				   "frac_sep_lt_0_25\tfrac_sep_lt_0_50\tfrac_sep_lt_1_00\t"
				   "final_contact_energy\tfinal_backbone_energy\tfinal_repulsion_energy\t"
				   "final_sep_energy\tfinal_total_energy\tfinal_force_l1\trepulsion_over_contact\t"
				   "contact_energy_per_wedge_k\tn_bad_iter\tn_relax_nonfinite_iter\t"
				   "n_coord_nonfinite\tstatus\telapsed_cpu_sec\n") < 0? -1 : 0;
#else
	return fprintf(fp,
				   "multiplier\tn_raw\tn_bpair\tn_beads\tn_iter\trelax_steps\trelax_step\t"
				   "temperature_start\ttemperature_end\trho_train_start\trho_train_end\t"
				   "k_rel_rep\tfinal_mean_entropy\tfinal_mean_pU\tmean_pU\tmedian_pU\t"
				   "frac_pU_lt_0_90\tfrac_pU_lt_0_75\tmean_pmax\tmedian_pmax\t"
				   "mean_margin\tmedian_margin\tmean_psame\tmean_pcross\t"
				   "sep_mean\tsep_median\tsep_min\tsep_max\tfrac_sep_lt_0_25\t"
				   "frac_sep_lt_0_50\tfrac_sep_lt_1_00\tfinal_contact_energy\t"
				   "final_backbone_energy\tfinal_repulsion_energy\tfinal_sep_energy\t"
				   "final_total_energy\tfinal_force_l1\tn_bad_iter\tn_relax_nonfinite_iter\t"
				   "n_coord_nonfinite\tstatus\telapsed_cpu_sec\n") < 0? -1 : 0;
#endif
}

static int append_summary_row(FILE *fp, const struct sweep_result *result)
{
	const struct sweep_config *c = result->config;
	const struct hk_blind_iter_loop_diag *d = &result->loop_diag;
	const struct bpair_dist_metrics *b = &result->bpair_metrics;
	const struct sep_dist_metrics *s = &result->sep_metrics;
	const char *status = result->status_ok? "OK" : "FAIL";

#ifdef HK_BLIND_REPULSION_DETAILED_SWEEP
	const double repulsion_over_contact = result->final_contact_energy > 0.0f?
		(double)result->final_repulsion_energy / (double)result->final_contact_energy : 0.0;
	const double contact_energy_per_wedge_k = d->final_sum_wedge_k > 0.0?
		(double)result->final_contact_energy / d->final_sum_wedge_k : 0.0;

	return fprintf(fp,
				   "%.9g\t%.9g\t%d\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%d\t"
				   "%d\t%s\t%.3f\n",
				   c->multiplier, result->k_rel_rep, result->n_raw, result->n_bpair,
				   result->n_beads, HK_BLIND_P9016_REP_SWEEP_N_ITER,
				   c->relax_steps, HK_BLIND_P9016_REP_SWEEP_RELAX_STEP,
				   HK_BLIND_P9016_REP_SWEEP_TEMPERATURE_START,
				   HK_BLIND_P9016_REP_SWEEP_TEMPERATURE_END,
				   d->final_mean_entropy, d->final_mean_pU,
				   b->mean_pU, b->median_pU, b->pU_p10, b->pU_p90,
				   b->frac_pU_lt_0_90, b->frac_pU_lt_0_75, b->frac_pU_lt_0_50,
				   b->mean_pmax, b->median_pmax, b->pmax_p10, b->pmax_p90,
				   b->mean_margin, b->median_margin, b->margin_p10, b->margin_p90,
				   b->mean_psame, b->mean_pcross, s->mean, s->median,
				   s->p10, s->p90, s->min, s->max, s->frac_lt_0_25,
				   s->frac_lt_0_50, s->frac_lt_1_00, result->final_contact_energy,
				   result->final_backbone_energy, result->final_repulsion_energy,
				   result->final_sep_energy, result->final_total_energy,
				   result->final_force_l1, repulsion_over_contact,
				   contact_energy_per_wedge_k, d->n_bad_iter,
				   d->n_relax_nonfinite_iter, d->n_coord_nonfinite, status,
				   result->elapsed_cpu_sec) < 0? -1 : 0;
#else
	return fprintf(fp,
				   "%.9g\t%d\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%d\t"
				   "%d\t%s\t%.3f\n",
				   c->multiplier, result->n_raw, result->n_bpair, result->n_beads,
				   HK_BLIND_P9016_REP_SWEEP_N_ITER, c->relax_steps,
				   HK_BLIND_P9016_REP_SWEEP_RELAX_STEP,
				   HK_BLIND_P9016_REP_SWEEP_TEMPERATURE_START,
				   HK_BLIND_P9016_REP_SWEEP_TEMPERATURE_END,
				   HK_BLIND_P9016_REP_SWEEP_RHO_TRAIN_START,
				   HK_BLIND_P9016_REP_SWEEP_RHO_TRAIN_END,
				   result->k_rel_rep, d->final_mean_entropy, d->final_mean_pU,
				   b->mean_pU, b->median_pU, b->frac_pU_lt_0_90,
				   b->frac_pU_lt_0_75, b->mean_pmax, b->median_pmax,
				   b->mean_margin, b->median_margin, b->mean_psame, b->mean_pcross,
				   s->mean, s->median, s->min, s->max,
				   s->frac_lt_0_25, s->frac_lt_0_50, s->frac_lt_1_00,
				   result->final_contact_energy, result->final_backbone_energy,
				   result->final_repulsion_energy, result->final_sep_energy, result->final_total_energy,
				   result->final_force_l1, d->n_bad_iter, d->n_relax_nonfinite_iter,
				   d->n_coord_nonfinite, status, result->elapsed_cpu_sec) < 0? -1 : 0;
#endif
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
			failed |= check_true("sweep summary header multiplier", strstr(line, "multiplier") != 0);
			failed |= check_true("sweep summary header mean pU", strstr(line, "mean_pU") != 0);
#ifdef HK_BLIND_REPULSION_DETAILED_SWEEP
			failed |= check_true("sweep summary header pU p10", strstr(line, "pU_p10") != 0);
#endif
			failed |= check_true("sweep summary header sep median", strstr(line, "sep_median") != 0);
			failed |= check_true("sweep summary header repulsion", strstr(line, "final_repulsion_energy") != 0);
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
						  const fvec3_t *haploid, FILE *summary_fp,
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

	fprintf(stderr, "%s multiplier=%.4g: n_iter=%d relax_steps=%d relax_step=%.4g\n",
			HK_BLIND_P9016_REP_SWEEP_LABEL,
			config->multiplier, HK_BLIND_P9016_REP_SWEEP_N_ITER,
			config->relax_steps, HK_BLIND_P9016_REP_SWEEP_RELAX_STEP);
	set = hk_blind_bpair_set_build(bmap, n_raw, raw);
	failed |= check_binned_set(bmap, set, n_raw);
	if (set) result->n_bpair = set->n_bpairs;
	if (failed) goto finish;

	diploid = (fvec3_t*)calloc((size_t)n_diploid, sizeof(*diploid));
	per_iter = (struct hk_blind_single_iter_diag*)calloc((size_t)HK_BLIND_P9016_REP_SWEEP_N_ITER, sizeof(*per_iter));
	if (diploid == 0 || per_iter == 0) {
		fprintf(stderr, "failed to allocate sweep coordinate/diagnostic data for multiplier %.4g\n",
				config->multiplier);
		failed = 1;
		goto finish;
	}
	failed |= check_i32("sweep init diploid",
						hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
																  diploid, HK_BLIND_P9016_REP_SWEEP_INIT_EPS,
																  HK_BLIND_P9016_REP_SWEEP_INIT_NOISE_SCALE,
																  HK_BLIND_P9016_REP_SWEEP_INIT_SEED), 0);
	failed |= check_coords_finite(diploid, n_diploid);
	if (failed) goto finish;

	hk_fdg_conf_init(&fdg_conf);
	fdg_conf.k_rel_rep *= config->multiplier;
	result->k_rel_rep = fdg_conf.k_rel_rep;
	set_schedule_conf(&schedule_conf, config);
	ret = hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &fdg_conf, diploid, 0,
											   &schedule_conf, per_iter, &result->loop_diag);
	failed |= check_i32("sweep scheduled ret", ret, 0);
	if (ret == 0) {
		failed |= validate_loop_health(config, &result->loop_diag, per_iter);
		failed |= check_coords_finite(diploid, n_diploid);
		failed |= validate_final_bpair_set(set);
		if (!failed) {
			const struct hk_blind_relax_diag *rd = &per_iter[HK_BLIND_P9016_REP_SWEEP_N_ITER - 1].relax_diag;
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
	result->status_ok = failed == 0;
	status = result->status_ok? "OK" : "FAIL";
	if (append_summary_row(summary_fp, result) != 0)
		failed = 1;
	result->status_ok = failed == 0;
	status = result->status_ok? "OK" : "FAIL";
	fflush(summary_fp);
	fprintf(stderr,
			"multiplier=%.4g k_rel_rep=%.6g n_iter=%d relax_steps=%d relax_step=%.4g mean_pU=%.6g "
			"frac_pU_lt_0.75=%.6g mean_pmax=%.6g mean_margin=%.6g "
			"sep_mean=%.6g repulsion_energy=%.6g status=%s\n",
			config->multiplier, result->k_rel_rep, HK_BLIND_P9016_REP_SWEEP_N_ITER,
			config->relax_steps, HK_BLIND_P9016_REP_SWEEP_RELAX_STEP,
			result->bpair_metrics.mean_pU, result->bpair_metrics.frac_pU_lt_0_75,
			result->bpair_metrics.mean_pmax, result->bpair_metrics.mean_margin,
			result->sep_metrics.mean, result->final_repulsion_energy, status);

	free(diploid);
	free(per_iter);
	hk_blind_bpair_set_destroy(set);
	return failed != 0;
}

static int run_one_config_to_row_file(size_t config_idx, const struct hk_bmap *bmap,
									  const struct hk_blind_pair *raw, int32_t n_raw,
									  const fvec3_t *haploid, const char *root_dir,
									  const char *row_path)
{
	FILE *fp = fopen(row_path, "w");
	struct sweep_result result;
	int ret;

	(void)root_dir;
	if (fp == 0) {
		fprintf(stderr, "failed to open sweep summary row %s\n", row_path);
		return 1;
	}
	ret = run_one_config(&sweep_configs[config_idx], bmap, raw, n_raw, haploid,
						 fp, &result);
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
								FILE *summary_fp)
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
		snprintf(row_paths[i], sizeof(row_paths[i]), "%s/multiplier_%.3g_rs%d.summary_row.tsv",
				 root_dir, sweep_configs[i].multiplier, sweep_configs[i].relax_steps);

	while (next < n_configs || running > 0) {
		while (next < n_configs && running < (size_t)jobs) {
			size_t idx = next++;
			pid_t pid = fork();
			if (pid < 0) {
				fprintf(stderr, "failed to fork sweep multiplier %.4g\n", sweep_configs[idx].multiplier);
				failed = 1;
				continue;
			}
			if (pid == 0) {
				int ret = run_one_config_to_row_file(idx, bmap, raw, n_raw, haploid,
													 root_dir, row_paths[idx]);
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
						fprintf(stderr, "sweep multiplier %.4g failed in child\n", sweep_configs[i].multiplier);
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
			fprintf(stderr, "sweep multiplier %.4g did not complete\n", sweep_configs[i].multiplier);
			failed = 1;
			continue;
		}
		if (append_summary_row_file(summary_fp, row_paths[i]) != 0) {
			fprintf(stderr, "failed to append summary row for multiplier %.4g\n", sweep_configs[i].multiplier);
			failed = 1;
		}
	}
	free(row_paths);
	free(pids);
	free(done);
	return failed != 0;
}

#ifdef HK_BLIND_REPULSION_DEPTH_SWEEP
struct depth_summary_row {
	float multiplier;
	int32_t relax_steps;
	float mean_pU;
	float mean_pmax;
	float mean_margin;
	float repulsion_over_contact;
	float contact_energy_per_wedge_k;
};

static int parse_summary_float(const char *token, float *out)
{
	char *end = 0;

	assert(token);
	assert(out);
	errno = 0;
	*out = strtof(token, &end);
	return errno == 0 && end != token && *end == 0 && isfinite(*out);
}

static int parse_summary_i32(const char *token, int32_t *out)
{
	char *end = 0;
	long v;

	assert(token);
	assert(out);
	errno = 0;
	v = strtol(token, &end, 10);
	if (errno != 0 || end == token || *end != 0)
		return 0;
	if (v < INT32_MIN || v > INT32_MAX)
		return 0;
	*out = (int32_t)v;
	return 1;
}

static int parse_depth_summary_row(char *line, struct depth_summary_row *out)
{
	char *save = 0;
	char *token;
	char status[16] = {0};
	int col = 0;
	int ok = 1;

	assert(line);
	assert(out);
	memset(out, 0, sizeof(*out));
	for (token = strtok_r(line, "\t\r\n", &save); token != 0;
		 token = strtok_r(0, "\t\r\n", &save), ++col) {
		switch (col) {
		case 0:
			ok &= parse_summary_float(token, &out->multiplier);
			break;
		case 6:
			ok &= parse_summary_i32(token, &out->relax_steps);
			break;
		case 12:
			ok &= parse_summary_float(token, &out->mean_pU);
			break;
		case 19:
			ok &= parse_summary_float(token, &out->mean_pmax);
			break;
		case 23:
			ok &= parse_summary_float(token, &out->mean_margin);
			break;
		case 44:
			ok &= parse_summary_float(token, &out->repulsion_over_contact);
			break;
		case 45:
			ok &= parse_summary_float(token, &out->contact_energy_per_wedge_k);
			break;
		case 49:
			strncpy(status, token, sizeof(status) - 1);
			break;
		default:
			break;
		}
	}
	if (!ok || col < 51)
		return 1;
	if (strcmp(status, "OK") != 0)
		return 1;
	return 0;
}

static int print_depth_interpretation_summary(const char *summary_path)
{
	FILE *fp = fopen(summary_path, "r");
	char line[8192];
	int row_idx = -1;
	int failed = 0;
	int has_best_pU = 0, has_best_pmax = 0, has_best_margin = 0;
	int has_low_contact = 0, has_low_ratio = 0;
	float best_pU = INFINITY, best_pmax = -INFINITY, best_margin = -INFINITY;
	float low_contact = INFINITY, low_ratio = INFINITY;
	float best_pU_mult = 0.0f, best_pmax_mult = 0.0f, best_margin_mult = 0.0f;
	float low_contact_mult = 0.0f, low_ratio_mult = 0.0f;
	int32_t best_pU_steps = 0, best_pmax_steps = 0, best_margin_steps = 0;
	int32_t low_contact_steps = 0, low_ratio_steps = 0;

	if (fp == 0) {
		fprintf(stderr, "failed to open depth summary for interpretation: %s\n", summary_path);
		return 1;
	}
	while (fgets(line, sizeof(line), fp)) {
		struct depth_summary_row row;
		if (row_idx < 0) {
			row_idx = 0;
			continue;
		}
		if (parse_depth_summary_row(line, &row) != 0) {
			fprintf(stderr, "failed to parse depth summary data row %d\n", row_idx);
			failed = 1;
			break;
		}
		if (!has_best_pU || row.mean_pU < best_pU) {
			has_best_pU = 1;
			best_pU = row.mean_pU;
			best_pU_mult = row.multiplier;
			best_pU_steps = row.relax_steps;
		}
		if (!has_best_pmax || row.mean_pmax > best_pmax) {
			has_best_pmax = 1;
			best_pmax = row.mean_pmax;
			best_pmax_mult = row.multiplier;
			best_pmax_steps = row.relax_steps;
		}
		if (!has_best_margin || row.mean_margin > best_margin) {
			has_best_margin = 1;
			best_margin = row.mean_margin;
			best_margin_mult = row.multiplier;
			best_margin_steps = row.relax_steps;
		}
		if (row.mean_pU < 0.75f &&
			(!has_low_contact || row.contact_energy_per_wedge_k < low_contact)) {
			has_low_contact = 1;
			low_contact = row.contact_energy_per_wedge_k;
			low_contact_mult = row.multiplier;
			low_contact_steps = row.relax_steps;
		}
		if (row.mean_margin > 0.20f &&
			(!has_low_ratio || row.repulsion_over_contact < low_ratio)) {
			has_low_ratio = 1;
			low_ratio = row.repulsion_over_contact;
			low_ratio_mult = row.multiplier;
			low_ratio_steps = row.relax_steps;
		}
		++row_idx;
	}
	if (ferror(fp))
		failed = 1;
	fclose(fp);
	if (failed)
		return 1;
	if (!has_best_pU || !has_best_pmax || !has_best_margin) {
		fprintf(stderr, "depth interpretation summary has no OK rows\n");
		return 1;
	}
	fprintf(stderr, "%s interpretation: best mean_pU multiplier=%.4g relax_steps=%d mean_pU=%.6g\n",
			HK_BLIND_P9016_REP_SWEEP_LABEL, best_pU_mult, best_pU_steps, best_pU);
	fprintf(stderr, "%s interpretation: best mean_pmax multiplier=%.4g relax_steps=%d mean_pmax=%.6g\n",
			HK_BLIND_P9016_REP_SWEEP_LABEL, best_pmax_mult, best_pmax_steps, best_pmax);
	fprintf(stderr, "%s interpretation: best mean_margin multiplier=%.4g relax_steps=%d mean_margin=%.6g\n",
			HK_BLIND_P9016_REP_SWEEP_LABEL, best_margin_mult, best_margin_steps, best_margin);
	if (has_low_contact)
		fprintf(stderr, "%s interpretation: lowest contact_energy_per_wedge_k with mean_pU<0.75 "
				"multiplier=%.4g relax_steps=%d contact_energy_per_wedge_k=%.6g\n",
				HK_BLIND_P9016_REP_SWEEP_LABEL, low_contact_mult, low_contact_steps, low_contact);
	else
		fprintf(stderr, "%s interpretation: no config reached mean_pU<0.75\n",
				HK_BLIND_P9016_REP_SWEEP_LABEL);
	if (has_low_ratio)
		fprintf(stderr, "%s interpretation: lowest repulsion_over_contact with mean_margin>0.20 "
				"multiplier=%.4g relax_steps=%d repulsion_over_contact=%.6g\n",
				HK_BLIND_P9016_REP_SWEEP_LABEL, low_ratio_mult, low_ratio_steps, low_ratio);
	else
		fprintf(stderr, "%s interpretation: no config reached mean_margin>0.20\n",
				HK_BLIND_P9016_REP_SWEEP_LABEL);
	return 0;
}
#endif

int main(void)
{
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	fvec3_t *haploid = 0;
	char root_dir[512] = {0};
	char summary_path[768];
	FILE *summary_fp = 0;
	int jobs = parse_sweep_jobs(sizeof(sweep_configs) / sizeof(sweep_configs[0]));
	int failed = 0;
	int32_t n_raw;
	int32_t i;
	size_t c;
	time_t t0 = time(0), t1;
	double elapsed_sec;

	if (jobs < 0)
		return 1;
	if (test_quantile_helper() != 0)
		return 1;
	if (!file_exists(HK_BLIND_P9016_REP_SWEEP_PATH)) {
		fprintf(stderr, "ERROR: %s not found; %s did not start\n",
				HK_BLIND_P9016_REP_SWEEP_PATH, HK_BLIND_P9016_REP_SWEEP_LABEL);
		return 1;
	}

	hk_verbose = 0;
	fprintf(stderr, "%s: reading %s\n", HK_BLIND_P9016_REP_SWEEP_LABEL,
			HK_BLIND_P9016_REP_SWEEP_PATH);
	m = hk_map_read(HK_BLIND_P9016_REP_SWEEP_PATH);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", HK_BLIND_P9016_REP_SWEEP_PATH);
		return 1;
	}
	n_raw = m->n_pairs;
	if (n_raw <= 0) {
		fprintf(stderr, "sweep input raw positive: predicate failed\n");
		failed = 1;
		goto cleanup;
	}

	fprintf(stderr, "%s: building 1Mb bmap from %d raw pairs\n",
			HK_BLIND_P9016_REP_SWEEP_LABEL, n_raw);
	bmap = hk_bmap_gen(m->d, n_raw, m->pairs, HK_BLIND_P9016_REP_SWEEP_RESOLUTION, 1);
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
	path_join(summary_path, sizeof(summary_path), root_dir, HK_BLIND_P9016_REP_SWEEP_SUMMARY_NAME);
	summary_fp = fopen(summary_path, "w");
	if (summary_fp == 0) {
		fprintf(stderr, "failed to open sweep summary %s\n", summary_path);
		failed = 1;
		goto cleanup;
	}
	failed |= check_i32("sweep summary header", write_summary_header(summary_fp), 0);
	if (failed) goto cleanup;
	fflush(summary_fp);

	fprintf(stderr, "%s: output root %s jobs=%d\n",
			HK_BLIND_P9016_REP_SWEEP_LABEL, root_dir, jobs);
	if (jobs == 1) {
		for (c = 0; c < sizeof(sweep_configs) / sizeof(sweep_configs[0]); ++c) {
			struct sweep_result result;
			int ret = run_one_config(&sweep_configs[c], bmap, raw, n_raw, haploid,
									 summary_fp, &result);
			if (ret != 0)
				failed = 1;
		}
	} else {
		failed |= run_configs_parallel(jobs, bmap, raw, n_raw, haploid,
									   root_dir, summary_fp);
	}
	if (fclose(summary_fp) != 0) {
		summary_fp = 0;
		failed = 1;
		goto cleanup;
	}
	summary_fp = 0;
	failed |= validate_summary_file(summary_path, (int)(sizeof(sweep_configs) / sizeof(sweep_configs[0])));
	if (failed) goto cleanup;
#ifdef HK_BLIND_REPULSION_DEPTH_SWEEP
	failed |= print_depth_interpretation_summary(summary_path);
	if (failed) goto cleanup;
#endif

	t1 = time(0);
	elapsed_sec = difftime(t1, t0);
	fprintf(stderr, "%s: output_root=%s elapsed_wall_sec=%.3f status=OK\n",
			HK_BLIND_P9016_REP_SWEEP_LABEL, root_dir, elapsed_sec);

cleanup:
	if (summary_fp) fclose(summary_fp);
	if (failed)
		fprintf(stderr, "%s failed%s%s\n", HK_BLIND_P9016_REP_SWEEP_LABEL,
				root_dir[0]? "; partial output remains in " : "",
				root_dir[0]? root_dir : "");
	free(haploid);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
