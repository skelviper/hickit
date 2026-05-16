#include <assert.h>
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

#define HK_BLIND_P9016_CURVE_PATH "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_CURVE_RESOLUTION 1000000
#define HK_BLIND_P9016_CURVE_N_ITER 3
#define HK_BLIND_P9016_CURVE_RELAX_STEPS 50
#define HK_BLIND_P9016_CURVE_RELAX_STEP 0.0005f
#define HK_BLIND_P9016_CURVE_TEMPERATURE_START 2.0f
#define HK_BLIND_P9016_CURVE_TEMPERATURE_END 1.0f
#define HK_BLIND_P9016_CURVE_RHO_TRAIN_START 1.0f
#define HK_BLIND_P9016_CURVE_RHO_TRAIN_END 1.0f
#define HK_BLIND_P9016_CURVE_UNIT 1.0f
#define HK_BLIND_P9016_CURVE_D_SCALE 1.0f
#define HK_BLIND_P9016_CURVE_BASE_K 2.0f
#define HK_BLIND_P9016_CURVE_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_CURVE_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_CURVE_INIT_EPS 0.5f
#define HK_BLIND_P9016_CURVE_INIT_NOISE_SCALE 0.0f
#define HK_BLIND_P9016_CURVE_INIT_SEED 17ULL

struct candidate_config {
	const char *name;
	float multiplier;
};

struct posterior_metrics {
	float mean_entropy;
	float mean_pU;
	float mean_pmax;
	float mean_margin;
	float mean_psame;
	float mean_pcross;
};

struct sep_metrics {
	float mean;
	float median;
	float min;
	float max;
	float frac_lt_0_25;
	float frac_lt_0_50;
	float frac_lt_1_00;
};

struct candidate_result {
	const struct candidate_config *config;
	char output_dir[768];
	int32_t n_raw;
	int32_t n_bpair;
	int32_t n_beads;
	float k_rel_rep;
	int32_t outer_rows;
	int32_t inner_rows;
	int32_t n_bad_iter;
	int32_t n_relax_nonfinite_iter;
	int32_t n_coord_nonfinite;
	int32_t total_chr_flipped;
	double final_sum_wedge_k;
	struct posterior_metrics posterior;
	struct sep_metrics sep;
	struct hk_blind_step_diag final_step;
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

static int make_output_root(char *dir, size_t dir_size)
{
	long pid = (long)getpid();
	long stamp = (long)time(0);
	int i;

	for (i = 0; i < 100; ++i) {
		snprintf(dir, dir_size, "/tmp/hk_blind_p9016_candidate_curves_%ld_%ld_%d",
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

static int count_nonfinite_coords(const fvec3_t *coords, int32_t n_diploid)
{
	int32_t i;
	int a;
	int n_bad = 0;

	assert(n_diploid >= 0);
	if (n_diploid == 0)
		return 0;
	assert(coords);
	for (i = 0; i < n_diploid; ++i)
		for (a = 0; a < 3; ++a)
			if (!isfinite(coords[i][a]))
				++n_bad;
	return n_bad;
}

static int step_diag_bad(const struct hk_blind_step_diag *diag)
{
	assert(diag);
	return diag->n_contact_nonfinite != 0 ||
		   diag->n_backbone_nonfinite != 0 ||
		   diag->n_repulsion_nonfinite != 0 ||
		   diag->n_force_nonfinite != 0 ||
		   !isfinite(diag->contact_energy) ||
		   !isfinite(diag->backbone_energy) ||
		   !isfinite(diag->repulsion_energy) ||
		   !isfinite(diag->sep_energy) ||
		   !isfinite(diag->total_energy) ||
		   !isfinite(diag->force_l1) ||
		   !isfinite(diag->backbone_force_l1) ||
		   !isfinite(diag->repulsion_force_l1);
}

static int iter_diag_bad(const struct hk_blind_iter_diag *diag)
{
	assert(diag);
	return diag->n_posterior_nonfinite != 0 ||
		   diag->n_posterior_bad_sum != 0 ||
		   diag->n_posterior_out_of_range != 0 ||
		   diag->n_uncertainty_nonfinite != 0 ||
		   diag->n_uncertainty_out_of_range != 0 ||
		   diag->n_five_state_bad_sum != 0 ||
		   diag->n_wedge_nonfinite != 0 ||
		   diag->n_wedge_bad_k != 0 ||
		   diag->n_wedge_bad_d_scale != 0 ||
		   diag->sep_force_nonfinite != 0;
}

static int compute_posterior_metrics(const struct hk_blind_bpair_set *set,
									 struct posterior_metrics *out)
{
	double sum_entropy = 0.0, sum_pU = 0.0, sum_pmax = 0.0, sum_margin = 0.0;
	double sum_psame = 0.0, sum_pcross = 0.0;
	int32_t i;
	int failed = 0;

	assert(set);
	assert(out);
	memset(out, 0, sizeof(*out));
	if (set->n_bpairs <= 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		float psame = bp->p4[HK_BLIND_STATE_00] + bp->p4[HK_BLIND_STATE_11];
		float pcross = bp->p4[HK_BLIND_STATE_01] + bp->p4[HK_BLIND_STATE_10];
		failed |= check_true("curve entropy finite", isfinite(bp->entropy));
		failed |= check_true("curve pU finite", isfinite(bp->pU));
		failed |= check_true("curve pmax finite", isfinite(bp->pmax));
		failed |= check_true("curve margin finite", isfinite(bp->margin));
		failed |= check_true("curve psame finite", isfinite(psame));
		failed |= check_true("curve pcross finite", isfinite(pcross));
		sum_entropy += bp->entropy;
		sum_pU += bp->pU;
		sum_pmax += bp->pmax;
		sum_margin += bp->margin;
		sum_psame += psame;
		sum_pcross += pcross;
	}
	if (failed)
		return failed;
	out->mean_entropy = (float)(sum_entropy / set->n_bpairs);
	out->mean_pU = (float)(sum_pU / set->n_bpairs);
	out->mean_pmax = (float)(sum_pmax / set->n_bpairs);
	out->mean_margin = (float)(sum_margin / set->n_bpairs);
	out->mean_psame = (float)(sum_psame / set->n_bpairs);
	out->mean_pcross = (float)(sum_pcross / set->n_bpairs);
	failed |= check_true("curve mean entropy finite", isfinite(out->mean_entropy));
	failed |= check_true("curve mean pU range", out->mean_pU >= -1e-6f && out->mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("curve mean pmax range", out->mean_pmax >= 0.25f - 1e-6f && out->mean_pmax <= 1.0f + 1e-6f);
	failed |= check_true("curve mean margin range", out->mean_margin >= -1e-6f && out->mean_margin <= 1.0f + 1e-6f);
	failed |= check_true("curve mean psame range", out->mean_psame >= -1e-6f && out->mean_psame <= 1.0f + 1e-6f);
	failed |= check_true("curve mean pcross range", out->mean_pcross >= -1e-6f && out->mean_pcross <= 1.0f + 1e-6f);
	return failed;
}

static int compute_sep_metrics(const struct hk_bmap *bmap, const fvec3_t *coords,
							   struct sep_metrics *out)
{
	float *sep = 0;
	double sum = 0.0;
	int32_t i;
	int failed = 0;

	assert(bmap);
	assert(coords);
	assert(out);
	memset(out, 0, sizeof(*out));
	if (bmap->n_beads <= 0)
		return 1;
	sep = (float*)malloc((size_t)bmap->n_beads * sizeof(*sep));
	if (sep == 0)
		return 1;
	out->min = INFINITY;
	for (i = 0; i < bmap->n_beads; ++i) {
		int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		sep[i] = dist3(coords[b0], coords[b1]);
		failed |= check_true("curve sep finite", isfinite(sep[i]));
		if (sep[i] < out->min) out->min = sep[i];
		if (sep[i] > out->max) out->max = sep[i];
		if (sep[i] < 0.25f) out->frac_lt_0_25 += 1.0f;
		if (sep[i] < 0.50f) out->frac_lt_0_50 += 1.0f;
		if (sep[i] < 1.00f) out->frac_lt_1_00 += 1.0f;
		sum += sep[i];
	}
	if (failed) {
		free(sep);
		return failed;
	}
	out->mean = (float)(sum / bmap->n_beads);
	out->median = median_float(sep, bmap->n_beads);
	out->frac_lt_0_25 /= (float)bmap->n_beads;
	out->frac_lt_0_50 /= (float)bmap->n_beads;
	out->frac_lt_1_00 /= (float)bmap->n_beads;
	failed |= check_true("curve sep mean finite", isfinite(out->mean));
	failed |= check_true("curve sep median finite", isfinite(out->median));
	failed |= check_true("curve sep min finite", isfinite(out->min));
	failed |= check_true("curve sep max finite", isfinite(out->max));
	free(sep);
	return failed;
}

static int validate_binned_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
							   int32_t n_raw)
{
	int failed = 0;
	int32_t i;

	failed |= check_true("curve set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_i32("curve raw count", set->n_raw, n_raw);
	failed |= check_true("curve bpair positive", set->n_bpairs > 0);
	for (i = 0; i < set->n_raw; ++i) {
		failed |= check_true("curve raw2binned id lower", set->raw2binned[i].bpair_id >= 0);
		failed |= check_true("curve raw2binned id upper", set->raw2binned[i].bpair_id < set->n_bpairs);
		failed |= check_true("curve raw2binned swapped",
							 set->raw2binned[i].swapped == 0 || set->raw2binned[i].swapped == 1);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		failed |= check_true("curve key sorted", bp->key.bid[0] <= bp->key.bid[1]);
		failed |= check_true("curve key bid0 valid", bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		failed |= check_true("curve key bid1 valid", bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
	}
	return failed;
}

static void set_schedule_conf(struct hk_blind_iter_schedule_conf *conf)
{
	memset(conf, 0, sizeof(*conf));
	conf->n_iter = HK_BLIND_P9016_CURVE_N_ITER;
	conf->base_conf.unit = HK_BLIND_P9016_CURVE_UNIT;
	conf->base_conf.d_scale = HK_BLIND_P9016_CURVE_D_SCALE;
	conf->base_conf.base_k = HK_BLIND_P9016_CURVE_BASE_K;
	conf->base_conf.temperature = HK_BLIND_P9016_CURVE_TEMPERATURE_START;
	conf->base_conf.rho_train = HK_BLIND_P9016_CURVE_RHO_TRAIN_START;
	conf->base_conf.min_sep_unit = HK_BLIND_P9016_CURVE_MIN_SEP_UNIT;
	conf->base_conf.lambda_sep = HK_BLIND_P9016_CURVE_LAMBDA_SEP;
	conf->base_conf.relax_step = HK_BLIND_P9016_CURVE_RELAX_STEP;
	conf->base_conf.relax_steps = HK_BLIND_P9016_CURVE_RELAX_STEPS;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->base_conf.d_scale_eps_count = 1e-6f;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = HK_BLIND_P9016_CURVE_TEMPERATURE_START;
	conf->temperature_end = HK_BLIND_P9016_CURVE_TEMPERATURE_END;
	conf->rho_train_start = HK_BLIND_P9016_CURVE_RHO_TRAIN_START;
	conf->rho_train_end = HK_BLIND_P9016_CURVE_RHO_TRAIN_END;
}

static int write_outer_header(FILE *fp)
{
	return fprintf(fp,
				   "candidate\titer\ttemperature\trho_train\trel_rep_multiplier\t"
				   "mean_entropy\tmean_pU\tmean_pmax\tmean_margin\tmean_psame\tmean_pcross\t"
				   "sep_mean\tsep_median\tsep_min\tsep_max\tfrac_sep_lt_0_25\t"
				   "frac_sep_lt_0_50\tfrac_sep_lt_1_00\tcontact_energy\tbackbone_energy\t"
				   "repulsion_energy\tsep_energy\ttotal_energy\tforce_l1\tn_chr_flipped\t"
				   "n_bad_iter\tn_relax_nonfinite_iter\tn_coord_nonfinite\n") < 0? -1 : 0;
}

static int write_inner_header(FILE *fp)
{
	return fprintf(fp,
				   "candidate\titer\trelax_step_index\trel_rep_k\tcontact_energy\t"
				   "backbone_energy\trepulsion_energy\tsep_energy\ttotal_energy\tforce_l1\t"
				   "repulsion_force_l1\tbackbone_force_l1\tsep_mean\tsep_median\t"
				   "frac_sep_lt_0_25\tfrac_sep_lt_0_50\tfrac_sep_lt_1_00\t"
				   "n_contact_nonfinite\tn_backbone_nonfinite\tn_repulsion_nonfinite\t"
				   "n_coord_nonfinite\n") < 0? -1 : 0;
}

static int write_final_header(FILE *fp)
{
	return fprintf(fp,
				   "candidate\toutput_dir\tmultiplier\tk_rel_rep\tn_raw\tn_bpair\tn_beads\t"
				   "n_iter\trelax_steps\trelax_step\touter_rows\tinner_rows\t"
				   "final_mean_entropy\tfinal_mean_pU\tfinal_mean_pmax\tfinal_mean_margin\t"
				   "final_mean_psame\tfinal_mean_pcross\tfinal_sep_mean\tfinal_sep_median\t"
				   "final_sep_min\tfinal_sep_max\tfinal_contact_energy\tfinal_backbone_energy\t"
				   "final_repulsion_energy\tfinal_sep_energy\tfinal_total_energy\tfinal_force_l1\t"
				   "final_sum_wedge_k\ttotal_chr_flipped\tn_bad_iter\tn_relax_nonfinite_iter\t"
				   "n_coord_nonfinite\tstatus\telapsed_cpu_sec\n") < 0? -1 : 0;
}

static int append_outer_row(FILE *fp, const struct candidate_config *config, int32_t iter,
							float temperature, float rho_train,
							const struct posterior_metrics *posterior,
							const struct sep_metrics *sep,
							const struct hk_blind_step_diag *step_diag,
							int32_t n_chr_flipped, int32_t n_bad_iter,
							int32_t n_relax_nonfinite_iter, int32_t n_coord_nonfinite)
{
	return fprintf(fp,
				   "%s\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%d\t%d\t%d\t%d\n",
				   config->name, iter, temperature, rho_train, config->multiplier,
				   posterior->mean_entropy, posterior->mean_pU, posterior->mean_pmax,
				   posterior->mean_margin, posterior->mean_psame, posterior->mean_pcross,
				   sep->mean, sep->median, sep->min, sep->max, sep->frac_lt_0_25,
				   sep->frac_lt_0_50, sep->frac_lt_1_00, step_diag->contact_energy,
				   step_diag->backbone_energy, step_diag->repulsion_energy, step_diag->sep_energy,
				   step_diag->total_energy, step_diag->force_l1, n_chr_flipped,
				   n_bad_iter, n_relax_nonfinite_iter, n_coord_nonfinite) < 0? -1 : 0;
}

static int append_inner_row(FILE *fp, const struct candidate_config *config, int32_t iter,
							int32_t step_index, float rel_rep_k,
							const struct hk_blind_step_diag *step_diag,
							const struct sep_metrics *sep, int32_t n_coord_nonfinite)
{
	return fprintf(fp,
				   "%s\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%d\t%d\t%d\n",
				   config->name, iter, step_index, rel_rep_k,
				   step_diag->contact_energy, step_diag->backbone_energy,
				   step_diag->repulsion_energy, step_diag->sep_energy,
				   step_diag->total_energy, step_diag->force_l1,
				   step_diag->repulsion_force_l1, step_diag->backbone_force_l1,
				   sep->mean, sep->median, sep->frac_lt_0_25, sep->frac_lt_0_50,
				   sep->frac_lt_1_00, step_diag->n_contact_nonfinite,
				   step_diag->n_backbone_nonfinite, step_diag->n_repulsion_nonfinite,
				   n_coord_nonfinite) < 0? -1 : 0;
}

static int append_final_row(FILE *fp, const struct candidate_result *result)
{
	const char *status = result->status_ok? "OK" : "FAIL";
	return fprintf(fp,
				   "%s\t%s\t%.9g\t%.9g\t%d\t%d\t%d\t%d\t%d\t%.9g\t%d\t%d\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.17g\t%d\t%d\t%d\t%d\t"
				   "%s\t%.3f\n",
				   result->config->name, result->output_dir, result->config->multiplier,
				   result->k_rel_rep, result->n_raw, result->n_bpair, result->n_beads,
				   HK_BLIND_P9016_CURVE_N_ITER, HK_BLIND_P9016_CURVE_RELAX_STEPS,
				   HK_BLIND_P9016_CURVE_RELAX_STEP, result->outer_rows, result->inner_rows,
				   result->posterior.mean_entropy, result->posterior.mean_pU,
				   result->posterior.mean_pmax, result->posterior.mean_margin,
				   result->posterior.mean_psame, result->posterior.mean_pcross,
				   result->sep.mean, result->sep.median, result->sep.min, result->sep.max,
				   result->final_step.contact_energy, result->final_step.backbone_energy,
				   result->final_step.repulsion_energy, result->final_step.sep_energy,
				   result->final_step.total_energy, result->final_step.force_l1,
				   result->final_sum_wedge_k, result->total_chr_flipped, result->n_bad_iter,
				   result->n_relax_nonfinite_iter, result->n_coord_nonfinite, status,
				   result->elapsed_cpu_sec) < 0? -1 : 0;
}

static double wedge_sum_k(const struct hk_blind_wedge_list *wedges)
{
	double sum = 0.0;
	int32_t i;
	assert(wedges);
	for (i = 0; i < wedges->n_edges; ++i)
		sum += wedges->edges[i].k;
	return sum;
}

static int run_one_candidate(const struct candidate_config *config, const struct hk_bmap *bmap,
							 const struct hk_blind_pair *raw, int32_t n_raw,
							 const fvec3_t *haploid, const char *root_dir,
							 struct candidate_result *result)
{
	struct hk_blind_bpair_set *set = 0;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_wedge_list wedges;
	fvec3_t *diploid = 0;
	float log_prior[HK_BLIND_N_STATE];
	char outer_path[1024], inner_path[1024], final_path[1024];
	FILE *outer_fp = 0, *inner_fp = 0, *final_fp = 0;
	int32_t n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	int failed = 0;
	clock_t t0 = clock(), t1;
	int32_t iter;

	memset(result, 0, sizeof(*result));
	result->config = config;
	hk_blind_step_diag_init(&result->final_step);
	hk_blind_wedge_list_init(&wedges);
	path_join(result->output_dir, sizeof(result->output_dir), root_dir, config->name);
	if (mkdir(result->output_dir, 0700) != 0) {
		fprintf(stderr, "failed to create curve directory %s\n", result->output_dir);
		failed = 1;
		goto cleanup;
	}
	path_join(outer_path, sizeof(outer_path), result->output_dir, "outer_iter_curve.tsv");
	path_join(inner_path, sizeof(inner_path), result->output_dir, "inner_relax_curve.tsv");
	path_join(final_path, sizeof(final_path), result->output_dir, "final_summary.tsv");

	outer_fp = fopen(outer_path, "w");
	inner_fp = fopen(inner_path, "w");
	final_fp = fopen(final_path, "w");
	if (outer_fp == 0 || inner_fp == 0 || final_fp == 0) {
		fprintf(stderr, "failed to open curve outputs in %s\n", result->output_dir);
		failed = 1;
		goto cleanup;
	}
	failed |= check_i32("curve outer header", write_outer_header(outer_fp), 0);
	failed |= check_i32("curve inner header", write_inner_header(inner_fp), 0);
	failed |= check_i32("curve final header", write_final_header(final_fp), 0);
	if (failed) goto cleanup;

	fprintf(stderr, "P9016 candidate curves %s: multiplier=%.4g relax_steps=%d\n",
			config->name, config->multiplier, HK_BLIND_P9016_CURVE_RELAX_STEPS);
	set = hk_blind_bpair_set_build(bmap, n_raw, raw);
	failed |= validate_binned_set(bmap, set, n_raw);
	if (set) {
		result->n_raw = set->n_raw;
		result->n_bpair = set->n_bpairs;
	}
	result->n_beads = bmap->n_beads;
	if (failed) goto cleanup;

	diploid = (fvec3_t*)calloc((size_t)n_diploid, sizeof(*diploid));
	if (diploid == 0) {
		fprintf(stderr, "failed to allocate curve coordinates for %s\n", config->name);
		failed = 1;
		goto cleanup;
	}
	failed |= check_i32("curve init diploid",
						hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
																  diploid, HK_BLIND_P9016_CURVE_INIT_EPS,
																  HK_BLIND_P9016_CURVE_INIT_NOISE_SCALE,
																  HK_BLIND_P9016_CURVE_INIT_SEED), 0);
	failed |= check_i32("curve initial coord finite", count_nonfinite_coords(diploid, n_diploid), 0);
	if (failed) goto cleanup;

	hk_fdg_conf_init(&fdg_conf);
	fdg_conf.k_rel_rep *= config->multiplier;
	result->k_rel_rep = fdg_conf.k_rel_rep;
	set_schedule_conf(&schedule_conf);
	hk_blind_init_uniform_log_prior(log_prior);

	for (iter = 0; iter < HK_BLIND_P9016_CURVE_N_ITER; ++iter) {
		struct hk_blind_iter_diag iter_diag;
		struct hk_blind_gauge_stats gauge_stats;
		struct posterior_metrics posterior;
		struct sep_metrics sep_after_step, sep_after_gauge;
		struct hk_blind_step_diag last_step_diag;
		fvec3_t *prev_coords = 0;
		float temperature = hk_blind_iter_schedule_temperature_at(&schedule_conf, iter);
		float rho_train = hk_blind_iter_schedule_rho_train_at(&schedule_conf, iter);
		int32_t iter_bad = 0;
		int32_t relax_bad = 0;
		int32_t coord_bad = 0;
		int32_t step;
		int ret;

		hk_blind_step_diag_init(&last_step_diag);
		prev_coords = (fvec3_t*)malloc((size_t)n_diploid * sizeof(*prev_coords));
		if (prev_coords == 0) {
			fprintf(stderr, "failed to allocate previous coords for %s iter %d\n",
					config->name, iter);
			failed = 1;
			break;
		}
		memcpy(prev_coords, diploid, (size_t)n_diploid * sizeof(*diploid));

		hk_blind_bpair_set_update_posterior_from_coords_params(set, &fdg_conf, diploid,
															   HK_BLIND_P9016_CURVE_UNIT,
															   log_prior, temperature);
		failed |= compute_posterior_metrics(set, &posterior);
		hk_blind_iter_diag_init(&iter_diag);
		iter_diag.n_bpair = set->n_bpairs;
		iter_diag.n_raw = set->n_raw;
		hk_blind_iter_diag_validate_bpair_set(set, &iter_diag);

		ret = hk_blind_wedge_list_build_from_bpair_set_params(&wedges, set, rho_train);
		failed |= check_i32("curve wedge build", ret, 0);
		if (ret == 0) {
			iter_diag.n_wedges_before_aggregation = wedges.n_edges_before_aggregation;
			iter_diag.n_skipped_self_edges = wedges.n_skipped_self_edges;
			ret = hk_blind_wedge_list_aggregate_exact(&wedges);
			failed |= check_i32("curve wedge aggregate", ret, 0);
			if (ret == 0) {
				iter_diag.n_wedges = wedges.n_edges;
				hk_blind_iter_diag_validate_wedge_list(&wedges, &iter_diag);
				result->final_sum_wedge_k = wedge_sum_k(&wedges);
			}
		}
		iter_bad = iter_diag_bad(&iter_diag);
		if (iter_bad)
			++result->n_bad_iter;
		if (failed || iter_bad) {
			free(prev_coords);
			break;
		}

		for (step = 0; step < HK_BLIND_P9016_CURVE_RELAX_STEPS; ++step) {
			struct hk_blind_step_diag step_diag;
			float rel_rep_k = hk_blind_rel_rep_schedule_at(step, HK_BLIND_P9016_CURVE_RELAX_STEPS);

			ret = hk_blind_relax_step_cpu(&fdg_conf, &wedges, bmap, bmap->n_beads, diploid,
										  HK_BLIND_P9016_CURVE_UNIT,
										  HK_BLIND_P9016_CURVE_RELAX_STEP,
										  HK_BLIND_P9016_CURVE_MIN_SEP_UNIT,
										  HK_BLIND_P9016_CURVE_LAMBDA_SEP,
										  1, HK_BLIND_REPULSION_CELL, rel_rep_k, &step_diag);
			coord_bad = count_nonfinite_coords(diploid, n_diploid);
			failed |= compute_sep_metrics(bmap, diploid, &sep_after_step);
			failed |= check_i32("curve relax step ret", ret, 0);
			failed |= check_i32("curve inner write",
								append_inner_row(inner_fp, config, iter, step, rel_rep_k,
												 &step_diag, &sep_after_step, coord_bad), 0);
			++result->inner_rows;
			last_step_diag = step_diag;
			if (ret != 0 || step_diag_bad(&step_diag) || coord_bad != 0) {
				relax_bad = 1;
				break;
			}
			if (failed)
				break;
		}
		if (relax_bad)
			++result->n_relax_nonfinite_iter;
		if (failed || relax_bad) {
			free(prev_coords);
			break;
		}

		hk_blind_gauge_stats_init(&gauge_stats);
		ret = hk_blind_temporal_gauge_stabilize_bmap(bmap, prev_coords, diploid, 0, &gauge_stats);
		free(prev_coords);
		failed |= check_i32("curve gauge stabilize", ret, 0);
		coord_bad = count_nonfinite_coords(diploid, n_diploid);
		result->n_coord_nonfinite = coord_bad;
		result->total_chr_flipped += gauge_stats.n_flipped;
		failed |= check_i32("curve post-gauge coord finite", coord_bad, 0);
		failed |= compute_sep_metrics(bmap, diploid, &sep_after_gauge);
		failed |= check_i32("curve outer write",
							append_outer_row(outer_fp, config, iter, temperature, rho_train,
											 &posterior, &sep_after_gauge, &last_step_diag,
											 gauge_stats.n_flipped, iter_bad, relax_bad,
											 coord_bad), 0);
		++result->outer_rows;
		result->posterior = posterior;
		result->sep = sep_after_gauge;
		result->final_step = last_step_diag;
		if (failed)
			break;
	}

	failed |= check_i32("curve outer row count", result->outer_rows, HK_BLIND_P9016_CURVE_N_ITER);
	failed |= check_i32("curve inner row count", result->inner_rows,
						HK_BLIND_P9016_CURVE_N_ITER * HK_BLIND_P9016_CURVE_RELAX_STEPS);
	failed |= check_i32("curve bad iter count", result->n_bad_iter, 0);
	failed |= check_i32("curve relax bad iter count", result->n_relax_nonfinite_iter, 0);
	failed |= check_i32("curve coord bad count", result->n_coord_nonfinite, 0);
	failed |= check_true("curve final mean pU finite", isfinite(result->posterior.mean_pU));
	failed |= check_true("curve final mean margin finite", isfinite(result->posterior.mean_margin));
	failed |= check_true("curve final force finite", isfinite(result->final_step.force_l1));
	failed |= check_true("curve final energy finite", isfinite(result->final_step.total_energy));

cleanup:
	t1 = clock();
	result->elapsed_cpu_sec = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
	if (outer_fp) {
		if (fclose(outer_fp) != 0) failed = 1;
		outer_fp = 0;
	}
	if (inner_fp) {
		if (fclose(inner_fp) != 0) failed = 1;
		inner_fp = 0;
	}
	result->status_ok = failed == 0 &&
		result->outer_rows == HK_BLIND_P9016_CURVE_N_ITER &&
		result->inner_rows == HK_BLIND_P9016_CURVE_N_ITER * HK_BLIND_P9016_CURVE_RELAX_STEPS &&
		result->n_bad_iter == 0 &&
		result->n_relax_nonfinite_iter == 0 &&
		result->n_coord_nonfinite == 0;
	if (final_fp) {
		if (append_final_row(final_fp, result) != 0)
			failed = 1;
		fflush(final_fp);
		if (fclose(final_fp) != 0)
			failed = 1;
		final_fp = 0;
	}
	if (failed)
		result->status_ok = 0;
	fprintf(stderr,
			"P9016 candidate curve summary: name=%s final_mean_pU=%.6g "
			"final_mean_margin=%.6g final_force_l1=%.6g final_total_energy=%.6g "
			"output_dir=%s status=%s\n",
			config->name, result->posterior.mean_pU, result->posterior.mean_margin,
			result->final_step.force_l1, result->final_step.total_energy,
			result->output_dir, result->status_ok? "OK" : "FAIL");
	if (failed)
		fprintf(stderr, "P9016 candidate curves %s failed; partial output remains in %s\n",
				config->name, result->output_dir);
	free(diploid);
	hk_blind_wedge_list_destroy(&wedges);
	hk_blind_bpair_set_destroy(set);
	return failed != 0 || !result->status_ok;
}

int main(void)
{
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	fvec3_t *haploid = 0;
	struct candidate_result results[sizeof(candidate_configs) / sizeof(candidate_configs[0])];
	char root_dir[512] = {0};
	int32_t n_raw;
	int32_t i;
	size_t c;
	int failed = 0;
	clock_t t0 = clock(), t1;
	double elapsed_sec;

	if (!file_exists(HK_BLIND_P9016_CURVE_PATH)) {
		fprintf(stderr, "ERROR: %s not found; candidate curve run did not start\n",
				HK_BLIND_P9016_CURVE_PATH);
		return 1;
	}

	hk_verbose = 0;
	fprintf(stderr, "P9016 candidate curves: reading %s\n", HK_BLIND_P9016_CURVE_PATH);
	m = hk_map_read(HK_BLIND_P9016_CURVE_PATH);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", HK_BLIND_P9016_CURVE_PATH);
		return 1;
	}
	n_raw = m->n_pairs;
	if (n_raw <= 0) {
		fprintf(stderr, "curve input raw positive: predicate failed\n");
		failed = 1;
		goto cleanup;
	}

	fprintf(stderr, "P9016 candidate curves: building 1Mb bmap from %d raw pairs\n", n_raw);
	bmap = hk_bmap_gen(m->d, n_raw, m->pairs, HK_BLIND_P9016_CURVE_RESOLUTION, 1);
	raw = (struct hk_blind_pair*)calloc((size_t)n_raw, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate curve bmap/raw data\n");
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n_raw; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);
	failed |= check_true("curve bead count positive", bmap->n_beads > 0);
	if (failed) goto cleanup;

	haploid = (fvec3_t*)calloc((size_t)bmap->n_beads, sizeof(*haploid));
	if (haploid == 0) {
		fprintf(stderr, "failed to allocate curve haploid scaffold\n");
		failed = 1;
		goto cleanup;
	}
	set_haploid_scaffold(bmap, haploid);

	failed |= check_i32("curve output root", make_output_root(root_dir, sizeof(root_dir)), 0);
	if (failed) goto cleanup;
	fprintf(stderr, "P9016 candidate curves: output root %s\n", root_dir);

	memset(results, 0, sizeof(results));
	for (c = 0; c < sizeof(candidate_configs) / sizeof(candidate_configs[0]); ++c) {
		failed |= run_one_candidate(&candidate_configs[c], bmap, raw, n_raw,
									haploid, root_dir, &results[c]);
	}
	if (failed) goto cleanup;

	t1 = clock();
	elapsed_sec = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
	fprintf(stderr, "P9016 candidate curves: output_root=%s elapsed_cpu_sec=%.3f status=OK\n",
			root_dir, elapsed_sec);
	for (c = 0; c < sizeof(candidate_configs) / sizeof(candidate_configs[0]); ++c) {
		fprintf(stderr,
				"P9016 candidate curve final: name=%s mean_pU=%.6g mean_margin=%.6g "
				"force_l1=%.6g total_energy=%.6g output_dir=%s\n",
				results[c].config->name, results[c].posterior.mean_pU,
				results[c].posterior.mean_margin, results[c].final_step.force_l1,
				results[c].final_step.total_energy, results[c].output_dir);
	}

cleanup:
	if (failed)
		fprintf(stderr, "P9016 candidate curves failed%s%s\n",
				root_dir[0]? "; partial output remains in " : "",
				root_dir[0]? root_dir : "");
	free(haploid);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
