#include <assert.h>
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
#include <zlib.h>
#include "hickit.h"

#define HK_BLIND_P9016_SCAN_PATH "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_SCAN_RESOLUTION 1000000
#define HK_BLIND_P9016_SCAN_RELAX_STEP 0.0005f
#define HK_BLIND_P9016_SCAN_TEMPERATURE_START 2.0f
#define HK_BLIND_P9016_SCAN_TEMPERATURE_END 1.0f
#define HK_BLIND_P9016_SCAN_RHO_TRAIN_START 1.0f
#define HK_BLIND_P9016_SCAN_RHO_TRAIN_END 1.0f
#define HK_BLIND_P9016_SCAN_UNIT 1.0f
#define HK_BLIND_P9016_SCAN_D_SCALE 1.0f
#define HK_BLIND_P9016_SCAN_BASE_K 2.0f
#define HK_BLIND_P9016_SCAN_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_SCAN_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_SCAN_INIT_EPS 0.5f
#define HK_BLIND_P9016_SCAN_INIT_NOISE_SCALE 0.0f
#define HK_BLIND_P9016_SCAN_RANDOM_SCALE 1.0f
#define HK_BLIND_P9016_SCAN_INIT_SEED 17ULL
#define HK_BLIND_P9016_SCAN_PRIOR_EPS 1e-6f
#define HK_BLIND_P9016_SCAN_SCAFFOLD_FDG_N_ITER 50
#define HK_BLIND_P9016_SCAN_ROOT_PREFIX "/tmp/hk_blind_p9016_rep1_fine_n_rs_grid"
#define HK_BLIND_P9016_SCAN_SUMMARY_NAME "scan_summary.tsv"
#define HK_BLIND_P9016_SCAN_COORDS_NAME "p9016_final.coords.tsv.gz"
#define HK_BLIND_P9016_SCAN_N_CONFIG 128

struct scan_config {
	char name[64];
	char mode_name[16];
	float multiplier;
	int32_t n_iter;
	int32_t relax_steps;
	int init_mode;
	int prior_mode;
	int rho_train_mode;
	int d_scale_mode;
	float init_eps;
	float init_noise_scale;
	float init_scale;
	uint64_t init_seed;
};

struct scan_mode_def {
	const char *name;
	int init_mode;
	int prior_mode;
	int rho_train_mode;
	int d_scale_mode;
	float init_eps;
	float init_noise_scale;
	float init_scale;
	uint64_t init_seed;
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

struct scan_result {
	const struct scan_config *config;
	char output_dir[1024];
	char coords_gz[1024];
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
	float final_repulsion_force_l1;
	float final_backbone_force_l1;
	int32_t n_repulsion_nonfinite_step;
	int same_bin_filter_enabled;
	int64_t n_raw_same_bin_excluded;
	int32_t n_bpair_same_bin_excluded;
	double repulsion_over_contact;
	double contact_energy_per_wedge_k;
	int status_ok;
	double elapsed_cpu_sec;
};

struct parsed_summary_row {
	char config_name[128];
	char coords_gz[1024];
	float mean_pU;
	float mean_pmax;
	float mean_margin;
	float sep_mean;
	float frac_sep_lt_0_25;
	double repulsion_over_contact;
	double contact_energy_per_wedge_k;
	int status_ok;
};

static int parse_summary_row(char *line, struct parsed_summary_row *out);

static const struct scan_mode_def scan_modes[] = {
	{
		"main",
		HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT,
		HK_BLIND_PRIOR_CIS_INTER_RATIO,
		HK_BLIND_RHO_TRAIN_ENTROPY,
		HK_BLIND_D_SCALE_RAW_COUNT,
		HK_BLIND_P9016_SCAN_INIT_EPS,
		HK_BLIND_P9016_SCAN_INIT_NOISE_SCALE,
		0.0f,
		HK_BLIND_P9016_SCAN_INIT_SEED
	},
	{
		"random",
		HK_BLIND_INIT_RANDOM_DIPLOID,
		HK_BLIND_PRIOR_UNIFORM,
		HK_BLIND_RHO_TRAIN_CONSTANT,
		HK_BLIND_D_SCALE_RAW_COUNT,
		0.0f,
		0.0f,
		HK_BLIND_P9016_SCAN_RANDOM_SCALE,
		HK_BLIND_P9016_SCAN_INIT_SEED
	}
};
static const float scan_multipliers[] = { 1.0f };
static const int32_t scan_n_iters[] = { 3, 5, 10, 20, 50, 80, 100, 120 };
static const int32_t scan_relax_steps[] = { 3, 5, 10, 25, 50, 80, 100, 120 };
static struct scan_config scan_configs[HK_BLIND_P9016_SCAN_N_CONFIG];
static size_t n_scan_configs = 0;

static int file_exists(const char *path)
{
	FILE *fp = fopen(path, "rb");
	if (fp == 0) return 0;
	fclose(fp);
	return 1;
}

static int parse_env_flag(const char *name)
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

static int validate_output_root_override(const char *dir)
{
	const char *prefix = HK_BLIND_P9016_SCAN_ROOT_PREFIX;
	size_t prefix_len = strlen(prefix);

	if (dir == 0 || *dir == 0)
		return 0;
	return strncmp(dir, prefix, prefix_len) == 0;
}

static int make_output_root(char *dir, size_t dir_size)
{
	const char *root_override = getenv("HK_BLIND_SCAN_ROOT");
	long pid = (long)getpid();
	long stamp = (long)time(0);
	int i;

	if (root_override != 0 && *root_override != 0) {
		if (!validate_output_root_override(root_override)) {
			fprintf(stderr, "ERROR: HK_BLIND_SCAN_ROOT must start with %s\n",
					HK_BLIND_P9016_SCAN_ROOT_PREFIX);
			return -1;
		}
		snprintf(dir, dir_size, "%s", root_override);
		if (mkdir(dir, 0700) == 0 || errno == EEXIST)
			return 0;
		return -1;
	}
	for (i = 0; i < 100; ++i) {
		snprintf(dir, dir_size, "%s_%ld_%ld_%d",
				 HK_BLIND_P9016_SCAN_ROOT_PREFIX, pid, stamp, i);
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

static int make_scan_configs(void)
{
	size_t idx = 0;
	size_t h, i, j, k;

	for (h = 0; h < sizeof(scan_modes) / sizeof(scan_modes[0]); ++h) {
		for (i = 0; i < sizeof(scan_multipliers) / sizeof(scan_multipliers[0]); ++i) {
			for (j = 0; j < sizeof(scan_n_iters) / sizeof(scan_n_iters[0]); ++j) {
				for (k = 0; k < sizeof(scan_relax_steps) / sizeof(scan_relax_steps[0]); ++k, ++idx) {
					struct scan_config *c;
					const struct scan_mode_def *mode = &scan_modes[h];
					if (idx >= HK_BLIND_P9016_SCAN_N_CONFIG)
						return 1;
					c = &scan_configs[idx];
					snprintf(c->mode_name, sizeof(c->mode_name), "%s", mode->name);
					c->multiplier = scan_multipliers[i];
					c->n_iter = scan_n_iters[j];
					c->relax_steps = scan_relax_steps[k];
					c->init_mode = mode->init_mode;
					c->prior_mode = mode->prior_mode;
					c->rho_train_mode = mode->rho_train_mode;
					c->d_scale_mode = mode->d_scale_mode;
					c->init_eps = mode->init_eps;
					c->init_noise_scale = mode->init_noise_scale;
					c->init_scale = mode->init_scale;
					c->init_seed = mode->init_seed;
					snprintf(c->name, sizeof(c->name), "%s_rep%d_n%d_rs%d",
							 mode->name, (int)c->multiplier, c->n_iter, c->relax_steps);
				}
			}
		}
	}
	n_scan_configs = idx;
	return idx == HK_BLIND_P9016_SCAN_N_CONFIG? 0 : 1;
}

static int test_scan_grid(void)
{
	size_t idx = 0;
	size_t h, i, j, k;
	int failed = 0;

	if (n_scan_configs != HK_BLIND_P9016_SCAN_N_CONFIG)
		return 1;
	failed |= check_i32("scan grid count",
						(int32_t)n_scan_configs, HK_BLIND_P9016_SCAN_N_CONFIG);
	for (h = 0; h < sizeof(scan_modes) / sizeof(scan_modes[0]); ++h) {
		for (i = 0; i < sizeof(scan_multipliers) / sizeof(scan_multipliers[0]); ++i) {
			for (j = 0; j < sizeof(scan_n_iters) / sizeof(scan_n_iters[0]); ++j) {
				for (k = 0; k < sizeof(scan_relax_steps) / sizeof(scan_relax_steps[0]); ++k, ++idx) {
					char expected_name[64];
					const struct scan_mode_def *mode = &scan_modes[h];
					snprintf(expected_name, sizeof(expected_name), "%s_rep%d_n%d_rs%d",
							 mode->name, (int)scan_multipliers[i], scan_n_iters[j], scan_relax_steps[k]);
					failed |= check_true("scan mode order", strcmp(scan_configs[idx].mode_name, mode->name) == 0);
					failed |= check_true("scan multiplier order",
										 fabsf(scan_configs[idx].multiplier - scan_multipliers[i]) < 1e-6f);
					failed |= check_i32("scan n_iter order", scan_configs[idx].n_iter, scan_n_iters[j]);
					failed |= check_i32("scan relax_steps order", scan_configs[idx].relax_steps, scan_relax_steps[k]);
					failed |= check_i32("scan init mode order", scan_configs[idx].init_mode, mode->init_mode);
					failed |= check_i32("scan prior mode order", scan_configs[idx].prior_mode, mode->prior_mode);
					failed |= check_i32("scan rho train mode order", scan_configs[idx].rho_train_mode, mode->rho_train_mode);
					failed |= check_i32("scan d scale mode order", scan_configs[idx].d_scale_mode, mode->d_scale_mode);
					failed |= check_true("scan config name order", strcmp(scan_configs[idx].name, expected_name) == 0);
				}
			}
		}
	}
	return failed;
}

static int test_job_count_clamping(void)
{
	const char *orig = getenv("HK_BLIND_SWEEP_JOBS");
	char *saved = 0;
	int failed = 0;

	if (orig != 0) {
		saved = (char*)malloc(strlen(orig) + 1);
		if (saved == 0)
			return 1;
		strcpy(saved, orig);
	}
	setenv("HK_BLIND_SWEEP_JOBS", "999", 1);
	failed |= check_i32("scan jobs clamp", parse_sweep_jobs(7), 7);
	setenv("HK_BLIND_SWEEP_JOBS", "3", 1);
	failed |= check_i32("scan jobs parse", parse_sweep_jobs(7), 3);
	if (saved != 0) {
		setenv("HK_BLIND_SWEEP_JOBS", saved, 1);
		free(saved);
	} else {
		unsetenv("HK_BLIND_SWEEP_JOBS");
	}
	return failed;
}

static int set_unphased_scaffold(struct hk_bmap *bmap, fvec3_t *haploid)
{
	struct hk_fdg_conf fdg_conf;

	hk_fdg_conf_init(&fdg_conf);
	fdg_conf.backend = HK_FDG_BACKEND_CPU;
	fdg_conf.n_iter = HK_BLIND_P9016_SCAN_SCAFFOLD_FDG_N_ITER;
	return hk_blind_init_haploid_scaffold_from_bmap_fdg(bmap, &fdg_conf, haploid,
														HK_BLIND_P9016_SCAN_INIT_SEED);
}

static float dist3(const fvec3_t a, const fvec3_t b)
{
	float dx = a[0] - b[0];
	float dy = a[1] - b[1];
	float dz = a[2] - b[2];
	return sqrtf(dx * dx + dy * dy + dz * dz);
}

static int check_coords_finite(const fvec3_t *coords, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("scan coord finite", isfinite(coords[i][a]));
	return failed;
}

static void set_schedule_conf(struct hk_blind_iter_schedule_conf *conf,
							  const struct scan_config *config)
{
	conf->n_iter = config->n_iter;
	conf->base_conf.unit = HK_BLIND_P9016_SCAN_UNIT;
	conf->base_conf.d_scale = HK_BLIND_P9016_SCAN_D_SCALE;
	conf->base_conf.base_k = HK_BLIND_P9016_SCAN_BASE_K;
	conf->base_conf.temperature = HK_BLIND_P9016_SCAN_TEMPERATURE_START;
	conf->base_conf.rho_train = HK_BLIND_P9016_SCAN_RHO_TRAIN_START;
	conf->base_conf.min_sep_unit = HK_BLIND_P9016_SCAN_MIN_SEP_UNIT;
	conf->base_conf.lambda_sep = HK_BLIND_P9016_SCAN_LAMBDA_SEP;
	conf->base_conf.relax_step = HK_BLIND_P9016_SCAN_RELAX_STEP;
	conf->base_conf.relax_steps = config->relax_steps;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.rho_train_mode = config->rho_train_mode;
	conf->base_conf.d_scale_mode = config->d_scale_mode;
	conf->base_conf.d_scale_eps_count = 1e-6f;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = HK_BLIND_P9016_SCAN_TEMPERATURE_START;
	conf->temperature_end = HK_BLIND_P9016_SCAN_TEMPERATURE_END;
	conf->rho_train_start = HK_BLIND_P9016_SCAN_RHO_TRAIN_START;
	conf->rho_train_end = HK_BLIND_P9016_SCAN_RHO_TRAIN_END;
}

static int check_binned_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
							int32_t n_raw)
{
	int failed = 0;
	int32_t i;

	failed |= check_true("scan set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_i32("scan raw count", set->n_raw, n_raw);
	failed |= check_true("scan bpair positive", set->n_bpairs > 0);
	for (i = 0; i < set->n_raw; ++i) {
		failed |= check_true("scan raw2binned id lower", set->raw2binned[i].bpair_id >= 0);
		failed |= check_true("scan raw2binned id upper", set->raw2binned[i].bpair_id < set->n_bpairs);
		failed |= check_true("scan raw2binned swapped",
							 set->raw2binned[i].swapped == 0 || set->raw2binned[i].swapped == 1);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		failed |= check_true("scan key sorted", bp->key.bid[0] <= bp->key.bid[1]);
		failed |= check_true("scan key bid0 valid", bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		failed |= check_true("scan key bid1 valid", bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
	}
	return failed;
}

static int validate_final_bpair_set(const struct hk_blind_bpair_set *set)
{
	struct hk_blind_iter_diag diag;
	int failed = 0;

	hk_blind_iter_diag_init(&diag);
	hk_blind_iter_diag_validate_bpair_set(set, &diag);
	failed |= check_i32("scan final posterior nonfinite", diag.n_posterior_nonfinite, 0);
	failed |= check_i32("scan final posterior bad sum", diag.n_posterior_bad_sum, 0);
	failed |= check_i32("scan final posterior out of range", diag.n_posterior_out_of_range, 0);
	failed |= check_i32("scan final uncertainty nonfinite", diag.n_uncertainty_nonfinite, 0);
	failed |= check_i32("scan final uncertainty out of range", diag.n_uncertainty_out_of_range, 0);
	failed |= check_i32("scan final five-state bad sum", diag.n_five_state_bad_sum, 0);
	return failed;
}

static int validate_loop_health(const struct scan_config *config,
								const struct hk_blind_iter_loop_diag *diag,
								const struct hk_blind_single_iter_diag *per_iter)
{
	const struct hk_blind_iter_diag *last_pre;
	int failed = 0;

	failed |= check_i32("scan completed", diag->n_completed, config->n_iter);
	failed |= check_i32("scan bad iter", diag->n_bad_iter, 0);
	failed |= check_i32("scan relax nonfinite iter", diag->n_relax_nonfinite_iter, 0);
	failed |= check_i32("scan coord nonfinite", diag->n_coord_nonfinite, 0);
	failed |= check_i32("scan repulsion nonfinite step total", diag->n_repulsion_nonfinite_step, 0);
	failed |= check_true("scan entropy finite", isfinite(diag->final_mean_entropy));
	failed |= check_true("scan pU finite", isfinite(diag->final_mean_pU));
	failed |= check_true("scan pU range", diag->final_mean_pU >= -1e-6f && diag->final_mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("scan sum wedge finite", isfinite(diag->final_sum_wedge_k));
	failed |= check_true("scan sum wedge nonnegative", diag->final_sum_wedge_k >= 0.0);
	failed |= check_true("scan final repulsion force finite", isfinite(diag->final_repulsion_force_l1));
	last_pre = &per_iter[config->n_iter - 1].pre_relax_diag;
	failed |= check_i32("scan final posterior nonfinite diag", last_pre->n_posterior_nonfinite, 0);
	failed |= check_i32("scan final posterior bad sum diag", last_pre->n_posterior_bad_sum, 0);
	failed |= check_i32("scan final posterior out range diag", last_pre->n_posterior_out_of_range, 0);
	failed |= check_i32("scan final uncertainty nonfinite diag", last_pre->n_uncertainty_nonfinite, 0);
	failed |= check_i32("scan final uncertainty out range diag", last_pre->n_uncertainty_out_of_range, 0);
	failed |= check_i32("scan final five-state bad sum diag", last_pre->n_five_state_bad_sum, 0);
	failed |= check_i32("scan final repulsion nonfinite step",
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
		failed |= check_true("scan metric pU finite", isfinite(pU[i]));
		failed |= check_true("scan metric pmax finite", isfinite(pmax[i]));
		failed |= check_true("scan metric margin finite", isfinite(margin[i]));
		failed |= check_true("scan metric psame finite", isfinite(same));
		failed |= check_true("scan metric pcross finite", isfinite(cross));
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
	failed |= check_true("scan mean pU range", out->mean_pU >= -1e-6f && out->mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("scan mean pmax range", out->mean_pmax >= 0.25f - 1e-6f && out->mean_pmax <= 1.0f + 1e-6f);
	failed |= check_true("scan mean margin range", out->mean_margin >= -1e-6f && out->mean_margin <= 1.0f + 1e-6f);
	failed |= check_true("scan mean psame range", out->mean_psame >= -1e-6f && out->mean_psame <= 1.0f + 1e-6f);
	failed |= check_true("scan mean pcross range", out->mean_pcross >= -1e-6f && out->mean_pcross <= 1.0f + 1e-6f);

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
		failed |= check_true("scan sep finite", isfinite(sep[i]));
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
	failed |= check_true("scan sep mean finite", isfinite(out->mean));
	failed |= check_true("scan sep median finite", isfinite(out->median));
	failed |= check_true("scan sep min finite", isfinite(out->min));
	failed |= check_true("scan sep max finite", isfinite(out->max));

cleanup:
	free(sep);
	return failed;
}

static int validate_result_metrics(const struct scan_result *result)
{
	const struct bpair_dist_metrics *b = &result->bpair_metrics;
	const struct sep_dist_metrics *s = &result->sep_metrics;
	int failed = 0;

	failed |= check_true("scan final contact energy finite", isfinite(result->final_contact_energy));
	failed |= check_true("scan final backbone energy finite", isfinite(result->final_backbone_energy));
	failed |= check_true("scan final repulsion energy finite", isfinite(result->final_repulsion_energy));
	failed |= check_true("scan final sep energy finite", isfinite(result->final_sep_energy));
	failed |= check_true("scan final total energy finite", isfinite(result->final_total_energy));
	failed |= check_true("scan final force finite", isfinite(result->final_force_l1));
	failed |= check_true("scan final repulsion force finite result", isfinite(result->final_repulsion_force_l1));
	failed |= check_true("scan final backbone force finite result", isfinite(result->final_backbone_force_l1));
	failed |= check_i32("scan result repulsion nonfinite step", result->n_repulsion_nonfinite_step, 0);
	failed |= check_true("scan mean pU finite", isfinite(b->mean_pU));
	failed |= check_true("scan median pU finite", isfinite(b->median_pU));
	failed |= check_true("scan pU p10 finite", isfinite(b->pU_p10));
	failed |= check_true("scan pU p90 finite", isfinite(b->pU_p90));
	failed |= check_true("scan mean pmax finite", isfinite(b->mean_pmax));
	failed |= check_true("scan pmax p10 finite", isfinite(b->pmax_p10));
	failed |= check_true("scan pmax p90 finite", isfinite(b->pmax_p90));
	failed |= check_true("scan mean margin finite", isfinite(b->mean_margin));
	failed |= check_true("scan margin p10 finite", isfinite(b->margin_p10));
	failed |= check_true("scan margin p90 finite", isfinite(b->margin_p90));
	failed |= check_true("scan sep mean finite result", isfinite(s->mean));
	failed |= check_true("scan sep median finite result", isfinite(s->median));
	failed |= check_true("scan sep p10 finite result", isfinite(s->p10));
	failed |= check_true("scan sep p90 finite result", isfinite(s->p90));
	failed |= check_true("scan sep min finite result", isfinite(s->min));
	failed |= check_true("scan sep max finite result", isfinite(s->max));
	failed |= check_true("scan repulsion/contact finite", isfinite(result->repulsion_over_contact));
	failed |= check_true("scan contact/wedge finite", isfinite(result->contact_energy_per_wedge_k));
	return failed;
}

static int validate_coords_gz_output(const char *path, const struct hk_bmap *bmap)
{
	struct stat st;
	gzFile fp = 0;
	char line[4096];
	uint8_t *copy_mask = 0;
	int32_t row_count = 0;
	int failed = 0;

	assert(path);
	assert(bmap);
	failed |= check_i32("scan coords gz stat", stat(path, &st), 0);
	if (failed)
		return 1;
	failed |= check_true("scan coords gz nonempty", st.st_size > 0);
	copy_mask = (uint8_t*)calloc((size_t)bmap->n_beads, sizeof(*copy_mask));
	if (copy_mask == 0)
		return 1;
	fp = gzopen(path, "rb");
	if (fp == 0) {
		free(copy_mask);
		return 1;
	}
	if (gzgets(fp, line, sizeof(line)) == 0) {
		failed = 1;
		goto cleanup;
	}
	failed |= check_true("scan coords gz header",
						 strcmp(line, "chr\tstart\tend\tbid\tcopy\tdiploid_bid\tx\ty\tz\n") == 0);
	while (gzgets(fp, line, sizeof(line)) != 0) {
		char chr[128];
		int st0 = 0, en = 0, bid = -1, copy = -1, diploid_bid = -1;
		double x = 0.0, y = 0.0, z = 0.0;
		int n_parsed;

		failed |= check_true("scan coords no phase0", strstr(line, "phase0") == 0);
		failed |= check_true("scan coords no phase1", strstr(line, "phase1") == 0);
		failed |= check_true("scan coords no truth", strstr(line, "truth") == 0);
		failed |= check_true("scan coords no oracle", strstr(line, "oracle") == 0);
		n_parsed = sscanf(line, "%127s\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf",
						  chr, &st0, &en, &bid, &copy, &diploid_bid, &x, &y, &z);
		failed |= check_i32("scan coords parsed", n_parsed, 9);
		failed |= check_true("scan coords bid range", bid >= 0 && bid < bmap->n_beads);
		failed |= check_true("scan coords copy range", copy == 0 || copy == 1);
		if (bid >= 0 && bid < bmap->n_beads && (copy == 0 || copy == 1)) {
			const struct hk_bead *bead = &bmap->beads[bid];
			copy_mask[bid] |= (uint8_t)(1u << copy);
			failed |= check_i32("scan coords start", st0, bead->st);
			failed |= check_i32("scan coords end", en, bead->en);
			if (bmap->d != 0 && bead->chr >= 0 && bead->chr < bmap->d->n && bmap->d->name != 0)
				failed |= check_true("scan coords chr", strcmp(chr, bmap->d->name[bead->chr]) == 0);
			failed |= check_i32("scan coords diploid bid", diploid_bid, hk_diploid_bid(bid, copy));
		}
		failed |= check_true("scan coords finite", isfinite(x) && isfinite(y) && isfinite(z));
		++row_count;
	}
	failed |= check_i32("scan coords row count", row_count, 2 * bmap->n_beads);
	if (!failed) {
		int32_t i;
		for (i = 0; i < bmap->n_beads; ++i)
			failed |= check_i32("scan coords copy mask", copy_mask[i], 3);
	}

cleanup:
	if (gzclose(fp) != Z_OK)
		failed = 1;
	free(copy_mask);
	return failed;
}

static const char *scan_scaffold_source(const struct scan_config *config)
{
	if (config->init_mode == HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT)
		return "unphased_fdg_1mb";
	if (config->init_mode == HK_BLIND_INIT_RANDOM_DIPLOID)
		return "random_diploid";
	if (config->init_mode == HK_BLIND_INIT_TOY_SPLIT)
		return "toy_split";
	return "unknown";
}

static int write_summary_header(FILE *fp)
{
	return fprintf(fp,
				   "config_name\toutput_dir\tcoords_gz\tmultiplier\tk_rel_rep\tn_iter\t"
				   "relax_steps\trelax_step\ttemperature_start\ttemperature_end\t"
				   "rho_train_start\trho_train_end\tunit\td_scale\tlegacy_base_k_unused\tmin_sep_unit\t"
				   "lambda_sep\tinit_eps\tinit_noise_scale\tinit_seed\tn_raw\tn_bpair\t"
				   "n_beads\tn_completed\tn_bad_iter\tn_relax_nonfinite_iter\t"
				   "n_coord_nonfinite\ttotal_chr_flipped\tstatus\telapsed_cpu_sec\t"
				   "final_mean_entropy\tfinal_mean_pU\tmean_pU\tmedian_pU\tpU_p10\t"
				   "pU_p90\tfrac_pU_lt_0_90\tfrac_pU_lt_0_75\tfrac_pU_lt_0_50\t"
				   "mean_pmax\tmedian_pmax\tpmax_p10\tpmax_p90\tfrac_pmax_gt_0_50\t"
				   "frac_pmax_gt_0_75\tmean_margin\tmedian_margin\tmargin_p10\t"
				   "margin_p90\tfrac_margin_gt_0_10\tfrac_margin_gt_0_25\t"
				   "mean_psame\tmean_pcross\tfrac_psame_gt_0_75\tfrac_pcross_gt_0_75\t"
				   "sep_mean\tsep_median\tsep_p10\tsep_p90\tsep_min\tsep_max\t"
				   "frac_sep_lt_0_25\tfrac_sep_lt_0_50\tfrac_sep_lt_1_00\t"
				   "final_contact_energy\tfinal_backbone_energy\tfinal_repulsion_energy\t"
				   "final_sep_energy\tfinal_total_energy\tfinal_force_l1\t"
				   "final_repulsion_force_l1\tfinal_backbone_force_l1\tn_repulsion_nonfinite_step\t"
				   "repulsion_over_contact\tcontact_energy_per_wedge_k\tfinal_sum_wedge_k\t"
				   "mode_group\tinit_mode\tprior_mode\trho_train_mode\td_scale_mode\t"
				   "scaffold_source\tinit_scale\tinit_seed_effective\t"
				   "same_bin_filter_enabled\tn_raw_same_bin_excluded\tn_bpair_same_bin_excluded\t"
				   "posterior_refreshed_after_final_relax\tfinal_mean_rho_train_bpair\t"
				   "final_min_rho_train_bpair\tfinal_max_rho_train_bpair\n") < 0? -1 : 0;
}

static int append_summary_row(FILE *fp, const struct scan_result *result)
{
	const struct scan_config *c = result->config;
	const struct hk_blind_iter_loop_diag *d = &result->loop_diag;
	const struct bpair_dist_metrics *b = &result->bpair_metrics;
	const struct sep_dist_metrics *s = &result->sep_metrics;
	const char *status = result->status_ok? "OK" : "FAIL";

	return fprintf(fp,
				   "%s\t%s\t%s\t%.9g\t%.9g\t%d\t%d\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%llu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%.3f\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				   "%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.17g\t%s\t%s\t%s\t%s\t%s\t"
				   "%s\t%.9g\t%llu\t%d\t%lld\t%d\t%d\t%.9g\t%.9g\t%.9g\n",
				   c->name, result->output_dir, result->coords_gz,
				   c->multiplier, result->k_rel_rep, c->n_iter, c->relax_steps,
				   HK_BLIND_P9016_SCAN_RELAX_STEP,
				   HK_BLIND_P9016_SCAN_TEMPERATURE_START,
				   HK_BLIND_P9016_SCAN_TEMPERATURE_END,
				   HK_BLIND_P9016_SCAN_RHO_TRAIN_START,
				   HK_BLIND_P9016_SCAN_RHO_TRAIN_END,
				   HK_BLIND_P9016_SCAN_UNIT,
				   HK_BLIND_P9016_SCAN_D_SCALE,
				   HK_BLIND_P9016_SCAN_BASE_K,
				   HK_BLIND_P9016_SCAN_MIN_SEP_UNIT,
				   HK_BLIND_P9016_SCAN_LAMBDA_SEP,
				   c->init_eps,
				   c->init_noise_scale,
				   (unsigned long long)c->init_seed,
				   result->n_raw, result->n_bpair, result->n_beads,
				   d->n_completed, d->n_bad_iter, d->n_relax_nonfinite_iter,
				   d->n_coord_nonfinite, d->total_chr_flipped, status,
				   result->elapsed_cpu_sec, d->final_mean_entropy, d->final_mean_pU,
				   b->mean_pU, b->median_pU, b->pU_p10, b->pU_p90,
				   b->frac_pU_lt_0_90, b->frac_pU_lt_0_75, b->frac_pU_lt_0_50,
				   b->mean_pmax, b->median_pmax, b->pmax_p10, b->pmax_p90,
				   b->frac_pmax_gt_0_50, b->frac_pmax_gt_0_75,
				   b->mean_margin, b->median_margin, b->margin_p10, b->margin_p90,
				   b->frac_margin_gt_0_10, b->frac_margin_gt_0_25,
				   b->mean_psame, b->mean_pcross,
				   b->frac_psame_gt_0_75, b->frac_pcross_gt_0_75,
				   s->mean, s->median, s->p10, s->p90, s->min, s->max,
				   s->frac_lt_0_25, s->frac_lt_0_50, s->frac_lt_1_00,
				   result->final_contact_energy, result->final_backbone_energy,
				   result->final_repulsion_energy, result->final_sep_energy,
				   result->final_total_energy, result->final_force_l1,
				   result->final_repulsion_force_l1, result->final_backbone_force_l1,
				   result->n_repulsion_nonfinite_step,
				   result->repulsion_over_contact, result->contact_energy_per_wedge_k,
				   d->final_sum_wedge_k,
				   c->mode_name,
				   hk_blind_init_mode_name(c->init_mode),
				   hk_blind_prior_mode_name(c->prior_mode),
				   hk_blind_rho_train_mode_name(c->rho_train_mode),
				   hk_blind_d_scale_mode_name(c->d_scale_mode),
				   scan_scaffold_source(c),
				   c->init_scale,
				   (unsigned long long)c->init_seed,
				   result->same_bin_filter_enabled,
				   (long long)result->n_raw_same_bin_excluded,
				   result->n_bpair_same_bin_excluded,
				   d->posterior_refreshed_after_final_relax,
				   d->final_mean_rho_train_bpair,
				   d->final_min_rho_train_bpair,
				   d->final_max_rho_train_bpair) < 0? -1 : 0;
}

static int run_one_config(const struct scan_config *config, const struct hk_bmap *bmap,
						  const struct hk_blind_pair *raw, int32_t n_raw,
						  const fvec3_t *haploid, const char *root_dir,
						  FILE *summary_fp, struct scan_result *result)
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
	path_join(result->output_dir, sizeof(result->output_dir), root_dir, config->name);
	path_join(result->coords_gz, sizeof(result->coords_gz), result->output_dir, HK_BLIND_P9016_SCAN_COORDS_NAME);
	if (mkdir(result->output_dir, 0700) != 0 && errno != EEXIST) {
		fprintf(stderr, "failed to create scan config dir %s\n", result->output_dir);
		failed = 1;
		goto finish;
	}

	fprintf(stderr, "P9016 main/random repel grid %s: mode=%s multiplier=%.4g n_iter=%d relax_steps=%d\n",
			config->name, config->mode_name, config->multiplier, config->n_iter, config->relax_steps);
	set = hk_blind_bpair_set_build(bmap, n_raw, raw);
	failed |= check_binned_set(bmap, set, n_raw);
	if (set) {
		struct hk_blind_prior_diag prior_diag;

		result->n_bpair = set->n_bpairs;
		result->same_bin_filter_enabled = set->same_bin_filter_enabled;
		result->n_raw_same_bin_excluded = set->n_raw_same_bin_excluded;
		result->n_bpair_same_bin_excluded = set->n_bpair_same_bin_excluded;
		if (config->prior_mode == HK_BLIND_PRIOR_CIS_INTER_RATIO) {
			failed |= check_i32("scan cis_inter_ratio prior",
								hk_blind_bpair_set_init_cis_inter_ratio_prior(bmap, set,
																			   HK_BLIND_P9016_SCAN_PRIOR_EPS,
																			   &prior_diag), 0);
		} else if (config->prior_mode == HK_BLIND_PRIOR_UNIFORM) {
			hk_blind_bpair_set_init_uniform_prior(set);
		} else {
			fprintf(stderr, "scan config %s has unsupported prior mode %d\n",
					config->name, config->prior_mode);
			failed = 1;
		}
	}
	if (failed) goto finish;

	diploid = (fvec3_t*)calloc((size_t)n_diploid, sizeof(*diploid));
	per_iter = (struct hk_blind_single_iter_diag*)calloc((size_t)config->n_iter, sizeof(*per_iter));
	if (diploid == 0 || per_iter == 0) {
		fprintf(stderr, "failed to allocate scan coordinate/diagnostic data for %s\n", config->name);
		failed = 1;
		goto finish;
	}
	if (config->init_mode == HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT) {
		failed |= check_i32("scan init unphased scaffold split",
							hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
																	  diploid, config->init_eps,
																	  config->init_noise_scale,
																	  config->init_seed), 0);
	} else if (config->init_mode == HK_BLIND_INIT_RANDOM_DIPLOID) {
		failed |= check_i32("scan init random diploid",
							hk_blind_init_random_diploid_coords(bmap, diploid,
																config->init_scale,
																config->init_seed), 0);
	} else {
		fprintf(stderr, "scan config %s has unsupported init mode %d\n",
				config->name, config->init_mode);
		failed = 1;
	}
	failed |= check_coords_finite(diploid, n_diploid);
	if (failed) goto finish;

	hk_fdg_conf_init(&fdg_conf);
	fdg_conf.k_rel_rep *= config->multiplier;
	result->k_rel_rep = fdg_conf.k_rel_rep;
	set_schedule_conf(&schedule_conf, config);
	ret = hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &fdg_conf, diploid, 0,
											   &schedule_conf, per_iter, &result->loop_diag);
	failed |= check_i32("scan scheduled ret", ret, 0);
	if (ret == 0) {
		const struct hk_blind_relax_diag *rd = &per_iter[config->n_iter - 1].relax_diag;
		result->final_contact_energy = rd->final_contact_energy;
		result->final_backbone_energy = rd->final_backbone_energy;
		result->final_repulsion_energy = rd->final_repulsion_energy;
		result->final_sep_energy = rd->final_sep_energy;
		result->final_total_energy = rd->final_total_energy;
		result->final_force_l1 = rd->final_force_l1;
		result->final_repulsion_force_l1 = rd->final_repulsion_force_l1;
		result->final_backbone_force_l1 = rd->final_backbone_force_l1;
		result->n_repulsion_nonfinite_step = result->loop_diag.n_repulsion_nonfinite_step;
		failed |= validate_loop_health(config, &result->loop_diag, per_iter);
		failed |= check_coords_finite(diploid, n_diploid);
		failed |= validate_final_bpair_set(set);
		failed |= compute_bpair_metrics(set, &result->bpair_metrics);
		failed |= compute_sep_metrics(bmap, diploid, &result->sep_metrics);
		result->repulsion_over_contact = result->final_contact_energy > 0.0f?
			(double)result->final_repulsion_energy / (double)result->final_contact_energy : 0.0;
		result->contact_energy_per_wedge_k = result->loop_diag.final_sum_wedge_k > 0.0?
			(double)result->final_contact_energy / result->loop_diag.final_sum_wedge_k : 0.0;
		failed |= validate_result_metrics(result);
	}
	if (diploid != 0) {
		if (hk_blind_write_diploid_coords_tsv_gz(result->coords_gz, bmap, diploid) != 0) {
			fprintf(stderr, "failed to write scan coords %s\n", result->coords_gz);
			failed = 1;
		} else if (validate_coords_gz_output(result->coords_gz, bmap) != 0) {
			fprintf(stderr, "failed to validate scan coords %s\n", result->coords_gz);
			failed = 1;
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
			"%s mean_pU=%.6g mean_pmax=%.6g mean_margin=%.6g sep_mean=%.6g "
			"contact_per_wedge=%.6g repulsion_over_contact=%.6g coords=%s status=%s\n",
			config->name, result->bpair_metrics.mean_pU, result->bpair_metrics.mean_pmax,
			result->bpair_metrics.mean_margin, result->sep_metrics.mean,
			result->contact_energy_per_wedge_k, result->repulsion_over_contact,
			result->coords_gz, status);

	free(diploid);
	free(per_iter);
	hk_blind_bpair_set_destroy(set);
	return failed != 0;
}

static void config_output_dir(char *dst, size_t dst_size, const char *root_dir,
							  const struct scan_config *config)
{
	path_join(dst, dst_size, root_dir, config->name);
}

static void config_summary_row_path(char *dst, size_t dst_size, const char *config_dir,
									const struct scan_config *config)
{
	char row_name[128];
	snprintf(row_name, sizeof(row_name), "%s.summary_row.tsv", config->name);
	path_join(dst, dst_size, config_dir, row_name);
}

static int summary_row_resume_ok(const char *row_path, const struct scan_config *config,
								 const struct hk_bmap *bmap)
{
	FILE *fp = fopen(row_path, "r");
	char line[8192];
	char parse_buf[8192];
	struct parsed_summary_row row;
	int ok = 0;

	if (fp == 0)
		return 0;
	if (fgets(line, sizeof(line), fp) == 0)
		goto finish;
	if (fgets(parse_buf, sizeof(parse_buf), fp) != 0)
		goto finish;
	strncpy(parse_buf, line, sizeof(parse_buf) - 1);
	parse_buf[sizeof(parse_buf) - 1] = 0;
	if (parse_summary_row(parse_buf, &row) != 0)
		goto finish;
	if (strcmp(row.config_name, config->name) != 0 || !row.status_ok)
		goto finish;
	if (validate_coords_gz_output(row.coords_gz, bmap) != 0)
		goto finish;
	ok = 1;

finish:
	fclose(fp);
	return ok;
}

static int run_one_config_to_row_file(size_t config_idx, const struct hk_bmap *bmap,
									  const struct hk_blind_pair *raw, int32_t n_raw,
									  const fvec3_t *haploid, const char *root_dir,
									  const char *row_path)
{
	const struct scan_config *config = &scan_configs[config_idx];
	char config_dir[1024];
	FILE *fp;
	struct scan_result result;
	int ret;

	if (!parse_env_flag("HK_BLIND_SCAN_FORCE") &&
		summary_row_resume_ok(row_path, config, bmap)) {
		fprintf(stderr, "P9016 large scan %s: resume skip from %s\n",
				config->name, row_path);
		return 0;
	}
	config_output_dir(config_dir, sizeof(config_dir), root_dir, config);
	if (mkdir(config_dir, 0700) != 0 && errno != EEXIST) {
		fprintf(stderr, "failed to create scan config dir %s\n", config_dir);
		return 1;
	}
	fp = fopen(row_path, "w");
	if (fp == 0) {
		fprintf(stderr, "failed to open scan summary row %s\n", row_path);
		return 1;
	}
	ret = run_one_config(config, bmap, raw, n_raw, haploid, root_dir, fp, &result);
	if (fclose(fp) != 0)
		ret = 1;
	return ret != 0;
}

static int append_summary_row_file(FILE *summary_fp, const char *row_path,
								   const char *expected_config_name,
								   const struct hk_bmap *bmap)
{
	FILE *fp = fopen(row_path, "r");
	char line[8192];
	char parse_buf[8192];
	struct parsed_summary_row row;
	int failed = 0;

	if (fp == 0) {
		fprintf(stderr, "failed to open completed scan summary row %s\n", row_path);
		return 1;
	}
	if (fgets(line, sizeof(line), fp) == 0) {
		fprintf(stderr, "scan summary row %s is empty\n", row_path);
		failed = 1;
		goto finish;
	}
	strncpy(parse_buf, line, sizeof(parse_buf) - 1);
	parse_buf[sizeof(parse_buf) - 1] = 0;
	if (parse_summary_row(parse_buf, &row) != 0 ||
		strcmp(row.config_name, expected_config_name) != 0) {
		fprintf(stderr, "scan summary row %s failed parse/config validation\n", row_path);
		failed = 1;
		goto finish;
	}
	if (!row.status_ok) {
		fprintf(stderr, "scan summary row %s has FAIL status\n", row_path);
		failed = 1;
		goto finish;
	}
	if (validate_coords_gz_output(row.coords_gz, bmap) != 0) {
		fprintf(stderr, "scan summary row %s has invalid coords_gz %s\n",
				row_path, row.coords_gz);
		failed = 1;
		goto finish;
	}
	if (fgets(parse_buf, sizeof(parse_buf), fp) != 0) {
		fprintf(stderr, "scan summary row %s contains more than one row\n", row_path);
		failed = 1;
		goto finish;
	}
	if (fputs(line, summary_fp) == EOF)
		failed = 1;

finish:
	if (ferror(fp)) failed = 1;
	fclose(fp);
	return failed;
}

static int run_configs_parallel(int jobs, const struct hk_bmap *bmap,
								const struct hk_blind_pair *raw, int32_t n_raw,
								const fvec3_t *haploid, const char *root_dir,
								FILE *summary_fp)
{
	size_t n_configs = n_scan_configs;
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
		fprintf(stderr, "failed to allocate scan parallel bookkeeping\n");
		free(row_paths);
		free(pids);
		free(done);
		return 1;
	}
	for (i = 0; i < n_configs; ++i) {
		char config_dir[1024];
		config_output_dir(config_dir, sizeof(config_dir), root_dir, &scan_configs[i]);
		config_summary_row_path(row_paths[i], sizeof(row_paths[i]), config_dir, &scan_configs[i]);
	}

	while (next < n_configs || running > 0) {
		while (next < n_configs && running < (size_t)jobs) {
			size_t idx = next++;
			pid_t pid = fork();
			if (pid < 0) {
				fprintf(stderr, "failed to fork scan config %s\n", scan_configs[idx].name);
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
				fprintf(stderr, "wait failed during scan parallel run\n");
				failed = 1;
				break;
			}
			--running;
			for (i = 0; i < n_configs; ++i) {
				if (pids[i] == pid) {
					done[i] = 1;
					if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
						fprintf(stderr, "scan config %s failed in child\n", scan_configs[i].name);
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
			fprintf(stderr, "scan config %s did not complete\n", scan_configs[i].name);
			failed = 1;
			continue;
		}
		if (append_summary_row_file(summary_fp, row_paths[i], scan_configs[i].name, bmap) != 0) {
			fprintf(stderr, "failed to append summary row for scan config %s\n", scan_configs[i].name);
			failed = 1;
		}
	}
	free(row_paths);
	free(pids);
	free(done);
	return failed != 0;
}

static int parse_summary_float(const char *token, float *out)
{
	char *end = 0;

	errno = 0;
	*out = strtof(token, &end);
	return errno == 0 && end != token && *end == 0 && isfinite(*out);
}

static int parse_summary_double(const char *token, double *out)
{
	char *end = 0;

	errno = 0;
	*out = strtod(token, &end);
	return errno == 0 && end != token && *end == 0 && isfinite(*out);
}

static int parse_summary_row(char *line, struct parsed_summary_row *out)
{
	char *save = 0;
	char *token;
	int col = 0;
	int ok = 1;

	memset(out, 0, sizeof(*out));
	for (token = strtok_r(line, "\t\r\n", &save); token != 0;
		 token = strtok_r(0, "\t\r\n", &save), ++col) {
		switch (col) {
		case 0:
			strncpy(out->config_name, token, sizeof(out->config_name) - 1);
			break;
		case 2:
			strncpy(out->coords_gz, token, sizeof(out->coords_gz) - 1);
			break;
		case 28:
			out->status_ok = strcmp(token, "OK") == 0;
			break;
		case 32:
			ok &= parse_summary_float(token, &out->mean_pU);
			break;
		case 39:
			ok &= parse_summary_float(token, &out->mean_pmax);
			break;
		case 45:
			ok &= parse_summary_float(token, &out->mean_margin);
			break;
		case 55:
			ok &= parse_summary_float(token, &out->sep_mean);
			break;
		case 61:
			ok &= parse_summary_float(token, &out->frac_sep_lt_0_25);
			break;
		case 73:
			ok &= parse_summary_double(token, &out->repulsion_over_contact);
			break;
		case 74:
			ok &= parse_summary_double(token, &out->contact_energy_per_wedge_k);
			break;
		default:
			break;
		}
	}
	return ok && col >= 76? 0 : 1;
}

static int validate_summary_file(const char *summary_path, int n_expected,
								 int *n_ok_out, int *n_fail_out)
{
	FILE *fp = fopen(summary_path, "r");
	char line[8192];
	int n_rows = -1;
	int n_ok = 0;
	int n_fail = 0;
	int failed = 0;

	if (fp == 0) return 1;
	while (fgets(line, sizeof(line), fp)) {
		if (n_rows < 0) {
			failed |= check_true("scan summary header config", strstr(line, "config_name") != 0);
			failed |= check_true("scan summary header coords", strstr(line, "coords_gz") != 0);
			failed |= check_true("scan summary header pU p10", strstr(line, "pU_p10") != 0);
			failed |= check_true("scan summary header repulsion", strstr(line, "final_repulsion_energy") != 0);
			failed |= check_true("scan summary header repulsion force", strstr(line, "final_repulsion_force_l1") != 0);
			failed |= check_true("scan summary header wedge", strstr(line, "final_sum_wedge_k") != 0);
			failed |= check_true("scan summary header mode", strstr(line, "mode_group") != 0);
			failed |= check_true("scan summary header init", strstr(line, "init_mode") != 0);
			failed |= check_true("scan summary header rho", strstr(line, "rho_train_mode") != 0);
		} else {
			struct parsed_summary_row row;
			if (parse_summary_row(line, &row) != 0) {
				fprintf(stderr, "failed to parse scan summary row %d\n", n_rows);
				failed = 1;
			} else if (row.status_ok) {
				++n_ok;
			} else {
				++n_fail;
			}
		}
		++n_rows;
	}
	if (ferror(fp)) failed = 1;
	fclose(fp);
	failed |= check_i32("scan summary row count", n_rows, n_expected);
	*n_ok_out = n_ok;
	*n_fail_out = n_fail;
	return failed;
}

static void print_compact_summary_row(const struct parsed_summary_row *row)
{
	fprintf(stderr, "%s\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%s\t%s\n",
			row->config_name, row->mean_pU, row->mean_pmax, row->mean_margin,
			row->sep_mean, row->contact_energy_per_wedge_k,
			row->repulsion_over_contact, row->coords_gz, row->status_ok? "OK" : "FAIL");
}

static int print_interpretation_summary(const char *summary_path, int n_ok, int n_fail)
{
	FILE *fp = fopen(summary_path, "r");
	char line[8192];
	struct parsed_summary_row *rows = 0, *next_rows;
	int n_rows = 0, m_rows = 0;
	int row_idx = -1, i;
	int failed = 0;
	int has_best_pU = 0, has_best_margin = 0, has_best_pmax = 0;
	int has_contact_1p0 = 0, has_contact_1p2 = 0, has_sep = 0;
	struct parsed_summary_row best_pU, best_margin, best_pmax;
	struct parsed_summary_row best_contact_1p0, best_contact_1p2, best_sep;

	if (fp == 0) {
		fprintf(stderr, "failed to open scan summary for interpretation: %s\n", summary_path);
		return 1;
	}
	fprintf(stderr, "scan OK=%d FAIL=%d\n", n_ok, n_fail);
	while (fgets(line, sizeof(line), fp)) {
		struct parsed_summary_row row;
		if (row_idx < 0) {
			row_idx = 0;
			continue;
		}
		if (parse_summary_row(line, &row) != 0) {
			fprintf(stderr, "failed to parse scan interpretation row %d\n", row_idx);
			failed = 1;
			break;
		}
		if (n_rows == m_rows) {
			int new_m = m_rows? m_rows * 2 : 256;
			next_rows = (struct parsed_summary_row*)realloc(rows, (size_t)new_m * sizeof(*rows));
			if (next_rows == 0) {
				failed = 1;
				break;
			}
			rows = next_rows;
			m_rows = new_m;
		}
		rows[n_rows++] = row;
		if (row.status_ok) {
			if (!has_best_pU || row.mean_pU < best_pU.mean_pU) {
				has_best_pU = 1;
				best_pU = row;
			}
			if (!has_best_margin || row.mean_margin > best_margin.mean_margin) {
				has_best_margin = 1;
				best_margin = row;
			}
			if (!has_best_pmax || row.mean_pmax > best_pmax.mean_pmax) {
				has_best_pmax = 1;
				best_pmax = row;
			}
			if (row.contact_energy_per_wedge_k <= 1.0 &&
				(!has_contact_1p0 || row.mean_pU < best_contact_1p0.mean_pU)) {
				has_contact_1p0 = 1;
				best_contact_1p0 = row;
			}
			if (row.contact_energy_per_wedge_k <= 1.2 &&
				(!has_contact_1p2 || row.mean_pU < best_contact_1p2.mean_pU)) {
				has_contact_1p2 = 1;
				best_contact_1p2 = row;
			}
			if (row.frac_sep_lt_0_25 == 0.0f &&
				(!has_sep || row.mean_pU < best_sep.mean_pU)) {
				has_sep = 1;
				best_sep = row;
			}
		}
		++row_idx;
	}
	if (ferror(fp))
		failed = 1;
	fclose(fp);
	if (failed) {
		free(rows);
		return 1;
	}
	if (!has_best_pU || !has_best_margin) {
		fprintf(stderr, "scan interpretation summary has no OK rows\n");
		free(rows);
		return 1;
	}
	fprintf(stderr, "scan parsed rows=%d summary_path=%s\n", n_rows, summary_path);
	fprintf(stderr, "best mean_pU: %s mean_pU=%.6g coords=%s\n",
			best_pU.config_name, best_pU.mean_pU, best_pU.coords_gz);
	fprintf(stderr, "best mean_margin: %s mean_margin=%.6g coords=%s\n",
			best_margin.config_name, best_margin.mean_margin, best_margin.coords_gz);
	fprintf(stderr, "best mean_pmax: %s mean_pmax=%.6g coords=%s\n",
			best_pmax.config_name, best_pmax.mean_pmax, best_pmax.coords_gz);
	if (has_sep)
		fprintf(stderr, "best frac_sep_lt_0_25==0: %s mean_pU=%.6g sep_mean=%.6g coords=%s\n",
				best_sep.config_name, best_sep.mean_pU, best_sep.sep_mean, best_sep.coords_gz);
	else
		fprintf(stderr, "best frac_sep_lt_0_25==0: none\n");
	if (has_contact_1p0)
		fprintf(stderr, "best contact_energy_per_wedge_k<=1.0: %s mean_pU=%.6g contact_per_wedge=%.6g coords=%s\n",
				best_contact_1p0.config_name, best_contact_1p0.mean_pU,
				best_contact_1p0.contact_energy_per_wedge_k, best_contact_1p0.coords_gz);
	else
		fprintf(stderr, "best contact_energy_per_wedge_k<=1.0: none\n");
	if (has_contact_1p2)
		fprintf(stderr, "best contact_energy_per_wedge_k<=1.2: %s mean_pU=%.6g contact_per_wedge=%.6g coords=%s\n",
				best_contact_1p2.config_name, best_contact_1p2.mean_pU,
				best_contact_1p2.contact_energy_per_wedge_k, best_contact_1p2.coords_gz);
	else
		fprintf(stderr, "best contact_energy_per_wedge_k<=1.2: none\n");
	fprintf(stderr, "top config compact table:\n");
	fprintf(stderr, "config_name\tmean_pU\tmean_pmax\tmean_margin\tsep_mean\t"
			"contact_energy_per_wedge_k\trepulsion_over_contact\tcoords_gz\tstatus\n");
	print_compact_summary_row(&best_pU);
	print_compact_summary_row(&best_margin);
	print_compact_summary_row(&best_pmax);
	if (has_sep) print_compact_summary_row(&best_sep);
	if (has_contact_1p0) print_compact_summary_row(&best_contact_1p0);
	if (has_contact_1p2) print_compact_summary_row(&best_contact_1p2);
	for (i = 0; i < n_rows; ++i) {
		if (!rows[i].status_ok)
			continue;
		if (strcmp(rows[i].config_name, best_pU.config_name) == 0 ||
			strcmp(rows[i].config_name, best_margin.config_name) == 0 ||
			strcmp(rows[i].config_name, best_pmax.config_name) == 0)
			continue;
		if (rows[i].mean_pU <= best_pU.mean_pU + 0.01f)
			print_compact_summary_row(&rows[i]);
	}
	free(rows);
	return 0;
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
	int jobs;
	int failed = 0;
	int n_ok = 0, n_fail = 0;
	int32_t n_raw;
	int32_t i;
	time_t t0 = time(0), t1;
	double elapsed_sec;

	if (make_scan_configs() != 0)
		return 1;
	if (test_scan_grid() != 0)
		return 1;
	if (test_job_count_clamping() != 0)
		return 1;
	jobs = parse_sweep_jobs(n_scan_configs);
	if (jobs < 0)
		return 1;
	if (!file_exists(HK_BLIND_P9016_SCAN_PATH)) {
		fprintf(stderr, "ERROR: %s not found; P9016 main/random repel grid did not start\n",
				HK_BLIND_P9016_SCAN_PATH);
		return 1;
	}

	hk_verbose = 0;
	fprintf(stderr, "P9016 main/random repel grid: reading %s\n", HK_BLIND_P9016_SCAN_PATH);
	m = hk_map_read(HK_BLIND_P9016_SCAN_PATH);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", HK_BLIND_P9016_SCAN_PATH);
		return 1;
	}
	n_raw = m->n_pairs;
	if (n_raw <= 0) {
		fprintf(stderr, "scan input raw positive: predicate failed\n");
		failed = 1;
		goto cleanup;
	}

	fprintf(stderr, "P9016 main/random repel grid: building 1Mb bmap from %d raw pairs\n", n_raw);
	bmap = hk_bmap_gen(m->d, n_raw, m->pairs, HK_BLIND_P9016_SCAN_RESOLUTION, 1);
	raw = (struct hk_blind_pair*)calloc((size_t)n_raw, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate scan bmap/raw data\n");
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n_raw; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);
	failed |= check_true("scan bead count positive", bmap->n_beads > 0);
	if (failed) goto cleanup;

	haploid = (fvec3_t*)calloc((size_t)bmap->n_beads, sizeof(*haploid));
	if (haploid == 0) {
		fprintf(stderr, "failed to allocate scan haploid scaffold\n");
		failed = 1;
		goto cleanup;
	}
	fprintf(stderr, "P9016 main/random repel grid: building unphased 1Mb FDG scaffold (%d iterations)\n",
			HK_BLIND_P9016_SCAN_SCAFFOLD_FDG_N_ITER);
	failed |= check_i32("scan unphased scaffold", set_unphased_scaffold(bmap, haploid), 0);
	failed |= check_coords_finite(haploid, bmap->n_beads);
	if (failed) goto cleanup;

	failed |= check_i32("scan output root", make_output_root(root_dir, sizeof(root_dir)), 0);
	if (failed) goto cleanup;
	path_join(summary_path, sizeof(summary_path), root_dir, HK_BLIND_P9016_SCAN_SUMMARY_NAME);
	summary_fp = fopen(summary_path, "w");
	if (summary_fp == 0) {
		fprintf(stderr, "failed to open scan summary %s\n", summary_path);
		failed = 1;
		goto cleanup;
	}
	failed |= check_i32("scan summary header", write_summary_header(summary_fp), 0);
	if (failed) goto cleanup;
	fflush(summary_fp);

	fprintf(stderr, "P9016 main/random repel grid: output root %s jobs=%d configs=%zu force=%d\n",
			root_dir, jobs, n_scan_configs, parse_env_flag("HK_BLIND_SCAN_FORCE"));
	failed |= run_configs_parallel(jobs, bmap, raw, n_raw, haploid, root_dir, summary_fp);
	if (fclose(summary_fp) != 0) {
		summary_fp = 0;
		failed = 1;
		goto cleanup;
	}
	summary_fp = 0;
	failed |= validate_summary_file(summary_path,
									(int)n_scan_configs,
									&n_ok, &n_fail);
	if (print_interpretation_summary(summary_path, n_ok, n_fail) != 0)
		failed = 1;
	if (n_fail != 0)
		failed = 1;
	if (failed) goto cleanup;

	t1 = time(0);
	elapsed_sec = difftime(t1, t0);
	fprintf(stderr, "P9016 main/random repel grid: output_root=%s summary=%s elapsed_wall_sec=%.3f status=OK\n",
			root_dir, summary_path, elapsed_sec);

cleanup:
	if (summary_fp) fclose(summary_fp);
	if (failed)
		fprintf(stderr, "P9016 main/random repel grid failed%s%s\n",
				root_dir[0]? "; partial output remains in " : "",
				root_dir[0]? root_dir : "");
	free(haploid);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
