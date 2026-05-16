#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include "hkpriv.h"

#define HK_BLIND_P9016_CLI_DEFAULT_INPUT "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_CLI_RESOLUTION 1000000
#define HK_BLIND_P9016_CLI_UNIT 1.0f
#define HK_BLIND_P9016_CLI_D_SCALE 1.0f
#define HK_BLIND_P9016_CLI_LEGACY_BASE_K_UNUSED 2.0f
#define HK_BLIND_P9016_CLI_TEMPERATURE_START 2.0f
#define HK_BLIND_P9016_CLI_TEMPERATURE_END 1.0f
#define HK_BLIND_P9016_CLI_RHO_TRAIN_START 1.0f
#define HK_BLIND_P9016_CLI_RHO_TRAIN_END 1.0f
#define HK_BLIND_P9016_CLI_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_CLI_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_CLI_RELAX_STEP 0.001f
#define HK_BLIND_P9016_CLI_RELAX_STEPS 5
#define HK_BLIND_P9016_CLI_N_ITER 3
#define HK_BLIND_P9016_CLI_INIT_EPS 0.5f
#define HK_BLIND_P9016_CLI_INIT_NOISE_SCALE 0.0f
#define HK_BLIND_P9016_CLI_INIT_SEED 17ULL
#define HK_BLIND_P9016_CLI_INIT_SCALE 1.0f
#define HK_BLIND_P9016_CLI_PRIOR_EPS 1e-6f
#define HK_BLIND_P9016_CLI_D_SCALE_EPS_COUNT 1e-6f
#define HK_BLIND_P9016_CLI_SCAFFOLD_FDG_N_ITER 50
#define HK_BLIND_P9016_CLI_BASE_K_NEIGHBOR_RADIUS 10000000
#define HK_BLIND_P9016_CLI_RUNNER_FAMILY "blind_p9016_cli"
#define HK_BLIND_P9016_CLI_RUNNER_VERSION "2026-04-30"
#define HK_BLIND_P9016_CLI_DEFAULT_PROFILE "p9016_1mb_cli_pairs_only_v1"
#define HK_BLIND_P9016_CLI_HASH_MOD 1000000ULL

struct hk_blind_p9016_cli_conf {
	const char *input_path;
	const char *out_dir;
	int n_iter;
	int relax_steps;
	float relax_step;
	float temperature_start;
	float temperature_end;
	float rho_train_start;
	float rho_train_end;
	int rho_train_mode;
	int base_k_mode;
	float heldout_fraction;
	uint64_t heldout_seed;
	int write_raw;
};

static void cli_conf_init(struct hk_blind_p9016_cli_conf *conf)
{
	memset(conf, 0, sizeof(*conf));
	conf->input_path = HK_BLIND_P9016_CLI_DEFAULT_INPUT;
	conf->out_dir = 0;
	conf->n_iter = HK_BLIND_P9016_CLI_N_ITER;
	conf->relax_steps = HK_BLIND_P9016_CLI_RELAX_STEPS;
	conf->relax_step = HK_BLIND_P9016_CLI_RELAX_STEP;
	conf->temperature_start = HK_BLIND_P9016_CLI_TEMPERATURE_START;
	conf->temperature_end = HK_BLIND_P9016_CLI_TEMPERATURE_END;
	conf->rho_train_start = HK_BLIND_P9016_CLI_RHO_TRAIN_START;
	conf->rho_train_end = HK_BLIND_P9016_CLI_RHO_TRAIN_END;
	conf->rho_train_mode = HK_BLIND_RHO_TRAIN_ENTROPY;
	conf->base_k_mode = HK_BLIND_BASE_K_NEIGHBOR_MEDIAN;
	conf->heldout_fraction = 0.0f;
	conf->heldout_seed = HK_BLIND_P9016_CLI_INIT_SEED;
	conf->write_raw = 0;
}

static void print_usage(FILE *fp)
{
	fprintf(fp, "Usage: hickit blind-p9016 [options]\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -i, --input FILE                 input P9016 .pairs[.gz] file [%s]\n",
			HK_BLIND_P9016_CLI_DEFAULT_INPUT);
	fprintf(fp, "  -o, --out-dir DIR                output directory [temporary directory]\n");
	fprintf(fp, "      --bd-iter INT                blind EM iterations [%d]\n",
			HK_BLIND_P9016_CLI_N_ITER);
	fprintf(fp, "      --bd-relax-steps INT         relax steps per iteration [%d]\n",
			HK_BLIND_P9016_CLI_RELAX_STEPS);
	fprintf(fp, "      --bd-relax-step FLOAT        relax step size [%.6g]\n",
			HK_BLIND_P9016_CLI_RELAX_STEP);
	fprintf(fp, "      --bd-rho-mode STR            constant|entropy|entropy_floor|entropy_cis_constant_trans|entropy_cis_floor_trans [entropy]\n");
	fprintf(fp, "      --bd-base-k-mode STR         uniform|neighbor_median [neighbor_median]\n");
	fprintf(fp, "      --bd-heldout-frac FLOAT      deterministic raw-contact holdout fraction [0]\n");
	fprintf(fp, "      --bd-heldout-seed INT        deterministic holdout salt [%llu]\n",
			(unsigned long long)HK_BLIND_P9016_CLI_INIT_SEED);
	fprintf(fp, "      --bd-write-raw               write raw posterior TSV\n");
	fprintf(fp, "      --help                       show this help\n");
}

static int parse_rho_mode(const char *s, int *mode)
{
	if (strcmp(s, "constant") == 0) *mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	else if (strcmp(s, "entropy") == 0) *mode = HK_BLIND_RHO_TRAIN_ENTROPY;
	else if (strcmp(s, "entropy_floor") == 0 || strcmp(s, "entropy_with_floor") == 0)
		*mode = HK_BLIND_RHO_TRAIN_ENTROPY_WITH_FLOOR;
	else if (strcmp(s, "entropy_cis_constant_trans") == 0)
		*mode = HK_BLIND_RHO_TRAIN_ENTROPY_CIS_CONSTANT_TRANS;
	else if (strcmp(s, "entropy_cis_floor_trans") == 0)
		*mode = HK_BLIND_RHO_TRAIN_ENTROPY_CIS_FLOOR_TRANS;
	else return -1;
	return 0;
}

static int parse_base_k_mode(const char *s, int *mode)
{
	if (strcmp(s, "uniform") == 0) *mode = HK_BLIND_BASE_K_UNIFORM;
	else if (strcmp(s, "neighbor_median") == 0) *mode = HK_BLIND_BASE_K_NEIGHBOR_MEDIAN;
	else return -1;
	return 0;
}

static int parse_positive_int(const char *s, int *out)
{
	char *end = 0;
	long v;
	errno = 0;
	v = strtol(s, &end, 10);
	if (errno != 0 || end == s || *end != 0 || v <= 0 || v > INT32_MAX)
		return -1;
	*out = (int)v;
	return 0;
}

static int parse_nonnegative_int(const char *s, int *out)
{
	char *end = 0;
	long v;
	errno = 0;
	v = strtol(s, &end, 10);
	if (errno != 0 || end == s || *end != 0 || v < 0 || v > INT32_MAX)
		return -1;
	*out = (int)v;
	return 0;
}

static int parse_float_range(const char *s, float lo, float hi, float *out)
{
	char *end = 0;
	float v;
	errno = 0;
	v = strtof(s, &end);
	if (errno != 0 || end == s || *end != 0 || !isfinite(v) || v < lo || v > hi)
		return -1;
	*out = v;
	return 0;
}

static int parse_u64_value(const char *s, uint64_t *out)
{
	char *end = 0;
	unsigned long long v;
	errno = 0;
	v = strtoull(s, &end, 10);
	if (errno != 0 || end == s || *end != 0)
		return -1;
	*out = (uint64_t)v;
	return 0;
}

static int parse_args(int argc, char *argv[], struct hk_blind_p9016_cli_conf *conf)
{
	enum {
		OPT_BD_ITER = 1000,
		OPT_BD_RELAX_STEPS,
		OPT_BD_RELAX_STEP,
		OPT_BD_RHO_MODE,
		OPT_BD_BASE_K_MODE,
		OPT_BD_HELDOUT_FRAC,
		OPT_BD_HELDOUT_SEED,
		OPT_BD_WRITE_RAW,
		OPT_HELP
	};
	static struct option long_options[] = {
		{ "input", required_argument, 0, 'i' },
		{ "out-dir", required_argument, 0, 'o' },
		{ "bd-iter", required_argument, 0, OPT_BD_ITER },
		{ "bd-relax-steps", required_argument, 0, OPT_BD_RELAX_STEPS },
		{ "bd-relax-step", required_argument, 0, OPT_BD_RELAX_STEP },
		{ "bd-rho-mode", required_argument, 0, OPT_BD_RHO_MODE },
		{ "bd-base-k-mode", required_argument, 0, OPT_BD_BASE_K_MODE },
		{ "bd-heldout-frac", required_argument, 0, OPT_BD_HELDOUT_FRAC },
		{ "bd-heldout-seed", required_argument, 0, OPT_BD_HELDOUT_SEED },
		{ "bd-write-raw", no_argument, 0, OPT_BD_WRITE_RAW },
		{ "help", no_argument, 0, OPT_HELP },
		{ 0, 0, 0, 0 }
	};
	int c, idx;

	optind = 1;
	while ((c = getopt_long(argc, argv, "i:o:", long_options, &idx)) >= 0) {
		switch (c) {
		case 'i':
			conf->input_path = optarg;
			break;
		case 'o':
			conf->out_dir = optarg;
			break;
		case OPT_BD_ITER:
			if (parse_positive_int(optarg, &conf->n_iter) != 0) return -1;
			break;
		case OPT_BD_RELAX_STEPS:
			if (parse_nonnegative_int(optarg, &conf->relax_steps) != 0) return -1;
			break;
		case OPT_BD_RELAX_STEP:
			if (parse_float_range(optarg, 0.0f, 1.0f, &conf->relax_step) != 0) return -1;
			break;
		case OPT_BD_RHO_MODE:
			if (parse_rho_mode(optarg, &conf->rho_train_mode) != 0) return -1;
			break;
		case OPT_BD_BASE_K_MODE:
			if (parse_base_k_mode(optarg, &conf->base_k_mode) != 0) return -1;
			break;
		case OPT_BD_HELDOUT_FRAC:
			if (parse_float_range(optarg, 0.0f, 0.95f, &conf->heldout_fraction) != 0) return -1;
			break;
		case OPT_BD_HELDOUT_SEED:
			if (parse_u64_value(optarg, &conf->heldout_seed) != 0) return -1;
			break;
		case OPT_BD_WRITE_RAW:
			conf->write_raw = 1;
			break;
		case OPT_HELP:
			print_usage(stdout);
			return 2;
		default:
			return -1;
		}
	}
	if (optind != argc)
		return -1;
	if (conf->n_iter <= 1) {
		conf->temperature_end = conf->temperature_start;
		conf->rho_train_end = conf->rho_train_start;
	}
	return 0;
}

static int path_exists_dir(const char *path)
{
	struct stat st;
	return stat(path, &st) == 0 && S_ISDIR(st.st_mode);
}

static int make_output_dir(const char *requested, char *dir, size_t dir_size)
{
	if (requested && *requested) {
		snprintf(dir, dir_size, "%s", requested);
		if (mkdir(dir, 0700) == 0 || (errno == EEXIST && path_exists_dir(dir)))
			return 0;
		return -1;
	}
	for (int i = 0; i < 100; ++i) {
		snprintf(dir, dir_size, "/tmp/hk_blind_p9016_cli_%ld_%ld_%d",
				 (long)getpid(), (long)time(0), i);
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

static uint64_t split_hash_pair(const struct hk_pair *p, uint64_t seed)
{
	uint64_t h = hash64(seed ^ 0x9e3779b97f4a7c15ULL);
	h ^= hash64(p->chr + 0x632be59bd9b4e019ULL);
	h ^= hash64((p->pos << 1) ^ (p->pos >> 63) ^ 0x85157af5ULL);
	h ^= hash64((uint64_t)(uint8_t)p->strand[0] << 8 | (uint8_t)p->strand[1]);
	return hash64(h);
}

static int split_pairs(const struct hk_blind_p9016_cli_conf *conf, const struct hk_map *m,
					   struct hk_pair **train_pairs_, int32_t *n_train_,
					   struct hk_blind_pair **heldout_raw_, int32_t *n_heldout_,
					   struct hk_blind_heldout_diag *heldout_diag)
{
	struct hk_pair *train_pairs = 0;
	struct hk_blind_pair *heldout_raw = 0;
	int32_t n_train = 0, n_heldout = 0;
	int32_t i;
	uint64_t threshold;

	assert(conf);
	assert(m && m->pairs);
	assert(train_pairs_ && n_train_ && heldout_raw_ && n_heldout_);
	assert(heldout_diag);
	hk_blind_heldout_diag_init(heldout_diag);
	heldout_diag->enabled = conf->heldout_fraction > 0.0f;
	heldout_diag->fraction = conf->heldout_fraction;
	heldout_diag->seed = conf->heldout_seed;
	heldout_diag->n_input_raw_total = m->n_pairs;

	train_pairs = (struct hk_pair*)calloc((size_t)m->n_pairs, sizeof(*train_pairs));
	heldout_raw = (struct hk_blind_pair*)calloc((size_t)m->n_pairs, sizeof(*heldout_raw));
	if (train_pairs == 0 || heldout_raw == 0) {
		free(train_pairs);
		free(heldout_raw);
		return -1;
	}

	threshold = (uint64_t)(conf->heldout_fraction * (float)HK_BLIND_P9016_CLI_HASH_MOD + 0.5f);
	for (i = 0; i < m->n_pairs; ++i) {
		int to_heldout = 0;
		if (threshold > 0) {
			uint64_t bucket = split_hash_pair(&m->pairs[i], conf->heldout_seed) % HK_BLIND_P9016_CLI_HASH_MOD;
			to_heldout = bucket < threshold;
		}
		if (to_heldout) {
			hk_blind_pair_from_pair(&heldout_raw[n_heldout++], &m->pairs[i]);
		} else {
			train_pairs[n_train++] = m->pairs[i];
		}
	}
	if (n_train <= 0) {
		free(train_pairs);
		free(heldout_raw);
		return -1;
	}

	heldout_diag->n_raw_train = n_train;
	heldout_diag->n_raw_heldout = n_heldout;
	*train_pairs_ = train_pairs;
	*n_train_ = n_train;
	*heldout_raw_ = heldout_raw;
	*n_heldout_ = n_heldout;
	return 0;
}

static void set_schedule_conf(const struct hk_blind_p9016_cli_conf *cli,
							  struct hk_blind_iter_schedule_conf *conf)
{
	conf->n_iter = cli->n_iter;
	conf->base_conf.unit = HK_BLIND_P9016_CLI_UNIT;
	conf->base_conf.d_scale = HK_BLIND_P9016_CLI_D_SCALE;
	conf->base_conf.base_k = HK_BLIND_P9016_CLI_LEGACY_BASE_K_UNUSED;
	conf->base_conf.temperature = cli->temperature_start;
	conf->base_conf.rho_train = cli->rho_train_start;
	conf->base_conf.min_sep_unit = HK_BLIND_P9016_CLI_MIN_SEP_UNIT;
	conf->base_conf.lambda_sep = HK_BLIND_P9016_CLI_LAMBDA_SEP;
	conf->base_conf.relax_step = cli->relax_step;
	conf->base_conf.relax_steps = cli->relax_steps;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.rho_train_mode = cli->rho_train_mode;
	conf->base_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->base_conf.d_scale_eps_count = HK_BLIND_P9016_CLI_D_SCALE_EPS_COUNT;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = cli->temperature_start;
	conf->temperature_end = cli->temperature_end;
	conf->rho_train_start = cli->rho_train_start;
	conf->rho_train_end = cli->rho_train_end;
}

static int init_prior(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
					  struct hk_blind_prior_diag *diag)
{
	return hk_blind_bpair_set_init_cis_inter_ratio_prior(bmap, set,
														 HK_BLIND_P9016_CLI_PRIOR_EPS,
														 diag);
}

static int init_coords(struct hk_bmap *bmap, fvec3_t *haploid, fvec3_t *diploid)
{
	struct hk_fdg_conf scaffold_conf;
	int ret;
	hk_fdg_conf_init(&scaffold_conf);
	scaffold_conf.backend = HK_FDG_BACKEND_CPU;
	scaffold_conf.n_iter = HK_BLIND_P9016_CLI_SCAFFOLD_FDG_N_ITER;
	ret = hk_blind_init_haploid_scaffold_from_bmap_fdg(bmap, &scaffold_conf, haploid,
													   HK_BLIND_P9016_CLI_INIT_SEED);
	if (ret != 0) return ret;
	return hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
													 diploid, HK_BLIND_P9016_CLI_INIT_EPS,
													 HK_BLIND_P9016_CLI_INIT_NOISE_SCALE,
													 HK_BLIND_P9016_CLI_INIT_SEED);
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
	FILE *fp = fopen(raw_path, "w");
	int ret;
	if (fp == 0) return -1;
	ret = hk_blind_write_raw_contact_posterior_tsv(fp, raw, n_raw, bmap, set);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int write_manifest(const char *manifest_path, const struct hk_blind_p9016_cli_conf *cli,
						  const char *out_dir, const char *posterior_path,
						  const char *coords_path, const char *diag_path,
						  const char *raw_path, const struct hk_bmap *bmap,
						  const struct hk_blind_bpair_set *set,
						  const struct hk_blind_prior_diag *prior_diag,
						  const struct hk_blind_base_k_stats *base_k_stats,
						  const struct hk_blind_iter_loop_diag *loop_diag,
						  const struct hk_blind_heldout_diag *heldout_diag)
{
	FILE *fp = fopen(manifest_path, "w");
	const char *status;
	if (fp == 0) return -1;
	status = (loop_diag->n_bad_iter == 0 &&
			  loop_diag->n_relax_nonfinite_iter == 0 &&
			  loop_diag->n_coord_nonfinite == 0)? "OK" : "WARN";
	if (fprintf(fp,
				"key\tvalue\n"
				"sample\tP9016\n"
				"runner_family\t%s\n"
				"runner_version\t%s\n"
				"default_profile\t%s\n"
				"input_path\t%s\n"
				"output_dir\t%s\n"
				"n_raw\t%d\n"
				"n_raw_input_total\t%lld\n"
				"n_raw_train\t%d\n"
				"n_raw_heldout\t%d\n"
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
				"rho_train_mode\t%s\n"
				"rho_train_floor\t%.9g\n"
				"contact_k_multiplier_cis\t1\n"
				"contact_k_multiplier_trans\t1\n"
				"d_scale_mode\traw_count\n"
				"d_scale_eps_count\t%.9g\n"
				"prior_mode\t%s\n"
				"prior_eps\t%.9g\n"
				"prior_inter_density\t%.17g\n"
				"prior_observed_inter\t%.17g\n"
				"prior_possible_inter\t%.17g\n"
				"prior_n_distance_bins\t%d\n"
				"prior_n_alpha\t%d\n"
				"prior_alpha_min\t%.9g\n"
				"prior_alpha_median\t%.9g\n"
				"prior_alpha_max\t%.9g\n"
				"prior_smoothing_method\tmoving_average_3bin_possible_weighted\n"
				"prior_alpha_clamp_min\t%.9g\n"
				"prior_alpha_clamp_max\t0.5\n"
				"same_bin_filter_enabled\t%d\n"
				"n_raw_same_bin_excluded\t%lld\n"
				"n_bpair_same_bin_excluded\t%d\n"
				"n_raw_cis\t%lld\n"
				"n_raw_trans\t%lld\n"
				"n_bpair_cis\t%d\n"
				"n_bpair_trans\t%d\n"
				"raw_posterior_same_bin_policy\tuniform_unknown_rows\n"
				"min_sep_unit\t%.9g\n"
				"lambda_sep\t%.9g\n"
				"relax_step\t%.9g\n"
				"relax_steps\t%d\n"
				"enable_repulsion\t1\n"
				"repulsion_mode\t%d\n"
				"repulsion_blocking_mode\tcurrent_edge_blocking\n"
				"temperature_start\t%.9g\n"
				"temperature_end\t%.9g\n"
				"rho_train_start\t%.9g\n"
				"rho_train_end\t%.9g\n"
				"init_mode\tunphased_scaffold_split\n"
				"init_scale\t%.9g\n"
				"init_eps_effective\t%.9g\n"
				"init_noise_scale_effective\t%.9g\n"
				"init_split_params_used\t1\n"
				"init_seed\t%llu\n"
				"scaffold_source\tunphased_bmap_fdg\n"
				"scaffold_fdg_n_iter\t%d\n"
				"heldout_enabled\t%d\n"
				"heldout_fraction\t%.9g\n"
				"heldout_seed\t%llu\n"
				"heldout_n_bpair\t%d\n"
				"heldout_n_bpair_eval\t%d\n"
				"heldout_n_raw_eval\t%lld\n"
				"heldout_mean_expected_energy\t%.17g\n"
				"heldout_mean_min_energy\t%.17g\n"
				"heldout_mean_entropy\t%.17g\n"
				"heldout_mean_pU\t%.17g\n"
				"heldout_mean_best_normalized_distance\t%.17g\n"
				"heldout_short_distance_frac\t%.17g\n"
				"write_raw_posterior\t%d\n"
				"output_bpair_posterior\t%s\n"
				"output_coords\t%s\n"
				"output_loop_diag\t%s\n",
				HK_BLIND_P9016_CLI_RUNNER_FAMILY,
				HK_BLIND_P9016_CLI_RUNNER_VERSION,
				HK_BLIND_P9016_CLI_DEFAULT_PROFILE,
				cli->input_path, out_dir, set->n_raw,
				(long long)heldout_diag->n_input_raw_total,
				heldout_diag->n_raw_train, heldout_diag->n_raw_heldout,
				set->n_bpairs, bmap->n_beads, HK_BLIND_P9016_CLI_RESOLUTION,
				loop_diag->n_completed, HK_BLIND_P9016_CLI_UNIT, HK_BLIND_P9016_CLI_D_SCALE,
				hk_blind_base_k_mode_name(cli->base_k_mode),
				base_k_stats->mean, base_k_stats->min, base_k_stats->mean, base_k_stats->max,
				base_k_stats->n_nonfinite, HK_BLIND_P9016_CLI_LEGACY_BASE_K_UNUSED,
				hk_blind_rho_train_mode_name(cli->rho_train_mode),
				HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
				HK_BLIND_P9016_CLI_D_SCALE_EPS_COUNT,
				hk_blind_prior_mode_name(prior_diag->prior_mode),
				prior_diag->eps_prior, prior_diag->inter_density,
				prior_diag->observed_inter, prior_diag->possible_inter,
				prior_diag->n_distance_bins, prior_diag->n_alpha,
				prior_diag->alpha_min, prior_diag->alpha_median, prior_diag->alpha_max,
				prior_diag->eps_prior,
				set->same_bin_filter_enabled, (long long)set->n_raw_same_bin_excluded,
				set->n_bpair_same_bin_excluded,
				(long long)set->n_raw_cis, (long long)set->n_raw_trans,
				set->n_bpair_cis, set->n_bpair_trans,
				HK_BLIND_P9016_CLI_MIN_SEP_UNIT, HK_BLIND_P9016_CLI_LAMBDA_SEP,
				cli->relax_step, cli->relax_steps, HK_BLIND_REPULSION_CELL,
				cli->temperature_start, cli->temperature_end,
				cli->rho_train_start, cli->rho_train_end,
				HK_BLIND_P9016_CLI_INIT_SCALE, HK_BLIND_P9016_CLI_INIT_EPS,
				HK_BLIND_P9016_CLI_INIT_NOISE_SCALE,
				(unsigned long long)HK_BLIND_P9016_CLI_INIT_SEED,
				HK_BLIND_P9016_CLI_SCAFFOLD_FDG_N_ITER,
				heldout_diag->enabled, heldout_diag->fraction,
				(unsigned long long)heldout_diag->seed,
				heldout_diag->n_bpair_heldout, heldout_diag->n_bpair_eval,
				(long long)heldout_diag->n_raw_eval,
				heldout_diag->mean_expected_energy, heldout_diag->mean_min_energy,
				heldout_diag->mean_entropy, heldout_diag->mean_pU,
				heldout_diag->mean_best_normalized_distance,
				heldout_diag->short_distance_frac,
				cli->write_raw, posterior_path, coords_path, diag_path) < 0) {
		fclose(fp);
		return -1;
	}
	if (cli->write_raw && fprintf(fp, "output_raw_posterior\t%s\n", raw_path) < 0) {
		fclose(fp);
		return -1;
	}
	if (fprintf(fp,
				"final_mean_entropy\t%.9g\n"
				"final_mean_pU\t%.9g\n"
				"final_mean_sep\t%.9g\n"
				"final_min_sep\t%.9g\n"
				"final_max_sep\t%.9g\n"
				"final_sum_wedge_k\t%.17g\n"
				"final_mean_rho_train_bpair\t%.9g\n"
				"final_min_rho_train_bpair\t%.9g\n"
				"final_max_rho_train_bpair\t%.9g\n"
				"final_n_skipped_same_bin_bpairs\t%lld\n"
				"final_repulsion_energy\t%.9g\n"
				"posterior_refreshed_after_final_relax\t%d\n"
				"posterior_refresh_temperature\t%.9g\n"
				"posterior_refresh_prior_mode\t%s\n"
				"posterior_refresh_mean_kl\t%.17g\n"
				"posterior_refresh_top_state_switch_frac\t%.9g\n"
				"posterior_refresh_mean_pU_before\t%.9g\n"
				"posterior_refresh_mean_pU_after\t%.9g\n"
				"n_bad_iter\t%d\n"
				"n_relax_nonfinite_iter\t%d\n"
				"n_coord_nonfinite\t%d\n"
				"status\t%s\n",
				loop_diag->final_mean_entropy, loop_diag->final_mean_pU,
				loop_diag->final_mean_sep, loop_diag->final_min_sep, loop_diag->final_max_sep,
				loop_diag->final_sum_wedge_k,
				loop_diag->final_mean_rho_train_bpair, loop_diag->final_min_rho_train_bpair,
				loop_diag->final_max_rho_train_bpair,
				(long long)loop_diag->final_n_skipped_same_bin_bpairs,
				loop_diag->final_repulsion_energy,
				loop_diag->posterior_refreshed_after_final_relax,
				loop_diag->posterior_refresh_temperature,
				hk_blind_prior_mode_name(prior_diag->prior_mode),
				loop_diag->posterior_refresh_mean_kl,
				loop_diag->posterior_refresh_top_state_switch_frac,
				loop_diag->posterior_refresh_mean_pU_before,
				loop_diag->posterior_refresh_mean_pU_after,
				loop_diag->n_bad_iter, loop_diag->n_relax_nonfinite_iter,
				loop_diag->n_coord_nonfinite, status) < 0) {
		fclose(fp);
		return -1;
	}
	if (fclose(fp) != 0) return -1;
	return 0;
}

int hk_blind_p9016_cli_main(int argc, char *argv[])
{
	struct hk_blind_p9016_cli_conf cli;
	struct hk_map *m = 0;
	struct hk_pair *train_pairs = 0;
	struct hk_blind_pair *raw = 0, *heldout_raw = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_bpair_set *set = 0, *heldout_set = 0;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	struct hk_blind_prior_diag prior_diag;
	struct hk_blind_base_k_stats base_k_stats;
	struct hk_blind_heldout_diag heldout_diag;
	struct hk_blind_single_iter_diag *per_iter = 0;
	fvec3_t *haploid = 0, *diploid = 0;
	char out_dir[1024], posterior_path[1200], coords_path[1200];
	char diag_path[1200], raw_path[1200], manifest_path[1200];
	int32_t n_train = 0, n_heldout = 0, n_diploid;
	int32_t i;
	int failed = 0;

	cli_conf_init(&cli);
	{
		int parse_ret = parse_args(argc, argv, &cli);
		if (parse_ret == 2) return 0;
		if (parse_ret != 0) {
			fprintf(stderr, "invalid blind-p9016 arguments\n");
			print_usage(stderr);
			return 1;
		}
	}

	hk_verbose = 0;
	m = hk_map_read(cli.input_path);
	if (m == 0 || m->pairs == 0 || m->n_pairs <= 0) {
		fprintf(stderr, "failed to read positive P9016 pairs input: %s\n", cli.input_path);
		failed = 1;
		goto cleanup;
	}
	if (split_pairs(&cli, m, &train_pairs, &n_train, &heldout_raw, &n_heldout,
					&heldout_diag) != 0) {
		fprintf(stderr, "failed to build deterministic train/heldout split\n");
		failed = 1;
		goto cleanup;
	}
	if (cli.base_k_mode == HK_BLIND_BASE_K_NEIGHBOR_MEDIAN)
		hk_pair_count_nei(n_train, train_pairs,
						  HK_BLIND_P9016_CLI_BASE_K_NEIGHBOR_RADIUS,
						  HK_BLIND_P9016_CLI_BASE_K_NEIGHBOR_RADIUS);
	bmap = hk_bmap_gen(m->d, n_train, train_pairs, HK_BLIND_P9016_CLI_RESOLUTION, 1);
	raw = (struct hk_blind_pair*)calloc((size_t)n_train, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate blind-p9016 bmap/raw arrays\n");
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n_train; ++i)
		hk_blind_pair_from_pair(&raw[i], &train_pairs[i]);

	set = hk_blind_bpair_set_build(bmap, n_train, raw);
	if (set == 0 || set->n_bpairs <= 0) {
		fprintf(stderr, "failed to build positive blind-p9016 bpair set\n");
		failed = 1;
		goto cleanup;
	}
	failed |= hk_blind_bpair_set_apply_base_k_mode(bmap, set, cli.base_k_mode) != 0;
	failed |= init_prior(bmap, set, &prior_diag) != 0;
	hk_blind_bpair_set_base_k_stats(set, &base_k_stats);
	if (failed) goto cleanup;

	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	haploid = (fvec3_t*)calloc((size_t)bmap->n_beads, sizeof(*haploid));
	diploid = (fvec3_t*)calloc((size_t)n_diploid, sizeof(*diploid));
	per_iter = (struct hk_blind_single_iter_diag*)calloc((size_t)cli.n_iter, sizeof(*per_iter));
	if (haploid == 0 || diploid == 0 || per_iter == 0) {
		fprintf(stderr, "failed to allocate blind-p9016 coordinate/diagnostic arrays\n");
		failed = 1;
		goto cleanup;
	}
	failed |= init_coords(bmap, haploid, diploid) != 0;
	if (failed) goto cleanup;

	hk_fdg_conf_init(&fdg_conf);
	fdg_conf.backend = HK_FDG_BACKEND_CPU;
	set_schedule_conf(&cli, &schedule_conf);
	failed |= hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &fdg_conf, diploid, 0,
												   &schedule_conf, per_iter, &loop_diag) != 0;
	if (loop_diag.n_bad_iter != 0 || loop_diag.n_relax_nonfinite_iter != 0 ||
		loop_diag.n_coord_nonfinite != 0) {
		fprintf(stderr, "blind-p9016 loop produced bad/nonfinite diagnostics\n");
		failed = 1;
	}
	if (failed) goto cleanup;

	if (n_heldout > 0) {
		heldout_set = hk_blind_bpair_set_build(bmap, n_heldout, heldout_raw);
		if (heldout_set == 0) {
			failed = 1;
			goto cleanup;
		}
		heldout_diag.n_bpair_heldout = heldout_set->n_bpairs;
		if (hk_blind_eval_heldout_bpair_set(bmap, heldout_set, &fdg_conf, diploid,
											 HK_BLIND_P9016_CLI_UNIT,
											 loop_diag.posterior_refresh_temperature > 0.0f?
												loop_diag.posterior_refresh_temperature : cli.temperature_end,
											 &heldout_diag) != 0) {
			failed = 1;
			goto cleanup;
		}
	}

	if (make_output_dir(cli.out_dir, out_dir, sizeof(out_dir)) != 0) {
		fprintf(stderr, "failed to create output directory\n");
		failed = 1;
		goto cleanup;
	}
	path_join(posterior_path, sizeof(posterior_path), out_dir, "p9016_full.bpair_posterior.tsv");
	path_join(coords_path, sizeof(coords_path), out_dir, "p9016_full.coords.tsv");
	path_join(diag_path, sizeof(diag_path), out_dir, "p9016_full.loop_diag.tsv");
	path_join(raw_path, sizeof(raw_path), out_dir, "p9016_full.raw_posterior.tsv");
	path_join(manifest_path, sizeof(manifest_path), out_dir, "p9016_full.manifest.tsv");

	failed |= write_required_outputs(posterior_path, coords_path, diag_path,
									 bmap, set, diploid, &loop_diag) != 0;
	if (cli.write_raw)
		failed |= write_optional_raw_output(raw_path, raw, n_train, bmap, set) != 0;
	failed |= write_manifest(manifest_path, &cli, out_dir, posterior_path, coords_path,
							 diag_path, raw_path, bmap, set, &prior_diag, &base_k_stats,
							 &loop_diag, &heldout_diag) != 0;
	if (failed) goto cleanup;

	fprintf(stderr,
			"blind-p9016: outdir=%s n_train=%d n_heldout=%d n_bpair=%d "
			"base_k_mode=%s rho_mode=%s status=OK\n",
			out_dir, n_train, n_heldout, set->n_bpairs,
			hk_blind_base_k_mode_name(cli.base_k_mode),
			hk_blind_rho_train_mode_name(cli.rho_train_mode));

cleanup:
	hk_blind_bpair_set_destroy(heldout_set);
	hk_blind_bpair_set_destroy(set);
	free(haploid);
	free(diploid);
	free(per_iter);
	free(raw);
	free(heldout_raw);
	free(train_pairs);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
