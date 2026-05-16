#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include "hkpriv.h"

#define HK_BLIND_P9016_FULL_PATH "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_FULL_RESOLUTION 1000000
#define HK_BLIND_P9016_FULL_N_ITER 3
#define HK_BLIND_P9016_FULL_UNIT 1.0f
#define HK_BLIND_P9016_FULL_D_SCALE 1.0f
#define HK_BLIND_P9016_FULL_LEGACY_BASE_K_UNUSED 2.0f
#define HK_BLIND_P9016_FULL_TEMPERATURE_START 2.0f
#define HK_BLIND_P9016_FULL_TEMPERATURE_END 1.0f
#define HK_BLIND_P9016_FULL_RHO_TRAIN_START 1.0f
#define HK_BLIND_P9016_FULL_RHO_TRAIN_END 1.0f
#define HK_BLIND_P9016_FULL_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_FULL_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_FULL_STEP 0.001f
#define HK_BLIND_P9016_FULL_RELAX_STEPS 5
#define HK_BLIND_P9016_FULL_INIT_EPS 0.5f
#define HK_BLIND_P9016_FULL_INIT_NOISE_SCALE 0.0f
#define HK_BLIND_P9016_FULL_INIT_SEED 17ULL
#define HK_BLIND_P9016_FULL_INIT_SCALE 1.0f
#define HK_BLIND_P9016_FULL_INIT_MODE HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT
#define HK_BLIND_P9016_FULL_PRIOR_MODE HK_BLIND_PRIOR_CIS_INTER_RATIO
#define HK_BLIND_P9016_FULL_PRIOR_EPS 1e-6f
#define HK_BLIND_P9016_FULL_RHO_TRAIN_MODE HK_BLIND_RHO_TRAIN_ENTROPY
#define HK_BLIND_P9016_FULL_D_SCALE_MODE HK_BLIND_D_SCALE_RAW_COUNT
#define HK_BLIND_P9016_FULL_D_SCALE_EPS_COUNT 1e-6f
#define HK_BLIND_P9016_FULL_SCAFFOLD_FDG_N_ITER 50
#define HK_BLIND_P9016_FULL_PREFIX_SCAN 8192
#define HK_BLIND_P9016_FULL_BASE_K_MODE HK_BLIND_BASE_K_NEIGHBOR_MEDIAN
#define HK_BLIND_P9016_FULL_BASE_K_NEIGHBOR_RADIUS 10000000
#define HK_BLIND_P9016_FULL_RUNNER_FAMILY "full_cpu"
#define HK_BLIND_P9016_FULL_RUNNER_VERSION "2026-04-30"
#define HK_BLIND_P9016_FULL_DEFAULT_PROFILE "p9016_1mb_pairs_only_v1"
#define HK_BLIND_P9016_FULL_HELDOUT_HASH_MOD 1000000ULL

struct p9016_full_cli_conf {
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
	int prior_mode;
	int base_k_mode;
	float heldout_fraction;
	int write_raw;
};

struct p9016_heldout_eval {
	int enabled;
	int64_t n_input_raw_total;
	int32_t n_train_raw;
	int32_t n_heldout_raw;
	int32_t n_eval_raw;
	float fraction;
	double mean_min_dist;
	double mean_min_energy;
	double mean_softmin_energy;
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

static void p9016_full_cli_conf_init(struct p9016_full_cli_conf *conf)
{
	memset(conf, 0, sizeof(*conf));
	conf->input_path = HK_BLIND_P9016_FULL_PATH;
	conf->out_dir = 0;
	conf->n_iter = HK_BLIND_P9016_FULL_N_ITER;
	conf->relax_steps = HK_BLIND_P9016_FULL_RELAX_STEPS;
	conf->relax_step = HK_BLIND_P9016_FULL_STEP;
	conf->temperature_start = HK_BLIND_P9016_FULL_TEMPERATURE_START;
	conf->temperature_end = HK_BLIND_P9016_FULL_TEMPERATURE_END;
	conf->rho_train_start = HK_BLIND_P9016_FULL_RHO_TRAIN_START;
	conf->rho_train_end = HK_BLIND_P9016_FULL_RHO_TRAIN_END;
	conf->rho_train_mode = HK_BLIND_P9016_FULL_RHO_TRAIN_MODE;
	conf->prior_mode = HK_BLIND_P9016_FULL_PRIOR_MODE;
	conf->base_k_mode = HK_BLIND_P9016_FULL_BASE_K_MODE;
	conf->heldout_fraction = 0.0f;
	conf->write_raw = 0;
}

static void print_usage(FILE *fp)
{
	fprintf(fp,
			"usage: run_blind_p9016_full_cpu.bin [options]\n"
			"  --input PATH                 P9016 pairs input (default: %s)\n"
			"  --out-dir DIR                output directory (default: create /tmp dir)\n"
			"  --bd-iter N                  EM iterations (default: %d)\n"
			"  --bd-relax-steps N           relax steps per iteration (default: %d)\n"
			"  --bd-relax-step X            Euler relax step (default: %.9g)\n"
			"  --bd-temp-start X            starting temperature (default: %.9g)\n"
			"  --bd-temp-end X              ending temperature (default: %.9g)\n"
			"  --bd-rho-start X             starting rho_train (default: %.9g)\n"
			"  --bd-rho-end X               ending rho_train (default: %.9g)\n"
			"  --bd-rho-mode MODE           constant|entropy|entropy_with_floor|entropy_cis_constant_trans|entropy_cis_floor_trans\n"
			"  --bd-prior-mode MODE         uniform|cis_inter_ratio\n"
			"  --bd-base-k-mode MODE        uniform|neighbor_median (default: %s)\n"
			"  --heldout-frac X             deterministic raw-contact heldout fraction in [0,1)\n"
			"  --write-raw                  write raw posterior\n"
			"  --no-write-raw               do not write raw posterior\n"
			"  --help                       show this help\n",
			HK_BLIND_P9016_FULL_PATH, HK_BLIND_P9016_FULL_N_ITER,
			HK_BLIND_P9016_FULL_RELAX_STEPS, HK_BLIND_P9016_FULL_STEP,
			HK_BLIND_P9016_FULL_TEMPERATURE_START, HK_BLIND_P9016_FULL_TEMPERATURE_END,
			HK_BLIND_P9016_FULL_RHO_TRAIN_START, HK_BLIND_P9016_FULL_RHO_TRAIN_END,
			hk_blind_base_k_mode_name(HK_BLIND_P9016_FULL_BASE_K_MODE));
}

static int parse_rho_train_mode(const char *s)
{
	if (strcmp(s, "constant") == 0) return HK_BLIND_RHO_TRAIN_CONSTANT;
	if (strcmp(s, "entropy") == 0) return HK_BLIND_RHO_TRAIN_ENTROPY;
	if (strcmp(s, "entropy_with_floor") == 0 || strcmp(s, "entropy_floor") == 0)
		return HK_BLIND_RHO_TRAIN_ENTROPY_WITH_FLOOR;
	if (strcmp(s, "entropy_cis_constant_trans") == 0) return HK_BLIND_RHO_TRAIN_ENTROPY_CIS_CONSTANT_TRANS;
	if (strcmp(s, "entropy_cis_floor_trans") == 0) return HK_BLIND_RHO_TRAIN_ENTROPY_CIS_FLOOR_TRANS;
	return -1;
}

static int parse_prior_mode(const char *s)
{
	if (strcmp(s, "uniform") == 0) return HK_BLIND_PRIOR_UNIFORM;
	if (strcmp(s, "cis_inter_ratio") == 0) return HK_BLIND_PRIOR_CIS_INTER_RATIO;
	return -1;
}

static int parse_base_k_mode(const char *s)
{
	if (strcmp(s, "uniform") == 0) return HK_BLIND_BASE_K_UNIFORM;
	if (strcmp(s, "neighbor_median") == 0 || strcmp(s, "neighbor") == 0)
		return HK_BLIND_BASE_K_NEIGHBOR_MEDIAN;
	return -1;
}

static int parse_cli_int(const char *label, const char *s, int min_value, int *out)
{
	char *end = 0;
	long v;
	errno = 0;
	v = strtol(s, &end, 10);
	if (errno != 0 || end == s || *end != 0 || v < min_value || v > INT32_MAX) {
		fprintf(stderr, "ERROR: bad %s: %s\n", label, s);
		return -1;
	}
	*out = (int)v;
	return 0;
}

static int parse_cli_float(const char *label, const char *s, float min_value, float max_value, float *out)
{
	char *end = 0;
	float v;
	errno = 0;
	v = strtof(s, &end);
	if (errno != 0 || end == s || *end != 0 || !isfinite(v) || v < min_value || v > max_value) {
		fprintf(stderr, "ERROR: bad %s: %s\n", label, s);
		return -1;
	}
	*out = v;
	return 0;
}

static int parse_full_cli(int argc, char **argv, struct p9016_full_cli_conf *conf)
{
	int i;
	p9016_full_cli_conf_init(conf);
	if (getenv("HK_BLIND_WRITE_RAW") != 0 && strcmp(getenv("HK_BLIND_WRITE_RAW"), "1") == 0)
		conf->write_raw = 1;
	for (i = 1; i < argc; ++i) {
		const char *arg = argv[i];
		const char *val = 0;
		if (strcmp(arg, "--help") == 0 || strcmp(arg, "-h") == 0) {
			print_usage(stdout);
			return 1;
		}
		if (strcmp(arg, "--write-raw") == 0) {
			conf->write_raw = 1;
			continue;
		}
		if (strcmp(arg, "--no-write-raw") == 0) {
			conf->write_raw = 0;
			continue;
		}
		if (i + 1 >= argc) {
			fprintf(stderr, "ERROR: option needs a value: %s\n", arg);
			return -1;
		}
		val = argv[++i];
		if (strcmp(arg, "--input") == 0) conf->input_path = val;
		else if (strcmp(arg, "--out-dir") == 0) conf->out_dir = val;
		else if (strcmp(arg, "--bd-iter") == 0) {
			if (parse_cli_int(arg, val, 1, &conf->n_iter) != 0) return -1;
		} else if (strcmp(arg, "--bd-relax-steps") == 0) {
			if (parse_cli_int(arg, val, 0, &conf->relax_steps) != 0) return -1;
		} else if (strcmp(arg, "--bd-relax-step") == 0) {
			if (parse_cli_float(arg, val, 0.0f, 1.0f, &conf->relax_step) != 0) return -1;
		} else if (strcmp(arg, "--bd-temp-start") == 0) {
			if (parse_cli_float(arg, val, 1e-6f, 1e6f, &conf->temperature_start) != 0) return -1;
		} else if (strcmp(arg, "--bd-temp-end") == 0) {
			if (parse_cli_float(arg, val, 1e-6f, 1e6f, &conf->temperature_end) != 0) return -1;
		} else if (strcmp(arg, "--bd-rho-start") == 0) {
			if (parse_cli_float(arg, val, 0.0f, 1e6f, &conf->rho_train_start) != 0) return -1;
		} else if (strcmp(arg, "--bd-rho-end") == 0) {
			if (parse_cli_float(arg, val, 0.0f, 1e6f, &conf->rho_train_end) != 0) return -1;
		} else if (strcmp(arg, "--bd-rho-mode") == 0) {
			conf->rho_train_mode = parse_rho_train_mode(val);
			if (!hk_blind_rho_train_mode_valid(conf->rho_train_mode)) {
				fprintf(stderr, "ERROR: bad --bd-rho-mode: %s\n", val);
				return -1;
			}
		} else if (strcmp(arg, "--bd-prior-mode") == 0) {
			conf->prior_mode = parse_prior_mode(val);
			if (!hk_blind_prior_mode_valid(conf->prior_mode)) {
				fprintf(stderr, "ERROR: bad --bd-prior-mode: %s\n", val);
				return -1;
			}
		} else if (strcmp(arg, "--bd-base-k-mode") == 0) {
			conf->base_k_mode = parse_base_k_mode(val);
			if (!hk_blind_base_k_mode_valid(conf->base_k_mode)) {
				fprintf(stderr, "ERROR: bad --bd-base-k-mode: %s\n", val);
				return -1;
			}
		} else if (strcmp(arg, "--heldout-frac") == 0) {
			if (parse_cli_float(arg, val, 0.0f, 0.999999f, &conf->heldout_fraction) != 0) return -1;
		} else {
			fprintf(stderr, "ERROR: unrecognized option: %s\n", arg);
			return -1;
		}
	}
	if (conf->n_iter == 1) {
		conf->temperature_end = conf->temperature_start;
		conf->rho_train_end = conf->rho_train_start;
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

static int check_no_forbidden_strings(const char *label, const char *s)
{
	int failed = 0;
	char buf[128];
	const char *forbidden[] = { "phase0", "phase1", "truth", "oracle" };
	int i;
	for (i = 0; i < 4; ++i) {
		snprintf(buf, sizeof(buf), "%s no %s", label, forbidden[i]);
		failed |= check_true(buf, strstr(s, forbidden[i]) == 0);
	}
	return failed;
}

static int make_output_dir(char *dir, size_t dir_size)
{
	long pid = (long)getpid();
	long stamp = (long)time(0);
	int i;

	for (i = 0; i < 100; ++i) {
		snprintf(dir, dir_size, "/tmp/hk_blind_p9016_full_cpu_%ld_%ld_%d", pid, stamp, i);
		if (mkdir(dir, 0700) == 0)
			return 0;
		if (errno != EEXIST)
			return -1;
	}
	return -1;
}

static int prepare_output_dir(const char *requested, char *dir, size_t dir_size)
{
	if (requested && requested[0]) {
		snprintf(dir, dir_size, "%s", requested);
		if (mkdir(dir, 0700) == 0 || errno == EEXIST)
			return 0;
		return -1;
	}
	return make_output_dir(dir, dir_size);
}

static const char *scaffold_source_for_init_mode(int init_mode)
{
	switch (init_mode) {
	case HK_BLIND_INIT_TOY_SPLIT: return "toy_deterministic";
	case HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT: return "unphased_bmap_fdg";
	case HK_BLIND_INIT_RANDOM_DIPLOID: return "random_diploid_direct";
	case HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT: return "random_haploid";
	default: return "unknown";
	}
}

static int init_coords_for_mode(struct hk_bmap *bmap, const struct hk_fdg_conf *scaffold_conf,
								int init_mode, fvec3_t *haploid, fvec3_t *diploid)
{
	int ret = -1;
	assert(bmap);
	assert(diploid);
	assert(hk_blind_init_mode_valid(init_mode));
	if (init_mode == HK_BLIND_INIT_RANDOM_DIPLOID) {
		return hk_blind_init_random_diploid_coords(bmap, diploid,
												   HK_BLIND_P9016_FULL_INIT_SCALE,
												   HK_BLIND_P9016_FULL_INIT_SEED);
	}
	assert(haploid);
	if (init_mode == HK_BLIND_INIT_TOY_SPLIT) {
		ret = hk_blind_init_toy_haploid_scaffold(bmap, haploid);
	} else if (init_mode == HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT) {
		ret = hk_blind_init_haploid_scaffold_from_bmap_fdg(bmap, scaffold_conf, haploid,
														   HK_BLIND_P9016_FULL_INIT_SEED);
	} else if (init_mode == HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT) {
		ret = hk_blind_init_random_haploid_scaffold(bmap, haploid,
													HK_BLIND_P9016_FULL_INIT_SCALE,
													HK_BLIND_P9016_FULL_INIT_SEED);
	}
	if (ret != 0) return ret;
	return hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
													 diploid, HK_BLIND_P9016_FULL_INIT_EPS,
													 HK_BLIND_P9016_FULL_INIT_NOISE_SCALE,
													 HK_BLIND_P9016_FULL_INIT_SEED);
}

static int check_coords_finite(const fvec3_t *coords, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("full coord finite", isfinite(coords[i][a]));
	return failed;
}

static void set_schedule_conf(struct hk_blind_iter_schedule_conf *conf,
							  const struct p9016_full_cli_conf *cli)
{
	conf->n_iter = cli->n_iter;
	conf->base_conf.unit = HK_BLIND_P9016_FULL_UNIT;
	conf->base_conf.d_scale = HK_BLIND_P9016_FULL_D_SCALE;
	conf->base_conf.base_k = HK_BLIND_P9016_FULL_LEGACY_BASE_K_UNUSED;
	conf->base_conf.temperature = cli->temperature_start;
	conf->base_conf.rho_train = cli->rho_train_start;
	conf->base_conf.min_sep_unit = HK_BLIND_P9016_FULL_MIN_SEP_UNIT;
	conf->base_conf.lambda_sep = HK_BLIND_P9016_FULL_LAMBDA_SEP;
	conf->base_conf.relax_step = cli->relax_step;
	conf->base_conf.relax_steps = cli->relax_steps;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.rho_train_mode = cli->rho_train_mode;
	conf->base_conf.d_scale_mode = HK_BLIND_P9016_FULL_D_SCALE_MODE;
	conf->base_conf.d_scale_eps_count = HK_BLIND_P9016_FULL_D_SCALE_EPS_COUNT;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = cli->temperature_start;
	conf->temperature_end = cli->temperature_end;
	conf->rho_train_start = cli->rho_train_start;
	conf->rho_train_end = cli->rho_train_end;
}

static int init_prior_for_mode(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
							   int prior_mode, struct hk_blind_prior_diag *diag)
{
	assert(diag);
	assert(hk_blind_prior_mode_valid(prior_mode));
	memset(diag, 0, sizeof(*diag));
	diag->prior_mode = prior_mode;
	diag->eps_prior = HK_BLIND_P9016_FULL_PRIOR_EPS;
	if (prior_mode == HK_BLIND_PRIOR_UNIFORM) {
		hk_blind_bpair_set_init_uniform_prior(set);
		return 0;
	}
	return hk_blind_bpair_set_init_cis_inter_ratio_prior(bmap, set,
														 HK_BLIND_P9016_FULL_PRIOR_EPS,
														 diag);
}

static uint64_t p9016_raw_contact_hash(const struct hk_pair *p, int32_t raw_id)
{
	uint64_t h;
	assert(p);
	h = hash64(p->chr);
	h ^= hash64(p->pos + 0x9e3779b97f4a7c15ULL);
	h ^= hash64((uint64_t)(uint8_t)(p->strand[0] + 2) << 8 | (uint8_t)(p->strand[1] + 2));
	h ^= hash64((uint64_t)(uint32_t)raw_id + 0xbf58476d1ce4e5b9ULL);
	return hash64(h);
}

static int p9016_is_heldout_contact(const struct hk_pair *p, int32_t raw_id, float heldout_fraction)
{
	uint64_t bucket;
	uint64_t cutoff;
	if (heldout_fraction <= 0.0f) return 0;
	bucket = p9016_raw_contact_hash(p, raw_id) % HK_BLIND_P9016_FULL_HELDOUT_HASH_MOD;
	cutoff = (uint64_t)((double)heldout_fraction * (double)HK_BLIND_P9016_FULL_HELDOUT_HASH_MOD + 0.5);
	if (cutoff >= HK_BLIND_P9016_FULL_HELDOUT_HASH_MOD)
		cutoff = HK_BLIND_P9016_FULL_HELDOUT_HASH_MOD - 1;
	return bucket < cutoff;
}

static int split_train_heldout_contacts(const struct hk_map *m, float heldout_fraction,
										struct hk_pair **train_pairs_out, int32_t *n_train_out,
										struct hk_blind_pair **heldout_raw_out,
										int32_t *n_heldout_out)
{
	struct hk_pair *train_pairs = 0;
	struct hk_blind_pair *heldout_raw = 0;
	int32_t i, n_train = 0, n_heldout = 0;

	assert(m);
	assert(train_pairs_out);
	assert(n_train_out);
	assert(heldout_raw_out);
	assert(n_heldout_out);
	*train_pairs_out = 0;
	*n_train_out = 0;
	*heldout_raw_out = 0;
	*n_heldout_out = 0;
	if (m->n_pairs <= 0) return -1;

	train_pairs = (struct hk_pair*)calloc((size_t)m->n_pairs, sizeof(*train_pairs));
	if (heldout_fraction > 0.0f)
		heldout_raw = (struct hk_blind_pair*)calloc((size_t)m->n_pairs, sizeof(*heldout_raw));
	if (train_pairs == 0 || (heldout_fraction > 0.0f && heldout_raw == 0)) {
		free(train_pairs);
		free(heldout_raw);
		return -1;
	}

	for (i = 0; i < m->n_pairs; ++i) {
		if (p9016_is_heldout_contact(&m->pairs[i], i, heldout_fraction)) {
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
	*train_pairs_out = train_pairs;
	*n_train_out = n_train;
	*heldout_raw_out = heldout_raw;
	*n_heldout_out = n_heldout;
	return 0;
}

static float p9016_min4(float a, float b, float c, float d)
{
	float m = a < b? a : b;
	m = m < c? m : c;
	return m < d? m : d;
}

static int evaluate_heldout_contacts(const struct hk_bmap *bmap, const struct hk_fdg_conf *fdg_conf,
									 const fvec3_t *coords, const struct hk_blind_pair *heldout_raw,
									 int32_t n_heldout, struct p9016_heldout_eval *eval)
{
	int32_t i;
	double sum_min_dist = 0.0, sum_min_energy = 0.0, sum_softmin_energy = 0.0;
	const float softmin_temperature = 1.0f;

	assert(eval);
	if (!eval->enabled || n_heldout <= 0 || heldout_raw == 0) {
		eval->n_eval_raw = 0;
		eval->mean_min_dist = 0.0;
		eval->mean_min_energy = 0.0;
		eval->mean_softmin_energy = 0.0;
		return 0;
	}
	for (i = 0; i < n_heldout; ++i) {
		int32_t bid[2];
		int32_t i0, i1, j0, j1;
		float d[HK_BLIND_N_STATE], e[HK_BLIND_N_STATE];
		float min_d, min_e, max_score;
		double z;
		int s;
		hk_blind_pair_to_bids(bmap, &heldout_raw[i], bid);
		i0 = hk_diploid_bid(bid[0], HK_DIPLOID_COPY0);
		i1 = hk_diploid_bid(bid[0], HK_DIPLOID_COPY1);
		j0 = hk_diploid_bid(bid[1], HK_DIPLOID_COPY0);
		j1 = hk_diploid_bid(bid[1], HK_DIPLOID_COPY1);
		d[HK_BLIND_STATE_00] = sqrtf((coords[i0][0] - coords[j0][0]) * (coords[i0][0] - coords[j0][0]) +
									  (coords[i0][1] - coords[j0][1]) * (coords[i0][1] - coords[j0][1]) +
									  (coords[i0][2] - coords[j0][2]) * (coords[i0][2] - coords[j0][2]));
		d[HK_BLIND_STATE_01] = sqrtf((coords[i0][0] - coords[j1][0]) * (coords[i0][0] - coords[j1][0]) +
									  (coords[i0][1] - coords[j1][1]) * (coords[i0][1] - coords[j1][1]) +
									  (coords[i0][2] - coords[j1][2]) * (coords[i0][2] - coords[j1][2]));
		d[HK_BLIND_STATE_10] = sqrtf((coords[i1][0] - coords[j0][0]) * (coords[i1][0] - coords[j0][0]) +
									  (coords[i1][1] - coords[j0][1]) * (coords[i1][1] - coords[j0][1]) +
									  (coords[i1][2] - coords[j0][2]) * (coords[i1][2] - coords[j0][2]));
		d[HK_BLIND_STATE_11] = sqrtf((coords[i1][0] - coords[j1][0]) * (coords[i1][0] - coords[j1][0]) +
									  (coords[i1][1] - coords[j1][1]) * (coords[i1][1] - coords[j1][1]) +
									  (coords[i1][2] - coords[j1][2]) * (coords[i1][2] - coords[j1][2]));
		for (s = 0; s < HK_BLIND_N_STATE; ++s)
			e[s] = hk_fdg_contact_energy_dist(fdg_conf, d[s], HK_BLIND_P9016_FULL_UNIT, 1.0f, 1.0f);
		min_d = p9016_min4(d[0], d[1], d[2], d[3]);
		min_e = p9016_min4(e[0], e[1], e[2], e[3]);
		max_score = -e[0] / softmin_temperature;
		for (s = 1; s < HK_BLIND_N_STATE; ++s) {
			float score = -e[s] / softmin_temperature;
			if (score > max_score) max_score = score;
		}
		z = 0.0;
		for (s = 0; s < HK_BLIND_N_STATE; ++s)
			z += exp((double)(-e[s] / softmin_temperature - max_score));
		sum_min_dist += min_d;
		sum_min_energy += min_e;
		sum_softmin_energy += -(double)softmin_temperature * (log(z / (double)HK_BLIND_N_STATE) + max_score);
	}
	eval->n_eval_raw = n_heldout;
	eval->mean_min_dist = sum_min_dist / n_heldout;
	eval->mean_min_energy = sum_min_energy / n_heldout;
	eval->mean_softmin_energy = sum_softmin_energy / n_heldout;
	return 0;
}

static int write_heldout_validation_tsv(const char *path, const struct p9016_heldout_eval *eval)
{
	FILE *fp;
	if (!eval->enabled) return 0;
	fp = fopen(path, "w");
	if (fp == 0) return -1;
	if (fprintf(fp,
				"metric\tvalue\n"
				"heldout_fraction\t%.9g\n"
				"n_input_raw_total\t%lld\n"
				"n_train_raw\t%d\n"
				"n_heldout_raw\t%d\n"
				"n_eval_raw\t%d\n"
				"mean_min_dist\t%.17g\n"
				"mean_min_energy\t%.17g\n"
				"mean_softmin_energy\t%.17g\n",
				eval->fraction, (long long)eval->n_input_raw_total, eval->n_train_raw,
				eval->n_heldout_raw, eval->n_eval_raw, eval->mean_min_dist,
				eval->mean_min_energy, eval->mean_softmin_energy) < 0) {
		fclose(fp);
		return -1;
	}
	if (fclose(fp) != 0) return -1;
	return 0;
}

static int check_binned_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set)
{
	int failed = 0;
	int32_t i;

	failed |= check_true("full set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_true("full raw positive", set->n_raw > 0);
	failed |= check_true("full bpair positive", set->n_bpairs > 0);
	for (i = 0; i < set->n_raw; ++i) {
		failed |= check_true("full raw2binned id lower", set->raw2binned[i].bpair_id >= 0);
		failed |= check_true("full raw2binned id upper", set->raw2binned[i].bpair_id < set->n_bpairs);
		failed |= check_true("full raw2binned swapped", set->raw2binned[i].swapped == 0 || set->raw2binned[i].swapped == 1);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		failed |= check_true("full key sorted", bp->key.bid[0] <= bp->key.bid[1]);
		failed |= check_true("full key bid0 valid", bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		failed |= check_true("full key bid1 valid", bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
	}
	return failed;
}

static int check_loop_health(const struct hk_blind_iter_loop_diag *diag,
							 const struct hk_blind_single_iter_diag *per_iter,
							 const struct p9016_full_cli_conf *cli)
{
	int failed = 0;
	int32_t i;

	failed |= check_i32("full loop n_iter", diag->n_iter, cli->n_iter);
	failed |= check_i32("full loop completed", diag->n_completed, cli->n_iter);
	failed |= check_i32("full loop bad iter", diag->n_bad_iter, 0);
	failed |= check_i32("full loop relax bad", diag->n_relax_nonfinite_iter, 0);
	failed |= check_i32("full loop coord bad", diag->n_coord_nonfinite, 0);
	failed |= check_close("full initial temperature", diag->initial_temperature, cli->temperature_start);
	failed |= check_close("full final temperature", diag->final_temperature, cli->temperature_end);
	failed |= check_close("full initial rho", diag->initial_rho_train, cli->rho_train_start);
	failed |= check_close("full final rho", diag->final_rho_train, cli->rho_train_end);
	failed |= check_true("full final entropy finite", isfinite(diag->final_mean_entropy));
	failed |= check_true("full final entropy range", diag->final_mean_entropy >= -1e-6f &&
						 diag->final_mean_entropy <= logf((float)HK_BLIND_N_STATE) + 1e-6f);
	failed |= check_true("full final pU finite", isfinite(diag->final_mean_pU));
	failed |= check_true("full final pU range", diag->final_mean_pU >= -1e-6f && diag->final_mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("full final sep finite", isfinite(diag->final_mean_sep));
	failed |= check_true("full final sum k finite", isfinite(diag->final_sum_wedge_k));
	failed |= check_true("full final sum k nonnegative", diag->final_sum_wedge_k >= 0.0);
	for (i = 0; i < cli->n_iter; ++i) {
		const struct hk_blind_iter_diag *pre = &per_iter[i].pre_relax_diag;
		failed |= check_i32("full posterior nonfinite", pre->n_posterior_nonfinite, 0);
		failed |= check_i32("full posterior bad sum", pre->n_posterior_bad_sum, 0);
		failed |= check_i32("full posterior out of range", pre->n_posterior_out_of_range, 0);
		failed |= check_i32("full uncertainty nonfinite", pre->n_uncertainty_nonfinite, 0);
		failed |= check_i32("full uncertainty out of range", pre->n_uncertainty_out_of_range, 0);
		failed |= check_i32("full five-state bad sum", pre->n_five_state_bad_sum, 0);
		failed |= check_i32("full wedge nonfinite", pre->n_wedge_nonfinite, 0);
		failed |= check_i32("full wedge bad k", pre->n_wedge_bad_k, 0);
		failed |= check_i32("full wedge bad d_scale", pre->n_wedge_bad_d_scale, 0);
		failed |= check_i32("full sep force nonfinite", pre->sep_force_nonfinite, 0);
		failed |= check_i32("full relax completed", per_iter[i].relax_diag.n_completed, cli->relax_steps);
		failed |= check_i32("full relax nonfinite step", per_iter[i].relax_diag.n_nonfinite_step, 0);
		failed |= check_i32("full relax coord bad", per_iter[i].relax_diag.n_coord_nonfinite, 0);
		failed |= check_i32("full backbone bad", per_iter[i].relax_diag.n_backbone_nonfinite_step, 0);
		failed |= check_i32("full repulsion bad", per_iter[i].relax_diag.n_repulsion_nonfinite_step, 0);
		failed |= check_true("full relax total finite", isfinite(per_iter[i].relax_diag.final_total_energy));
		failed |= check_true("full repulsion finite", isfinite(per_iter[i].relax_diag.final_repulsion_energy));
	}
	return failed;
}

static int validate_final_bpair_set(const struct hk_blind_bpair_set *set)
{
	struct hk_blind_iter_diag diag;
	int failed = 0;

	hk_blind_iter_diag_init(&diag);
	hk_blind_iter_diag_validate_bpair_set(set, &diag);
	failed |= check_i32("full final posterior nonfinite", diag.n_posterior_nonfinite, 0);
	failed |= check_i32("full final posterior bad sum", diag.n_posterior_bad_sum, 0);
	failed |= check_i32("full final posterior out of range", diag.n_posterior_out_of_range, 0);
	failed |= check_i32("full final uncertainty nonfinite", diag.n_uncertainty_nonfinite, 0);
	failed |= check_i32("full final uncertainty out of range", diag.n_uncertainty_out_of_range, 0);
	failed |= check_i32("full final five-state bad sum", diag.n_five_state_bad_sum, 0);
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
						  const char *heldout_path,
						  const struct p9016_full_cli_conf *cli,
						  const struct p9016_heldout_eval *heldout_eval,
						  const struct hk_bmap *bmap,
						  const struct hk_blind_bpair_set *set,
						  const struct hk_blind_prior_diag *prior_diag,
						  const struct hk_blind_base_k_stats *base_k_stats,
						  const struct hk_blind_iter_loop_diag *loop_diag)
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
				"n_input_raw_total\t%lld\n"
				"n_train_raw\t%d\n"
				"n_heldout_raw\t%d\n"
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
				"rho_train_mode\t%s\n"
				"d_scale_mode\t%s\n"
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
				"prior_smoothing_method\t%s\n"
				"prior_alpha_clamp_min\t%.9g\n"
				"prior_alpha_clamp_max\t%.9g\n"
				"same_bin_filter_enabled\t%d\n"
				"n_raw_same_bin_excluded\t%lld\n"
				"n_bpair_same_bin_excluded\t%d\n"
				"raw_posterior_same_bin_policy\t%s\n"
				"min_sep_unit\t%.9g\n"
				"lambda_sep\t%.9g\n"
				"relax_step\t%.9g\n"
				"relax_steps\t%d\n"
				"enable_repulsion\t%d\n"
				"repulsion_mode\t%d\n"
				"repulsion_blocking_mode\t%s\n"
				"temperature_start\t%.9g\n"
				"temperature_end\t%.9g\n"
				"rho_train_start\t%.9g\n"
				"rho_train_end\t%.9g\n"
				"init_mode\t%s\n"
				"init_scale\t%.9g\n"
				"init_eps_effective\t%.9g\n"
				"init_noise_scale_effective\t%.9g\n"
				"init_split_params_used\t%d\n"
				"init_seed\t%llu\n"
				"scaffold_source\t%s\n"
				"scaffold_fdg_n_iter\t%d\n"
				"write_raw_posterior\t%d\n"
				"heldout_enabled\t%d\n"
				"heldout_fraction\t%.9g\n"
				"heldout_split_key\t%s\n"
				"heldout_n_eval_raw\t%d\n"
				"heldout_mean_min_dist\t%.17g\n"
				"heldout_mean_min_energy\t%.17g\n"
				"heldout_mean_softmin_energy\t%.17g\n"
				"output_bpair_posterior\t%s\n"
				"output_coords\t%s\n"
				"output_loop_diag\t%s\n"
				,
				HK_BLIND_P9016_FULL_RUNNER_FAMILY,
				HK_BLIND_P9016_FULL_RUNNER_VERSION,
				HK_BLIND_P9016_FULL_DEFAULT_PROFILE,
				cli->input_path, out_dir, (long long)heldout_eval->n_input_raw_total,
				heldout_eval->n_train_raw, heldout_eval->n_heldout_raw,
				set->n_raw, set->n_bpairs, bmap->n_beads,
				HK_BLIND_P9016_FULL_RESOLUTION, cli->n_iter,
				HK_BLIND_P9016_FULL_UNIT, HK_BLIND_P9016_FULL_D_SCALE,
				hk_blind_base_k_mode_name(cli->base_k_mode),
				base_k_stats->mean, base_k_stats->min, base_k_stats->mean, base_k_stats->max,
				base_k_stats->n_nonfinite, HK_BLIND_P9016_FULL_LEGACY_BASE_K_UNUSED,
				hk_blind_rho_train_mode_name(cli->rho_train_mode),
				hk_blind_d_scale_mode_name(HK_BLIND_P9016_FULL_D_SCALE_MODE),
				HK_BLIND_P9016_FULL_D_SCALE_EPS_COUNT,
				hk_blind_prior_mode_name(prior_diag->prior_mode),
				prior_diag->eps_prior, prior_diag->inter_density,
				prior_diag->observed_inter, prior_diag->possible_inter,
				prior_diag->n_distance_bins, prior_diag->n_alpha,
				prior_diag->alpha_min, prior_diag->alpha_median, prior_diag->alpha_max,
				prior_diag->prior_mode == HK_BLIND_PRIOR_CIS_INTER_RATIO?
					"moving_average_3bin_possible_weighted" : "none",
				prior_diag->eps_prior, 0.5,
				set->same_bin_filter_enabled, (long long)set->n_raw_same_bin_excluded,
				set->n_bpair_same_bin_excluded, "uniform_unknown_rows",
				HK_BLIND_P9016_FULL_MIN_SEP_UNIT,
				HK_BLIND_P9016_FULL_LAMBDA_SEP, cli->relax_step,
				cli->relax_steps, 1, HK_BLIND_REPULSION_CELL,
				"current_edge_blocking",
				cli->temperature_start,
				cli->temperature_end, cli->rho_train_start,
				cli->rho_train_end,
				hk_blind_init_mode_name(HK_BLIND_P9016_FULL_INIT_MODE),
				HK_BLIND_P9016_FULL_INIT_SCALE,
				HK_BLIND_P9016_FULL_INIT_MODE == HK_BLIND_INIT_RANDOM_DIPLOID? 0.0 : HK_BLIND_P9016_FULL_INIT_EPS,
				HK_BLIND_P9016_FULL_INIT_MODE == HK_BLIND_INIT_RANDOM_DIPLOID? 0.0 : HK_BLIND_P9016_FULL_INIT_NOISE_SCALE,
				HK_BLIND_P9016_FULL_INIT_MODE == HK_BLIND_INIT_RANDOM_DIPLOID? 0 : 1,
				(unsigned long long)HK_BLIND_P9016_FULL_INIT_SEED,
				scaffold_source_for_init_mode(HK_BLIND_P9016_FULL_INIT_MODE),
				HK_BLIND_P9016_FULL_SCAFFOLD_FDG_N_ITER, cli->write_raw,
				heldout_eval->enabled, heldout_eval->fraction,
				heldout_eval->enabled? "raw_contact_fields_and_index" : "none",
				heldout_eval->n_eval_raw, heldout_eval->mean_min_dist,
				heldout_eval->mean_min_energy, heldout_eval->mean_softmin_energy,
				posterior_path, coords_path, diag_path) < 0) {
		fclose(fp);
		return -1;
	}
	if (cli->write_raw && fprintf(fp, "output_raw_posterior\t%s\n", raw_path) < 0) {
		fclose(fp);
		return -1;
	}
	if (heldout_eval->enabled && fprintf(fp, "output_heldout_validation\t%s\n", heldout_path) < 0) {
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
	failed |= check_true("full output size positive", file_size_or_negative(path) > 0);
	prefix = read_file_prefix(path, HK_BLIND_P9016_FULL_PREFIX_SCAN);
	failed |= check_true("full output prefix readable", prefix != 0);
	if (prefix == 0) return 1;
	failed |= check_true("full output header token", strstr(prefix, header_token) != 0);
	failed |= check_no_forbidden_strings(label, prefix);
	free(prefix);
	return failed;
}

static int validate_manifest_prefix(const char *path, const struct p9016_full_cli_conf *cli,
									const struct p9016_heldout_eval *heldout_eval)
{
	char *prefix = read_file_prefix(path, HK_BLIND_P9016_FULL_PREFIX_SCAN);
	int failed = 0;

	failed |= check_true("full manifest exists", file_exists(path));
	failed |= check_true("full manifest size positive", file_size_or_negative(path) > 0);
	failed |= check_true("full manifest readable", prefix != 0);
	if (prefix == 0) return 1;
	failed |= check_true("full manifest sample", strstr(prefix, "sample\tP9016") != 0);
	failed |= check_true("full manifest runner family",
						 strstr(prefix, "runner_family\t" HK_BLIND_P9016_FULL_RUNNER_FAMILY) != 0);
	failed |= check_true("full manifest runner version",
						 strstr(prefix, "runner_version\t" HK_BLIND_P9016_FULL_RUNNER_VERSION) != 0);
	failed |= check_true("full manifest default profile",
						 strstr(prefix, "default_profile\t" HK_BLIND_P9016_FULL_DEFAULT_PROFILE) != 0);
	failed |= check_true("full manifest n_raw", strstr(prefix, "n_raw") != 0);
	failed |= check_true("full manifest n input raw total", strstr(prefix, "n_input_raw_total") != 0);
	failed |= check_true("full manifest n train raw", strstr(prefix, "n_train_raw") != 0);
	failed |= check_true("full manifest n_bpair", strstr(prefix, "n_bpair") != 0);
	failed |= check_true("full manifest n_beads", strstr(prefix, "n_beads") != 0);
	failed |= check_true("full manifest init mode", strstr(prefix, "init_mode\tunphased_scaffold_split") != 0);
	failed |= check_true("full manifest scaffold source", strstr(prefix, "scaffold_source\tunphased_bmap_fdg") != 0);
	failed |= check_true("full manifest prior mode", strstr(prefix, hk_blind_prior_mode_name(cli->prior_mode)) != 0);
	failed |= check_true("full manifest rho train mode", strstr(prefix, hk_blind_rho_train_mode_name(cli->rho_train_mode)) != 0);
	failed |= check_true("full manifest d scale mode", strstr(prefix, "d_scale_mode\traw_count") != 0);
	failed |= check_true("full manifest base k effective", strstr(prefix, "base_k_effective") != 0);
	failed |= check_true("full manifest base k mode",
						 strstr(prefix, hk_blind_base_k_mode_name(cli->base_k_mode)) != 0);
	failed |= check_true("full manifest legacy base k", strstr(prefix, "legacy_base_k_unused") != 0);
	failed |= check_true("full manifest same-bin filter", strstr(prefix, "same_bin_filter_enabled\t1") != 0);
	failed |= check_true("full manifest posterior refresh",
						 strstr(prefix, "posterior_refreshed_after_final_relax\t1") != 0);
	failed |= check_true("full manifest status OK", strstr(prefix, "status\tOK") != 0);
	failed |= check_true("full manifest raw flag",
						 cli->write_raw? strstr(prefix, "write_raw_posterior\t1") != 0
									   : strstr(prefix, "write_raw_posterior\t0") != 0);
	failed |= check_true("full manifest heldout flag",
						 heldout_eval->enabled? strstr(prefix, "heldout_enabled\t1") != 0
											  : strstr(prefix, "heldout_enabled\t0") != 0);
	failed |= check_no_forbidden_strings("full manifest", prefix);
	free(prefix);
	return failed;
}

int hk_blind_p9016_full_cpu_main(int argc, char **argv)
{
	struct p9016_full_cli_conf cli;
	struct p9016_heldout_eval heldout_eval;
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_pair *train_pairs = 0;
	struct hk_blind_pair *raw = 0;
	struct hk_blind_pair *heldout_raw = 0;
	struct hk_blind_bpair_set *set = 0;
	struct hk_fdg_conf fdg_conf;
	struct hk_fdg_conf scaffold_conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	struct hk_blind_prior_diag prior_diag;
	struct hk_blind_base_k_stats base_k_stats;
	struct hk_blind_single_iter_diag *per_iter = 0;
	fvec3_t *haploid = 0, *diploid = 0;
	char out_dir[512] = {0};
	char posterior_path[768], coords_path[768], diag_path[768], raw_path[768], heldout_path[768], manifest_path[768];
	int32_t n_input_raw, n_train_raw = 0, n_heldout_raw = 0, n_diploid;
	int32_t i;
	int failed = 0;
	int parse_ret;
	clock_t t0 = clock(), t1;
	double elapsed_sec;

	memset(&heldout_eval, 0, sizeof(heldout_eval));
	parse_ret = parse_full_cli(argc, argv, &cli);
	if (parse_ret > 0) return 0;
	if (parse_ret < 0) {
		print_usage(stderr);
		return 1;
	}

	if (!file_exists(cli.input_path)) {
		fprintf(stderr, "ERROR: %s not found; full P9016 blind CPU run did not start\n",
				cli.input_path);
		return 1;
	}

	hk_verbose = 0;
	fprintf(stderr, "P9016 full CPU blind run: reading %s\n", cli.input_path);
	m = hk_map_read(cli.input_path);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", cli.input_path);
		return 1;
	}
	n_input_raw = m->n_pairs;
	if (n_input_raw <= 0) {
		fprintf(stderr, "full input raw positive: predicate failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_i32("full split train/heldout",
						split_train_heldout_contacts(m, cli.heldout_fraction,
													 &train_pairs, &n_train_raw,
													 &heldout_raw, &n_heldout_raw), 0);
	if (failed) goto cleanup;
	heldout_eval.enabled = cli.heldout_fraction > 0.0f;
	heldout_eval.fraction = cli.heldout_fraction;
	heldout_eval.n_input_raw_total = n_input_raw;
	heldout_eval.n_train_raw = n_train_raw;
	heldout_eval.n_heldout_raw = n_heldout_raw;

	fprintf(stderr, "P9016 full CPU blind run: building 1Mb blind contacts from %d train raw pairs",
			n_train_raw);
	if (n_heldout_raw > 0)
		fprintf(stderr, " (%d held out)", n_heldout_raw);
	fprintf(stderr, "\n");
	if (cli.base_k_mode == HK_BLIND_BASE_K_NEIGHBOR_MEDIAN)
		hk_pair_count_nei(n_train_raw, train_pairs,
						  HK_BLIND_P9016_FULL_BASE_K_NEIGHBOR_RADIUS,
						  HK_BLIND_P9016_FULL_BASE_K_NEIGHBOR_RADIUS);
	bmap = hk_bmap_gen(m->d, n_train_raw, train_pairs, HK_BLIND_P9016_FULL_RESOLUTION, 1);
	raw = (struct hk_blind_pair*)calloc((size_t)n_train_raw, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate full P9016 bmap/raw data\n");
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n_train_raw; ++i)
		hk_blind_pair_from_pair(&raw[i], &train_pairs[i]);

	set = hk_blind_bpair_set_build(bmap, n_train_raw, raw);
	failed |= check_binned_set(bmap, set);
	failed |= check_true("full bead count positive", bmap && bmap->n_beads > 0);
	failed |= check_i32("full base_k mode",
						hk_blind_bpair_set_apply_base_k_mode(bmap, set,
															 cli.base_k_mode), 0);
	failed |= check_i32("full prior init",
						init_prior_for_mode(bmap, set, cli.prior_mode, &prior_diag), 0);
	hk_blind_bpair_set_base_k_stats(set, &base_k_stats);
	if (failed) goto cleanup;

	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	haploid = (fvec3_t*)calloc(bmap->n_beads, sizeof(*haploid));
	diploid = (fvec3_t*)calloc(n_diploid, sizeof(*diploid));
	per_iter = (struct hk_blind_single_iter_diag*)calloc((size_t)cli.n_iter, sizeof(*per_iter));
	if (haploid == 0 || diploid == 0 || per_iter == 0) {
		fprintf(stderr, "failed to allocate full P9016 coordinate/diagnostic data\n");
		failed = 1;
		goto cleanup;
	}

	hk_fdg_conf_init(&scaffold_conf);
	scaffold_conf.backend = HK_FDG_BACKEND_CPU;
	scaffold_conf.n_iter = HK_BLIND_P9016_FULL_SCAFFOLD_FDG_N_ITER;
	failed |= check_i32("full init coords",
						init_coords_for_mode(bmap, &scaffold_conf, HK_BLIND_P9016_FULL_INIT_MODE,
											 haploid, diploid), 0);
	failed |= check_coords_finite(diploid, n_diploid);
	if (failed) goto cleanup;

	fprintf(stderr, "P9016 full CPU blind run: scheduled loop n_iter=%d relax_steps=%d\n",
			cli.n_iter, cli.relax_steps);
	hk_fdg_conf_init(&fdg_conf);
	fdg_conf.backend = HK_FDG_BACKEND_CPU;
	set_schedule_conf(&schedule_conf, &cli);
	failed |= check_i32("full scheduled ret",
						hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &fdg_conf, diploid, 0,
															 &schedule_conf, per_iter, &loop_diag), 0);
	failed |= check_loop_health(&loop_diag, per_iter, &cli);
	failed |= check_coords_finite(diploid, n_diploid);
	failed |= validate_final_bpair_set(set);
	if (failed) goto cleanup;

	if (heldout_eval.enabled)
		failed |= check_i32("full heldout eval",
							evaluate_heldout_contacts(bmap, &fdg_conf, diploid, heldout_raw,
													  n_heldout_raw, &heldout_eval), 0);
	if (failed) goto cleanup;

	failed |= check_i32("full output dir", prepare_output_dir(cli.out_dir, out_dir, sizeof(out_dir)), 0);
	if (failed) goto cleanup;
	snprintf(posterior_path, sizeof(posterior_path), "%s/p9016_full.bpair_posterior.tsv", out_dir);
	snprintf(coords_path, sizeof(coords_path), "%s/p9016_full.coords.tsv", out_dir);
	snprintf(diag_path, sizeof(diag_path), "%s/p9016_full.loop_diag.tsv", out_dir);
	snprintf(raw_path, sizeof(raw_path), "%s/p9016_full.raw_posterior.tsv", out_dir);
	snprintf(heldout_path, sizeof(heldout_path), "%s/p9016_full.heldout_contact_validation.tsv", out_dir);
	snprintf(manifest_path, sizeof(manifest_path), "%s/p9016_full.manifest.tsv", out_dir);

	fprintf(stderr, "P9016 full CPU blind run: writing outputs to %s\n", out_dir);
	failed |= check_i32("full write required outputs",
						write_required_outputs(posterior_path, coords_path, diag_path,
											   bmap, set, diploid, &loop_diag), 0);
	if (cli.write_raw) {
		fprintf(stderr, "P9016 full CPU blind run: writing raw posterior\n");
		failed |= check_i32("full write raw output",
							write_optional_raw_output(raw_path, raw, n_train_raw, bmap, set), 0);
	}
	if (heldout_eval.enabled)
		failed |= check_i32("full write heldout validation",
							write_heldout_validation_tsv(heldout_path, &heldout_eval), 0);
	failed |= check_i32("full write manifest",
						write_manifest(manifest_path, out_dir, posterior_path, coords_path,
									   diag_path, raw_path, heldout_path, &cli, &heldout_eval,
									   bmap, set,
									   &prior_diag, &base_k_stats, &loop_diag), 0);
	if (failed) goto cleanup;

	failed |= validate_file_prefix("full bpair posterior", posterior_path, "p00");
	failed |= validate_file_prefix("full coords", coords_path, "diploid_bid");
	failed |= validate_file_prefix("full loop diag", diag_path, "final_mean_entropy");
	failed |= validate_manifest_prefix(manifest_path, &cli, &heldout_eval);
	if (cli.write_raw)
		failed |= validate_file_prefix("full raw posterior", raw_path, "raw_id");
	if (heldout_eval.enabled)
		failed |= validate_file_prefix("full heldout validation", heldout_path, "mean_min_dist");
	if (failed) goto cleanup;

	t1 = clock();
	elapsed_sec = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
	fprintf(stderr,
			"P9016 full CPU blind run: n_raw=%d n_bpair=%d n_beads=%d "
			"final_entropy=%.8g final_pU=%.8g final_sum_wedge_k=%.8g "
			"heldout=%d outdir=%s elapsed_cpu_sec=%.3f status=OK\n",
			set->n_raw, set->n_bpairs, bmap->n_beads,
			loop_diag.final_mean_entropy, loop_diag.final_mean_pU,
			loop_diag.final_sum_wedge_k, n_heldout_raw, out_dir, elapsed_sec);

cleanup:
	if (failed)
		fprintf(stderr, "P9016 full CPU blind run failed%s%s\n",
				out_dir[0]? "; partial outputs remain in " : "",
				out_dir[0]? out_dir : "");
	free(haploid);
	free(diploid);
	free(per_iter);
	hk_blind_bpair_set_destroy(set);
	free(heldout_raw);
	free(raw);
	free(train_pairs);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}

int main(int argc, char **argv)
{
	return hk_blind_p9016_full_cpu_main(argc, argv);
}
