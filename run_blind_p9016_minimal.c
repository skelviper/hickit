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
#include "fdg_gpu.h"

#define HK_P9016_DEFAULT_PAIRS "/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz"
#define HK_P9016_DEFAULT_BIN_SIZE_BP 1000000
#define HK_P9016_DEFAULT_N_ITER 100
#define HK_P9016_DEFAULT_RELAX_STEPS 100
#define HK_P9016_DEFAULT_RELAX_STEP 0.012f
#define HK_P9016_UNIT 1.0f
#define HK_P9016_D_SCALE 1.0f
#define HK_P9016_D_SCALE_EPS_COUNT 1e-6f
#define HK_P9016_PRIOR_EPS 1e-6f
#define HK_P9016_PRIOR_ALPHA_CLAMP_MIN 1e-6f
#define HK_P9016_BASE_K_NEIGHBOR_RADIUS 10000000
#define HK_P9016_INIT_EPS 0.5f
#define HK_P9016_INIT_NOISE_SCALE 0.0f
#define HK_P9016_INIT_SCALE 10.0f
#define HK_P9016_INIT_SEED 17ULL
#define HK_P9016_SCAFFOLD_FDG_N_ITER 50
#define HK_P9016_CONFIG_NAME "minimal_soft_sep_off"
#define HK_P9016_REPULSION_MULTIPLIER 1.0f
#define HK_P9016_DEFAULT_CHAIN "4000000,1000000,200000"
#define HK_P9016_DEFAULT_CHILD_OFFSET_STEP 0.05f
#define HK_P9016_DEFAULT_PARENT_ANCHOR_K 0.05f
#define HK_P9016_MAX_STAGES 8

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

struct minimal_result {
	char output_dir[1024];
	char stage_name[32];
	char init_mode[128];
	int32_t stage_index;
	int32_t n_stages;
	int32_t bin_size_bp;
	int32_t parent_bin_size_bp;
	char init_source[128];
	char relax_backend[32];
	int32_t n_raw;
	int32_t n_bpair;
	int32_t n_beads;
	int64_t n_raw_cis;
	int64_t n_raw_trans;
	int32_t n_bpair_cis;
	int32_t n_bpair_trans;
	int64_t n_same_bin_excluded;
	float final_mean_entropy;
	float final_mean_pU;
	float final_mean_sep;
	float final_min_sep;
	float final_max_sep;
	double final_sum_wedge_k;
	double final_refreshed_sum_wedge_k;
	int64_t final_refreshed_n_wedges;
	float final_mean_rho_train_bpair;
	float final_min_rho_train_bpair;
	float final_max_rho_train_bpair;
	int64_t final_n_expanded_edges;
	int64_t final_n_softall_candidate_pairs;
	int64_t final_n_softall_filter1_pairs;
	int64_t final_n_softall_filter2_pairs;
	int64_t final_n_softall_bmap_pairs;
	int64_t final_n_softall_selected_raw;
	int64_t final_n_softall_gate_skip_raw;
	int64_t final_n_softall_same_bin_skip_raw;
	int64_t final_n_softall_state_raw_count[HK_BLIND_N_STATE];
	float final_contact_energy;
	float final_repulsion_energy;
	float final_backbone_energy;
	int64_t final_n_repulsion_pairs_considered;
	int64_t final_n_repulsion_pairs_blocked;
	int64_t final_n_repulsion_pairs_active;
	float k_rel_rep_effective;
	int posterior_refreshed_after_final_relax;
	int32_t n_bad_iter;
	int32_t n_relax_nonfinite_iter;
	int32_t n_coord_nonfinite;
	int status_ok;
};

struct stage_state {
	struct hk_map *map;
	struct hk_bmap *bmap;
	struct hk_blind_pair *raw;
	fvec3_t *coords;
	int32_t n_raw;
	int bin_size_bp;
};

struct hk_p9016_minimal_params {
	float init_eps;
	float init_noise_scale;
	float init_scale;
	float min_sep_unit;
	float lambda_sep;
	int d_scale_mode;
	float d_scale_eps_count;
	uint64_t init_seed;
	char config_name[128];
};

struct hk_p9016_sep_distribution {
	int32_t n_beads;
	int32_t n_sep_below_min;
	float min_sep_unit;
	float lambda_sep;
	float sep_min;
	float sep_p01;
	float sep_p05;
	float sep_p10;
	float sep_p25;
	float sep_median;
	float sep_mean;
	float sep_p75;
	float sep_p90;
	float sep_p95;
	float sep_p99;
	float sep_max;
	float frac_sep_below_min;
};

static int check_true(const char *label, int ok)
{
	if (!ok) {
		fprintf(stderr, "check failed: %s\n", label);
		return 1;
	}
	return 0;
}

static int check_i32(const char *label, int32_t got, int32_t expected)
{
	if (got != expected) {
		fprintf(stderr, "check failed: %s got %d expected %d\n", label, got, expected);
		return 1;
	}
	return 0;
}

static const char *env_or_default(const char *name, const char *fallback)
{
	const char *s = getenv(name);
	return s && s[0]? s : fallback;
}

static int env_int_or_default(const char *name, int fallback, int min_value)
{
	const char *s = getenv(name);
	char *end = 0;
	long v;
	if (s == 0 || s[0] == 0)
		return fallback;
	errno = 0;
	v = strtol(s, &end, 10);
	if (errno || end == s || *end != 0 || v < min_value || v > INT32_MAX) {
		fprintf(stderr, "invalid %s=%s; using %d\n", name, s, fallback);
		return fallback;
	}
	return (int)v;
}

static uint64_t env_u64_or_default(const char *name, uint64_t fallback)
{
	const char *s = getenv(name);
	char *end = 0;
	unsigned long long v;
	if (s == 0 || s[0] == 0)
		return fallback;
	errno = 0;
	v = strtoull(s, &end, 10);
	if (errno || end == s || *end != 0) {
		fprintf(stderr, "invalid %s=%s; using %llu\n",
				name, s, (unsigned long long)fallback);
		return fallback;
	}
	return (uint64_t)v;
}

static float env_float_or_default(const char *name, float fallback, float min_value)
{
	const char *s = getenv(name);
	char *end = 0;
	float v;
	if (s == 0 || s[0] == 0)
		return fallback;
	errno = 0;
	v = strtof(s, &end);
	if (errno || end == s || *end != 0 || !isfinite(v) || v < min_value) {
		fprintf(stderr, "invalid %s=%s; using %.9g\n", name, s, fallback);
		return fallback;
	}
	return v;
}

static int env_float_checked(const char *name, float fallback, float min_value, float *out)
{
	const char *s = getenv(name);
	char *end = 0;
	float v;
	assert(out);
	if (s == 0 || s[0] == 0) {
		*out = fallback;
		return 0;
	}
	errno = 0;
	v = strtof(s, &end);
	if (errno || end == s || *end != 0 || !isfinite(v) || v < min_value) {
		fprintf(stderr, "invalid %s=%s; expected finite value >= %.9g\n",
				name, s, min_value);
		return -1;
	}
	*out = v;
	return 0;
}

static int config_name_is_safe_slug(const char *s)
{
	const unsigned char *p;
	if (s == 0 || s[0] == 0 || strcmp(s, ".") == 0 || strcmp(s, "..") == 0)
		return 0;
	if (strstr(s, "..") != 0)
		return 0;
	for (p = (const unsigned char*)s; *p; ++p) {
		if ((*p >= 'A' && *p <= 'Z') || (*p >= 'a' && *p <= 'z') ||
			(*p >= '0' && *p <= '9') || *p == '_' || *p == '-' || *p == '.')
			continue;
		return 0;
	}
	return 1;
}

static int load_minimal_params(struct hk_p9016_minimal_params *p)
{
	const char *config_name;
	const char *config_env;
	int n;
	assert(p);
	memset(p, 0, sizeof(*p));
	if (env_float_checked("HK_BLIND_P9016_INIT_EPS",
						  HK_P9016_INIT_EPS, 0.0f, &p->init_eps) != 0 ||
		env_float_checked("HK_BLIND_P9016_INIT_NOISE_SCALE",
						  HK_P9016_INIT_NOISE_SCALE, 0.0f,
						  &p->init_noise_scale) != 0 ||
		env_float_checked("HK_BLIND_P9016_INIT_SCALE",
						  HK_P9016_INIT_SCALE, 0.0f, &p->init_scale) != 0 ||
		env_float_checked("HK_BLIND_P9016_MIN_SEP_UNIT",
						  0.0f, 0.0f, &p->min_sep_unit) != 0 ||
		env_float_checked("HK_BLIND_P9016_LAMBDA_SEP",
						  0.0f, 0.0f, &p->lambda_sep) != 0 ||
		env_float_checked("HK_BLIND_P9016_D_SCALE_EPS_COUNT",
						  HK_P9016_D_SCALE_EPS_COUNT, 0.0f,
						  &p->d_scale_eps_count) != 0)
		return -1;
	if (p->d_scale_eps_count <= 0.0f) {
		fprintf(stderr, "HK_BLIND_P9016_D_SCALE_EPS_COUNT must be > 0\n");
		return -1;
	}
	p->init_seed = env_u64_or_default("HK_BLIND_P9016_INIT_SEED",
									  HK_P9016_INIT_SEED);
	config_env = getenv("HK_BLIND_P9016_D_SCALE_MODE");
	if (config_env == 0 || config_env[0] == 0 || strcmp(config_env, "raw_count") == 0)
		p->d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	else if (strcmp(config_env, "expected_count") == 0)
		p->d_scale_mode = HK_BLIND_D_SCALE_EXPECTED_COUNT;
	else {
		fprintf(stderr, "invalid HK_BLIND_P9016_D_SCALE_MODE=%s; expected raw_count or expected_count\n",
				config_env);
		return -1;
	}
	config_env = getenv("HK_BLIND_P9016_CONFIG_NAME");
	config_name = config_env? config_env : HK_P9016_CONFIG_NAME;
	if (!config_name_is_safe_slug(config_name)) {
		fprintf(stderr, "HK_BLIND_P9016_CONFIG_NAME must be a non-empty safe slug using [A-Za-z0-9._-], without '..'\n");
		return -1;
	}
	n = snprintf(p->config_name, sizeof(p->config_name), "%s", config_name);
	if (n < 0) {
		fprintf(stderr, "failed to copy HK_BLIND_P9016_CONFIG_NAME\n");
		return -1;
	}
	if ((size_t)n >= sizeof(p->config_name)) {
		fprintf(stderr, "HK_BLIND_P9016_CONFIG_NAME is too long; max length is %zu\n",
				sizeof(p->config_name) - 1);
		return -1;
	}
	return 0;
}

static int parse_init_mode(int *ok)
{
	const char *s = getenv("HK_BLIND_P9016_INIT_MODE");
	if (ok) *ok = 1;
	if (s == 0 || s[0] == 0 || strcmp(s, "unphased_scaffold_split") == 0)
		return HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT;
	if (strcmp(s, "random_diploid") == 0)
		return HK_BLIND_INIT_RANDOM_DIPLOID;
	if (strcmp(s, "random_haploid_split") == 0)
		return HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT;
	if (strcmp(s, "toy_split") == 0)
		return HK_BLIND_INIT_TOY_SPLIT;
	fprintf(stderr, "invalid HK_BLIND_P9016_INIT_MODE=%s\n", s);
	if (ok) *ok = 0;
	return HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT;
}

static int env_flag_enabled(const char *name)
{
	const char *s = getenv(name);
	return s && s[0] && strcmp(s, "0") != 0 &&
		strcmp(s, "false") != 0 && strcmp(s, "FALSE") != 0;
}

static int parse_relax_backend(void)
{
	const char *s = getenv("HK_BLIND_P9016_RELAX_BACKEND");
	if (s == 0 || s[0] == 0 || strcmp(s, "cpu") == 0)
		return HK_FDG_BACKEND_CPU;
	if (strcmp(s, "gpu") == 0)
		return HK_FDG_BACKEND_GPU;
	if (strcmp(s, "auto") == 0)
		return HK_FDG_BACKEND_AUTO;
	fprintf(stderr, "invalid HK_BLIND_P9016_RELAX_BACKEND=%s; using cpu\n", s);
	return HK_FDG_BACKEND_CPU;
}

static const char *relax_backend_name(int backend)
{
	if (backend == HK_FDG_BACKEND_GPU) return "gpu";
	if (backend == HK_FDG_BACKEND_AUTO) return "auto";
	return "cpu";
}

static const char *input_contact_source_from_path(const char *path)
{
	if (path == 0)
		return "unknown";
	if (strstr(path, "impute") || strstr(path, "positive_control") ||
		strstr(path, "native_hickit"))
		return "imputed_or_positive_control_pairs";
	return "raw_pairs";
}

static int approved_pairs_realpath_match(const char *pairs_path)
{
	char approved[PATH_MAX];
	char got[PATH_MAX];
	if (pairs_path == 0 || pairs_path[0] == 0)
		return 0;
	if (realpath(HK_P9016_DEFAULT_PAIRS, approved) == 0) {
		fprintf(stderr, "failed to resolve approved P9016 pairs path %s: %s\n",
				HK_P9016_DEFAULT_PAIRS, strerror(errno));
		return 0;
	}
	if (realpath(pairs_path, got) == 0) {
		fprintf(stderr, "failed to resolve P9016 pairs path %s: %s\n",
				pairs_path, strerror(errno));
		return 0;
	}
	return strcmp(approved, got) == 0;
}

static void path_join(char *dst, size_t dst_size, const char *a, const char *b)
{
	size_t n;
	assert(dst);
	assert(dst_size > 0);
	assert(a);
	assert(b);
	n = strlen(a);
	if (n > 0 && a[n - 1] == '/')
		snprintf(dst, dst_size, "%s%s", a, b);
	else
		snprintf(dst, dst_size, "%s/%s", a, b);
}

static int mkdir_if_missing(const char *path)
{
	if (mkdir(path, 0700) == 0)
		return 0;
	if (errno == EEXIST)
		return 0;
	fprintf(stderr, "failed to create directory %s: %s\n", path, strerror(errno));
	return -1;
}

static int make_output_root(char *root, size_t root_size)
{
	const char *requested = getenv("HK_BLIND_P9016_OUTPUT_ROOT");
	time_t t;
	struct tm tmv;
	char stamp[32];
	if (requested && requested[0]) {
		snprintf(root, root_size, "%s", requested);
		return mkdir_if_missing(root);
	}
	t = time(0);
	localtime_r(&t, &tmv);
	strftime(stamp, sizeof(stamp), "%Y%m%d_%H%M%S", &tmv);
	snprintf(root, root_size, "/tmp/hk_blind_p9016_minimal_%s_%ld", stamp, (long)getpid());
	return mkdir_if_missing(root);
}

static void stage_state_destroy(struct stage_state *state)
{
	if (state == 0) return;
	free(state->coords);
	free(state->raw);
	if (state->bmap) hk_bmap_destroy(state->bmap);
	if (state->map) hk_map_destroy(state->map);
	memset(state, 0, sizeof(*state));
}

static const char *resolution_label(int bin_size_bp, char label[32])
{
	if (bin_size_bp % 1000000 == 0)
		snprintf(label, 32, "%dMb", bin_size_bp / 1000000);
	else if (bin_size_bp % 1000 == 0)
		snprintf(label, 32, "%dkb", bin_size_bp / 1000);
	else
		snprintf(label, 32, "%dbp", bin_size_bp);
	return label;
}

static const char *resolution_slug(int bin_size_bp, char label[32])
{
	if (bin_size_bp % 1000000 == 0)
		snprintf(label, 32, "%dm", bin_size_bp / 1000000);
	else if (bin_size_bp % 1000 == 0)
		snprintf(label, 32, "%dk", bin_size_bp / 1000);
	else
		snprintf(label, 32, "%dbp", bin_size_bp);
	return label;
}

static int parse_resolution_chain(const char *s, int *bins, int max_bins)
{
	char *buf, *tok, *saveptr = 0;
	int n = 0;
	if (s == 0 || s[0] == 0)
		return 0;
	buf = strdup(s);
	if (buf == 0)
		return -1;
	for (tok = strtok_r(buf, ",", &saveptr); tok; tok = strtok_r(0, ",", &saveptr)) {
		char *end = 0;
		long v;
		if (n >= max_bins) {
			free(buf);
			return -1;
		}
		errno = 0;
		v = strtol(tok, &end, 10);
		if (errno || end == tok || *end != 0 || v <= 0 || v > INT32_MAX) {
			free(buf);
			return -1;
		}
		bins[n++] = (int)v;
	}
	free(buf);
	return n;
}

static int write_coarse_to_fine_map(const char *path, const struct hk_bmap *coarse,
									const struct hk_bmap *fine,
									const struct hk_blind_coarse_to_fine_map *map)
{
	FILE *fp = fopen(path, "w");
	int ret;
	if (fp == 0) return -1;
	ret = hk_blind_write_coarse_to_fine_map_tsv(fp, coarse, fine, map);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int write_bmap_summary(const char *path, const struct hk_bmap *bmap)
{
	FILE *fp = fopen(path, "w");
	int ret;
	if (fp == 0) return -1;
	ret = hk_blind_write_bmap_summary_tsv(fp, bmap, 0);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int read_blind_pairs(const char *pairs_path, int bin_size_bp,
							struct hk_map **map_out, struct hk_bmap **bmap_out,
							struct hk_blind_pair **raw_out, int32_t *n_raw_out)
{
	struct hk_map *map = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	int32_t i;

	assert(map_out);
	assert(bmap_out);
	assert(raw_out);
	assert(n_raw_out);
	map = hk_map_read(pairs_path);
	if (map == 0 || map->pairs == 0 || map->n_pairs <= 0) {
		fprintf(stderr, "failed to read P9016 pairs from %s\n", pairs_path);
		goto fail;
	}
	hk_pair_count_nei(map->n_pairs, map->pairs,
					  HK_P9016_BASE_K_NEIGHBOR_RADIUS,
					  HK_P9016_BASE_K_NEIGHBOR_RADIUS);
	bmap = hk_bmap_gen(map->d, map->n_pairs, map->pairs, bin_size_bp, 1);
	if (bmap == 0 || bmap->n_beads <= 0) {
		fprintf(stderr, "failed to build %d bp bin map\n", bin_size_bp);
		goto fail;
	}
	raw = (struct hk_blind_pair*)calloc((size_t)map->n_pairs, sizeof(*raw));
	if (raw == 0)
		goto fail;
	for (i = 0; i < map->n_pairs; ++i)
		hk_blind_pair_from_pair(&raw[i], &map->pairs[i]);
	*map_out = map;
	*bmap_out = bmap;
	*raw_out = raw;
	*n_raw_out = map->n_pairs;
	return 0;
fail:
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (map) hk_map_destroy(map);
	return -1;
}

static const char *init_scaffold_source_name(int init_mode)
{
	switch (init_mode) {
	case HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT: return "unphased_fdg";
	case HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT: return "random_haploid";
	case HK_BLIND_INIT_TOY_SPLIT: return "toy_haploid";
	case HK_BLIND_INIT_RANDOM_DIPLOID: return "none";
	default: return "unknown";
	}
}

static int init_scaffold_fdg_n_iter(int init_mode)
{
	return init_mode == HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT?
		HK_P9016_SCAFFOLD_FDG_N_ITER : 0;
}

static float init_eps_effective(int init_mode, const struct hk_p9016_minimal_params *p)
{
	return init_mode == HK_BLIND_INIT_RANDOM_DIPLOID? 0.0f : p->init_eps;
}

static float init_noise_scale_effective(int init_mode,
										const struct hk_p9016_minimal_params *p)
{
	return init_mode == HK_BLIND_INIT_RANDOM_DIPLOID? 0.0f : p->init_noise_scale;
}

static float init_scale_effective(int init_mode, float init_scale)
{
	return (init_mode == HK_BLIND_INIT_RANDOM_DIPLOID ||
			init_mode == HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT)? init_scale : 0.0f;
}

static int init_minimal_coords(struct hk_bmap *bmap, fvec3_t *haploid, fvec3_t *diploid,
							   int init_mode, float init_scale,
							   const struct hk_p9016_minimal_params *p)
{
	struct hk_fdg_conf scaffold_conf;
	assert(bmap);
	assert(haploid);
	assert(diploid);
	assert(p);
	switch (init_mode) {
	case HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT:
		hk_fdg_conf_init(&scaffold_conf);
		scaffold_conf.backend = HK_FDG_BACKEND_CPU;
		scaffold_conf.n_iter = HK_P9016_SCAFFOLD_FDG_N_ITER;
		if (hk_blind_init_haploid_scaffold_from_bmap_fdg(bmap, &scaffold_conf, haploid,
														 p->init_seed) != 0)
			return -1;
		return hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
														 diploid, p->init_eps,
														 p->init_noise_scale,
														 p->init_seed);
	case HK_BLIND_INIT_RANDOM_DIPLOID:
		return hk_blind_init_random_diploid_coords(bmap, diploid, init_scale, p->init_seed);
	case HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT:
		if (hk_blind_init_random_haploid_scaffold(bmap, haploid, init_scale, p->init_seed) != 0)
			return -1;
		return hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
														 diploid, p->init_eps,
														 p->init_noise_scale,
														 p->init_seed);
	case HK_BLIND_INIT_TOY_SPLIT:
		if (hk_blind_init_toy_haploid_scaffold(bmap, haploid) != 0)
			return -1;
		return hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
														 diploid, p->init_eps,
														 p->init_noise_scale,
														 p->init_seed);
	default:
		fprintf(stderr, "unsupported init mode %d\n", init_mode);
		return -1;
	}
}

static void set_minimal_schedule(struct hk_blind_iter_schedule_conf *conf, int n_iter,
								 int relax_steps, float relax_step,
								 const struct hk_p9016_minimal_params *p)
{
	memset(conf, 0, sizeof(*conf));
	conf->n_iter = n_iter;
	conf->base_conf.unit = HK_P9016_UNIT;
	conf->base_conf.d_scale = HK_P9016_D_SCALE;
	conf->base_conf.base_k = 1.0f;
	conf->base_conf.temperature = 1.0f;
	conf->base_conf.rho_train = 1.0f;
	conf->base_conf.rho_train_floor = 0.0f;
	conf->base_conf.min_sep_unit = p->min_sep_unit;
	conf->base_conf.lambda_sep = p->lambda_sep;
	conf->base_conf.chr_sep_unit = 0.0f;
	conf->base_conf.lambda_chr_sep = 0.0f;
	conf->base_conf.relax_step = relax_step;
	conf->base_conf.relax_steps = relax_steps;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.repulsion_block_k_min = 0.0f;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = p->d_scale_mode;
	conf->base_conf.d_scale_eps_count = p->d_scale_eps_count;
	conf->base_conf.estep_score_mode = HK_BLIND_ESTEP_SCORE_FDG_FLAT;
	conf->temperature_start = 1.0f;
	conf->temperature_end = 1.0f;
	conf->rho_train_start = 1.0f;
	conf->rho_train_end = 1.0f;
}

static int float_cmp(const void *a_, const void *b_)
{
	const float a = *(const float*)a_;
	const float b = *(const float*)b_;
	return (a > b) - (a < b);
}

static float percentile_sorted(const float *values, int32_t n, double q)
{
	double pos, lo_f, hi_f, t;
	int32_t lo, hi;
	if (n <= 0 || values == 0)
		return NAN;
	if (q <= 0.0)
		return values[0];
	if (q >= 1.0)
		return values[n - 1];
	pos = q * (double)(n - 1);
	lo_f = floor(pos);
	hi_f = ceil(pos);
	lo = (int32_t)lo_f;
	hi = (int32_t)hi_f;
	t = pos - lo_f;
	return (float)((1.0 - t) * (double)values[lo] + t * (double)values[hi]);
}

static int compute_sep_distribution(const struct hk_bmap *bmap, const fvec3_t *coords,
									float min_sep_unit, float lambda_sep,
									struct hk_p9016_sep_distribution *diag)
{
	float *values;
	double sum = 0.0;
	int32_t i, n = 0, below = 0;
	assert(bmap);
	assert(coords);
	assert(diag);
	memset(diag, 0, sizeof(*diag));
	diag->n_beads = bmap->n_beads;
	diag->min_sep_unit = min_sep_unit;
	diag->lambda_sep = lambda_sep;
	if (bmap->n_beads <= 0)
		return 0;
	values = (float*)calloc((size_t)bmap->n_beads, sizeof(*values));
	if (values == 0)
		return -1;
	for (i = 0; i < bmap->n_beads; ++i) {
		int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		float dx = coords[b0][0] - coords[b1][0];
		float dy = coords[b0][1] - coords[b1][1];
		float dz = coords[b0][2] - coords[b1][2];
		float d = sqrtf(dx * dx + dy * dy + dz * dz);
		if (!isfinite(d))
			continue;
		values[n++] = d;
		sum += d;
		if (min_sep_unit > 0.0f && d < min_sep_unit)
			++below;
	}
	if (n <= 0) {
		free(values);
		return -1;
	}
	qsort(values, (size_t)n, sizeof(*values), float_cmp);
	diag->n_beads = n;
	diag->n_sep_below_min = below;
	diag->sep_min = values[0];
	diag->sep_p01 = percentile_sorted(values, n, 0.01);
	diag->sep_p05 = percentile_sorted(values, n, 0.05);
	diag->sep_p10 = percentile_sorted(values, n, 0.10);
	diag->sep_p25 = percentile_sorted(values, n, 0.25);
	diag->sep_median = percentile_sorted(values, n, 0.50);
	diag->sep_mean = (float)(sum / n);
	diag->sep_p75 = percentile_sorted(values, n, 0.75);
	diag->sep_p90 = percentile_sorted(values, n, 0.90);
	diag->sep_p95 = percentile_sorted(values, n, 0.95);
	diag->sep_p99 = percentile_sorted(values, n, 0.99);
	diag->sep_max = values[n - 1];
	diag->frac_sep_below_min = n > 0? (float)below / (float)n : 0.0f;
	free(values);
	return 0;
}

static int write_sep_diag(const char *path, const struct hk_bmap *bmap,
						  const fvec3_t *coords, float min_sep_unit,
						  float lambda_sep,
						  struct hk_p9016_sep_distribution *diag_out)
{
	struct hk_p9016_sep_distribution diag;
	FILE *fp;
	if (compute_sep_distribution(bmap, coords, min_sep_unit, lambda_sep, &diag) != 0)
		return -1;
	if (diag_out)
		*diag_out = diag;
	fp = fopen(path, "w");
	if (fp == 0)
		return -1;
	if (fprintf(fp,
				"n_beads\tmin_sep_unit\tlambda_sep\tsep_min\tsep_p01\tsep_p05\t"
				"sep_p10\tsep_p25\tsep_median\tsep_mean\tsep_p75\tsep_p90\t"
				"sep_p95\tsep_p99\tsep_max\tn_sep_below_min\tfrac_sep_below_min\n"
				"%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				"%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\n",
				diag.n_beads, diag.min_sep_unit, diag.lambda_sep,
				diag.sep_min, diag.sep_p01, diag.sep_p05, diag.sep_p10,
				diag.sep_p25, diag.sep_median, diag.sep_mean, diag.sep_p75,
				diag.sep_p90, diag.sep_p95, diag.sep_p99, diag.sep_max,
				diag.n_sep_below_min, diag.frac_sep_below_min) < 0) {
		fclose(fp);
		return -1;
	}
	return fclose(fp) == 0? 0 : -1;
}

static int check_coords_finite(const fvec3_t *coords, int32_t n)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			if (!isfinite(coords[i][a])) {
				fprintf(stderr, "non-finite coordinate at %d,%d\n", i, a);
				return 1;
			}
	return 0;
}

static int validate_final_bpair_set(const struct hk_blind_bpair_set *set)
{
	struct hk_blind_iter_diag diag;
	int failed = 0;
	hk_blind_iter_diag_init(&diag);
	hk_blind_iter_diag_validate_bpair_set(set, &diag);
	failed |= check_i32("posterior nonfinite", diag.n_posterior_nonfinite, 0);
	failed |= check_i32("posterior bad sum", diag.n_posterior_bad_sum, 0);
	failed |= check_i32("posterior out of range", diag.n_posterior_out_of_range, 0);
	failed |= check_i32("uncertainty nonfinite", diag.n_uncertainty_nonfinite, 0);
	failed |= check_i32("uncertainty out of range", diag.n_uncertainty_out_of_range, 0);
	return failed? -1 : 0;
}

static int write_required_outputs(const char *posterior_path, const char *coords_path,
								  const char *coords_gz_path, const char *diag_path,
								  const struct hk_bmap *bmap,
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
	ret = hk_blind_write_diploid_coords_tsv_gz(coords_gz_path, bmap, coords);
	if (ret != 0) return -1;
	fp = fopen(diag_path, "w");
	if (fp == 0) return -1;
	ret = hk_blind_write_iter_loop_diag_tsv(fp, loop_diag);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int write_raw_output(const char *raw_path, const struct hk_blind_pair *raw, int32_t n_raw,
							const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set)
{
	FILE *fp = fopen(raw_path, "w");
	int ret;
	if (fp == 0) return -1;
	ret = hk_blind_write_raw_contact_posterior_tsv(fp, raw, n_raw, bmap, set);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int write_force_class_diag(const char *force_path, const struct hk_fdg_conf *fdg_conf,
								  const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
								  const fvec3_t *coords,
								  const struct hk_blind_relax_diag *relax_diag,
								  int d_scale_mode, float d_scale_eps_count,
								  struct hk_blind_contact_class_diag *diag_out)
{
	struct hk_blind_wedge_list wedges;
	struct hk_blind_contact_class_diag diag;
	FILE *fp;
	int ret;
	hk_blind_wedge_list_init(&wedges);
	ret = hk_blind_wedge_list_build_softall_mode(&wedges, bmap, set,
												 d_scale_mode, d_scale_eps_count);
	if (ret == 0)
		ret = hk_blind_wedge_list_aggregate_exact(&wedges);
	if (ret == 0)
		ret = hk_blind_contact_class_diag_accumulate(fdg_conf, bmap, &wedges, coords,
													 bmap->n_beads, HK_P9016_UNIT, &diag);
	hk_blind_wedge_list_destroy(&wedges);
	if (ret != 0)
		return ret;
	if (diag_out)
		*diag_out = diag;
	fp = fopen(force_path, "w");
	if (fp == 0)
		return -1;
	if (fprintf(fp,
				"class\tn_wedges\tsum_wedge_k\tcontact_energy\tcontact_force_l1\t"
				"backbone_force_l1\trepulsion_force_l1\thomolog_sep_force_l1\tn_nonfinite\n"
				"cis\t%lld\t%.17g\t%.9g\t%.9g\t0\t0\t0\t%d\n"
				"trans\t%lld\t%.17g\t%.9g\t%.9g\t0\t0\t0\t%d\n"
				"total\t%lld\t%.17g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\n",
				(long long)diag.n_wedges_cis, diag.sum_wedge_k_cis,
				diag.contact_energy_cis, diag.contact_force_l1_cis, diag.n_nonfinite,
				(long long)diag.n_wedges_trans, diag.sum_wedge_k_trans,
				diag.contact_energy_trans, diag.contact_force_l1_trans, diag.n_nonfinite,
				(long long)(diag.n_wedges_cis + diag.n_wedges_trans),
				diag.sum_wedge_k_cis + diag.sum_wedge_k_trans,
				diag.contact_energy_cis + diag.contact_energy_trans,
				diag.contact_force_l1_cis + diag.contact_force_l1_trans,
				relax_diag? relax_diag->final_backbone_force_l1 : 0.0f,
				relax_diag? relax_diag->final_repulsion_force_l1 : 0.0f,
				relax_diag? relax_diag->final_sep_force_l1 : 0.0f,
				diag.n_nonfinite) < 0) {
		fclose(fp);
		return -1;
	}
	return fclose(fp) == 0? 0 : -1;
}

static int write_manifest(const char *manifest_path, const char *pairs_path, const char *out_dir,
						  const char *posterior_path, const char *coords_path,
						  const char *coords_gz_path, const char *diag_path,
						  const char *force_diag_path, const char *sep_diag_path,
						  const char *raw_path,
						  const char *map_path,
						  int bin_size_bp, int n_iter, int relax_steps, int write_raw,
						  float relax_step, const char *input_contact_source,
						  const char *stage_name, int stage_index, int n_stages,
						  int parent_bin_size_bp, const char *init_source,
						  float child_offset_step, float parent_anchor_k,
						  int init_mode, float init_scale,
						  const struct hk_p9016_minimal_params *p,
						  const char *relax_backend, int approved_pairs_realpath,
						  const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
						  const struct hk_blind_base_k_stats *base_k_stats,
						  const struct hk_blind_iter_loop_diag *loop_diag,
						  const struct hk_blind_contact_class_diag *final_graph_diag,
						  const struct hk_blind_relax_diag *final_relax_diag,
						  float k_rel_rep_effective)
{
	FILE *fp = fopen(manifest_path, "w");
	const char *status;
	double final_refreshed_sum_wedge_k;
	int64_t final_refreshed_n_wedges;
	char label[32];
	if (fp == 0) return -1;
	assert(p);
	assert(final_graph_diag);
	final_refreshed_sum_wedge_k = final_graph_diag->sum_wedge_k_cis +
		final_graph_diag->sum_wedge_k_trans;
	final_refreshed_n_wedges = final_graph_diag->n_wedges_cis +
		final_graph_diag->n_wedges_trans;
	resolution_label(bin_size_bp, label);
	status = (loop_diag->n_bad_iter == 0 &&
			  loop_diag->n_relax_nonfinite_iter == 0 &&
			  loop_diag->n_coord_nonfinite == 0)? "OK" : "WARN";
	if (fprintf(fp,
				"key\tvalue\n"
				"sample\tP9016\n"
				"runner_family\tp9016_minimal\n"
				"runner_version\t2026-06-15\n"
				"default_profile\tp9016_softall_minimal_v2\n"
				"config_name\t%s\n"
				"input_path\t%s\n"
				"input_contact_source\t%s\n"
				"output_dir\t%s\n"
				"n_raw\t%d\n"
				"n_bpair\t%d\n"
					"n_beads\t%d\n"
					"resolution\t%d\n"
					"n_iter\t%d\n"
					"bin_size_bp\t%d\n"
					"resolution_label\t%s\n"
					"refinement_stage\t%s\n"
					"refinement_stage_index\t%d\n"
					"refinement_n_stages\t%d\n"
					"parent_bin_size_bp\t%d\n"
					"refinement_init_source\t%s\n"
					"refinement_child_offset_step\t%.9g\n"
					"refinement_parent_anchor_k\t%.9g\n"
					"relax_backend\t%s\n"
					"approved_p9016_raw_pairs_realpath\t%d\n"
					"unit\t%.9g\n"
					"d_scale\t%.9g\n"
					"base_k_mode\tneighbor_median\n"
				"base_k_effective\t%.9g\n"
				"base_k_min\t%.9g\n"
				"base_k_mean\t%.9g\n"
				"base_k_max\t%.9g\n"
				"base_k_n_nonfinite\t%d\n"
				"init_mode\t%s\n"
				"init_eps\t%.9g\n"
				"init_noise_scale\t%.9g\n"
				"init_scale\t%.9g\n"
				"init_eps_effective\t%.9g\n"
				"init_noise_scale_effective\t%.9g\n"
				"init_scale_effective\t%.9g\n"
				"init_seed\t%llu\n"
				"scaffold_source\t%s\n"
				"scaffold_fdg_n_iter\t%d\n"
				"prior_mode\tuniform\n"
				"prior_eps\t%.9g\n"
				"prior_inter_density\t0\n"
				"prior_observed_inter\t0\n"
				"prior_possible_inter\t0\n"
				"prior_n_distance_bins\t0\n"
				"prior_n_alpha\t0\n"
				"prior_alpha_min\t0\n"
				"prior_alpha_median\t0\n"
				"prior_alpha_max\t0\n"
				"prior_smoothing_method\tnone\n"
				"prior_alpha_clamp_min\t%.9g\n"
				"prior_alpha_clamp_max\t0.5\n"
				"rho_train_mode\tconstant\n"
				"rho_train\t1\n"
				"rho_train_floor\t0\n"
				"baseline\tsoftall\n"
				"mstep_graph_mode\traw_expected_soft_all\n"
				"training_graph_mode\traw_expected_soft_all\n"
				"softall_mode\traw_expected_soft_all\n"
				"training_graph_weighted_filter\t1\n"
				"training_graph_probability_weighted\t1\n"
				"training_graph_dscale_probability_weighted\t1\n"
				"estep_score_mode\tfdg_flat\n"
				"d_scale_mode\t%s\n"
				"d_scale_eps_count\t%.9g\n"
				"same_bin_filter_enabled\t%d\n"
				"n_raw_same_bin_excluded\t%lld\n"
				"n_bpair_same_bin_excluded\t%d\n"
				"n_raw_cis\t%lld\n"
				"n_raw_trans\t%lld\n"
				"n_bpair_cis\t%d\n"
				"n_bpair_trans\t%d\n"
				"raw_posterior_same_bin_policy\tuniform_unknown_rows\n"
				"copy_labels_are_gauge_only\t1\n"
				"uses_phase_labels\t0\n"
				"uses_charm_or_reference\t0\n"
				"uses_charm_for_training\t0\n"
				"min_sep_unit\t%.9g\n"
				"lambda_sep\t%.9g\n"
				"chr_sep_unit\t%.9g\n"
				"lambda_chr_sep\t%.9g\n"
				"relax_step\t%.9g\n"
				"relax_steps\t%d\n"
				"temperature_start\t1\n"
				"temperature_end\t1\n"
				"rho_train_start\t1\n"
				"rho_train_end\t1\n"
				"enable_repulsion\t1\n"
				"repulsion_mode\t%d\n"
				"repulsion_blocking_mode\tcurrent_edge_blocking\n"
				"repulsion_multiplier\t%.9g\n"
				"k_rel_rep_effective\t%.9g\n"
				"write_raw_posterior\t%d\n"
				"output_bpair_posterior\t%s\n"
				"output_coords\t%s\n"
				"output_coords_gz\t%s\n"
				"output_loop_diag\t%s\n"
				"output_force_class_diag\t%s\n"
				"output_sep_diag\t%s\n",
					p->config_name, pairs_path, input_contact_source, out_dir,
					set->n_raw, set->n_bpairs, bmap->n_beads,
					bin_size_bp, n_iter, bin_size_bp, label,
						stage_name, stage_index, n_stages, parent_bin_size_bp,
						init_source, child_offset_step, parent_anchor_k,
						relax_backend, approved_pairs_realpath,
					HK_P9016_UNIT, HK_P9016_D_SCALE,
				base_k_stats->mean, base_k_stats->min, base_k_stats->mean,
				base_k_stats->max, base_k_stats->n_nonfinite,
				hk_blind_init_mode_name(init_mode),
				p->init_eps, p->init_noise_scale, p->init_scale,
				init_eps_effective(init_mode, p), init_noise_scale_effective(init_mode, p),
				init_scale_effective(init_mode, init_scale),
				(unsigned long long)p->init_seed, init_scaffold_source_name(init_mode),
				init_scaffold_fdg_n_iter(init_mode),
				HK_P9016_PRIOR_EPS, HK_P9016_PRIOR_ALPHA_CLAMP_MIN,
				hk_blind_d_scale_mode_name(p->d_scale_mode),
				p->d_scale_eps_count, set->same_bin_filter_enabled,
				(long long)set->n_raw_same_bin_excluded,
				set->n_bpair_same_bin_excluded, (long long)set->n_raw_cis,
				(long long)set->n_raw_trans, set->n_bpair_cis, set->n_bpair_trans,
				p->min_sep_unit, p->lambda_sep, 0.0f, 0.0f,
				relax_step, relax_steps, HK_BLIND_REPULSION_CELL,
				HK_P9016_REPULSION_MULTIPLIER, k_rel_rep_effective, write_raw,
				posterior_path, coords_path, coords_gz_path, diag_path,
				force_diag_path, sep_diag_path) < 0) {
			fclose(fp);
			return -1;
		}
		if (map_path && map_path[0] && fprintf(fp, "output_coarse_to_fine_map\t%s\n", map_path) < 0) {
			fclose(fp);
			return -1;
		}
		if (write_raw && fprintf(fp, "output_raw_posterior\t%s\n", raw_path) < 0) {
			fclose(fp);
			return -1;
	}
	if (fprintf(fp,
				"status\t%s\n"
				"final_mean_entropy\t%.9g\n"
				"final_mean_pU\t%.9g\n"
				"final_mean_sep\t%.9g\n"
				"final_min_sep\t%.9g\n"
				"final_max_sep\t%.9g\n"
				"final_sum_wedge_k\t%.17g\n"
				"last_training_sum_wedge_k\t%.17g\n"
				"final_refreshed_sum_wedge_k\t%.17g\n"
				"final_refreshed_n_wedges\t%lld\n"
				"final_mean_rho_train_bpair\t%.9g\n"
				"final_min_rho_train_bpair\t%.9g\n"
				"final_max_rho_train_bpair\t%.9g\n"
				"final_contact_energy\t%.9g\n"
				"final_repulsion_energy\t%.9g\n"
				"final_backbone_energy\t%.9g\n"
				"final_sep_energy\t%.9g\n"
				"final_sep_force_l1\t%.9g\n"
				"posterior_refreshed_after_final_relax\t%d\n"
				"posterior_refresh_temperature\t%.9g\n"
				"posterior_refresh_prior_mode\tuniform\n"
				"posterior_refresh_mean_kl\t%.17g\n"
				"posterior_refresh_top_state_switch_frac\t%.9g\n"
				"posterior_refresh_mean_pU_before\t%.9g\n"
				"posterior_refresh_mean_pU_after\t%.9g\n"
				"n_bad_iter\t%d\n"
				"n_relax_nonfinite_iter\t%d\n"
				"n_coord_nonfinite\t%d\n",
				status, loop_diag->final_mean_entropy, loop_diag->final_mean_pU,
				loop_diag->final_mean_sep, loop_diag->final_min_sep,
				loop_diag->final_max_sep, loop_diag->final_sum_wedge_k,
				loop_diag->final_sum_wedge_k, final_refreshed_sum_wedge_k,
				(long long)final_refreshed_n_wedges,
				loop_diag->final_mean_rho_train_bpair,
				loop_diag->final_min_rho_train_bpair,
				loop_diag->final_max_rho_train_bpair,
				final_relax_diag? final_relax_diag->final_contact_energy : 0.0f,
				final_relax_diag? final_relax_diag->final_repulsion_energy : loop_diag->final_repulsion_energy,
				final_relax_diag? final_relax_diag->final_backbone_energy : 0.0f,
				final_relax_diag? final_relax_diag->final_sep_energy : 0.0f,
				final_relax_diag? final_relax_diag->final_sep_force_l1 : 0.0f,
				loop_diag->posterior_refreshed_after_final_relax,
				loop_diag->posterior_refresh_temperature,
				loop_diag->posterior_refresh_mean_kl,
				loop_diag->posterior_refresh_top_state_switch_frac,
				loop_diag->posterior_refresh_mean_pU_before,
				loop_diag->posterior_refresh_mean_pU_after,
				loop_diag->n_bad_iter, loop_diag->n_relax_nonfinite_iter,
				loop_diag->n_coord_nonfinite) < 0) {
		fclose(fp);
		return -1;
	}
	return fclose(fp) == 0? 0 : -1;
}

static int write_summary_header(FILE *fp)
{
	return fprintf(fp,
				   "stage_name\tstage_index\tn_stages\tbin_size_bp\tparent_bin_size_bp\t"
				   "config_name\toutput_dir\tinput_contact_source\tbaseline\t"
				   "init_mode\tn_iter\trelax_steps\trelax_step\tn_raw\tn_bpair\tn_beads\t"
				   "n_raw_cis\tn_raw_trans\tn_bpair_cis\tn_bpair_trans\t"
				   "final_mean_entropy\tfinal_mean_pU\tfinal_mean_sep\tfinal_min_sep\t"
				   "final_max_sep\tfinal_sum_wedge_k\tfinal_refreshed_sum_wedge_k\t"
				   "final_refreshed_n_wedges\tfinal_mean_rho_train_bpair\t"
				   "final_min_rho_train_bpair\tfinal_max_rho_train_bpair\t"
				   "final_n_expanded_edges\tfinal_n_softall_candidate_pairs\t"
				   "final_n_softall_bmap_pairs\tfinal_n_softall_selected_raw\t"
				   "final_n_softall_same_bin_skip_raw\t"
				   "final_n_softall_state00_raw\tfinal_n_softall_state01_raw\t"
				   "final_n_softall_state10_raw\tfinal_n_softall_state11_raw\t"
				   "final_contact_energy\tfinal_repulsion_energy\tfinal_backbone_energy\t"
				   "final_n_repulsion_pairs_considered\tfinal_n_repulsion_pairs_blocked\t"
				   "final_n_repulsion_pairs_active\tn_same_bin_excluded\t"
				   "posterior_refreshed_after_final_relax\tn_bad_iter\t"
				   "n_relax_nonfinite_iter\tn_coord_nonfinite\taudit_status\t"
				   "write_raw_posterior\tinit_source\trelax_backend\n") < 0? -1 : 0;
}

static int append_summary_row(FILE *fp, const struct minimal_result *result, int n_iter,
							  int relax_steps, float relax_step,
							  const char *input_contact_source, int write_raw,
							  const char *config_name)
{
	return fprintf(fp,
				   "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\tsoftall\t%s\t"
				   "%d\t%d\t%.9g\t%d\t%d\t%d\t%lld\t%lld\t%d\t%d\t"
				   "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.17g\t%.17g\t%lld\t"
				   "%.9g\t%.9g\t%.9g\t%lld\t%lld\t%lld\t%lld\t%lld\t"
				   "%lld\t%lld\t%lld\t%lld\t%.9g\t%.9g\t%.9g\t%lld\t%lld\t%lld\t"
				   "%lld\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\n",
				   result->stage_name, result->stage_index, result->n_stages,
				   result->bin_size_bp, result->parent_bin_size_bp,
				   config_name, result->output_dir, input_contact_source,
				   result->init_mode, n_iter, relax_steps, relax_step,
				   result->n_raw, result->n_bpair, result->n_beads,
				   (long long)result->n_raw_cis, (long long)result->n_raw_trans,
				   result->n_bpair_cis, result->n_bpair_trans,
				   result->final_mean_entropy, result->final_mean_pU,
				   result->final_mean_sep, result->final_min_sep, result->final_max_sep,
				   result->final_sum_wedge_k, result->final_refreshed_sum_wedge_k,
				   (long long)result->final_refreshed_n_wedges,
				   result->final_mean_rho_train_bpair,
				   result->final_min_rho_train_bpair, result->final_max_rho_train_bpair,
				   (long long)result->final_n_expanded_edges,
				   (long long)result->final_n_softall_candidate_pairs,
				   (long long)result->final_n_softall_bmap_pairs,
				   (long long)result->final_n_softall_selected_raw,
				   (long long)result->final_n_softall_same_bin_skip_raw,
				   (long long)result->final_n_softall_state_raw_count[HK_BLIND_STATE_00],
				   (long long)result->final_n_softall_state_raw_count[HK_BLIND_STATE_01],
				   (long long)result->final_n_softall_state_raw_count[HK_BLIND_STATE_10],
				   (long long)result->final_n_softall_state_raw_count[HK_BLIND_STATE_11],
				   result->final_contact_energy, result->final_repulsion_energy,
				   result->final_backbone_energy,
				   (long long)result->final_n_repulsion_pairs_considered,
				   (long long)result->final_n_repulsion_pairs_blocked,
				   (long long)result->final_n_repulsion_pairs_active,
				   (long long)result->n_same_bin_excluded,
				   result->posterior_refreshed_after_final_relax,
				   result->n_bad_iter, result->n_relax_nonfinite_iter,
				   result->n_coord_nonfinite,
				   result->status_ok? "OK" : "WARN", write_raw,
				   result->init_source, result->relax_backend) < 0? -1 : 0;
}

static int run_minimal_config(struct hk_bmap *bmap, const struct hk_blind_pair *raw,
							  int32_t n_raw, const char *pairs_path, const char *root_dir,
							  int bin_size_bp, int n_iter, int relax_steps, int write_raw,
							  float relax_step, const char *input_contact_source,
							  const char *stage_name, int stage_index, int n_stages,
							  const struct hk_bmap *parent_bmap,
							  const fvec3_t *parent_coords,
							  int parent_bin_size_bp,
							  float child_offset_step, float parent_anchor_k,
							  int init_mode, float init_scale,
							  const struct hk_p9016_minimal_params *p,
							  int relax_backend, int approved_pairs_realpath,
							  FILE *summary_fp, struct minimal_result *result,
							  fvec3_t **coords_out)
{
	struct hk_blind_bpair_set *set = 0;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	struct hk_blind_base_k_stats base_k_stats;
	struct hk_blind_contact_class_diag final_graph_diag;
	struct hk_blind_single_iter_diag *per_iter = 0;
	struct hk_blind_coarse_to_fine_map *anchor_map = 0;
	fvec3_t *haploid = 0;
	fvec3_t *diploid = 0;
	int32_t n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	char posterior_path[1024], coords_path[1024], coords_gz_path[1024], diag_path[1024];
	char force_diag_path[1024], sep_diag_path[1024], raw_path[1024], manifest_path[1024], map_path[1024];
	char init_source[128];
	char backend_name[32];
	int failed = 0;

	assert(stage_name);
	assert(p);
	assert(coords_out);
	*coords_out = 0;
	memset(result, 0, sizeof(*result));
	hk_blind_contact_class_diag_init(&final_graph_diag);
	path_join(result->output_dir, sizeof(result->output_dir), root_dir, p->config_name);
	if (mkdir_if_missing(result->output_dir) != 0)
		return 1;
	path_join(posterior_path, sizeof(posterior_path), result->output_dir, "p9016_full.bpair_posterior.tsv");
	path_join(coords_path, sizeof(coords_path), result->output_dir, "p9016_full.coords.tsv");
	path_join(coords_gz_path, sizeof(coords_gz_path), result->output_dir, "p9016_full.coords.tsv.gz");
	path_join(diag_path, sizeof(diag_path), result->output_dir, "p9016_full.loop_diag.tsv");
	path_join(force_diag_path, sizeof(force_diag_path), result->output_dir, "p9016_full.force_class_diag.tsv");
	path_join(sep_diag_path, sizeof(sep_diag_path), result->output_dir, "p9016_full.sep_diag.tsv");
	path_join(raw_path, sizeof(raw_path), result->output_dir, "p9016_full.raw_posterior.tsv");
	path_join(manifest_path, sizeof(manifest_path), result->output_dir, "p9016_full.manifest.tsv");
	path_join(map_path, sizeof(map_path), result->output_dir, "p9016_full.coarse_to_fine_map.tsv");
	snprintf(backend_name, sizeof(backend_name), "%s", relax_backend_name(relax_backend));
	snprintf(result->stage_name, sizeof(result->stage_name), "%s", stage_name);
	result->stage_index = stage_index;
	result->n_stages = n_stages;
	result->bin_size_bp = bin_size_bp;
	result->parent_bin_size_bp = parent_bin_size_bp;
	snprintf(result->init_mode, sizeof(result->init_mode), "%s",
			 hk_blind_init_mode_name(init_mode));
	snprintf(result->relax_backend, sizeof(result->relax_backend), "%s", backend_name);

	fprintf(stderr, "minimal P9016[%s]: build binned contacts\n", stage_name);
	set = hk_blind_bpair_set_build(bmap, n_raw, raw);
	failed |= check_true("bpair set", set != 0);
	if (failed) goto cleanup;
	failed |= check_i32("base-k mode", hk_blind_bpair_set_apply_base_k_mode(bmap, set,
																			 HK_BLIND_BASE_K_NEIGHBOR_MEDIAN), 0);
	hk_blind_bpair_set_init_uniform_prior(set);
	hk_blind_bpair_set_base_k_stats(set, &base_k_stats);
	if (failed) goto cleanup;

	diploid = (fvec3_t*)calloc((size_t)n_diploid, sizeof(*diploid));
	per_iter = (struct hk_blind_single_iter_diag*)calloc((size_t)n_iter, sizeof(*per_iter));
	if (diploid == 0 || per_iter == 0) {
		failed = 1;
		goto cleanup;
	}
	if (parent_bmap && parent_coords) {
		anchor_map = hk_blind_coarse_to_fine_map_build(parent_bmap, bmap);
		failed |= check_true("coarse-to-fine map", anchor_map != 0);
		if (failed) goto cleanup;
		failed |= check_i32("write coarse-to-fine map",
							write_coarse_to_fine_map(map_path, parent_bmap, bmap, anchor_map), 0);
		fprintf(stderr, "minimal P9016[%s]: lift coords from parent %d bp using child_offset_step=%.8g anchor_k=%.8g\n",
				stage_name, parent_bin_size_bp, child_offset_step, parent_anchor_k);
		failed |= check_i32("lift coords",
							hk_blind_lift_from_4mb(parent_bmap, anchor_map, parent_coords,
												   diploid, child_offset_step), 0);
		snprintf(init_source, sizeof(init_source), "previous_resolution_lift");
	} else {
		haploid = (fvec3_t*)calloc((size_t)bmap->n_beads, sizeof(*haploid));
		if (haploid == 0) {
			failed = 1;
			goto cleanup;
		}
		fprintf(stderr, "minimal P9016[%s]: init mode=%s seed=%llu scale=%.8g\n",
				stage_name, hk_blind_init_mode_name(init_mode),
				(unsigned long long)p->init_seed, init_scale);
		failed |= check_i32("init coords",
							init_minimal_coords(bmap, haploid, diploid,
												init_mode, init_scale, p), 0);
		snprintf(init_source, sizeof(init_source), "%s", hk_blind_init_mode_name(init_mode));
	}
	snprintf(result->init_source, sizeof(result->init_source), "%s", init_source);
	failed |= check_coords_finite(diploid, n_diploid);
	if (failed) goto cleanup;

	fprintf(stderr, "minimal P9016[%s]: run softall, n_iter=%d relax_steps=%d relax_step=%.8g relax_backend=%s\n",
			stage_name, n_iter, relax_steps, relax_step, backend_name);
	hk_fdg_conf_init(&fdg_conf);
	fdg_conf.backend = (enum hk_fdg_backend)relax_backend;
	fdg_conf.k_rel_rep *= HK_P9016_REPULSION_MULTIPLIER;
	set_minimal_schedule(&schedule_conf, n_iter, relax_steps, relax_step, p);
	if (anchor_map)
		failed |= check_i32("scheduled parent-anchor run",
							hk_blind_run_iter_loop_scheduled_parent_anchor_cpu(bmap, set, &fdg_conf, diploid,
																			   0, &schedule_conf, anchor_map,
																			   parent_coords, parent_anchor_k,
																			   per_iter, &loop_diag), 0);
	else
		failed |= check_i32("scheduled run",
							hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &fdg_conf, diploid,
																 0, &schedule_conf, per_iter,
																 &loop_diag), 0);
	failed |= check_coords_finite(diploid, n_diploid);
	failed |= validate_final_bpair_set(set);
	failed |= check_i32("bad iter", loop_diag.n_bad_iter, 0);
	failed |= check_i32("relax nonfinite iter", loop_diag.n_relax_nonfinite_iter, 0);
	failed |= check_i32("coord nonfinite", loop_diag.n_coord_nonfinite, 0);
	if (failed) goto cleanup;

	fprintf(stderr, "minimal P9016[%s]: write outputs to %s\n", stage_name, result->output_dir);
	failed |= check_i32("write required outputs",
						write_required_outputs(posterior_path, coords_path, coords_gz_path, diag_path,
											   bmap, set, diploid, &loop_diag), 0);
	failed |= check_i32("write force class diag",
						write_force_class_diag(force_diag_path, &fdg_conf, bmap, set, diploid,
											   n_iter > 0? &per_iter[n_iter - 1].relax_diag : 0,
											   p->d_scale_mode, p->d_scale_eps_count,
											   &final_graph_diag), 0);
	failed |= check_i32("write sep diag",
						write_sep_diag(sep_diag_path, bmap, diploid, p->min_sep_unit,
									   p->lambda_sep, 0), 0);
	if (write_raw)
		failed |= check_i32("write raw posterior", write_raw_output(raw_path, raw, n_raw, bmap, set), 0);
	failed |= check_i32("write manifest",
						write_manifest(manifest_path, pairs_path, result->output_dir,
									   posterior_path, coords_path, coords_gz_path, diag_path,
									   force_diag_path, sep_diag_path, raw_path,
									   anchor_map? map_path : 0,
									   bin_size_bp, n_iter, relax_steps,
									   write_raw, relax_step, input_contact_source,
									   stage_name, stage_index, n_stages, parent_bin_size_bp,
									   init_source, child_offset_step, parent_anchor_k,
									   init_mode, init_scale, p, backend_name,
									   approved_pairs_realpath,
									   bmap, set,
									   &base_k_stats, &loop_diag, &final_graph_diag,
									   n_iter > 0? &per_iter[n_iter - 1].relax_diag : 0,
									   fdg_conf.k_rel_rep), 0);
	if (failed) goto cleanup;

	result->n_raw = set->n_raw;
	result->n_bpair = set->n_bpairs;
	result->n_beads = bmap->n_beads;
	result->n_raw_cis = set->n_raw_cis;
	result->n_raw_trans = set->n_raw_trans;
	result->n_bpair_cis = set->n_bpair_cis;
	result->n_bpair_trans = set->n_bpair_trans;
	result->n_same_bin_excluded = set->n_raw_same_bin_excluded;
	result->final_mean_entropy = loop_diag.final_mean_entropy;
	result->final_mean_pU = loop_diag.final_mean_pU;
	result->final_mean_sep = loop_diag.final_mean_sep;
	result->final_min_sep = loop_diag.final_min_sep;
	result->final_max_sep = loop_diag.final_max_sep;
	result->final_sum_wedge_k = loop_diag.final_sum_wedge_k;
	result->final_refreshed_sum_wedge_k = final_graph_diag.sum_wedge_k_cis +
		final_graph_diag.sum_wedge_k_trans;
	result->final_refreshed_n_wedges = final_graph_diag.n_wedges_cis +
		final_graph_diag.n_wedges_trans;
	result->final_mean_rho_train_bpair = loop_diag.final_mean_rho_train_bpair;
	result->final_min_rho_train_bpair = loop_diag.final_min_rho_train_bpair;
	result->final_max_rho_train_bpair = loop_diag.final_max_rho_train_bpair;
	result->final_n_expanded_edges = loop_diag.final_n_expanded_edges;
	result->final_n_softall_candidate_pairs = loop_diag.final_n_softall_candidate_pairs;
	result->final_n_softall_filter1_pairs = loop_diag.final_n_softall_filter1_pairs;
	result->final_n_softall_filter2_pairs = loop_diag.final_n_softall_filter2_pairs;
	result->final_n_softall_bmap_pairs = loop_diag.final_n_softall_bmap_pairs;
	result->final_n_softall_selected_raw = loop_diag.final_n_softall_selected_raw;
	result->final_n_softall_gate_skip_raw = loop_diag.final_n_softall_gate_skip_raw;
	result->final_n_softall_same_bin_skip_raw = loop_diag.final_n_softall_same_bin_skip_raw;
	memcpy(result->final_n_softall_state_raw_count, loop_diag.final_n_softall_state_raw_count,
		   sizeof(result->final_n_softall_state_raw_count));
	result->k_rel_rep_effective = fdg_conf.k_rel_rep;
	result->final_contact_energy = per_iter[n_iter - 1].relax_diag.final_contact_energy;
	result->final_repulsion_energy = per_iter[n_iter - 1].relax_diag.final_repulsion_energy;
	result->final_backbone_energy = per_iter[n_iter - 1].relax_diag.final_backbone_energy;
	result->final_n_repulsion_pairs_considered = per_iter[n_iter - 1].relax_diag.final_n_repulsion_pairs_considered;
	result->final_n_repulsion_pairs_blocked = per_iter[n_iter - 1].relax_diag.final_n_repulsion_pairs_blocked;
	result->final_n_repulsion_pairs_active = per_iter[n_iter - 1].relax_diag.final_n_repulsion_pairs_active;
	result->posterior_refreshed_after_final_relax = loop_diag.posterior_refreshed_after_final_relax;
	result->n_bad_iter = loop_diag.n_bad_iter;
	result->n_relax_nonfinite_iter = loop_diag.n_relax_nonfinite_iter;
	result->n_coord_nonfinite = loop_diag.n_coord_nonfinite;
	result->status_ok = loop_diag.n_bad_iter == 0 &&
		loop_diag.n_relax_nonfinite_iter == 0 &&
		loop_diag.n_coord_nonfinite == 0;
	failed |= check_i32("append summary", append_summary_row(summary_fp, result, n_iter,
															 relax_steps, relax_step,
															 input_contact_source, write_raw,
															 p->config_name), 0);
	fflush(summary_fp);
	*coords_out = diploid;
	diploid = 0;
	fprintf(stderr, "minimal P9016[%s]: status=%s final_entropy=%.8g final_pU=%.8g\n",
			stage_name, result->status_ok? "OK" : "WARN", result->final_mean_entropy,
			result->final_mean_pU);

cleanup:
	if (anchor_map) hk_blind_coarse_to_fine_map_destroy(anchor_map);
	free(per_iter);
	free(diploid);
	free(haploid);
	if (set) hk_blind_bpair_set_destroy(set);
	return failed? 1 : 0;
}

int main(void)
{
	const char *pairs_path = env_or_default("HK_BLIND_P9016_PAIRS", HK_P9016_DEFAULT_PAIRS);
	int bin_size_bp = env_int_or_default("HK_BLIND_P9016_BIN_SIZE_BP",
									 HK_P9016_DEFAULT_BIN_SIZE_BP, 1);
	int run_chain = env_flag_enabled("HK_BLIND_P9016_CHAIN");
	const char *chain_spec = env_or_default("HK_BLIND_P9016_RESOLUTION_CHAIN",
											HK_P9016_DEFAULT_CHAIN);
	int n_iter = env_int_or_default("HK_BLIND_P9016_MINIMAL_N_ITER",
									HK_P9016_DEFAULT_N_ITER, 1);
	int relax_steps = env_int_or_default("HK_BLIND_P9016_MINIMAL_RELAX_STEPS",
										 HK_P9016_DEFAULT_RELAX_STEPS, 0);
	float relax_step = env_float_or_default("HK_BLIND_P9016_RELAX_STEP",
											HK_P9016_DEFAULT_RELAX_STEP, 0.0f);
	float child_offset_step = env_float_or_default("HK_BLIND_P9016_CHILD_OFFSET_STEP",
												   HK_P9016_DEFAULT_CHILD_OFFSET_STEP, 0.0f);
	float parent_anchor_k = env_float_or_default("HK_BLIND_P9016_PARENT_ANCHOR_K",
												 HK_P9016_DEFAULT_PARENT_ANCHOR_K, 0.0f);
	int init_mode_ok = 1;
	int init_mode = parse_init_mode(&init_mode_ok);
	int relax_backend = parse_relax_backend();
	int approved_pairs_realpath;
	int write_raw = env_flag_enabled("HK_BLIND_WRITE_RAW");
	struct hk_p9016_minimal_params params;
	struct stage_state prev_state;
	int bins[HK_P9016_MAX_STAGES];
	int n_bins = 1;
	const char *input_contact_source;
	char root_dir[1024], summary_path[1024];
	FILE *summary_fp = 0;
	int stage_idx;
	int failed = 0;

	memset(&prev_state, 0, sizeof(prev_state));
	if (load_minimal_params(&params) != 0)
		return 1;
	if (!init_mode_ok)
		return 1;
	if ((init_mode == HK_BLIND_INIT_RANDOM_DIPLOID ||
		 init_mode == HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT) &&
		params.init_scale <= 0.0f) {
		fprintf(stderr, "HK_BLIND_P9016_INIT_SCALE must be > 0 for init_mode=%s\n",
				hk_blind_init_mode_name(init_mode));
		return 1;
	}
	if (make_output_root(root_dir, sizeof(root_dir)) != 0)
		return 1;
	input_contact_source = input_contact_source_from_path(pairs_path);
	if (strcmp(input_contact_source, "raw_pairs") != 0) {
		fprintf(stderr, "minimal P9016 requires raw P9016 pairs, got source=%s path=%s\n",
				input_contact_source, pairs_path);
		return 1;
	}
	approved_pairs_realpath = approved_pairs_realpath_match(pairs_path);
	if (!approved_pairs_realpath &&
		!env_flag_enabled("HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS")) {
		fprintf(stderr, "minimal P9016 baseline requires the approved raw pairs realpath: %s\n",
				HK_P9016_DEFAULT_PAIRS);
		fprintf(stderr, "got pairs path: %s\n", pairs_path);
		fprintf(stderr, "set HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS=1 only for explicit non-baseline diagnostics\n");
		return 1;
	}
	if (run_chain) {
		n_bins = parse_resolution_chain(chain_spec, bins, HK_P9016_MAX_STAGES);
		if (n_bins <= 0) {
			fprintf(stderr, "invalid HK_BLIND_P9016_RESOLUTION_CHAIN=%s\n", chain_spec);
			return 1;
		}
	} else {
		bins[0] = bin_size_bp;
		n_bins = 1;
	}
	if (relax_backend == HK_FDG_BACKEND_AUTO) {
		relax_backend = hk_fdg_gpu_is_available()? HK_FDG_BACKEND_GPU : HK_FDG_BACKEND_CPU;
	}
	fprintf(stderr,
			"minimal P9016: input=%s output_root=%s config=%s baseline=softall n_iter=%d relax_steps=%d relax_step=%.8g chain=%d init_mode=%s init_eps=%.8g init_noise=%.8g init_scale=%.8g init_seed=%llu min_sep=%.8g lambda_sep=%.8g d_scale_mode=%s d_scale_eps=%.8g relax_backend=%s\n",
			pairs_path, root_dir, params.config_name, n_iter, relax_steps, relax_step,
			run_chain, hk_blind_init_mode_name(init_mode),
			params.init_eps, params.init_noise_scale, params.init_scale,
			(unsigned long long)params.init_seed, params.min_sep_unit,
			params.lambda_sep, hk_blind_d_scale_mode_name(params.d_scale_mode),
			params.d_scale_eps_count, relax_backend_name(relax_backend));
	path_join(summary_path, sizeof(summary_path), root_dir, "matrix_summary.tsv");
	summary_fp = fopen(summary_path, "w");
	if (summary_fp == 0) {
		fprintf(stderr, "failed to open summary %s\n", summary_path);
		failed = 1;
		goto cleanup;
	}
	failed |= check_i32("write summary header", write_summary_header(summary_fp), 0);
	for (stage_idx = 0; !failed && stage_idx < n_bins; ++stage_idx) {
		struct stage_state cur_state;
		struct minimal_result result;
		char stage_name[32];
		char stage_root[1024];
		char bmap_summary_path[1024];
		fvec3_t *coords = 0;
		int parent_bin = stage_idx > 0? prev_state.bin_size_bp : 0;

		memset(&cur_state, 0, sizeof(cur_state));
		resolution_slug(bins[stage_idx], stage_name);
		if (run_chain)
			path_join(stage_root, sizeof(stage_root), root_dir, stage_name);
		else
			snprintf(stage_root, sizeof(stage_root), "%s", root_dir);
		if (mkdir_if_missing(stage_root) != 0) {
			failed = 1;
			break;
		}
		fprintf(stderr, "minimal P9016: stage %d/%d name=%s bin_size_bp=%d parent=%d\n",
				stage_idx + 1, n_bins, stage_name, bins[stage_idx], parent_bin);
		if (read_blind_pairs(pairs_path, bins[stage_idx], &cur_state.map,
							 &cur_state.bmap, &cur_state.raw, &cur_state.n_raw) != 0) {
			failed = 1;
			stage_state_destroy(&cur_state);
			break;
		}
		cur_state.bin_size_bp = bins[stage_idx];
		path_join(bmap_summary_path, sizeof(bmap_summary_path), stage_root, "bmap_summary.tsv");
		failed |= check_i32("write bmap summary", write_bmap_summary(bmap_summary_path, cur_state.bmap), 0);
		if (!failed)
			failed |= run_minimal_config(cur_state.bmap, cur_state.raw, cur_state.n_raw,
										 pairs_path, stage_root, bins[stage_idx],
										 n_iter, relax_steps, write_raw, relax_step,
											 input_contact_source, stage_name, stage_idx + 1, n_bins,
												 stage_idx > 0? prev_state.bmap : 0,
												 stage_idx > 0? (const fvec3_t*)prev_state.coords : 0,
												 parent_bin, child_offset_step, parent_anchor_k,
												 init_mode, params.init_scale, &params,
												 relax_backend, approved_pairs_realpath,
												 summary_fp, &result, &coords);
		if (failed) {
			free(coords);
			stage_state_destroy(&cur_state);
			break;
		}
		cur_state.coords = coords;
		stage_state_destroy(&prev_state);
		prev_state = cur_state;
		memset(&cur_state, 0, sizeof(cur_state));
	}
cleanup:
	if (summary_fp && fclose(summary_fp) != 0)
		failed = 1;
	stage_state_destroy(&prev_state);
	if (failed) {
		fprintf(stderr, "minimal P9016: FAILED\n");
		return 1;
	}
	fprintf(stderr, "minimal P9016: summary %s\n", summary_path);
	return 0;
}
