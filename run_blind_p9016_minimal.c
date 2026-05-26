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

#define HK_P9016_DEFAULT_PAIRS "../pairs/P9016.pairs.gz"
#define HK_P9016_DEFAULT_BIN_SIZE_BP 1000000
#define HK_P9016_DEFAULT_N_ITER 100
#define HK_P9016_DEFAULT_RELAX_STEPS 100
#define HK_P9016_UNIT 1.0f
#define HK_P9016_D_SCALE 1.0f
#define HK_P9016_D_SCALE_EPS_COUNT 1e-6f
#define HK_P9016_BASE_K_NEIGHBOR_RADIUS 10000000
#define HK_P9016_INIT_EPS 0.5f
#define HK_P9016_INIT_NOISE_SCALE 0.0f
#define HK_P9016_INIT_SEED 17ULL
#define HK_P9016_INIT_SCALE 1.0f
#define HK_P9016_SCAFFOLD_FDG_N_ITER 50
#define HK_P9016_DEFAULT_RELAX_STEP 0.012f
#define HK_P9016_REPULSION_MULTIPLIER 1.0f
#define HK_P9016_CONFIG_NAME "minimal_soft_sep_off"
#define HK_P9016_HARD_PC_SEED 17ULL
#define HK_P9016_HARD_PC_KIND_HARD_STATE 0
#define HK_P9016_HARD_PC_KIND_STATE_FREQ 1
#define HK_P9016_HARD_PC_KIND_IMPUTED_P4_TOP 2
#define HK_P9016_IMPUTED_P4_THRESHOLD 0.75f
#define HK_P9016_STRICT_RAW_PAIRS_REALPATH "/mnt/shared/zliu/CHARM/CHARM_mesc/data/pairs/P9016.pairs.gz"
#define HK_P9016_STRICT_SMOKE_PAIRS_REALPATH "/mnt/ssd/zliu/phase3/hickit/testdata/p9016_blind_smoke.pairs"

struct minimal_init_config {
	int mode;
	uint64_t seed;
	float scale;
	float eps;
	float noise_scale;
	int32_t scaffold_fdg_n_iter;
	char config_name[256];
};

struct minimal_run_config {
	struct minimal_init_config init;
	int hard_pc_sample_size_percent;
	int hard_pc_kind;
	int d_scale_mode;
	int state_weight_mode;
	int mstep_graph_mode;
	float trans_pmax_min;
	float trans_margin_min;
	float trans_posterior_power_gamma;
	float imputed_p4_threshold;
	int raw_split_confidence_mode;
	float raw_split_confidence_floor;
	float raw_split_trans_scale;
	float raw_split_trans_confidence_power;
	float raw_split_state_p_min;
	float raw_split_posterior_power_gamma;
	int raw_outlier_enable;
	float raw_outlier_beta_cis;
	float raw_outlier_beta_trans;
	float raw_outlier_prior_cis;
	float raw_outlier_prior_trans;
	float raw_soft_min_q;
	uint64_t raw_split_sample_salt;
	uint64_t hard_pc_seed;
	float min_sep_unit;
	float lambda_sep;
	float chr_sep_unit;
	float lambda_chr_sep;
	float contact_k_multiplier_cis;
	float contact_k_multiplier_trans;
	float trans_d_scale_multiplier;
	float temperature_start;
	float temperature_end;
	float repulsion_block_k_min;
	int estep_score_mode;
	int strict_no_prior_guard;
	char input_contact_source[64];
	char config_name[256];
};

struct minimal_result {
	char output_dir[1024];
	int32_t n_raw;
	int32_t n_bpair;
	int32_t n_beads;
	struct hk_blind_phase_lock_diag phase_lock_diag;
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
	double last_training_sum_wedge_k;
	double final_refreshed_sum_wedge_k;
	int64_t final_refreshed_n_wedges;
	float final_mean_rho_train_bpair;
	float final_min_rho_train_bpair;
	float final_max_rho_train_bpair;
	int64_t final_n_expanded_edges;
	int64_t final_n_split_candidate_pairs;
	int64_t final_n_split_filter1_pairs;
	int64_t final_n_split_filter2_pairs;
	int64_t final_n_split_bmap_pairs;
	int64_t final_n_split_selected_raw;
	int64_t final_n_split_gate_skip_raw;
	int64_t final_n_split_same_bin_skip_raw;
	int64_t final_n_split_locked_skip_raw;
	int64_t final_n_split_state_raw_count[HK_BLIND_N_STATE];
	struct hk_blind_raw_outlier_diag final_raw_outlier_diag;
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
		fprintf(stderr, "invalid %s=%s; using %llu\n", name, s,
				(unsigned long long)fallback);
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

static int env_flag_enabled(const char *name)
{
	const char *s = getenv(name);
	return s && s[0] && strcmp(s, "0") != 0 && strcmp(s, "false") != 0 && strcmp(s, "FALSE") != 0;
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

static int strict_no_prior_pairs_path_ok(const char *path)
{
	char resolved[4096];
	if (path == 0 || path[0] == 0)
		return 0;
	if (realpath(path, resolved) == 0)
		return 0;
	return strcmp(resolved, HK_P9016_STRICT_RAW_PAIRS_REALPATH) == 0 ||
		   strcmp(resolved, HK_P9016_STRICT_SMOKE_PAIRS_REALPATH) == 0;
}

static int strict_no_prior_source_ok(const char *source, const char *path)
{
	return source && strcmp(source, "raw_pairs") == 0 &&
		   strict_no_prior_pairs_path_ok(path);
}

static int hard_pc_kind_from_env(void)
{
	const char *s = getenv("HK_BLIND_P9016_HARD_PC_KIND");
	if (s == 0 || s[0] == 0 || strcmp(s, "hard_state") == 0 ||
		strcmp(s, "raw_phase_lock") == 0 || strcmp(s, "majority") == 0)
		return HK_P9016_HARD_PC_KIND_HARD_STATE;
	if (strcmp(s, "state_freq") == 0 || strcmp(s, "state_frequency") == 0 ||
		strcmp(s, "fixed_p4") == 0)
		return HK_P9016_HARD_PC_KIND_STATE_FREQ;
	if (strcmp(s, "imputed_p4_top") == 0 || strcmp(s, "pairs_p4_top") == 0 ||
		strcmp(s, "native_p4_top") == 0)
		return HK_P9016_HARD_PC_KIND_IMPUTED_P4_TOP;
	fprintf(stderr, "invalid HK_BLIND_P9016_HARD_PC_KIND=%s\n", s);
	return -1;
}

static const char *hard_pc_kind_name(int kind)
{
	switch (kind) {
	case HK_P9016_HARD_PC_KIND_HARD_STATE: return "raw_phase_lock";
	case HK_P9016_HARD_PC_KIND_STATE_FREQ: return "raw_phase_state_freq";
	case HK_P9016_HARD_PC_KIND_IMPUTED_P4_TOP: return "imputed_p4_top";
	default: return "unknown";
	}
}

static int d_scale_mode_from_env(void)
{
	const char *s = getenv("HK_BLIND_P9016_D_SCALE_MODE");
	if (s == 0 || s[0] == 0 || strcmp(s, "raw") == 0 || strcmp(s, "raw_count") == 0)
		return HK_BLIND_D_SCALE_RAW_COUNT;
	if (strcmp(s, "expected") == 0 || strcmp(s, "expected_count") == 0)
		return HK_BLIND_D_SCALE_EXPECTED_COUNT;
	if (strcmp(s, "density") == 0 || strcmp(s, "density_normalized_raw_count") == 0)
		return HK_BLIND_D_SCALE_DENSITY_NORMALIZED_RAW_COUNT;
	if (strcmp(s, "capped_density") == 0 || strcmp(s, "capped_density_raw_count") == 0)
		return HK_BLIND_D_SCALE_CAPPED_DENSITY_RAW_COUNT;
	fprintf(stderr, "invalid HK_BLIND_P9016_D_SCALE_MODE=%s\n", s);
	return -1;
}

static const char *d_scale_mode_config_suffix(int mode)
{
	switch (mode) {
	case HK_BLIND_D_SCALE_RAW_COUNT: return "";
	case HK_BLIND_D_SCALE_EXPECTED_COUNT: return "_dsexpected";
	case HK_BLIND_D_SCALE_DENSITY_NORMALIZED_RAW_COUNT: return "_dsdensity";
	case HK_BLIND_D_SCALE_CAPPED_DENSITY_RAW_COUNT: return "_dscappeddensity";
	default: return "_dsunknown";
	}
}

static int state_weight_mode_from_env(void)
{
	const char *s = getenv("HK_BLIND_P9016_STATE_WEIGHT_MODE");
	if (s == 0 || s[0] == 0 || strcmp(s, "posterior") == 0)
		return HK_BLIND_STATE_WEIGHT_POSTERIOR;
	if (strcmp(s, "binary_support") == 0 || strcmp(s, "support") == 0 ||
		strcmp(s, "selected") == 0)
		return HK_BLIND_STATE_WEIGHT_BINARY_SUPPORT;
	if (strcmp(s, "top_only") == 0 || strcmp(s, "hard_top") == 0 ||
		strcmp(s, "max_state") == 0)
		return HK_BLIND_STATE_WEIGHT_TOP_ONLY;
	fprintf(stderr, "invalid HK_BLIND_P9016_STATE_WEIGHT_MODE=%s\n", s);
	return -1;
}

static const char *state_weight_mode_config_suffix(int mode, float pmax_min, float margin_min)
{
	static char buf[80];
	switch (mode) {
	case HK_BLIND_STATE_WEIGHT_POSTERIOR:
		return "";
	case HK_BLIND_STATE_WEIGHT_BINARY_SUPPORT:
		snprintf(buf, sizeof(buf), "_swbinary_p%.3g", pmax_min);
		return buf;
	case HK_BLIND_STATE_WEIGHT_TOP_ONLY:
		snprintf(buf, sizeof(buf), "_swtop_p%.3g_m%.3g", pmax_min, margin_min);
		return buf;
	default:
		return "_swunknown";
	}
}

static int mstep_graph_mode_uses_final_phased_prob_weight(int mode);
static int mstep_graph_mode_uses_stochastic_selection(int mode);

static const char *posterior_gamma_config_suffix(int state_weight_mode, float gamma)
{
	static char buf[32];
	if (state_weight_mode != HK_BLIND_STATE_WEIGHT_POSTERIOR || gamma <= 1.0f)
		return "";
	snprintf(buf, sizeof(buf), "_g%.3g", gamma);
	return buf;
}

static const char *raw_split_posterior_gamma_config_suffix(int mstep_graph_mode, float gamma)
{
	static char buf[32];
	if (!mstep_graph_mode_uses_final_phased_prob_weight(mstep_graph_mode) || gamma <= 1.0f)
		return "";
	snprintf(buf, sizeof(buf), "_rpg%.3g", gamma);
	return buf;
}

static const char *raw_split_sample_salt_config_suffix(int mstep_graph_mode, uint64_t sample_salt)
{
	static char buf[64];
	if (!mstep_graph_mode_uses_stochastic_selection(mstep_graph_mode) || sample_salt == 0)
		return "";
	snprintf(buf, sizeof(buf), "_ss%llu", (unsigned long long)sample_salt);
	return buf;
}

static int mstep_graph_mode_uses_top_threshold_suffix(int mode)
{
	return mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_TOP ||
	   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_TOP1 ||
	   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_TOP2 ||
	   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP1 ||
	   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP2 ||
	   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1 ||
	   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1_CONF_KWEIGHT ||
	   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_SAMPLE1 ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF_WEIGHTED;
}

static int mstep_graph_mode_uses_stochastic_selection(int mode)
{
	return mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1 ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1_CONF_KWEIGHT ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_SAMPLE1 ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF_WEIGHTED;
}

static int mstep_graph_mode_uses_weighted_filter(int mode)
{
	return mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_FULL ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDWEIGHT ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDHALF ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_COUNT ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_LOCKED_ONLY ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_CIS ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_ALL ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_ORACLE_ALL ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_ALL ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER;
}

static int mstep_graph_mode_uses_final_phased_prob_weight(int mode)
{
	return mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP1 ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP2 ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1_CONF_KWEIGHT ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_SAMPLE1 ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF_WEIGHTED ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1 ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1_RENORM ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF ||
			   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER ||
				   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT ||
				   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDWEIGHT ||
				   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDHALF ||
				   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_COUNT ||
				   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_ALL ||
				   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER ||
			   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_LOCKED_ONLY ||
			   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT ||
			   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_CIS ||
			   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_ALL ||
			   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_ORACLE_ALL;
}

static int mstep_graph_mode_uses_probability_dscale(int mode)
{
	return mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP1 ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP2 ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_SAMPLE1 ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF_WEIGHTED ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL ||
			   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1 ||
			   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1_RENORM ||
				   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF ||
				   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER ||
				   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDWEIGHT ||
				   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDHALF ||
					   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_COUNT ||
					   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_ALL ||
					   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER ||
					   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_LOCKED_ONLY ||
					   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT ||
					   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_CIS ||
					   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_ALL ||
					   mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_ORACLE_ALL;
}

static int mstep_graph_mode_has_posterior_squared_risk(int mode)
{
	return mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1_CONF_KWEIGHT ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_SAMPLE1 ||
		   mode == HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF_WEIGHTED;
}

static const char *mstep_graph_threshold_config_suffix(int mode, float pmax_min,
													   float margin_min)
{
	static char buf[80];
	if (!mstep_graph_mode_uses_top_threshold_suffix(mode) ||
		(pmax_min <= 0.0f && margin_min <= 0.0f))
		return "";
	snprintf(buf, sizeof(buf), "_p%.3g_m%.3g", pmax_min, margin_min);
	return buf;
}

static int mstep_graph_mode_from_env(void)
{
	const char *s = getenv("HK_BLIND_P9016_MSTEP_GRAPH_MODE");
	if (s == 0 || s[0] == 0 || strcmp(s, "bpair") == 0 ||
		strcmp(s, "bpair_four_state") == 0 || strcmp(s, "four_state") == 0)
		return HK_BLIND_MSTEP_GRAPH_BPAIR;
	if (strcmp(s, "raw_split_top") == 0 || strcmp(s, "raw_split") == 0 ||
		strcmp(s, "split_then_bin") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_TOP;
	if (strcmp(s, "raw_expected_count") == 0 || strcmp(s, "raw_expected") == 0 ||
		strcmp(s, "expected_count_raw") == 0 || strcmp(s, "sufficient_stat") == 0 ||
		strcmp(s, "raw_sufficient_stat") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_COUNT;
	if (strcmp(s, "raw_expected_soft_all") == 0 ||
		strcmp(s, "raw_expected_all") == 0 ||
		strcmp(s, "soft_all") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_ALL;
	if (strcmp(s, "raw_expected_soft_outlier") == 0 ||
		strcmp(s, "soft_outlier") == 0 ||
		strcmp(s, "raw_soft_outlier") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER;
	if (strcmp(s, "raw_expected_locked_only") == 0 ||
		strcmp(s, "raw_expected_locked") == 0 ||
		strcmp(s, "locked_only") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_LOCKED_ONLY;
	if (strcmp(s, "raw_expected_pcut") == 0 ||
		strcmp(s, "raw_expected_gate") == 0 ||
		strcmp(s, "expected_pcut") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT;
	if (strcmp(s, "raw_expected_pcut_oracle_cis") == 0 ||
		strcmp(s, "raw_expected_pcut_cis_oracle") == 0 ||
		strcmp(s, "pcut_oracle_cis") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_CIS;
	if (strcmp(s, "raw_expected_pcut_oracle_all") == 0 ||
		strcmp(s, "raw_expected_pcut_all_oracle") == 0 ||
		strcmp(s, "pcut_oracle_all") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_ALL;
	if (strcmp(s, "raw_expected_oracle_all") == 0 ||
		strcmp(s, "raw_oracle_all") == 0 ||
		strcmp(s, "oracle_all") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_ORACLE_ALL;
	if (strcmp(s, "raw_split_soft") == 0 || strcmp(s, "raw_soft") == 0 ||
		strcmp(s, "soft_raw_split") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT;
	if (strcmp(s, "raw_split_soft_top1") == 0 || strcmp(s, "raw_soft_top1") == 0 ||
		strcmp(s, "soft_raw_split_top1") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_TOP1;
	if (strcmp(s, "raw_split_soft_top2") == 0 || strcmp(s, "raw_soft_top2") == 0 ||
		strcmp(s, "soft_raw_split_top2") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_TOP2;
	if (strcmp(s, "raw_split_soft_filtered_top1") == 0 ||
		strcmp(s, "raw_soft_filtered_top1") == 0 ||
		strcmp(s, "soft_filtered_top1") == 0 ||
		strcmp(s, "filtered_soft_top1") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_TOP1;
	if (strcmp(s, "raw_split_soft_filtered_top2") == 0 ||
		strcmp(s, "raw_soft_filtered_top2") == 0 ||
		strcmp(s, "soft_filtered_top2") == 0 ||
		strcmp(s, "filtered_soft_top2") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_TOP2;
	if (strcmp(s, "raw_split_soft_filtered_weighted_top1") == 0 ||
		strcmp(s, "raw_soft_filtered_weighted_top1") == 0 ||
		strcmp(s, "soft_filtered_weighted_top1") == 0 ||
		strcmp(s, "filtered_weighted_top1") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP1;
	if (strcmp(s, "raw_split_soft_filtered_weighted_top2") == 0 ||
		strcmp(s, "raw_soft_filtered_weighted_top2") == 0 ||
		strcmp(s, "soft_filtered_weighted_top2") == 0 ||
		strcmp(s, "filtered_weighted_top2") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP2;
	if (strcmp(s, "raw_split_soft_filtered_sample1") == 0 ||
		strcmp(s, "raw_soft_filtered_sample1") == 0 ||
		strcmp(s, "soft_filtered_sample1") == 0 ||
		strcmp(s, "filtered_sample1") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1;
	if (strcmp(s, "raw_split_soft_filtered_weighted_sample1") == 0 ||
		strcmp(s, "raw_soft_filtered_weighted_sample1") == 0 ||
		strcmp(s, "soft_filtered_weighted_sample1") == 0 ||
		strcmp(s, "filtered_weighted_sample1") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_SAMPLE1;
	if (strcmp(s, "raw_split_soft_filtered_sample1_conf_kweight") == 0 ||
		strcmp(s, "raw_soft_filtered_sample1_conf_kweight") == 0 ||
		strcmp(s, "soft_filtered_sample1_conf_kweight") == 0 ||
		strcmp(s, "filtered_sample1_conf_kweight") == 0 ||
		strcmp(s, "sample1_conf_kweight") == 0 ||
		strcmp(s, "sample1_kweight") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1_CONF_KWEIGHT;
	if (strcmp(s, "raw_split_soft_filtered_bernoulli") == 0 ||
		strcmp(s, "raw_soft_filtered_bernoulli") == 0 ||
		strcmp(s, "soft_filtered_bernoulli") == 0 ||
		strcmp(s, "filtered_bernoulli") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI;
	if (strcmp(s, "raw_split_soft_filtered_bernoulli_conf") == 0 ||
		strcmp(s, "raw_soft_filtered_bernoulli_conf") == 0 ||
		strcmp(s, "soft_filtered_bernoulli_conf") == 0 ||
		strcmp(s, "filtered_bernoulli_conf") == 0 ||
		strcmp(s, "bernoulli_conf") == 0 ||
		strcmp(s, "posterior_thin_conf") == 0 ||
		strcmp(s, "filtered_thin_conf") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF;
	if (strcmp(s, "raw_split_soft_filtered_bernoulli_conf_weighted") == 0 ||
		strcmp(s, "raw_soft_filtered_bernoulli_conf_weighted") == 0 ||
		strcmp(s, "soft_filtered_bernoulli_conf_weighted") == 0 ||
		strcmp(s, "filtered_bernoulli_conf_weighted") == 0 ||
		strcmp(s, "bernoulli_conf_weighted") == 0 ||
		strcmp(s, "posterior_thin_conf_weighted") == 0 ||
		strcmp(s, "filtered_thin_conf_weighted") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF_WEIGHTED;
	if (strcmp(s, "raw_split_soft_filtered_all") == 0 ||
		strcmp(s, "raw_soft_filtered_all") == 0 ||
		strcmp(s, "soft_filtered_all") == 0 ||
		strcmp(s, "filtered_all") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL;
	if (strcmp(s, "raw_split_soft_filtered_all_conf") == 0 ||
		strcmp(s, "raw_soft_filtered_all_conf") == 0 ||
		strcmp(s, "soft_filtered_all_conf") == 0 ||
		strcmp(s, "filtered_all_conf") == 0 ||
		strcmp(s, "filtered_entropy_all") == 0 ||
		strcmp(s, "filtered_entropy") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF;
	if (strcmp(s, "raw_split_soft_filtered_all_conf_wfilter") == 0 ||
		strcmp(s, "raw_soft_filtered_all_conf_wfilter") == 0 ||
		strcmp(s, "soft_filtered_all_conf_wfilter") == 0 ||
		strcmp(s, "filtered_all_conf_wfilter") == 0 ||
		strcmp(s, "weighted_filter_all_conf") == 0 ||
		strcmp(s, "wfilter_all_conf") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER;
	if (strcmp(s, "raw_split_soft_filtered_all_conf_wfilter_full") == 0 ||
		strcmp(s, "raw_soft_filtered_all_conf_wfilter_full") == 0 ||
		strcmp(s, "soft_filtered_all_conf_wfilter_full") == 0 ||
		strcmp(s, "filtered_all_conf_wfilter_full") == 0 ||
		strcmp(s, "weighted_filter_all_conf_full") == 0 ||
		strcmp(s, "wfilter_all_conf_full") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_FULL;
	if (strcmp(s, "raw_split_soft_filtered_all_conf_wfilter_kweight") == 0 ||
		strcmp(s, "raw_soft_filtered_all_conf_wfilter_kweight") == 0 ||
		strcmp(s, "soft_filtered_all_conf_wfilter_kweight") == 0 ||
		strcmp(s, "filtered_all_conf_wfilter_kweight") == 0 ||
		strcmp(s, "weighted_filter_all_conf_kweight") == 0 ||
		strcmp(s, "wfilter_all_conf_kweight") == 0 ||
		strcmp(s, "wfilter_kweight") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT;
	if (strcmp(s, "raw_split_soft_filtered_all_conf_wfilter_kdweight") == 0 ||
		strcmp(s, "raw_soft_filtered_all_conf_wfilter_kdweight") == 0 ||
		strcmp(s, "soft_filtered_all_conf_wfilter_kdweight") == 0 ||
		strcmp(s, "filtered_all_conf_wfilter_kdweight") == 0 ||
		strcmp(s, "weighted_filter_all_conf_kdweight") == 0 ||
		strcmp(s, "wfilter_all_conf_kdweight") == 0 ||
		strcmp(s, "wfilter_kdweight") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDWEIGHT;
	if (strcmp(s, "raw_split_soft_filtered_all_conf_wfilter_kdhalf") == 0 ||
		strcmp(s, "raw_soft_filtered_all_conf_wfilter_kdhalf") == 0 ||
		strcmp(s, "soft_filtered_all_conf_wfilter_kdhalf") == 0 ||
		strcmp(s, "filtered_all_conf_wfilter_kdhalf") == 0 ||
		strcmp(s, "weighted_filter_all_conf_kdhalf") == 0 ||
		strcmp(s, "wfilter_all_conf_kdhalf") == 0 ||
		strcmp(s, "wfilter_kdhalf") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDHALF;
	if (strcmp(s, "raw_split_soft_filtered_pcut_top1") == 0 ||
		strcmp(s, "raw_soft_filtered_pcut_top1") == 0 ||
		strcmp(s, "soft_filtered_pcut_top1") == 0 ||
		strcmp(s, "filtered_pcut_top1") == 0 ||
		strcmp(s, "pcut_top1") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1;
	if (strcmp(s, "raw_split_soft_filtered_pcut_top1_renorm") == 0 ||
		strcmp(s, "raw_soft_filtered_pcut_top1_renorm") == 0 ||
		strcmp(s, "soft_filtered_pcut_top1_renorm") == 0 ||
		strcmp(s, "filtered_pcut_top1_renorm") == 0 ||
		strcmp(s, "pcut_top1_renorm") == 0 ||
		strcmp(s, "pcut_renorm") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1_RENORM;
	fprintf(stderr, "invalid HK_BLIND_P9016_MSTEP_GRAPH_MODE=%s\n", s);
	return -1;
}

static const char *mstep_graph_mode_config_suffix(int mode)
{
	switch (mode) {
	case HK_BLIND_MSTEP_GRAPH_BPAIR:
		return "";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_TOP:
		return "_mrawsplit";
	case HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_COUNT:
		return "_mrawexpected";
	case HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_LOCKED_ONLY:
		return "_mrawexpectedlock";
	case HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT:
		return "_mrawexpectedpcut";
	case HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_CIS:
		return "_mrawexpectedpcutocis";
	case HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_ALL:
		return "_mrawexpectedpcutoall";
	case HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_ORACLE_ALL:
		return "_mrawexpectedoall";
	case HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_ALL:
		return "_mrawexpectedsoftall";
	case HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER:
		return "_mrawexpectedoutlier";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT:
		return "_mrawsoft";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_TOP1:
		return "_mrawsofttop1";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_TOP2:
		return "_mrawsofttop2";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_TOP1:
		return "_mrawsoftfilttop1";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_TOP2:
		return "_mrawsoftfilttop2";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP1:
		return "_mrawsoftfiltwtop1";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP2:
		return "_mrawsoftfiltwtop2";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1:
		return "_mrawsoftfiltsample1";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1_CONF_KWEIGHT:
		return "_mrawsoftfiltsample1confk";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_SAMPLE1:
		return "_mrawsoftfiltwsample1";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI:
		return "_mrawsoftfiltbernoulli";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF:
		return "_mrawsoftfiltbernoulliconf";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF_WEIGHTED:
		return "_mrawsoftfiltbernoulliconfw";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL:
		return "_mrawsoftfiltall";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF:
		return "_mrawsoftfiltallconf";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER:
		return "_mrawsoftfiltallconfwf";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_FULL:
		return "_mrawsoftfiltallconfwffull";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT:
		return "_mrawsoftfiltallconfwfk";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDWEIGHT:
		return "_mrawsoftfiltallconfwfkd";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDHALF:
		return "_mrawsoftfiltallconfwfkdhalf";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1:
		return "_mrawsoftfiltpcuttop1";
	case HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1_RENORM:
		return "_mrawsoftfiltpcuttop1renorm";
	default:
		return "_munknown";
	}
}

static const char *raw_split_state_p_min_config_suffix(int mstep_graph_mode, float p_min)
{
	static char buf[32];
	if ((mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1 &&
		 mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1_RENORM &&
		 mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT &&
		 mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_CIS &&
		 mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_ALL &&
		 mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_ORACLE_ALL) ||
		p_min <= 0.0f)
		return "";
	snprintf(buf, sizeof(buf), "_sp%.3g", p_min);
	return buf;
}

static const char *raw_outlier_config_suffix(int mstep_graph_mode, int enable,
											 float beta_cis, float beta_trans,
											 float prior_cis, float prior_trans,
											 float min_q)
{
	static char buf[128];
	if (mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER)
		return "";
	if (!enable)
		return "_uoff";
	if (beta_cis == HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS &&
		beta_trans == HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS &&
		prior_cis == HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS &&
		prior_trans == HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS &&
		min_q == HK_BLIND_RAW_SOFT_DEFAULT_MIN_Q)
		return "_u";
	snprintf(buf, sizeof(buf), "_u_bc%.3g_bt%.3g_pc%.3g_pt%.3g_mq%.3g",
			 beta_cis, beta_trans, prior_cis, prior_trans, min_q);
	return buf;
}

static int raw_split_confidence_mode_from_env(void)
{
	const char *s = getenv("HK_BLIND_P9016_RAW_SPLIT_CONFIDENCE_MODE");
	if (s == 0 || s[0] == 0 || strcmp(s, "entropy_linear") == 0 ||
		strcmp(s, "entropy") == 0)
		return HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR;
	if (strcmp(s, "pmax_margin") == 0 || strcmp(s, "pmax") == 0)
		return HK_BLIND_RAW_SPLIT_CONF_PMAX_MARGIN;
	if (strcmp(s, "floor_entropy") == 0 || strcmp(s, "entropy_floor") == 0)
		return HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY;
	fprintf(stderr, "invalid HK_BLIND_P9016_RAW_SPLIT_CONFIDENCE_MODE=%s\n", s);
	return -1;
}

static int estep_score_mode_from_env(void)
{
	const char *s = getenv("HK_BLIND_P9016_ESTEP_SCORE_MODE");
	if (s == 0 || s[0] == 0 || strcmp(s, "fdg_flat") == 0 ||
		strcmp(s, "fdg") == 0 || strcmp(s, "flat") == 0)
		return HK_BLIND_ESTEP_SCORE_FDG_FLAT;
	if (strcmp(s, "logdist2") == 0 || strcmp(s, "log_distance2") == 0 ||
		strcmp(s, "log_distance_squared") == 0)
		return HK_BLIND_ESTEP_SCORE_LOGDIST2;
	if (strcmp(s, "dist2") == 0 || strcmp(s, "distance2") == 0 ||
		strcmp(s, "distance_squared") == 0)
		return HK_BLIND_ESTEP_SCORE_DIST2;
	fprintf(stderr, "invalid HK_BLIND_P9016_ESTEP_SCORE_MODE=%s\n", s);
	return -1;
}

static const char *raw_split_confidence_config_suffix(int mstep_graph_mode, int confidence_mode,
													 float confidence_floor)
{
	static char buf[80];
	if (mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_TOP1 &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_TOP2 &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1_CONF_KWEIGHT &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF &&
				mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF_WEIGHTED &&
					mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER &&
					mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_FULL &&
					mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT &&
					mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDWEIGHT &&
					mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDHALF)
		return "";
	switch (confidence_mode) {
	case HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR:
		return "_csentropy";
	case HK_BLIND_RAW_SPLIT_CONF_PMAX_MARGIN:
		return "_cspmax";
	case HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY:
		snprintf(buf, sizeof(buf), "_csfloor%.3g", confidence_floor);
		return buf;
	default:
		return "_csunknown";
	}
}

static const char *raw_split_trans_config_suffix(int mstep_graph_mode, float trans_scale,
												 float trans_confidence_power)
{
	static char buf[80];
	if (mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_TOP1 &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_TOP2 &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1_CONF_KWEIGHT &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF &&
				mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF_WEIGHTED &&
					mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER &&
					mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_FULL &&
					mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT &&
					mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDWEIGHT &&
					mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDHALF)
		return "";
	if (trans_scale == 1.0f && trans_confidence_power == 1.0f)
		return "";
	snprintf(buf, sizeof(buf), "_tcs%.3g_tcp%.3g", trans_scale, trans_confidence_power);
	return buf;
}

static const char *temperature_config_suffix(float temperature_start, float temperature_end)
{
	static char buf[64];
	if (temperature_start == 1.0f && temperature_end == 1.0f)
		return "";
	if (temperature_start == temperature_end)
		snprintf(buf, sizeof(buf), "_T%.3g", temperature_start);
	else
		snprintf(buf, sizeof(buf), "_T%.3gto%.3g", temperature_start, temperature_end);
	return buf;
}

static const char *estep_score_config_suffix(int mode)
{
	switch (mode) {
	case HK_BLIND_ESTEP_SCORE_FDG_FLAT:
		return "";
	case HK_BLIND_ESTEP_SCORE_LOGDIST2:
		return "_elogdist2";
	case HK_BLIND_ESTEP_SCORE_DIST2:
		return "_edist2";
	default:
		return "_eunknown";
	}
}

static const char *homolog_sep_config_suffix(float min_sep_unit, float lambda_sep)
{
	static char buf[64];
	if (min_sep_unit <= 0.0f && lambda_sep <= 0.0f)
		return "";
	snprintf(buf, sizeof(buf), "_sep%.3g_lam%.3g", min_sep_unit, lambda_sep);
	return buf;
}

static const char *chr_sep_config_suffix(float chr_sep_unit, float lambda_chr_sep)
{
	static char buf[64];
	if (chr_sep_unit <= 0.0f && lambda_chr_sep <= 0.0f)
		return "";
	snprintf(buf, sizeof(buf), "_chrsep%.3g_chrlam%.3g", chr_sep_unit, lambda_chr_sep);
	return buf;
}

static const char *contact_k_multiplier_config_suffix(float cis_multiplier,
													  float trans_multiplier)
{
	static char buf[80];
	if (cis_multiplier == 1.0f && trans_multiplier == 1.0f)
		return "";
	if (cis_multiplier == trans_multiplier)
		snprintf(buf, sizeof(buf), "_ck%.3g", cis_multiplier);
	else
		snprintf(buf, sizeof(buf), "_ckc%.3g_ckt%.3g", cis_multiplier, trans_multiplier);
	return buf;
}

static const char *repulsion_block_config_suffix(float k_min)
{
	static char buf[64];
	if (k_min <= 0.0f)
		return "";
	snprintf(buf, sizeof(buf), "_rbk%.3g", k_min);
	return buf;
}

static int env_percent_list_or_default(const char *name, int values[8], int max_values)
{
	const char *s = getenv(name);
	char buf[256];
	char *p;
	int n = 0;
	assert(values);
	assert(max_values > 0);
	if (s == 0 || s[0] == 0) {
		values[0] = 10;
		values[1] = 20;
		values[2] = 50;
		values[3] = 100;
		return 4;
	}
	snprintf(buf, sizeof(buf), "%s", s);
	p = strtok(buf, ",");
	while (p && n < max_values) {
		char *end = 0;
		long v;
		errno = 0;
		v = strtol(p, &end, 10);
		if (errno || end == p || *end != 0 || v < 0 || v > 100) {
			fprintf(stderr, "invalid %s entry=%s; expected integer percent 0..100\n", name, p);
			return -1;
		}
		values[n++] = (int)v;
		p = strtok(0, ",");
	}
	if (p != 0 || n == 0) {
		fprintf(stderr, "invalid %s=%s; too many or empty entries\n", name, s);
		return -1;
	}
	return n;
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

static int init_mode_from_env(const char *name, int fallback)
{
	const char *s = getenv(name);
	if (s == 0 || s[0] == 0)
		return fallback;
	if (strcmp(s, "unphased_scaffold_split") == 0 ||
		strcmp(s, "unphased") == 0 ||
		strcmp(s, "scaffold") == 0)
		return HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT;
	if (strcmp(s, "random_diploid") == 0 ||
		strcmp(s, "random") == 0)
		return HK_BLIND_INIT_RANDOM_DIPLOID;
	if (strcmp(s, "random_haploid_split") == 0 ||
		strcmp(s, "random_haploid") == 0)
		return HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT;
	fprintf(stderr, "invalid %s=%s\n", name, s);
	return -1;
}

static const char *scaffold_source_for_init_mode(int mode)
{
	switch (mode) {
	case HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT: return "unphased_fdg";
	case HK_BLIND_INIT_RANDOM_DIPLOID: return "random_diploid";
	case HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT: return "random_haploid";
	default: return "unknown";
	}
}

static int init_mode_uses_split(int mode)
{
	return mode != HK_BLIND_INIT_RANDOM_DIPLOID;
}

static void minimal_init_config_from_env(struct minimal_init_config *conf)
{
	int default_mode = HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT;
	assert(conf);
	memset(conf, 0, sizeof(*conf));
	conf->mode = init_mode_from_env("HK_BLIND_P9016_INIT_MODE", default_mode);
	conf->seed = env_u64_or_default("HK_BLIND_P9016_INIT_SEED", HK_P9016_INIT_SEED);
	conf->scale = env_float_or_default("HK_BLIND_P9016_INIT_SCALE", HK_P9016_INIT_SCALE, 0.0f);
	conf->eps = env_float_or_default("HK_BLIND_P9016_INIT_EPS", HK_P9016_INIT_EPS, 0.0f);
	conf->noise_scale = env_float_or_default("HK_BLIND_P9016_INIT_NOISE_SCALE",
											 HK_P9016_INIT_NOISE_SCALE, 0.0f);
	conf->scaffold_fdg_n_iter = HK_P9016_SCAFFOLD_FDG_N_ITER;
	if (conf->mode == default_mode &&
		conf->seed == HK_P9016_INIT_SEED &&
		conf->eps == HK_P9016_INIT_EPS &&
		conf->noise_scale == HK_P9016_INIT_NOISE_SCALE) {
		snprintf(conf->config_name, sizeof(conf->config_name), "%s", HK_P9016_CONFIG_NAME);
	} else {
		snprintf(conf->config_name, sizeof(conf->config_name), "%s_%s_seed%llu",
				 HK_P9016_CONFIG_NAME, hk_blind_init_mode_name(conf->mode),
				 (unsigned long long)conf->seed);
	}
}

static void minimal_run_config_init(struct minimal_run_config *run_conf,
									const struct minimal_init_config *init_conf,
									int hard_pc_sample_size_percent, int hard_pc_kind,
									int d_scale_mode,
									int state_weight_mode, int mstep_graph_mode,
									float trans_pmax_min,
									float trans_margin_min, float trans_posterior_power_gamma,
									float imputed_p4_threshold,
										int raw_split_confidence_mode,
											float raw_split_confidence_floor,
											float raw_split_trans_scale,
											float raw_split_trans_confidence_power,
											float raw_split_state_p_min,
											float raw_split_posterior_power_gamma,
											int raw_outlier_enable,
											float raw_outlier_beta_cis,
											float raw_outlier_beta_trans,
											float raw_outlier_prior_cis,
											float raw_outlier_prior_trans,
											float raw_soft_min_q,
											uint64_t raw_split_sample_salt,
											uint64_t hard_pc_seed, float min_sep_unit,
											float lambda_sep,
											float chr_sep_unit, float lambda_chr_sep,
												float contact_k_multiplier_cis,
												float contact_k_multiplier_trans,
												float trans_d_scale_multiplier,
												float temperature_start, float temperature_end,
											float repulsion_block_k_min,
											int estep_score_mode,
											int strict_no_prior_guard,
											const char *input_contact_source)
{
	char base_name[256];
	assert(run_conf);
	assert(init_conf);
	assert(hard_pc_sample_size_percent >= 0 && hard_pc_sample_size_percent <= 100);
	memset(run_conf, 0, sizeof(*run_conf));
	run_conf->init = *init_conf;
	run_conf->hard_pc_sample_size_percent = hard_pc_sample_size_percent;
	run_conf->hard_pc_kind = hard_pc_kind;
	run_conf->d_scale_mode = d_scale_mode;
	run_conf->state_weight_mode = state_weight_mode;
	run_conf->mstep_graph_mode = mstep_graph_mode;
	run_conf->trans_pmax_min = trans_pmax_min;
	run_conf->trans_margin_min = trans_margin_min;
	run_conf->trans_posterior_power_gamma = trans_posterior_power_gamma;
	run_conf->imputed_p4_threshold = imputed_p4_threshold;
	run_conf->raw_split_confidence_mode = raw_split_confidence_mode;
	run_conf->raw_split_confidence_floor = raw_split_confidence_floor;
	run_conf->raw_split_trans_scale = raw_split_trans_scale;
	run_conf->raw_split_trans_confidence_power = raw_split_trans_confidence_power;
	run_conf->raw_split_state_p_min = raw_split_state_p_min;
	run_conf->raw_split_posterior_power_gamma = raw_split_posterior_power_gamma;
	run_conf->raw_outlier_enable = raw_outlier_enable;
	run_conf->raw_outlier_beta_cis = raw_outlier_beta_cis;
	run_conf->raw_outlier_beta_trans = raw_outlier_beta_trans;
	run_conf->raw_outlier_prior_cis = raw_outlier_prior_cis;
	run_conf->raw_outlier_prior_trans = raw_outlier_prior_trans;
	run_conf->raw_soft_min_q = raw_soft_min_q;
	run_conf->raw_split_sample_salt = raw_split_sample_salt;
	run_conf->hard_pc_seed = hard_pc_seed;
	run_conf->min_sep_unit = min_sep_unit;
	run_conf->lambda_sep = lambda_sep;
	run_conf->chr_sep_unit = chr_sep_unit;
	run_conf->lambda_chr_sep = lambda_chr_sep;
	run_conf->contact_k_multiplier_cis = contact_k_multiplier_cis;
	run_conf->contact_k_multiplier_trans = contact_k_multiplier_trans;
	run_conf->trans_d_scale_multiplier = trans_d_scale_multiplier;
	run_conf->temperature_start = temperature_start;
	run_conf->temperature_end = temperature_end;
	run_conf->repulsion_block_k_min = repulsion_block_k_min;
	run_conf->estep_score_mode = estep_score_mode;
	run_conf->strict_no_prior_guard = strict_no_prior_guard;
	snprintf(run_conf->input_contact_source, sizeof(run_conf->input_contact_source),
			 "%s", input_contact_source? input_contact_source : "unknown");
	snprintf(base_name, sizeof(base_name), "%.80s", init_conf->config_name);
	strncat(base_name, d_scale_mode_config_suffix(d_scale_mode),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, state_weight_mode_config_suffix(state_weight_mode, trans_pmax_min,
													   trans_margin_min),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, posterior_gamma_config_suffix(state_weight_mode, trans_posterior_power_gamma),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, mstep_graph_mode_config_suffix(mstep_graph_mode),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, mstep_graph_threshold_config_suffix(mstep_graph_mode,
														  trans_pmax_min,
														  trans_margin_min),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, raw_split_confidence_config_suffix(mstep_graph_mode, raw_split_confidence_mode,
													 raw_split_confidence_floor),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, raw_split_trans_config_suffix(mstep_graph_mode, raw_split_trans_scale,
													 raw_split_trans_confidence_power),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, raw_split_state_p_min_config_suffix(mstep_graph_mode,
														   raw_split_state_p_min),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, raw_split_posterior_gamma_config_suffix(mstep_graph_mode,
															   raw_split_posterior_power_gamma),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, raw_outlier_config_suffix(mstep_graph_mode, raw_outlier_enable,
												 raw_outlier_beta_cis, raw_outlier_beta_trans,
												 raw_outlier_prior_cis, raw_outlier_prior_trans,
												 raw_soft_min_q),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, raw_split_sample_salt_config_suffix(mstep_graph_mode,
														   raw_split_sample_salt),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, contact_k_multiplier_config_suffix(contact_k_multiplier_cis,
														  contact_k_multiplier_trans),
			sizeof(base_name) - strlen(base_name) - 1);
	if (trans_d_scale_multiplier != 1.0f) {
		char buf[64];
		snprintf(buf, sizeof(buf), "_tds%.3g", trans_d_scale_multiplier);
		strncat(base_name, buf, sizeof(base_name) - strlen(base_name) - 1);
	}
	strncat(base_name, temperature_config_suffix(temperature_start, temperature_end),
			sizeof(base_name) - strlen(base_name) - 1);
	strncat(base_name, estep_score_config_suffix(estep_score_mode),
			sizeof(base_name) - strlen(base_name) - 1);
		strncat(base_name, repulsion_block_config_suffix(repulsion_block_k_min),
				sizeof(base_name) - strlen(base_name) - 1);
		strncat(base_name, chr_sep_config_suffix(chr_sep_unit, lambda_chr_sep),
				sizeof(base_name) - strlen(base_name) - 1);
	if (hard_pc_sample_size_percent > 0) {
		const char *kind_suffix = hard_pc_kind == HK_P9016_HARD_PC_KIND_STATE_FREQ? "_statefreq" :
			hard_pc_kind == HK_P9016_HARD_PC_KIND_IMPUTED_P4_TOP? "_impp4top" : "";
		if (min_sep_unit > 0.0f || lambda_sep > 0.0f)
			snprintf(run_conf->config_name, sizeof(run_conf->config_name),
					 "%.190s_hardpc%03d%s_sep%.3g_lam%.3g",
					 base_name, hard_pc_sample_size_percent, kind_suffix,
					 min_sep_unit, lambda_sep);
		else
			snprintf(run_conf->config_name, sizeof(run_conf->config_name), "%.220s_hardpc%03d%s",
					 base_name, hard_pc_sample_size_percent, kind_suffix);
		}
		else
			snprintf(run_conf->config_name, sizeof(run_conf->config_name), "%.190s%s",
					 base_name, homolog_sep_config_suffix(min_sep_unit, lambda_sep));
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

static int init_minimal_coords(struct hk_bmap *bmap, const struct minimal_init_config *init_conf,
							   fvec3_t *haploid, fvec3_t *diploid)
{
	struct hk_fdg_conf scaffold_conf;
	int ret;
	assert(init_conf);
	assert(hk_blind_init_mode_valid(init_conf->mode));
	if (init_conf->mode == HK_BLIND_INIT_RANDOM_DIPLOID)
		return hk_blind_init_random_diploid_coords(bmap, diploid,
												   init_conf->scale,
												   init_conf->seed);
	assert(haploid);
	if (init_conf->mode == HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT) {
		ret = hk_blind_init_random_haploid_scaffold(bmap, haploid,
													init_conf->scale,
													init_conf->seed);
		if (ret != 0) return ret;
		return hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
														 diploid, init_conf->eps,
														 init_conf->noise_scale,
														 init_conf->seed);
	}
	hk_fdg_conf_init(&scaffold_conf);
	scaffold_conf.backend = HK_FDG_BACKEND_CPU;
	scaffold_conf.n_iter = init_conf->scaffold_fdg_n_iter;
	if (hk_blind_init_haploid_scaffold_from_bmap_fdg(bmap, &scaffold_conf, haploid,
													 init_conf->seed) != 0)
		return -1;
	return hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
													 diploid, init_conf->eps,
													 init_conf->noise_scale,
													 init_conf->seed);
}

static void set_minimal_schedule(struct hk_blind_iter_schedule_conf *conf, int n_iter, int relax_steps,
								 float relax_step, float min_sep_unit, float lambda_sep,
								 float chr_sep_unit, float lambda_chr_sep,
								 float contact_k_multiplier_cis,
								 float contact_k_multiplier_trans,
								 float trans_d_scale_multiplier,
								 int d_scale_mode, int state_weight_mode,
								 int mstep_graph_mode,
								 float trans_pmax_min, float trans_margin_min,
								 float trans_posterior_power_gamma,
								 int raw_split_confidence_mode,
								 float raw_split_confidence_floor,
									 float raw_split_trans_scale,
									 float raw_split_trans_confidence_power,
									 float raw_split_state_p_min,
									 float raw_split_posterior_power_gamma,
									 int raw_outlier_enable,
									 float raw_outlier_beta_cis,
									 float raw_outlier_beta_trans,
									 float raw_outlier_prior_cis,
									 float raw_outlier_prior_trans,
									 float raw_soft_min_q,
									 float temperature_start, float temperature_end,
									 float repulsion_block_k_min,
									 int estep_score_mode)
{
	memset(conf, 0, sizeof(*conf));
	conf->n_iter = n_iter;
	conf->base_conf.unit = HK_P9016_UNIT;
	conf->base_conf.d_scale = HK_P9016_D_SCALE;
	conf->base_conf.base_k = 1.0f;
	conf->base_conf.temperature = temperature_start;
	conf->base_conf.rho_train = 1.0f;
	conf->base_conf.rho_train_floor = 0.0f;
	conf->base_conf.min_sep_unit = min_sep_unit;
	conf->base_conf.lambda_sep = lambda_sep;
	conf->base_conf.chr_sep_unit = chr_sep_unit;
	conf->base_conf.lambda_chr_sep = lambda_chr_sep;
	conf->base_conf.relax_step = relax_step;
	conf->base_conf.relax_steps = relax_steps;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.repulsion_block_k_min = repulsion_block_k_min;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = d_scale_mode;
	conf->base_conf.d_scale_eps_count = HK_P9016_D_SCALE_EPS_COUNT;
	conf->base_conf.contact_k_multiplier_cis = contact_k_multiplier_cis;
	conf->base_conf.contact_k_multiplier_trans = contact_k_multiplier_trans;
	conf->base_conf.trans_d_scale_multiplier = trans_d_scale_multiplier;
	conf->base_conf.state_weight_mode = state_weight_mode;
	conf->base_conf.mstep_graph_mode = mstep_graph_mode;
	conf->base_conf.trans_margin_min = trans_margin_min;
	conf->base_conf.trans_pmax_min = trans_pmax_min;
	conf->base_conf.trans_posterior_power_gamma = trans_posterior_power_gamma;
	conf->base_conf.raw_split_confidence_mode = raw_split_confidence_mode;
	conf->base_conf.raw_split_confidence_floor = raw_split_confidence_floor;
	conf->base_conf.raw_split_trans_scale = raw_split_trans_scale;
	conf->base_conf.raw_split_trans_confidence_power = raw_split_trans_confidence_power;
	conf->base_conf.raw_split_state_p_min = raw_split_state_p_min;
	conf->base_conf.raw_split_posterior_power_gamma = raw_split_posterior_power_gamma;
	conf->base_conf.raw_outlier_enable = raw_outlier_enable;
	conf->base_conf.raw_outlier_beta_cis = raw_outlier_beta_cis;
	conf->base_conf.raw_outlier_beta_trans = raw_outlier_beta_trans;
	conf->base_conf.raw_outlier_prior_cis = raw_outlier_prior_cis;
	conf->base_conf.raw_outlier_prior_trans = raw_outlier_prior_trans;
	conf->base_conf.raw_soft_min_q = raw_soft_min_q;
	conf->base_conf.estep_score_mode = estep_score_mode;
	conf->temperature_start = temperature_start;
	conf->temperature_end = temperature_end;
	conf->rho_train_start = 1.0f;
	conf->rho_train_end = 1.0f;
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
								  const struct minimal_run_config *run_conf,
								  struct hk_blind_contact_class_diag *diag_out)
{
	struct hk_blind_wedge_list wedges;
	struct hk_blind_contact_class_diag diag;
	FILE *fp;
	int ret;
	if (run_conf == 0)
		return -1;
	hk_blind_wedge_list_init(&wedges);
	ret = hk_blind_wedge_list_build_mstep_graph_ex(
		&wedges, bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT, 0.0f,
		run_conf->d_scale_mode, HK_P9016_D_SCALE_EPS_COUNT,
		run_conf->contact_k_multiplier_cis, run_conf->contact_k_multiplier_trans,
		run_conf->state_weight_mode,
		run_conf->mstep_graph_mode, run_conf->trans_margin_min,
		run_conf->trans_pmax_min, run_conf->trans_posterior_power_gamma,
		run_conf->raw_split_confidence_mode,
		run_conf->raw_split_confidence_floor,
		run_conf->raw_split_trans_scale,
		run_conf->raw_split_trans_confidence_power,
		run_conf->raw_split_state_p_min,
		run_conf->raw_split_posterior_power_gamma,
		run_conf->raw_outlier_enable,
		run_conf->raw_outlier_beta_cis,
		run_conf->raw_outlier_beta_trans,
		run_conf->raw_outlier_prior_cis,
		run_conf->raw_outlier_prior_trans,
		run_conf->raw_soft_min_q,
		run_conf->temperature_end);
	if (ret == 0)
		hk_blind_wedge_list_apply_trans_d_scale_multiplier(bmap, &wedges,
														   run_conf->trans_d_scale_multiplier);
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
						  const char *force_diag_path, const char *raw_path,
						  int bin_size_bp, int n_iter, int relax_steps, int write_raw,
						  float relax_step,
						  const struct minimal_run_config *run_conf,
						  const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
						  const struct hk_blind_base_k_stats *base_k_stats,
						  const struct hk_blind_iter_loop_diag *loop_diag,
						  const struct hk_blind_contact_class_diag *final_graph_diag,
						  float k_rel_rep_effective)
{
	FILE *fp = fopen(manifest_path, "w");
	const char *status;
	const struct hk_blind_phase_lock_diag *lock_diag;
	const char *init_mode_name;
	const char *scaffold_source;
	double final_refreshed_sum_wedge_k = 0.0;
	int64_t final_refreshed_n_wedges = 0;
	float init_scale_effective;
	float init_eps_effective;
	float init_noise_effective;
	int init_split_params_used;
	int32_t scaffold_fdg_n_iter;
	char label[32];
	if (fp == 0) return -1;
	assert(run_conf);
	assert(final_graph_diag);
	lock_diag = &set->phase_lock_diag;
	final_refreshed_sum_wedge_k = final_graph_diag->sum_wedge_k_cis +
		final_graph_diag->sum_wedge_k_trans;
	final_refreshed_n_wedges = final_graph_diag->n_wedges_cis +
		final_graph_diag->n_wedges_trans;
	resolution_label(bin_size_bp, label);
	{
		const struct minimal_init_config *init_conf = &run_conf->init;
	init_mode_name = hk_blind_init_mode_name(init_conf->mode);
	scaffold_source = scaffold_source_for_init_mode(init_conf->mode);
	init_scale_effective = init_conf->mode == HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT? 0.0f : init_conf->scale;
	init_eps_effective = init_mode_uses_split(init_conf->mode)? init_conf->eps : 0.0f;
	init_noise_effective = init_mode_uses_split(init_conf->mode)? init_conf->noise_scale : 0.0f;
	init_split_params_used = init_mode_uses_split(init_conf->mode)? 1 : 0;
	scaffold_fdg_n_iter = init_conf->mode == HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT?
		init_conf->scaffold_fdg_n_iter : 0;
	}
	status = (loop_diag->n_bad_iter == 0 &&
			  loop_diag->n_relax_nonfinite_iter == 0 &&
			  loop_diag->n_coord_nonfinite == 0)? "OK" : "WARN";
	if (fprintf(fp,
				"key\tvalue\n"
				"sample\tP9016\n"
				"runner_family\tp9016_minimal\n"
				"runner_version\t2026-05-17\n"
				"default_profile\tp9016_minimal_soft_sep_off_v1\n"
				"input_path\t%s\n"
				"input_contact_source\t%s\n"
				"strict_no_prior_guard\t%d\n"
				"output_dir\t%s\n"
				"n_raw\t%d\n"
				"n_bpair\t%d\n"
				"n_beads\t%d\n"
				"resolution\t%d\n"
				"bin_size_bp\t%d\n"
				"resolution_label\t%s\n"
				"bmap_merge_mode\tskip_merge_uniform\n"
				"hk_bmap_gen_skip_merge_arg\t1\n"
				"n_haploid_beads\t%d\n"
				"n_diploid_beads\t%d\n"
				"scan_config\t%s\n"
				"n_iter\t%d\n"
				"unit\t%.9g\n"
				"d_scale\t%.9g\n"
				"base_k_mode\tneighbor_median\n"
				"base_k_effective\t%.9g\n"
				"base_k_min\t%.9g\n"
				"base_k_mean\t%.9g\n"
				"base_k_max\t%.9g\n"
				"base_k_n_nonfinite\t%d\n"
				"legacy_base_k_unused\t1\n"
				"init_mode\t%s\n"
				"init_scale\t%.9g\n"
				"init_eps_effective\t%.9g\n"
				"init_noise_scale_effective\t%.9g\n"
				"init_split_params_used\t%d\n"
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
				"rho_train_floor\t0\n"
				"contact_k_multiplier_cis\t%.9g\n"
				"contact_k_multiplier_trans\t%.9g\n"
				"trans_d_scale_multiplier\t%.9g\n"
				"state_weight_mode\t%s\n"
				"mstep_graph_mode\t%s\n"
				"mstep_uses_stochastic_selection\t%d\n"
				"mstep_uses_weighted_filter\t%d\n"
				"mstep_final_force_probability_weighted\t%d\n"
				"mstep_dscale_probability_weighted\t%d\n"
				"mstep_posterior_squared_risk\t%d\n"
				"estep_score_mode\t%s\n"
				"trans_margin_min\t%.9g\n"
				"trans_pmax_min\t%.9g\n"
				"trans_posterior_power_gamma\t%.9g\n"
				"imputed_p4_threshold\t%.9g\n"
				"raw_split_confidence_mode\t%s\n"
				"raw_split_confidence_floor\t%.9g\n"
				"raw_split_trans_scale\t%.9g\n"
				"raw_split_trans_confidence_power\t%.9g\n"
				"raw_split_state_p_min\t%.9g\n"
				"raw_split_posterior_power_gamma\t%.9g\n"
				"raw_outlier_enable\t%d\n"
				"raw_outlier_beta_cis\t%.9g\n"
				"raw_outlier_beta_trans\t%.9g\n"
				"raw_outlier_prior_cis\t%.9g\n"
				"raw_outlier_prior_trans\t%.9g\n"
				"raw_soft_min_q\t%.9g\n"
				"raw_split_sample_salt\t%llu\n"
				"trans_ambiguous_action\tnone\n"
				"d_scale_mode\t%s\n"
				"d_scale_eps_count\t%.9g\n"
				"count_exposure_mode\tnone\n"
				"count_cap\t%.9g\n"
				"count_exponent\t%.9g\n"
				"coarse_to_fine_map_available\t0\n"
				"lift_from_4mb_available\t0\n"
				"inherited_copy_labels_are_gauge_only\t%d\n"
				"copy_labels_are_gauge_only\t%d\n"
				"uses_phase_labels\t%d\n"
				"hard_positive_control_enabled\t%d\n"
				"hard_positive_control_kind\t%s\n"
				"hard_positive_control_sample_size_requested\t%d\n"
				"hard_positive_control_sample_seed\t%llu\n"
				"hard_positive_control_sample_universe\tall_raw_contacts\n"
				"hard_positive_control_sampled_raw\t%lld\n"
				"hard_positive_control_full_phase_raw\t%lld\n"
				"hard_positive_control_partial_phase_raw\t%lld\n"
				"hard_positive_control_unphased_raw\t%lld\n"
				"hard_positive_control_sampled_full_phase_raw\t%lld\n"
				"hard_positive_control_sampled_partial_phase_raw\t%lld\n"
				"hard_positive_control_sampled_unphased_raw\t%lld\n"
				"hard_positive_control_locked_raw\t%lld\n"
				"hard_positive_control_locked_bpair\t%d\n"
				"hard_positive_control_conflict_raw\t%lld\n"
				"hard_positive_control_conflict_bpair\t%d\n"
				"hard_positive_control_same_bin_locked_raw\t%lld\n"
				"hard_positive_control_same_bin_locked_bpair\t%d\n"
				"hard_positive_control_locked_state_00_bpair\t%d\n"
				"hard_positive_control_locked_state_01_bpair\t%d\n"
				"hard_positive_control_locked_state_10_bpair\t%d\n"
				"hard_positive_control_locked_state_11_bpair\t%d\n"
				"hard_positive_control_locked_state_00_raw\t%lld\n"
				"hard_positive_control_locked_state_01_raw\t%lld\n"
				"hard_positive_control_locked_state_10_raw\t%lld\n"
				"hard_positive_control_locked_state_11_raw\t%lld\n"
				"uses_charm_or_reference\t0\n"
				"uses_charm_for_training\t0\n"
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
				"chr_sep_unit\t%.9g\n"
				"lambda_chr_sep\t%.9g\n"
				"relax_step\t%.9g\n"
					"relax_steps\t%d\n"
					"enable_repulsion\t1\n"
					"repulsion_mode\t%d\n"
					"repulsion_blocking_mode\t%s\n"
					"repulsion_block_k_min\t%.9g\n"
					"repulsion_multiplier\t%.9g\n"
					"k_rel_rep_effective\t%.9g\n"
					"temperature_start\t%.9g\n"
					"temperature_end\t%.9g\n"
					"rho_train_start\t1\n"
					"rho_train_end\t1\n"
				"write_raw_posterior\t%d\n"
				"output_bpair_posterior\t%s\n"
				"output_coords\t%s\n"
				"output_coords_gz\t%s\n"
				"output_loop_diag\t%s\n"
				"output_force_class_diag\t%s\n",
				pairs_path, run_conf->input_contact_source,
				run_conf->strict_no_prior_guard,
				out_dir, set->n_raw, set->n_bpairs, bmap->n_beads,
				bin_size_bp, bin_size_bp, label, bmap->n_beads,
				bmap->n_beads * HK_DIPLOID_N_COPY, run_conf->config_name,
				n_iter, HK_P9016_UNIT,
				HK_P9016_D_SCALE, base_k_stats->mean, base_k_stats->min,
				base_k_stats->mean, base_k_stats->max, base_k_stats->n_nonfinite,
				init_mode_name, init_scale_effective,
				init_eps_effective, init_noise_effective, init_split_params_used,
				(unsigned long long)run_conf->init.seed, scaffold_source, scaffold_fdg_n_iter,
				HK_P9016_D_SCALE_EPS_COUNT, HK_P9016_D_SCALE_EPS_COUNT,
				run_conf->contact_k_multiplier_cis,
				run_conf->contact_k_multiplier_trans,
				run_conf->trans_d_scale_multiplier,
				hk_blind_state_weight_mode_name(run_conf->state_weight_mode),
				hk_blind_mstep_graph_mode_name(run_conf->mstep_graph_mode),
				mstep_graph_mode_uses_stochastic_selection(run_conf->mstep_graph_mode),
				mstep_graph_mode_uses_weighted_filter(run_conf->mstep_graph_mode),
				mstep_graph_mode_uses_final_phased_prob_weight(run_conf->mstep_graph_mode),
				mstep_graph_mode_uses_probability_dscale(run_conf->mstep_graph_mode),
				mstep_graph_mode_has_posterior_squared_risk(run_conf->mstep_graph_mode),
				hk_blind_estep_score_mode_name(run_conf->estep_score_mode),
				run_conf->trans_margin_min, run_conf->trans_pmax_min,
				run_conf->trans_posterior_power_gamma,
				run_conf->imputed_p4_threshold,
				hk_blind_raw_split_confidence_mode_name(run_conf->raw_split_confidence_mode),
				run_conf->raw_split_confidence_floor,
				run_conf->raw_split_trans_scale,
				run_conf->raw_split_trans_confidence_power,
				run_conf->raw_split_state_p_min,
				run_conf->raw_split_posterior_power_gamma,
				run_conf->raw_outlier_enable,
				run_conf->raw_outlier_beta_cis,
				run_conf->raw_outlier_beta_trans,
				run_conf->raw_outlier_prior_cis,
				run_conf->raw_outlier_prior_trans,
				run_conf->raw_soft_min_q,
				(unsigned long long)run_conf->raw_split_sample_salt,
				hk_blind_d_scale_mode_name(run_conf->d_scale_mode),
				HK_P9016_D_SCALE_EPS_COUNT, HK_BLIND_D_SCALE_DEFAULT_COUNT_CAP,
				1.0f / 3.0f,
				run_conf->hard_pc_sample_size_percent > 0? 0 : 1,
				run_conf->hard_pc_sample_size_percent > 0? 0 : 1,
				run_conf->hard_pc_sample_size_percent > 0? 1 : 0,
				run_conf->hard_pc_sample_size_percent > 0? 1 : 0,
				run_conf->hard_pc_sample_size_percent > 0? hard_pc_kind_name(run_conf->hard_pc_kind) : "none",
				run_conf->hard_pc_sample_size_percent,
				(unsigned long long)run_conf->hard_pc_seed,
				(long long)lock_diag->n_sampled_raw,
				(long long)lock_diag->n_full_phase_raw,
				(long long)lock_diag->n_partial_phase_raw,
				(long long)lock_diag->n_unphased_raw,
				(long long)lock_diag->n_sampled_full_phase_raw,
				(long long)lock_diag->n_sampled_partial_phase_raw,
				(long long)lock_diag->n_sampled_unphased_raw,
				(long long)lock_diag->n_locked_raw,
				lock_diag->n_locked_bpair,
				(long long)lock_diag->n_conflict_raw,
				lock_diag->n_conflict_bpair,
				(long long)lock_diag->n_same_bin_locked_raw,
				lock_diag->n_same_bin_locked_bpair,
				lock_diag->state_bpair_count[HK_BLIND_STATE_00],
				lock_diag->state_bpair_count[HK_BLIND_STATE_01],
				lock_diag->state_bpair_count[HK_BLIND_STATE_10],
				lock_diag->state_bpair_count[HK_BLIND_STATE_11],
				(long long)lock_diag->state_raw_count[HK_BLIND_STATE_00],
				(long long)lock_diag->state_raw_count[HK_BLIND_STATE_01],
				(long long)lock_diag->state_raw_count[HK_BLIND_STATE_10],
				(long long)lock_diag->state_raw_count[HK_BLIND_STATE_11],
				set->same_bin_filter_enabled,
				(long long)set->n_raw_same_bin_excluded,
				set->n_bpair_same_bin_excluded, (long long)set->n_raw_cis,
						(long long)set->n_raw_trans, set->n_bpair_cis, set->n_bpair_trans,
							run_conf->min_sep_unit, run_conf->lambda_sep,
							run_conf->chr_sep_unit, run_conf->lambda_chr_sep,
							relax_step, relax_steps, HK_BLIND_REPULSION_CELL,
						run_conf->repulsion_block_k_min > 0.0f? "contact_k_threshold" : "current_edge_blocking",
						run_conf->repulsion_block_k_min,
						HK_P9016_REPULSION_MULTIPLIER, k_rel_rep_effective,
					run_conf->temperature_start, run_conf->temperature_end, write_raw,
					posterior_path, coords_path, coords_gz_path, diag_path, force_diag_path) < 0) {
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
				"posterior_refreshed_after_final_relax\t%d\n"
				"posterior_refresh_temperature\t%.9g\n"
				"posterior_refresh_prior_mode\tuniform\n"
				"posterior_refresh_mean_kl\t%.17g\n"
				"posterior_refresh_top_state_switch_frac\t%.9g\n"
				"posterior_refresh_mean_pU_before\t%.9g\n"
				"posterior_refresh_mean_pU_after\t%.9g\n"
				"final_raw_outlier_n_seen\t%lld\n"
				"final_raw_outlier_n_locked\t%lld\n"
				"final_raw_outlier_n_unlocked\t%lld\n"
				"final_raw_outlier_n_cis\t%lld\n"
				"final_raw_outlier_n_trans\t%lld\n"
				"final_raw_outlier_real_mass_all\t%.17g\n"
				"final_raw_outlier_outlier_mass_all\t%.17g\n"
				"final_raw_outlier_emitted_real_mass_all\t%.17g\n"
				"final_raw_outlier_truncated_real_mass_all\t%.17g\n"
				"n_bad_iter\t%d\n"
				"n_relax_nonfinite_iter\t%d\n"
				"n_coord_nonfinite\t%d\n",
				status, loop_diag->final_mean_entropy, loop_diag->final_mean_pU,
				loop_diag->final_mean_sep, loop_diag->final_min_sep,
				loop_diag->final_max_sep, loop_diag->final_sum_wedge_k,
				loop_diag->final_sum_wedge_k,
				final_refreshed_sum_wedge_k,
				(long long)final_refreshed_n_wedges,
				loop_diag->final_mean_rho_train_bpair,
				loop_diag->final_min_rho_train_bpair,
				loop_diag->final_max_rho_train_bpair,
				loop_diag->posterior_refreshed_after_final_relax,
				loop_diag->posterior_refresh_temperature,
				loop_diag->posterior_refresh_mean_kl,
				loop_diag->posterior_refresh_top_state_switch_frac,
				loop_diag->posterior_refresh_mean_pU_before,
				loop_diag->posterior_refresh_mean_pU_after,
				(long long)loop_diag->final_raw_outlier_diag.n_seen,
				(long long)loop_diag->final_raw_outlier_diag.n_locked,
				(long long)loop_diag->final_raw_outlier_diag.n_unlocked,
				(long long)loop_diag->final_raw_outlier_diag.n_cis,
				(long long)loop_diag->final_raw_outlier_diag.n_trans,
				loop_diag->final_raw_outlier_diag.real_mass_all,
				loop_diag->final_raw_outlier_diag.outlier_mass_all,
				loop_diag->final_raw_outlier_diag.emitted_real_mass_all,
				loop_diag->final_raw_outlier_diag.truncated_real_mass_all,
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
				   "config_name\toutput_dir\tinit_mode\tprior_mode\trho_train_mode\t"
				   "input_contact_source\tstrict_no_prior_guard\t"
				   "rho_train_floor\tcontact_k_multiplier_cis\tcontact_k_multiplier_trans\t"
				   "trans_d_scale_multiplier\t"
				   "state_weight_mode\tmstep_graph_mode\testep_score_mode\t"
				   "mstep_uses_stochastic_selection\tmstep_uses_weighted_filter\t"
				   "mstep_final_force_probability_weighted\t"
				   "mstep_dscale_probability_weighted\t"
				   "mstep_posterior_squared_risk\t"
				   "trans_margin_min\ttrans_pmax_min\t"
				   "trans_posterior_power_gamma\timputed_p4_threshold\t"
						   "raw_split_confidence_mode\traw_split_confidence_floor\t"
						   "raw_split_trans_scale\traw_split_trans_confidence_power\t"
						   "raw_split_state_p_min\traw_split_posterior_power_gamma\t"
						   "raw_outlier_enable\traw_outlier_beta_cis\traw_outlier_beta_trans\t"
						   "raw_outlier_prior_cis\traw_outlier_prior_trans\traw_soft_min_q\t"
						   "raw_split_sample_salt\t"
					   "enable_repulsion\td_scale_mode\tn_iter\t"
				   "relax_steps\tinit_seed\trelax_step\tscan_config\trepulsion_multiplier\t"
				   "k_rel_rep_effective\trepulsion_block_k_min\ttemperature_start\ttemperature_end\t"
				   "rho_train_start\trho_train_end\tn_raw\tn_bpair\tn_beads\tn_raw_cis\t"
				   "n_raw_trans\tn_bpair_cis\tn_bpair_trans\tfinal_mean_entropy\t"
				   "final_mean_pU\tfinal_mean_sep\tfinal_min_sep\tfinal_max_sep\t"
				   "final_sum_wedge_k\tlast_training_sum_wedge_k\t"
					   "final_refreshed_sum_wedge_k\tfinal_refreshed_n_wedges\t"
					   "final_mean_rho_train_bpair\t"
					   "final_min_rho_train_bpair\tfinal_max_rho_train_bpair\t"
						   "final_n_expanded_edges\tfinal_n_split_candidate_pairs\t"
						   "final_n_split_filter1_pairs\tfinal_n_split_filter2_pairs\t"
						   "final_n_split_bmap_pairs\t"
						   "final_n_split_selected_raw\tfinal_n_split_gate_skip_raw\t"
						   "final_n_split_same_bin_skip_raw\tfinal_n_split_locked_skip_raw\t"
						   "final_n_split_state00_raw\tfinal_n_split_state01_raw\t"
						   "final_n_split_state10_raw\tfinal_n_split_state11_raw\t"
						   "final_raw_outlier_n_seen\tfinal_raw_outlier_n_locked\t"
						   "final_raw_outlier_n_unlocked\tfinal_raw_outlier_n_cis\t"
						   "final_raw_outlier_n_trans\tfinal_raw_outlier_n_bad_log_norm\t"
						   "final_raw_outlier_n_all_real_truncated\t"
						   "final_raw_outlier_real_mass_all\tfinal_raw_outlier_outlier_mass_all\t"
						   "final_raw_outlier_emitted_real_mass_all\t"
						   "final_raw_outlier_truncated_real_mass_all\t"
						   "final_raw_outlier_real_mass_cis\tfinal_raw_outlier_outlier_mass_cis\t"
						   "final_raw_outlier_real_mass_trans\tfinal_raw_outlier_outlier_mass_trans\t"
							   "final_contact_energy\tfinal_repulsion_energy\tfinal_backbone_energy\t"
							   "final_n_repulsion_pairs_considered\t"
							   "final_n_repulsion_pairs_blocked\tfinal_n_repulsion_pairs_active\t"
					   "n_same_bin_excluded\tposterior_refreshed_after_final_relax\tn_bad_iter\t"
				   "n_relax_nonfinite_iter\tn_coord_nonfinite\trepulsion_mode\t"
				   "audit_status\twrite_raw_posterior\thard_pc_sample_size\t"
				   "hard_pc_seed\thard_pc_sampled_raw\thard_pc_locked_raw\t"
				   "hard_pc_locked_bpair\thard_pc_conflict_bpair\t"
				   "hard_pc_same_bin_locked_bpair\thard_pc_kind\t"
				   "min_sep_unit\tlambda_sep\tchr_sep_unit\tlambda_chr_sep\n") < 0? -1 : 0;
}

static int append_summary_row(FILE *fp, const struct minimal_result *result, int n_iter,
							  int relax_steps, const struct minimal_run_config *run_conf,
							  float relax_step,
							  int write_raw)
{
	const struct hk_blind_phase_lock_diag *lock_diag;
	assert(run_conf);
	lock_diag = &result->phase_lock_diag;
	return fprintf(fp,
				   "%s\t%s\t%s\tuniform\tconstant\t%s\t%d\t0\t%.9g\t%.9g\t%.9g\t"
					   "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%s\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%llu\t1\t%s\t%d\t%d\t%llu\t%.9g\t"
						   "%s\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t1\t1\t%d\t%d\t%d\t"
						   "%lld\t%lld\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
						   "%.17g\t%.17g\t%.17g\t%lld\t%.9g\t%.9g\t%.9g\t"
						   "%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t"
							   "%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t"
							   "%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t"
							   "%.9g\t%.9g\t%.9g\t%lld\t%lld\t%lld\t%lld\t%d\t%d\t"
							   "%d\t%d\t%d\t%s\t%d\t%d\t%llu\t%lld\t%lld\t%d\t%d\t%d\t%s\t%.9g\t%.9g\t%.9g\t%.9g\n",
				   run_conf->config_name, result->output_dir,
				   hk_blind_init_mode_name(run_conf->init.mode),
				   run_conf->input_contact_source,
				   run_conf->strict_no_prior_guard,
				   run_conf->contact_k_multiplier_cis,
				   run_conf->contact_k_multiplier_trans,
				   run_conf->trans_d_scale_multiplier,
				   hk_blind_state_weight_mode_name(run_conf->state_weight_mode),
				   hk_blind_mstep_graph_mode_name(run_conf->mstep_graph_mode),
				   hk_blind_estep_score_mode_name(run_conf->estep_score_mode),
					   mstep_graph_mode_uses_stochastic_selection(run_conf->mstep_graph_mode),
					   mstep_graph_mode_uses_weighted_filter(run_conf->mstep_graph_mode),
					   mstep_graph_mode_uses_final_phased_prob_weight(run_conf->mstep_graph_mode),
					   mstep_graph_mode_uses_probability_dscale(run_conf->mstep_graph_mode),
					   mstep_graph_mode_has_posterior_squared_risk(run_conf->mstep_graph_mode),
				   run_conf->trans_margin_min, run_conf->trans_pmax_min,
				   run_conf->trans_posterior_power_gamma,
				   run_conf->imputed_p4_threshold,
				   hk_blind_raw_split_confidence_mode_name(run_conf->raw_split_confidence_mode),
				   run_conf->raw_split_confidence_floor,
				   run_conf->raw_split_trans_scale,
				   run_conf->raw_split_trans_confidence_power,
				   run_conf->raw_split_state_p_min,
				   run_conf->raw_split_posterior_power_gamma,
				   run_conf->raw_outlier_enable,
				   run_conf->raw_outlier_beta_cis,
				   run_conf->raw_outlier_beta_trans,
				   run_conf->raw_outlier_prior_cis,
				   run_conf->raw_outlier_prior_trans,
				   run_conf->raw_soft_min_q,
				   (unsigned long long)run_conf->raw_split_sample_salt,
				   hk_blind_d_scale_mode_name(run_conf->d_scale_mode),
				   n_iter, relax_steps,
				   (unsigned long long)run_conf->init.seed, relax_step,
					   run_conf->config_name,
					   HK_P9016_REPULSION_MULTIPLIER, result->k_rel_rep_effective,
					   run_conf->repulsion_block_k_min,
					   run_conf->temperature_start, run_conf->temperature_end,
				   result->n_raw, result->n_bpair, result->n_beads,
				   (long long)result->n_raw_cis, (long long)result->n_raw_trans,
				   result->n_bpair_cis, result->n_bpair_trans,
				   result->final_mean_entropy, result->final_mean_pU,
				   result->final_mean_sep, result->final_min_sep, result->final_max_sep,
				   result->final_sum_wedge_k, result->last_training_sum_wedge_k,
				   result->final_refreshed_sum_wedge_k,
					   (long long)result->final_refreshed_n_wedges,
					   result->final_mean_rho_train_bpair,
					   result->final_min_rho_train_bpair, result->final_max_rho_train_bpair,
					   (long long)result->final_n_expanded_edges,
					   (long long)result->final_n_split_candidate_pairs,
						   (long long)result->final_n_split_filter1_pairs,
						   (long long)result->final_n_split_filter2_pairs,
						   (long long)result->final_n_split_bmap_pairs,
						   (long long)result->final_n_split_selected_raw,
						   (long long)result->final_n_split_gate_skip_raw,
						   (long long)result->final_n_split_same_bin_skip_raw,
						   (long long)result->final_n_split_locked_skip_raw,
						   (long long)result->final_n_split_state_raw_count[HK_BLIND_STATE_00],
						   (long long)result->final_n_split_state_raw_count[HK_BLIND_STATE_01],
						   (long long)result->final_n_split_state_raw_count[HK_BLIND_STATE_10],
							   (long long)result->final_n_split_state_raw_count[HK_BLIND_STATE_11],
							   (long long)result->final_raw_outlier_diag.n_seen,
							   (long long)result->final_raw_outlier_diag.n_locked,
							   (long long)result->final_raw_outlier_diag.n_unlocked,
							   (long long)result->final_raw_outlier_diag.n_cis,
							   (long long)result->final_raw_outlier_diag.n_trans,
							   (long long)result->final_raw_outlier_diag.n_bad_log_norm,
							   (long long)result->final_raw_outlier_diag.n_all_real_truncated,
							   result->final_raw_outlier_diag.real_mass_all,
							   result->final_raw_outlier_diag.outlier_mass_all,
							   result->final_raw_outlier_diag.emitted_real_mass_all,
							   result->final_raw_outlier_diag.truncated_real_mass_all,
							   result->final_raw_outlier_diag.real_mass_cis,
							   result->final_raw_outlier_diag.outlier_mass_cis,
							   result->final_raw_outlier_diag.real_mass_trans,
							   result->final_raw_outlier_diag.outlier_mass_trans,
							   result->final_contact_energy, result->final_repulsion_energy,
					   result->final_backbone_energy,
					   (long long)result->final_n_repulsion_pairs_considered,
					   (long long)result->final_n_repulsion_pairs_blocked,
					   (long long)result->final_n_repulsion_pairs_active,
					   (long long)result->n_same_bin_excluded,
					   result->posterior_refreshed_after_final_relax,
				   result->n_bad_iter, result->n_relax_nonfinite_iter,
				   result->n_coord_nonfinite, HK_BLIND_REPULSION_CELL,
				   result->status_ok? "OK" : "WARN", write_raw,
				   run_conf->hard_pc_sample_size_percent,
				   (unsigned long long)run_conf->hard_pc_seed,
				   (long long)lock_diag->n_sampled_raw,
				   (long long)lock_diag->n_locked_raw,
				   lock_diag->n_locked_bpair,
				   lock_diag->n_conflict_bpair,
				   lock_diag->n_same_bin_locked_bpair,
					   hard_pc_kind_name(run_conf->hard_pc_kind),
					   run_conf->min_sep_unit,
					   run_conf->lambda_sep,
					   run_conf->chr_sep_unit,
					   run_conf->lambda_chr_sep) < 0? -1 : 0;
}

static int run_minimal_config(struct hk_bmap *bmap, const struct hk_blind_pair *raw,
							  const struct hk_pair *pairs, int32_t n_raw,
							  const char *pairs_path, const char *root_dir,
							  int bin_size_bp, int n_iter, int relax_steps, int write_raw,
							  float relax_step,
							  const struct minimal_run_config *run_conf,
							  FILE *summary_fp, struct minimal_result *result)
{
	struct hk_blind_bpair_set *set = 0;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	struct hk_blind_base_k_stats base_k_stats;
	struct hk_blind_contact_class_diag final_graph_diag;
	struct hk_blind_single_iter_diag *per_iter = 0;
	fvec3_t *haploid = 0;
	fvec3_t *diploid = 0;
	int32_t n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	char posterior_path[1024], coords_path[1024], coords_gz_path[1024], diag_path[1024];
	char force_diag_path[1024], raw_path[1024], manifest_path[1024];
	int failed = 0;

	assert(run_conf);
	memset(result, 0, sizeof(*result));
	hk_blind_contact_class_diag_init(&final_graph_diag);
	path_join(result->output_dir, sizeof(result->output_dir), root_dir, run_conf->config_name);
	if (mkdir_if_missing(result->output_dir) != 0)
		return 1;
	path_join(posterior_path, sizeof(posterior_path), result->output_dir, "p9016_full.bpair_posterior.tsv");
	path_join(coords_path, sizeof(coords_path), result->output_dir, "p9016_full.coords.tsv");
	path_join(coords_gz_path, sizeof(coords_gz_path), result->output_dir, "p9016_full.coords.tsv.gz");
	path_join(diag_path, sizeof(diag_path), result->output_dir, "p9016_full.loop_diag.tsv");
	path_join(force_diag_path, sizeof(force_diag_path), result->output_dir, "p9016_full.force_class_diag.tsv");
	path_join(raw_path, sizeof(raw_path), result->output_dir, "p9016_full.raw_posterior.tsv");
	path_join(manifest_path, sizeof(manifest_path), result->output_dir, "p9016_full.manifest.tsv");

	fprintf(stderr, "minimal P9016: build binned contacts\n");
	set = hk_blind_bpair_set_build(bmap, n_raw, raw);
	failed |= check_true("bpair set", set != 0);
	if (failed) goto cleanup;
	failed |= check_i32("base-k mode", hk_blind_bpair_set_apply_base_k_mode(bmap, set,
																			 HK_BLIND_BASE_K_NEIGHBOR_MEDIAN), 0);
	hk_blind_bpair_set_init_uniform_prior(set);
	if (run_conf->hard_pc_kind == HK_P9016_HARD_PC_KIND_STATE_FREQ)
		failed |= check_i32("apply state-frequency positive control locks",
							hk_blind_bpair_set_apply_raw_phase_state_freq_locks(set, pairs, n_raw,
																				run_conf->hard_pc_sample_size_percent,
																				run_conf->hard_pc_seed,
																				&set->phase_lock_diag), 0);
	else if (run_conf->hard_pc_kind == HK_P9016_HARD_PC_KIND_IMPUTED_P4_TOP)
		failed |= check_i32("apply imputed-p4-top positive control locks",
							hk_blind_bpair_set_apply_imputed_p4_top_locks(set, pairs, n_raw,
																		 run_conf->hard_pc_sample_size_percent,
																		 run_conf->hard_pc_seed,
																		 run_conf->imputed_p4_threshold,
																		 &set->phase_lock_diag), 0);
	else
		failed |= check_i32("apply hard positive control locks",
							hk_blind_bpair_set_apply_raw_phase_locks(set, pairs, n_raw,
																	 run_conf->hard_pc_sample_size_percent,
																	 run_conf->hard_pc_seed,
																	 &set->phase_lock_diag), 0);
	hk_blind_bpair_set_base_k_stats(set, &base_k_stats);
	if (failed) goto cleanup;

	haploid = (fvec3_t*)calloc((size_t)bmap->n_beads, sizeof(*haploid));
	diploid = (fvec3_t*)calloc((size_t)n_diploid, sizeof(*diploid));
	per_iter = (struct hk_blind_single_iter_diag*)calloc((size_t)n_iter, sizeof(*per_iter));
	if (haploid == 0 || diploid == 0 || per_iter == 0) {
		failed = 1;
		goto cleanup;
	}
	fprintf(stderr, "minimal P9016: init %s seed=%llu scale=%.9g\n",
			hk_blind_init_mode_name(run_conf->init.mode),
			(unsigned long long)run_conf->init.seed, run_conf->init.scale);
	if (run_conf->hard_pc_sample_size_percent > 0) {
		fprintf(stderr,
				"minimal P9016: hard PC sample=%d%% sampled=%lld full=%lld locked_bpair=%d conflict_bpair=%d\n",
				run_conf->hard_pc_sample_size_percent,
				(long long)set->phase_lock_diag.n_sampled_raw,
				(long long)set->phase_lock_diag.n_sampled_full_phase_raw,
				set->phase_lock_diag.n_locked_bpair,
				set->phase_lock_diag.n_conflict_bpair);
	}
	failed |= check_i32("init coords", init_minimal_coords(bmap, &run_conf->init, haploid, diploid), 0);
	failed |= check_coords_finite(diploid, n_diploid);
	if (failed) goto cleanup;

		fprintf(stderr, "minimal P9016: run soft posterior, min_sep=%.8g lambda_sep=%.8g chr_sep=%.8g lambda_chr_sep=%.8g, n_iter=%d relax_steps=%d relax_step=%.8g\n",
				run_conf->min_sep_unit, run_conf->lambda_sep,
				run_conf->chr_sep_unit, run_conf->lambda_chr_sep,
				n_iter, relax_steps, relax_step);
		fprintf(stderr, "minimal P9016: temperature %.8g -> %.8g\n",
				run_conf->temperature_start, run_conf->temperature_end);
	hk_fdg_conf_init(&fdg_conf);
	fdg_conf.backend = HK_FDG_BACKEND_CPU;
	fdg_conf.k_rel_rep *= HK_P9016_REPULSION_MULTIPLIER;
	set_minimal_schedule(&schedule_conf, n_iter, relax_steps, relax_step,
						 run_conf->min_sep_unit, run_conf->lambda_sep,
						 run_conf->chr_sep_unit, run_conf->lambda_chr_sep,
						 run_conf->contact_k_multiplier_cis,
						 run_conf->contact_k_multiplier_trans,
						 run_conf->trans_d_scale_multiplier,
						 run_conf->d_scale_mode, run_conf->state_weight_mode,
						 run_conf->mstep_graph_mode,
						 run_conf->trans_pmax_min, run_conf->trans_margin_min,
							 run_conf->trans_posterior_power_gamma,
								 run_conf->raw_split_confidence_mode,
						 run_conf->raw_split_confidence_floor,
						 run_conf->raw_split_trans_scale,
						 run_conf->raw_split_trans_confidence_power,
						 run_conf->raw_split_state_p_min,
						 run_conf->raw_split_posterior_power_gamma,
						 run_conf->raw_outlier_enable,
						 run_conf->raw_outlier_beta_cis,
						 run_conf->raw_outlier_beta_trans,
						 run_conf->raw_outlier_prior_cis,
						 run_conf->raw_outlier_prior_trans,
						 run_conf->raw_soft_min_q,
						 run_conf->temperature_start, run_conf->temperature_end,
						 run_conf->repulsion_block_k_min,
								 run_conf->estep_score_mode);
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

	fprintf(stderr, "minimal P9016: write outputs to %s\n", result->output_dir);
	failed |= check_i32("write required outputs",
						write_required_outputs(posterior_path, coords_path, coords_gz_path, diag_path,
											   bmap, set, diploid, &loop_diag), 0);
	failed |= check_i32("write force class diag",
						write_force_class_diag(force_diag_path, &fdg_conf, bmap, set, diploid,
											   n_iter > 0? &per_iter[n_iter - 1].relax_diag : 0,
											   run_conf, &final_graph_diag), 0);
	if (write_raw)
		failed |= check_i32("write raw posterior", write_raw_output(raw_path, raw, n_raw, bmap, set), 0);
	failed |= check_i32("write manifest",
						write_manifest(manifest_path, pairs_path, result->output_dir,
									   posterior_path, coords_path, coords_gz_path, diag_path,
									   force_diag_path, raw_path, bin_size_bp, n_iter, relax_steps,
									   write_raw, relax_step, run_conf, bmap, set, &base_k_stats, &loop_diag,
									   &final_graph_diag,
									   fdg_conf.k_rel_rep), 0);
	if (failed) goto cleanup;

	result->n_raw = set->n_raw;
	result->n_bpair = set->n_bpairs;
	result->n_beads = bmap->n_beads;
	result->phase_lock_diag = set->phase_lock_diag;
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
	result->last_training_sum_wedge_k = loop_diag.final_sum_wedge_k;
	result->final_refreshed_sum_wedge_k = final_graph_diag.sum_wedge_k_cis +
		final_graph_diag.sum_wedge_k_trans;
	result->final_refreshed_n_wedges = final_graph_diag.n_wedges_cis +
		final_graph_diag.n_wedges_trans;
	result->final_mean_rho_train_bpair = loop_diag.final_mean_rho_train_bpair;
	result->final_min_rho_train_bpair = loop_diag.final_min_rho_train_bpair;
	result->final_max_rho_train_bpair = loop_diag.final_max_rho_train_bpair;
	result->final_n_expanded_edges = loop_diag.final_n_expanded_edges;
	result->final_n_split_candidate_pairs = loop_diag.final_n_split_candidate_pairs;
	result->final_n_split_filter1_pairs = loop_diag.final_n_split_filter1_pairs;
	result->final_n_split_filter2_pairs = loop_diag.final_n_split_filter2_pairs;
	result->final_n_split_bmap_pairs = loop_diag.final_n_split_bmap_pairs;
	result->final_n_split_selected_raw = loop_diag.final_n_split_selected_raw;
	result->final_n_split_gate_skip_raw = loop_diag.final_n_split_gate_skip_raw;
	result->final_n_split_same_bin_skip_raw = loop_diag.final_n_split_same_bin_skip_raw;
	result->final_n_split_locked_skip_raw = loop_diag.final_n_split_locked_skip_raw;
	memcpy(result->final_n_split_state_raw_count, loop_diag.final_n_split_state_raw_count,
		   sizeof(result->final_n_split_state_raw_count));
	result->final_raw_outlier_diag = loop_diag.final_raw_outlier_diag;
	result->k_rel_rep_effective = fdg_conf.k_rel_rep;
	if (n_iter > 0) {
		result->final_contact_energy = per_iter[n_iter - 1].relax_diag.final_contact_energy;
		result->final_repulsion_energy = per_iter[n_iter - 1].relax_diag.final_repulsion_energy;
		result->final_backbone_energy = per_iter[n_iter - 1].relax_diag.final_backbone_energy;
		result->final_n_repulsion_pairs_considered = per_iter[n_iter - 1].relax_diag.final_n_repulsion_pairs_considered;
		result->final_n_repulsion_pairs_blocked = per_iter[n_iter - 1].relax_diag.final_n_repulsion_pairs_blocked;
		result->final_n_repulsion_pairs_active = per_iter[n_iter - 1].relax_diag.final_n_repulsion_pairs_active;
	}
	result->posterior_refreshed_after_final_relax = loop_diag.posterior_refreshed_after_final_relax;
	result->n_bad_iter = loop_diag.n_bad_iter;
	result->n_relax_nonfinite_iter = loop_diag.n_relax_nonfinite_iter;
	result->n_coord_nonfinite = loop_diag.n_coord_nonfinite;
	result->status_ok = loop_diag.n_bad_iter == 0 &&
		loop_diag.n_relax_nonfinite_iter == 0 &&
		loop_diag.n_coord_nonfinite == 0;
	failed |= check_i32("append summary", append_summary_row(summary_fp, result, n_iter,
															 relax_steps, run_conf, relax_step, write_raw), 0);
	fflush(summary_fp);
	fprintf(stderr, "minimal P9016: status=%s final_entropy=%.8g final_pU=%.8g\n",
			result->status_ok? "OK" : "WARN", result->final_mean_entropy,
			result->final_mean_pU);

cleanup:
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
	int n_iter = env_int_or_default("HK_BLIND_P9016_MINIMAL_N_ITER",
									HK_P9016_DEFAULT_N_ITER, 1);
	int relax_steps = env_int_or_default("HK_BLIND_P9016_MINIMAL_RELAX_STEPS",
										 HK_P9016_DEFAULT_RELAX_STEPS, 0);
	float relax_step = env_float_or_default("HK_BLIND_P9016_RELAX_STEP",
											HK_P9016_DEFAULT_RELAX_STEP, 0.0f);
	int write_raw = env_flag_enabled("HK_BLIND_WRITE_RAW");
	struct hk_map *map = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	struct minimal_init_config init_conf;
	struct minimal_run_config run_conf;
	uint64_t hard_pc_seed;
	int hard_pc_kind;
	int d_scale_mode;
	int state_weight_mode;
	int mstep_graph_mode;
	float trans_pmax_min, trans_margin_min, trans_posterior_power_gamma;
	float imputed_p4_threshold;
	int raw_split_confidence_mode;
	float raw_split_confidence_floor;
	float raw_split_trans_scale;
	float raw_split_trans_confidence_power;
	float raw_split_state_p_min;
	float raw_split_posterior_power_gamma;
	int raw_outlier_enable;
	float raw_outlier_beta_cis;
	float raw_outlier_beta_trans;
	float raw_outlier_prior_cis;
	float raw_outlier_prior_trans;
	float raw_soft_min_q;
	uint64_t raw_split_sample_salt;
	float min_sep_unit, lambda_sep;
	float chr_sep_unit, lambda_chr_sep;
	float contact_k_multiplier_cis, contact_k_multiplier_trans;
	float trans_d_scale_multiplier;
	float temperature_start, temperature_end;
	float repulsion_block_k_min;
	int estep_score_mode;
	int strict_no_prior_guard;
	const char *input_contact_source;
	int hard_pc_sizes[8];
	int n_hard_pc_sizes, i;
	int32_t n_raw = 0;
	char root_dir[1024], summary_path[1024], bmap_summary_path[1024];
	FILE *summary_fp = 0;
	struct minimal_result result;
	int failed = 0;

	if (make_output_root(root_dir, sizeof(root_dir)) != 0)
		return 1;
	minimal_init_config_from_env(&init_conf);
	hard_pc_seed = env_u64_or_default("HK_BLIND_P9016_HARD_PC_SEED", HK_P9016_HARD_PC_SEED);
	hard_pc_kind = hard_pc_kind_from_env();
	d_scale_mode = d_scale_mode_from_env();
	state_weight_mode = state_weight_mode_from_env();
	mstep_graph_mode = mstep_graph_mode_from_env();
	if (hard_pc_kind < 0 || d_scale_mode < 0 || state_weight_mode < 0 ||
		mstep_graph_mode < 0)
		return 1;
	trans_pmax_min = env_float_or_default("HK_BLIND_P9016_TRANS_PMAX_MIN", 0.0f, 0.0f);
	if (trans_pmax_min > 1.0f) trans_pmax_min = 1.0f;
	trans_margin_min = env_float_or_default("HK_BLIND_P9016_TRANS_MARGIN_MIN", 0.0f, 0.0f);
	if (trans_margin_min > 1.0f) trans_margin_min = 1.0f;
	trans_posterior_power_gamma = env_float_or_default("HK_BLIND_P9016_TRANS_POSTERIOR_POWER_GAMMA", 1.0f, 0.0f);
	if (trans_posterior_power_gamma > 0.0f && trans_posterior_power_gamma < 1.0f)
		trans_posterior_power_gamma = 1.0f;
	imputed_p4_threshold = env_float_or_default("HK_BLIND_P9016_IMPUTED_P4_THRESHOLD",
												HK_P9016_IMPUTED_P4_THRESHOLD, 0.0f);
	if (imputed_p4_threshold > 1.0f) imputed_p4_threshold = 1.0f;
	raw_split_confidence_mode = raw_split_confidence_mode_from_env();
	if (raw_split_confidence_mode < 0)
		return 1;
	raw_split_confidence_floor = env_float_or_default("HK_BLIND_P9016_RAW_SPLIT_CONFIDENCE_FLOOR",
													  0.0f, 0.0f);
	if (raw_split_confidence_floor > 1.0f) raw_split_confidence_floor = 1.0f;
	raw_split_trans_scale = env_float_or_default("HK_BLIND_P9016_RAW_SPLIT_TRANS_SCALE",
												 1.0f, 0.0f);
	raw_split_trans_confidence_power =
		env_float_or_default("HK_BLIND_P9016_RAW_SPLIT_TRANS_CONFIDENCE_POWER",
							 1.0f, 1.0f);
	raw_split_state_p_min = env_float_or_default("HK_BLIND_P9016_RAW_SPLIT_STATE_P_MIN",
												 0.0f, 0.0f);
	if (raw_split_state_p_min > 1.0f) raw_split_state_p_min = 1.0f;
	raw_split_posterior_power_gamma =
		env_float_or_default("HK_BLIND_P9016_RAW_SPLIT_POSTERIOR_POWER_GAMMA",
							 1.0f, 0.0f);
	if (raw_split_posterior_power_gamma > 0.0f &&
		raw_split_posterior_power_gamma < 1.0f)
		raw_split_posterior_power_gamma = 1.0f;
	if (raw_split_posterior_power_gamma > 1.0f &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDWEIGHT &&
		mstep_graph_mode != HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDHALF) {
		fprintf(stderr, "HK_BLIND_P9016_RAW_SPLIT_POSTERIOR_POWER_GAMMA>1 is only supported for mstep_graph_mode=wfilter_kweight, wfilter_kdweight or wfilter_kdhalf\n");
		return 1;
	}
	raw_outlier_enable =
		mstep_graph_mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER;
	if (getenv("HK_BLIND_P9016_RAW_OUTLIER_ENABLE"))
		raw_outlier_enable = env_flag_enabled("HK_BLIND_P9016_RAW_OUTLIER_ENABLE");
	if (getenv("HK_BLIND_RAW_OUTLIER_ENABLE"))
		raw_outlier_enable = env_flag_enabled("HK_BLIND_RAW_OUTLIER_ENABLE");
	raw_outlier_beta_cis =
		env_float_or_default("HK_BLIND_P9016_RAW_OUTLIER_BETA_CIS",
							 HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS, 0.0f);
	raw_outlier_beta_trans =
		env_float_or_default("HK_BLIND_P9016_RAW_OUTLIER_BETA_TRANS",
							 HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS, 0.0f);
	raw_outlier_prior_cis =
		env_float_or_default("HK_BLIND_P9016_RAW_OUTLIER_PRIOR_CIS",
							 HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS, 0.0f);
	raw_outlier_prior_trans =
		env_float_or_default("HK_BLIND_P9016_RAW_OUTLIER_PRIOR_TRANS",
							 HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS, 0.0f);
	raw_soft_min_q =
		env_float_or_default("HK_BLIND_P9016_RAW_SOFT_MIN_Q",
							 HK_BLIND_RAW_SOFT_DEFAULT_MIN_Q, 0.0f);
	if (raw_soft_min_q > 1.0f)
		raw_soft_min_q = 1.0f;
	raw_split_sample_salt = env_u64_or_default("HK_BLIND_P9016_RAW_SPLIT_SAMPLE_SALT", 0);
	if (mstep_graph_mode_uses_stochastic_selection(mstep_graph_mode)) {
		char salt_buf[32];
		snprintf(salt_buf, sizeof(salt_buf), "%llu",
				 (unsigned long long)raw_split_sample_salt);
		setenv("HK_BLIND_RAW_SPLIT_SAMPLE_SALT", salt_buf, 1);
	}
	min_sep_unit = env_float_or_default("HK_BLIND_P9016_MIN_SEP_UNIT", 0.0f, 0.0f);
	lambda_sep = env_float_or_default("HK_BLIND_P9016_LAMBDA_SEP", 0.0f, 0.0f);
	chr_sep_unit = env_float_or_default("HK_BLIND_P9016_CHR_SEP_UNIT", 0.0f, 0.0f);
	lambda_chr_sep = env_float_or_default("HK_BLIND_P9016_LAMBDA_CHR_SEP", 0.0f, 0.0f);
	contact_k_multiplier_cis = env_float_or_default("HK_BLIND_P9016_CONTACT_K_MULTIPLIER_CIS",
													 1.0f, 0.0f);
	contact_k_multiplier_trans = env_float_or_default("HK_BLIND_P9016_CONTACT_K_MULTIPLIER_TRANS",
													   1.0f, 0.0f);
	trans_d_scale_multiplier = env_float_or_default("HK_BLIND_P9016_TRANS_D_SCALE_MULTIPLIER",
													 1.0f, 0.0f);
	if (trans_d_scale_multiplier <= 0.0f)
		trans_d_scale_multiplier = 1.0f;
	temperature_start = env_float_or_default("HK_BLIND_P9016_TEMPERATURE_START", 1.0f, 1e-6f);
	temperature_end = env_float_or_default("HK_BLIND_P9016_TEMPERATURE_END", 1.0f, 1e-6f);
	repulsion_block_k_min = env_float_or_default("HK_BLIND_P9016_REPULSION_BLOCK_K_MIN", 0.0f, 0.0f);
	estep_score_mode = estep_score_mode_from_env();
	if (estep_score_mode < 0 || init_conf.mode < 0)
		return 1;
	strict_no_prior_guard = env_flag_enabled("HK_BLIND_P9016_STRICT_NO_PRIOR");
	input_contact_source = input_contact_source_from_path(pairs_path);
	n_hard_pc_sizes = env_percent_list_or_default("HK_BLIND_P9016_HARD_PC_SIZES",
												  hard_pc_sizes, 8);
	if (n_hard_pc_sizes < 0)
		return 1;
	if (strict_no_prior_guard) {
		if (init_conf.mode != HK_BLIND_INIT_RANDOM_DIPLOID) {
			fprintf(stderr, "strict no-prior requires HK_BLIND_P9016_INIT_MODE=random_diploid\n");
			return 1;
		}
		if (mstep_graph_mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_CIS ||
			mstep_graph_mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_ALL ||
			mstep_graph_mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_ORACLE_ALL) {
			fprintf(stderr, "strict no-prior rejects oracle addback modes because they use imputed p4 addback labels\n");
			return 1;
		}
		if (!strict_no_prior_source_ok(input_contact_source, pairs_path)) {
			fprintf(stderr, "strict no-prior rejects input_contact_source=%s path=%s; only the P9016 raw pairs path or smoke fixture is allowed\n",
					input_contact_source, pairs_path);
			return 1;
		}
		if (n_hard_pc_sizes != 1 || hard_pc_sizes[0] != 0) {
			fprintf(stderr, "strict no-prior requires HK_BLIND_P9016_HARD_PC_SIZES=0\n");
			return 1;
		}
	}
	fprintf(stderr, "minimal P9016: input=%s bin_size_bp=%d output_root=%s init=%s seed=%llu d_scale_mode=%s state_weight_mode=%s mstep_graph_mode=%s estep_score_mode=%s pmax_min=%.4g margin_min=%.4g gamma=%.4g state_p_min=%.4g raw_split_gamma=%.4g outlier=%d beta=%.4g/%.4g prior=%.4g/%.4g minq=%.4g contact_k=%.4g/%.4g trans_d_scale=%.4g temperature=%.4g->%.4g repulsion_block_k_min=%.4g chr_sep=%.4g lambda_chr_sep=%.4g\n",
			pairs_path, bin_size_bp, root_dir, hk_blind_init_mode_name(init_conf.mode),
			(unsigned long long)init_conf.seed, hk_blind_d_scale_mode_name(d_scale_mode),
			hk_blind_state_weight_mode_name(state_weight_mode),
			hk_blind_mstep_graph_mode_name(mstep_graph_mode),
			hk_blind_estep_score_mode_name(estep_score_mode), trans_pmax_min,
			trans_margin_min, trans_posterior_power_gamma,
			raw_split_state_p_min, raw_split_posterior_power_gamma,
			raw_outlier_enable, raw_outlier_beta_cis, raw_outlier_beta_trans,
			raw_outlier_prior_cis, raw_outlier_prior_trans, raw_soft_min_q,
			contact_k_multiplier_cis,
			contact_k_multiplier_trans, trans_d_scale_multiplier,
			temperature_start, temperature_end,
			repulsion_block_k_min, chr_sep_unit, lambda_chr_sep);
	if (read_blind_pairs(pairs_path, bin_size_bp, &map, &bmap, &raw, &n_raw) != 0)
		return 1;
	path_join(bmap_summary_path, sizeof(bmap_summary_path), root_dir, "bmap_summary.tsv");
	failed |= check_i32("write bmap summary", write_bmap_summary(bmap_summary_path, bmap), 0);
	path_join(summary_path, sizeof(summary_path), root_dir, "matrix_summary.tsv");
	summary_fp = fopen(summary_path, "w");
	if (summary_fp == 0) {
		fprintf(stderr, "failed to open summary %s\n", summary_path);
		failed = 1;
		goto cleanup;
	}
	failed |= check_i32("write summary header", write_summary_header(summary_fp), 0);
	for (i = 0; !failed && i < n_hard_pc_sizes; ++i) {
		minimal_run_config_init(&run_conf, &init_conf, hard_pc_sizes[i], hard_pc_kind,
								d_scale_mode, state_weight_mode, mstep_graph_mode,
								trans_pmax_min,
								trans_margin_min, trans_posterior_power_gamma,
								imputed_p4_threshold,
								raw_split_confidence_mode,
								raw_split_confidence_floor,
								raw_split_trans_scale,
								raw_split_trans_confidence_power,
								raw_split_state_p_min,
								raw_split_posterior_power_gamma,
								raw_outlier_enable,
								raw_outlier_beta_cis,
								raw_outlier_beta_trans,
								raw_outlier_prior_cis,
								raw_outlier_prior_trans,
								raw_soft_min_q,
								raw_split_sample_salt,
								hard_pc_seed, min_sep_unit, lambda_sep,
								chr_sep_unit, lambda_chr_sep,
								contact_k_multiplier_cis,
								contact_k_multiplier_trans,
								trans_d_scale_multiplier,
								temperature_start, temperature_end,
								repulsion_block_k_min,
								estep_score_mode,
								strict_no_prior_guard,
								input_contact_source);
		failed |= run_minimal_config(bmap, raw, map->pairs, n_raw, pairs_path, root_dir,
									 bin_size_bp, n_iter, relax_steps, write_raw,
									 relax_step, &run_conf, summary_fp, &result);
	}
cleanup:
	if (summary_fp && fclose(summary_fp) != 0)
		failed = 1;
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (map) hk_map_destroy(map);
	if (failed) {
		fprintf(stderr, "minimal P9016: FAILED\n");
		return 1;
	}
	fprintf(stderr, "minimal P9016: summary %s\n", summary_path);
	return 0;
}
