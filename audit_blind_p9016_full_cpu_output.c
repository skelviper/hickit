#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "hickit.h"

#define HK_AUDIT_LINE_MAX 16384
#define HK_AUDIT_EXAMPLE_MAX 5
#define HK_AUDIT_TOL 1e-4
#define HK_AUDIT_MAX_FIELDS 256

enum manifest_key {
	MK_SAMPLE = 0,
	MK_RUNNER_FAMILY,
	MK_RUNNER_VERSION,
	MK_DEFAULT_PROFILE,
	MK_INPUT_PATH,
	MK_OUTPUT_DIR,
	MK_N_RAW,
	MK_N_BPAIR,
	MK_N_BEADS,
	MK_RESOLUTION,
	MK_N_ITER,
	MK_UNIT,
	MK_D_SCALE,
	MK_BASE_K_MODE,
	MK_BASE_K_EFFECTIVE,
	MK_BASE_K_MIN,
	MK_BASE_K_MEAN,
	MK_BASE_K_MAX,
	MK_BASE_K_N_NONFINITE,
	MK_INIT_MODE,
	MK_PRIOR_MODE,
	MK_PRIOR_EPS,
	MK_PRIOR_INTER_DENSITY,
	MK_PRIOR_OBSERVED_INTER,
	MK_PRIOR_POSSIBLE_INTER,
	MK_PRIOR_N_DISTANCE_BINS,
	MK_PRIOR_N_ALPHA,
	MK_PRIOR_ALPHA_MIN,
	MK_PRIOR_ALPHA_MEDIAN,
	MK_PRIOR_ALPHA_MAX,
	MK_PRIOR_SMOOTHING_METHOD,
	MK_PRIOR_ALPHA_CLAMP_MIN,
	MK_PRIOR_ALPHA_CLAMP_MAX,
	MK_RHO_TRAIN_MODE,
	MK_RHO_TRAIN,
	MK_D_SCALE_MODE,
	MK_D_SCALE_EPS_COUNT,
	MK_D_SCALE_POSTERIOR_GAMMA,
	MK_SAME_BIN_FILTER_ENABLED,
	MK_HELDOUT_ENABLED,
	MK_HELDOUT_FRACTION,
	MK_HELDOUT_SEED,
	MK_HELDOUT_N_INPUT_RAW_TOTAL,
	MK_HELDOUT_N_RAW_TRAIN,
	MK_HELDOUT_N_RAW_HELDOUT,
	MK_HELDOUT_N_BPAIR_HELDOUT,
	MK_HELDOUT_N_BPAIR_EVAL,
	MK_HELDOUT_N_RAW_EVAL,
	MK_HELDOUT_MEAN_EXPECTED_ENERGY,
	MK_HELDOUT_MEAN_MIN_ENERGY,
	MK_HELDOUT_MEAN_ENTROPY,
	MK_HELDOUT_MEAN_PU,
	MK_HELDOUT_MEAN_BEST_NORMALIZED_DISTANCE,
	MK_HELDOUT_SHORT_DISTANCE_FRAC,
	MK_N_RAW_SAME_BIN_EXCLUDED,
	MK_N_BPAIR_SAME_BIN_EXCLUDED,
	MK_RAW_POSTERIOR_SAME_BIN_POLICY,
	MK_MIN_SEP_UNIT,
	MK_LAMBDA_SEP,
	MK_RELAX_STEP,
	MK_RELAX_STEPS,
	MK_TEMPERATURE_START,
	MK_TEMPERATURE_END,
	MK_RHO_TRAIN_START,
	MK_RHO_TRAIN_END,
	MK_INIT_EPS_EFFECTIVE,
	MK_INIT_NOISE_SCALE_EFFECTIVE,
	MK_INIT_SEED,
	MK_ENABLE_REPULSION,
	MK_REPULSION_MODE,
	MK_REPULSION_BLOCKING_MODE,
	MK_WRITE_RAW_POSTERIOR,
	MK_OUTPUT_BPAIR_POSTERIOR,
	MK_OUTPUT_COORDS,
	MK_OUTPUT_LOOP_DIAG,
	MK_FINAL_MEAN_ENTROPY,
	MK_FINAL_MEAN_PU,
	MK_FINAL_SUM_WEDGE_K,
	MK_POSTERIOR_REFRESHED_AFTER_FINAL_RELAX,
	MK_POSTERIOR_REFRESH_TEMPERATURE,
	MK_POSTERIOR_REFRESH_PRIOR_MODE,
	MK_POSTERIOR_REFRESH_MEAN_KL,
	MK_POSTERIOR_REFRESH_TOP_STATE_SWITCH_FRAC,
	MK_POSTERIOR_REFRESH_MEAN_PU_BEFORE,
	MK_POSTERIOR_REFRESH_MEAN_PU_AFTER,
	MK_N_BAD_ITER,
	MK_N_RELAX_NONFINITE_ITER,
	MK_N_COORD_NONFINITE,
	MK_STATUS,
	MK_BASELINE,
	MK_MSTEP_GRAPH_MODE,
	MK_BIN_SIZE_BP,
	MK_N_REQUIRED,
	MK_RHO_TRAIN_FLOOR = MK_N_REQUIRED,
	MK_N_RAW_CIS,
	MK_N_RAW_TRANS,
	MK_N_BPAIR_CIS,
	MK_N_BPAIR_TRANS,
	MK_USES_PHASE_LABELS,
	MK_INIT_SCALE_EFFECTIVE,
	MK_CONFIG_NAME,
	MK_INPUT_CONTACT_SOURCE,
	MK_APPROVED_P9016_RAW_PAIRS_REALPATH,
	MK_INIT_EPS,
	MK_INIT_NOISE_SCALE,
	MK_INIT_SCALE,
	MK_CHR_SEP_UNIT,
	MK_LAMBDA_CHR_SEP,
	MK_LAMBDA_COPYTRACK,
	MK_COPYTRACK_USES_PHASE_LABELS,
	MK_COPYTRACK_USES_CHARM_OR_REFERENCE,
	MK_LAMBDA_GLOBAL_COPYTRACK,
	MK_GLOBAL_COPYTRACK_USES_PHASE_LABELS,
	MK_GLOBAL_COPYTRACK_USES_CHARM_OR_REFERENCE,
	MK_LAMBDA_NORMDIR_COPYTRACK,
	MK_NORMDIR_COPYTRACK_EPS_UNIT,
	MK_NORMDIR_COPYTRACK_USES_PHASE_LABELS,
	MK_NORMDIR_COPYTRACK_USES_CHARM_OR_REFERENCE,
	MK_USES_CHARM_OR_REFERENCE,
	MK_USES_CHARM_FOR_TRAINING,
	MK_REFERENCE_DERIVED_POSITIVE_CONTROL,
	MK_TRANS_CHR_PAIR_PRIOR_MODE,
	MK_TRANS_CHR_PAIR_PRIOR_LAMBDA,
	MK_TRANS_CHR_PAIR_PRIOR_EPS,
	MK_TRANS_CHR_PAIR_PRIOR_POWER,
	MK_TRANS_CHR_PAIR_PRIOR_WARMUP_ITER,
	MK_TRANS_CHR_PAIR_PRIOR_USES_PHASE_LABELS,
	MK_TRANS_CHR_PAIR_PRIOR_USES_CHARM_OR_REFERENCE,
	MK_TRANS_CHR_PAIR_MSTEP_MODE,
	MK_TRANS_CHR_PAIR_MSTEP_LAMBDA,
	MK_TRANS_CHR_PAIR_MSTEP_EPS,
	MK_TRANS_CHR_PAIR_MSTEP_POWER,
	MK_TRANS_CHR_PAIR_MSTEP_WARMUP_ITER,
	MK_TRANS_CHR_PAIR_MSTEP_SCOPE,
	MK_TRANS_CHR_PAIR_MSTEP_APPLICATION_POINT,
	MK_TRANS_CHR_PAIR_MSTEP_USES_PHASE_LABELS,
	MK_TRANS_CHR_PAIR_MSTEP_USES_CHARM_OR_REFERENCE,
	MK_TRANS_CONTACT_SCALING_MODE,
	MK_TRANS_K_MULTIPLIER,
	MK_TRANS_DSCALE_MULTIPLIER,
	MK_TRANS_D_SCALE_POSTERIOR_GAMMA,
	MK_TRANS_CONTACT_SCALING_SCOPE,
	MK_TRANS_CONTACT_SCALING_USES_PHASE_LABELS,
	MK_TRANS_CONTACT_SCALING_USES_CHARM_OR_REFERENCE,
	MK_TRANS_TOP1_MSTEP_MODE,
	MK_TRANS_TOP1_MSTEP_MIN_PMAX,
	MK_TRANS_TOP1_MSTEP_MIN_MARGIN,
	MK_TRANS_TOP1_MSTEP_MIX_WEIGHT,
	MK_TRANS_TOP1_MSTEP_SCOPE,
	MK_TRANS_TOP1_MSTEP_APPLICATION_POINT,
	MK_TRANS_TOP1_MSTEP_USES_PHASE_LABELS,
	MK_TRANS_TOP1_MSTEP_USES_CHARM_OR_REFERENCE,
	MK_TRANS_CALLABLE_ANCHOR_MODE,
	MK_TRANS_CALLABLE_ANCHOR_TOP_FRAC,
	MK_TRANS_CALLABLE_ANCHOR_MIX_WEIGHT,
	MK_TRANS_CALLABLE_ANCHOR_MIN_N_RAW,
	MK_TRANS_CALLABLE_ANCHOR_WARMUP_ITER,
	MK_TRANS_CALLABLE_ANCHOR_USES_PHASE_LABELS,
	MK_TRANS_CALLABLE_ANCHOR_USES_CHARM_OR_REFERENCE,
	MK_TRANS_GATE_MODE,
	MK_TRANS_GATE_MIN_PMAX,
	MK_TRANS_GATE_MIN_MARGIN,
	MK_TRANS_GATE_MIN_NEG_ENTROPY,
	MK_TRANS_GATE_SCOPE,
	MK_TRANS_GATE_APPLICATION_POINT,
	MK_TRANS_GATE_SOURCE,
	MK_TRANS_GATE_USES_PHASE_LABELS,
	MK_TRANS_GATE_USES_CHARM_OR_REFERENCE,
	MK_N_KEYS
};

struct manifest_info {
	uint8_t seen[MK_N_KEYS];
	int has_output_raw_posterior;
	char sample[64];
	char status[64];
	char runner_family[128];
	char runner_version[128];
	char default_profile[128];
	char output_bpair_posterior[1024];
	char output_coords[1024];
	char output_loop_diag[1024];
	char output_raw_posterior[1024];
	char mstep_graph_mode[128];
	int64_t n_raw;
	int64_t n_bpair;
	int64_t n_beads;
	int64_t resolution;
	int64_t bin_size_bp;
	int64_t n_iter;
	int64_t relax_steps;
	int64_t init_seed;
	int64_t enable_repulsion;
	int64_t repulsion_mode;
	int64_t write_raw_posterior;
	int64_t base_k_n_nonfinite;
	int64_t same_bin_filter_enabled;
	int64_t heldout_enabled;
	int64_t heldout_seed;
	int64_t heldout_n_input_raw_total;
	int64_t heldout_n_raw_train;
	int64_t heldout_n_raw_heldout;
	int64_t heldout_n_bpair_heldout;
	int64_t heldout_n_bpair_eval;
	int64_t heldout_n_raw_eval;
	int64_t n_raw_same_bin_excluded;
	int64_t n_bpair_same_bin_excluded;
	int64_t n_raw_cis;
	int64_t n_raw_trans;
	int64_t n_bpair_cis;
	int64_t n_bpair_trans;
	int64_t uses_phase_labels;
	int64_t uses_charm_or_reference;
	int64_t uses_charm_for_training;
	int64_t reference_derived_positive_control;
	int64_t approved_p9016_raw_pairs_realpath;
	int64_t posterior_refreshed_after_final_relax;
	int64_t n_bad_iter;
	int64_t n_relax_nonfinite_iter;
	int64_t n_coord_nonfinite;
	double unit;
	double d_scale;
	double base_k_effective;
	double base_k_min;
	double base_k_mean;
	double base_k_max;
	double d_scale_eps_count;
	double d_scale_posterior_gamma;
	double heldout_fraction;
	double heldout_mean_expected_energy;
	double heldout_mean_min_energy;
	double heldout_mean_entropy;
	double heldout_mean_pU;
	double heldout_mean_best_normalized_distance;
	double heldout_short_distance_frac;
	double min_sep_unit;
	double lambda_sep;
	double lambda_copytrack;
	double relax_step;
	double temperature_start;
	double temperature_end;
	double rho_train_start;
	double rho_train_end;
	double rho_train;
	double rho_train_floor;
	double init_eps_effective;
	double init_noise_scale_effective;
	double init_eps;
	double init_noise_scale;
	double init_scale;
	double init_scale_effective;
	double chr_sep_unit;
	double lambda_chr_sep;
	double lambda_global_copytrack;
	double lambda_normdir_copytrack;
	double normdir_copytrack_eps_unit;
	double prior_eps;
	double prior_inter_density;
	double prior_observed_inter;
	double prior_possible_inter;
	double prior_alpha_min;
	double prior_alpha_median;
	double prior_alpha_max;
	double prior_alpha_clamp_min;
	double prior_alpha_clamp_max;
	double final_mean_entropy;
	double final_mean_pU;
	double final_sum_wedge_k;
	double posterior_refresh_temperature;
	double posterior_refresh_mean_kl;
	double posterior_refresh_top_state_switch_frac;
	double posterior_refresh_mean_pU_before;
	double posterior_refresh_mean_pU_after;
	double trans_chr_pair_prior_lambda;
	double trans_chr_pair_prior_eps;
	double trans_chr_pair_prior_power;
	double trans_chr_pair_mstep_lambda;
	double trans_chr_pair_mstep_eps;
	double trans_chr_pair_mstep_power;
	double trans_k_multiplier;
	double trans_dscale_multiplier;
	double trans_d_scale_posterior_gamma;
	double trans_top1_mstep_min_pmax;
	double trans_top1_mstep_min_margin;
	double trans_top1_mstep_mix_weight;
	double trans_callable_anchor_top_frac;
	double trans_callable_anchor_mix_weight;
	double trans_gate_min_pmax;
	double trans_gate_min_margin;
	double trans_gate_min_neg_entropy;
	int64_t prior_n_distance_bins;
	int64_t prior_n_alpha;
	int64_t trans_chr_pair_prior_warmup_iter;
	int64_t trans_chr_pair_prior_uses_phase_labels;
	int64_t trans_chr_pair_prior_uses_charm_or_reference;
	int64_t trans_chr_pair_mstep_warmup_iter;
	int64_t trans_chr_pair_mstep_uses_phase_labels;
	int64_t trans_chr_pair_mstep_uses_charm_or_reference;
	int64_t trans_contact_scaling_uses_phase_labels;
	int64_t trans_contact_scaling_uses_charm_or_reference;
	int64_t trans_top1_mstep_uses_phase_labels;
	int64_t trans_top1_mstep_uses_charm_or_reference;
	int64_t trans_callable_anchor_min_n_raw;
	int64_t trans_callable_anchor_warmup_iter;
	int64_t trans_callable_anchor_uses_phase_labels;
	int64_t trans_callable_anchor_uses_charm_or_reference;
	int64_t trans_gate_uses_phase_labels;
	int64_t trans_gate_uses_charm_or_reference;
	int64_t copytrack_uses_phase_labels;
	int64_t copytrack_uses_charm_or_reference;
	int64_t global_copytrack_uses_phase_labels;
	int64_t global_copytrack_uses_charm_or_reference;
	int64_t normdir_copytrack_uses_phase_labels;
	int64_t normdir_copytrack_uses_charm_or_reference;
	char base_k_mode[64];
	char init_mode[128];
	char prior_mode[128];
	char trans_chr_pair_prior_mode[128];
	char trans_chr_pair_mstep_mode[128];
	char trans_chr_pair_mstep_scope[128];
	char trans_chr_pair_mstep_application_point[128];
	char trans_contact_scaling_mode[128];
	char trans_contact_scaling_scope[128];
	char trans_top1_mstep_mode[128];
	char trans_top1_mstep_scope[128];
	char trans_top1_mstep_application_point[128];
	char trans_callable_anchor_mode[128];
	char trans_gate_mode[128];
	char trans_gate_scope[128];
	char trans_gate_application_point[128];
	char trans_gate_source[128];
	char prior_smoothing_method[128];
	char rho_train_mode[128];
	char d_scale_mode[128];
	char repulsion_blocking_mode[128];
	char raw_posterior_same_bin_policy[128];
	char posterior_refresh_prior_mode[128];
	char baseline[64];
	char config_name[128];
	char input_contact_source[128];
};

struct file_audit {
	int64_t rows;
	int64_t same_bin_rows;
	int64_t bad_rows;
	int64_t forbidden_hits;
	int n_examples;
	char examples[HK_AUDIT_EXAMPLE_MAX][256];
};

static void audit_init(struct file_audit *audit)
{
	memset(audit, 0, sizeof(*audit));
}

static void trim_line(char *s)
{
	size_t n = strlen(s);
	while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r')) {
		s[--n] = 0;
	}
}

static int check_close(double got, double expected)
{
	double scale = fabs(expected) > 1.0? fabs(expected) : 1.0;
	return fabs(got - expected) <= HK_AUDIT_TOL * scale;
}

static int in_unit_range(double x)
{
	return isfinite(x) && x >= -HK_AUDIT_TOL && x <= 1.0 + HK_AUDIT_TOL;
}

static double expected_base_d_scale_from_n_raw(int n_raw)
{
	return pow((double)n_raw, -1.0 / 3.0);
}

static int recognized_init_mode(const char *mode)
{
	return strcmp(mode, "unphased_scaffold_split") == 0 ||
		   strcmp(mode, "random_diploid") == 0 ||
		   strcmp(mode, "random_haploid_split") == 0 ||
		   strcmp(mode, "toy_split") == 0;
}

static int recognized_rho_train_mode(const char *mode)
{
	return strcmp(mode, "constant") == 0 ||
		   strcmp(mode, "entropy") == 0 ||
		   strcmp(mode, "entropy_with_floor") == 0 ||
		   strcmp(mode, "entropy_cis_constant_trans") == 0 ||
		   strcmp(mode, "entropy_cis_floor_trans") == 0 ||
		   strcmp(mode, "trans_entropy") == 0 ||
		   strcmp(mode, "trans_entropy_with_floor") == 0;
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

static void path_join(char *dst, size_t dst_size, const char *dir, const char *name)
{
	size_t n = strlen(dir);
	if (n > 0 && dir[n - 1] == '/')
		snprintf(dst, dst_size, "%s%s", dir, name);
	else
		snprintf(dst, dst_size, "%s/%s", dir, name);
}

static void add_example(struct file_audit *audit, const char *msg)
{
	if (audit->n_examples >= HK_AUDIT_EXAMPLE_MAX) return;
	snprintf(audit->examples[audit->n_examples], sizeof(audit->examples[audit->n_examples]), "%s", msg);
	++audit->n_examples;
}

static int env_flag_enabled(const char *name)
{
	const char *v = getenv(name);
	return v && strcmp(v, "1") == 0;
}

static void add_example_fmt(struct file_audit *audit, const char *prefix, int64_t row, const char *detail, double value)
{
	char buf[256];
	if (isfinite(value))
		snprintf(buf, sizeof(buf), "%s row %lld %s %.12g", prefix, (long long)row, detail, value);
	else
		snprintf(buf, sizeof(buf), "%s row %lld %s", prefix, (long long)row, detail);
	add_example(audit, buf);
}

static int scan_forbidden_line(struct file_audit *audit, const char *file_label, int64_t row, const char *line)
{
	static const char *forbidden[] = { "phase0", "phase1", "truth" };
	int i, hit = 0;
	for (i = 0; i < (int)(sizeof(forbidden) / sizeof(forbidden[0])); ++i) {
		if (strstr(line, forbidden[i]) != 0) {
			char buf[256];
			++audit->forbidden_hits;
			++hit;
			snprintf(buf, sizeof(buf), "%s row %lld forbidden string %s",
					 file_label, (long long)row, forbidden[i]);
			add_example(audit, buf);
		}
	}
	return hit;
}

static int header_has_token(const char *header, const char *token)
{
	size_t len = strlen(token);
	const char *p = header;
	while ((p = strstr(p, token)) != 0) {
		char before = p == header? '\t' : p[-1];
		char after = p[len];
		if ((p == header || before == '\t') && (after == '\t' || after == '\n' || after == '\r' || after == 0))
			return 1;
		++p;
	}
	return 0;
}

static int check_header_tokens(struct file_audit *audit, const char *file_label, const char *header,
							   const char *const *tokens, int n_tokens)
{
	int failed = 0;
	int i;
	for (i = 0; i < n_tokens; ++i) {
		if (!header_has_token(header, tokens[i])) {
			char buf[256];
			snprintf(buf, sizeof(buf), "%s header missing %s", file_label, tokens[i]);
			add_example(audit, buf);
			failed = 1;
		}
	}
	if (failed) ++audit->bad_rows;
	return failed;
}

static int manifest_key_index(const char *key)
{
	if (strcmp(key, "sample") == 0) return MK_SAMPLE;
	if (strcmp(key, "runner_family") == 0) return MK_RUNNER_FAMILY;
	if (strcmp(key, "runner_version") == 0) return MK_RUNNER_VERSION;
	if (strcmp(key, "default_profile") == 0) return MK_DEFAULT_PROFILE;
	if (strcmp(key, "input_path") == 0) return MK_INPUT_PATH;
	if (strcmp(key, "output_dir") == 0) return MK_OUTPUT_DIR;
	if (strcmp(key, "n_raw") == 0) return MK_N_RAW;
	if (strcmp(key, "n_bpair") == 0) return MK_N_BPAIR;
	if (strcmp(key, "n_beads") == 0) return MK_N_BEADS;
	if (strcmp(key, "resolution") == 0) return MK_RESOLUTION;
	if (strcmp(key, "bin_size_bp") == 0) return MK_BIN_SIZE_BP;
	if (strcmp(key, "n_iter") == 0) return MK_N_ITER;
	if (strcmp(key, "unit") == 0) return MK_UNIT;
	if (strcmp(key, "d_scale") == 0) return MK_D_SCALE;
	if (strcmp(key, "base_k_mode") == 0) return MK_BASE_K_MODE;
	if (strcmp(key, "base_k_effective") == 0) return MK_BASE_K_EFFECTIVE;
	if (strcmp(key, "base_k_min") == 0) return MK_BASE_K_MIN;
	if (strcmp(key, "base_k_mean") == 0) return MK_BASE_K_MEAN;
	if (strcmp(key, "base_k_max") == 0) return MK_BASE_K_MAX;
	if (strcmp(key, "base_k_n_nonfinite") == 0) return MK_BASE_K_N_NONFINITE;
	if (strcmp(key, "init_mode") == 0) return MK_INIT_MODE;
	if (strcmp(key, "prior_mode") == 0) return MK_PRIOR_MODE;
	if (strcmp(key, "prior_eps") == 0) return MK_PRIOR_EPS;
	if (strcmp(key, "prior_inter_density") == 0) return MK_PRIOR_INTER_DENSITY;
	if (strcmp(key, "prior_observed_inter") == 0) return MK_PRIOR_OBSERVED_INTER;
	if (strcmp(key, "prior_possible_inter") == 0) return MK_PRIOR_POSSIBLE_INTER;
	if (strcmp(key, "prior_n_distance_bins") == 0) return MK_PRIOR_N_DISTANCE_BINS;
	if (strcmp(key, "prior_n_alpha") == 0) return MK_PRIOR_N_ALPHA;
	if (strcmp(key, "prior_alpha_min") == 0) return MK_PRIOR_ALPHA_MIN;
	if (strcmp(key, "prior_alpha_median") == 0) return MK_PRIOR_ALPHA_MEDIAN;
	if (strcmp(key, "prior_alpha_max") == 0) return MK_PRIOR_ALPHA_MAX;
	if (strcmp(key, "prior_smoothing_method") == 0) return MK_PRIOR_SMOOTHING_METHOD;
	if (strcmp(key, "prior_alpha_clamp_min") == 0) return MK_PRIOR_ALPHA_CLAMP_MIN;
	if (strcmp(key, "prior_alpha_clamp_max") == 0) return MK_PRIOR_ALPHA_CLAMP_MAX;
	if (strcmp(key, "rho_train_mode") == 0) return MK_RHO_TRAIN_MODE;
	if (strcmp(key, "rho_train") == 0) return MK_RHO_TRAIN;
	if (strcmp(key, "d_scale_mode") == 0) return MK_D_SCALE_MODE;
	if (strcmp(key, "d_scale_eps_count") == 0) return MK_D_SCALE_EPS_COUNT;
	if (strcmp(key, "d_scale_posterior_gamma") == 0) return MK_D_SCALE_POSTERIOR_GAMMA;
	if (strcmp(key, "dscale_posterior_gamma") == 0) return MK_D_SCALE_POSTERIOR_GAMMA;
	if (strcmp(key, "same_bin_filter_enabled") == 0) return MK_SAME_BIN_FILTER_ENABLED;
	if (strcmp(key, "heldout_enabled") == 0) return MK_HELDOUT_ENABLED;
	if (strcmp(key, "heldout_fraction") == 0) return MK_HELDOUT_FRACTION;
	if (strcmp(key, "heldout_seed") == 0) return MK_HELDOUT_SEED;
	if (strcmp(key, "heldout_n_input_raw_total") == 0) return MK_HELDOUT_N_INPUT_RAW_TOTAL;
	if (strcmp(key, "heldout_n_raw_train") == 0) return MK_HELDOUT_N_RAW_TRAIN;
	if (strcmp(key, "heldout_n_raw_heldout") == 0) return MK_HELDOUT_N_RAW_HELDOUT;
	if (strcmp(key, "heldout_n_bpair_heldout") == 0) return MK_HELDOUT_N_BPAIR_HELDOUT;
	if (strcmp(key, "heldout_n_bpair_eval") == 0) return MK_HELDOUT_N_BPAIR_EVAL;
	if (strcmp(key, "heldout_n_raw_eval") == 0) return MK_HELDOUT_N_RAW_EVAL;
	if (strcmp(key, "heldout_mean_expected_energy") == 0) return MK_HELDOUT_MEAN_EXPECTED_ENERGY;
	if (strcmp(key, "heldout_mean_min_energy") == 0) return MK_HELDOUT_MEAN_MIN_ENERGY;
	if (strcmp(key, "heldout_mean_entropy") == 0) return MK_HELDOUT_MEAN_ENTROPY;
	if (strcmp(key, "heldout_mean_pU") == 0) return MK_HELDOUT_MEAN_PU;
	if (strcmp(key, "heldout_mean_best_normalized_distance") == 0) return MK_HELDOUT_MEAN_BEST_NORMALIZED_DISTANCE;
	if (strcmp(key, "heldout_short_distance_frac") == 0) return MK_HELDOUT_SHORT_DISTANCE_FRAC;
	if (strcmp(key, "n_raw_same_bin_excluded") == 0) return MK_N_RAW_SAME_BIN_EXCLUDED;
	if (strcmp(key, "n_bpair_same_bin_excluded") == 0) return MK_N_BPAIR_SAME_BIN_EXCLUDED;
	if (strcmp(key, "raw_posterior_same_bin_policy") == 0) return MK_RAW_POSTERIOR_SAME_BIN_POLICY;
	if (strcmp(key, "min_sep_unit") == 0) return MK_MIN_SEP_UNIT;
	if (strcmp(key, "lambda_sep") == 0) return MK_LAMBDA_SEP;
	if (strcmp(key, "relax_step") == 0) return MK_RELAX_STEP;
	if (strcmp(key, "relax_steps") == 0) return MK_RELAX_STEPS;
	if (strcmp(key, "temperature_start") == 0) return MK_TEMPERATURE_START;
	if (strcmp(key, "temperature_end") == 0) return MK_TEMPERATURE_END;
	if (strcmp(key, "rho_train_start") == 0) return MK_RHO_TRAIN_START;
	if (strcmp(key, "rho_train_end") == 0) return MK_RHO_TRAIN_END;
	if (strcmp(key, "rho_train_floor") == 0) return MK_RHO_TRAIN_FLOOR;
	if (strcmp(key, "baseline") == 0) return MK_BASELINE;
	if (strcmp(key, "mstep_graph_mode") == 0) return MK_MSTEP_GRAPH_MODE;
	if (strcmp(key, "n_raw_cis") == 0) return MK_N_RAW_CIS;
	if (strcmp(key, "n_raw_trans") == 0) return MK_N_RAW_TRANS;
	if (strcmp(key, "n_bpair_cis") == 0) return MK_N_BPAIR_CIS;
	if (strcmp(key, "n_bpair_trans") == 0) return MK_N_BPAIR_TRANS;
	if (strcmp(key, "uses_phase_labels") == 0) return MK_USES_PHASE_LABELS;
	if (strcmp(key, "init_eps_effective") == 0) return MK_INIT_EPS_EFFECTIVE;
	if (strcmp(key, "init_noise_scale_effective") == 0) return MK_INIT_NOISE_SCALE_EFFECTIVE;
	if (strcmp(key, "init_scale_effective") == 0) return MK_INIT_SCALE_EFFECTIVE;
	if (strcmp(key, "config_name") == 0) return MK_CONFIG_NAME;
	if (strcmp(key, "input_contact_source") == 0) return MK_INPUT_CONTACT_SOURCE;
	if (strcmp(key, "approved_p9016_raw_pairs_realpath") == 0) return MK_APPROVED_P9016_RAW_PAIRS_REALPATH;
	if (strcmp(key, "init_eps") == 0) return MK_INIT_EPS;
	if (strcmp(key, "init_noise_scale") == 0) return MK_INIT_NOISE_SCALE;
	if (strcmp(key, "init_scale") == 0) return MK_INIT_SCALE;
	if (strcmp(key, "chr_sep_unit") == 0) return MK_CHR_SEP_UNIT;
	if (strcmp(key, "lambda_chr_sep") == 0) return MK_LAMBDA_CHR_SEP;
	if (strcmp(key, "lambda_copytrack") == 0) return MK_LAMBDA_COPYTRACK;
	if (strcmp(key, "copytrack_uses_phase_labels") == 0) return MK_COPYTRACK_USES_PHASE_LABELS;
	if (strcmp(key, "copytrack_uses_charm_or_reference") == 0) return MK_COPYTRACK_USES_CHARM_OR_REFERENCE;
	if (strcmp(key, "lambda_global_copytrack") == 0) return MK_LAMBDA_GLOBAL_COPYTRACK;
	if (strcmp(key, "global_copytrack_uses_phase_labels") == 0) return MK_GLOBAL_COPYTRACK_USES_PHASE_LABELS;
	if (strcmp(key, "global_copytrack_uses_charm_or_reference") == 0) return MK_GLOBAL_COPYTRACK_USES_CHARM_OR_REFERENCE;
	if (strcmp(key, "lambda_normdir_copytrack") == 0) return MK_LAMBDA_NORMDIR_COPYTRACK;
	if (strcmp(key, "normdir_copytrack_eps_unit") == 0) return MK_NORMDIR_COPYTRACK_EPS_UNIT;
	if (strcmp(key, "normdir_copytrack_uses_phase_labels") == 0) return MK_NORMDIR_COPYTRACK_USES_PHASE_LABELS;
	if (strcmp(key, "normdir_copytrack_uses_charm_or_reference") == 0) return MK_NORMDIR_COPYTRACK_USES_CHARM_OR_REFERENCE;
	if (strcmp(key, "uses_charm_or_reference") == 0) return MK_USES_CHARM_OR_REFERENCE;
	if (strcmp(key, "uses_charm_for_training") == 0) return MK_USES_CHARM_FOR_TRAINING;
	if (strcmp(key, "reference_derived_positive_control") == 0) return MK_REFERENCE_DERIVED_POSITIVE_CONTROL;
	if (strcmp(key, "trans_chr_pair_prior_mode") == 0) return MK_TRANS_CHR_PAIR_PRIOR_MODE;
	if (strcmp(key, "trans_chr_pair_prior_lambda") == 0) return MK_TRANS_CHR_PAIR_PRIOR_LAMBDA;
	if (strcmp(key, "trans_chr_pair_prior_eps") == 0) return MK_TRANS_CHR_PAIR_PRIOR_EPS;
	if (strcmp(key, "trans_chr_pair_prior_power") == 0) return MK_TRANS_CHR_PAIR_PRIOR_POWER;
	if (strcmp(key, "trans_chr_pair_prior_warmup_iter") == 0) return MK_TRANS_CHR_PAIR_PRIOR_WARMUP_ITER;
	if (strcmp(key, "trans_chr_pair_prior_uses_phase_labels") == 0) return MK_TRANS_CHR_PAIR_PRIOR_USES_PHASE_LABELS;
	if (strcmp(key, "trans_chr_pair_prior_uses_charm_or_reference") == 0) return MK_TRANS_CHR_PAIR_PRIOR_USES_CHARM_OR_REFERENCE;
	if (strcmp(key, "trans_chr_pair_mstep_mode") == 0) return MK_TRANS_CHR_PAIR_MSTEP_MODE;
	if (strcmp(key, "trans_chr_pair_mstep_lambda") == 0) return MK_TRANS_CHR_PAIR_MSTEP_LAMBDA;
	if (strcmp(key, "trans_chr_pair_mstep_eps") == 0) return MK_TRANS_CHR_PAIR_MSTEP_EPS;
	if (strcmp(key, "trans_chr_pair_mstep_power") == 0) return MK_TRANS_CHR_PAIR_MSTEP_POWER;
	if (strcmp(key, "trans_chr_pair_mstep_warmup_iter") == 0) return MK_TRANS_CHR_PAIR_MSTEP_WARMUP_ITER;
	if (strcmp(key, "trans_chr_pair_mstep_scope") == 0) return MK_TRANS_CHR_PAIR_MSTEP_SCOPE;
	if (strcmp(key, "trans_chr_pair_mstep_application_point") == 0) return MK_TRANS_CHR_PAIR_MSTEP_APPLICATION_POINT;
	if (strcmp(key, "trans_chr_pair_mstep_uses_phase_labels") == 0) return MK_TRANS_CHR_PAIR_MSTEP_USES_PHASE_LABELS;
	if (strcmp(key, "trans_chr_pair_mstep_uses_charm_or_reference") == 0) return MK_TRANS_CHR_PAIR_MSTEP_USES_CHARM_OR_REFERENCE;
	if (strcmp(key, "trans_contact_scaling_mode") == 0) return MK_TRANS_CONTACT_SCALING_MODE;
	if (strcmp(key, "trans_k_multiplier") == 0) return MK_TRANS_K_MULTIPLIER;
	if (strcmp(key, "trans_dscale_multiplier") == 0) return MK_TRANS_DSCALE_MULTIPLIER;
	if (strcmp(key, "trans_d_scale_posterior_gamma") == 0) return MK_TRANS_D_SCALE_POSTERIOR_GAMMA;
	if (strcmp(key, "trans_contact_scaling_scope") == 0) return MK_TRANS_CONTACT_SCALING_SCOPE;
	if (strcmp(key, "trans_contact_scaling_uses_phase_labels") == 0) return MK_TRANS_CONTACT_SCALING_USES_PHASE_LABELS;
	if (strcmp(key, "trans_contact_scaling_uses_charm_or_reference") == 0) return MK_TRANS_CONTACT_SCALING_USES_CHARM_OR_REFERENCE;
	if (strcmp(key, "trans_top1_mstep_mode") == 0) return MK_TRANS_TOP1_MSTEP_MODE;
	if (strcmp(key, "trans_top1_mstep_min_pmax") == 0) return MK_TRANS_TOP1_MSTEP_MIN_PMAX;
	if (strcmp(key, "trans_top1_mstep_min_margin") == 0) return MK_TRANS_TOP1_MSTEP_MIN_MARGIN;
	if (strcmp(key, "trans_top1_mstep_mix_weight") == 0) return MK_TRANS_TOP1_MSTEP_MIX_WEIGHT;
	if (strcmp(key, "trans_top1_mstep_scope") == 0) return MK_TRANS_TOP1_MSTEP_SCOPE;
	if (strcmp(key, "trans_top1_mstep_application_point") == 0) return MK_TRANS_TOP1_MSTEP_APPLICATION_POINT;
	if (strcmp(key, "trans_top1_mstep_uses_phase_labels") == 0) return MK_TRANS_TOP1_MSTEP_USES_PHASE_LABELS;
	if (strcmp(key, "trans_top1_mstep_uses_charm_or_reference") == 0) return MK_TRANS_TOP1_MSTEP_USES_CHARM_OR_REFERENCE;
	if (strcmp(key, "trans_callable_anchor_mode") == 0) return MK_TRANS_CALLABLE_ANCHOR_MODE;
	if (strcmp(key, "trans_callable_anchor_top_frac") == 0) return MK_TRANS_CALLABLE_ANCHOR_TOP_FRAC;
	if (strcmp(key, "trans_callable_anchor_mix_weight") == 0) return MK_TRANS_CALLABLE_ANCHOR_MIX_WEIGHT;
	if (strcmp(key, "trans_callable_anchor_min_n_raw") == 0) return MK_TRANS_CALLABLE_ANCHOR_MIN_N_RAW;
	if (strcmp(key, "trans_callable_anchor_warmup_iter") == 0) return MK_TRANS_CALLABLE_ANCHOR_WARMUP_ITER;
	if (strcmp(key, "trans_callable_anchor_uses_phase_labels") == 0) return MK_TRANS_CALLABLE_ANCHOR_USES_PHASE_LABELS;
	if (strcmp(key, "trans_callable_anchor_uses_charm_or_reference") == 0) return MK_TRANS_CALLABLE_ANCHOR_USES_CHARM_OR_REFERENCE;
	if (strcmp(key, "trans_gate_mode") == 0) return MK_TRANS_GATE_MODE;
	if (strcmp(key, "trans_gate_min_pmax") == 0) return MK_TRANS_GATE_MIN_PMAX;
	if (strcmp(key, "trans_gate_min_margin") == 0) return MK_TRANS_GATE_MIN_MARGIN;
	if (strcmp(key, "trans_gate_min_neg_entropy") == 0) return MK_TRANS_GATE_MIN_NEG_ENTROPY;
	if (strcmp(key, "trans_gate_scope") == 0) return MK_TRANS_GATE_SCOPE;
	if (strcmp(key, "trans_gate_application_point") == 0) return MK_TRANS_GATE_APPLICATION_POINT;
	if (strcmp(key, "trans_gate_source") == 0) return MK_TRANS_GATE_SOURCE;
	if (strcmp(key, "trans_gate_uses_phase_labels") == 0) return MK_TRANS_GATE_USES_PHASE_LABELS;
	if (strcmp(key, "trans_gate_uses_charm_or_reference") == 0) return MK_TRANS_GATE_USES_CHARM_OR_REFERENCE;
	if (strcmp(key, "init_seed") == 0) return MK_INIT_SEED;
	if (strcmp(key, "enable_repulsion") == 0) return MK_ENABLE_REPULSION;
	if (strcmp(key, "repulsion_mode") == 0) return MK_REPULSION_MODE;
	if (strcmp(key, "repulsion_blocking_mode") == 0) return MK_REPULSION_BLOCKING_MODE;
	if (strcmp(key, "write_raw_posterior") == 0) return MK_WRITE_RAW_POSTERIOR;
	if (strcmp(key, "output_bpair_posterior") == 0) return MK_OUTPUT_BPAIR_POSTERIOR;
	if (strcmp(key, "output_coords") == 0) return MK_OUTPUT_COORDS;
	if (strcmp(key, "output_loop_diag") == 0) return MK_OUTPUT_LOOP_DIAG;
	if (strcmp(key, "final_mean_entropy") == 0) return MK_FINAL_MEAN_ENTROPY;
	if (strcmp(key, "final_mean_pU") == 0) return MK_FINAL_MEAN_PU;
	if (strcmp(key, "final_sum_wedge_k") == 0) return MK_FINAL_SUM_WEDGE_K;
	if (strcmp(key, "posterior_refreshed_after_final_relax") == 0) return MK_POSTERIOR_REFRESHED_AFTER_FINAL_RELAX;
	if (strcmp(key, "posterior_refresh_temperature") == 0) return MK_POSTERIOR_REFRESH_TEMPERATURE;
	if (strcmp(key, "posterior_refresh_prior_mode") == 0) return MK_POSTERIOR_REFRESH_PRIOR_MODE;
	if (strcmp(key, "posterior_refresh_mean_kl") == 0) return MK_POSTERIOR_REFRESH_MEAN_KL;
	if (strcmp(key, "posterior_refresh_top_state_switch_frac") == 0) return MK_POSTERIOR_REFRESH_TOP_STATE_SWITCH_FRAC;
	if (strcmp(key, "posterior_refresh_mean_pU_before") == 0) return MK_POSTERIOR_REFRESH_MEAN_PU_BEFORE;
	if (strcmp(key, "posterior_refresh_mean_pU_after") == 0) return MK_POSTERIOR_REFRESH_MEAN_PU_AFTER;
	if (strcmp(key, "n_bad_iter") == 0) return MK_N_BAD_ITER;
	if (strcmp(key, "n_relax_nonfinite_iter") == 0) return MK_N_RELAX_NONFINITE_ITER;
	if (strcmp(key, "n_coord_nonfinite") == 0) return MK_N_COORD_NONFINITE;
	if (strcmp(key, "status") == 0) return MK_STATUS;
	return -1;
}

static int manifest_optional_current_key(const char *key)
{
	return strcmp(key, "resolution_label") == 0 ||
		   strcmp(key, "refinement_stage") == 0 ||
		   strcmp(key, "refinement_stage_index") == 0 ||
		   strcmp(key, "refinement_n_stages") == 0 ||
		   strcmp(key, "parent_bin_size_bp") == 0 ||
		   strcmp(key, "refinement_init_source") == 0 ||
		   strcmp(key, "refinement_child_offset_step") == 0 ||
			   strcmp(key, "refinement_parent_anchor_k") == 0 ||
			   strcmp(key, "relax_backend") == 0 ||
			   strcmp(key, "scaffold_source") == 0 ||
		   strcmp(key, "scaffold_fdg_n_iter") == 0 ||
			   strcmp(key, "trans_chr_pair_prior_source") == 0 ||
			   strcmp(key, "trans_chr_pair_prior_scope") == 0 ||
				   strcmp(key, "trans_chr_pair_mstep_source") == 0 ||
				   strcmp(key, "trans_chr_pair_mstep_update_rule") == 0 ||
				   strcmp(key, "trans_chr_pair_mstep_aggregate_mode") == 0 ||
				   strcmp(key, "trans_callable_anchor_scope") == 0 ||
				   strcmp(key, "trans_callable_anchor_application_point") == 0 ||
				   strcmp(key, "trans_callable_anchor_score_source") == 0 ||
				   strcmp(key, "trans_callable_anchor_training_rule") == 0 ||
				   strcmp(key, "training_graph_weighted_filter") == 0 ||
			   strcmp(key, "training_graph_probability_weighted") == 0 ||
				   strcmp(key, "edge_k_probability_weighted") == 0 ||
					   strcmp(key, "dscale_mode") == 0 ||
					   strcmp(key, "d_scale_mode_input_string") == 0 ||
					   strcmp(key, "dscale_effective_count_formula") == 0 ||
				   strcmp(key, "dscale_probability_weighted") == 0 ||
				   strcmp(key, "trans_dscale_effective_count_formula") == 0 ||
				   strcmp(key, "trans_dscale_probability_weighted") == 0 ||
				   strcmp(key, "legacy_expected_count_alias_used") == 0 ||
			   strcmp(key, "training_graph_dscale_probability_weighted") == 0 ||
			   strcmp(key, "trans_gate_confidence_unit") == 0 ||
			   strcmp(key, "trans_gate_neg_entropy_formula") == 0 ||
			   strcmp(key, "trans_gate_active_criteria") == 0 ||
			   strcmp(key, "readgroup_mode") == 0 ||
				   strcmp(key, "readgroup_max_segments") == 0 ||
				   strcmp(key, "readgroup_eps") == 0 ||
				   strcmp(key, "readgroup_entries") == 0 ||
				   strcmp(key, "readgroup_groups") == 0 ||
				   strcmp(key, "readgroup_groups_used") == 0 ||
				   strcmp(key, "readgroup_raw_decoded") == 0 ||
				   strcmp(key, "readgroup_groups_skipped_too_large") == 0 ||
				   strcmp(key, "readgroup_groups_skipped_bad") == 0 ||
				   strcmp(key, "readgroup_uses_phase_labels") == 0 ||
				   strcmp(key, "readgroup_uses_charm_or_reference") == 0 ||
				   strcmp(key, "training_graph_mode") == 0 ||
		   strcmp(key, "softall_mode") == 0 ||
		   strcmp(key, "estep_score_mode") == 0 ||
		   strcmp(key, "copy_labels_are_gauge_only") == 0 ||
		   strcmp(key, "reference_training_source_3dg") == 0 ||
		   strcmp(key, "repulsion_multiplier") == 0 ||
		   strcmp(key, "k_rel_rep_effective") == 0 ||
		   strcmp(key, "output_coords_gz") == 0 ||
			   strcmp(key, "output_force_class_diag") == 0 ||
			   strcmp(key, "output_sep_diag") == 0 ||
			   strcmp(key, "output_heldout_diag") == 0 ||
			   strcmp(key, "output_coarse_to_fine_map") == 0 ||
		   strcmp(key, "final_mean_sep") == 0 ||
		   strcmp(key, "final_min_sep") == 0 ||
		   strcmp(key, "final_max_sep") == 0 ||
		   strcmp(key, "last_training_sum_wedge_k") == 0 ||
		   strcmp(key, "final_refreshed_sum_wedge_k") == 0 ||
		   strcmp(key, "final_refreshed_n_wedges") == 0 ||
		   strcmp(key, "final_mean_rho_train_bpair") == 0 ||
		   strcmp(key, "final_min_rho_train_bpair") == 0 ||
		   strcmp(key, "final_max_rho_train_bpair") == 0 ||
		   strcmp(key, "final_contact_energy") == 0 ||
		   strcmp(key, "final_repulsion_energy") == 0 ||
			   strcmp(key, "final_backbone_energy") == 0 ||
			   strcmp(key, "final_sep_energy") == 0 ||
			   strcmp(key, "final_sep_force_l1") == 0 ||
				   strcmp(key, "copytrack_prior_mode") == 0 ||
				   strcmp(key, "copytrack_scope") == 0 ||
				   strcmp(key, "copytrack_vector_definition") == 0 ||
				   strcmp(key, "init_coord_anchor_enabled") == 0 ||
				   strcmp(key, "init_coord_anchor_k") == 0 ||
				   strcmp(key, "init_coord_anchor_mode") == 0 ||
				   strcmp(key, "init_coord_anchor_source") == 0 ||
				   strcmp(key, "init_coord_anchor_map_type") == 0 ||
				   strcmp(key, "init_coord_anchor_uses_phase_labels") == 0 ||
				   strcmp(key, "init_coord_anchor_uses_charm_or_reference") == 0 ||
				   strcmp(key, "global_copytrack_prior_mode") == 0 ||
			   strcmp(key, "global_copytrack_scope") == 0 ||
			   strcmp(key, "global_copytrack_vector_definition") == 0 ||
			   strcmp(key, "global_copytrack_reference") == 0 ||
			   strcmp(key, "global_copytrack_weight_mode") == 0 ||
			   strcmp(key, "normdir_copytrack_prior_mode") == 0 ||
			   strcmp(key, "normdir_copytrack_scope") == 0 ||
			   strcmp(key, "normdir_copytrack_vector_definition") == 0 ||
			   strcmp(key, "normdir_copytrack_short_vector_skip_eps_unit") == 0 ||
			   strcmp(key, "final_copytrack_energy") == 0 ||
			   strcmp(key, "final_copytrack_force_l1") == 0 ||
			   strcmp(key, "final_global_copytrack_energy") == 0 ||
				   strcmp(key, "final_global_copytrack_force_l1") == 0 ||
				   strcmp(key, "final_normdir_copytrack_energy") == 0 ||
				   strcmp(key, "final_normdir_copytrack_force_l1") == 0 ||
				   strcmp(key, "final_anchor_energy") == 0 ||
				   strcmp(key, "final_anchor_force_l1") == 0 ||
				   strcmp(key, "n_copytrack_nonfinite_step") == 0 ||
				   strcmp(key, "n_global_copytrack_nonfinite_step") == 0 ||
				   strcmp(key, "n_normdir_copytrack_nonfinite_step") == 0 ||
				   strcmp(key, "n_anchor_nonfinite_step") == 0 ||
			   strcmp(key, "git_commit") == 0 ||
			   strcmp(key, "git_dirty_count") == 0 ||
			   strcmp(key, "binary_hash") == 0;
		}

static int parse_i64_value(const char *s, int64_t *out)
{
	char *end = 0;
	long long v;
	errno = 0;
	v = strtoll(s, &end, 10);
	if (errno != 0 || end == s || (*end != 0 && *end != '\n' && *end != '\r')) return -1;
	*out = (int64_t)v;
	return 0;
}

static int parse_double_value(const char *s, double *out)
{
	char *end = 0;
	double v;
	errno = 0;
	v = strtod(s, &end);
	if (errno != 0 || end == s || (*end != 0 && *end != '\n' && *end != '\r')) return -1;
	*out = v;
	return 0;
}

static int split_tsv_line(char *line, char **fields, int max_fields)
{
	int n = 0;
	char *p = line;
	assert(fields);
	assert(max_fields > 0);
	while (n < max_fields && p) {
		char *tab = strchr(p, '\t');
		fields[n++] = p;
		if (tab == 0)
			break;
		*tab = 0;
		p = tab + 1;
	}
	return n;
}

static const char *tsv_get_value(char **header_fields, char **value_fields, int n_header, int n_value,
								 const char *name)
{
	int i;
	assert(header_fields);
	assert(value_fields);
	assert(name);
	for (i = 0; i < n_header && i < n_value; ++i)
		if (strcmp(header_fields[i], name) == 0)
			return value_fields[i];
	return 0;
}

static int tsv_get_i32(char **header_fields, char **value_fields, int n_header, int n_value,
					   const char *name, int *out)
{
	int64_t v;
	const char *s = tsv_get_value(header_fields, value_fields, n_header, n_value, name);
	if (s == 0 || parse_i64_value(s, &v) != 0 || v < INT32_MIN || v > INT32_MAX)
		return -1;
	*out = (int)v;
	return 0;
}

static int tsv_get_i64(char **header_fields, char **value_fields, int n_header, int n_value,
					   const char *name, long long *out)
{
	int64_t v;
	const char *s = tsv_get_value(header_fields, value_fields, n_header, n_value, name);
	if (s == 0 || parse_i64_value(s, &v) != 0)
		return -1;
	*out = (long long)v;
	return 0;
}

static int tsv_get_double(char **header_fields, char **value_fields, int n_header, int n_value,
						  const char *name, double *out)
{
	const char *s = tsv_get_value(header_fields, value_fields, n_header, n_value, name);
	if (s == 0 || parse_double_value(s, out) != 0)
		return -1;
	return 0;
}

static int manifest_set_value(struct manifest_info *info, const char *key, const char *value)
{
	int idx = manifest_key_index(key);
	if (idx >= 0) info->seen[idx] = 1;
	switch (idx) {
	case MK_SAMPLE: snprintf(info->sample, sizeof(info->sample), "%s", value); return 0;
	case MK_RUNNER_FAMILY: snprintf(info->runner_family, sizeof(info->runner_family), "%s", value); return 0;
	case MK_RUNNER_VERSION: snprintf(info->runner_version, sizeof(info->runner_version), "%s", value); return 0;
	case MK_DEFAULT_PROFILE: snprintf(info->default_profile, sizeof(info->default_profile), "%s", value); return 0;
	case MK_INPUT_PATH: return 0;
	case MK_OUTPUT_DIR: return 0;
	case MK_STATUS: snprintf(info->status, sizeof(info->status), "%s", value); return 0;
	case MK_BASE_K_MODE: snprintf(info->base_k_mode, sizeof(info->base_k_mode), "%s", value); return 0;
	case MK_INIT_MODE: snprintf(info->init_mode, sizeof(info->init_mode), "%s", value); return 0;
	case MK_PRIOR_MODE: snprintf(info->prior_mode, sizeof(info->prior_mode), "%s", value); return 0;
	case MK_TRANS_CHR_PAIR_PRIOR_MODE: snprintf(info->trans_chr_pair_prior_mode, sizeof(info->trans_chr_pair_prior_mode), "%s", value); return 0;
	case MK_TRANS_CHR_PAIR_MSTEP_MODE: snprintf(info->trans_chr_pair_mstep_mode, sizeof(info->trans_chr_pair_mstep_mode), "%s", value); return 0;
	case MK_TRANS_CHR_PAIR_MSTEP_SCOPE: snprintf(info->trans_chr_pair_mstep_scope, sizeof(info->trans_chr_pair_mstep_scope), "%s", value); return 0;
	case MK_TRANS_CHR_PAIR_MSTEP_APPLICATION_POINT: snprintf(info->trans_chr_pair_mstep_application_point, sizeof(info->trans_chr_pair_mstep_application_point), "%s", value); return 0;
	case MK_TRANS_CONTACT_SCALING_MODE: snprintf(info->trans_contact_scaling_mode, sizeof(info->trans_contact_scaling_mode), "%s", value); return 0;
	case MK_TRANS_CONTACT_SCALING_SCOPE: snprintf(info->trans_contact_scaling_scope, sizeof(info->trans_contact_scaling_scope), "%s", value); return 0;
	case MK_TRANS_TOP1_MSTEP_MODE: snprintf(info->trans_top1_mstep_mode, sizeof(info->trans_top1_mstep_mode), "%s", value); return 0;
		case MK_TRANS_TOP1_MSTEP_SCOPE: snprintf(info->trans_top1_mstep_scope, sizeof(info->trans_top1_mstep_scope), "%s", value); return 0;
		case MK_TRANS_TOP1_MSTEP_APPLICATION_POINT: snprintf(info->trans_top1_mstep_application_point, sizeof(info->trans_top1_mstep_application_point), "%s", value); return 0;
		case MK_TRANS_CALLABLE_ANCHOR_MODE: snprintf(info->trans_callable_anchor_mode, sizeof(info->trans_callable_anchor_mode), "%s", value); return 0;
		case MK_TRANS_GATE_MODE: snprintf(info->trans_gate_mode, sizeof(info->trans_gate_mode), "%s", value); return 0;
	case MK_TRANS_GATE_SCOPE: snprintf(info->trans_gate_scope, sizeof(info->trans_gate_scope), "%s", value); return 0;
	case MK_TRANS_GATE_APPLICATION_POINT: snprintf(info->trans_gate_application_point, sizeof(info->trans_gate_application_point), "%s", value); return 0;
	case MK_TRANS_GATE_SOURCE: snprintf(info->trans_gate_source, sizeof(info->trans_gate_source), "%s", value); return 0;
	case MK_PRIOR_SMOOTHING_METHOD: snprintf(info->prior_smoothing_method, sizeof(info->prior_smoothing_method), "%s", value); return 0;
	case MK_RHO_TRAIN_MODE: snprintf(info->rho_train_mode, sizeof(info->rho_train_mode), "%s", value); return 0;
	case MK_D_SCALE_MODE: snprintf(info->d_scale_mode, sizeof(info->d_scale_mode), "%s", value); return 0;
	case MK_REPULSION_BLOCKING_MODE: snprintf(info->repulsion_blocking_mode, sizeof(info->repulsion_blocking_mode), "%s", value); return 0;
	case MK_RAW_POSTERIOR_SAME_BIN_POLICY: snprintf(info->raw_posterior_same_bin_policy, sizeof(info->raw_posterior_same_bin_policy), "%s", value); return 0;
	case MK_POSTERIOR_REFRESH_PRIOR_MODE: snprintf(info->posterior_refresh_prior_mode, sizeof(info->posterior_refresh_prior_mode), "%s", value); return 0;
	case MK_BASELINE: snprintf(info->baseline, sizeof(info->baseline), "%s", value); return 0;
	case MK_MSTEP_GRAPH_MODE: snprintf(info->mstep_graph_mode, sizeof(info->mstep_graph_mode), "%s", value); return 0;
	case MK_CONFIG_NAME: snprintf(info->config_name, sizeof(info->config_name), "%s", value); return 0;
	case MK_INPUT_CONTACT_SOURCE: snprintf(info->input_contact_source, sizeof(info->input_contact_source), "%s", value); return 0;
	case MK_N_RAW: return parse_i64_value(value, &info->n_raw);
	case MK_N_BPAIR: return parse_i64_value(value, &info->n_bpair);
	case MK_N_BEADS: return parse_i64_value(value, &info->n_beads);
	case MK_RESOLUTION: return parse_i64_value(value, &info->resolution);
	case MK_BIN_SIZE_BP: return parse_i64_value(value, &info->bin_size_bp);
	case MK_N_ITER: return parse_i64_value(value, &info->n_iter);
	case MK_PRIOR_N_DISTANCE_BINS: return parse_i64_value(value, &info->prior_n_distance_bins);
	case MK_PRIOR_N_ALPHA: return parse_i64_value(value, &info->prior_n_alpha);
	case MK_RELAX_STEPS: return parse_i64_value(value, &info->relax_steps);
	case MK_INIT_SEED: return parse_i64_value(value, &info->init_seed);
	case MK_ENABLE_REPULSION: return parse_i64_value(value, &info->enable_repulsion);
	case MK_REPULSION_MODE: return parse_i64_value(value, &info->repulsion_mode);
	case MK_WRITE_RAW_POSTERIOR: return parse_i64_value(value, &info->write_raw_posterior);
		case MK_BASE_K_N_NONFINITE: return parse_i64_value(value, &info->base_k_n_nonfinite);
		case MK_SAME_BIN_FILTER_ENABLED: return parse_i64_value(value, &info->same_bin_filter_enabled);
		case MK_HELDOUT_ENABLED: return parse_i64_value(value, &info->heldout_enabled);
		case MK_HELDOUT_SEED: return parse_i64_value(value, &info->heldout_seed);
		case MK_HELDOUT_N_INPUT_RAW_TOTAL: return parse_i64_value(value, &info->heldout_n_input_raw_total);
		case MK_HELDOUT_N_RAW_TRAIN: return parse_i64_value(value, &info->heldout_n_raw_train);
		case MK_HELDOUT_N_RAW_HELDOUT: return parse_i64_value(value, &info->heldout_n_raw_heldout);
		case MK_HELDOUT_N_BPAIR_HELDOUT: return parse_i64_value(value, &info->heldout_n_bpair_heldout);
		case MK_HELDOUT_N_BPAIR_EVAL: return parse_i64_value(value, &info->heldout_n_bpair_eval);
		case MK_HELDOUT_N_RAW_EVAL: return parse_i64_value(value, &info->heldout_n_raw_eval);
		case MK_N_RAW_SAME_BIN_EXCLUDED: return parse_i64_value(value, &info->n_raw_same_bin_excluded);
	case MK_N_BPAIR_SAME_BIN_EXCLUDED: return parse_i64_value(value, &info->n_bpair_same_bin_excluded);
	case MK_N_RAW_CIS: return parse_i64_value(value, &info->n_raw_cis);
	case MK_N_RAW_TRANS: return parse_i64_value(value, &info->n_raw_trans);
	case MK_N_BPAIR_CIS: return parse_i64_value(value, &info->n_bpair_cis);
	case MK_N_BPAIR_TRANS: return parse_i64_value(value, &info->n_bpair_trans);
	case MK_USES_PHASE_LABELS: return parse_i64_value(value, &info->uses_phase_labels);
	case MK_APPROVED_P9016_RAW_PAIRS_REALPATH: return parse_i64_value(value, &info->approved_p9016_raw_pairs_realpath);
	case MK_USES_CHARM_OR_REFERENCE: return parse_i64_value(value, &info->uses_charm_or_reference);
	case MK_USES_CHARM_FOR_TRAINING: return parse_i64_value(value, &info->uses_charm_for_training);
	case MK_REFERENCE_DERIVED_POSITIVE_CONTROL: return parse_i64_value(value, &info->reference_derived_positive_control);
	case MK_TRANS_CHR_PAIR_PRIOR_USES_PHASE_LABELS: return parse_i64_value(value, &info->trans_chr_pair_prior_uses_phase_labels);
	case MK_TRANS_CHR_PAIR_PRIOR_USES_CHARM_OR_REFERENCE: return parse_i64_value(value, &info->trans_chr_pair_prior_uses_charm_or_reference);
	case MK_TRANS_CHR_PAIR_PRIOR_WARMUP_ITER: return parse_i64_value(value, &info->trans_chr_pair_prior_warmup_iter);
	case MK_TRANS_CHR_PAIR_MSTEP_USES_PHASE_LABELS: return parse_i64_value(value, &info->trans_chr_pair_mstep_uses_phase_labels);
	case MK_TRANS_CHR_PAIR_MSTEP_USES_CHARM_OR_REFERENCE: return parse_i64_value(value, &info->trans_chr_pair_mstep_uses_charm_or_reference);
	case MK_TRANS_CHR_PAIR_MSTEP_WARMUP_ITER: return parse_i64_value(value, &info->trans_chr_pair_mstep_warmup_iter);
	case MK_TRANS_CONTACT_SCALING_USES_PHASE_LABELS: return parse_i64_value(value, &info->trans_contact_scaling_uses_phase_labels);
	case MK_TRANS_CONTACT_SCALING_USES_CHARM_OR_REFERENCE: return parse_i64_value(value, &info->trans_contact_scaling_uses_charm_or_reference);
		case MK_TRANS_TOP1_MSTEP_USES_PHASE_LABELS: return parse_i64_value(value, &info->trans_top1_mstep_uses_phase_labels);
		case MK_TRANS_TOP1_MSTEP_USES_CHARM_OR_REFERENCE: return parse_i64_value(value, &info->trans_top1_mstep_uses_charm_or_reference);
		case MK_TRANS_CALLABLE_ANCHOR_USES_PHASE_LABELS: return parse_i64_value(value, &info->trans_callable_anchor_uses_phase_labels);
		case MK_TRANS_CALLABLE_ANCHOR_USES_CHARM_OR_REFERENCE: return parse_i64_value(value, &info->trans_callable_anchor_uses_charm_or_reference);
		case MK_TRANS_CALLABLE_ANCHOR_MIN_N_RAW: return parse_i64_value(value, &info->trans_callable_anchor_min_n_raw);
		case MK_TRANS_CALLABLE_ANCHOR_WARMUP_ITER: return parse_i64_value(value, &info->trans_callable_anchor_warmup_iter);
		case MK_TRANS_GATE_USES_PHASE_LABELS: return parse_i64_value(value, &info->trans_gate_uses_phase_labels);
	case MK_TRANS_GATE_USES_CHARM_OR_REFERENCE: return parse_i64_value(value, &info->trans_gate_uses_charm_or_reference);
	case MK_POSTERIOR_REFRESHED_AFTER_FINAL_RELAX: return parse_i64_value(value, &info->posterior_refreshed_after_final_relax);
	case MK_N_BAD_ITER: return parse_i64_value(value, &info->n_bad_iter);
	case MK_N_RELAX_NONFINITE_ITER: return parse_i64_value(value, &info->n_relax_nonfinite_iter);
	case MK_N_COORD_NONFINITE: return parse_i64_value(value, &info->n_coord_nonfinite);
	case MK_UNIT: return parse_double_value(value, &info->unit);
	case MK_D_SCALE: return parse_double_value(value, &info->d_scale);
	case MK_BASE_K_EFFECTIVE: return parse_double_value(value, &info->base_k_effective);
	case MK_BASE_K_MIN: return parse_double_value(value, &info->base_k_min);
	case MK_BASE_K_MEAN: return parse_double_value(value, &info->base_k_mean);
	case MK_BASE_K_MAX: return parse_double_value(value, &info->base_k_max);
	case MK_D_SCALE_EPS_COUNT: return parse_double_value(value, &info->d_scale_eps_count);
	case MK_D_SCALE_POSTERIOR_GAMMA: return parse_double_value(value, &info->d_scale_posterior_gamma);
	case MK_HELDOUT_FRACTION: return parse_double_value(value, &info->heldout_fraction);
		case MK_HELDOUT_MEAN_EXPECTED_ENERGY: return parse_double_value(value, &info->heldout_mean_expected_energy);
		case MK_HELDOUT_MEAN_MIN_ENERGY: return parse_double_value(value, &info->heldout_mean_min_energy);
		case MK_HELDOUT_MEAN_ENTROPY: return parse_double_value(value, &info->heldout_mean_entropy);
		case MK_HELDOUT_MEAN_PU: return parse_double_value(value, &info->heldout_mean_pU);
		case MK_HELDOUT_MEAN_BEST_NORMALIZED_DISTANCE: return parse_double_value(value, &info->heldout_mean_best_normalized_distance);
		case MK_HELDOUT_SHORT_DISTANCE_FRAC: return parse_double_value(value, &info->heldout_short_distance_frac);
		case MK_MIN_SEP_UNIT: return parse_double_value(value, &info->min_sep_unit);
	case MK_LAMBDA_SEP: return parse_double_value(value, &info->lambda_sep);
	case MK_LAMBDA_COPYTRACK: return parse_double_value(value, &info->lambda_copytrack);
	case MK_COPYTRACK_USES_PHASE_LABELS: return parse_i64_value(value, &info->copytrack_uses_phase_labels);
	case MK_COPYTRACK_USES_CHARM_OR_REFERENCE: return parse_i64_value(value, &info->copytrack_uses_charm_or_reference);
	case MK_LAMBDA_GLOBAL_COPYTRACK: return parse_double_value(value, &info->lambda_global_copytrack);
	case MK_GLOBAL_COPYTRACK_USES_PHASE_LABELS: return parse_i64_value(value, &info->global_copytrack_uses_phase_labels);
	case MK_GLOBAL_COPYTRACK_USES_CHARM_OR_REFERENCE: return parse_i64_value(value, &info->global_copytrack_uses_charm_or_reference);
	case MK_LAMBDA_NORMDIR_COPYTRACK: return parse_double_value(value, &info->lambda_normdir_copytrack);
	case MK_NORMDIR_COPYTRACK_EPS_UNIT: return parse_double_value(value, &info->normdir_copytrack_eps_unit);
	case MK_NORMDIR_COPYTRACK_USES_PHASE_LABELS: return parse_i64_value(value, &info->normdir_copytrack_uses_phase_labels);
	case MK_NORMDIR_COPYTRACK_USES_CHARM_OR_REFERENCE: return parse_i64_value(value, &info->normdir_copytrack_uses_charm_or_reference);
	case MK_TRANS_CHR_PAIR_PRIOR_LAMBDA: return parse_double_value(value, &info->trans_chr_pair_prior_lambda);
	case MK_TRANS_CHR_PAIR_PRIOR_EPS: return parse_double_value(value, &info->trans_chr_pair_prior_eps);
	case MK_TRANS_CHR_PAIR_PRIOR_POWER: return parse_double_value(value, &info->trans_chr_pair_prior_power);
	case MK_TRANS_CHR_PAIR_MSTEP_LAMBDA: return parse_double_value(value, &info->trans_chr_pair_mstep_lambda);
	case MK_TRANS_CHR_PAIR_MSTEP_EPS: return parse_double_value(value, &info->trans_chr_pair_mstep_eps);
	case MK_TRANS_CHR_PAIR_MSTEP_POWER: return parse_double_value(value, &info->trans_chr_pair_mstep_power);
	case MK_TRANS_K_MULTIPLIER: return parse_double_value(value, &info->trans_k_multiplier);
	case MK_TRANS_DSCALE_MULTIPLIER: return parse_double_value(value, &info->trans_dscale_multiplier);
	case MK_TRANS_D_SCALE_POSTERIOR_GAMMA: return parse_double_value(value, &info->trans_d_scale_posterior_gamma);
	case MK_TRANS_TOP1_MSTEP_MIN_PMAX: return parse_double_value(value, &info->trans_top1_mstep_min_pmax);
		case MK_TRANS_TOP1_MSTEP_MIN_MARGIN: return parse_double_value(value, &info->trans_top1_mstep_min_margin);
		case MK_TRANS_TOP1_MSTEP_MIX_WEIGHT: return parse_double_value(value, &info->trans_top1_mstep_mix_weight);
		case MK_TRANS_CALLABLE_ANCHOR_TOP_FRAC: return parse_double_value(value, &info->trans_callable_anchor_top_frac);
		case MK_TRANS_CALLABLE_ANCHOR_MIX_WEIGHT: return parse_double_value(value, &info->trans_callable_anchor_mix_weight);
		case MK_TRANS_GATE_MIN_PMAX: return parse_double_value(value, &info->trans_gate_min_pmax);
	case MK_TRANS_GATE_MIN_MARGIN: return parse_double_value(value, &info->trans_gate_min_margin);
	case MK_TRANS_GATE_MIN_NEG_ENTROPY: return parse_double_value(value, &info->trans_gate_min_neg_entropy);
	case MK_RELAX_STEP: return parse_double_value(value, &info->relax_step);
	case MK_TEMPERATURE_START: return parse_double_value(value, &info->temperature_start);
	case MK_TEMPERATURE_END: return parse_double_value(value, &info->temperature_end);
	case MK_RHO_TRAIN_START: return parse_double_value(value, &info->rho_train_start);
	case MK_RHO_TRAIN_END: return parse_double_value(value, &info->rho_train_end);
	case MK_RHO_TRAIN: return parse_double_value(value, &info->rho_train);
	case MK_RHO_TRAIN_FLOOR: return parse_double_value(value, &info->rho_train_floor);
	case MK_INIT_EPS_EFFECTIVE: return parse_double_value(value, &info->init_eps_effective);
	case MK_INIT_NOISE_SCALE_EFFECTIVE: return parse_double_value(value, &info->init_noise_scale_effective);
	case MK_INIT_SCALE_EFFECTIVE: return parse_double_value(value, &info->init_scale_effective);
	case MK_INIT_EPS: return parse_double_value(value, &info->init_eps);
	case MK_INIT_NOISE_SCALE: return parse_double_value(value, &info->init_noise_scale);
	case MK_INIT_SCALE: return parse_double_value(value, &info->init_scale);
	case MK_CHR_SEP_UNIT: return parse_double_value(value, &info->chr_sep_unit);
	case MK_LAMBDA_CHR_SEP: return parse_double_value(value, &info->lambda_chr_sep);
	case MK_PRIOR_EPS: return parse_double_value(value, &info->prior_eps);
	case MK_PRIOR_INTER_DENSITY: return parse_double_value(value, &info->prior_inter_density);
	case MK_PRIOR_OBSERVED_INTER: return parse_double_value(value, &info->prior_observed_inter);
	case MK_PRIOR_POSSIBLE_INTER: return parse_double_value(value, &info->prior_possible_inter);
	case MK_PRIOR_ALPHA_MIN: return parse_double_value(value, &info->prior_alpha_min);
	case MK_PRIOR_ALPHA_MEDIAN: return parse_double_value(value, &info->prior_alpha_median);
	case MK_PRIOR_ALPHA_MAX: return parse_double_value(value, &info->prior_alpha_max);
	case MK_PRIOR_ALPHA_CLAMP_MIN: return parse_double_value(value, &info->prior_alpha_clamp_min);
	case MK_PRIOR_ALPHA_CLAMP_MAX: return parse_double_value(value, &info->prior_alpha_clamp_max);
	case MK_FINAL_MEAN_ENTROPY: return parse_double_value(value, &info->final_mean_entropy);
	case MK_FINAL_MEAN_PU: return parse_double_value(value, &info->final_mean_pU);
	case MK_FINAL_SUM_WEDGE_K: return parse_double_value(value, &info->final_sum_wedge_k);
	case MK_POSTERIOR_REFRESH_TEMPERATURE: return parse_double_value(value, &info->posterior_refresh_temperature);
	case MK_POSTERIOR_REFRESH_MEAN_KL: return parse_double_value(value, &info->posterior_refresh_mean_kl);
	case MK_POSTERIOR_REFRESH_TOP_STATE_SWITCH_FRAC: return parse_double_value(value, &info->posterior_refresh_top_state_switch_frac);
	case MK_POSTERIOR_REFRESH_MEAN_PU_BEFORE: return parse_double_value(value, &info->posterior_refresh_mean_pU_before);
	case MK_POSTERIOR_REFRESH_MEAN_PU_AFTER: return parse_double_value(value, &info->posterior_refresh_mean_pU_after);
	case MK_OUTPUT_BPAIR_POSTERIOR: snprintf(info->output_bpair_posterior, sizeof(info->output_bpair_posterior), "%s", value); return 0;
	case MK_OUTPUT_COORDS: snprintf(info->output_coords, sizeof(info->output_coords), "%s", value); return 0;
	case MK_OUTPUT_LOOP_DIAG: snprintf(info->output_loop_diag, sizeof(info->output_loop_diag), "%s", value); return 0;
	default:
		if (strcmp(key, "output_raw_posterior") == 0) {
			info->has_output_raw_posterior = 1;
			snprintf(info->output_raw_posterior, sizeof(info->output_raw_posterior), "%s", value);
			return 0;
		}
		if (manifest_optional_current_key(key))
			return 0;
		return -1;
	}
}

static int audit_manifest(const char *path, struct manifest_info *info, struct file_audit *audit)
{
	FILE *fp = fopen(path, "r");
	char line[HK_AUDIT_LINE_MAX];
	int failed = 0;
	int64_t line_no = 0;
	int i;
	const char *expected_sample;

	memset(info, 0, sizeof(*info));
	info->rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	info->rho_train = NAN;
	info->init_eps = NAN;
	info->init_noise_scale = NAN;
	info->init_scale = NAN;
	info->chr_sep_unit = NAN;
	info->lambda_chr_sep = NAN;
	info->lambda_copytrack = 0.0;
	info->lambda_global_copytrack = 0.0;
	info->lambda_normdir_copytrack = 0.0;
	info->normdir_copytrack_eps_unit = 0.0;
	info->copytrack_uses_phase_labels = 0;
	info->copytrack_uses_charm_or_reference = 0;
	info->normdir_copytrack_uses_phase_labels = 0;
	info->normdir_copytrack_uses_charm_or_reference = 0;
	snprintf(info->trans_gate_mode, sizeof(info->trans_gate_mode), "off");
	snprintf(info->trans_gate_scope, sizeof(info->trans_gate_scope), "trans_raw_contacts_by_bpair_posterior");
	snprintf(info->trans_gate_application_point, sizeof(info->trans_gate_application_point), "post_estep_pre_softall_mstep");
	snprintf(info->trans_gate_source, sizeof(info->trans_gate_source), "blind_current_posterior_confidence");
	info->trans_gate_min_pmax = 1.0;
	info->trans_gate_min_margin = 1.0;
	info->trans_gate_min_neg_entropy = -log((double)HK_BLIND_N_STATE);
	info->trans_gate_uses_phase_labels = 0;
	info->trans_gate_uses_charm_or_reference = 0;
	audit_init(audit);
	if (fp == 0) {
		add_example(audit, "manifest open failed");
		audit->bad_rows++;
		return 1;
	}
	while (fgets(line, sizeof(line), fp)) {
		char *tab;
		++line_no;
		scan_forbidden_line(audit, "manifest", line_no, line);
		trim_line(line);
		if (line_no == 1) {
			if (strcmp(line, "key\tvalue") != 0) {
				add_example(audit, "manifest header is not key/value");
				++audit->bad_rows;
				failed = 1;
			}
			continue;
		}
		tab = strchr(line, '\t');
		if (tab == 0) {
			add_example_fmt(audit, "manifest", line_no, "missing tab", NAN);
			++audit->bad_rows;
			failed = 1;
			continue;
		}
		*tab++ = 0;
		if (manifest_set_value(info, line, tab) != 0) {
			char buf[256];
			snprintf(buf, sizeof(buf), "manifest row %lld bad numeric value for %s",
					 (long long)line_no, line);
			add_example(audit, buf);
			++audit->bad_rows;
			failed = 1;
		}
		++audit->rows;
	}
	fclose(fp);
	for (i = 0; i < MK_N_REQUIRED; ++i) {
		if (!info->seen[i]) {
			char buf[256];
			snprintf(buf, sizeof(buf), "manifest missing required key index %d", i);
			add_example(audit, buf);
			failed = 1;
		}
	}
	expected_sample = getenv("HK_BLIND_SAMPLE");
	if (expected_sample == 0 || expected_sample[0] == 0)
		expected_sample = "P9016";
	if (strcmp(info->sample, expected_sample) != 0) {
		add_example(audit, "manifest sample differs from HK_BLIND_SAMPLE/P9016 expectation");
		failed = 1;
	}
	if (info->runner_family[0] == 0 || info->runner_version[0] == 0 ||
		info->default_profile[0] == 0 || info->baseline[0] == 0) {
		add_example(audit, "manifest runner identity fields are empty");
		failed = 1;
	}
	if (strcmp(info->runner_family, "p9016_minimal") != 0 &&
		strcmp(info->runner_family, "test_fixture") != 0) {
		add_example(audit, "manifest runner_family is not recognized");
		failed = 1;
	}
	if (strcmp(info->baseline, "softall") != 0) {
		add_example(audit, "manifest baseline is not softall");
		failed = 1;
	}
	if (strcmp(info->mstep_graph_mode, "raw_expected_soft_all") != 0) {
		add_example(audit, "manifest mstep_graph_mode is not raw_expected_soft_all");
		failed = 1;
	}
	if (strcmp(info->status, "OK") != 0) {
		add_example(audit, "manifest status is not OK");
		failed = 1;
	}
	if (info->n_raw <= 0 || info->n_bpair <= 0 || info->n_beads <= 0) {
		add_example(audit, "manifest count is nonpositive");
		failed = 1;
	}
	if (info->resolution <= 0) {
		add_example(audit, "manifest resolution is nonpositive");
		failed = 1;
	}
	if (info->bin_size_bp <= 0) {
		add_example(audit, "manifest bin_size_bp is nonpositive");
		failed = 1;
	}
	if (info->resolution != info->bin_size_bp) {
		add_example(audit, "manifest resolution disagrees with bin_size_bp");
		failed = 1;
	}
	if (info->n_iter <= 0 || info->relax_steps < 0 || info->init_seed < 0) {
		add_example(audit, "manifest iteration/init integer value out of range");
		failed = 1;
	}
	if (!isfinite(info->unit) || info->unit <= 0.0 ||
		!isfinite(info->d_scale) || info->d_scale <= 0.0 ||
		!isfinite(info->base_k_effective) || info->base_k_effective < 0.0 ||
		!isfinite(info->base_k_min) || !isfinite(info->base_k_mean) ||
		!isfinite(info->base_k_max) || info->base_k_min < 0.0 ||
		info->base_k_mean < 0.0 || info->base_k_max < 0.0 ||
		info->base_k_n_nonfinite != 0) {
		add_example(audit, "manifest unit/d_scale/base_k_effective out of range");
		failed = 1;
	}
	if (strcmp(info->base_k_mode, "uniform") == 0 &&
		(!check_close(info->base_k_effective, info->base_k_min) ||
		 !check_close(info->base_k_effective, info->base_k_mean) ||
		 !check_close(info->base_k_effective, info->base_k_max))) {
		add_example(audit, "manifest uniform base_k stats disagree");
		failed = 1;
	}
	if (strcmp(info->base_k_mode, "uniform") != 0 &&
		strcmp(info->base_k_mode, "neighbor_median") != 0) {
		add_example(audit, "manifest base_k_mode is not recognized");
		failed = 1;
	}
	if (!recognized_init_mode(info->init_mode)) {
		add_example(audit, "manifest init_mode is not recognized");
		failed = 1;
	}
	if (strcmp(info->prior_mode, "uniform") != 0) {
		add_example(audit, "manifest prior_mode is not recognized");
		failed = 1;
	}
	if (!isfinite(info->prior_eps) || info->prior_eps <= 0.0 || info->prior_eps > 0.5 ||
		!isfinite(info->prior_inter_density) || info->prior_inter_density < 0.0 ||
		!isfinite(info->prior_observed_inter) || info->prior_observed_inter < 0.0 ||
		!isfinite(info->prior_possible_inter) || info->prior_possible_inter < 0.0 ||
		info->prior_n_distance_bins < 0 || info->prior_n_alpha < 0 ||
		!isfinite(info->prior_alpha_min) || !isfinite(info->prior_alpha_median) ||
		!isfinite(info->prior_alpha_max) ||
		!isfinite(info->prior_alpha_clamp_min) || info->prior_alpha_clamp_min <= 0.0 ||
		info->prior_alpha_clamp_min > 0.5 ||
		!isfinite(info->prior_alpha_clamp_max) ||
		!check_close(info->prior_alpha_clamp_max, 0.5)) {
		add_example(audit, "manifest prior metadata out of range");
		failed = 1;
	}
	if (info->prior_n_alpha != 0 || strcmp(info->prior_smoothing_method, "none") != 0) {
		add_example(audit, "manifest uniform prior smoothing method is not none");
		failed = 1;
	}
	if (!recognized_rho_train_mode(info->rho_train_mode)) {
		add_example(audit, "manifest rho_train_mode is not recognized");
		failed = 1;
	}
	if (info->seen[MK_RHO_TRAIN] &&
		(!isfinite(info->rho_train) || info->rho_train < 0.0)) {
		add_example(audit, "manifest rho_train is out of range");
		failed = 1;
	}
	if (strcmp(info->d_scale_mode, "raw_count") != 0 &&
		strcmp(info->d_scale_mode, "posterior_count") != 0 &&
		strcmp(info->d_scale_mode, "expected_count") != 0 &&
		strcmp(info->d_scale_mode, "tempered_posterior_count") != 0) {
		add_example(audit, "manifest d_scale_mode is not recognized");
		failed = 1;
	}
	if (info->same_bin_filter_enabled != 1 || info->n_raw_same_bin_excluded < 0 ||
		info->n_bpair_same_bin_excluded < 0) {
		add_example(audit, "manifest same-bin filter fields out of range");
		failed = 1;
	}
	if (info->seen[MK_HELDOUT_ENABLED] || info->seen[MK_HELDOUT_FRACTION] ||
		info->seen[MK_HELDOUT_N_RAW_TRAIN] || info->seen[MK_HELDOUT_N_RAW_HELDOUT]) {
		if (!(info->seen[MK_HELDOUT_ENABLED] && info->seen[MK_HELDOUT_FRACTION] &&
			  info->seen[MK_HELDOUT_SEED] && info->seen[MK_HELDOUT_N_INPUT_RAW_TOTAL] &&
			  info->seen[MK_HELDOUT_N_RAW_TRAIN] && info->seen[MK_HELDOUT_N_RAW_HELDOUT] &&
			  info->seen[MK_HELDOUT_N_BPAIR_HELDOUT] && info->seen[MK_HELDOUT_N_BPAIR_EVAL] &&
			  info->seen[MK_HELDOUT_N_RAW_EVAL] &&
			  info->seen[MK_HELDOUT_MEAN_EXPECTED_ENERGY] &&
			  info->seen[MK_HELDOUT_MEAN_MIN_ENERGY] &&
			  info->seen[MK_HELDOUT_MEAN_ENTROPY] &&
			  info->seen[MK_HELDOUT_MEAN_PU] &&
			  info->seen[MK_HELDOUT_MEAN_BEST_NORMALIZED_DISTANCE] &&
			  info->seen[MK_HELDOUT_SHORT_DISTANCE_FRAC])) {
			add_example(audit, "manifest heldout fields are incomplete");
			failed = 1;
		} else if ((info->heldout_enabled != 0 && info->heldout_enabled != 1) ||
				   !isfinite(info->heldout_fraction) || info->heldout_fraction < 0.0 ||
				   info->heldout_fraction >= 1.0 || info->heldout_seed < 0 ||
				   info->heldout_n_input_raw_total < 0 ||
				   info->heldout_n_raw_train < 0 ||
				   info->heldout_n_raw_heldout < 0 ||
				   info->heldout_n_bpair_heldout < 0 ||
				   info->heldout_n_bpair_eval < 0 ||
				   info->heldout_n_raw_eval < 0 ||
				   info->heldout_n_raw_train + info->heldout_n_raw_heldout !=
				   info->heldout_n_input_raw_total ||
				   info->heldout_n_raw_train != info->n_raw ||
				   info->heldout_n_bpair_eval > info->heldout_n_bpair_heldout ||
				   info->heldout_n_raw_eval > info->heldout_n_raw_heldout ||
				   !isfinite(info->heldout_mean_expected_energy) ||
				   !isfinite(info->heldout_mean_min_energy) ||
				   !isfinite(info->heldout_mean_entropy) ||
				   !isfinite(info->heldout_mean_pU) ||
				   !isfinite(info->heldout_mean_best_normalized_distance) ||
				   !isfinite(info->heldout_short_distance_frac) ||
				   info->heldout_mean_entropy < 0.0 ||
				   info->heldout_mean_pU < 0.0 || info->heldout_mean_pU > 1.0 ||
				   info->heldout_short_distance_frac < 0.0 ||
				   info->heldout_short_distance_frac > 1.0) {
			add_example(audit, "manifest heldout fields are out of range");
			failed = 1;
		} else if (info->heldout_enabled == 0 &&
				   (info->heldout_n_raw_heldout != 0 ||
					info->heldout_n_bpair_eval != 0 ||
					info->heldout_n_raw_eval != 0 ||
					info->heldout_fraction != 0.0)) {
			add_example(audit, "manifest heldout disabled fields are inconsistent");
			failed = 1;
		} else if (info->heldout_enabled == 1 &&
				   (info->heldout_fraction <= 0.0 ||
					info->heldout_n_raw_heldout <= 0 ||
					info->heldout_n_bpair_heldout <= 0)) {
			add_example(audit, "manifest heldout enabled fields are inconsistent");
			failed = 1;
		}
	}
	if (!info->seen[MK_USES_PHASE_LABELS]) {
		add_example(audit, "manifest missing uses_phase_labels");
		failed = 1;
	} else if (info->uses_phase_labels != 0) {
		add_example(audit, "manifest top-level uses_phase_labels must be 0");
		failed = 1;
	}
	if (!(info->seen[MK_INPUT_CONTACT_SOURCE] &&
		  info->seen[MK_USES_CHARM_OR_REFERENCE] &&
		  info->seen[MK_USES_CHARM_FOR_TRAINING] &&
		  info->seen[MK_REFERENCE_DERIVED_POSITIVE_CONTROL])) {
		add_example(audit, "manifest missing top-level training-boundary fields");
		failed = 1;
		} else {
			const int allow_reference = env_flag_enabled("HK_BLIND_AUDIT_ALLOW_REFERENCE_DERIVED_TRAINING");
			const int allow_nonstandard = env_flag_enabled("HK_BLIND_AUDIT_ALLOW_NONSTANDARD_PAIRS");
			if (info->uses_charm_or_reference != 0 ||
				info->uses_charm_for_training != 0 ||
				info->reference_derived_positive_control != 0) {
			if (!allow_reference) {
				add_example(audit, "manifest top-level reference/CHARM training flags are nonzero");
				failed = 1;
				}
			} else if (strcmp(info->input_contact_source, "raw_pairs") != 0) {
				if (!(allow_nonstandard &&
					  strcmp(info->input_contact_source, "contacts_seg_derived_pairs") == 0)) {
					add_example(audit, "manifest input_contact_source is not raw_pairs for blind training");
					failed = 1;
				}
			}
		}
	if (!info->seen[MK_APPROVED_P9016_RAW_PAIRS_REALPATH]) {
		add_example(audit, "manifest missing approved_p9016_raw_pairs_realpath");
		failed = 1;
		} else if (info->approved_p9016_raw_pairs_realpath != 1) {
			const int allow_custom = env_flag_enabled("HK_BLIND_AUDIT_ALLOW_CUSTOM_RAW_PAIRS");
			const int allow_nonstandard = env_flag_enabled("HK_BLIND_AUDIT_ALLOW_NONSTANDARD_PAIRS");
			const int allow_reference = env_flag_enabled("HK_BLIND_AUDIT_ALLOW_REFERENCE_DERIVED_TRAINING");
			if (!allow_custom && !allow_nonstandard && !allow_reference) {
				add_example(audit, "manifest approved_p9016_raw_pairs_realpath is not 1");
				failed = 1;
			}
	}
	if (strcmp(info->raw_posterior_same_bin_policy, "uniform_unknown_rows") != 0) {
		add_example(audit, "manifest raw posterior same-bin policy is not uniform_unknown_rows");
		failed = 1;
	}
	if (!isfinite(info->min_sep_unit) || info->min_sep_unit < 0.0 ||
		!isfinite(info->lambda_sep) || info->lambda_sep < 0.0 ||
		!isfinite(info->lambda_copytrack) || info->lambda_copytrack < 0.0 ||
		!isfinite(info->lambda_global_copytrack) || info->lambda_global_copytrack < 0.0 ||
		!isfinite(info->lambda_normdir_copytrack) || info->lambda_normdir_copytrack < 0.0 ||
		!isfinite(info->normdir_copytrack_eps_unit) || info->normdir_copytrack_eps_unit < 0.0 ||
		!isfinite(info->relax_step) || info->relax_step < 0.0) {
		add_example(audit, "manifest separation/relax value out of range");
		failed = 1;
	}
	if (info->lambda_normdir_copytrack > 0.0 && info->normdir_copytrack_eps_unit <= 0.0) {
		add_example(audit, "normdir copytrack eps must be positive when normdir force is enabled");
		failed = 1;
	}
	if (info->copytrack_uses_phase_labels != 0 ||
		info->copytrack_uses_charm_or_reference != 0) {
		add_example(audit, "copytrack manifest is not blind-safe");
		failed = 1;
	}
	if (info->global_copytrack_uses_phase_labels != 0 ||
		info->global_copytrack_uses_charm_or_reference != 0) {
		add_example(audit, "global copytrack manifest is not blind-safe");
		failed = 1;
	}
	if (info->normdir_copytrack_uses_phase_labels != 0 ||
		info->normdir_copytrack_uses_charm_or_reference != 0) {
		add_example(audit, "normdir copytrack manifest is not blind-safe");
		failed = 1;
	}
	if (!isfinite(info->temperature_start) || info->temperature_start <= 0.0 ||
		!isfinite(info->temperature_end) || info->temperature_end <= 0.0 ||
		!isfinite(info->rho_train_start) || info->rho_train_start < 0.0 ||
		!isfinite(info->rho_train_end) || info->rho_train_end < 0.0) {
		add_example(audit, "manifest schedule value out of range");
		failed = 1;
	}
	if (!isfinite(info->rho_train_floor) || info->rho_train_floor < 0.0 ||
		info->rho_train_floor > 1.0) {
		add_example(audit, "manifest rho floor out of range");
		failed = 1;
	}
	if (info->seen[MK_N_RAW_CIS] || info->seen[MK_N_RAW_TRANS] ||
		info->seen[MK_N_BPAIR_CIS] || info->seen[MK_N_BPAIR_TRANS]) {
		if (!(info->seen[MK_N_RAW_CIS] && info->seen[MK_N_RAW_TRANS] &&
			  info->seen[MK_N_BPAIR_CIS] && info->seen[MK_N_BPAIR_TRANS]) ||
			info->n_raw_cis < 0 || info->n_raw_trans < 0 ||
			info->n_bpair_cis < 0 || info->n_bpair_trans < 0 ||
			info->n_raw_cis + info->n_raw_trans != info->n_raw ||
			info->n_bpair_cis + info->n_bpair_trans != info->n_bpair) {
			add_example(audit, "manifest cis/trans count fields are inconsistent");
			failed = 1;
		}
	}
	if (info->n_iter <= 1 &&
		(!check_close(info->temperature_start, info->temperature_end) ||
		 !check_close(info->rho_train_start, info->rho_train_end))) {
		add_example(audit, "manifest single-iteration schedule endpoints differ from effective value");
		failed = 1;
	}
	if (!isfinite(info->init_eps_effective) || info->init_eps_effective < 0.0 ||
		!isfinite(info->init_noise_scale_effective) || info->init_noise_scale_effective < 0.0 ||
		!isfinite(info->init_scale_effective) || info->init_scale_effective < 0.0) {
		add_example(audit, "manifest init value out of range");
		failed = 1;
	}
	if (info->seen[MK_INIT_EPS] &&
		(!isfinite(info->init_eps) || info->init_eps < 0.0)) {
		add_example(audit, "manifest init_eps out of range");
		failed = 1;
	}
	if (info->seen[MK_INIT_NOISE_SCALE] &&
		(!isfinite(info->init_noise_scale) || info->init_noise_scale < 0.0)) {
		add_example(audit, "manifest init_noise_scale out of range");
		failed = 1;
	}
	if (info->seen[MK_INIT_SCALE] &&
		(!isfinite(info->init_scale) || info->init_scale < 0.0)) {
		add_example(audit, "manifest init_scale out of range");
		failed = 1;
	}
	if (info->seen[MK_CHR_SEP_UNIT] &&
		(!isfinite(info->chr_sep_unit) || !check_close(info->chr_sep_unit, 0.0))) {
		add_example(audit, "manifest chr_sep_unit must remain zero");
		failed = 1;
	}
	if (info->seen[MK_LAMBDA_CHR_SEP] &&
		(!isfinite(info->lambda_chr_sep) || !check_close(info->lambda_chr_sep, 0.0))) {
		add_example(audit, "manifest lambda_chr_sep must remain zero");
		failed = 1;
	}
	if (info->seen[MK_TRANS_CHR_PAIR_PRIOR_LAMBDA] ||
		info->seen[MK_TRANS_CHR_PAIR_PRIOR_EPS] ||
		info->seen[MK_TRANS_CHR_PAIR_PRIOR_POWER] ||
		info->seen[MK_TRANS_CHR_PAIR_PRIOR_WARMUP_ITER] ||
		info->seen[MK_TRANS_CHR_PAIR_PRIOR_MODE]) {
		if (!(info->seen[MK_TRANS_CHR_PAIR_PRIOR_LAMBDA] &&
			  info->seen[MK_TRANS_CHR_PAIR_PRIOR_EPS] &&
			  info->seen[MK_TRANS_CHR_PAIR_PRIOR_POWER] &&
			  info->seen[MK_TRANS_CHR_PAIR_PRIOR_WARMUP_ITER] &&
			  info->seen[MK_TRANS_CHR_PAIR_PRIOR_MODE]) ||
			!isfinite(info->trans_chr_pair_prior_lambda) ||
			info->trans_chr_pair_prior_lambda < 0.0 ||
			info->trans_chr_pair_prior_lambda > 1.0 ||
			!isfinite(info->trans_chr_pair_prior_eps) ||
			info->trans_chr_pair_prior_eps < 0.0 ||
			!isfinite(info->trans_chr_pair_prior_power) ||
			info->trans_chr_pair_prior_power < 0.0 ||
			info->trans_chr_pair_prior_warmup_iter < 0) {
			add_example(audit, "manifest trans chromosome-pair prior metadata out of range");
			failed = 1;
		}
		if (info->trans_chr_pair_prior_lambda > 0.0) {
			if (strcmp(info->trans_chr_pair_prior_mode, "blind_posterior_chrom_pair") != 0) {
				add_example(audit, "manifest trans chromosome-pair prior mode is inconsistent");
				failed = 1;
			}
		} else if (strcmp(info->trans_chr_pair_prior_mode, "off") != 0) {
			add_example(audit, "manifest trans chromosome-pair prior mode should be off");
			failed = 1;
		}
		if (info->seen[MK_TRANS_CHR_PAIR_PRIOR_USES_PHASE_LABELS] &&
			info->trans_chr_pair_prior_uses_phase_labels != 0) {
			add_example(audit, "manifest trans chromosome-pair prior uses phase labels");
			failed = 1;
		}
		if (info->seen[MK_TRANS_CHR_PAIR_PRIOR_USES_CHARM_OR_REFERENCE] &&
			info->trans_chr_pair_prior_uses_charm_or_reference != 0) {
			add_example(audit, "manifest trans chromosome-pair prior uses CHARM/reference");
			failed = 1;
		}
	}
	if (info->seen[MK_TRANS_CHR_PAIR_MSTEP_LAMBDA] ||
		info->seen[MK_TRANS_CHR_PAIR_MSTEP_EPS] ||
		info->seen[MK_TRANS_CHR_PAIR_MSTEP_POWER] ||
		info->seen[MK_TRANS_CHR_PAIR_MSTEP_WARMUP_ITER] ||
		info->seen[MK_TRANS_CHR_PAIR_MSTEP_MODE]) {
		if (!(info->seen[MK_TRANS_CHR_PAIR_MSTEP_LAMBDA] &&
			  info->seen[MK_TRANS_CHR_PAIR_MSTEP_EPS] &&
				  info->seen[MK_TRANS_CHR_PAIR_MSTEP_POWER] &&
				  info->seen[MK_TRANS_CHR_PAIR_MSTEP_WARMUP_ITER] &&
				  info->seen[MK_TRANS_CHR_PAIR_MSTEP_MODE] &&
				  info->seen[MK_TRANS_CHR_PAIR_MSTEP_SCOPE] &&
				  info->seen[MK_TRANS_CHR_PAIR_MSTEP_APPLICATION_POINT] &&
				  info->seen[MK_TRANS_CHR_PAIR_MSTEP_USES_PHASE_LABELS] &&
				  info->seen[MK_TRANS_CHR_PAIR_MSTEP_USES_CHARM_OR_REFERENCE]) ||
			!isfinite(info->trans_chr_pair_mstep_lambda) ||
			info->trans_chr_pair_mstep_lambda < 0.0 ||
			info->trans_chr_pair_mstep_lambda > 1.0 ||
			!isfinite(info->trans_chr_pair_mstep_eps) ||
			info->trans_chr_pair_mstep_eps < 0.0 ||
			!isfinite(info->trans_chr_pair_mstep_power) ||
			info->trans_chr_pair_mstep_power < 0.0 ||
			info->trans_chr_pair_mstep_warmup_iter < 0) {
			add_example(audit, "manifest trans chromosome-pair M-step metadata out of range");
			failed = 1;
		}
		if (info->trans_chr_pair_mstep_lambda > 0.0) {
			if (strcmp(info->trans_chr_pair_mstep_mode, "blind_posterior_chrom_pair_mstep") != 0) {
				add_example(audit, "manifest trans chromosome-pair M-step mode is inconsistent");
				failed = 1;
			}
		} else if (strcmp(info->trans_chr_pair_mstep_mode, "off") != 0) {
			add_example(audit, "manifest trans chromosome-pair M-step mode should be off");
			failed = 1;
		}
		if (strcmp(info->trans_chr_pair_mstep_scope, "trans_edges_only") != 0) {
			add_example(audit, "manifest trans chromosome-pair M-step scope is wrong");
			failed = 1;
		}
		if (strcmp(info->trans_chr_pair_mstep_application_point, "post_estep_pre_softall_mstep") != 0) {
			add_example(audit, "manifest trans chromosome-pair M-step application point is wrong");
			failed = 1;
		}
			if (info->trans_chr_pair_mstep_uses_phase_labels != 0) {
				add_example(audit, "manifest trans chromosome-pair M-step uses phase labels");
				failed = 1;
			}
			if (info->trans_chr_pair_mstep_uses_charm_or_reference != 0) {
				add_example(audit, "manifest trans chromosome-pair M-step uses CHARM/reference");
				failed = 1;
			}
	}
		if (info->seen[MK_TRANS_CONTACT_SCALING_MODE] ||
			info->seen[MK_TRANS_K_MULTIPLIER] ||
			info->seen[MK_TRANS_DSCALE_MULTIPLIER] ||
			info->seen[MK_TRANS_D_SCALE_POSTERIOR_GAMMA] ||
			info->seen[MK_TRANS_CONTACT_SCALING_SCOPE]) {
			if (!(info->seen[MK_TRANS_CONTACT_SCALING_MODE] &&
				  info->seen[MK_TRANS_K_MULTIPLIER] &&
				  info->seen[MK_TRANS_DSCALE_MULTIPLIER] &&
				  info->seen[MK_D_SCALE_POSTERIOR_GAMMA] &&
				  info->seen[MK_TRANS_D_SCALE_POSTERIOR_GAMMA] &&
				  info->seen[MK_TRANS_CONTACT_SCALING_SCOPE]) ||
				!isfinite(info->trans_k_multiplier) ||
				info->trans_k_multiplier < 0.0 ||
				!isfinite(info->trans_dscale_multiplier) ||
				info->trans_dscale_multiplier <= 0.0 ||
				!isfinite(info->trans_d_scale_posterior_gamma) ||
				info->trans_d_scale_posterior_gamma < 0.0) {
			add_example(audit, "manifest trans contact scaling metadata out of range");
			failed = 1;
		}
		if (strcmp(info->trans_contact_scaling_scope, "trans_edges_only") != 0) {
			add_example(audit, "manifest trans contact scaling scope is not trans_edges_only");
			failed = 1;
		}
			if (check_close(info->trans_k_multiplier, 1.0) &&
				check_close(info->trans_dscale_multiplier, 1.0) &&
				check_close(info->trans_d_scale_posterior_gamma,
							info->d_scale_posterior_gamma)) {
			if (strcmp(info->trans_contact_scaling_mode, "off") != 0) {
				add_example(audit, "manifest trans contact scaling mode should be off");
				failed = 1;
			}
		} else if (strcmp(info->trans_contact_scaling_mode, "trans_edge_scaling") != 0) {
			add_example(audit, "manifest trans contact scaling mode is inconsistent");
			failed = 1;
		}
		if (info->seen[MK_TRANS_CONTACT_SCALING_USES_PHASE_LABELS] &&
			info->trans_contact_scaling_uses_phase_labels != 0) {
			add_example(audit, "manifest trans contact scaling uses phase labels");
			failed = 1;
		}
		if (info->seen[MK_TRANS_CONTACT_SCALING_USES_CHARM_OR_REFERENCE] &&
			info->trans_contact_scaling_uses_charm_or_reference != 0) {
			add_example(audit, "manifest trans contact scaling uses CHARM/reference");
			failed = 1;
		}
	}
		if (info->seen[MK_TRANS_TOP1_MSTEP_MODE] ||
			info->seen[MK_TRANS_TOP1_MSTEP_MIN_PMAX] ||
			info->seen[MK_TRANS_TOP1_MSTEP_MIN_MARGIN] ||
		info->seen[MK_TRANS_TOP1_MSTEP_MIX_WEIGHT] ||
		info->seen[MK_TRANS_TOP1_MSTEP_SCOPE] ||
		info->seen[MK_TRANS_TOP1_MSTEP_APPLICATION_POINT]) {
		if (!(info->seen[MK_TRANS_TOP1_MSTEP_MODE] &&
			  info->seen[MK_TRANS_TOP1_MSTEP_MIN_PMAX] &&
			  info->seen[MK_TRANS_TOP1_MSTEP_MIN_MARGIN] &&
			  info->seen[MK_TRANS_TOP1_MSTEP_MIX_WEIGHT] &&
			  info->seen[MK_TRANS_TOP1_MSTEP_SCOPE] &&
			  info->seen[MK_TRANS_TOP1_MSTEP_APPLICATION_POINT]) ||
			!isfinite(info->trans_top1_mstep_min_pmax) ||
			info->trans_top1_mstep_min_pmax < 0.0 ||
			info->trans_top1_mstep_min_pmax > 1.0 ||
			!isfinite(info->trans_top1_mstep_min_margin) ||
			info->trans_top1_mstep_min_margin < 0.0 ||
			info->trans_top1_mstep_min_margin > 1.0 ||
			!isfinite(info->trans_top1_mstep_mix_weight) ||
			info->trans_top1_mstep_mix_weight < 0.0 ||
			info->trans_top1_mstep_mix_weight > 1.0) {
			add_example(audit, "manifest trans top1 M-step metadata out of range");
			failed = 1;
		}
		if (strcmp(info->trans_top1_mstep_mode, "off") != 0 &&
			strcmp(info->trans_top1_mstep_mode, "hard") != 0 &&
			strcmp(info->trans_top1_mstep_mode, "mix") != 0) {
			add_example(audit, "manifest trans top1 M-step mode is invalid");
			failed = 1;
		}
		if (strcmp(info->trans_top1_mstep_scope, "trans_edges_only") != 0) {
			add_example(audit, "manifest trans top1 M-step scope is not trans_edges_only");
			failed = 1;
		}
		if (strcmp(info->trans_top1_mstep_application_point, "post_estep_pre_softall_mstep") != 0) {
			add_example(audit, "manifest trans top1 M-step application point is wrong");
			failed = 1;
		}
		if (strcmp(info->trans_top1_mstep_mode, "off") == 0 &&
			(!check_close(info->trans_top1_mstep_min_pmax, 1.0) ||
			 !check_close(info->trans_top1_mstep_min_margin, 1.0) ||
			 !check_close(info->trans_top1_mstep_mix_weight, 0.0))) {
			add_example(audit, "manifest trans top1 M-step off defaults are inconsistent");
			failed = 1;
		}
		if (strcmp(info->trans_top1_mstep_mode, "hard") == 0 &&
			!check_close(info->trans_top1_mstep_mix_weight, 1.0)) {
			add_example(audit, "manifest trans top1 hard mode should have mix weight 1");
			failed = 1;
		}
		if (info->seen[MK_TRANS_TOP1_MSTEP_USES_PHASE_LABELS] &&
			info->trans_top1_mstep_uses_phase_labels != 0) {
			add_example(audit, "manifest trans top1 M-step uses phase labels");
			failed = 1;
		}
			if (info->seen[MK_TRANS_TOP1_MSTEP_USES_CHARM_OR_REFERENCE] &&
				info->trans_top1_mstep_uses_charm_or_reference != 0) {
				add_example(audit, "manifest trans top1 M-step uses CHARM/reference");
					failed = 1;
				}
			}
		if (info->seen[MK_TRANS_CALLABLE_ANCHOR_MODE] ||
			info->seen[MK_TRANS_CALLABLE_ANCHOR_TOP_FRAC] ||
			info->seen[MK_TRANS_CALLABLE_ANCHOR_MIX_WEIGHT] ||
			info->seen[MK_TRANS_CALLABLE_ANCHOR_MIN_N_RAW] ||
			info->seen[MK_TRANS_CALLABLE_ANCHOR_WARMUP_ITER]) {
			if (!(info->seen[MK_TRANS_CALLABLE_ANCHOR_MODE] &&
					  info->seen[MK_TRANS_CALLABLE_ANCHOR_TOP_FRAC] &&
					  info->seen[MK_TRANS_CALLABLE_ANCHOR_MIX_WEIGHT] &&
					  info->seen[MK_TRANS_CALLABLE_ANCHOR_MIN_N_RAW] &&
					  info->seen[MK_TRANS_CALLABLE_ANCHOR_WARMUP_ITER] &&
					  info->seen[MK_TRANS_CALLABLE_ANCHOR_USES_PHASE_LABELS] &&
					  info->seen[MK_TRANS_CALLABLE_ANCHOR_USES_CHARM_OR_REFERENCE]) ||
				!isfinite(info->trans_callable_anchor_top_frac) ||
				info->trans_callable_anchor_top_frac < 0.0 ||
				info->trans_callable_anchor_top_frac > 1.0 ||
				!isfinite(info->trans_callable_anchor_mix_weight) ||
				info->trans_callable_anchor_mix_weight < 0.0 ||
				info->trans_callable_anchor_mix_weight > 1.0 ||
				info->trans_callable_anchor_min_n_raw < 0 ||
				info->trans_callable_anchor_warmup_iter < 0) {
				add_example(audit, "manifest trans callable-anchor metadata out of range");
				failed = 1;
			}
			if (strcmp(info->trans_callable_anchor_mode, "off") != 0 &&
				strcmp(info->trans_callable_anchor_mode, "neg_entropy") != 0 &&
				strcmp(info->trans_callable_anchor_mode, "pmax") != 0 &&
				strcmp(info->trans_callable_anchor_mode, "margin") != 0) {
				add_example(audit, "manifest trans callable-anchor mode is invalid");
				failed = 1;
			}
				if (info->trans_callable_anchor_uses_phase_labels != 0) {
					add_example(audit, "manifest trans callable-anchor uses phase labels");
					failed = 1;
				}
				if (info->trans_callable_anchor_uses_charm_or_reference != 0) {
					add_example(audit, "manifest trans callable-anchor uses CHARM/reference");
					failed = 1;
				}
		}
			if (info->seen[MK_TRANS_GATE_MODE] ||
				info->seen[MK_TRANS_GATE_MIN_PMAX] ||
			info->seen[MK_TRANS_GATE_MIN_MARGIN] ||
			info->seen[MK_TRANS_GATE_MIN_NEG_ENTROPY] ||
			info->seen[MK_TRANS_GATE_SCOPE] ||
			info->seen[MK_TRANS_GATE_APPLICATION_POINT] ||
			info->seen[MK_TRANS_GATE_SOURCE]) {
			if (!(info->seen[MK_TRANS_GATE_MODE] &&
				  info->seen[MK_TRANS_GATE_MIN_PMAX] &&
				  info->seen[MK_TRANS_GATE_MIN_MARGIN] &&
				  info->seen[MK_TRANS_GATE_MIN_NEG_ENTROPY] &&
				  info->seen[MK_TRANS_GATE_SCOPE] &&
				  info->seen[MK_TRANS_GATE_APPLICATION_POINT] &&
				  info->seen[MK_TRANS_GATE_SOURCE]) ||
				!isfinite(info->trans_gate_min_pmax) ||
				info->trans_gate_min_pmax < 0.0 ||
				info->trans_gate_min_pmax > 1.0 ||
				!isfinite(info->trans_gate_min_margin) ||
				info->trans_gate_min_margin < 0.0 ||
				info->trans_gate_min_margin > 1.0 ||
				!isfinite(info->trans_gate_min_neg_entropy) ||
				info->trans_gate_min_neg_entropy < -log((double)HK_BLIND_N_STATE) - HK_AUDIT_TOL ||
				info->trans_gate_min_neg_entropy > 0.0) {
				add_example(audit, "manifest trans gate metadata out of range");
				failed = 1;
			}
			if (strcmp(info->trans_gate_mode, "off") != 0 &&
				strcmp(info->trans_gate_mode, "pmax") != 0 &&
				strcmp(info->trans_gate_mode, "margin") != 0 &&
				strcmp(info->trans_gate_mode, "neg_entropy") != 0 &&
				strcmp(info->trans_gate_mode, "pmax_margin") != 0) {
				add_example(audit, "manifest trans gate mode is invalid");
				failed = 1;
			}
		if (strcmp(info->trans_gate_scope, "trans_raw_contacts_by_bpair_posterior") != 0) {
			add_example(audit, "manifest trans gate scope is not trans_raw_contacts_by_bpair_posterior");
			failed = 1;
		}
			if (strcmp(info->trans_gate_application_point, "post_estep_pre_softall_mstep") != 0) {
				add_example(audit, "manifest trans gate application point is wrong");
				failed = 1;
			}
			if (strcmp(info->trans_gate_source, "blind_current_posterior_confidence") != 0) {
				add_example(audit, "manifest trans gate source is not blind posterior confidence");
				failed = 1;
			}
			if (strcmp(info->trans_gate_mode, "off") == 0 &&
				(!check_close(info->trans_gate_min_pmax, 1.0) ||
				 !check_close(info->trans_gate_min_margin, 1.0) ||
				 !check_close(info->trans_gate_min_neg_entropy, -log((double)HK_BLIND_N_STATE)))) {
				add_example(audit, "manifest trans gate off defaults are inconsistent");
				failed = 1;
			}
			if (info->seen[MK_TRANS_GATE_USES_PHASE_LABELS] &&
				info->trans_gate_uses_phase_labels != 0) {
				add_example(audit, "manifest trans gate uses phase labels");
				failed = 1;
			}
			if (info->seen[MK_TRANS_GATE_USES_CHARM_OR_REFERENCE] &&
				info->trans_gate_uses_charm_or_reference != 0) {
				add_example(audit, "manifest trans gate uses CHARM/reference");
				failed = 1;
			}
		}
		if (strcmp(info->init_mode, "random_diploid") == 0) {
			if (!check_close(info->init_eps_effective, 0.0) ||
				info->init_scale_effective <= 0.0) {
			add_example(audit, "manifest random_diploid init metadata is inconsistent");
			failed = 1;
		}
	} else if (strcmp(info->init_mode, "random_haploid_split") == 0) {
		if (info->init_eps_effective <= 0.0 || info->init_scale_effective <= 0.0) {
			add_example(audit, "manifest random_haploid_split init metadata is inconsistent");
			failed = 1;
		}
	} else if (strcmp(info->init_mode, "unphased_scaffold_split") == 0 ||
			   strcmp(info->init_mode, "toy_split") == 0) {
		if (info->init_eps_effective <= 0.0 ||
			!check_close(info->init_scale_effective, 0.0)) {
			add_example(audit, "manifest split init metadata is inconsistent");
			failed = 1;
		}
	}
	if (!isfinite(info->d_scale_eps_count) || info->d_scale_eps_count <= 0.0) {
		add_example(audit, "manifest d_scale_eps_count out of range");
		failed = 1;
	}
	if (!isfinite(info->final_mean_entropy) || info->final_mean_entropy < -HK_AUDIT_TOL ||
		info->final_mean_entropy > log((double)HK_BLIND_N_STATE) + HK_AUDIT_TOL) {
		add_example_fmt(audit, "manifest", 0, "bad final_mean_entropy", info->final_mean_entropy);
		failed = 1;
	}
	if (info->n_bad_iter != 0 || info->n_relax_nonfinite_iter != 0 || info->n_coord_nonfinite != 0) {
		add_example(audit, "manifest bad/nonfinite counters are nonzero");
		failed = 1;
	}
	if (!in_unit_range(info->final_mean_pU)) {
		add_example_fmt(audit, "manifest", 0, "bad final_mean_pU", info->final_mean_pU);
		failed = 1;
	}
	if (!isfinite(info->final_sum_wedge_k) || info->final_sum_wedge_k < 0.0) {
		add_example_fmt(audit, "manifest", 0, "bad final_sum_wedge_k", info->final_sum_wedge_k);
		failed = 1;
	}
	const int trans_pair_prior_on = info->seen[MK_TRANS_CHR_PAIR_PRIOR_LAMBDA] &&
		info->trans_chr_pair_prior_lambda > 0.0;
	const char *expected_refresh_prior_mode = trans_pair_prior_on?
		"blind_posterior_chrom_pair" : info->prior_mode;
	if (info->posterior_refreshed_after_final_relax != 1 ||
		!isfinite(info->posterior_refresh_temperature) ||
		info->posterior_refresh_temperature <= 0.0 ||
		strcmp(info->posterior_refresh_prior_mode, expected_refresh_prior_mode) != 0 ||
		!isfinite(info->posterior_refresh_mean_kl) ||
		info->posterior_refresh_mean_kl < -HK_AUDIT_TOL ||
		!in_unit_range(info->posterior_refresh_top_state_switch_frac) ||
		!in_unit_range(info->posterior_refresh_mean_pU_before) ||
		!in_unit_range(info->posterior_refresh_mean_pU_after)) {
		add_example(audit, "manifest posterior refresh fields are invalid");
		failed = 1;
	}
	if (info->write_raw_posterior != 0 && info->write_raw_posterior != 1) {
		add_example(audit, "manifest write_raw_posterior is not 0/1");
		failed = 1;
	}
	if (info->enable_repulsion != 1) {
		add_example(audit, "manifest enable_repulsion is not 1");
		failed = 1;
	}
	if (info->repulsion_mode != HK_BLIND_REPULSION_N2 &&
		info->repulsion_mode != HK_BLIND_REPULSION_CELL) {
		add_example(audit, "manifest repulsion_mode is not N2/CELL");
		failed = 1;
	}
	if (strcmp(info->repulsion_blocking_mode, "current_edge_blocking") != 0) {
		add_example(audit, "manifest repulsion_blocking_mode is not recognized");
		failed = 1;
	}
	if (info->write_raw_posterior == 1 && !info->has_output_raw_posterior) {
		add_example(audit, "manifest raw output flag set but output_raw_posterior missing");
		failed = 1;
	}
	if (audit->forbidden_hits > 0) failed = 1;
	if (failed && audit->bad_rows == 0) audit->bad_rows = 1;
	return failed;
}

static int manifest_has_heldout_fields(const struct manifest_info *info)
{
	assert(info);
	return info->seen[MK_HELDOUT_ENABLED] ||
		   info->seen[MK_HELDOUT_FRACTION] ||
		   info->seen[MK_HELDOUT_SEED] ||
		   info->seen[MK_HELDOUT_N_INPUT_RAW_TOTAL] ||
		   info->seen[MK_HELDOUT_N_RAW_TRAIN] ||
		   info->seen[MK_HELDOUT_N_RAW_HELDOUT] ||
		   info->seen[MK_HELDOUT_N_BPAIR_HELDOUT] ||
		   info->seen[MK_HELDOUT_N_BPAIR_EVAL] ||
		   info->seen[MK_HELDOUT_N_RAW_EVAL] ||
		   info->seen[MK_HELDOUT_MEAN_EXPECTED_ENERGY] ||
		   info->seen[MK_HELDOUT_MEAN_MIN_ENERGY] ||
		   info->seen[MK_HELDOUT_MEAN_ENTROPY] ||
		   info->seen[MK_HELDOUT_MEAN_PU] ||
		   info->seen[MK_HELDOUT_MEAN_BEST_NORMALIZED_DISTANCE] ||
		   info->seen[MK_HELDOUT_SHORT_DISTANCE_FRAC];
}

static void mark_bad_row(struct file_audit *audit, int *row_bad)
{
	if (!*row_bad) {
		++audit->bad_rows;
		*row_bad = 1;
	}
}

static int base_k_matches_manifest(const struct manifest_info *info, double base_k)
{
	if (strcmp(info->base_k_mode, "uniform") == 0)
		return check_close(base_k, info->base_k_effective);
	return base_k >= info->base_k_min - HK_AUDIT_TOL &&
		   base_k <= info->base_k_max + HK_AUDIT_TOL;
}

static int audit_bpair_file(const char *path, const struct manifest_info *info, struct file_audit *audit)
{
	static const char *const header_tokens[] = {
		"chr1", "start1", "end1", "chr2", "start2", "end2", "bid1", "bid2",
		"n_raw", "base_d_scale", "base_k", "p00", "p01", "p10", "p11", "pU",
		"psame_raw", "pcross_raw", "entropy", "margin", "pmax", "rho_output"
	};
	FILE *fp = fopen(path, "r");
	char line[HK_AUDIT_LINE_MAX];
	int failed = 0;
	int64_t row = 0;

	audit_init(audit);
	if (fp == 0) {
		add_example(audit, "bpair posterior open failed");
		audit->bad_rows++;
		return 1;
	}
	if (fgets(line, sizeof(line), fp) == 0) {
		add_example(audit, "bpair posterior missing header");
		audit->bad_rows++;
		fclose(fp);
		return 1;
	}
	scan_forbidden_line(audit, "bpair", 0, line);
	failed |= check_header_tokens(audit, "bpair", line, header_tokens, (int)(sizeof(header_tokens) / sizeof(header_tokens[0])));
	while (fgets(line, sizeof(line), fp)) {
		char chr1[128], chr2[128];
		char contact_class[64];
		int st1, en1, st2, en2, bid1, bid2, n_raw;
		double base_d_scale, base_k, expected_d_scale;
		double p00, p01, p10, p11, pU, psame, pcross, entropy, margin, pmax, rho;
		double sum_p4;
		int row_bad = 0;
		int n, offset = 0, has_contact_class = 0;
		++row;
		scan_forbidden_line(audit, "bpair", row, line);
		n = sscanf(line, "%127[^\t]\t%d\t%d\t%127[^\t]\t%d\t%d\t%d\t%d\t%d\t"
				   "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf%n",
				   chr1, &st1, &en1, chr2, &st2, &en2, &bid1, &bid2, &n_raw,
				   &base_d_scale, &base_k, &p00, &p01, &p10, &p11, &pU, &psame, &pcross,
				   &entropy, &margin, &pmax, &rho, &offset);
		if (n != 22) {
			add_example_fmt(audit, "bpair", row, "parse failed", (double)n);
			++audit->bad_rows;
			continue;
		}
		if (offset > 0 && line[offset] == '\t') {
			if (sscanf(line + offset + 1, "%63[^\t\r\n]", contact_class) == 1) {
				has_contact_class = 1;
			} else {
				add_example_fmt(audit, "bpair", row, "bad contact_class", NAN);
				mark_bad_row(audit, &row_bad);
			}
		}
		sum_p4 = p00 + p01 + p10 + p11;
		if (bid1 < 0 || bid1 >= info->n_beads || bid2 < 0 || bid2 >= info->n_beads) {
			add_example_fmt(audit, "bpair", row, "bad bid", (double)bid1);
			mark_bad_row(audit, &row_bad);
		}
		if (n_raw <= 0) {
			add_example_fmt(audit, "bpair", row, "bad n_raw", (double)n_raw);
			mark_bad_row(audit, &row_bad);
		}
		if (!isfinite(base_d_scale) || base_d_scale <= 0.0) {
			add_example_fmt(audit, "bpair", row, "bad base_d_scale", base_d_scale);
			mark_bad_row(audit, &row_bad);
		}
		if (!isfinite(base_k) || base_k < 0.0) {
			add_example_fmt(audit, "bpair", row, "bad base_k", base_k);
			mark_bad_row(audit, &row_bad);
		}
		if (n_raw > 0 && isfinite(base_d_scale) && base_d_scale > 0.0) {
			expected_d_scale = expected_base_d_scale_from_n_raw(n_raw);
			if (!check_close(base_d_scale, expected_d_scale)) {
				add_example_fmt(audit, "bpair", row, "bad base_d_scale", base_d_scale);
				mark_bad_row(audit, &row_bad);
			}
		}
		if (isfinite(base_k) && !base_k_matches_manifest(info, base_k)) {
			add_example_fmt(audit, "bpair", row, "bad base_k", base_k);
			mark_bad_row(audit, &row_bad);
		}
		if (!isfinite(p00) || !isfinite(p01) || !isfinite(p10) || !isfinite(p11) ||
			!isfinite(pU) || !isfinite(entropy) || !isfinite(margin) ||
			!isfinite(pmax) || !isfinite(rho)) {
			add_example_fmt(audit, "bpair", row, "nonfinite posterior field", NAN);
			mark_bad_row(audit, &row_bad);
		}
		if (!in_unit_range(p00) || !in_unit_range(p01) || !in_unit_range(p10) ||
			!in_unit_range(p11) || !in_unit_range(pU) || !in_unit_range(rho)) {
			add_example_fmt(audit, "bpair", row, "probability out of range", sum_p4);
			mark_bad_row(audit, &row_bad);
		}
		if (!check_close(sum_p4, 1.0)) {
			add_example_fmt(audit, "bpair", row, "bad p4 sum", sum_p4);
			mark_bad_row(audit, &row_bad);
		}
		if (!check_close(rho * sum_p4 + pU, 1.0)) {
			add_example_fmt(audit, "bpair", row, "bad five-state sum", rho * sum_p4 + pU);
			mark_bad_row(audit, &row_bad);
		}
		if (!check_close(psame, p00 + p11)) {
			add_example_fmt(audit, "bpair", row, "bad psame_raw", psame);
			mark_bad_row(audit, &row_bad);
		}
		if (!check_close(pcross, p01 + p10)) {
			add_example_fmt(audit, "bpair", row, "bad pcross_raw", pcross);
			mark_bad_row(audit, &row_bad);
		}
		if (bid1 == bid2) {
			++audit->same_bin_rows;
			if (info->uses_phase_labels == 0 &&
				(!check_close(p00, 0.25) || !check_close(p01, 0.25) ||
				!check_close(p10, 0.25) || !check_close(p11, 0.25) ||
				 !check_close(pU, 1.0) || !check_close(rho, 0.0))) {
				add_example_fmt(audit, "bpair", row, "same-bin posterior not uniform unknown", pU);
				mark_bad_row(audit, &row_bad);
			}
		}
		if (has_contact_class) {
			const char *expected_class = strcmp(chr1, chr2) == 0? "cis" : "trans";
			if (strcmp(contact_class, "cis") != 0 && strcmp(contact_class, "trans") != 0) {
				add_example_fmt(audit, "bpair", row, "bad contact_class", NAN);
				mark_bad_row(audit, &row_bad);
			} else if (strcmp(contact_class, expected_class) != 0) {
				add_example_fmt(audit, "bpair", row, "contact_class disagrees with chr", NAN);
				mark_bad_row(audit, &row_bad);
			}
		}
	}
	fclose(fp);
	audit->rows = row;
	if (row != info->n_bpair) {
		add_example_fmt(audit, "bpair", row, "row count mismatch", (double)row);
		failed = 1;
	}
	if (audit->same_bin_rows != info->n_bpair_same_bin_excluded) {
		add_example_fmt(audit, "bpair", row, "same-bin row count mismatch",
						(double)audit->same_bin_rows);
		failed = 1;
	}
	if (audit->bad_rows > 0 || audit->forbidden_hits > 0) failed = 1;
	return failed;
}

static int audit_coords_file(const char *path, const struct manifest_info *info, struct file_audit *audit)
{
	static const char *const header_tokens[] = {
		"chr", "start", "end", "bid", "copy", "diploid_bid", "x", "y", "z"
	};
	FILE *fp = fopen(path, "r");
	char line[HK_AUDIT_LINE_MAX];
	uint8_t *copy_mask = 0;
	int failed = 0;
	int64_t row = 0;
	int64_t i;

	audit_init(audit);
	if (fp == 0) {
		add_example(audit, "coords open failed");
		audit->bad_rows++;
		return 1;
	}
	copy_mask = (uint8_t*)calloc((size_t)info->n_beads, sizeof(*copy_mask));
	if (copy_mask == 0) {
		add_example(audit, "coords copy mask allocation failed");
		audit->bad_rows++;
		fclose(fp);
		return 1;
	}
	if (fgets(line, sizeof(line), fp) == 0) {
		add_example(audit, "coords missing header");
		audit->bad_rows++;
		free(copy_mask);
		fclose(fp);
		return 1;
	}
	scan_forbidden_line(audit, "coords", 0, line);
	failed |= check_header_tokens(audit, "coords", line, header_tokens, (int)(sizeof(header_tokens) / sizeof(header_tokens[0])));
	while (fgets(line, sizeof(line), fp)) {
		char chr[128];
		int st, en, bid, copy, diploid_bid;
		double x, y, z;
		int row_bad = 0;
		int n;
		++row;
		scan_forbidden_line(audit, "coords", row, line);
		n = sscanf(line, "%127[^\t]\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf",
				   chr, &st, &en, &bid, &copy, &diploid_bid, &x, &y, &z);
		if (n != 9) {
			add_example_fmt(audit, "coords", row, "parse failed", (double)n);
			++audit->bad_rows;
			continue;
		}
		if (bid < 0 || bid >= info->n_beads) {
			add_example_fmt(audit, "coords", row, "bad bid", (double)bid);
			mark_bad_row(audit, &row_bad);
		}
		if (copy != 0 && copy != 1) {
			add_example_fmt(audit, "coords", row, "bad copy", (double)copy);
			mark_bad_row(audit, &row_bad);
		}
		if (bid >= 0 && bid < info->n_beads && (copy == 0 || copy == 1)) {
			uint8_t bit = (uint8_t)(1u << copy);
			if (copy_mask[bid] & bit) {
				add_example_fmt(audit, "coords", row, "duplicate bid/copy", (double)bid);
				mark_bad_row(audit, &row_bad);
			}
			copy_mask[bid] |= bit;
			if (diploid_bid != hk_diploid_bid(bid, copy)) {
				add_example_fmt(audit, "coords", row, "bad diploid_bid", (double)diploid_bid);
				mark_bad_row(audit, &row_bad);
			}
		}
		if (!isfinite(x) || !isfinite(y) || !isfinite(z)) {
			add_example_fmt(audit, "coords", row, "nonfinite coordinate", NAN);
			mark_bad_row(audit, &row_bad);
		}
	}
	fclose(fp);
	audit->rows = row;
	if (row != 2 * info->n_beads) {
		add_example_fmt(audit, "coords", row, "row count mismatch", (double)row);
		failed = 1;
	}
	for (i = 0; i < info->n_beads; ++i) {
		if (copy_mask[i] != 3) {
			char buf[256];
			snprintf(buf, sizeof(buf), "coords bid %lld missing copy mask %u", (long long)i, (unsigned)copy_mask[i]);
			add_example(audit, buf);
			++audit->bad_rows;
			failed = 1;
			break;
		}
	}
	free(copy_mask);
	if (audit->bad_rows > 0 || audit->forbidden_hits > 0) failed = 1;
	return failed;
}

static int audit_loop_diag_file(const char *path, const struct manifest_info *info, struct file_audit *audit)
{
	static const char *const header_tokens[] = {
		"n_iter", "n_completed", "final_mean_entropy", "final_mean_pU",
		"n_bad_iter", "n_relax_nonfinite_iter", "n_coord_nonfinite", "final_sum_wedge_k",
		"final_mean_rho_train_bpair", "final_n_skipped_same_bin_bpairs",
		"final_repulsion_energy", "final_repulsion_force_l1", "n_repulsion_nonfinite_step",
		"repulsion_mode", "posterior_refreshed_after_final_relax", "posterior_refresh_temperature",
		"final_copytrack_energy", "final_copytrack_force_l1", "n_copytrack_nonfinite_step",
		"final_global_copytrack_energy", "final_global_copytrack_force_l1", "n_global_copytrack_nonfinite_step",
		"final_normdir_copytrack_energy", "final_normdir_copytrack_force_l1", "n_normdir_copytrack_nonfinite_step"
	};
	FILE *fp = fopen(path, "r");
	char header_line[HK_AUDIT_LINE_MAX];
	char value_line[HK_AUDIT_LINE_MAX];
	char *header_fields[HK_AUDIT_MAX_FIELDS];
	char *value_fields[HK_AUDIT_MAX_FIELDS];
	int n_header, n_value;
	int failed = 0;
	int row_bad = 0;
	int n_iter, n_completed, total_chr_flipped, n_bad_iter, n_relax_bad, n_coord_bad;
	double expected_final_temperature, expected_final_rho;
	double initial_entropy, final_entropy, initial_pU, final_pU;
	double initial_temperature, final_temperature, initial_rho, final_rho;
	double mean_sep, min_sep, max_sep, sum_k;
	double mean_rho_train, min_rho_train, max_rho_train;
	long long n_skipped_same_bin;
	long long n_rep_considered, n_rep_blocked, n_rep_active;
	double rep_energy, rep_force;
	int n_rep_bad, rep_mode, posterior_refreshed;
	double copytrack_energy, copytrack_force;
	double global_copytrack_energy, global_copytrack_force;
	double normdir_copytrack_energy, normdir_copytrack_force;
	int n_copytrack_bad, n_global_copytrack_bad, n_normdir_copytrack_bad;
	double refresh_temperature, refresh_kl, refresh_switch_frac, refresh_pu_before, refresh_pu_after;

	audit_init(audit);
	if (fp == 0) {
		add_example(audit, "loop diag open failed");
		audit->bad_rows++;
		return 1;
	}
	if (fgets(header_line, sizeof(header_line), fp) == 0) {
		add_example(audit, "loop diag missing header");
		audit->bad_rows++;
		fclose(fp);
		return 1;
	}
	scan_forbidden_line(audit, "loop_diag", 0, header_line);
	failed |= check_header_tokens(audit, "loop_diag", header_line, header_tokens, (int)(sizeof(header_tokens) / sizeof(header_tokens[0])));
	if (fgets(value_line, sizeof(value_line), fp) == 0) {
		add_example(audit, "loop diag missing data row");
		audit->bad_rows++;
		fclose(fp);
		return 1;
	}
	scan_forbidden_line(audit, "loop_diag", 1, value_line);
	audit->rows = 1;
	trim_line(header_line);
	trim_line(value_line);
	n_header = split_tsv_line(header_line, header_fields, HK_AUDIT_MAX_FIELDS);
	n_value = split_tsv_line(value_line, value_fields, HK_AUDIT_MAX_FIELDS);
	if (tsv_get_i32(header_fields, value_fields, n_header, n_value, "n_iter", &n_iter) != 0 ||
		tsv_get_i32(header_fields, value_fields, n_header, n_value, "n_completed", &n_completed) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "initial_mean_entropy", &initial_entropy) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_mean_entropy", &final_entropy) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "initial_mean_pU", &initial_pU) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_mean_pU", &final_pU) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "initial_temperature", &initial_temperature) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_temperature", &final_temperature) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "initial_rho_train", &initial_rho) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_rho_train", &final_rho) != 0 ||
		tsv_get_i32(header_fields, value_fields, n_header, n_value, "total_chr_flipped", &total_chr_flipped) != 0 ||
		tsv_get_i32(header_fields, value_fields, n_header, n_value, "n_bad_iter", &n_bad_iter) != 0 ||
		tsv_get_i32(header_fields, value_fields, n_header, n_value, "n_relax_nonfinite_iter", &n_relax_bad) != 0 ||
		tsv_get_i32(header_fields, value_fields, n_header, n_value, "n_coord_nonfinite", &n_coord_bad) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_mean_sep", &mean_sep) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_min_sep", &min_sep) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_max_sep", &max_sep) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_sum_wedge_k", &sum_k) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_mean_rho_train_bpair", &mean_rho_train) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_min_rho_train_bpair", &min_rho_train) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_max_rho_train_bpair", &max_rho_train) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "final_n_skipped_same_bin_bpairs", &n_skipped_same_bin) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_repulsion_energy", &rep_energy) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_repulsion_force_l1", &rep_force) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "final_n_repulsion_pairs_considered", &n_rep_considered) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "final_n_repulsion_pairs_blocked", &n_rep_blocked) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "final_n_repulsion_pairs_active", &n_rep_active) != 0 ||
		tsv_get_i32(header_fields, value_fields, n_header, n_value, "n_repulsion_nonfinite_step", &n_rep_bad) != 0 ||
		tsv_get_i32(header_fields, value_fields, n_header, n_value, "repulsion_mode", &rep_mode) != 0 ||
		tsv_get_i32(header_fields, value_fields, n_header, n_value, "posterior_refreshed_after_final_relax", &posterior_refreshed) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "posterior_refresh_temperature", &refresh_temperature) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "posterior_refresh_mean_kl", &refresh_kl) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "posterior_refresh_top_state_switch_frac", &refresh_switch_frac) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "posterior_refresh_mean_pU_before", &refresh_pu_before) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "posterior_refresh_mean_pU_after", &refresh_pu_after) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_copytrack_energy", &copytrack_energy) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_copytrack_force_l1", &copytrack_force) != 0 ||
		tsv_get_i32(header_fields, value_fields, n_header, n_value, "n_copytrack_nonfinite_step", &n_copytrack_bad) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_global_copytrack_energy", &global_copytrack_energy) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_global_copytrack_force_l1", &global_copytrack_force) != 0 ||
		tsv_get_i32(header_fields, value_fields, n_header, n_value, "n_global_copytrack_nonfinite_step", &n_global_copytrack_bad) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_normdir_copytrack_energy", &normdir_copytrack_energy) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "final_normdir_copytrack_force_l1", &normdir_copytrack_force) != 0 ||
		tsv_get_i32(header_fields, value_fields, n_header, n_value, "n_normdir_copytrack_nonfinite_step", &n_normdir_copytrack_bad) != 0) {
		add_example(audit, "loop diag parse failed");
		++audit->bad_rows;
		fclose(fp);
		return 1;
	}
	if (fgets(value_line, sizeof(value_line), fp) != 0) {
		add_example(audit, "loop diag has more than one data row");
		++audit->bad_rows;
		failed = 1;
	}
	fclose(fp);
	expected_final_temperature = info->n_iter <= 1? info->temperature_start : info->temperature_end;
	expected_final_rho = info->n_iter <= 1? info->rho_train_start : info->rho_train_end;
	if (n_iter != info->n_iter || n_completed != info->n_iter) {
		add_example_fmt(audit, "loop_diag", 1, "bad n_completed", (double)n_completed);
		mark_bad_row(audit, &row_bad);
	}
	if (!check_close(initial_temperature, info->temperature_start) ||
		!check_close(final_temperature, expected_final_temperature) ||
		!check_close(initial_rho, info->rho_train_start) ||
		!check_close(final_rho, expected_final_rho)) {
		add_example(audit, "loop diag schedule values differ from manifest");
		mark_bad_row(audit, &row_bad);
	}
	if (n_bad_iter != 0 || n_relax_bad != 0 || n_coord_bad != 0) {
		add_example(audit, "loop diag bad counters are nonzero");
		mark_bad_row(audit, &row_bad);
	}
	if (n_rep_bad != 0) {
		add_example(audit, "loop diag repulsion bad counter is nonzero");
		mark_bad_row(audit, &row_bad);
	}
	if (n_copytrack_bad != 0) {
		add_example(audit, "loop diag copytrack bad counter is nonzero");
		mark_bad_row(audit, &row_bad);
	}
	if (n_global_copytrack_bad != 0) {
		add_example(audit, "loop diag global copytrack bad counter is nonzero");
		mark_bad_row(audit, &row_bad);
	}
	if (n_normdir_copytrack_bad != 0) {
		add_example(audit, "loop diag normdir copytrack bad counter is nonzero");
		mark_bad_row(audit, &row_bad);
	}
	if (rep_mode != HK_BLIND_REPULSION_N2 && rep_mode != HK_BLIND_REPULSION_CELL) {
		add_example(audit, "loop diag repulsion_mode is not N2/CELL");
		mark_bad_row(audit, &row_bad);
	}
	if (rep_mode != info->repulsion_mode) {
		add_example(audit, "loop diag repulsion_mode differs from manifest");
		mark_bad_row(audit, &row_bad);
	}
	if (!in_unit_range(final_pU)) {
		add_example_fmt(audit, "loop_diag", 1, "bad final_mean_pU", final_pU);
		mark_bad_row(audit, &row_bad);
	}
	if (!isfinite(sum_k) || sum_k < 0.0) {
		add_example_fmt(audit, "loop_diag", 1, "bad final_sum_wedge_k", sum_k);
		mark_bad_row(audit, &row_bad);
	}
	if (!isfinite(final_entropy) || !isfinite(mean_sep) || !isfinite(min_sep) || !isfinite(max_sep)) {
		add_example(audit, "loop diag nonfinite summary value");
		mark_bad_row(audit, &row_bad);
	}
	if (!isfinite(rep_energy) || rep_energy < 0.0 || !isfinite(rep_force) || rep_force < 0.0) {
		add_example(audit, "loop diag bad repulsion summary value");
		mark_bad_row(audit, &row_bad);
	}
	if (!isfinite(copytrack_energy) || copytrack_energy < 0.0 ||
		!isfinite(copytrack_force) || copytrack_force < 0.0) {
		add_example(audit, "loop diag bad copytrack summary value");
		mark_bad_row(audit, &row_bad);
	}
	if (!isfinite(global_copytrack_energy) || global_copytrack_energy < 0.0 ||
		!isfinite(global_copytrack_force) || global_copytrack_force < 0.0) {
		add_example(audit, "loop diag bad global copytrack summary value");
		mark_bad_row(audit, &row_bad);
	}
	if (!isfinite(normdir_copytrack_energy) || normdir_copytrack_energy < 0.0 ||
		!isfinite(normdir_copytrack_force) || normdir_copytrack_force < 0.0) {
		add_example(audit, "loop diag bad normdir copytrack summary value");
		mark_bad_row(audit, &row_bad);
	}
	if (n_rep_considered < 0 || n_rep_blocked < 0 || n_rep_active < 0 ||
		n_rep_blocked > n_rep_considered || n_rep_active > n_rep_considered) {
		add_example(audit, "loop diag bad repulsion pair counters");
		mark_bad_row(audit, &row_bad);
	}
	if (!isfinite(mean_rho_train) || !isfinite(min_rho_train) || !isfinite(max_rho_train) ||
		min_rho_train < -HK_AUDIT_TOL || max_rho_train < min_rho_train - HK_AUDIT_TOL) {
		add_example(audit, "loop diag bad rho_train_bpair summary value");
		mark_bad_row(audit, &row_bad);
	}
	if (n_skipped_same_bin < 0) {
		add_example(audit, "loop diag bad same-bin skipped count");
		mark_bad_row(audit, &row_bad);
	}
	if (posterior_refreshed != 1 || !isfinite(refresh_temperature) ||
		refresh_temperature <= 0.0 || !isfinite(refresh_kl) || refresh_kl < -HK_AUDIT_TOL ||
		!in_unit_range(refresh_switch_frac) || !in_unit_range(refresh_pu_before) ||
		!in_unit_range(refresh_pu_after)) {
		add_example(audit, "loop diag bad posterior refresh fields");
		mark_bad_row(audit, &row_bad);
	}
	if (!check_close(refresh_temperature, info->posterior_refresh_temperature) ||
		!check_close(refresh_kl, info->posterior_refresh_mean_kl) ||
		!check_close(refresh_switch_frac, info->posterior_refresh_top_state_switch_frac) ||
		!check_close(refresh_pu_before, info->posterior_refresh_mean_pU_before) ||
		!check_close(refresh_pu_after, info->posterior_refresh_mean_pU_after)) {
		add_example(audit, "loop diag posterior refresh fields differ from manifest");
		mark_bad_row(audit, &row_bad);
	}
	if (audit->bad_rows > 0 || audit->forbidden_hits > 0) failed = 1;
	return failed;
}

static int audit_raw_file(const char *path, const struct manifest_info *info, struct file_audit *audit)
{
	static const char *const header_tokens[] = {
		"raw_id", "chr1", "pos1", "chr2", "pos2", "bid1_raw", "bid2_raw",
		"bpair_id", "bid1_canonical", "bid2_canonical", "swapped", "p00", "p01",
		"p10", "p11", "pU", "psame_raw", "pcross_raw", "entropy", "margin",
		"pmax", "rho_output"
	};
	FILE *fp = fopen(path, "r");
	char line[HK_AUDIT_LINE_MAX];
	int failed = 0;
	int64_t row = 0;

	audit_init(audit);
	if (fp == 0) {
		add_example(audit, "raw posterior open failed");
		audit->bad_rows++;
		return 1;
	}
	if (fgets(line, sizeof(line), fp) == 0) {
		add_example(audit, "raw posterior missing header");
		audit->bad_rows++;
		fclose(fp);
		return 1;
	}
	scan_forbidden_line(audit, "raw", 0, line);
	failed |= check_header_tokens(audit, "raw", line, header_tokens, (int)(sizeof(header_tokens) / sizeof(header_tokens[0])));
	while (fgets(line, sizeof(line), fp)) {
		char chr1[128], chr2[128];
		int raw_id, pos1, pos2, bid1_raw, bid2_raw, bpair_id, bid1_can, bid2_can, swapped;
		double p00, p01, p10, p11, pU, psame, pcross, entropy, margin, pmax, rho;
		double sum_p4;
		int row_bad = 0;
		int n;
		++row;
		scan_forbidden_line(audit, "raw", row, line);
		n = sscanf(line, "%d\t%127[^\t]\t%d\t%127[^\t]\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t"
				   "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				   &raw_id, chr1, &pos1, chr2, &pos2, &bid1_raw, &bid2_raw,
				   &bpair_id, &bid1_can, &bid2_can, &swapped,
				   &p00, &p01, &p10, &p11, &pU, &psame, &pcross,
				   &entropy, &margin, &pmax, &rho);
		if (n != 22) {
			add_example_fmt(audit, "raw", row, "parse failed", (double)n);
			++audit->bad_rows;
			continue;
		}
		sum_p4 = p00 + p01 + p10 + p11;
		if (raw_id < 0 || raw_id >= info->n_raw) {
			add_example_fmt(audit, "raw", row, "bad raw_id", (double)raw_id);
			mark_bad_row(audit, &row_bad);
		}
		if (bpair_id < 0 || bpair_id >= info->n_bpair) {
			add_example_fmt(audit, "raw", row, "bad bpair_id", (double)bpair_id);
			mark_bad_row(audit, &row_bad);
		}
		if (swapped != 0 && swapped != 1) {
			add_example_fmt(audit, "raw", row, "bad swapped", (double)swapped);
			mark_bad_row(audit, &row_bad);
		}
		if (bid1_can < 0 || bid1_can >= info->n_beads || bid2_can < 0 ||
			bid2_can >= info->n_beads || bid1_can > bid2_can) {
			add_example_fmt(audit, "raw", row, "bad canonical bid", (double)bid1_can);
			mark_bad_row(audit, &row_bad);
		}
		if (!isfinite(p00) || !isfinite(p01) || !isfinite(p10) || !isfinite(p11) ||
			!isfinite(pU) || !isfinite(entropy) || !isfinite(margin) ||
			!isfinite(pmax) || !isfinite(rho)) {
			add_example_fmt(audit, "raw", row, "nonfinite posterior field", NAN);
			mark_bad_row(audit, &row_bad);
		}
		if (!in_unit_range(p00) || !in_unit_range(p01) || !in_unit_range(p10) ||
			!in_unit_range(p11) || !in_unit_range(pU) || !in_unit_range(rho)) {
			add_example_fmt(audit, "raw", row, "probability out of range", sum_p4);
			mark_bad_row(audit, &row_bad);
		}
		if (!check_close(sum_p4, 1.0)) {
			add_example_fmt(audit, "raw", row, "bad p4 sum", sum_p4);
			mark_bad_row(audit, &row_bad);
		}
		if (!check_close(rho * sum_p4 + pU, 1.0)) {
			add_example_fmt(audit, "raw", row, "bad five-state sum", rho * sum_p4 + pU);
			mark_bad_row(audit, &row_bad);
		}
		if (!check_close(psame, p00 + p11)) {
			add_example_fmt(audit, "raw", row, "bad psame_raw", psame);
			mark_bad_row(audit, &row_bad);
		}
		if (!check_close(pcross, p01 + p10)) {
			add_example_fmt(audit, "raw", row, "bad pcross_raw", pcross);
			mark_bad_row(audit, &row_bad);
		}
		if (bid1_can == bid2_can) {
			++audit->same_bin_rows;
			if (info->uses_phase_labels == 0 &&
				(!check_close(p00, 0.25) || !check_close(p01, 0.25) ||
				!check_close(p10, 0.25) || !check_close(p11, 0.25) ||
				 !check_close(pU, 1.0) || !check_close(rho, 0.0))) {
				add_example_fmt(audit, "raw", row, "same-bin posterior not uniform unknown", pU);
				mark_bad_row(audit, &row_bad);
			}
		}
	}
	fclose(fp);
	audit->rows = row;
	if (row != info->n_raw) {
		add_example_fmt(audit, "raw", row, "row count mismatch", (double)row);
		failed = 1;
	}
	if (audit->same_bin_rows != info->n_raw_same_bin_excluded) {
		add_example_fmt(audit, "raw", row, "same-bin row count mismatch",
						(double)audit->same_bin_rows);
		failed = 1;
	}
	if (audit->bad_rows > 0 || audit->forbidden_hits > 0) failed = 1;
	return failed;
}

static int audit_heldout_diag_file(const char *path, const struct manifest_info *info, struct file_audit *audit)
{
	static const char *const header_tokens[] = {
		"enabled", "fraction", "seed", "n_input_raw_total", "n_raw_train", "n_raw_heldout",
		"n_bpair_heldout", "n_bpair_eval", "n_raw_eval", "n_bpair_same_bin_skipped",
		"n_raw_same_bin_skipped", "mean_expected_energy", "mean_min_energy", "mean_entropy",
		"mean_pU", "mean_best_normalized_distance", "short_distance_frac"
	};
	FILE *fp = fopen(path, "r");
	char header_line[HK_AUDIT_LINE_MAX];
	char value_line[HK_AUDIT_LINE_MAX];
	char *header_fields[HK_AUDIT_MAX_FIELDS];
	char *value_fields[HK_AUDIT_MAX_FIELDS];
	int n_header, n_value;
	int failed = 0;
	int row_bad = 0;
	int enabled;
	long long seed, n_input_raw_total, n_raw_train, n_raw_heldout;
	long long n_bpair_heldout, n_bpair_eval, n_raw_eval;
	long long n_bpair_same_bin_skipped, n_raw_same_bin_skipped;
	double fraction, mean_expected_energy, mean_min_energy, mean_entropy;
	double mean_pU, mean_best_normalized_distance, short_distance_frac;

	audit_init(audit);
	if (fp == 0) {
		add_example(audit, "heldout diag open failed");
		audit->bad_rows++;
		return 1;
	}
	if (fgets(header_line, sizeof(header_line), fp) == 0) {
		add_example(audit, "heldout diag missing header");
		audit->bad_rows++;
		fclose(fp);
		return 1;
	}
	scan_forbidden_line(audit, "heldout_diag", 0, header_line);
	failed |= check_header_tokens(audit, "heldout_diag", header_line, header_tokens,
								  (int)(sizeof(header_tokens) / sizeof(header_tokens[0])));
	if (fgets(value_line, sizeof(value_line), fp) == 0) {
		add_example(audit, "heldout diag missing data row");
		audit->bad_rows++;
		fclose(fp);
		return 1;
	}
	scan_forbidden_line(audit, "heldout_diag", 1, value_line);
	audit->rows = 1;
	trim_line(header_line);
	trim_line(value_line);
	n_header = split_tsv_line(header_line, header_fields, HK_AUDIT_MAX_FIELDS);
	n_value = split_tsv_line(value_line, value_fields, HK_AUDIT_MAX_FIELDS);
	if (tsv_get_i32(header_fields, value_fields, n_header, n_value, "enabled", &enabled) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "fraction", &fraction) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "seed", &seed) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "n_input_raw_total", &n_input_raw_total) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "n_raw_train", &n_raw_train) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "n_raw_heldout", &n_raw_heldout) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "n_bpair_heldout", &n_bpair_heldout) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "n_bpair_eval", &n_bpair_eval) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "n_raw_eval", &n_raw_eval) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "n_bpair_same_bin_skipped", &n_bpair_same_bin_skipped) != 0 ||
		tsv_get_i64(header_fields, value_fields, n_header, n_value, "n_raw_same_bin_skipped", &n_raw_same_bin_skipped) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "mean_expected_energy", &mean_expected_energy) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "mean_min_energy", &mean_min_energy) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "mean_entropy", &mean_entropy) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "mean_pU", &mean_pU) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "mean_best_normalized_distance", &mean_best_normalized_distance) != 0 ||
		tsv_get_double(header_fields, value_fields, n_header, n_value, "short_distance_frac", &short_distance_frac) != 0) {
		add_example(audit, "heldout diag parse failed");
		++audit->bad_rows;
		fclose(fp);
		return 1;
	}
	if (fgets(value_line, sizeof(value_line), fp) != 0) {
		add_example(audit, "heldout diag has more than one data row");
		++audit->bad_rows;
		failed = 1;
	}
	fclose(fp);
	if (enabled != info->heldout_enabled ||
		!check_close(fraction, info->heldout_fraction) ||
		seed != info->heldout_seed ||
		n_input_raw_total != info->heldout_n_input_raw_total ||
		n_raw_train != info->heldout_n_raw_train ||
		n_raw_heldout != info->heldout_n_raw_heldout ||
		n_bpair_heldout != info->heldout_n_bpair_heldout ||
		n_bpair_eval != info->heldout_n_bpair_eval ||
		n_raw_eval != info->heldout_n_raw_eval ||
		!check_close(mean_expected_energy, info->heldout_mean_expected_energy) ||
		!check_close(mean_min_energy, info->heldout_mean_min_energy) ||
		!check_close(mean_entropy, info->heldout_mean_entropy) ||
		!check_close(mean_pU, info->heldout_mean_pU) ||
		!check_close(mean_best_normalized_distance, info->heldout_mean_best_normalized_distance) ||
		!check_close(short_distance_frac, info->heldout_short_distance_frac)) {
		add_example(audit, "heldout diag disagrees with manifest");
		mark_bad_row(audit, &row_bad);
	}
	if ((enabled != 0 && enabled != 1) || !isfinite(fraction) || fraction < 0.0 ||
		fraction >= 1.0 || n_input_raw_total < 0 || n_raw_train < 0 ||
		n_raw_heldout < 0 || n_bpair_heldout < 0 || n_bpair_eval < 0 ||
		n_raw_eval < 0 || n_bpair_same_bin_skipped < 0 ||
		n_raw_same_bin_skipped < 0 ||
		n_raw_train + n_raw_heldout != n_input_raw_total ||
		n_bpair_eval + n_bpair_same_bin_skipped > n_bpair_heldout ||
		n_raw_eval + n_raw_same_bin_skipped > n_raw_heldout ||
		!isfinite(mean_expected_energy) || !isfinite(mean_min_energy) ||
		!isfinite(mean_entropy) || !isfinite(mean_pU) ||
		!isfinite(mean_best_normalized_distance) || !isfinite(short_distance_frac) ||
		mean_entropy < 0.0 || mean_pU < 0.0 || mean_pU > 1.0 ||
		short_distance_frac < 0.0 || short_distance_frac > 1.0) {
		add_example(audit, "heldout diag fields are out of range");
		mark_bad_row(audit, &row_bad);
	}
	if (row_bad || audit->bad_rows > 0 || audit->forbidden_hits > 0) failed = 1;
	return failed;
}

static void print_examples(const char *label, const struct file_audit *audit)
{
	int i;
	for (i = 0; i < audit->n_examples; ++i)
		fprintf(stderr, "%s example: %s\n", label, audit->examples[i]);
}

static int require_file(const char *label, const char *path)
{
	if (!file_exists(path)) {
		fprintf(stderr, "%s missing: %s\n", label, path);
		return 1;
	}
	if (file_size_or_negative(path) <= 0) {
		fprintf(stderr, "%s empty: %s\n", label, path);
		return 1;
	}
	return 0;
}

int main(int argc, char **argv)
{
	struct manifest_info info;
	struct file_audit manifest_audit, bpair_audit, coords_audit, diag_audit, raw_audit, heldout_audit;
	char manifest_path[1024], bpair_path[1024], coords_path[1024], diag_path[1024];
	char raw_path[1024], heldout_path[1024];
	const char *outdir;
	int failed = 0;
	int audit_failed = 0;
	int audit_heldout = 0;
	int64_t raw_rows = 0, raw_bad_rows = 0, raw_forbidden_hits = 0;
	int64_t forbidden_hits;

	memset(&info, 0, sizeof(info));
	audit_init(&manifest_audit);
	audit_init(&bpair_audit);
	audit_init(&coords_audit);
	audit_init(&diag_audit);
	audit_init(&raw_audit);
	audit_init(&heldout_audit);

	if (argc != 2 || argv[1][0] == 0) {
		fprintf(stderr, "usage: %s OUTDIR\n", argv[0]);
		return 1;
	}
	outdir = argv[1];
	path_join(manifest_path, sizeof(manifest_path), outdir, "p9016_full.manifest.tsv");
	path_join(bpair_path, sizeof(bpair_path), outdir, "p9016_full.bpair_posterior.tsv");
	path_join(coords_path, sizeof(coords_path), outdir, "p9016_full.coords.tsv");
	path_join(diag_path, sizeof(diag_path), outdir, "p9016_full.loop_diag.tsv");
	path_join(raw_path, sizeof(raw_path), outdir, "p9016_full.raw_posterior.tsv");
	path_join(heldout_path, sizeof(heldout_path), outdir, "p9016_full.heldout_diag.tsv");

	failed |= require_file("manifest", manifest_path);
	failed |= require_file("bpair posterior", bpair_path);
	failed |= require_file("coords", coords_path);
	failed |= require_file("loop diag", diag_path);
	if (failed) return 1;

	audit_failed |= audit_manifest(manifest_path, &info, &manifest_audit);
	if (audit_failed) {
		goto summary;
	}
	audit_heldout = manifest_has_heldout_fields(&info) || file_exists(heldout_path);
	if (audit_heldout)
		failed |= require_file("heldout diag", heldout_path);
	if (info.write_raw_posterior == 1)
		failed |= require_file("raw posterior", raw_path);
	if (failed) {
		audit_failed = 1;
		goto summary;
	}

	audit_failed |= audit_bpair_file(bpair_path, &info, &bpair_audit);
	audit_failed |= audit_coords_file(coords_path, &info, &coords_audit);
	audit_failed |= audit_loop_diag_file(diag_path, &info, &diag_audit);
	if (audit_heldout)
		audit_failed |= audit_heldout_diag_file(heldout_path, &info, &heldout_audit);
	if (info.write_raw_posterior == 1) {
		audit_failed |= audit_raw_file(raw_path, &info, &raw_audit);
		raw_rows = raw_audit.rows;
		raw_bad_rows = raw_audit.bad_rows;
		raw_forbidden_hits = raw_audit.forbidden_hits;
	} else {
		audit_init(&raw_audit);
	}

summary:
	forbidden_hits = manifest_audit.forbidden_hits + bpair_audit.forbidden_hits +
		coords_audit.forbidden_hits + diag_audit.forbidden_hits +
		heldout_audit.forbidden_hits + raw_forbidden_hits;
	if (audit_failed) {
		print_examples("manifest", &manifest_audit);
		print_examples("bpair", &bpair_audit);
		print_examples("coords", &coords_audit);
		print_examples("loop_diag", &diag_audit);
		print_examples("heldout_diag", &heldout_audit);
		print_examples("raw", &raw_audit);
	}
	fprintf(stderr,
			"audit status: %s n_raw=%lld n_bpair=%lld n_beads=%lld "
			"bpair_rows=%lld coord_rows=%lld raw_rows=%lld "
			"posterior_bad_rows=%lld coord_bad_rows=%lld raw_bad_rows=%lld "
			"forbidden_string_hits=%lld\n",
			audit_failed? "FAIL" : "OK",
			(long long)info.n_raw, (long long)info.n_bpair, (long long)info.n_beads,
			(long long)bpair_audit.rows, (long long)coords_audit.rows, (long long)raw_rows,
			(long long)bpair_audit.bad_rows, (long long)coords_audit.bad_rows,
			(long long)raw_bad_rows, (long long)forbidden_hits);
	return audit_failed? 1 : 0;
}
