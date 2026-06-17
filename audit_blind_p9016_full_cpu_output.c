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
	MK_SAME_BIN_FILTER_ENABLED,
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
	MK_INIT_EPS,
	MK_INIT_NOISE_SCALE,
	MK_INIT_SCALE,
	MK_CHR_SEP_UNIT,
	MK_LAMBDA_CHR_SEP,
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
	int64_t n_raw_same_bin_excluded;
	int64_t n_bpair_same_bin_excluded;
	int64_t n_raw_cis;
	int64_t n_raw_trans;
	int64_t n_bpair_cis;
	int64_t n_bpair_trans;
	int64_t uses_phase_labels;
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
	double min_sep_unit;
	double lambda_sep;
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
	int64_t prior_n_distance_bins;
	int64_t prior_n_alpha;
	char base_k_mode[64];
	char init_mode[128];
	char prior_mode[128];
	char prior_smoothing_method[128];
	char rho_train_mode[128];
	char d_scale_mode[128];
	char repulsion_blocking_mode[128];
	char raw_posterior_same_bin_policy[128];
	char posterior_refresh_prior_mode[128];
	char baseline[64];
	char config_name[128];
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
	if (strcmp(key, "same_bin_filter_enabled") == 0) return MK_SAME_BIN_FILTER_ENABLED;
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
	if (strcmp(key, "init_eps") == 0) return MK_INIT_EPS;
	if (strcmp(key, "init_noise_scale") == 0) return MK_INIT_NOISE_SCALE;
	if (strcmp(key, "init_scale") == 0) return MK_INIT_SCALE;
	if (strcmp(key, "chr_sep_unit") == 0) return MK_CHR_SEP_UNIT;
	if (strcmp(key, "lambda_chr_sep") == 0) return MK_LAMBDA_CHR_SEP;
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
	return strcmp(key, "input_contact_source") == 0 ||
		   strcmp(key, "resolution_label") == 0 ||
		   strcmp(key, "refinement_stage") == 0 ||
		   strcmp(key, "refinement_stage_index") == 0 ||
		   strcmp(key, "refinement_n_stages") == 0 ||
		   strcmp(key, "parent_bin_size_bp") == 0 ||
		   strcmp(key, "refinement_init_source") == 0 ||
		   strcmp(key, "refinement_child_offset_step") == 0 ||
			   strcmp(key, "refinement_parent_anchor_k") == 0 ||
			   strcmp(key, "relax_backend") == 0 ||
			   strcmp(key, "approved_p9016_raw_pairs_realpath") == 0 ||
			   strcmp(key, "scaffold_source") == 0 ||
		   strcmp(key, "scaffold_fdg_n_iter") == 0 ||
			   strcmp(key, "training_graph_weighted_filter") == 0 ||
			   strcmp(key, "training_graph_probability_weighted") == 0 ||
				   strcmp(key, "edge_k_probability_weighted") == 0 ||
				   strcmp(key, "dscale_mode") == 0 ||
				   strcmp(key, "d_scale_mode_input_string") == 0 ||
				   strcmp(key, "d_scale_posterior_gamma") == 0 ||
				   strcmp(key, "dscale_posterior_gamma") == 0 ||
				   strcmp(key, "dscale_effective_count_formula") == 0 ||
				   strcmp(key, "dscale_probability_weighted") == 0 ||
				   strcmp(key, "legacy_expected_count_alias_used") == 0 ||
				   strcmp(key, "training_graph_dscale_probability_weighted") == 0 ||
				   strcmp(key, "training_graph_mode") == 0 ||
		   strcmp(key, "softall_mode") == 0 ||
		   strcmp(key, "estep_score_mode") == 0 ||
		   strcmp(key, "copy_labels_are_gauge_only") == 0 ||
		   strcmp(key, "uses_charm_or_reference") == 0 ||
		   strcmp(key, "uses_charm_for_training") == 0 ||
		   strcmp(key, "repulsion_multiplier") == 0 ||
		   strcmp(key, "k_rel_rep_effective") == 0 ||
		   strcmp(key, "output_coords_gz") == 0 ||
		   strcmp(key, "output_force_class_diag") == 0 ||
		   strcmp(key, "output_sep_diag") == 0 ||
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
	case MK_PRIOR_SMOOTHING_METHOD: snprintf(info->prior_smoothing_method, sizeof(info->prior_smoothing_method), "%s", value); return 0;
	case MK_RHO_TRAIN_MODE: snprintf(info->rho_train_mode, sizeof(info->rho_train_mode), "%s", value); return 0;
	case MK_D_SCALE_MODE: snprintf(info->d_scale_mode, sizeof(info->d_scale_mode), "%s", value); return 0;
	case MK_REPULSION_BLOCKING_MODE: snprintf(info->repulsion_blocking_mode, sizeof(info->repulsion_blocking_mode), "%s", value); return 0;
	case MK_RAW_POSTERIOR_SAME_BIN_POLICY: snprintf(info->raw_posterior_same_bin_policy, sizeof(info->raw_posterior_same_bin_policy), "%s", value); return 0;
	case MK_POSTERIOR_REFRESH_PRIOR_MODE: snprintf(info->posterior_refresh_prior_mode, sizeof(info->posterior_refresh_prior_mode), "%s", value); return 0;
	case MK_BASELINE: snprintf(info->baseline, sizeof(info->baseline), "%s", value); return 0;
	case MK_MSTEP_GRAPH_MODE: snprintf(info->mstep_graph_mode, sizeof(info->mstep_graph_mode), "%s", value); return 0;
	case MK_CONFIG_NAME: snprintf(info->config_name, sizeof(info->config_name), "%s", value); return 0;
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
	case MK_N_RAW_SAME_BIN_EXCLUDED: return parse_i64_value(value, &info->n_raw_same_bin_excluded);
	case MK_N_BPAIR_SAME_BIN_EXCLUDED: return parse_i64_value(value, &info->n_bpair_same_bin_excluded);
	case MK_N_RAW_CIS: return parse_i64_value(value, &info->n_raw_cis);
	case MK_N_RAW_TRANS: return parse_i64_value(value, &info->n_raw_trans);
	case MK_N_BPAIR_CIS: return parse_i64_value(value, &info->n_bpair_cis);
	case MK_N_BPAIR_TRANS: return parse_i64_value(value, &info->n_bpair_trans);
	case MK_USES_PHASE_LABELS: return parse_i64_value(value, &info->uses_phase_labels);
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
	case MK_MIN_SEP_UNIT: return parse_double_value(value, &info->min_sep_unit);
	case MK_LAMBDA_SEP: return parse_double_value(value, &info->lambda_sep);
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
	if (strcmp(info->rho_train_mode, "constant") != 0) {
		add_example(audit, "manifest rho_train_mode is not recognized");
		failed = 1;
	}
	if (info->seen[MK_RHO_TRAIN] &&
		(!isfinite(info->rho_train) || !check_close(info->rho_train, 1.0))) {
		add_example(audit, "manifest rho_train must remain one");
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
	if (strcmp(info->raw_posterior_same_bin_policy, "uniform_unknown_rows") != 0) {
		add_example(audit, "manifest raw posterior same-bin policy is not uniform_unknown_rows");
		failed = 1;
	}
	if (!isfinite(info->min_sep_unit) || info->min_sep_unit < 0.0 ||
		!isfinite(info->lambda_sep) || info->lambda_sep < 0.0 ||
		!isfinite(info->relax_step) || info->relax_step < 0.0) {
		add_example(audit, "manifest separation/relax value out of range");
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
	if (info->posterior_refreshed_after_final_relax != 1 ||
		!isfinite(info->posterior_refresh_temperature) ||
		info->posterior_refresh_temperature <= 0.0 ||
		strcmp(info->posterior_refresh_prior_mode, info->prior_mode) != 0 ||
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
		"repulsion_mode", "posterior_refreshed_after_final_relax", "posterior_refresh_temperature"
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
		tsv_get_double(header_fields, value_fields, n_header, n_value, "posterior_refresh_mean_pU_after", &refresh_pu_after) != 0) {
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
	struct file_audit manifest_audit, bpair_audit, coords_audit, diag_audit, raw_audit;
	char manifest_path[1024], bpair_path[1024], coords_path[1024], diag_path[1024], raw_path[1024];
	const char *outdir;
	int failed = 0;
	int audit_failed = 0;
	int64_t raw_rows = 0, raw_bad_rows = 0, raw_forbidden_hits = 0;
	int64_t forbidden_hits;

	memset(&info, 0, sizeof(info));
	audit_init(&manifest_audit);
	audit_init(&bpair_audit);
	audit_init(&coords_audit);
	audit_init(&diag_audit);
	audit_init(&raw_audit);

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

	failed |= require_file("manifest", manifest_path);
	failed |= require_file("bpair posterior", bpair_path);
	failed |= require_file("coords", coords_path);
	failed |= require_file("loop diag", diag_path);
	if (failed) return 1;

	audit_failed |= audit_manifest(manifest_path, &info, &manifest_audit);
	if (audit_failed) {
		goto summary;
	}
	if (info.write_raw_posterior == 1)
		failed |= require_file("raw posterior", raw_path);
	if (failed) {
		audit_failed = 1;
		goto summary;
	}

	audit_failed |= audit_bpair_file(bpair_path, &info, &bpair_audit);
	audit_failed |= audit_coords_file(coords_path, &info, &coords_audit);
	audit_failed |= audit_loop_diag_file(diag_path, &info, &diag_audit);
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
		coords_audit.forbidden_hits + diag_audit.forbidden_hits + raw_forbidden_hits;
	if (audit_failed) {
		print_examples("manifest", &manifest_audit);
		print_examples("bpair", &bpair_audit);
		print_examples("coords", &coords_audit);
		print_examples("loop_diag", &diag_audit);
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
