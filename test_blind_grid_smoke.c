#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "hickit.h"

#define GRID_SMOKE_RESOLUTION 1000000

struct grid_smoke_config {
	const char *id;
	int init_mode;
	int prior_mode;
	int rho_train_mode;
	int d_scale_mode;
	uint64_t seed;
};

static const struct grid_smoke_config configs[] = {
	{ "random_uniform_constant_raw", HK_BLIND_INIT_RANDOM_DIPLOID, HK_BLIND_PRIOR_UNIFORM,
	  HK_BLIND_RHO_TRAIN_CONSTANT, HK_BLIND_D_SCALE_RAW_COUNT, 11ULL },
	{ "toy_cis_entropy_expected", HK_BLIND_INIT_TOY_SPLIT, HK_BLIND_PRIOR_CIS_INTER_RATIO,
	  HK_BLIND_RHO_TRAIN_ENTROPY, HK_BLIND_D_SCALE_EXPECTED_COUNT, 12ULL }
};

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

static void path_join(char *dst, size_t dst_size, const char *dir, const char *name)
{
	size_t n = strlen(dir);
	if (n > 0 && dir[n - 1] == '/')
		snprintf(dst, dst_size, "%s%s", dir, name);
	else
		snprintf(dst, dst_size, "%s/%s", dir, name);
}

static int make_root(char *dir, size_t dir_size)
{
	int i;
	for (i = 0; i < 100; ++i) {
		snprintf(dir, dir_size, "/tmp/hk_blind_grid_smoke_%ld_%d", (long)getpid(), i);
		if (mkdir(dir, 0700) == 0)
			return 0;
		if (errno != EEXIST)
			return -1;
	}
	return -1;
}

static void remove_run_dir(const char *dir)
{
	char path[1024];
	if (dir == 0 || dir[0] == 0) return;
	path_join(path, sizeof(path), dir, "p9016_full.manifest.tsv"); remove(path);
	path_join(path, sizeof(path), dir, "p9016_full.bpair_posterior.tsv"); remove(path);
	path_join(path, sizeof(path), dir, "p9016_full.coords.tsv"); remove(path);
	path_join(path, sizeof(path), dir, "p9016_full.loop_diag.tsv"); remove(path);
	rmdir(dir);
}

static void remove_root(const char *root, const char run_dirs[][1024], size_t n_run)
{
	char path[1024];
	size_t i;
	for (i = 0; i < n_run; ++i)
		remove_run_dir(run_dirs[i]);
	path_join(path, sizeof(path), root, "grid_summary.tsv");
	remove(path);
	rmdir(root);
}

static int init_coords(struct hk_bmap *bmap, const struct grid_smoke_config *config,
					   fvec3_t *haploid, fvec3_t *diploid)
{
	if (config->init_mode == HK_BLIND_INIT_RANDOM_DIPLOID)
		return hk_blind_init_random_diploid_coords(bmap, diploid, 1.0f, config->seed);
	if (config->init_mode == HK_BLIND_INIT_TOY_SPLIT) {
		if (hk_blind_init_toy_haploid_scaffold(bmap, haploid) != 0)
			return -1;
		return hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
														 diploid, 0.5f, 0.0f, config->seed);
	}
	return -1;
}

static int init_prior(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
					  const struct grid_smoke_config *config,
					  struct hk_blind_prior_diag *diag)
{
	memset(diag, 0, sizeof(*diag));
	diag->prior_mode = config->prior_mode;
	diag->eps_prior = 1e-6f;
	if (config->prior_mode == HK_BLIND_PRIOR_UNIFORM) {
		hk_blind_bpair_set_init_uniform_prior(set);
		return 0;
	}
	return hk_blind_bpair_set_init_cis_inter_ratio_prior(bmap, set, 1e-6f, diag);
}

static int write_manifest(const char *path, const char *out_dir, const struct grid_smoke_config *config,
						  const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
						  const struct hk_blind_prior_diag *prior_diag,
						  const struct hk_blind_iter_loop_diag *loop_diag)
{
	FILE *fp = fopen(path, "w");
	struct hk_blind_base_k_stats stats;
	if (fp == 0) return -1;
	hk_blind_bpair_set_base_k_stats(set, &stats);
	if (fprintf(fp,
				"key\tvalue\n"
				"sample\tP9016\n"
				"runner_family\ttest_fixture\n"
				"runner_version\t2026-04-30\n"
				"default_profile\tp9016_grid_smoke_fixture_v1\n"
				"input_path\ttestdata/p9016_blind_fixture.pairs\n"
				"output_dir\t%s\n"
				"n_raw\t%d\n"
				"n_bpair\t%d\n"
				"n_beads\t%d\n"
				"resolution\t%d\n"
				"n_iter\t1\n"
				"unit\t1\n"
				"d_scale\t1\n"
				"base_k_mode\tuniform\n"
				"base_k_effective\t%.9g\n"
				"base_k_min\t%.9g\n"
				"base_k_mean\t%.9g\n"
				"base_k_max\t%.9g\n"
				"base_k_n_nonfinite\t%d\n"
				"legacy_base_k_unused\t2\n"
				"init_mode\t%s\n"
				"init_scale\t1\n"
				"init_eps_effective\t%.9g\n"
				"init_noise_scale_effective\t%.9g\n"
				"init_split_params_used\t%d\n"
				"init_seed\t%llu\n"
				"scaffold_source\t%s\n"
				"scaffold_fdg_n_iter\t0\n"
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
				"rho_train_mode\t%s\n"
				"d_scale_mode\t%s\n"
				"d_scale_eps_count\t1e-06\n"
				"same_bin_filter_enabled\t%d\n"
				"n_raw_same_bin_excluded\t%lld\n"
				"n_bpair_same_bin_excluded\t%d\n"
				"raw_posterior_same_bin_policy\tuniform_unknown_rows\n"
				"min_sep_unit\t0.25\n"
				"lambda_sep\t0.05\n"
				"relax_step\t0.001\n"
				"relax_steps\t1\n"
				"enable_repulsion\t1\n"
				"repulsion_mode\t2\n"
				"repulsion_blocking_mode\tcurrent_edge_blocking\n"
				"temperature_start\t2\n"
				"temperature_end\t2\n"
				"rho_train_start\t1\n"
				"rho_train_end\t1\n"
				"write_raw_posterior\t0\n"
				"output_bpair_posterior\t%s/p9016_full.bpair_posterior.tsv\n"
				"output_coords\t%s/p9016_full.coords.tsv\n"
				"output_loop_diag\t%s/p9016_full.loop_diag.tsv\n"
				"final_mean_entropy\t%.9g\n"
				"final_mean_pU\t%.9g\n"
				"final_sum_wedge_k\t%.17g\n"
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
				"status\tOK\n",
				out_dir, set->n_raw, set->n_bpairs, bmap->n_beads, GRID_SMOKE_RESOLUTION,
				stats.mean, stats.min, stats.mean, stats.max, stats.n_nonfinite,
				hk_blind_init_mode_name(config->init_mode),
				config->init_mode == HK_BLIND_INIT_RANDOM_DIPLOID? 0.0 : 0.5,
				config->init_mode == HK_BLIND_INIT_RANDOM_DIPLOID? 0.0 : 0.0,
				config->init_mode == HK_BLIND_INIT_RANDOM_DIPLOID? 0 : 1,
				(unsigned long long)config->seed,
				config->init_mode == HK_BLIND_INIT_RANDOM_DIPLOID? "random_diploid_direct" : "toy_deterministic",
				hk_blind_prior_mode_name(prior_diag->prior_mode), prior_diag->eps_prior,
				prior_diag->inter_density, prior_diag->observed_inter, prior_diag->possible_inter,
				prior_diag->n_distance_bins, prior_diag->n_alpha, prior_diag->alpha_min,
				prior_diag->alpha_median, prior_diag->alpha_max,
				prior_diag->prior_mode == HK_BLIND_PRIOR_CIS_INTER_RATIO?
					"moving_average_3bin_possible_weighted" : "none",
				prior_diag->eps_prior, 0.5,
				hk_blind_rho_train_mode_name(config->rho_train_mode),
				hk_blind_d_scale_mode_name(config->d_scale_mode),
				set->same_bin_filter_enabled, (long long)set->n_raw_same_bin_excluded,
				set->n_bpair_same_bin_excluded, out_dir, out_dir, out_dir,
				loop_diag->final_mean_entropy, loop_diag->final_mean_pU,
				loop_diag->final_sum_wedge_k, loop_diag->final_repulsion_energy,
				loop_diag->posterior_refreshed_after_final_relax,
				loop_diag->posterior_refresh_temperature,
				hk_blind_prior_mode_name(prior_diag->prior_mode),
				loop_diag->posterior_refresh_mean_kl,
				loop_diag->posterior_refresh_top_state_switch_frac,
				loop_diag->posterior_refresh_mean_pU_before,
				loop_diag->posterior_refresh_mean_pU_after,
				loop_diag->n_bad_iter,
				loop_diag->n_relax_nonfinite_iter, loop_diag->n_coord_nonfinite) < 0) {
		fclose(fp);
		return -1;
	}
	if (fclose(fp) != 0) return -1;
	return 0;
}

static int run_auditor(const char *out_dir)
{
	char cmd[4096];
	int ret;
	snprintf(cmd, sizeof(cmd), "./audit_blind_p9016_full_cpu_output.bin \"%s\"", out_dir);
	ret = system(cmd);
	return ret == 0? 0 : -1;
}

static int run_one(const struct grid_smoke_config *config, const char *root, char *run_dir,
				   size_t run_dir_size, const struct hk_bmap *bmap_template,
				   const struct hk_blind_pair *raw, int32_t n_raw, FILE *summary)
{
	struct hk_blind_bpair_set *set = 0;
	struct hk_blind_prior_diag prior_diag;
	struct hk_blind_iter_schedule_conf loop_conf;
	struct hk_blind_single_iter_diag iter_diag[1];
	struct hk_blind_iter_loop_diag loop_diag;
	struct hk_fdg_conf fdg_conf;
	fvec3_t *haploid = 0, *diploid = 0;
	struct hk_bmap *bmap = hk_bmap_bead_dup(bmap_template);
	int32_t n_diploid;
	int failed = 0;
	char path[1024];

	if (bmap == 0) return 1;
	path_join(run_dir, run_dir_size, root, config->id);
	if (mkdir(run_dir, 0700) != 0) {
		hk_bmap_destroy(bmap);
		return 1;
	}
	set = hk_blind_bpair_set_build(bmap, n_raw, raw);
	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	haploid = (fvec3_t*)calloc((size_t)bmap->n_beads, sizeof(*haploid));
	diploid = (fvec3_t*)calloc((size_t)n_diploid, sizeof(*diploid));
	if (set == 0 || haploid == 0 || diploid == 0) {
		failed = 1;
		goto cleanup;
	}
	failed |= init_prior(bmap, set, config, &prior_diag) != 0;
	failed |= init_coords(bmap, config, haploid, diploid) != 0;
	hk_fdg_conf_init(&fdg_conf);
	fdg_conf.backend = HK_FDG_BACKEND_CPU;
	memset(&loop_conf, 0, sizeof(loop_conf));
	loop_conf.n_iter = 1;
	loop_conf.base_conf.unit = 1.0f;
	loop_conf.base_conf.d_scale = 1.0f;
	loop_conf.base_conf.base_k = 2.0f;
	loop_conf.base_conf.temperature = 2.0f;
	loop_conf.base_conf.rho_train = 1.0f;
	loop_conf.base_conf.min_sep_unit = 0.25f;
	loop_conf.base_conf.lambda_sep = 0.05f;
	loop_conf.base_conf.relax_step = 0.001f;
	loop_conf.base_conf.relax_steps = 1;
	loop_conf.base_conf.enable_repulsion = 1;
	loop_conf.base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	loop_conf.base_conf.rho_train_mode = config->rho_train_mode;
	loop_conf.base_conf.d_scale_mode = config->d_scale_mode;
	loop_conf.base_conf.d_scale_eps_count = 1e-6f;
	loop_conf.base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	loop_conf.base_conf.contact_k_multiplier_cis = 1.0f;
	loop_conf.base_conf.contact_k_multiplier_trans = 1.0f;
	loop_conf.temperature_start = 2.0f;
	loop_conf.temperature_end = 2.0f;
	loop_conf.rho_train_start = 1.0f;
	loop_conf.rho_train_end = 1.0f;
	failed |= hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &fdg_conf, diploid, 0,
												   &loop_conf, iter_diag, &loop_diag) != 0;
	if (failed) goto cleanup;
	path_join(path, sizeof(path), run_dir, "p9016_full.bpair_posterior.tsv");
	{
		FILE *fp = fopen(path, "w");
		if (fp == 0) { failed = 1; goto cleanup; }
		failed |= hk_blind_write_bpair_posterior_tsv(fp, bmap, set) != 0;
		fclose(fp);
	}
	path_join(path, sizeof(path), run_dir, "p9016_full.coords.tsv");
	{
		FILE *fp = fopen(path, "w");
		if (fp == 0) { failed = 1; goto cleanup; }
		failed |= hk_blind_write_diploid_coords_tsv(fp, bmap, diploid) != 0;
		fclose(fp);
	}
	path_join(path, sizeof(path), run_dir, "p9016_full.loop_diag.tsv");
	{
		FILE *fp = fopen(path, "w");
		if (fp == 0) { failed = 1; goto cleanup; }
		failed |= hk_blind_write_iter_loop_diag_tsv(fp, &loop_diag) != 0;
		fclose(fp);
	}
	path_join(path, sizeof(path), run_dir, "p9016_full.manifest.tsv");
	failed |= write_manifest(path, run_dir, config, bmap, set, &prior_diag, &loop_diag) != 0;
	failed |= run_auditor(run_dir) != 0;
	if (!failed) {
		fprintf(summary, "%s\t%s\t%s\t%s\t%s\t%s\t%.9g\t%.9g\tOK\n",
				config->id, run_dir, hk_blind_init_mode_name(config->init_mode),
				hk_blind_prior_mode_name(config->prior_mode),
				hk_blind_rho_train_mode_name(config->rho_train_mode),
				hk_blind_d_scale_mode_name(config->d_scale_mode),
				loop_diag.final_mean_entropy, loop_diag.final_mean_pU);
	}

cleanup:
	free(haploid);
	free(diploid);
	hk_blind_bpair_set_destroy(set);
	hk_bmap_destroy(bmap);
	return failed;
}

static int count_summary_rows(const char *path)
{
	FILE *fp = fopen(path, "r");
	char line[2048];
	int rows = -1;
	if (fp == 0) return -1;
	while (fgets(line, sizeof(line), fp))
		++rows;
	fclose(fp);
	return rows;
}

int main(void)
{
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	char root[512] = {0}, summary_path[1024], run_dirs[2][1024] = {{0}};
	FILE *summary = 0;
	int failed = 0;
	int32_t i;
	size_t c;

	if (!file_exists("./audit_blind_p9016_full_cpu_output.bin")) {
		fprintf(stderr, "grid smoke auditor binary missing\n");
		return 1;
	}
	m = hk_map_read("testdata/p9016_blind_fixture.pairs");
	if (m == 0) return 1;
	bmap = hk_bmap_gen(m->d, m->n_pairs, m->pairs, GRID_SMOKE_RESOLUTION, 1);
	raw = (struct hk_blind_pair*)calloc((size_t)m->n_pairs, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < m->n_pairs; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);
	failed |= check_true("grid root", make_root(root, sizeof(root)) == 0);
	path_join(summary_path, sizeof(summary_path), root, "grid_summary.tsv");
	summary = fopen(summary_path, "w");
	failed |= check_true("grid summary open", summary != 0);
	if (failed) goto cleanup;
	fprintf(summary, "config_id\toutput_dir\tinit_mode\tprior_mode\trho_train_mode\td_scale_mode\tfinal_mean_entropy\tfinal_mean_pU\taudit_status\n");
	for (c = 0; c < sizeof(configs) / sizeof(configs[0]); ++c) {
		failed |= run_one(&configs[c], root, run_dirs[c], sizeof(run_dirs[c]), bmap, raw,
						  m->n_pairs, summary);
	}
	fclose(summary);
	summary = 0;
	failed |= check_true("grid summary row count",
						 count_summary_rows(summary_path) == (int)(sizeof(configs) / sizeof(configs[0])));
	failed |= check_true("grid manifest distinct", strcmp(run_dirs[0], run_dirs[1]) != 0);

cleanup:
	if (summary) fclose(summary);
	remove_root(root, run_dirs, sizeof(configs) / sizeof(configs[0]));
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
