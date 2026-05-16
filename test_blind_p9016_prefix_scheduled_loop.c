#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "hickit.h"

#define HK_BLIND_P9016_SCHEDULED_PREFIX_N 1000
#define HK_BLIND_P9016_SCHEDULED_PATH "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_SCHEDULED_N_ITER 3
#define HK_BLIND_P9016_SCHEDULED_UNIT 1.0f
#define HK_BLIND_P9016_SCHEDULED_D_SCALE 1.0f
#define HK_BLIND_P9016_SCHEDULED_BASE_K 2.0f
#define HK_BLIND_P9016_SCHEDULED_TEMPERATURE_START 2.0f
#define HK_BLIND_P9016_SCHEDULED_TEMPERATURE_END 1.0f
#define HK_BLIND_P9016_SCHEDULED_RHO_TRAIN_START 1.0f
#define HK_BLIND_P9016_SCHEDULED_RHO_TRAIN_END 1.0f
#define HK_BLIND_P9016_SCHEDULED_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_SCHEDULED_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_SCHEDULED_STEP 0.001f
#define HK_BLIND_P9016_SCHEDULED_RELAX_STEPS 5

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

static int file_exists(const char *path)
{
	FILE *fp = fopen(path, "rb");
	if (fp == 0) return 0;
	fclose(fp);
	return 1;
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

static void copy_coords(fvec3_t *dst, const fvec3_t *src, int32_t n)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			dst[i][a] = src[i][a];
}

static int coords_any_changed(const fvec3_t *before, const fvec3_t *after, int32_t n, float tol)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			if (fabsf(after[i][a] - before[i][a]) > tol)
				return 1;
	return 0;
}

static int check_coords_finite(const fvec3_t *coords, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("scheduled-loop coord finite", isfinite(coords[i][a]));
	return failed;
}

static int check_binned_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set, int32_t max_raw)
{
	int failed = 0;
	int32_t i;

	failed |= check_true("scheduled-loop set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_true("scheduled-loop raw positive", set->n_raw > 0);
	failed |= check_true("scheduled-loop raw bounded", set->n_raw <= max_raw);
	failed |= check_true("scheduled-loop bpair positive", set->n_bpairs > 0);
	for (i = 0; i < set->n_raw; ++i) {
		failed |= check_true("scheduled-loop raw2binned id lower", set->raw2binned[i].bpair_id >= 0);
		failed |= check_true("scheduled-loop raw2binned id upper", set->raw2binned[i].bpair_id < set->n_bpairs);
		failed |= check_true("scheduled-loop raw2binned swapped", set->raw2binned[i].swapped == 0 || set->raw2binned[i].swapped == 1);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		failed |= check_true("scheduled-loop key sorted", bp->key.bid[0] <= bp->key.bid[1]);
		failed |= check_true("scheduled-loop key bid0 valid", bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		failed |= check_true("scheduled-loop key bid1 valid", bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
	}
	return failed;
}

static void set_schedule_conf(struct hk_blind_iter_schedule_conf *conf)
{
	conf->n_iter = HK_BLIND_P9016_SCHEDULED_N_ITER;
	conf->base_conf.unit = HK_BLIND_P9016_SCHEDULED_UNIT;
	conf->base_conf.d_scale = HK_BLIND_P9016_SCHEDULED_D_SCALE;
	conf->base_conf.base_k = HK_BLIND_P9016_SCHEDULED_BASE_K;
	conf->base_conf.temperature = HK_BLIND_P9016_SCHEDULED_TEMPERATURE_START;
	conf->base_conf.rho_train = HK_BLIND_P9016_SCHEDULED_RHO_TRAIN_START;
	conf->base_conf.min_sep_unit = HK_BLIND_P9016_SCHEDULED_MIN_SEP_UNIT;
	conf->base_conf.lambda_sep = HK_BLIND_P9016_SCHEDULED_LAMBDA_SEP;
	conf->base_conf.relax_step = HK_BLIND_P9016_SCHEDULED_STEP;
	conf->base_conf.relax_steps = HK_BLIND_P9016_SCHEDULED_RELAX_STEPS;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->base_conf.d_scale_eps_count = 1e-6f;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = HK_BLIND_P9016_SCHEDULED_TEMPERATURE_START;
	conf->temperature_end = HK_BLIND_P9016_SCHEDULED_TEMPERATURE_END;
	conf->rho_train_start = HK_BLIND_P9016_SCHEDULED_RHO_TRAIN_START;
	conf->rho_train_end = HK_BLIND_P9016_SCHEDULED_RHO_TRAIN_END;
}

static int check_clean_pre_relax(const char *label, const struct hk_blind_iter_diag *diag)
{
	int failed = 0;
	char buf[160];

	snprintf(buf, sizeof(buf), "%s posterior nonfinite", label);
	failed |= check_i32(buf, diag->n_posterior_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s posterior bad sum", label);
	failed |= check_i32(buf, diag->n_posterior_bad_sum, 0);
	snprintf(buf, sizeof(buf), "%s posterior out of range", label);
	failed |= check_i32(buf, diag->n_posterior_out_of_range, 0);
	snprintf(buf, sizeof(buf), "%s uncertainty nonfinite", label);
	failed |= check_i32(buf, diag->n_uncertainty_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s uncertainty out of range", label);
	failed |= check_i32(buf, diag->n_uncertainty_out_of_range, 0);
	snprintf(buf, sizeof(buf), "%s five-state bad sum", label);
	failed |= check_i32(buf, diag->n_five_state_bad_sum, 0);
	snprintf(buf, sizeof(buf), "%s wedge nonfinite", label);
	failed |= check_i32(buf, diag->n_wedge_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s wedge bad k", label);
	failed |= check_i32(buf, diag->n_wedge_bad_k, 0);
	snprintf(buf, sizeof(buf), "%s wedge bad d_scale", label);
	failed |= check_i32(buf, diag->n_wedge_bad_d_scale, 0);
	snprintf(buf, sizeof(buf), "%s sep force nonfinite", label);
	failed |= check_i32(buf, diag->sep_force_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s mean entropy finite", label);
	failed |= check_true(buf, isfinite(diag->mean_entropy));
	snprintf(buf, sizeof(buf), "%s mean entropy range", label);
	failed |= check_true(buf, diag->mean_entropy >= -1e-6f && diag->mean_entropy <= logf((float)HK_BLIND_N_STATE) + 1e-6f);
	snprintf(buf, sizeof(buf), "%s mean pU finite", label);
	failed |= check_true(buf, isfinite(diag->mean_pU));
	snprintf(buf, sizeof(buf), "%s mean pU range", label);
	failed |= check_true(buf, diag->mean_pU >= -1e-6f && diag->mean_pU <= 1.0f + 1e-6f);
	return failed;
}

static int check_relax_diag(const char *label, const struct hk_blind_relax_diag *diag)
{
	int failed = 0;
	char buf[160];

	snprintf(buf, sizeof(buf), "%s relax completed", label);
	failed |= check_i32(buf, diag->n_completed, HK_BLIND_P9016_SCHEDULED_RELAX_STEPS);
	snprintf(buf, sizeof(buf), "%s relax nonfinite step", label);
	failed |= check_i32(buf, diag->n_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s relax coord nonfinite", label);
	failed |= check_i32(buf, diag->n_coord_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s backbone nonfinite step", label);
	failed |= check_i32(buf, diag->n_backbone_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s repulsion nonfinite step", label);
	failed |= check_i32(buf, diag->n_repulsion_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s total finite", label);
	failed |= check_true(buf, isfinite(diag->final_total_energy));
	snprintf(buf, sizeof(buf), "%s contact finite", label);
	failed |= check_true(buf, isfinite(diag->final_contact_energy));
	snprintf(buf, sizeof(buf), "%s backbone finite", label);
	failed |= check_true(buf, isfinite(diag->final_backbone_energy));
	snprintf(buf, sizeof(buf), "%s repulsion finite", label);
	failed |= check_true(buf, isfinite(diag->final_repulsion_energy));
	snprintf(buf, sizeof(buf), "%s sep finite", label);
	failed |= check_true(buf, isfinite(diag->final_sep_energy));
	snprintf(buf, sizeof(buf), "%s force finite", label);
	failed |= check_true(buf, isfinite(diag->final_force_l1));
	return failed;
}

static int check_gauge_diag(const char *label, const struct hk_blind_gauge_stats *stats)
{
	int failed = 0;
	char buf[160];

	snprintf(buf, sizeof(buf), "%s gauge no-flip finite", label);
	failed |= check_true(buf, isfinite(stats->mismatch_no_flip));
	snprintf(buf, sizeof(buf), "%s gauge flip finite", label);
	failed |= check_true(buf, isfinite(stats->mismatch_flip));
	snprintf(buf, sizeof(buf), "%s gauge chosen finite", label);
	failed |= check_true(buf, isfinite(stats->mismatch_chosen));
	snprintf(buf, sizeof(buf), "%s gauge flipped nonnegative", label);
	failed |= check_true(buf, stats->n_flipped >= 0);
	return failed;
}

static int check_loop_diag(const struct hk_blind_iter_loop_diag *diag)
{
	int failed = 0;

	failed |= check_i32("scheduled-loop n_iter", diag->n_iter, HK_BLIND_P9016_SCHEDULED_N_ITER);
	failed |= check_i32("scheduled-loop completed", diag->n_completed, HK_BLIND_P9016_SCHEDULED_N_ITER);
	failed |= check_i32("scheduled-loop bad iter", diag->n_bad_iter, 0);
	failed |= check_i32("scheduled-loop relax nonfinite iter", diag->n_relax_nonfinite_iter, 0);
	failed |= check_i32("scheduled-loop coord nonfinite", diag->n_coord_nonfinite, 0);
	failed |= check_close("scheduled-loop initial temperature", diag->initial_temperature, HK_BLIND_P9016_SCHEDULED_TEMPERATURE_START);
	failed |= check_close("scheduled-loop final temperature", diag->final_temperature, HK_BLIND_P9016_SCHEDULED_TEMPERATURE_END);
	failed |= check_close("scheduled-loop initial rho", diag->initial_rho_train, HK_BLIND_P9016_SCHEDULED_RHO_TRAIN_START);
	failed |= check_close("scheduled-loop final rho", diag->final_rho_train, HK_BLIND_P9016_SCHEDULED_RHO_TRAIN_END);
	failed |= check_true("scheduled-loop final entropy finite", isfinite(diag->final_mean_entropy));
	failed |= check_true("scheduled-loop final entropy range",
						 diag->final_mean_entropy >= -1e-6f &&
						 diag->final_mean_entropy <= logf((float)HK_BLIND_N_STATE) + 1e-6f);
	failed |= check_true("scheduled-loop final pU finite", isfinite(diag->final_mean_pU));
	failed |= check_true("scheduled-loop final pU range", diag->final_mean_pU >= -1e-6f && diag->final_mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("scheduled-loop final sep finite", isfinite(diag->final_mean_sep));
	failed |= check_true("scheduled-loop final min sep finite", isfinite(diag->final_min_sep));
	failed |= check_true("scheduled-loop final max sep finite", isfinite(diag->final_max_sep));
	failed |= check_true("scheduled-loop final sum k finite", isfinite(diag->final_sum_wedge_k));
	failed |= check_true("scheduled-loop final sum k nonnegative", diag->final_sum_wedge_k >= 0.0);
	failed |= check_true("scheduled-loop total chr flipped nonnegative", diag->total_chr_flipped >= 0);
	return failed;
}

int main(void)
{
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	struct hk_blind_bpair_set *set = 0;
	struct hk_fdg_conf conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	struct hk_blind_single_iter_diag *per_iter = 0;
	fvec3_t *haploid = 0, *diploid = 0, *diploid_before = 0;
	int32_t n_prefix, n_diploid;
	int32_t i;
	float max_force_l1 = 0.0f;
	int moved = 0;
	int failed = 0;

	if (!file_exists(HK_BLIND_P9016_SCHEDULED_PATH)) {
		fprintf(stderr, "SKIP: %s not found\n", HK_BLIND_P9016_SCHEDULED_PATH);
		return 0;
	}

	hk_verbose = 0;
	m = hk_map_read(HK_BLIND_P9016_SCHEDULED_PATH);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", HK_BLIND_P9016_SCHEDULED_PATH);
		return 1;
	}
	n_prefix = m->n_pairs < HK_BLIND_P9016_SCHEDULED_PREFIX_N? m->n_pairs : HK_BLIND_P9016_SCHEDULED_PREFIX_N;
	if (n_prefix <= 0) {
		fprintf(stderr, "no pairs in %s\n", HK_BLIND_P9016_SCHEDULED_PATH);
		failed = 1;
		goto cleanup;
	}

	bmap = hk_bmap_gen(m->d, n_prefix, m->pairs, 1000000, 1);
	raw = (struct hk_blind_pair*)calloc(n_prefix, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate P9016 scheduled-loop bmap/raw data\n");
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n_prefix; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);

	set = hk_blind_bpair_set_build(bmap, n_prefix, raw);
	failed |= check_binned_set(bmap, set, HK_BLIND_P9016_SCHEDULED_PREFIX_N);
	if (failed) goto cleanup;

	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	haploid = (fvec3_t*)calloc(bmap->n_beads, sizeof(*haploid));
	diploid = (fvec3_t*)calloc(n_diploid, sizeof(*diploid));
	diploid_before = (fvec3_t*)calloc(n_diploid, sizeof(*diploid_before));
	per_iter = (struct hk_blind_single_iter_diag*)calloc(HK_BLIND_P9016_SCHEDULED_N_ITER, sizeof(*per_iter));
	if (haploid == 0 || diploid == 0 || diploid_before == 0 || per_iter == 0) {
		fprintf(stderr, "failed to allocate P9016 scheduled-loop data\n");
		failed = 1;
		goto cleanup;
	}
	set_haploid_scaffold(bmap, haploid);
	failed |= check_i32("scheduled-loop init diploid",
						hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
																  diploid, 0.5f, 0.0f, 17), 0);
	failed |= check_coords_finite(diploid, n_diploid);
	if (failed) goto cleanup;
	copy_coords(diploid_before, diploid, n_diploid);

	hk_fdg_conf_init(&conf);
	set_schedule_conf(&schedule_conf);
	failed |= check_i32("scheduled-loop ret",
						hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &conf, diploid, 0,
															 &schedule_conf, per_iter, &loop_diag), 0);
	failed |= check_loop_diag(&loop_diag);
	for (i = 0; i < HK_BLIND_P9016_SCHEDULED_N_ITER; ++i) {
		char label[64];
		snprintf(label, sizeof(label), "scheduled-loop iter%d", (int)i);
		failed |= check_clean_pre_relax(label, &per_iter[i].pre_relax_diag);
		failed |= check_relax_diag(label, &per_iter[i].relax_diag);
		failed |= check_gauge_diag(label, &per_iter[i].gauge_stats);
		failed |= check_true("scheduled-loop per-iter wedges nonnegative", per_iter[i].n_wedges >= 0);
		failed |= check_true("scheduled-loop per-iter sum k finite", isfinite(per_iter[i].sum_wedge_k));
		failed |= check_true("scheduled-loop per-iter sum k nonnegative", per_iter[i].sum_wedge_k >= 0.0);
		if (per_iter[i].relax_diag.max_force_l1 > max_force_l1)
			max_force_l1 = per_iter[i].relax_diag.max_force_l1;
	}
	failed |= check_coords_finite(diploid, n_diploid);
	moved = coords_any_changed(diploid_before, diploid, n_diploid, 1e-8f);
	if (max_force_l1 > 1e-12f)
		failed |= check_true("scheduled-loop coordinate movement", moved);
	else
		fprintf(stderr, "NOTE: P9016 scheduled-loop force_l1 is zero; no movement expected\n");

	if (!failed) {
		fprintf(stderr,
				"P9016 scheduled-loop: n_raw=%d n_bpair=%d n_iter=%d T=%.3g->%.3g "
				"final_entropy=%.8g final_pU=%.8g final_sum_wedge_k=%.8g total_chr_flipped=%d moved=%d\n",
				set->n_raw, set->n_bpairs, loop_diag.n_completed,
				loop_diag.initial_temperature, loop_diag.final_temperature,
				loop_diag.final_mean_entropy, loop_diag.final_mean_pU, loop_diag.final_sum_wedge_k,
				loop_diag.total_chr_flipped, moved);
	}

cleanup:
	free(haploid);
	free(diploid);
	free(diploid_before);
	free(per_iter);
	hk_blind_bpair_set_destroy(set);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
