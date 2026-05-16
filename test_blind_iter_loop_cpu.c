#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "hickit.h"

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

static void set_coord(fvec3_t x, float x0, float x1, float x2)
{
	x[0] = x0;
	x[1] = x1;
	x[2] = x2;
}

static void copy_coords(fvec3_t *dst, const fvec3_t *src, int32_t n)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			dst[i][a] = src[i][a];
}

static int coords_changed(const fvec3_t *a, const fvec3_t *b, int32_t n)
{
	int32_t i;
	int axis;
	for (i = 0; i < n; ++i)
		for (axis = 0; axis < 3; ++axis)
			if (fabsf(a[i][axis] - b[i][axis]) > 1e-6f)
				return 1;
	return 0;
}

static int check_coords_close(const char *label, const fvec3_t *got, const fvec3_t *expected, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	char buf[96];
	for (i = 0; i < n; ++i) {
		for (a = 0; a < 3; ++a) {
			snprintf(buf, sizeof(buf), "%s bead%d axis%d", label, (int)i, a);
			failed |= check_close(buf, got[i][a], expected[i][a]);
		}
	}
	return failed;
}

static int check_coords_finite(const fvec3_t *coords, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("coord finite", isfinite(coords[i][a]));
	return failed;
}

static void set_bmap(struct hk_bmap *bmap, struct hk_bead *beads, int32_t n_beads)
{
	int32_t i;
	memset(bmap, 0, sizeof(*bmap));
	bmap->n_beads = n_beads;
	bmap->unit = 1.0f;
	bmap->beads = beads;
	for (i = 0; i < n_beads; ++i) {
		beads[i].chr = 0;
		beads[i].st = i * 1000000;
		beads[i].en = beads[i].st + 1000000;
	}
}

static void init_bpair(struct hk_blind_bpair *bp, int32_t bid0, int32_t bid1)
{
	int s;
	bp->key.bid[0] = bid0;
	bp->key.bid[1] = bid1;
	bp->n_raw = 1;
	bp->base_d_scale = 1.0f;
	bp->base_k = 1.0f;
	for (s = 0; s < HK_BLIND_N_STATE; ++s) {
		bp->log_prior[s] = -logf((float)HK_BLIND_N_STATE);
		bp->p4[s] = 0.25f;
	}
	bp->entropy = logf((float)HK_BLIND_N_STATE);
	bp->pmax = 0.25f;
	bp->margin = 0.0f;
	bp->rho_output = 0.0f;
	bp->pU = 1.0f;
}

static void set_bpair_set(struct hk_blind_bpair_set *set, struct hk_blind_bpair *bpairs,
						  struct hk_blind_raw2binned *raw2binned)
{
	init_bpair(&bpairs[0], 0, 1);
	init_bpair(&bpairs[1], 1, 1);
	hk_blind_raw2binned_set(&raw2binned[0], 0, 0);
	hk_blind_raw2binned_set(&raw2binned[1], 1, 0);
	set->bpairs = bpairs;
	set->n_bpairs = 2;
	set->raw2binned = raw2binned;
	set->n_raw = 2;
}

static void set_coords(fvec3_t coords[4])
{
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 2.5f, 0.0f, 0.0f);
	set_coord(coords[3], 2.5f, 1.0f, 0.0f);
}

static void set_loop_conf(struct hk_blind_iter_loop_conf *conf, int32_t n_iter)
{
	conf->n_iter = n_iter;
	conf->single_iter_conf.unit = 1.0f;
	conf->single_iter_conf.d_scale = 1.0f;
	conf->single_iter_conf.base_k = 2.0f;
	conf->single_iter_conf.temperature = 1.0f;
	conf->single_iter_conf.rho_train = 1.0f;
	conf->single_iter_conf.min_sep_unit = 0.0f;
	conf->single_iter_conf.lambda_sep = 0.0f;
	conf->single_iter_conf.relax_step = 0.01f;
	conf->single_iter_conf.relax_steps = 2;
	conf->single_iter_conf.enable_repulsion = 0;
	conf->single_iter_conf.repulsion_mode = HK_BLIND_REPULSION_NONE;
	conf->single_iter_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->single_iter_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->single_iter_conf.d_scale_eps_count = 1e-6f;
	conf->single_iter_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->single_iter_conf.contact_k_multiplier_cis = 1.0f;
	conf->single_iter_conf.contact_k_multiplier_trans = 1.0f;
}

static void set_schedule_conf(struct hk_blind_iter_schedule_conf *conf, int32_t n_iter,
							  float temperature_start, float temperature_end,
							  float rho_train_start, float rho_train_end)
{
	conf->n_iter = n_iter;
	conf->base_conf.unit = 1.0f;
	conf->base_conf.d_scale = 1.0f;
	conf->base_conf.base_k = 2.0f;
	conf->base_conf.temperature = temperature_start;
	conf->base_conf.rho_train = rho_train_start;
	conf->base_conf.min_sep_unit = 0.0f;
	conf->base_conf.lambda_sep = 0.0f;
	conf->base_conf.relax_step = 0.01f;
	conf->base_conf.relax_steps = 2;
	conf->base_conf.enable_repulsion = 0;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_NONE;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->base_conf.d_scale_eps_count = 1e-6f;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = temperature_start;
	conf->temperature_end = temperature_end;
	conf->rho_train_start = rho_train_start;
	conf->rho_train_end = rho_train_end;
}

static int check_posterior_valid(const struct hk_blind_bpair_set *set)
{
	int failed = 0;
	int32_t i;
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		double sum = 0.0;
		int s;
		char buf[96];
		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			snprintf(buf, sizeof(buf), "posterior finite %d/%d", (int)i, s);
			failed |= check_true(buf, isfinite(bp->p4[s]));
			snprintf(buf, sizeof(buf), "posterior range %d/%d", (int)i, s);
			failed |= check_true(buf, bp->p4[s] >= -1e-6f && bp->p4[s] <= 1.0f + 1e-6f);
			sum += bp->p4[s];
		}
		snprintf(buf, sizeof(buf), "posterior sum %d", (int)i);
		failed |= check_close(buf, sum, 1.0);
		snprintf(buf, sizeof(buf), "five-state sum %d", (int)i);
		failed |= check_close(buf, bp->rho_output * sum + bp->pU, 1.0);
	}
	return failed;
}

static int check_final_posterior_matches_coords(const char *label, const struct hk_blind_bpair_set *set,
												const struct hk_fdg_conf *fdg_conf, const fvec3_t *coords,
												float unit, float temperature)
{
	struct hk_blind_bpair expected_bpairs[2];
	struct hk_blind_raw2binned expected_raw2binned[2];
	struct hk_blind_bpair_set expected;
	int failed = 0;
	int32_t i;
	int s;
	char buf[128];

	failed |= check_i32("posterior match fixture size", set->n_bpairs, 2);
	if (failed) return failed;
	memcpy(expected_bpairs, set->bpairs, sizeof(expected_bpairs));
	memcpy(expected_raw2binned, set->raw2binned, sizeof(expected_raw2binned));
	expected = *set;
	expected.bpairs = expected_bpairs;
	expected.raw2binned = expected_raw2binned;
	hk_blind_bpair_set_update_posterior_from_coords_params(&expected, fdg_conf, coords, unit, 0, temperature);
	for (i = 0; i < set->n_bpairs; ++i) {
		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			snprintf(buf, sizeof(buf), "%s p4 bpair%d state%d", label, (int)i, s);
			failed |= check_close(buf, set->bpairs[i].p4[s], expected.bpairs[i].p4[s]);
		}
		snprintf(buf, sizeof(buf), "%s pU bpair%d", label, (int)i);
		failed |= check_close(buf, set->bpairs[i].pU, expected.bpairs[i].pU);
		snprintf(buf, sizeof(buf), "%s rho bpair%d", label, (int)i);
		failed |= check_close(buf, set->bpairs[i].rho_output, expected.bpairs[i].rho_output);
	}
	return failed;
}

static int check_clean_single_iter_diag(const char *label, const struct hk_blind_single_iter_diag *diag, int32_t relax_steps)
{
	int failed = 0;
	char buf[128];
	snprintf(buf, sizeof(buf), "%s posterior nonfinite", label);
	failed |= check_i32(buf, diag->pre_relax_diag.n_posterior_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s posterior bad sum", label);
	failed |= check_i32(buf, diag->pre_relax_diag.n_posterior_bad_sum, 0);
	snprintf(buf, sizeof(buf), "%s posterior out of range", label);
	failed |= check_i32(buf, diag->pre_relax_diag.n_posterior_out_of_range, 0);
	snprintf(buf, sizeof(buf), "%s uncertainty nonfinite", label);
	failed |= check_i32(buf, diag->pre_relax_diag.n_uncertainty_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s uncertainty out of range", label);
	failed |= check_i32(buf, diag->pre_relax_diag.n_uncertainty_out_of_range, 0);
	snprintf(buf, sizeof(buf), "%s five-state bad sum", label);
	failed |= check_i32(buf, diag->pre_relax_diag.n_five_state_bad_sum, 0);
	snprintf(buf, sizeof(buf), "%s wedge nonfinite", label);
	failed |= check_i32(buf, diag->pre_relax_diag.n_wedge_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s wedge bad k", label);
	failed |= check_i32(buf, diag->pre_relax_diag.n_wedge_bad_k, 0);
	snprintf(buf, sizeof(buf), "%s wedge bad d_scale", label);
	failed |= check_i32(buf, diag->pre_relax_diag.n_wedge_bad_d_scale, 0);
	snprintf(buf, sizeof(buf), "%s sep force nonfinite", label);
	failed |= check_i32(buf, diag->pre_relax_diag.sep_force_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s relax completed", label);
	failed |= check_i32(buf, diag->relax_diag.n_completed, relax_steps);
	snprintf(buf, sizeof(buf), "%s relax nonfinite", label);
	failed |= check_i32(buf, diag->relax_diag.n_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s relax coord nonfinite", label);
	failed |= check_i32(buf, diag->relax_diag.n_coord_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s backbone nonfinite", label);
	failed |= check_i32(buf, diag->relax_diag.n_backbone_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s repulsion nonfinite", label);
	failed |= check_i32(buf, diag->relax_diag.n_repulsion_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s mean entropy finite", label);
	failed |= check_true(buf, isfinite(diag->pre_relax_diag.mean_entropy));
	snprintf(buf, sizeof(buf), "%s mean pU finite", label);
	failed |= check_true(buf, isfinite(diag->pre_relax_diag.mean_pU));
	snprintf(buf, sizeof(buf), "%s backbone finite", label);
	failed |= check_true(buf, isfinite(diag->relax_diag.final_backbone_energy));
	snprintf(buf, sizeof(buf), "%s repulsion finite", label);
	failed |= check_true(buf, isfinite(diag->relax_diag.final_repulsion_energy));
	snprintf(buf, sizeof(buf), "%s backbone active", label);
	failed |= check_true(buf, diag->relax_diag.max_backbone_force_l1 > 0.0f);
	return failed;
}

static int check_loop_diag_clean(const struct hk_blind_iter_loop_diag *diag, int32_t n_iter)
{
	int failed = 0;
	failed |= check_i32("loop n_iter", diag->n_iter, n_iter);
	failed |= check_i32("loop completed", diag->n_completed, n_iter);
	failed |= check_i32("loop bad iter", diag->n_bad_iter, 0);
	failed |= check_i32("loop relax bad iter", diag->n_relax_nonfinite_iter, 0);
	failed |= check_i32("loop coord nonfinite", diag->n_coord_nonfinite, 0);
	failed |= check_true("loop initial entropy finite", isfinite(diag->initial_mean_entropy));
	failed |= check_true("loop final entropy finite", isfinite(diag->final_mean_entropy));
	failed |= check_true("loop initial pU finite", isfinite(diag->initial_mean_pU));
	failed |= check_true("loop final pU finite", isfinite(diag->final_mean_pU));
	failed |= check_true("loop entropy range", diag->final_mean_entropy >= -1e-6f && diag->final_mean_entropy <= logf((float)HK_BLIND_N_STATE) + 1e-6f);
	failed |= check_true("loop pU range", diag->final_mean_pU >= -1e-6f && diag->final_mean_pU <= 1.0f + 1e-6f);
	if (n_iter > 0) {
		failed |= check_true("loop initial temperature finite", isfinite(diag->initial_temperature));
		failed |= check_true("loop final temperature finite", isfinite(diag->final_temperature));
		failed |= check_true("loop initial rho finite", isfinite(diag->initial_rho_train));
		failed |= check_true("loop final rho finite", isfinite(diag->final_rho_train));
	}
	failed |= check_true("loop chr flipped nonnegative", diag->total_chr_flipped >= 0);
	failed |= check_true("loop sep mean finite", isfinite(diag->final_mean_sep));
	failed |= check_true("loop sep min finite", isfinite(diag->final_min_sep));
	failed |= check_true("loop sep max finite", isfinite(diag->final_max_sep));
	failed |= check_true("loop sum k finite", isfinite(diag->final_sum_wedge_k));
	failed |= check_true("loop sum k nonnegative", diag->final_sum_wedge_k >= 0.0);
	return failed;
}

static int check_zero_iter(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_raw2binned raw2binned[2];
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_loop_conf loop_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	fvec3_t coords[4], before[4];
	int failed = 0;

	set_bmap(&bmap, beads, 2);
	set_bpair_set(&set, bpairs, raw2binned);
	set_coords(coords);
	copy_coords(before, coords, 4);
	hk_fdg_conf_init(&fdg_conf);
	set_loop_conf(&loop_conf, 0);

	failed |= check_i32("zero loop ret", hk_blind_run_iter_loop_cpu(&bmap, &set, &fdg_conf, coords,
																	0, &loop_conf, 0, &loop_diag), 0);
	failed |= check_coords_close("zero loop coords", coords, before, 4);
	failed |= check_i32("zero loop completed", loop_diag.n_completed, 0);
	failed |= check_i32("zero loop coord nonfinite", loop_diag.n_coord_nonfinite, 0);
	failed |= check_close("zero loop initial entropy", loop_diag.initial_mean_entropy, 0.0);
	failed |= check_true("zero loop sep finite", isfinite(loop_diag.final_mean_sep));
	return failed;
}

static int check_fixed_k_loop(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_raw2binned raw2binned[2];
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_loop_conf loop_conf;
	struct hk_blind_single_iter_diag per_iter[3];
	struct hk_blind_iter_loop_diag loop_diag;
	fvec3_t coords[4], before[4];
	int failed = 0;
	int i;

	set_bmap(&bmap, beads, 2);
	set_bpair_set(&set, bpairs, raw2binned);
	set_coords(coords);
	copy_coords(before, coords, 4);
	hk_fdg_conf_init(&fdg_conf);
	set_loop_conf(&loop_conf, 3);

	failed |= check_i32("fixed loop ret", hk_blind_run_iter_loop_cpu(&bmap, &set, &fdg_conf, coords,
																	 0, &loop_conf, per_iter, &loop_diag), 0);
	failed |= check_loop_diag_clean(&loop_diag, 3);
	failed |= check_true("fixed loop coords changed", coords_changed(coords, before, 4));
	failed |= check_coords_finite(coords, 4);
	failed |= check_posterior_valid(&set);
	for (i = 0; i < 3; ++i) {
		char buf[64];
		snprintf(buf, sizeof(buf), "iter%d", i);
		failed |= check_clean_single_iter_diag(buf, &per_iter[i], loop_conf.single_iter_conf.relax_steps);
	}
	failed |= check_close("fixed initial entropy", loop_diag.initial_mean_entropy, per_iter[0].pre_relax_diag.mean_entropy);
	failed |= check_close("fixed initial pU", loop_diag.initial_mean_pU, per_iter[0].pre_relax_diag.mean_pU);
	failed |= check_close("fixed final sum k", loop_diag.final_sum_wedge_k, per_iter[2].sum_wedge_k);
	failed |= check_i32("fixed posterior refreshed", loop_diag.posterior_refreshed_after_final_relax, 1);
	failed |= check_true("fixed refresh kl finite", isfinite(loop_diag.posterior_refresh_mean_kl));
	failed |= check_true("fixed refresh changed posterior",
						  loop_diag.posterior_refresh_mean_kl > 1e-12 ||
						  loop_diag.posterior_refresh_top_state_switch_frac > 0.0f ||
						  fabsf(loop_diag.posterior_refresh_mean_pU_before -
								 loop_diag.posterior_refresh_mean_pU_after) > 1e-7f);
	failed |= check_final_posterior_matches_coords("fixed final refreshed", &set, &fdg_conf, coords,
												   loop_conf.single_iter_conf.unit,
												   loop_conf.single_iter_conf.temperature);
	failed |= check_close("fixed initial temperature", loop_diag.initial_temperature, loop_conf.single_iter_conf.temperature);
	failed |= check_close("fixed final temperature", loop_diag.final_temperature, loop_conf.single_iter_conf.temperature);
	failed |= check_close("fixed initial rho", loop_diag.initial_rho_train, loop_conf.single_iter_conf.rho_train);
	failed |= check_close("fixed final rho", loop_diag.final_rho_train, loop_conf.single_iter_conf.rho_train);
	return failed;
}

static int check_loop_without_per_iter_buffer(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_raw2binned raw2binned[2];
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_loop_conf loop_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	fvec3_t coords[4];
	int failed = 0;

	set_bmap(&bmap, beads, 2);
	set_bpair_set(&set, bpairs, raw2binned);
	set_coords(coords);
	hk_fdg_conf_init(&fdg_conf);
	set_loop_conf(&loop_conf, 2);

	failed |= check_i32("null per-iter ret", hk_blind_run_iter_loop_cpu(&bmap, &set, &fdg_conf, coords,
																		0, &loop_conf, 0, &loop_diag), 0);
	failed |= check_loop_diag_clean(&loop_diag, 2);
	failed |= check_coords_finite(coords, 4);
	return failed;
}

static int check_schedule_values(void)
{
	struct hk_blind_iter_schedule_conf conf;
	int failed = 0;

	set_schedule_conf(&conf, 3, 2.0f, 1.0f, 1.0f, 0.5f);
	failed |= check_close("schedule temp t0", hk_blind_iter_schedule_temperature_at(&conf, 0), 2.0);
	failed |= check_close("schedule temp t1", hk_blind_iter_schedule_temperature_at(&conf, 1), sqrt(2.0));
	failed |= check_close("schedule temp t2", hk_blind_iter_schedule_temperature_at(&conf, 2), 1.0);
	failed |= check_close("schedule rho t0", hk_blind_iter_schedule_rho_train_at(&conf, 0), 1.0);
	failed |= check_close("schedule rho t1", hk_blind_iter_schedule_rho_train_at(&conf, 1), 0.75);
	failed |= check_close("schedule rho t2", hk_blind_iter_schedule_rho_train_at(&conf, 2), 0.5);

	set_schedule_conf(&conf, 1, 2.0f, 1.0f, 0.25f, 0.75f);
	failed |= check_close("single schedule temp", hk_blind_iter_schedule_temperature_at(&conf, 0), 2.0);
	failed |= check_close("single schedule rho", hk_blind_iter_schedule_rho_train_at(&conf, 0), 0.25);
	return failed;
}

static int check_scheduled_fixed_equivalence(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair fixed_bpairs[2], scheduled_bpairs[2];
	struct hk_blind_raw2binned fixed_raw2binned[2], scheduled_raw2binned[2];
	struct hk_blind_bpair_set fixed_set, scheduled_set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_loop_conf fixed_conf;
	struct hk_blind_iter_schedule_conf scheduled_conf;
	struct hk_blind_iter_loop_diag fixed_diag, scheduled_diag;
	fvec3_t fixed_coords[4], scheduled_coords[4];
	int failed = 0;

	set_bmap(&bmap, beads, 2);
	set_bpair_set(&fixed_set, fixed_bpairs, fixed_raw2binned);
	set_bpair_set(&scheduled_set, scheduled_bpairs, scheduled_raw2binned);
	set_coords(fixed_coords);
	copy_coords(scheduled_coords, fixed_coords, 4);
	hk_fdg_conf_init(&fdg_conf);
	set_loop_conf(&fixed_conf, 2);
	fixed_conf.single_iter_conf.temperature = 1.25f;
	fixed_conf.single_iter_conf.rho_train = 0.75f;
	set_schedule_conf(&scheduled_conf, 2, 1.25f, 1.25f, 0.75f, 0.75f);

	failed |= check_i32("fixed eq ret", hk_blind_run_iter_loop_cpu(&bmap, &fixed_set, &fdg_conf, fixed_coords,
																	0, &fixed_conf, 0, &fixed_diag), 0);
	failed |= check_i32("scheduled eq ret", hk_blind_run_iter_loop_scheduled_cpu(&bmap, &scheduled_set, &fdg_conf,
																				  scheduled_coords, 0, &scheduled_conf,
																				  0, &scheduled_diag), 0);
	failed |= check_loop_diag_clean(&fixed_diag, 2);
	failed |= check_loop_diag_clean(&scheduled_diag, 2);
	failed |= check_coords_close("scheduled fixed coords", scheduled_coords, fixed_coords, 4);
	failed |= check_close("scheduled fixed sum k", scheduled_diag.final_sum_wedge_k, fixed_diag.final_sum_wedge_k);
	failed |= check_close("scheduled fixed final entropy", scheduled_diag.final_mean_entropy, fixed_diag.final_mean_entropy);
	failed |= check_close("scheduled initial temperature", scheduled_diag.initial_temperature, 1.25);
	failed |= check_close("scheduled final temperature", scheduled_diag.final_temperature, 1.25);
	failed |= check_close("scheduled initial rho", scheduled_diag.initial_rho_train, 0.75);
	failed |= check_close("scheduled final rho", scheduled_diag.final_rho_train, 0.75);
	failed |= check_final_posterior_matches_coords("scheduled final refreshed", &scheduled_set, &fdg_conf,
												   scheduled_coords, scheduled_conf.base_conf.unit,
												   scheduled_conf.temperature_end);
	return failed;
}

static int check_rho_train_schedule_sum_k(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair rho1_bpairs[2], rho05_bpairs[2];
	struct hk_blind_raw2binned rho1_raw2binned[2], rho05_raw2binned[2];
	struct hk_blind_bpair_set rho1_set, rho05_set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_schedule_conf rho1_conf, rho05_conf;
	struct hk_blind_iter_loop_diag rho1_diag, rho05_diag;
	fvec3_t rho1_coords[4], rho05_coords[4];
	int failed = 0;

	set_bmap(&bmap, beads, 2);
	set_bpair_set(&rho1_set, rho1_bpairs, rho1_raw2binned);
	set_bpair_set(&rho05_set, rho05_bpairs, rho05_raw2binned);
	set_coords(rho1_coords);
	copy_coords(rho05_coords, rho1_coords, 4);
	hk_fdg_conf_init(&fdg_conf);
	set_schedule_conf(&rho1_conf, 1, 1.0f, 1.0f, 1.0f, 1.0f);
	set_schedule_conf(&rho05_conf, 1, 1.0f, 1.0f, 0.5f, 0.5f);
	rho1_conf.base_conf.relax_steps = 0;
	rho05_conf.base_conf.relax_steps = 0;

	failed |= check_i32("rho1 scheduled ret", hk_blind_run_iter_loop_scheduled_cpu(&bmap, &rho1_set, &fdg_conf,
																				   rho1_coords, 0, &rho1_conf, 0,
																				   &rho1_diag), 0);
	failed |= check_i32("rho05 scheduled ret", hk_blind_run_iter_loop_scheduled_cpu(&bmap, &rho05_set, &fdg_conf,
																					rho05_coords, 0, &rho05_conf, 0,
																					&rho05_diag), 0);
	failed |= check_loop_diag_clean(&rho1_diag, 1);
	failed |= check_loop_diag_clean(&rho05_diag, 1);
	failed |= check_true("rho1 sum positive", rho1_diag.final_sum_wedge_k > 0.0);
	failed |= check_close("rho half sum k", rho05_diag.final_sum_wedge_k, 0.5 * rho1_diag.final_sum_wedge_k);
	failed |= check_close("rho05 final rho", rho05_diag.final_rho_train, 0.5);
	return failed;
}

static int check_scheduled_loop_runs(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_raw2binned raw2binned[2];
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_single_iter_diag per_iter[3];
	struct hk_blind_iter_loop_diag loop_diag;
	fvec3_t coords[4], before[4];
	int failed = 0;
	int i;

	set_bmap(&bmap, beads, 2);
	set_bpair_set(&set, bpairs, raw2binned);
	set_coords(coords);
	copy_coords(before, coords, 4);
	hk_fdg_conf_init(&fdg_conf);
	set_schedule_conf(&schedule_conf, 3, 2.0f, 1.0f, 1.0f, 1.0f);

	failed |= check_i32("scheduled loop ret", hk_blind_run_iter_loop_scheduled_cpu(&bmap, &set, &fdg_conf, coords,
																				   0, &schedule_conf, per_iter,
																				   &loop_diag), 0);
	failed |= check_loop_diag_clean(&loop_diag, 3);
	failed |= check_close("scheduled loop initial temp", loop_diag.initial_temperature, 2.0);
	failed |= check_close("scheduled loop final temp", loop_diag.final_temperature, 1.0);
	failed |= check_close("scheduled loop initial rho", loop_diag.initial_rho_train, 1.0);
	failed |= check_close("scheduled loop final rho", loop_diag.final_rho_train, 1.0);
	failed |= check_true("scheduled loop coords changed", coords_changed(coords, before, 4));
	failed |= check_coords_finite(coords, 4);
	failed |= check_final_posterior_matches_coords("scheduled loop final refreshed", &set, &fdg_conf,
												   coords, schedule_conf.base_conf.unit,
												   schedule_conf.temperature_end);
	for (i = 0; i < 3; ++i) {
		char buf[64];
		snprintf(buf, sizeof(buf), "scheduled iter%d", i);
		failed |= check_clean_single_iter_diag(buf, &per_iter[i], schedule_conf.base_conf.relax_steps);
	}
	return failed;
}

static int run_invalid_temperature(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_raw2binned raw2binned[2];
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	fvec3_t coords[4];

	set_bmap(&bmap, beads, 2);
	set_bpair_set(&set, bpairs, raw2binned);
	set_coords(coords);
	hk_fdg_conf_init(&fdg_conf);
	set_schedule_conf(&schedule_conf, 1, 0.0f, 1.0f, 1.0f, 1.0f);
	return hk_blind_run_iter_loop_scheduled_cpu(&bmap, &set, &fdg_conf, coords, 0, &schedule_conf, 0, &loop_diag);
}

static int run_invalid_rho_train(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_raw2binned raw2binned[2];
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	fvec3_t coords[4];

	set_bmap(&bmap, beads, 2);
	set_bpair_set(&set, bpairs, raw2binned);
	set_coords(coords);
	hk_fdg_conf_init(&fdg_conf);
	set_schedule_conf(&schedule_conf, 1, 1.0f, 1.0f, -0.5f, 1.0f);
	return hk_blind_run_iter_loop_scheduled_cpu(&bmap, &set, &fdg_conf, coords, 0, &schedule_conf, 0, &loop_diag);
}

int main(int argc, char **argv)
{
	int failed = 0;
	if (argc > 1 && strcmp(argv[1], "--invalid-schedule-temperature") == 0)
		return run_invalid_temperature();
	if (argc > 1 && strcmp(argv[1], "--invalid-schedule-rho") == 0)
		return run_invalid_rho_train();
	failed |= check_zero_iter();
	failed |= check_fixed_k_loop();
	failed |= check_loop_without_per_iter_buffer();
	failed |= check_schedule_values();
	failed |= check_scheduled_fixed_equivalence();
	failed |= check_rho_train_schedule_sum_k();
	failed |= check_scheduled_loop_runs();
	return failed != 0;
}
