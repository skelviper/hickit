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

static int check_i64(const char *label, int64_t got, int64_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %lld, expected %lld\n", label, (long long)got, (long long)expected);
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

static int check_contains(const char *label, const char *s, const char *needle)
{
	if (strstr(s, needle) == 0) {
		fprintf(stderr, "%s: missing substring '%s'\n", label, needle);
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

static float dist3(const fvec3_t x, const fvec3_t y)
{
	float dx = x[0] - y[0];
	float dy = x[1] - y[1];
	float dz = x[2] - y[2];
	return sqrtf(dx * dx + dy * dy + dz * dz);
}

static int coords_changed(const fvec3_t *a, const fvec3_t *b, int32_t n)
{
	int32_t i;
	int k;
	for (i = 0; i < n; ++i)
		for (k = 0; k < 3; ++k)
			if (fabsf(a[i][k] - b[i][k]) > 1e-6f)
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
			failed |= check_true("finite coord", isfinite(coords[i][a]));
	return failed;
}

static int p4_close(const float a[HK_BLIND_N_STATE], const float b[HK_BLIND_N_STATE])
{
	int s;
	for (s = 0; s < HK_BLIND_N_STATE; ++s)
		if (fabsf(a[s] - b[s]) > 1e-5f)
			return 0;
	return 1;
}

static int check_p4_close(const char *label, const float got[HK_BLIND_N_STATE], const float expected[HK_BLIND_N_STATE])
{
	int failed = 0;
	int s;
	char buf[96];
	for (s = 0; s < HK_BLIND_N_STATE; ++s) {
		snprintf(buf, sizeof(buf), "%s state%d", label, s);
		failed |= check_close(buf, got[s], expected[s]);
	}
	return failed;
}

static int check_p4_differs(const char *label, const float a[HK_BLIND_N_STATE], const float b[HK_BLIND_N_STATE])
{
	if (p4_close(a, b)) {
		fprintf(stderr, "%s: p4 unexpectedly matched\n", label);
		return 1;
	}
	return 0;
}

static void set_bmap(struct hk_bmap *bmap, struct hk_bead *beads, int32_t n_beads, int two_chr)
{
	int32_t i;
	memset(bmap, 0, sizeof(*bmap));
	bmap->n_beads = n_beads;
	bmap->unit = 1.0f;
	bmap->beads = beads;
	for (i = 0; i < n_beads; ++i) {
		beads[i].chr = two_chr && i >= n_beads / 2? 1 : 0;
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

static void set_one_bpair_set(struct hk_blind_bpair_set *set, struct hk_blind_bpair *bp,
							  struct hk_blind_raw2binned *raw2binned)
{
	init_bpair(bp, 0, 1);
	hk_blind_raw2binned_set(raw2binned, 0, 0);
	set->bpairs = bp;
	set->n_bpairs = 1;
	set->raw2binned = raw2binned;
	set->n_raw = 1;
}

static void set_empty_bpair_set(struct hk_blind_bpair_set *set)
{
	set->bpairs = 0;
	set->n_bpairs = 0;
	set->raw2binned = 0;
	set->n_raw = 0;
}

static void set_iter_conf(struct hk_blind_single_iter_conf *conf, int32_t relax_steps)
{
	conf->unit = 1.0f;
	conf->d_scale = 1.0f;
	conf->base_k = 2.0f;
	conf->temperature = 1.0f;
	conf->rho_train = 1.0f;
	conf->min_sep_unit = 0.0f;
	conf->lambda_sep = 0.0f;
	conf->relax_step = 0.01f;
	conf->relax_steps = relax_steps;
	conf->enable_repulsion = 0;
	conf->repulsion_mode = HK_BLIND_REPULSION_NONE;
	conf->rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->d_scale_eps_count = 1e-6f;
	conf->rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->contact_k_multiplier_cis = 1.0f;
	conf->contact_k_multiplier_trans = 1.0f;
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
			snprintf(buf, sizeof(buf), "p4 finite %d/%d", (int)i, s);
			failed |= check_true(buf, isfinite(bp->p4[s]));
			snprintf(buf, sizeof(buf), "p4 range %d/%d", (int)i, s);
			failed |= check_true(buf, bp->p4[s] >= -1e-6f && bp->p4[s] <= 1.0f + 1e-6f);
			sum += bp->p4[s];
		}
		snprintf(buf, sizeof(buf), "p4 sum %d", (int)i);
		failed |= check_close(buf, sum, 1.0);
		snprintf(buf, sizeof(buf), "five-state sum %d", (int)i);
		failed |= check_close(buf, bp->rho_output * sum + bp->pU, 1.0);
	}
	return failed;
}

static int check_clean_pre_diag(const char *label, const struct hk_blind_iter_diag *diag)
{
	int failed = 0;
	char buf[128];
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
	snprintf(buf, sizeof(buf), "%s sep nonfinite", label);
	failed |= check_i32(buf, diag->sep_stats.n_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s sep force nonfinite", label);
	failed |= check_i32(buf, diag->sep_force_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s wedge nonfinite", label);
	failed |= check_i32(buf, diag->n_wedge_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s wedge bad k", label);
	failed |= check_i32(buf, diag->n_wedge_bad_k, 0);
	snprintf(buf, sizeof(buf), "%s wedge bad d_scale", label);
	failed |= check_i32(buf, diag->n_wedge_bad_d_scale, 0);
	return failed;
}

static int check_relax_diag(const char *label, const struct hk_blind_relax_diag *diag, int32_t n_steps)
{
	int failed = 0;
	char buf[128];
	snprintf(buf, sizeof(buf), "%s steps", label);
	failed |= check_i32(buf, diag->n_steps, n_steps);
	snprintf(buf, sizeof(buf), "%s completed", label);
	failed |= check_i32(buf, diag->n_completed, n_steps);
	snprintf(buf, sizeof(buf), "%s nonfinite step", label);
	failed |= check_i32(buf, diag->n_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s backbone nonfinite", label);
	failed |= check_i32(buf, diag->n_backbone_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s repulsion nonfinite", label);
	failed |= check_i32(buf, diag->n_repulsion_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s coord nonfinite", label);
	failed |= check_i32(buf, diag->n_coord_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s initial total finite", label);
	failed |= check_true(buf, isfinite(diag->initial_total_energy));
	snprintf(buf, sizeof(buf), "%s final total finite", label);
	failed |= check_true(buf, isfinite(diag->final_total_energy));
	snprintf(buf, sizeof(buf), "%s total sum initial", label);
	failed |= check_close(buf, diag->initial_total_energy,
						  diag->initial_contact_energy + diag->initial_backbone_energy +
						  diag->initial_repulsion_energy + diag->initial_sep_energy);
	snprintf(buf, sizeof(buf), "%s total sum final", label);
	failed |= check_close(buf, diag->final_total_energy,
						  diag->final_contact_energy + diag->final_backbone_energy +
						  diag->final_repulsion_energy + diag->final_sep_energy);
	return failed;
}

static float compute_contact_energy_with_params(const struct hk_fdg_conf *conf, const struct hk_blind_bpair_set *set,
												const fvec3_t *coords, int32_t n_diploid)
{
	struct hk_blind_wedge_list list;
	fvec3_t force[4] = {{0}};
	float energy;
	int32_t n_nonfinite = 0;

	hk_blind_wedge_list_init(&list);
	hk_blind_wedge_list_build_from_bpair_set_params(&list, set, 1.0f);
	hk_blind_wedge_list_aggregate_exact(&list);
	energy = hk_blind_wedge_list_accumulate_contact_force_cpu(conf, &list, coords, n_diploid, 1.0f, force, &n_nonfinite);
	hk_blind_wedge_list_destroy(&list);
	return n_nonfinite == 0? energy : NAN;
}

static void compute_posterior_params(const struct hk_fdg_conf *conf, const struct hk_blind_bpair *bp, const fvec3_t *coords,
									 float p4[HK_BLIND_N_STATE])
{
	float log_prior[HK_BLIND_N_STATE];
	float entropy, pmax, margin, rho_output, pU;

	hk_blind_init_uniform_log_prior(log_prior);
	hk_blind_bpair_posterior_from_coords(conf, bp->key.bid[0], bp->key.bid[1], coords,
										 1.0f, bp->base_d_scale, bp->base_k, log_prior, 1.0f,
										 0, p4, &entropy, &pmax, &margin, &rho_output, &pU);
}

static void compute_posterior_global(const struct hk_fdg_conf *conf, const struct hk_blind_bpair *bp, const fvec3_t *coords,
									 float d_scale, float k, float p4[HK_BLIND_N_STATE])
{
	float log_prior[HK_BLIND_N_STATE];
	float entropy, pmax, margin, rho_output, pU;

	hk_blind_init_uniform_log_prior(log_prior);
	hk_blind_bpair_posterior_from_coords(conf, bp->key.bid[0], bp->key.bid[1], coords,
										 1.0f, d_scale, k, log_prior, 1.0f,
										 0, p4, &entropy, &pmax, &margin, &rho_output, &pU);
}

static int check_per_bpair_posterior_update_uses_contact_params(void)
{
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf fdg_conf;
	float log_prior[HK_BLIND_N_STATE];
	fvec3_t coords[4];
	float global0[HK_BLIND_N_STATE], global1[HK_BLIND_N_STATE];
	float param0[HK_BLIND_N_STATE], param1[HK_BLIND_N_STATE];
	int failed = 0;

	init_bpair(&bpairs[0], 0, 1);
	init_bpair(&bpairs[1], 0, 1);
	bpairs[0].base_d_scale = 1.0f;
	bpairs[1].base_d_scale = 0.5f;
	set.bpairs = bpairs;
	set.n_bpairs = 2;
	set.raw2binned = 0;
	set.n_raw = 2;

	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 2.5f, 0.0f, 0.0f);
	set_coord(coords[3], 2.5f, 1.0f, 0.0f);
	hk_fdg_conf_init(&fdg_conf);
	hk_blind_init_uniform_log_prior(log_prior);

	hk_blind_bpair_set_update_posterior_from_coords(&set, &fdg_conf, coords, 1.0f, 1.0f, 1.0f, log_prior, 1.0f);
	memcpy(global0, bpairs[0].p4, sizeof(global0));
	memcpy(global1, bpairs[1].p4, sizeof(global1));
	failed |= check_p4_close("global helper ignores bp0 params", global0, global1);

	hk_blind_bpair_set_update_posterior_from_coords_params(&set, &fdg_conf, coords, 1.0f, log_prior, 1.0f);
	memcpy(param0, bpairs[0].p4, sizeof(param0));
	memcpy(param1, bpairs[1].p4, sizeof(param1));
	failed |= check_p4_differs("param helper uses base_d_scale", param0, param1);

	bpairs[0].base_d_scale = 1.0f;
	bpairs[1].base_d_scale = 1.0f;
	bpairs[0].base_k = 1.0f;
	bpairs[1].base_k = 3.0f;
	hk_blind_bpair_set_update_posterior_from_coords(&set, &fdg_conf, coords, 1.0f, 1.0f, 1.0f, log_prior, 1.0f);
	memcpy(global0, bpairs[0].p4, sizeof(global0));
	memcpy(global1, bpairs[1].p4, sizeof(global1));
	failed |= check_p4_close("global helper ignores bp base_k", global0, global1);

	hk_blind_bpair_set_update_posterior_from_coords_params(&set, &fdg_conf, coords, 1.0f, log_prior, 1.0f);
	memcpy(param0, bpairs[0].p4, sizeof(param0));
	memcpy(param1, bpairs[1].p4, sizeof(param1));
	failed |= check_p4_differs("param helper uses base_k", param0, param1);
	failed |= check_posterior_valid(&set);
	return failed;
}

static int check_basic_single_iteration(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bp;
	struct hk_blind_raw2binned raw2binned;
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_single_iter_conf iter_conf;
	struct hk_blind_single_iter_diag diag;
	fvec3_t coords[4], before[4];
	char summary[2048];
	int failed = 0;

	set_bmap(&bmap, beads, 2, 0);
	set_one_bpair_set(&set, &bp, &raw2binned);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 2.5f, 0.0f, 0.0f);
	set_coord(coords[3], 2.5f, 1.0f, 0.0f);
	copy_coords(before, coords, 4);
	hk_fdg_conf_init(&fdg_conf);
	set_iter_conf(&iter_conf, 5);

	failed |= check_i32("single iter ret", hk_blind_run_single_iter_cpu(&bmap, &set, &fdg_conf, coords,
																		0, &iter_conf, &diag), 0);
	failed |= check_i32("single n_raw", diag.pre_relax_diag.n_raw, 1);
	failed |= check_i32("single n_bpair", diag.pre_relax_diag.n_bpair, 1);
	failed |= check_true("single coords changed", coords_changed(coords, before, 4));
	failed |= check_coords_finite(coords, 4);
	failed |= check_posterior_valid(&set);
	failed |= check_clean_pre_diag("single", &diag.pre_relax_diag);
	failed |= check_relax_diag("single", &diag.relax_diag, 5);
	failed |= check_true("single wedges", diag.n_wedges > 0);
	failed |= check_i32("single top wedges", diag.n_wedges, diag.pre_relax_diag.n_wedges);
	failed |= check_i32("single wedges before", diag.n_wedges_before_aggregation, diag.pre_relax_diag.n_wedges_before_aggregation);
	failed |= check_i64("single skipped self", diag.n_skipped_self_edges, diag.pre_relax_diag.n_skipped_self_edges);
	failed |= check_close("single sum k", diag.sum_wedge_k, diag.pre_relax_diag.sum_wedge_k);
	failed |= check_true("single backbone energy", diag.relax_diag.initial_backbone_energy > 0.0f);
	failed |= check_true("single backbone force", diag.relax_diag.max_backbone_force_l1 > 0.0f);
	failed |= check_i32("single gauge flips", diag.gauge_stats.n_flipped, 0);
	failed |= check_i32("single n chr flipped", diag.n_chr_flipped, 0);
	hk_blind_iter_diag_snprintf(summary, sizeof(summary), &diag.pre_relax_diag);
	failed |= check_contains("single summary", summary, "status: OK");
	return failed;
}

static int check_single_iteration_uses_bpair_contact_params(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bp;
	struct hk_blind_raw2binned raw2binned;
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_single_iter_conf iter_conf;
	struct hk_blind_single_iter_diag diag;
	fvec3_t coords[4], before[4];
	float expected_contact_energy;
	float expected_p4[HK_BLIND_N_STATE], global_p4[HK_BLIND_N_STATE];
	int failed = 0;

	set_bmap(&bmap, beads, 2, 0);
	set_one_bpair_set(&set, &bp, &raw2binned);
	bp.base_d_scale = 0.5f;
	bp.base_k = 1.0f;
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 2.5f, 0.0f, 0.0f);
	set_coord(coords[3], 2.5f, 1.0f, 0.0f);
	copy_coords(before, coords, 4);
	hk_fdg_conf_init(&fdg_conf);
	set_iter_conf(&iter_conf, 1);
	iter_conf.base_k = 7.0f;
	iter_conf.d_scale = 3.0f;
	iter_conf.relax_step = 0.0f;

	compute_posterior_params(&fdg_conf, &bp, before, expected_p4);
	compute_posterior_global(&fdg_conf, &bp, before, iter_conf.d_scale, iter_conf.base_k, global_p4);
	failed |= check_i32("single params ret", hk_blind_run_single_iter_cpu(&bmap, &set, &fdg_conf, coords,
																		 0, &iter_conf, &diag), 0);
	expected_contact_energy = compute_contact_energy_with_params(&fdg_conf, &set, before, 4);
	failed |= check_coords_close("single params coords", coords, before, 4);
	failed |= check_p4_close("single params p4", bp.p4, expected_p4);
	failed |= check_p4_differs("single params not global p4", bp.p4, global_p4);
	failed |= check_close("single params sum k", diag.sum_wedge_k, 1.0);
	failed |= check_close("single params contact energy", diag.relax_diag.initial_contact_energy, expected_contact_energy);
	failed |= check_true("single params contact differs from zero", diag.relax_diag.initial_contact_energy > 0.0f);
	failed |= check_clean_pre_diag("single params", &diag.pre_relax_diag);
	failed |= check_relax_diag("single params", &diag.relax_diag, 1);
	return failed;
}

static int check_zero_relax_steps(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bp;
	struct hk_blind_raw2binned raw2binned;
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_single_iter_conf iter_conf;
	struct hk_blind_single_iter_diag diag;
	fvec3_t coords[4], before[4];
	int failed = 0;

	set_bmap(&bmap, beads, 2, 0);
	set_one_bpair_set(&set, &bp, &raw2binned);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 2.5f, 0.0f, 0.0f);
	set_coord(coords[3], 2.5f, 1.0f, 0.0f);
	copy_coords(before, coords, 4);
	hk_fdg_conf_init(&fdg_conf);
	set_iter_conf(&iter_conf, 0);

	failed |= check_i32("zero single ret", hk_blind_run_single_iter_cpu(&bmap, &set, &fdg_conf, coords,
																		0, &iter_conf, &diag), 0);
	failed |= check_coords_close("zero coords", coords, before, 4);
	failed |= check_posterior_valid(&set);
	failed |= check_clean_pre_diag("zero", &diag.pre_relax_diag);
	failed |= check_relax_diag("zero", &diag.relax_diag, 0);
	failed |= check_i32("zero completed", diag.relax_diag.n_completed, 0);
	failed |= check_true("zero wedges built", diag.n_wedges > 0);
	failed |= check_i32("zero gauge flips", diag.gauge_stats.n_flipped, 0);
	return failed;
}

static int check_backbone_only_single_iteration(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_single_iter_conf iter_conf;
	struct hk_blind_single_iter_diag diag;
	fvec3_t coords[4];
	float before, after;
	int failed = 0;

	set_bmap(&bmap, beads, 2, 0);
	set_empty_bpair_set(&set);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 2.5f, 0.0f, 0.0f);
	set_coord(coords[3], 2.5f, 1.0f, 0.0f);
	before = dist3(coords[0], coords[2]);
	hk_fdg_conf_init(&fdg_conf);
	set_iter_conf(&iter_conf, 4);

	failed |= check_i32("backbone only ret", hk_blind_run_single_iter_cpu(&bmap, &set, &fdg_conf, coords,
																		 0, &iter_conf, &diag), 0);
	after = dist3(coords[0], coords[2]);
	failed |= check_true("backbone only distance decreases", after < before);
	failed |= check_i32("backbone only bpair", diag.pre_relax_diag.n_bpair, 0);
	failed |= check_i32("backbone only raw", diag.pre_relax_diag.n_raw, 0);
	failed |= check_i32("backbone only wedges", diag.n_wedges, 0);
	failed |= check_i32("backbone only wedges pre", diag.pre_relax_diag.n_wedges, 0);
	failed |= check_true("backbone only energy", diag.relax_diag.initial_backbone_energy > 0.0f);
	failed |= check_true("backbone only force", diag.relax_diag.max_backbone_force_l1 > 0.0f);
	failed |= check_relax_diag("backbone only", &diag.relax_diag, 4);
	failed |= check_clean_pre_diag("backbone only", &diag.pre_relax_diag);
	failed |= check_coords_finite(coords, 4);
	return failed;
}

static void swap_copies(fvec3_t *coords, int32_t haploid_bid)
{
	int32_t b0 = hk_diploid_bid(haploid_bid, HK_DIPLOID_COPY0);
	int32_t b1 = hk_diploid_bid(haploid_bid, HK_DIPLOID_COPY1);
	int a;
	for (a = 0; a < 3; ++a) {
		float t = coords[b0][a];
		coords[b0][a] = coords[b1][a];
		coords[b1][a] = t;
	}
}

static int check_gauge_lower_level_reference(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	struct hk_blind_gauge_stats stats;
	fvec3_t prev[8], cur[8];
	int failed = 0;

	set_bmap(&bmap, beads, 4, 1);
	set_coord(prev[0], 0.0f, 0.0f, 0.0f);
	set_coord(prev[1], 0.0f, 1.0f, 0.0f);
	set_coord(prev[2], 1.0f, 0.0f, 0.0f);
	set_coord(prev[3], 1.0f, 1.0f, 0.0f);
	set_coord(prev[4], 10.0f, 0.0f, 0.0f);
	set_coord(prev[5], 10.0f, 1.0f, 0.0f);
	set_coord(prev[6], 11.0f, 0.0f, 0.0f);
	set_coord(prev[7], 11.0f, 1.0f, 0.0f);
	copy_coords(cur, prev, 8);
	swap_copies(cur, 0);
	swap_copies(cur, 1);

	failed |= check_i32("gauge ref ret", hk_blind_temporal_gauge_stabilize_bmap(&bmap, prev, cur, 0, &stats), 0);
	failed |= check_i32("gauge ref flipped", stats.n_flipped, 1);
	failed |= check_coords_close("gauge ref coords", cur, prev, 8);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_per_bpair_posterior_update_uses_contact_params();
	failed |= check_basic_single_iteration();
	failed |= check_single_iteration_uses_bpair_contact_params();
	failed |= check_zero_relax_steps();
	failed |= check_backbone_only_single_iteration();
	failed |= check_gauge_lower_level_reference();
	return failed != 0;
}
