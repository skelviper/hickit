#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "hickit.h"

static int check_i32(const char *label, int32_t got, int32_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %d, expected %d\n", label, got, expected);
		return 1;
	}
	return 0;
}

static int check_close(const char *label, float got, float expected)
{
	float tol = 1e-5f;
	float scale = fabsf(expected) > 1.0f? fabsf(expected) : 1.0f;
	if (fabsf(got - expected) > tol * scale) {
		fprintf(stderr, "%s: got %.8g, expected %.8g\n", label, got, expected);
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

static float dist3(const fvec3_t a, const fvec3_t b)
{
	float dx = a[0] - b[0];
	float dy = a[1] - b[1];
	float dz = a[2] - b[2];
	return sqrtf(dx * dx + dy * dy + dz * dz);
}

static void set_test_bmap(struct hk_bmap *bmap, struct hk_bead beads[5])
{
	int i;
	for (i = 0; i < 5; ++i) {
		beads[i].chr = i < 3? 0 : 1;
		beads[i].st = 1000000 * i;
		beads[i].en = 1000000 * (i + 1);
	}
	bmap->n_beads = 5;
	bmap->n_pairs = 0;
	bmap->unit = 1.0f;
	bmap->d = 0;
	bmap->beads = beads;
	bmap->offcnt = 0;
	bmap->pairs = 0;
	bmap->x = 0;
	bmap->feat = 0;
	bmap->cpg = 0;
	bmap->gc_bias = 0;
	bmap->gc_corrected = 0;
}

static void set_haploid_coords(fvec3_t x[5])
{
	set_coord(x[0], 1.0f, 2.0f, 3.0f);
	set_coord(x[1], 4.0f, 5.0f, 6.0f);
	set_coord(x[2], -1.0f, 7.0f, 0.5f);
	set_coord(x[3], 2.5f, -3.0f, 9.0f);
	set_coord(x[4], -4.0f, 1.5f, -2.0f);
}

static void zero_force(fvec3_t *force, int32_t n_diploid)
{
	int32_t i;
	int a;
	for (i = 0; i < n_diploid; ++i)
		for (a = 0; a < 3; ++a)
			force[i][a] = 0.0f;
}

static void set_sentinel_force(fvec3_t *force, int32_t n_diploid)
{
	int32_t i;
	int a;
	for (i = 0; i < n_diploid; ++i)
		for (a = 0; a < 3; ++a)
			force[i][a] = 10.0f + (float)i + 0.125f * (float)a;
}

static void copy_force(fvec3_t *dst, const fvec3_t *src, int32_t n_diploid)
{
	int32_t i;
	int a;
	for (i = 0; i < n_diploid; ++i)
		for (a = 0; a < 3; ++a)
			dst[i][a] = src[i][a];
}

static int coords_equal(const fvec3_t *a, const fvec3_t *b, int32_t n)
{
	int32_t i;
	int k;
	for (i = 0; i < n; ++i)
		for (k = 0; k < 3; ++k)
			if (a[i][k] != b[i][k])
				return 0;
	return 1;
}

static int check_all_coords_finite(const fvec3_t *x, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("finite initialized coordinate", isfinite(x[i][a]));
	return failed;
}

static int check_initialized_geometry(const fvec3_t *diploid, const fvec3_t *haploid, int32_t n_haploid, float eps)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n_haploid; ++i) {
		int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		for (a = 0; a < 3; ++a)
			failed |= check_close("initialized midpoint", 0.5f * (diploid[b0][a] + diploid[b1][a]), haploid[i][a]);
		failed |= check_close("initialized homolog separation", dist3(diploid[b0], diploid[b1]), 2.0f * eps);
	}
	return failed;
}

static int check_force_finite_nonzero_and_net_zero(const fvec3_t *force, int32_t n_diploid)
{
	fvec3_t net;
	float max_abs = 0.0f;
	int failed = 0;
	int32_t i;
	int a;

	set_coord(net, 0.0f, 0.0f, 0.0f);
	for (i = 0; i < n_diploid; ++i) {
		for (a = 0; a < 3; ++a) {
			failed |= check_true("finite active force", isfinite(force[i][a]));
			net[a] += force[i][a];
			if (fabsf(force[i][a]) > max_abs) max_abs = fabsf(force[i][a]);
		}
	}
	failed |= check_close("active net force x", net[0], 0.0f);
	failed |= check_close("active net force y", net[1], 0.0f);
	failed |= check_close("active net force z", net[2], 0.0f);
	failed |= check_true("active force has nonzero component", max_abs > 0.0f);
	return failed;
}

static int check_force_accumulation(const fvec3_t *sentinel, const fvec3_t *added, const fvec3_t *final_force, int32_t n_diploid)
{
	int failed = 0;
	int32_t i;
	int a;
	char buf[96];

	for (i = 0; i < n_diploid; ++i) {
		for (a = 0; a < 3; ++a) {
			snprintf(buf, sizeof(buf), "sentinel plus added bead%d axis%d", (int)i, a);
			failed |= check_close(buf, final_force[i][a], sentinel[i][a] + added[i][a]);
		}
	}
	return failed;
}

static int check_stats_equal(const struct hk_blind_homolog_sep_stats *a, const struct hk_blind_homolog_sep_stats *b)
{
	int failed = 0;
	failed |= check_i32("det stats n_haploid", a->n_haploid, b->n_haploid);
	failed |= check_i32("det stats n_finite", a->n_finite, b->n_finite);
	failed |= check_i32("det stats n_nonfinite", a->n_nonfinite, b->n_nonfinite);
	failed |= check_i32("det stats n_collapsed", a->n_collapsed, b->n_collapsed);
	failed |= check_close("det stats min", a->min_sep, b->min_sep);
	failed |= check_close("det stats max", a->max_sep, b->max_sep);
	failed |= check_close("det stats mean", a->mean_sep, b->mean_sep);
	failed |= check_close("det stats threshold", a->collapse_threshold, b->collapse_threshold);
	return failed;
}

static int check_init_diagnostics_and_forces(void)
{
	const int32_t n_haploid = 5, n_diploid = 10;
	struct hk_bmap bmap;
	struct hk_bead beads[5];
	struct hk_blind_homolog_sep_stats stats, stats_repeat;
	fvec3_t haploid[5], diploid[10], diploid_repeat[10];
	fvec3_t force_inactive[10], force_inactive_initial[10];
	fvec3_t force_active[10], force_sentinel[10], force_final[10];
	float eps = 0.35f;
	float unit = 1.0f;
	float lambda_sep = 0.05f;
	float inactive_min_sep_unit = 0.6f;
	float active_min_sep_unit = 0.8f;
	float collapse_threshold = 0.6f;
	float expected_active_e = (float)n_haploid * lambda_sep * (active_min_sep_unit - 2.0f * eps / unit) * (active_min_sep_unit - 2.0f * eps / unit);
	float inactive_e, active_e, active_e_sentinel;
	int failed = 0;

	set_test_bmap(&bmap, beads);
	set_haploid_coords(haploid);
	failed |= check_i32("init ret", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, n_haploid, diploid, eps, 0.0f, 12345), 0);
	failed |= check_i32("repeat init ret", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, n_haploid, diploid_repeat, eps, 0.0f, 12345), 0);
	failed |= check_true("same seed identical coords", coords_equal(diploid, diploid_repeat, n_diploid));
	failed |= check_all_coords_finite(diploid, n_diploid);
	failed |= check_initialized_geometry(diploid, haploid, n_haploid, eps);

	failed |= check_i32("stats ret", hk_blind_homolog_sep_compute_stats(n_haploid, diploid, collapse_threshold, &stats), 0);
	failed |= check_i32("repeat stats ret", hk_blind_homolog_sep_compute_stats(n_haploid, diploid_repeat, collapse_threshold, &stats_repeat), 0);
	failed |= check_stats_equal(&stats, &stats_repeat);
	failed |= check_i32("stats n_haploid", stats.n_haploid, n_haploid);
	failed |= check_i32("stats n_finite", stats.n_finite, n_haploid);
	failed |= check_i32("stats n_nonfinite", stats.n_nonfinite, 0);
	failed |= check_i32("stats n_collapsed", stats.n_collapsed, 0);
	failed |= check_close("stats min sep", stats.min_sep, 2.0f * eps);
	failed |= check_close("stats max sep", stats.max_sep, 2.0f * eps);
	failed |= check_close("stats mean sep", stats.mean_sep, 2.0f * eps);

	set_sentinel_force(force_inactive, n_diploid);
	copy_force(force_inactive_initial, force_inactive, n_diploid);
	inactive_e = hk_blind_homolog_sep_accumulate_force(n_haploid, diploid, force_inactive, unit, inactive_min_sep_unit, lambda_sep);
	failed |= check_close("inactive energy", inactive_e, 0.0f);
	failed |= check_true("inactive min threshold below sep", inactive_min_sep_unit * unit < 2.0f * eps);
	failed |= check_true("inactive force unchanged", coords_equal(force_inactive, force_inactive_initial, n_diploid));

	zero_force(force_active, n_diploid);
	active_e = hk_blind_homolog_sep_accumulate_force(n_haploid, diploid, force_active, unit, active_min_sep_unit, lambda_sep);
	failed |= check_true("active min threshold above sep", active_min_sep_unit * unit > 2.0f * eps);
	failed |= check_true("active energy positive", active_e > 0.0f);
	failed |= check_close("active expected energy", active_e, expected_active_e);
	failed |= check_force_finite_nonzero_and_net_zero(force_active, n_diploid);

	set_sentinel_force(force_sentinel, n_diploid);
	copy_force(force_final, force_sentinel, n_diploid);
	active_e_sentinel = hk_blind_homolog_sep_accumulate_force(n_haploid, diploid, force_final, unit, active_min_sep_unit, lambda_sep);
	failed |= check_close("active sentinel energy", active_e_sentinel, active_e);
	failed |= check_force_accumulation(force_sentinel, force_active, force_final, n_diploid);
	return failed;
}

int main(void)
{
	return check_init_diagnostics_and_forces() != 0;
}
