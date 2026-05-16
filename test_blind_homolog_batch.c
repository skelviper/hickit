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

static float make_qnan(void)
{
	union { uint32_t u; float f; } v;
	v.u = 0x7fc00000u;
	return v.f;
}

static void set_coord(fvec3_t x, float x0, float x1, float x2)
{
	x[0] = x0;
	x[1] = x1;
	x[2] = x2;
}

static void set_pair(fvec3_t *coords, int32_t i, float x00, float x01, float x02, float x10, float x11, float x12)
{
	set_coord(coords[hk_diploid_bid(i, HK_DIPLOID_COPY0)], x00, x01, x02);
	set_coord(coords[hk_diploid_bid(i, HK_DIPLOID_COPY1)], x10, x11, x12);
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

static int check_force_eq(const char *label, const fvec3_t *got, const fvec3_t *expected, int32_t n_diploid)
{
	int failed = 0;
	int32_t i;
	int a;
	char buf[96];
	for (i = 0; i < n_diploid; ++i) {
		for (a = 0; a < 3; ++a) {
			snprintf(buf, sizeof(buf), "%s bead%d axis%d", label, (int)i, a);
			failed |= check_close(buf, got[i][a], expected[i][a]);
		}
	}
	return failed;
}

static void add_vec(fvec3_t dst, const fvec3_t src)
{
	int a;
	for (a = 0; a < 3; ++a)
		dst[a] += src[a];
}

static int check_force_finite(const fvec3_t *force, int32_t n_diploid)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n_diploid; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("finite force", isfinite(force[i][a]));
	return failed;
}

static int check_no_beads_accumulator(void)
{
	float e = hk_blind_homolog_sep_accumulate_force(0, 0, 0, 1.0f, 0.5f, 0.05f);
	return check_close("no-beads accumulator energy", e, 0.0f);
}

static int check_beyond_min_sep_unchanged(void)
{
	fvec3_t coords[4], force[4], initial[4];
	float e;
	int failed = 0;

	set_pair(coords, 0, 0.0f, 0.0f, 0.0f, 0.6f, 0.0f, 0.0f);
	set_pair(coords, 1, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f);
	set_sentinel_force(force, 4);
	copy_force(initial, force, 4);
	e = hk_blind_homolog_sep_accumulate_force(2, coords, force, 1.0f, 0.5f, 0.05f);
	failed |= check_close("beyond batch energy", e, 0.0f);
	failed |= check_force_eq("beyond force unchanged", force, initial, 4);
	return failed;
}

static int check_one_collapsed_pair(void)
{
	fvec3_t coords[4], force[4], expected[4], f0, f1;
	float e, single_e;
	int failed = 0;

	set_pair(coords, 0, 0.0f, 0.0f, 0.0f, 0.25f, 0.0f, 0.0f);
	set_pair(coords, 1, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f);
	set_sentinel_force(force, 4);
	copy_force(expected, force, 4);
	single_e = hk_blind_homolog_sep_energy_force(coords[0], coords[1], 1.0f, 0.5f, 0.05f, f0, f1);
	add_vec(expected[0], f0);
	add_vec(expected[1], f1);
	e = hk_blind_homolog_sep_accumulate_force(2, coords, force, 1.0f, 0.5f, 0.05f);
	failed |= check_close("one collapsed energy", e, single_e);
	failed |= check_force_eq("one collapsed accumulated force", force, expected, 4);
	return failed;
}

static int check_multiple_collapsed_pairs_and_net_zero(void)
{
	fvec3_t coords[6], force[6], expected[6], f0, f1, net;
	float e, single_e, expected_e = 0.0f;
	int failed = 0, i, a;

	set_pair(coords, 0, 0.0f, 0.0f, 0.0f, 0.25f, 0.0f, 0.0f);
	set_pair(coords, 1, 2.0f, 0.0f, 0.0f, 2.0f, 0.25f, 0.0f);
	set_pair(coords, 2, 0.0f, 3.0f, 0.0f, 1.0f, 3.0f, 0.0f);
	zero_force(force, 6);
	zero_force(expected, 6);

	single_e = hk_blind_homolog_sep_energy_force(coords[0], coords[1], 1.0f, 0.5f, 0.05f, f0, f1);
	expected_e += single_e;
	add_vec(expected[0], f0);
	add_vec(expected[1], f1);
	single_e = hk_blind_homolog_sep_energy_force(coords[2], coords[3], 1.0f, 0.5f, 0.05f, f0, f1);
	expected_e += single_e;
	add_vec(expected[2], f0);
	add_vec(expected[3], f1);

	e = hk_blind_homolog_sep_accumulate_force(3, coords, force, 1.0f, 0.5f, 0.05f);
	failed |= check_close("multiple collapsed energy", e, expected_e);
	failed |= check_force_eq("multiple collapsed force", force, expected, 6);
	set_coord(net, 0.0f, 0.0f, 0.0f);
	for (i = 0; i < 6; ++i)
		for (a = 0; a < 3; ++a)
			net[a] += force[i][a];
	failed |= check_close("net force x", net[0], 0.0f);
	failed |= check_close("net force y", net[1], 0.0f);
	failed |= check_close("net force z", net[2], 0.0f);
	return failed;
}

static int check_force_accumulation(void)
{
	fvec3_t coords[2], force[2], expected[2], f0, f1;
	float e, single_e;
	int failed = 0;

	set_pair(coords, 0, 0.0f, 0.0f, 0.0f, 0.25f, 0.0f, 0.0f);
	set_sentinel_force(force, 2);
	copy_force(expected, force, 2);
	single_e = hk_blind_homolog_sep_energy_force(coords[0], coords[1], 1.0f, 0.5f, 0.05f, f0, f1);
	add_vec(expected[0], f0);
	add_vec(expected[1], f1);
	e = hk_blind_homolog_sep_accumulate_force(1, coords, force, 1.0f, 0.5f, 0.05f);
	failed |= check_close("accumulation energy", e, single_e);
	failed |= check_force_eq("accumulation force", force, expected, 2);
	return failed;
}

static int check_zero_distance_fallback_batch(void)
{
	fvec3_t coords[2], force[2];
	float e, expected_e = 0.05f * 0.5f * 0.5f;
	float expected_f = 2.0f * 0.05f * 0.5f;
	int failed = 0;

	set_pair(coords, 0, 1.0f, 2.0f, 3.0f, 1.0f, 2.0f, 3.0f);
	zero_force(force, 2);
	e = hk_blind_homolog_sep_accumulate_force(1, coords, force, 1.0f, 0.5f, 0.05f);
	failed |= check_close("zero batch energy", e, expected_e);
	failed |= check_force_finite(force, 2);
	failed |= check_close("zero batch f0 x", force[0][0], expected_f);
	failed |= check_close("zero batch f1 x", force[1][0], -expected_f);
	failed |= check_close("zero batch f0 y", force[0][1], 0.0f);
	failed |= check_close("zero batch f0 z", force[0][2], 0.0f);
	failed |= check_close("zero batch f1 y", force[1][1], 0.0f);
	failed |= check_close("zero batch f1 z", force[1][2], 0.0f);
	return failed;
}

static int check_lambda_zero_batch(void)
{
	fvec3_t coords[2], force[2], initial[2];
	float e;
	int failed = 0;

	set_pair(coords, 0, 0.0f, 0.0f, 0.0f, 0.25f, 0.0f, 0.0f);
	set_sentinel_force(force, 2);
	copy_force(initial, force, 2);
	e = hk_blind_homolog_sep_accumulate_force(1, coords, force, 1.0f, 0.5f, 0.0f);
	failed |= check_close("lambda-zero batch energy", e, 0.0f);
	failed |= check_force_eq("lambda-zero force unchanged", force, initial, 2);
	return failed;
}

static int check_basic_stats(void)
{
	struct hk_blind_homolog_sep_stats stats;
	fvec3_t coords[6];
	int failed = 0;

	set_pair(coords, 0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
	set_pair(coords, 1, 0.0f, 0.0f, 0.0f, 0.5f, 0.0f, 0.0f);
	set_pair(coords, 2, 0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f);
	failed |= check_i32("basic stats ret", hk_blind_homolog_sep_compute_stats(3, coords, 0.6f, &stats), 0);
	failed |= check_i32("basic n_haploid", stats.n_haploid, 3);
	failed |= check_i32("basic n_finite", stats.n_finite, 3);
	failed |= check_i32("basic n_nonfinite", stats.n_nonfinite, 0);
	failed |= check_i32("basic n_collapsed", stats.n_collapsed, 2);
	failed |= check_close("basic min", stats.min_sep, 0.0f);
	failed |= check_close("basic max", stats.max_sep, 2.0f);
	failed |= check_close("basic mean", stats.mean_sep, (0.0f + 0.5f + 2.0f) / 3.0f);
	failed |= check_close("basic threshold", stats.collapse_threshold, 0.6f);
	return failed;
}

static int check_no_beads_stats(void)
{
	struct hk_blind_homolog_sep_stats stats;
	int failed = 0;

	failed |= check_i32("empty stats ret", hk_blind_homolog_sep_compute_stats(0, 0, 0.5f, &stats), 0);
	failed |= check_i32("empty n_haploid", stats.n_haploid, 0);
	failed |= check_i32("empty n_finite", stats.n_finite, 0);
	failed |= check_i32("empty n_nonfinite", stats.n_nonfinite, 0);
	failed |= check_i32("empty n_collapsed", stats.n_collapsed, 0);
	failed |= check_close("empty min", stats.min_sep, 0.0f);
	failed |= check_close("empty max", stats.max_sep, 0.0f);
	failed |= check_close("empty mean", stats.mean_sep, 0.0f);
	failed |= check_close("empty threshold", stats.collapse_threshold, 0.5f);
	return failed;
}

static int check_nonfinite_stats(void)
{
	struct hk_blind_homolog_sep_stats stats;
	fvec3_t coords[6];
	int failed = 0;

	set_pair(coords, 0, 0.0f, 0.0f, 0.0f, 0.5f, 0.0f, 0.0f);
	set_pair(coords, 1, 1.0f, 0.0f, 0.0f, make_qnan(), 0.0f, 0.0f);
	set_pair(coords, 2, 0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f);
	failed |= check_i32("nonfinite stats ret", hk_blind_homolog_sep_compute_stats(3, coords, 1.0f, &stats), 0);
	failed |= check_i32("nonfinite n_haploid", stats.n_haploid, 3);
	failed |= check_i32("nonfinite n_finite", stats.n_finite, 2);
	failed |= check_i32("nonfinite n_nonfinite", stats.n_nonfinite, 1);
	failed |= check_i32("nonfinite n_collapsed", stats.n_collapsed, 1);
	failed |= check_close("nonfinite min", stats.min_sep, 0.5f);
	failed |= check_close("nonfinite max", stats.max_sep, 2.0f);
	failed |= check_close("nonfinite mean", stats.mean_sep, (0.5f + 2.0f) / 2.0f);
	return failed;
}

static int check_threshold_boundary_stats(void)
{
	struct hk_blind_homolog_sep_stats stats;
	fvec3_t coords[2];
	int failed = 0;

	set_pair(coords, 0, 0.0f, 0.0f, 0.0f, 0.6f, 0.0f, 0.0f);
	failed |= check_i32("boundary stats ret", hk_blind_homolog_sep_compute_stats(1, coords, 0.6f, &stats), 0);
	failed |= check_i32("boundary collapsed", stats.n_collapsed, 0);
	failed |= check_close("boundary min", stats.min_sep, 0.6f);
	failed |= check_close("boundary max", stats.max_sep, 0.6f);
	failed |= check_close("boundary mean", stats.mean_sep, 0.6f);
	return failed;
}

static void set_test_bmap(struct hk_bmap *bmap, struct hk_bead beads[4])
{
	int i;
	for (i = 0; i < 4; ++i) {
		beads[i].chr = i < 2? 0 : 1;
		beads[i].st = 1000000 * i;
		beads[i].en = 1000000 * (i + 1);
	}
	bmap->n_beads = 4;
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

static void set_haploid_coords(fvec3_t x[4])
{
	set_coord(x[0], 1.0f, 2.0f, 3.0f);
	set_coord(x[1], 4.0f, 5.0f, 6.0f);
	set_coord(x[2], -1.0f, 7.0f, 0.5f);
	set_coord(x[3], 2.5f, -3.0f, 9.0f);
}

static int check_initialized_coords_stats(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	struct hk_blind_homolog_sep_stats stats;
	fvec3_t haploid[4], diploid[8];
	float eps = 0.35f;
	int failed = 0;

	set_test_bmap(&bmap, beads);
	set_haploid_coords(haploid);
	failed |= check_i32("init coords ret", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, 4, diploid, eps, 0.0f, 101), 0);
	failed |= check_i32("init stats ret", hk_blind_homolog_sep_compute_stats(4, diploid, 0.6f, &stats), 0);
	failed |= check_i32("init n_finite", stats.n_finite, 4);
	failed |= check_i32("init n_nonfinite", stats.n_nonfinite, 0);
	failed |= check_i32("init n_collapsed", stats.n_collapsed, 0);
	failed |= check_close("init min sep", stats.min_sep, 2.0f * eps);
	failed |= check_close("init max sep", stats.max_sep, 2.0f * eps);
	failed |= check_close("init mean sep", stats.mean_sep, 2.0f * eps);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_no_beads_accumulator();
	failed |= check_beyond_min_sep_unchanged();
	failed |= check_one_collapsed_pair();
	failed |= check_multiple_collapsed_pairs_and_net_zero();
	failed |= check_force_accumulation();
	failed |= check_zero_distance_fallback_batch();
	failed |= check_lambda_zero_batch();
	failed |= check_basic_stats();
	failed |= check_no_beads_stats();
	failed |= check_nonfinite_stats();
	failed |= check_threshold_boundary_stats();
	failed |= check_initialized_coords_stats();
	return failed != 0;
}
