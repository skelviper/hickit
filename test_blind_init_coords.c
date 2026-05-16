#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
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
	x[0][0] = 1.0f;  x[0][1] = 2.0f;  x[0][2] = 3.0f;
	x[1][0] = 4.0f;  x[1][1] = 5.0f;  x[1][2] = 6.0f;
	x[2][0] = -1.0f; x[2][1] = 7.0f;  x[2][2] = 0.5f;
	x[3][0] = 2.5f;  x[3][1] = -3.0f; x[3][2] = 9.0f;
}

static float dist3(const fvec3_t a, const fvec3_t b)
{
	float dx = a[0] - b[0];
	float dy = a[1] - b[1];
	float dz = a[2] - b[2];
	return sqrtf(dx * dx + dy * dy + dz * dz);
}

static int check_coord(const char *label, const fvec3_t got, const fvec3_t expected)
{
	int failed = 0;
	char buf[64];
	snprintf(buf, sizeof(buf), "%s x", label);
	failed |= check_close(buf, got[0], expected[0]);
	snprintf(buf, sizeof(buf), "%s y", label);
	failed |= check_close(buf, got[1], expected[1]);
	snprintf(buf, sizeof(buf), "%s z", label);
	failed |= check_close(buf, got[2], expected[2]);
	return failed;
}

static int check_zero_separation_zero_noise(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t haploid[4], diploid[8];
	int failed = 0, i;

	set_test_bmap(&bmap, beads);
	set_haploid_coords(haploid);
	failed |= check_i32("zero init ret", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, 4, diploid, 0.0f, 0.0f, 17), 0);
	for (i = 0; i < 4; ++i) {
		failed |= check_coord("zero copy0", diploid[hk_diploid_bid(i, HK_DIPLOID_COPY0)], haploid[i]);
		failed |= check_coord("zero copy1", diploid[hk_diploid_bid(i, HK_DIPLOID_COPY1)], haploid[i]);
	}
	return failed;
}

static int check_separation_no_noise(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t haploid[4], diploid[8];
	float eps = 1.25f;
	int failed = 0, i, a;

	set_test_bmap(&bmap, beads);
	set_haploid_coords(haploid);
	failed |= check_i32("separation init ret", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, 4, diploid, eps, 0.0f, 91), 0);
	for (i = 0; i < 4; ++i) {
		int32_t d0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t d1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		for (a = 0; a < 3; ++a)
			failed |= check_close("separation midpoint", 0.5f * (diploid[d0][a] + diploid[d1][a]), haploid[i][a]);
		failed |= check_close("separation magnitude", dist3(diploid[d0], diploid[d1]), 2.0f * eps);
	}
	return failed;
}

static int check_same_chromosome_direction(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t haploid[4], diploid[8];
	float off0[3], off1[3], off2[3];
	int failed = 0, a;

	set_test_bmap(&bmap, beads);
	set_haploid_coords(haploid);
	failed |= check_i32("same-chr init ret", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, 4, diploid, 0.75f, 0.0f, 123), 0);
	for (a = 0; a < 3; ++a) {
		off0[a] = diploid[hk_diploid_bid(0, HK_DIPLOID_COPY0)][a] - haploid[0][a];
		off1[a] = diploid[hk_diploid_bid(1, HK_DIPLOID_COPY0)][a] - haploid[1][a];
		off2[a] = diploid[hk_diploid_bid(2, HK_DIPLOID_COPY0)][a] - haploid[2][a];
		failed |= check_close("same-chr offset", off0[a], off1[a]);
	}
	failed |= check_true("different-chr offset differs", fabsf(off0[0] - off2[0]) > 1e-6f ||
													 fabsf(off0[1] - off2[1]) > 1e-6f ||
													 fabsf(off0[2] - off2[2]) > 1e-6f);
	return failed;
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

static int coords_any_differ(const fvec3_t *a, const fvec3_t *b, int32_t n)
{
	return !coords_equal(a, b, n);
}

static int coords_all_finite(const fvec3_t *x, int32_t n)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			if (!isfinite(x[i][a]))
				return 0;
	return 1;
}

static int check_seed_determinism(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t haploid[4], diploid_a[8], diploid_b[8], diploid_c[8];
	int failed = 0;

	set_test_bmap(&bmap, beads);
	set_haploid_coords(haploid);
	failed |= check_i32("determinism init a", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, 4, diploid_a, 0.5f, 0.1f, 777), 0);
	failed |= check_i32("determinism init b", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, 4, diploid_b, 0.5f, 0.1f, 777), 0);
	failed |= check_i32("determinism init c", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, 4, diploid_c, 0.5f, 0.1f, 778), 0);
	failed |= check_true("same seed identical", coords_equal(diploid_a, diploid_b, 8));
	failed |= check_true("different seed differs", coords_any_differ(diploid_a, diploid_c, 8));
	return failed;
}

static int check_noise_bound(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t haploid[4], diploid[8];
	float noise_scale = 0.125f;
	int failed = 0, i, a;

	set_test_bmap(&bmap, beads);
	set_haploid_coords(haploid);
	failed |= check_i32("noise-bound init ret", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, 4, diploid, 0.0f, noise_scale, 31337), 0);
	for (i = 0; i < 4; ++i) {
		int32_t d0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t d1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		for (a = 0; a < 3; ++a) {
			failed |= check_true("noise copy0 bounded", fabsf(diploid[d0][a] - haploid[i][a]) <= noise_scale + 1e-6f);
			failed |= check_true("noise copy1 bounded", fabsf(diploid[d1][a] - haploid[i][a]) <= noise_scale + 1e-6f);
		}
	}
	return failed;
}

static int check_toy_haploid_scaffold_available(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t haploid[4];
	int failed = 0;

	set_test_bmap(&bmap, beads);
	failed |= check_i32("toy scaffold ret", hk_blind_init_toy_haploid_scaffold(&bmap, haploid), 0);
	failed |= check_close("toy scaffold bead0 x", haploid[0][0], 0.0f);
	failed |= check_close("toy scaffold bead1 x", haploid[1][0], 0.10f);
	failed |= check_close("toy scaffold bead2 z", haploid[2][2], 0.06f);
	failed |= check_true("toy scaffold finite", coords_all_finite(haploid, 4));
	return failed;
}

static int check_random_haploid_scaffold_seed(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t a[4], b[4], c[4];
	int failed = 0;

	set_test_bmap(&bmap, beads);
	failed |= check_i32("random haploid a", hk_blind_init_random_haploid_scaffold(&bmap, a, 3.0f, 1001), 0);
	failed |= check_i32("random haploid b", hk_blind_init_random_haploid_scaffold(&bmap, b, 3.0f, 1001), 0);
	failed |= check_i32("random haploid c", hk_blind_init_random_haploid_scaffold(&bmap, c, 3.0f, 1002), 0);
	failed |= check_true("random haploid same seed identical", coords_equal(a, b, 4));
	failed |= check_true("random haploid different seed differs", coords_any_differ(a, c, 4));
	failed |= check_true("random haploid finite", coords_all_finite(a, 4));
	return failed;
}

static int check_random_diploid_seed(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t a[8], b[8], c[8];
	int failed = 0;

	set_test_bmap(&bmap, beads);
	failed |= check_i32("random diploid a", hk_blind_init_random_diploid_coords(&bmap, a, 3.0f, 2001), 0);
	failed |= check_i32("random diploid b", hk_blind_init_random_diploid_coords(&bmap, b, 3.0f, 2001), 0);
	failed |= check_i32("random diploid c", hk_blind_init_random_diploid_coords(&bmap, c, 3.0f, 2002), 0);
	failed |= check_true("random diploid same seed identical", coords_equal(a, b, 8));
	failed |= check_true("random diploid different seed differs", coords_any_differ(a, c, 8));
	failed |= check_true("random diploid finite", coords_all_finite(a, 8));
	return failed;
}

static int check_unphased_fdg_scaffold_seed(void)
{
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_fdg_conf conf;
	fvec3_t *a = 0, *b = 0, *c = 0;
	int failed = 0;

	m = hk_map_read("testdata/p9016_blind_fixture.pairs");
	if (m == 0) {
		fprintf(stderr, "unphased fdg fixture: failed to read map\n");
		return 1;
	}
	bmap = hk_bmap_gen(m->d, m->n_pairs, m->pairs, 1000000, 1);
	if (bmap == 0 || bmap->n_beads <= 0) {
		fprintf(stderr, "unphased fdg fixture: bad bmap\n");
		if (bmap) hk_bmap_destroy(bmap);
		hk_map_destroy(m);
		return 1;
	}
	a = (fvec3_t*)malloc((size_t)bmap->n_beads * sizeof(fvec3_t));
	b = (fvec3_t*)malloc((size_t)bmap->n_beads * sizeof(fvec3_t));
	c = (fvec3_t*)malloc((size_t)bmap->n_beads * sizeof(fvec3_t));
	if (a == 0 || b == 0 || c == 0) {
		fprintf(stderr, "unphased fdg fixture: failed to allocate coords\n");
		free(a); free(b); free(c);
		hk_bmap_destroy(bmap);
		hk_map_destroy(m);
		return 1;
	}
	hk_fdg_conf_init(&conf);
	conf.backend = HK_FDG_BACKEND_CPU;
	conf.n_iter = 2;
	failed |= check_i32("unphased fdg scaffold a", hk_blind_init_haploid_scaffold_from_bmap_fdg(bmap, &conf, a, 3001), 0);
	failed |= check_i32("unphased fdg scaffold b", hk_blind_init_haploid_scaffold_from_bmap_fdg(bmap, &conf, b, 3001), 0);
	failed |= check_i32("unphased fdg scaffold c", hk_blind_init_haploid_scaffold_from_bmap_fdg(bmap, &conf, c, 3002), 0);
	failed |= check_true("unphased fdg same seed identical", coords_equal(a, b, bmap->n_beads));
	failed |= check_true("unphased fdg different seed differs", coords_any_differ(a, c, bmap->n_beads));
	failed |= check_true("unphased fdg finite", coords_all_finite(a, bmap->n_beads));
	free(a); free(b); free(c);
	hk_bmap_destroy(bmap);
	hk_map_destroy(m);
	return failed;
}

static int check_diploid_indexing_and_finite(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t haploid[4], diploid[8];
	int failed = 0, i, a;

	set_test_bmap(&bmap, beads);
	set_haploid_coords(haploid);
	failed |= check_i32("indexing init ret", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, 4, diploid, 0.25f, 0.05f, 99), 0);
	failed |= check_i32("index copy0", hk_diploid_bid(2, HK_DIPLOID_COPY0), 4);
	failed |= check_i32("index copy1", hk_diploid_bid(2, HK_DIPLOID_COPY1), 5);
	for (i = 0; i < 8; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("finite diploid coordinate", isfinite(diploid[i][a]));
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_zero_separation_zero_noise();
	failed |= check_separation_no_noise();
	failed |= check_same_chromosome_direction();
	failed |= check_seed_determinism();
	failed |= check_noise_bound();
	failed |= check_toy_haploid_scaffold_available();
	failed |= check_random_haploid_scaffold_seed();
	failed |= check_random_diploid_seed();
	failed |= check_unphased_fdg_scaffold_seed();
	failed |= check_diploid_indexing_and_finite();
	return failed != 0;
}
