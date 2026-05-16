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

static int check_u8(const char *label, uint8_t got, uint8_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %u, expected %u\n", label, got, expected);
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

static int check_range(const char *label, float x, float lo, float hi)
{
	float tol = 1e-5f;
	if (!(x >= lo - tol && x <= hi + tol)) {
		fprintf(stderr, "%s: expected %.8g in [%.8g, %.8g]\n", label, x, lo, hi);
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

static void init_synthetic_bmap(struct hk_bmap *bmap, struct hk_sdict *dict, int32_t len[1],
								struct hk_bead beads[6], uint64_t offcnt[2])
{
	int i;
	len[0] = 6000000;
	dict->n = 1;
	dict->m = 1;
	dict->name = 0;
	dict->len = len;
	dict->h = 0;
	for (i = 0; i < 6; ++i) {
		beads[i].chr = 0;
		beads[i].st = 1000000 * i;
		beads[i].en = 1000000 * (i + 1);
	}
	offcnt[0] = (uint64_t)6;
	offcnt[1] = (uint64_t)6 << 32;
	bmap->n_beads = 6;
	bmap->n_pairs = 0;
	bmap->unit = 1.0f;
	bmap->d = dict;
	bmap->beads = beads;
	bmap->offcnt = offcnt;
	bmap->pairs = 0;
	bmap->x = 0;
	bmap->feat = 0;
	bmap->cpg = 0;
	bmap->gc_bias = 0;
	bmap->gc_corrected = 0;
}

static void init_haploid_scaffold(fvec3_t haploid[6])
{
	set_coord(haploid[0], 0.0f, 0.0f, 0.0f);
	set_coord(haploid[1], 1.0f, 0.0f, 0.0f);
	set_coord(haploid[2], 2.0f, 0.5f, 0.0f);
	set_coord(haploid[3], 4.0f, 1.0f, 0.0f);
	set_coord(haploid[4], 7.0f, 1.5f, 0.0f);
	set_coord(haploid[5], 8.0f, 2.0f, 0.0f);
}

static struct hk_blind_pair make_blind_pair(int32_t bid0, int32_t bid1)
{
	struct hk_blind_pair p;
	p.chr[0] = 0;
	p.chr[1] = 0;
	p.pos[0] = bid0 * 1000000 + 1234;
	p.pos[1] = bid1 * 1000000 + 5678;
	p.strand[0] = 1;
	p.strand[1] = -1;
	return p;
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

static int check_diploid_separation(const fvec3_t *diploid, int32_t n_haploid, float eps)
{
	int failed = 0;
	int32_t i;
	for (i = 0; i < n_haploid; ++i) {
		int32_t d0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t d1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		failed |= check_close("init homolog separation", dist3(diploid[d0], diploid[d1]), 2.0f * eps);
	}
	return failed;
}

static int check_midpoints(const fvec3_t *haploid, const fvec3_t *diploid, int32_t n_haploid)
{
	int failed = 0;
	int32_t i;
	int k;
	for (i = 0; i < n_haploid; ++i) {
		int32_t d0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t d1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		for (k = 0; k < 3; ++k)
			failed |= check_close("init midpoint", 0.5f * (diploid[d0][k] + diploid[d1][k]), haploid[i][k]);
	}
	return failed;
}

static int check_bpair_diagnostics(const struct hk_blind_bpair *p)
{
	float p_sum = 0.0f, out_sum;
	int failed = 0, i;
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		failed |= check_true("posterior p4 finite", isfinite(p->p4[i]));
		failed |= check_range("posterior p4 range", p->p4[i], 0.0f, 1.0f);
		p_sum += p->p4[i];
	}
	out_sum = p->rho_output * p_sum + p->pU;
	failed |= check_true("posterior entropy finite", isfinite(p->entropy));
	failed |= check_true("posterior pmax finite", isfinite(p->pmax));
	failed |= check_true("posterior margin finite", isfinite(p->margin));
	failed |= check_true("posterior rho finite", isfinite(p->rho_output));
	failed |= check_true("posterior pU finite", isfinite(p->pU));
	failed |= check_close("posterior p4 sum", p_sum, 1.0f);
	failed |= check_close("posterior output-with-U sum", out_sum, 1.0f);
	failed |= check_range("posterior entropy range", p->entropy, 0.0f, logf(4.0f));
	failed |= check_range("posterior pmax range", p->pmax, 0.25f, 1.0f);
	failed |= check_range("posterior margin range", p->margin, 0.0f, 1.0f);
	failed |= check_range("posterior rho range", p->rho_output, 0.0f, 1.0f);
	failed |= check_range("posterior pU range", p->pU, 0.0f, 1.0f);
	return failed;
}

static int check_init_to_posterior_smoke(void)
{
	struct hk_bmap bmap;
	struct hk_sdict dict;
	struct hk_bead beads[6];
	struct hk_blind_pair raw[4];
	struct hk_blind_bpair_set *set;
	struct hk_fdg_conf conf;
	uint64_t offcnt[2];
	int32_t len[1];
	fvec3_t haploid[6], diploid_a[12], diploid_b[12];
	float log_prior[HK_BLIND_N_STATE];
	float eps = 0.75f;
	int failed = 0, i, n_nonuniform = 0;

	init_synthetic_bmap(&bmap, &dict, len, beads, offcnt);
	init_haploid_scaffold(haploid);
	failed |= check_i32("init smoke ret a", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, 6, diploid_a, eps, 0.0f, 20240517), 0);
	failed |= check_i32("init smoke ret b", hk_blind_init_diploid_coords_from_haploid(&bmap, haploid, 6, diploid_b, eps, 0.0f, 20240517), 0);
	failed |= check_true("init smoke deterministic", coords_equal(diploid_a, diploid_b, 12));
	failed |= check_diploid_separation(diploid_a, 6, eps);
	failed |= check_midpoints(haploid, diploid_a, 6);

	raw[0] = make_blind_pair(2, 5);
	raw[1] = make_blind_pair(5, 2);
	raw[2] = make_blind_pair(1, 4);
	raw[3] = make_blind_pair(4, 1);
	set = hk_blind_bpair_set_build(&bmap, 4, raw);
	failed |= check_i32("init smoke n_raw", set->n_raw, 4);
	failed |= check_i32("init smoke n_bpairs", set->n_bpairs, 2);
	failed |= check_i32("init smoke key0 bid0", set->bpairs[0].key.bid[0], 1);
	failed |= check_i32("init smoke key0 bid1", set->bpairs[0].key.bid[1], 4);
	failed |= check_i32("init smoke key1 bid0", set->bpairs[1].key.bid[0], 2);
	failed |= check_i32("init smoke key1 bid1", set->bpairs[1].key.bid[1], 5);
	failed |= check_i32("init smoke raw0 bpair", set->raw2binned[0].bpair_id, 1);
	failed |= check_u8("init smoke raw0 swapped", set->raw2binned[0].swapped, 0);
	failed |= check_i32("init smoke raw1 bpair", set->raw2binned[1].bpair_id, 1);
	failed |= check_u8("init smoke raw1 swapped", set->raw2binned[1].swapped, 1);
	failed |= check_i32("init smoke raw2 bpair", set->raw2binned[2].bpair_id, 0);
	failed |= check_u8("init smoke raw2 swapped", set->raw2binned[2].swapped, 0);
	failed |= check_i32("init smoke raw3 bpair", set->raw2binned[3].bpair_id, 0);
	failed |= check_u8("init smoke raw3 swapped", set->raw2binned[3].swapped, 1);

	hk_fdg_conf_init(&conf);
	hk_blind_init_uniform_log_prior(log_prior);
	hk_blind_bpair_set_update_posterior_from_coords(set, &conf, diploid_a, 1.0f, 1.0f, 2.0f, log_prior, 1.0f);
	for (i = 0; i < set->n_bpairs; ++i) {
		failed |= check_bpair_diagnostics(&set->bpairs[i]);
		if (set->bpairs[i].pmax > 0.25f + 1e-5f && set->bpairs[i].pU < 1.0f - 1e-5f)
			++n_nonuniform;
	}
	failed |= check_true("init smoke nonuniform posterior", n_nonuniform > 0);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

int main(void)
{
	return check_init_to_posterior_smoke() != 0;
}
