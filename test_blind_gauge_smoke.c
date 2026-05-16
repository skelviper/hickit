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
		fprintf(stderr, "%s: got %u, expected %u\n", label, (unsigned)got, (unsigned)expected);
		return 1;
	}
	return 0;
}

static int check_close(const char *label, double got, double expected)
{
	double tol = 1e-6;
	double scale = fabs(expected) > 1.0? fabs(expected) : 1.0;
	if (fabs(got - expected) > tol * scale) {
		fprintf(stderr, "%s: got %.12g, expected %.12g\n", label, got, expected);
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

static void set_pair(fvec3_t *coords, int32_t i, float x00, float x01, float x02, float x10, float x11, float x12)
{
	set_coord(coords[hk_diploid_bid(i, HK_DIPLOID_COPY0)], x00, x01, x02);
	set_coord(coords[hk_diploid_bid(i, HK_DIPLOID_COPY1)], x10, x11, x12);
}

static void copy_coords(fvec3_t *dst, const fvec3_t *src, int32_t n_diploid)
{
	int32_t i;
	int a;
	for (i = 0; i < n_diploid; ++i)
		for (a = 0; a < 3; ++a)
			dst[i][a] = src[i][a];
}

static int check_coords_close(const char *label, const fvec3_t *got, const fvec3_t *expected, int32_t n_diploid)
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

static int check_chr_coords_close(const char *label, const struct hk_bmap *bmap, int32_t chr, const fvec3_t *got, const fvec3_t *expected)
{
	int failed = 0;
	int32_t i;
	for (i = 0; i < bmap->n_beads; ++i) {
		if (bmap->beads[i].chr == chr) {
			int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
			int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
			failed |= check_coords_close(label, &got[b0], &expected[b0], 1);
			failed |= check_coords_close(label, &got[b1], &expected[b1], 1);
		}
	}
	return failed;
}

static int check_p4_close(const char *label, const float got[HK_BLIND_N_STATE], const float expected[HK_BLIND_N_STATE])
{
	int failed = 0;
	int i;
	char buf[64];
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		snprintf(buf, sizeof(buf), "%s state%d", label, i);
		failed |= check_close(buf, got[i], expected[i]);
	}
	return failed;
}

static void set_test_bmap(struct hk_bmap *bmap, struct hk_bead beads[4])
{
	int32_t chr[4] = {0, 1, 1, 0};
	int32_t i;
	for (i = 0; i < 4; ++i) {
		beads[i].chr = chr[i];
		beads[i].st = i * 1000000;
		beads[i].en = beads[i].st + 1000000;
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

static void set_prev_coords(fvec3_t coords[8])
{
	set_pair(coords, 0, 0.0f, 0.0f, 0.0f, 10.0f, 0.0f, 0.0f);
	set_pair(coords, 1, 100.0f, 0.0f, 0.0f, 110.0f, 0.0f, 0.0f);
	set_pair(coords, 2, 100.0f, 1.0f, 0.0f, 110.0f, 1.0f, 0.0f);
	set_pair(coords, 3, 0.0f, 1.0f, 0.0f, 10.0f, 1.0f, 0.0f);
}

static void swap_pair_copies(fvec3_t *coords, int32_t i)
{
	int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
	int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
	float t;
	int a;
	for (a = 0; a < 3; ++a) {
		t = coords[b0][a];
		coords[b0][a] = coords[b1][a];
		coords[b1][a] = t;
	}
}

static void calc_uncertainty(const float p4[HK_BLIND_N_STATE], float *entropy, float *pmax, float *margin,
							 float *rho_output, float *pU)
{
	float top1 = 0.0f, top2 = 0.0f;
	int i;
	*entropy = 0.0f;
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		if (p4[i] > 0.0f)
			*entropy -= p4[i] * logf(p4[i]);
		if (p4[i] > top1) {
			top2 = top1;
			top1 = p4[i];
		} else if (p4[i] > top2) {
			top2 = p4[i];
		}
	}
	*pmax = top1;
	*margin = top1 - top2;
	*rho_output = 1.0f - *entropy / logf((float)HK_BLIND_N_STATE);
	if (*rho_output < 0.0f) *rho_output = 0.0f;
	if (*rho_output > 1.0f) *rho_output = 1.0f;
	*pU = 1.0f - *rho_output;
}

static void set_bpair(struct hk_blind_bpair *bp, int32_t bid0, int32_t bid1, const float p4[HK_BLIND_N_STATE])
{
	int i;
	bp->key.bid[0] = bid0;
	bp->key.bid[1] = bid1;
	bp->n_raw = 1;
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		bp->p4[i] = p4[i];
	calc_uncertainty(bp->p4, &bp->entropy, &bp->pmax, &bp->margin, &bp->rho_output, &bp->pU);
}

static int check_derived_fields_invariant(const char *label, const struct hk_blind_bpair *bp,
										  float entropy, float pmax, float margin, float rho_output, float pU)
{
	int failed = 0;
	char buf[96];
	snprintf(buf, sizeof(buf), "%s entropy", label);
	failed |= check_close(buf, bp->entropy, entropy);
	snprintf(buf, sizeof(buf), "%s pmax", label);
	failed |= check_close(buf, bp->pmax, pmax);
	snprintf(buf, sizeof(buf), "%s margin", label);
	failed |= check_close(buf, bp->margin, margin);
	snprintf(buf, sizeof(buf), "%s rho", label);
	failed |= check_close(buf, bp->rho_output, rho_output);
	snprintf(buf, sizeof(buf), "%s pU", label);
	failed |= check_close(buf, bp->pU, pU);
	return failed;
}

static int check_gauge_bookkeeping_smoke(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t prev[8], cur[8], cur_before[8];
	struct hk_blind_bpair bpairs[4];
	struct hk_blind_raw2binned raw2binned[2];
	struct hk_blind_bpair_set set;
	struct hk_blind_gauge_stats stats;
	uint8_t chr_flipped[2] = {0, 0};
	float p4[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	float entropy0[4], pmax0[4], margin0[4], rho0[4], pU0[4];
	float expected_a[HK_BLIND_N_STATE] = {0.3f, 0.4f, 0.1f, 0.2f};
	float expected_b[HK_BLIND_N_STATE] = {0.2f, 0.1f, 0.4f, 0.3f};
	float expected_c[HK_BLIND_N_STATE] = {0.4f, 0.3f, 0.2f, 0.1f};
	float expected_d[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	float expected_raw_swapped[HK_BLIND_N_STATE] = {0.3f, 0.1f, 0.4f, 0.2f};
	float expected_raw_ordered[HK_BLIND_N_STATE] = {0.3f, 0.4f, 0.1f, 0.2f};
	float raw_p4[HK_BLIND_N_STATE];
	int failed = 0;
	int i;

	set_test_bmap(&bmap, beads);
	set_prev_coords(prev);
	copy_coords(cur, prev, 8);
	swap_pair_copies(cur, 0);
	swap_pair_copies(cur, 3);
	copy_coords(cur_before, cur, 8);

	failed |= check_i32("gauge ret", hk_blind_temporal_gauge_stabilize_bmap(&bmap, prev, cur, chr_flipped, &stats), 0);
	failed |= check_u8("chr0 flipped", chr_flipped[0], 1);
	failed |= check_u8("chr1 not flipped", chr_flipped[1], 0);
	failed |= check_i32("stats n_chr", stats.n_chr, 2);
	failed |= check_i32("stats n_flipped", stats.n_flipped, 1);
	failed |= check_chr_coords_close("chr0 stabilized", &bmap, 0, cur, prev);
	failed |= check_chr_coords_close("chr1 unchanged", &bmap, 1, cur, cur_before);
	failed |= check_coords_close("all stabilized", cur, prev, 8);

	// The bead order is interleaved to cover both endpoint0-only and endpoint1-only canonical flips.
	set_bpair(&bpairs[0], 0, 1, p4); // A: chr0, chr1
	set_bpair(&bpairs[1], 2, 3, p4); // B: chr1, chr0
	set_bpair(&bpairs[2], 0, 3, p4); // C: chr0, chr0
	set_bpair(&bpairs[3], 1, 2, p4); // D: chr1, chr1
	for (i = 0; i < 4; ++i) {
		entropy0[i] = bpairs[i].entropy;
		pmax0[i] = bpairs[i].pmax;
		margin0[i] = bpairs[i].margin;
		rho0[i] = bpairs[i].rho_output;
		pU0[i] = bpairs[i].pU;
	}
	hk_blind_raw2binned_set(&raw2binned[0], 0, 1);
	hk_blind_raw2binned_set(&raw2binned[1], 0, 0);
	set.bpairs = bpairs;
	set.n_bpairs = 4;
	set.raw2binned = raw2binned;
	set.n_raw = 2;

	failed |= check_i32("apply p4 flips ret", hk_blind_bpair_set_apply_chr_flips(&set, &bmap, chr_flipped, 2), 0);
	failed |= check_p4_close("bpair A endpoint0 flipped", bpairs[0].p4, expected_a);
	failed |= check_p4_close("bpair B endpoint1 flipped", bpairs[1].p4, expected_b);
	failed |= check_p4_close("bpair C both flipped", bpairs[2].p4, expected_c);
	failed |= check_p4_close("bpair D no flip", bpairs[3].p4, expected_d);
	for (i = 0; i < 4; ++i)
		failed |= check_derived_fields_invariant("derived invariant", &bpairs[i], entropy0[i], pmax0[i], margin0[i], rho0[i], pU0[i]);

	hk_blind_p4_to_raw_order(bpairs[raw2binned[0].bpair_id].p4, raw2binned[0].swapped, raw_p4);
	failed |= check_p4_close("raw swapped after gauge", raw_p4, expected_raw_swapped);
	hk_blind_p4_to_raw_order(bpairs[raw2binned[1].bpair_id].p4, raw2binned[1].swapped, raw_p4);
	failed |= check_p4_close("raw ordered after gauge", raw_p4, expected_raw_ordered);

	return failed;
}

int main(void)
{
	return check_gauge_bookkeeping_smoke() != 0;
}
