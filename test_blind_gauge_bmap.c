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

static int check_stats_close(const char *label, const struct hk_blind_gauge_stats *got, const struct hk_blind_gauge_stats *expected)
{
	int failed = 0;
	char buf[96];
	snprintf(buf, sizeof(buf), "%s n_chr", label);
	failed |= check_i32(buf, got->n_chr, expected->n_chr);
	snprintf(buf, sizeof(buf), "%s n_flipped", label);
	failed |= check_i32(buf, got->n_flipped, expected->n_flipped);
	snprintf(buf, sizeof(buf), "%s mismatch no", label);
	failed |= check_close(buf, got->mismatch_no_flip, expected->mismatch_no_flip);
	snprintf(buf, sizeof(buf), "%s mismatch flip", label);
	failed |= check_close(buf, got->mismatch_flip, expected->mismatch_flip);
	snprintf(buf, sizeof(buf), "%s mismatch chosen", label);
	failed |= check_close(buf, got->mismatch_chosen, expected->mismatch_chosen);
	return failed;
}

static void set_prev_coords(fvec3_t coords[8])
{
	set_pair(coords, 0, 0.0f, 0.0f, 0.0f, 10.0f, 0.0f, 0.0f);
	set_pair(coords, 1, 0.0f, 1.0f, 0.0f, 10.0f, 1.0f, 0.0f);
	set_pair(coords, 2, 100.0f, 0.0f, 0.0f, 110.0f, 0.0f, 0.0f);
	set_pair(coords, 3, 100.0f, 1.0f, 0.0f, 110.0f, 1.0f, 0.0f);
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

static void set_gauge_bmap(struct hk_bmap *bmap, struct hk_bead beads[4])
{
	int32_t i;
	for (i = 0; i < 4; ++i) {
		beads[i].chr = i < 2? 0 : 1;
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

static void set_empty_bmap(struct hk_bmap *bmap)
{
	bmap->n_beads = 0;
	bmap->n_pairs = 0;
	bmap->unit = 1.0f;
	bmap->d = 0;
	bmap->beads = 0;
	bmap->offcnt = 0;
	bmap->pairs = 0;
	bmap->x = 0;
	bmap->feat = 0;
	bmap->cpg = 0;
	bmap->gc_bias = 0;
	bmap->gc_corrected = 0;
}

static void set_p4_bmap(struct hk_bmap *bmap, struct hk_bead beads[4])
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

static int check_bmap_wrapper_matches_low_level(void)
{
	int32_t chr_by_haploid[4] = {0, 0, 1, 1};
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t prev[8], cur_low[8], cur_bmap[8];
	uint8_t flipped_low[2] = {0, 0}, flipped_bmap[2] = {0, 0};
	struct hk_blind_gauge_stats stats_low, stats_bmap;
	int failed = 0;

	set_gauge_bmap(&bmap, beads);
	set_prev_coords(prev);
	copy_coords(cur_low, prev, 8);
	copy_coords(cur_bmap, prev, 8);
	swap_pair_copies(cur_low, 0);
	swap_pair_copies(cur_low, 1);
	swap_pair_copies(cur_bmap, 0);
	swap_pair_copies(cur_bmap, 1);

	failed |= check_i32("low ret", hk_blind_temporal_gauge_stabilize_coords(4, chr_by_haploid, 2, prev, cur_low, flipped_low, &stats_low), 0);
	failed |= check_i32("bmap ret", hk_blind_temporal_gauge_stabilize_bmap(&bmap, prev, cur_bmap, flipped_bmap, &stats_bmap), 0);
	failed |= check_u8("bmap match flip0", flipped_bmap[0], flipped_low[0]);
	failed |= check_u8("bmap match flip1", flipped_bmap[1], flipped_low[1]);
	failed |= check_coords_close("bmap match coords", cur_bmap, cur_low, 8);
	failed |= check_stats_close("bmap match stats", &stats_bmap, &stats_low);
	return failed;
}

static int check_bmap_wrapper_null_outputs(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t prev[8], cur[8];
	uint8_t flipped[2] = {0, 0};
	struct hk_blind_gauge_stats stats;
	int failed = 0;

	set_gauge_bmap(&bmap, beads);
	set_prev_coords(prev);
	copy_coords(cur, prev, 8);
	swap_pair_copies(cur, 0);
	swap_pair_copies(cur, 1);
	failed |= check_i32("bmap null flipped ret", hk_blind_temporal_gauge_stabilize_bmap(&bmap, prev, cur, 0, &stats), 0);
	failed |= check_coords_close("bmap null flipped coords", cur, prev, 8);
	failed |= check_i32("bmap null flipped stats chr", stats.n_chr, 2);
	failed |= check_i32("bmap null flipped stats flips", stats.n_flipped, 1);

	copy_coords(cur, prev, 8);
	swap_pair_copies(cur, 0);
	swap_pair_copies(cur, 1);
	failed |= check_i32("bmap null stats ret", hk_blind_temporal_gauge_stabilize_bmap(&bmap, prev, cur, flipped, 0), 0);
	failed |= check_u8("bmap null stats flip0", flipped[0], 1);
	failed |= check_u8("bmap null stats flip1", flipped[1], 0);
	failed |= check_coords_close("bmap null stats coords", cur, prev, 8);
	return failed;
}

static int check_bmap_wrapper_no_beads(void)
{
	struct hk_bmap bmap;
	struct hk_blind_gauge_stats stats;
	int failed = 0;

	set_empty_bmap(&bmap);
	failed |= check_i32("empty bmap ret", hk_blind_temporal_gauge_stabilize_bmap(&bmap, 0, 0, 0, &stats), 0);
	failed |= check_i32("empty bmap stats chr", stats.n_chr, 0);
	failed |= check_i32("empty bmap stats flips", stats.n_flipped, 0);
	failed |= check_close("empty bmap no", stats.mismatch_no_flip, 0.0);
	failed |= check_close("empty bmap flip", stats.mismatch_flip, 0.0);
	failed |= check_close("empty bmap chosen", stats.mismatch_chosen, 0.0);
	return failed;
}

static int check_bmap_wrapper_uses_dict_chr_count(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	struct hk_sdict dict;
	fvec3_t prev[8], cur[8];
	uint8_t flipped[3] = {7, 7, 7};
	struct hk_blind_gauge_stats stats;
	int failed = 0;

	set_gauge_bmap(&bmap, beads);
	dict.n = 3;
	dict.m = 3;
	dict.name = 0;
	dict.len = 0;
	dict.h = 0;
	bmap.d = &dict;
	set_prev_coords(prev);
	copy_coords(cur, prev, 8);
	failed |= check_i32("dict n_chr ret", hk_blind_temporal_gauge_stabilize_bmap(&bmap, prev, cur, flipped, &stats), 0);
	failed |= check_i32("dict n_chr stats", stats.n_chr, 3);
	failed |= check_u8("dict n_chr flip0", flipped[0], 0);
	failed |= check_u8("dict n_chr flip1", flipped[1], 0);
	failed |= check_u8("dict n_chr flip2", flipped[2], 0);
	return failed;
}

static int run_p4_case(const char *label, int32_t bid0, int32_t bid1, const uint8_t chr_flipped[2],
					   const float expected_p4[HK_BLIND_N_STATE])
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	struct hk_blind_bpair bpair;
	struct hk_blind_bpair_set set;
	float p4[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	float entropy0, pmax0, margin0, rho0, pU0;
	int failed = 0;

	set_p4_bmap(&bmap, beads);
	set_bpair(&bpair, bid0, bid1, p4);
	entropy0 = bpair.entropy;
	pmax0 = bpair.pmax;
	margin0 = bpair.margin;
	rho0 = bpair.rho_output;
	pU0 = bpair.pU;
	set.bpairs = &bpair;
	set.n_bpairs = 1;
	set.raw2binned = 0;
	set.n_raw = 0;

	failed |= check_i32(label, hk_blind_bpair_set_apply_chr_flips(&set, &bmap, chr_flipped, 2), 0);
	failed |= check_p4_close(label, bpair.p4, expected_p4);
	failed |= check_close("p4 entropy invariant", bpair.entropy, entropy0);
	failed |= check_close("p4 pmax invariant", bpair.pmax, pmax0);
	failed |= check_close("p4 margin invariant", bpair.margin, margin0);
	failed |= check_close("p4 rho invariant", bpair.rho_output, rho0);
	failed |= check_close("p4 pU invariant", bpair.pU, pU0);
	return failed;
}

static int check_bpair_p4_chr_flips(void)
{
	uint8_t no_flip[2] = {0, 0};
	uint8_t endpoint0_flip[2] = {1, 0};
	uint8_t endpoint1_flip[2] = {0, 1};
	uint8_t both_flip[2] = {1, 1};
	float expected_no[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	float expected_e0[HK_BLIND_N_STATE] = {0.3f, 0.4f, 0.1f, 0.2f};
	float expected_e1[HK_BLIND_N_STATE] = {0.2f, 0.1f, 0.4f, 0.3f};
	float expected_both[HK_BLIND_N_STATE] = {0.4f, 0.3f, 0.2f, 0.1f};
	int failed = 0;

	failed |= run_p4_case("p4 no chr flip", 0, 1, no_flip, expected_no);
	failed |= run_p4_case("p4 endpoint0 chr flip", 0, 1, endpoint0_flip, expected_e0);
	failed |= run_p4_case("p4 endpoint1 chr flip", 0, 1, endpoint1_flip, expected_e1);
	failed |= run_p4_case("p4 both chr flip", 0, 1, both_flip, expected_both);
	failed |= run_p4_case("p4 same chr both endpoints", 0, 3, endpoint0_flip, expected_both);
	return failed;
}

static int check_raw_order_composition(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	struct hk_blind_bpair bpair;
	struct hk_blind_bpair_set set;
	uint8_t endpoint0_flip[2] = {1, 0};
	float p4[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	float raw[HK_BLIND_N_STATE];
	float expected_raw[HK_BLIND_N_STATE] = {0.3f, 0.1f, 0.4f, 0.2f};
	int failed = 0;

	set_p4_bmap(&bmap, beads);
	set_bpair(&bpair, 0, 1, p4);
	set.bpairs = &bpair;
	set.n_bpairs = 1;
	set.raw2binned = 0;
	set.n_raw = 0;
	failed |= check_i32("raw composition flip ret", hk_blind_bpair_set_apply_chr_flips(&set, &bmap, endpoint0_flip, 2), 0);
	hk_blind_p4_to_raw_order(bpair.p4, 1, raw);
	failed |= check_p4_close("raw composition", raw, expected_raw);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_bmap_wrapper_matches_low_level();
	failed |= check_bmap_wrapper_null_outputs();
	failed |= check_bmap_wrapper_no_beads();
	failed |= check_bmap_wrapper_uses_dict_chr_count();
	failed |= check_bpair_p4_chr_flips();
	failed |= check_raw_order_composition();
	return failed != 0;
}
