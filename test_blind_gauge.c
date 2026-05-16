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

static int check_stats_common(const struct hk_blind_gauge_stats *stats, int32_t n_chr, int32_t n_flipped)
{
	int failed = 0;
	failed |= check_i32("stats n_chr", stats->n_chr, n_chr);
	failed |= check_i32("stats n_flipped", stats->n_flipped, n_flipped);
	failed |= check_true("stats no finite", isfinite(stats->mismatch_no_flip));
	failed |= check_true("stats flip finite", isfinite(stats->mismatch_flip));
	failed |= check_true("stats chosen finite", isfinite(stats->mismatch_chosen));
	failed |= check_true("stats chosen <= no", stats->mismatch_chosen <= stats->mismatch_no_flip + 1e-9);
	failed |= check_true("stats chosen <= flip", stats->mismatch_chosen <= stats->mismatch_flip + 1e-9);
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

static int check_no_flip_needed(void)
{
	int32_t chr_by_haploid[4] = {0, 0, 1, 1};
	fvec3_t prev[8], cur[8], expected[8];
	uint8_t flipped[2] = {9, 9};
	struct hk_blind_gauge_stats stats;
	int failed = 0;

	set_prev_coords(prev);
	copy_coords(cur, prev, 8);
	copy_coords(expected, cur, 8);
	failed |= check_i32("no flip ret", hk_blind_temporal_gauge_stabilize_coords(4, chr_by_haploid, 2, prev, cur, flipped, &stats), 0);
	failed |= check_u8("no flip chr0", flipped[0], 0);
	failed |= check_u8("no flip chr1", flipped[1], 0);
	failed |= check_coords_close("no flip coords", cur, expected, 8);
	failed |= check_stats_common(&stats, 2, 0);
	failed |= check_close("no flip chosen", stats.mismatch_chosen, 0.0);
	return failed;
}

static int check_one_chromosome_flipped(void)
{
	int32_t chr_by_haploid[4] = {0, 0, 1, 1};
	fvec3_t prev[8], cur[8];
	uint8_t flipped[2] = {0, 0};
	struct hk_blind_gauge_stats stats;
	int failed = 0;

	set_prev_coords(prev);
	copy_coords(cur, prev, 8);
	swap_pair_copies(cur, 0);
	swap_pair_copies(cur, 1);
	failed |= check_i32("one flip ret", hk_blind_temporal_gauge_stabilize_coords(4, chr_by_haploid, 2, prev, cur, flipped, &stats), 0);
	failed |= check_u8("one flip chr0", flipped[0], 1);
	failed |= check_u8("one flip chr1", flipped[1], 0);
	failed |= check_coords_close("one flip coords", cur, prev, 8);
	failed |= check_stats_common(&stats, 2, 1);
	failed |= check_close("one flip chosen", stats.mismatch_chosen, 0.0);
	return failed;
}

static int check_multiple_chromosomes_flipped(void)
{
	int32_t chr_by_haploid[4] = {0, 0, 1, 1};
	fvec3_t prev[8], cur[8];
	uint8_t flipped[2] = {0, 0};
	struct hk_blind_gauge_stats stats;
	int failed = 0, i;

	set_prev_coords(prev);
	copy_coords(cur, prev, 8);
	for (i = 0; i < 4; ++i)
		swap_pair_copies(cur, i);
	failed |= check_i32("multi flip ret", hk_blind_temporal_gauge_stabilize_coords(4, chr_by_haploid, 2, prev, cur, flipped, &stats), 0);
	failed |= check_u8("multi flip chr0", flipped[0], 1);
	failed |= check_u8("multi flip chr1", flipped[1], 1);
	failed |= check_coords_close("multi flip coords", cur, prev, 8);
	failed |= check_stats_common(&stats, 2, 2);
	failed |= check_close("multi flip chosen", stats.mismatch_chosen, 0.0);
	return failed;
}

static int check_tie_policy_no_flip(void)
{
	int32_t chr_by_haploid[1] = {0};
	fvec3_t prev[2], cur[2], expected[2];
	uint8_t flipped[1] = {7};
	struct hk_blind_gauge_stats stats;
	int failed = 0;

	set_pair(prev, 0, -1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f);
	set_pair(cur, 0, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
	copy_coords(expected, cur, 2);
	failed |= check_i32("tie ret", hk_blind_temporal_gauge_stabilize_coords(1, chr_by_haploid, 1, prev, cur, flipped, &stats), 0);
	failed |= check_u8("tie flipped", flipped[0], 0);
	failed |= check_coords_close("tie coords unchanged", cur, expected, 2);
	failed |= check_stats_common(&stats, 1, 0);
	failed |= check_close("tie no vs flip", stats.mismatch_no_flip, stats.mismatch_flip);
	failed |= check_close("tie chosen", stats.mismatch_chosen, stats.mismatch_no_flip);
	return failed;
}

static int check_no_beads(void)
{
	uint8_t flipped[3] = {5, 5, 5};
	struct hk_blind_gauge_stats stats;
	int failed = 0;

	failed |= check_i32("empty ret", hk_blind_temporal_gauge_stabilize_coords(0, 0, 3, 0, 0, flipped, &stats), 0);
	failed |= check_u8("empty chr0", flipped[0], 0);
	failed |= check_u8("empty chr1", flipped[1], 0);
	failed |= check_u8("empty chr2", flipped[2], 0);
	failed |= check_i32("empty stats n_chr", stats.n_chr, 3);
	failed |= check_i32("empty stats n_flipped", stats.n_flipped, 0);
	failed |= check_close("empty no", stats.mismatch_no_flip, 0.0);
	failed |= check_close("empty flip", stats.mismatch_flip, 0.0);
	failed |= check_close("empty chosen", stats.mismatch_chosen, 0.0);
	return failed;
}

static int check_p4_permutations(void)
{
	float in[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	float out[HK_BLIND_N_STATE];
	float expected_no[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	float expected_e0[HK_BLIND_N_STATE] = {0.3f, 0.4f, 0.1f, 0.2f};
	float expected_e1[HK_BLIND_N_STATE] = {0.2f, 0.1f, 0.4f, 0.3f};
	float expected_both[HK_BLIND_N_STATE] = {0.4f, 0.3f, 0.2f, 0.1f};
	float inplace[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	int failed = 0;

	hk_blind_p4_apply_endpoint_flips(in, 0, 0, out);
	failed |= check_p4_close("p4 no flip", out, expected_no);
	hk_blind_p4_apply_endpoint_flips(in, 1, 0, out);
	failed |= check_p4_close("p4 endpoint0", out, expected_e0);
	hk_blind_p4_apply_endpoint_flips(in, 0, 1, out);
	failed |= check_p4_close("p4 endpoint1", out, expected_e1);
	hk_blind_p4_apply_endpoint_flips(in, 1, 1, out);
	failed |= check_p4_close("p4 both", out, expected_both);
	hk_blind_p4_apply_endpoint_flips(inplace, 1, 0, inplace);
	failed |= check_p4_close("p4 inplace", inplace, expected_e0);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_no_flip_needed();
	failed |= check_one_chromosome_flipped();
	failed |= check_multiple_chromosomes_flipped();
	failed |= check_tie_policy_no_flip();
	failed |= check_no_beads();
	failed |= check_p4_permutations();
	return failed != 0;
}
