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

static int check_greater(const char *label, float a, float b)
{
	if (!(a > b)) {
		fprintf(stderr, "%s: expected %.8g > %.8g\n", label, a, b);
		return 1;
	}
	return 0;
}

static int check_range(const char *label, float x, float lo, float hi)
{
	if (!(x >= lo && x <= hi)) {
		fprintf(stderr, "%s: expected %.8g in [%.8g, %.8g]\n", label, x, lo, hi);
		return 1;
	}
	return 0;
}

static int check_key(const char *label, const struct hk_blind_bpair_key *key, int32_t bid0, int32_t bid1)
{
	int failed = 0;
	char buf[64];
	snprintf(buf, sizeof(buf), "%s bid0", label);
	failed |= check_i32(buf, key->bid[0], bid0);
	snprintf(buf, sizeof(buf), "%s bid1", label);
	failed |= check_i32(buf, key->bid[1], bid1);
	return failed;
}

static int check_p4(const char *label, const float got[HK_BLIND_N_STATE], const float expected[HK_BLIND_N_STATE])
{
	int failed = 0, i;
	char buf[64];
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		snprintf(buf, sizeof(buf), "%s p4[%d]", label, i);
		failed |= check_close(buf, got[i], expected[i]);
	}
	return failed;
}

static int check_raw_p4_inheritance(const char *label, const float raw[HK_BLIND_N_STATE],
									const float canonical[HK_BLIND_N_STATE], uint8_t swapped)
{
	int failed = 0;
	char buf[64];
	if (swapped) {
		snprintf(buf, sizeof(buf), "%s p00", label);
		failed |= check_close(buf, raw[HK_BLIND_STATE_00], canonical[HK_BLIND_STATE_00]);
		snprintf(buf, sizeof(buf), "%s p01", label);
		failed |= check_close(buf, raw[HK_BLIND_STATE_01], canonical[HK_BLIND_STATE_10]);
		snprintf(buf, sizeof(buf), "%s p10", label);
		failed |= check_close(buf, raw[HK_BLIND_STATE_10], canonical[HK_BLIND_STATE_01]);
		snprintf(buf, sizeof(buf), "%s p11", label);
		failed |= check_close(buf, raw[HK_BLIND_STATE_11], canonical[HK_BLIND_STATE_11]);
	} else {
		failed |= check_p4(label, raw, canonical);
	}
	return failed;
}

static int check_finite_bpair(const struct hk_blind_bpair *p)
{
	int failed = 0, i;
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		if (!isfinite(p->p4[i])) {
			fprintf(stderr, "p4[%d] is not finite: %.8g\n", i, p->p4[i]);
			failed = 1;
		}
	}
	if (!isfinite(p->entropy)) {
		fprintf(stderr, "entropy is not finite: %.8g\n", p->entropy);
		failed = 1;
	}
	if (!isfinite(p->pmax)) {
		fprintf(stderr, "pmax is not finite: %.8g\n", p->pmax);
		failed = 1;
	}
	if (!isfinite(p->margin)) {
		fprintf(stderr, "margin is not finite: %.8g\n", p->margin);
		failed = 1;
	}
	if (!isfinite(p->rho_output)) {
		fprintf(stderr, "rho_output is not finite: %.8g\n", p->rho_output);
		failed = 1;
	}
	if (!isfinite(p->pU)) {
		fprintf(stderr, "pU is not finite: %.8g\n", p->pU);
		failed = 1;
	}
	return failed;
}

static int check_probability_sums(const struct hk_blind_bpair *p)
{
	float p_sum = 0.0f, out_sum;
	int failed = 0, i;
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		p_sum += p->p4[i];
	out_sum = p->rho_output * p_sum + p->pU;
	failed |= check_close("posterior p4 sum", p_sum, 1.0f);
	failed |= check_close("posterior output-with-U sum", out_sum, 1.0f);
	failed |= check_range("posterior entropy", p->entropy, 0.0f, logf(4.0f));
	failed |= check_range("posterior pmax", p->pmax, 0.25f, 1.0f);
	failed |= check_range("posterior margin", p->margin, 0.0f, 1.0f);
	failed |= check_range("posterior rho_output", p->rho_output, 0.0f, 1.0f);
	failed |= check_range("posterior pU", p->pU, 0.0f, 1.0f);
	return failed;
}

static void set_coord(fvec3_t x, float x0, float x1, float x2)
{
	x[0] = x0;
	x[1] = x1;
	x[2] = x2;
}

static void init_coords_for_p00_best(fvec3_t *coords, int32_t n)
{
	int32_t i;
	for (i = 0; i < n; ++i)
		set_coord(coords[i], 100.0f + 10.0f * i, 100.0f, 100.0f);

	set_coord(coords[hk_diploid_bid(2, HK_DIPLOID_COPY0)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(5, HK_DIPLOID_COPY0)], 1.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(5, HK_DIPLOID_COPY1)], 3.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(2, HK_DIPLOID_COPY1)], -0.25f, sqrtf(14.4375f), 0.0f);
}

static int check_smoke_path(void)
{
	struct hk_map *m;
	struct hk_bmap *b;
	struct hk_blind_pair raw[8];
	struct hk_blind_bpair_set *set;
	struct hk_fdg_conf conf;
	fvec3_t coords[12];
	float log_prior[HK_BLIND_N_STATE];
	float raw_p4[HK_BLIND_N_STATE];
	float first_forward[HK_BLIND_N_STATE] = {0.0f, 0.0f, 0.0f, 0.0f};
	float first_reverse[HK_BLIND_N_STATE] = {0.0f, 0.0f, 0.0f, 0.0f};
	int have_forward = 0, have_reverse = 0;
	int32_t n_forward = 0, n_reverse = 0;
	int failed = 0;
	int32_t i;

	hk_verbose = 0;
	m = hk_map_read("testdata/p9016_blind_smoke.pairs");
	if (m == 0) {
		fprintf(stderr, "failed to read smoke fixture\n");
		return 1;
	}
	failed |= check_i32("smoke fixture n_pairs", m->n_pairs, 8);
	if (m->n_pairs != 8) {
		hk_map_destroy(m);
		return 1;
	}

	b = hk_bmap_gen(m->d, m->n_pairs, m->pairs, 1000000, 1);
	if (b == 0) {
		fprintf(stderr, "failed to build smoke bmap\n");
		hk_map_destroy(m);
		return 1;
	}

	for (i = 0; i < m->n_pairs; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);
	set = hk_blind_bpair_set_build(b, m->n_pairs, raw);

	failed |= check_i32("smoke set n_raw", set->n_raw, 8);
	failed |= check_i32("smoke set n_bpairs", set->n_bpairs, 1);
	failed |= check_key("smoke canonical key", &set->bpairs[0].key, 2, 5);
	failed |= check_i32("smoke canonical n_raw", set->bpairs[0].n_raw, 8);

	hk_fdg_conf_init(&conf);
	hk_blind_init_uniform_log_prior(log_prior);
	init_coords_for_p00_best(coords, 12);
	hk_blind_bpair_set_update_posterior_from_coords(set, &conf, coords, 1.0f, 1.0f, 2.0f, log_prior, 1.0f);

	failed |= check_finite_bpair(&set->bpairs[0]);
	failed |= check_probability_sums(&set->bpairs[0]);
	failed |= check_greater("smoke p00 > p01", set->bpairs[0].p4[HK_BLIND_STATE_00], set->bpairs[0].p4[HK_BLIND_STATE_01]);
	failed |= check_greater("smoke p00 > p10", set->bpairs[0].p4[HK_BLIND_STATE_00], set->bpairs[0].p4[HK_BLIND_STATE_10]);
	failed |= check_greater("smoke p00 > p11", set->bpairs[0].p4[HK_BLIND_STATE_00], set->bpairs[0].p4[HK_BLIND_STATE_11]);
	failed |= check_greater("smoke p01 > p10", set->bpairs[0].p4[HK_BLIND_STATE_01], set->bpairs[0].p4[HK_BLIND_STATE_10]);

	for (i = 0; i < set->n_raw; ++i) {
		int32_t bid[2];
		hk_blind_pair_to_bids(b, &raw[i], bid);
		failed |= check_i32("smoke raw bpair id", set->raw2binned[i].bpair_id, 0);
		if (bid[0] == 2 && bid[1] == 5) {
			++n_forward;
			failed |= check_u8("smoke forward swapped", set->raw2binned[i].swapped, 0);
		} else if (bid[0] == 5 && bid[1] == 2) {
			++n_reverse;
			failed |= check_u8("smoke reverse swapped", set->raw2binned[i].swapped, 1);
		} else {
			fprintf(stderr, "unexpected smoke raw bids: %d,%d\n", bid[0], bid[1]);
			failed = 1;
		}

		hk_blind_p4_to_raw_order(set->bpairs[0].p4, set->raw2binned[i].swapped, raw_p4);
		failed |= check_raw_p4_inheritance("smoke raw inherited p4", raw_p4, set->bpairs[0].p4, set->raw2binned[i].swapped);
		if (set->raw2binned[i].swapped) {
			if (!have_reverse) {
				int j;
				for (j = 0; j < HK_BLIND_N_STATE; ++j)
					first_reverse[j] = raw_p4[j];
				have_reverse = 1;
			} else {
				failed |= check_p4("smoke same reverse phase variants", raw_p4, first_reverse);
			}
		} else {
			if (!have_forward) {
				int j;
				for (j = 0; j < HK_BLIND_N_STATE; ++j)
					first_forward[j] = raw_p4[j];
				have_forward = 1;
			} else {
				failed |= check_p4("smoke same forward phase variants", raw_p4, first_forward);
			}
		}
	}

	failed |= check_i32("smoke forward count", n_forward, 4);
	failed |= check_i32("smoke reverse count", n_reverse, 4);
	if (have_forward)
		failed |= check_p4("smoke forward equals canonical", first_forward, set->bpairs[0].p4);
	if (have_reverse)
		failed |= check_raw_p4_inheritance("smoke reverse swaps p01 p10", first_reverse, set->bpairs[0].p4, 1);

	hk_blind_bpair_set_destroy(set);
	hk_bmap_destroy(b);
	hk_map_destroy(m);
	return failed;
}

int main(void)
{
	return check_smoke_path() != 0;
}
