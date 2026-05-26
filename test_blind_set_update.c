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

static int check_u8(const char *label, uint8_t got, uint8_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %u, expected %u\n", label, got, expected);
		return 1;
	}
	return 0;
}

static int check_true(const char *label, int ok)
{
	if (!ok) {
		fprintf(stderr, "%s: false\n", label);
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

static int check_less(const char *label, float a, float b)
{
	if (!(a < b)) {
		fprintf(stderr, "%s: expected %.8g < %.8g\n", label, a, b);
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

static int check_sums(const char *label, const struct hk_blind_bpair *p)
{
	float p_sum = 0.0f, out_sum;
	int i, failed = 0;
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		p_sum += p->p4[i];
	out_sum = p->rho_output * p_sum + p->pU;
	failed |= check_close(label, p_sum, 1.0f);
	failed |= check_close("output-with-U sum", out_sum, 1.0f);
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

static struct hk_blind_pair make_blind_pair(int32_t chr0, int32_t pos0, int32_t chr1, int32_t pos1)
{
	struct hk_blind_pair p;
	p.chr[0] = chr0;
	p.chr[1] = chr1;
	p.pos[0] = pos0;
	p.pos[1] = pos1;
	p.strand[0] = 1;
	p.strand[1] = -1;
	return p;
}

static void set_coord(fvec3_t x, float x0, float x1, float x2)
{
	x[0] = x0;
	x[1] = x1;
	x[2] = x2;
}

static void init_far_coords(fvec3_t *coords, int32_t n)
{
	int32_t i;
	for (i = 0; i < n; ++i)
		set_coord(coords[i], 1000.0f + 10.0f * i, 1000.0f, 1000.0f);
}

static int check_updated_uncertainty(const char *label, const struct hk_blind_bpair *p)
{
	int failed = 0;
	failed |= check_finite_bpair(p);
	failed |= check_sums(label, p);
	failed |= check_less("updated entropy", p->entropy, logf(4.0f));
	failed |= check_greater("updated pmax", p->pmax, 0.25f);
	failed |= check_greater("updated margin", p->margin, 0.0f);
	failed |= check_greater("updated rho_output", p->rho_output, 0.0f);
	failed |= check_less("updated pU", p->pU, 1.0f);
	return failed;
}

static int check_bpair_set_update(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[3];
	struct hk_blind_bpair_set *set;
	struct hk_fdg_conf conf;
	fvec3_t coords[16];
	float log_prior[HK_BLIND_N_STATE];
	int failed = 0;

	raw[0] = make_blind_pair(0, 2000000, 0, 5000000);
	raw[1] = make_blind_pair(0, 5000000, 0, 2000000);
	raw[2] = make_blind_pair(0, 3000000, 0, 7000000);
	set = hk_blind_bpair_set_build(b, 3, raw);

	failed |= check_i32("update n_bpairs", set->n_bpairs, 2);
	failed |= check_key("update key0", &set->bpairs[0].key, 2, 5);
	failed |= check_key("update key1", &set->bpairs[1].key, 3, 7);
	failed |= check_i32("update key0 n_raw", set->bpairs[0].n_raw, 2);
	failed |= check_i32("update key1 n_raw", set->bpairs[1].n_raw, 1);
	failed |= check_i32("update raw0 bpair", set->raw2binned[0].bpair_id, 0);
	failed |= check_i32("update raw1 bpair", set->raw2binned[1].bpair_id, 0);
	failed |= check_i32("update raw2 bpair", set->raw2binned[2].bpair_id, 1);
	failed |= check_u8("update raw0 swapped", set->raw2binned[0].swapped, 0);
	failed |= check_u8("update raw1 swapped", set->raw2binned[1].swapped, 1);
	failed |= check_u8("update raw2 swapped", set->raw2binned[2].swapped, 0);

	hk_fdg_conf_init(&conf);
	hk_blind_init_uniform_log_prior(log_prior);
	init_far_coords(coords, 16);

	set_coord(coords[hk_diploid_bid(2, HK_DIPLOID_COPY0)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(2, HK_DIPLOID_COPY1)], 100.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(5, HK_DIPLOID_COPY0)], 1.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(5, HK_DIPLOID_COPY1)], 200.0f, 0.0f, 0.0f);

	set_coord(coords[hk_diploid_bid(3, HK_DIPLOID_COPY0)], 300.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(3, HK_DIPLOID_COPY1)], 400.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(7, HK_DIPLOID_COPY0)], 500.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(7, HK_DIPLOID_COPY1)], 401.0f, 0.0f, 0.0f);

	hk_blind_bpair_set_update_posterior_from_coords(set, &conf, coords, 1.0f, 1.0f, 2.0f, log_prior, 1.0f);

	failed |= check_updated_uncertainty("update key0 prob sum", &set->bpairs[0]);
	failed |= check_greater("key0 p00 > p01", set->bpairs[0].p4[HK_BLIND_STATE_00], set->bpairs[0].p4[HK_BLIND_STATE_01]);
	failed |= check_greater("key0 p00 > p10", set->bpairs[0].p4[HK_BLIND_STATE_00], set->bpairs[0].p4[HK_BLIND_STATE_10]);
	failed |= check_greater("key0 p00 > p11", set->bpairs[0].p4[HK_BLIND_STATE_00], set->bpairs[0].p4[HK_BLIND_STATE_11]);

	failed |= check_updated_uncertainty("update key1 prob sum", &set->bpairs[1]);
	failed |= check_greater("key1 p11 > p00", set->bpairs[1].p4[HK_BLIND_STATE_11], set->bpairs[1].p4[HK_BLIND_STATE_00]);
	failed |= check_greater("key1 p11 > p01", set->bpairs[1].p4[HK_BLIND_STATE_11], set->bpairs[1].p4[HK_BLIND_STATE_01]);
	failed |= check_greater("key1 p11 > p10", set->bpairs[1].p4[HK_BLIND_STATE_11], set->bpairs[1].p4[HK_BLIND_STATE_10]);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_posterior_independent_of_bpair_params(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[1];
	struct hk_blind_bpair_set *set;
	struct hk_fdg_conf conf;
	fvec3_t coords[12];
	float log_prior[HK_BLIND_N_STATE];
	float first_p4[HK_BLIND_N_STATE];
	float first_entropy, first_pmax, first_margin, first_rho_output, first_pU;
	int failed = 0, i;

	raw[0] = make_blind_pair(0, 2000000, 0, 5000000);
	set = hk_blind_bpair_set_build(b, 1, raw);
	failed |= check_i32("posterior params n_bpairs", set->n_bpairs, 1);

	hk_fdg_conf_init(&conf);
	hk_blind_init_uniform_log_prior(log_prior);
	init_far_coords(coords, 12);
	set_coord(coords[hk_diploid_bid(2, HK_DIPLOID_COPY0)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(2, HK_DIPLOID_COPY1)], 100.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(5, HK_DIPLOID_COPY0)], 1.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(5, HK_DIPLOID_COPY1)], 200.0f, 0.0f, 0.0f);

	set->bpairs[0].base_d_scale = 0.25f;
	set->bpairs[0].base_k = 99.0f;
	hk_blind_bpair_set_update_posterior_from_coords(set, &conf, coords, 1.0f, 1.0f, 2.0f, log_prior, 1.0f);
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		first_p4[i] = set->bpairs[0].p4[i];
	first_entropy = set->bpairs[0].entropy;
	first_pmax = set->bpairs[0].pmax;
	first_margin = set->bpairs[0].margin;
	first_rho_output = set->bpairs[0].rho_output;
	first_pU = set->bpairs[0].pU;

	set->bpairs[0].base_d_scale = 4.0f;
	set->bpairs[0].base_k = 0.125f;
	hk_blind_bpair_set_update_posterior_from_coords(set, &conf, coords, 1.0f, 1.0f, 2.0f, log_prior, 1.0f);
	failed |= check_p4("posterior params p4 unchanged", set->bpairs[0].p4, first_p4);
	failed |= check_close("posterior params entropy unchanged", set->bpairs[0].entropy, first_entropy);
	failed |= check_close("posterior params pmax unchanged", set->bpairs[0].pmax, first_pmax);
	failed |= check_close("posterior params margin unchanged", set->bpairs[0].margin, first_margin);
	failed |= check_close("posterior params rho unchanged", set->bpairs[0].rho_output, first_rho_output);
	failed |= check_close("posterior params pU unchanged", set->bpairs[0].pU, first_pU);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_estep_logdist2_breaks_fdg_flat_ties(void)
{
	struct hk_fdg_conf conf;
	fvec3_t coords[4];
	float log_prior[HK_BLIND_N_STATE];
	float energy[HK_BLIND_N_STATE], p4[HK_BLIND_N_STATE];
	float entropy, pmax, margin, rho_output, pU;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	hk_blind_init_uniform_log_prior(log_prior);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 10.0f, 0.0f);
	set_coord(coords[2], 0.75f, 0.0f, 0.0f);
	set_coord(coords[3], 1.25f, 0.0f, 0.0f);
	hk_blind_bpair_posterior_from_coords(&conf, 0, 1, coords, 1.0f, 1.0f, 1.0f,
										 log_prior, 1.0f, energy, p4,
										 &entropy, &pmax, &margin, &rho_output, &pU);
	failed |= check_close("flat tie p00 p01", p4[HK_BLIND_STATE_00], p4[HK_BLIND_STATE_01]);
	hk_blind_bpair_posterior_from_coords_score_mode(&conf, 0, 1, coords, 1.0f, 1.0f, 1.0f,
													log_prior, 1.0f,
													HK_BLIND_ESTEP_SCORE_LOGDIST2,
													energy, p4, &entropy, &pmax,
													&margin, &rho_output, &pU);
	failed |= check_greater("logdist2 nearer state higher", p4[HK_BLIND_STATE_00],
							p4[HK_BLIND_STATE_01]);
	return failed;
}

static int check_raw_order_helper(void)
{
	float canonical[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	float expected_same[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	float expected_swapped[HK_BLIND_N_STATE] = {0.1f, 0.3f, 0.2f, 0.4f};
	float raw[HK_BLIND_N_STATE];
	float inplace[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	int failed = 0;

	hk_blind_p4_to_raw_order(canonical, 0, raw);
	failed |= check_p4("raw order same", raw, expected_same);

	hk_blind_p4_to_raw_order(canonical, 1, raw);
	failed |= check_p4("raw order swapped", raw, expected_swapped);
	hk_blind_p4_to_raw_order(inplace, 1, inplace);
	failed |= check_p4("raw order swapped inplace", inplace, expected_swapped);
	failed |= check_i32("state 00 same", hk_blind_raw_state_to_canonical_state(HK_BLIND_STATE_00, 0), HK_BLIND_STATE_00);
	failed |= check_i32("state 01 same", hk_blind_raw_state_to_canonical_state(HK_BLIND_STATE_01, 0), HK_BLIND_STATE_01);
	failed |= check_i32("state 10 same", hk_blind_raw_state_to_canonical_state(HK_BLIND_STATE_10, 0), HK_BLIND_STATE_10);
	failed |= check_i32("state 11 same", hk_blind_raw_state_to_canonical_state(HK_BLIND_STATE_11, 0), HK_BLIND_STATE_11);
	failed |= check_i32("state 00 swapped", hk_blind_raw_state_to_canonical_state(HK_BLIND_STATE_00, 1), HK_BLIND_STATE_00);
	failed |= check_i32("state 01 swapped", hk_blind_raw_state_to_canonical_state(HK_BLIND_STATE_01, 1), HK_BLIND_STATE_10);
	failed |= check_i32("state 10 swapped", hk_blind_raw_state_to_canonical_state(HK_BLIND_STATE_10, 1), HK_BLIND_STATE_01);
	failed |= check_i32("state 11 swapped", hk_blind_raw_state_to_canonical_state(HK_BLIND_STATE_11, 1), HK_BLIND_STATE_11);
	return failed;
}

static int check_raw2binned_inheritance(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[2];
	struct hk_blind_bpair_set *set;
	float canonical[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	float expected_same[HK_BLIND_N_STATE] = {0.1f, 0.2f, 0.3f, 0.4f};
	float expected_swapped[HK_BLIND_N_STATE] = {0.1f, 0.3f, 0.2f, 0.4f};
	float raw_p4[HK_BLIND_N_STATE];
	int failed = 0, i;

	raw[0] = make_blind_pair(0, 2000000, 0, 5000000);
	raw[1] = make_blind_pair(0, 5000000, 0, 2000000);
	set = hk_blind_bpair_set_build(b, 2, raw);

	failed |= check_i32("inherit n_bpairs", set->n_bpairs, 1);
	failed |= check_key("inherit key", &set->bpairs[0].key, 2, 5);
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		set->bpairs[0].p4[i] = canonical[i];

	hk_blind_p4_to_raw_order(set->bpairs[set->raw2binned[0].bpair_id].p4, set->raw2binned[0].swapped, raw_p4);
	failed |= check_p4("inherit raw0", raw_p4, expected_same);
	hk_blind_p4_to_raw_order(set->bpairs[set->raw2binned[1].bpair_id].p4, set->raw2binned[1].swapped, raw_p4);
	failed |= check_p4("inherit raw1", raw_p4, expected_swapped);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_raw_phase_locks_record_raw_state(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[2];
	struct hk_pair pairs[2];
	struct hk_blind_bpair_set *set;
	int failed = 0;

	memset(pairs, 0, sizeof(pairs));
	raw[0] = make_blind_pair(0, 0, 0, 1000000);
	raw[1] = make_blind_pair(0, 1000000, 0, 0);
	pairs[0].chr = (uint64_t)0 << 32 | 0u;
	pairs[0].pos = (uint64_t)0 << 32 | 1000000u;
	pairs[0].phase[0] = 0;
	pairs[0].phase[1] = 1;
	pairs[1].chr = (uint64_t)0 << 32 | 0u;
	pairs[1].pos = (uint64_t)1000000 << 32 | 0u;
	pairs[1].phase[0] = 0;
	pairs[1].phase[1] = 1;

	set = hk_blind_bpair_set_build(b, 2, raw);
	failed |= check_true("raw lock set exists", set != 0);
	if (set == 0)
		return 1;
	failed |= check_i32("raw lock init0", set->raw_locked_state[0], -1);
	failed |= check_i32("raw lock init1", set->raw_locked_state[1], -1);
	failed |= check_i32("apply raw locks for raw states",
						hk_blind_bpair_set_apply_raw_phase_locks(set, pairs, 2, 100,
																 17, &set->phase_lock_diag), 0);
	failed |= check_i32("raw lock forward 01 canonical", set->raw_locked_state[0],
						HK_BLIND_STATE_01);
	failed |= check_i32("raw lock swapped 01 canonical", set->raw_locked_state[1],
						HK_BLIND_STATE_10);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_imputed_p4_top_locks_record_raw_state(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[3];
	struct hk_pair pairs[3];
	struct hk_blind_bpair_set *set;
	int failed = 0;

	memset(pairs, 0, sizeof(pairs));
	raw[0] = make_blind_pair(0, 0, 0, 1000000);
	raw[1] = make_blind_pair(0, 1000000, 0, 0);
	raw[2] = make_blind_pair(0, 2000000, 0, 3000000);

	pairs[0]._.p4[HK_BLIND_STATE_00] = 0.01f;
	pairs[0]._.p4[HK_BLIND_STATE_01] = 0.80f;
	pairs[0]._.p4[HK_BLIND_STATE_10] = 0.10f;
	pairs[0]._.p4[HK_BLIND_STATE_11] = 0.09f;
	pairs[1]._.p4[HK_BLIND_STATE_00] = 0.01f;
	pairs[1]._.p4[HK_BLIND_STATE_01] = 0.80f;
	pairs[1]._.p4[HK_BLIND_STATE_10] = 0.10f;
	pairs[1]._.p4[HK_BLIND_STATE_11] = 0.09f;
	pairs[2]._.p4[HK_BLIND_STATE_00] = 0.30f;
	pairs[2]._.p4[HK_BLIND_STATE_01] = 0.30f;
	pairs[2]._.p4[HK_BLIND_STATE_10] = 0.20f;
	pairs[2]._.p4[HK_BLIND_STATE_11] = 0.20f;

	set = hk_blind_bpair_set_build(b, 3, raw);
	failed |= check_true("imputed p4 lock set exists", set != 0);
	if (set == 0)
		return 1;
	failed |= check_i32("apply imputed p4 top locks",
						hk_blind_bpair_set_apply_imputed_p4_top_locks(set, pairs, 3, 100,
																	 17, 0.75f,
																	 &set->phase_lock_diag), 0);
	failed |= check_i32("imputed p4 lock forward 01 canonical", set->raw_locked_state[0],
						HK_BLIND_STATE_01);
	failed |= check_i32("imputed p4 lock swapped 01 canonical", set->raw_locked_state[1],
						HK_BLIND_STATE_10);
	failed |= check_i32("imputed p4 lock low confidence skip", set->raw_locked_state[2],
						HK_BLIND_RAW_LOCKED_STATE_SKIP);
	failed |= check_i32("imputed p4 top selected raw", (int32_t)set->phase_lock_diag.n_locked_raw, 2);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_same_bin_posterior_is_unknown(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[1];
	struct hk_blind_bpair_set *set;
	struct hk_fdg_conf conf;
	fvec3_t coords[8];
	float log_prior[HK_BLIND_N_STATE];
	int failed = 0;

	raw[0] = make_blind_pair(0, 3000000, 0, 3000100);
	set = hk_blind_bpair_set_build(b, 1, raw);
	failed |= check_i32("same-bin unknown n_bpairs", set->n_bpairs, 1);
	failed |= check_key("same-bin unknown key", &set->bpairs[0].key, 3, 3);
	hk_fdg_conf_init(&conf);
	hk_blind_init_uniform_log_prior(log_prior);
	init_far_coords(coords, 8);
	set_coord(coords[hk_diploid_bid(3, HK_DIPLOID_COPY0)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(3, HK_DIPLOID_COPY1)], 100.0f, 0.0f, 0.0f);
	set->bpairs[0].p4[HK_BLIND_STATE_00] = 0.9f;
	set->bpairs[0].p4[HK_BLIND_STATE_01] = 0.05f;
	set->bpairs[0].p4[HK_BLIND_STATE_10] = 0.03f;
	set->bpairs[0].p4[HK_BLIND_STATE_11] = 0.02f;
	hk_blind_bpair_set_update_posterior_from_coords(set, &conf, coords, 1.0f, 1.0f, 2.0f, log_prior, 1.0f);
	failed |= check_close("same-bin p00 unknown", set->bpairs[0].p4[HK_BLIND_STATE_00], 0.25f);
	failed |= check_close("same-bin p01 unknown", set->bpairs[0].p4[HK_BLIND_STATE_01], 0.25f);
	failed |= check_close("same-bin p10 unknown", set->bpairs[0].p4[HK_BLIND_STATE_10], 0.25f);
	failed |= check_close("same-bin p11 unknown", set->bpairs[0].p4[HK_BLIND_STATE_11], 0.25f);
	failed |= check_close("same-bin pU unknown", set->bpairs[0].pU, 1.0f);
	failed |= check_close("same-bin rho unknown", set->bpairs[0].rho_output, 0.0f);
	hk_blind_bpair_set_update_posterior_from_coords_params(set, &conf, coords, 1.0f, log_prior, 1.0f);
	failed |= check_close("same-bin params p00 unknown", set->bpairs[0].p4[HK_BLIND_STATE_00], 0.25f);
	failed |= check_close("same-bin params pU unknown", set->bpairs[0].pU, 1.0f);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_phase_lock_skips_posterior_update(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[2];
	struct hk_pair pairs[2];
	struct hk_blind_bpair_set *set;
	struct hk_fdg_conf conf;
	struct hk_blind_phase_lock_diag diag;
	fvec3_t coords[12];
	float log_prior[HK_BLIND_N_STATE];
	int failed = 0;

	raw[0] = make_blind_pair(0, 2000000, 0, 5000000);
	raw[1] = make_blind_pair(0, 5000000, 0, 2000000);
	memset(pairs, 0, sizeof(pairs));
	pairs[0].phase[0] = 0;
	pairs[0].phase[1] = 1;
	pairs[1].phase[0] = -1;
	pairs[1].phase[1] = -1;
	set = hk_blind_bpair_set_build(b, 2, raw);
	hk_fdg_conf_init(&conf);
	hk_blind_init_uniform_log_prior(log_prior);
	init_far_coords(coords, 12);
	set_coord(coords[hk_diploid_bid(2, HK_DIPLOID_COPY0)], 100.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(2, HK_DIPLOID_COPY1)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(5, HK_DIPLOID_COPY0)], 1.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(5, HK_DIPLOID_COPY1)], 200.0f, 0.0f, 0.0f);

	failed |= check_i32("phase lock apply", hk_blind_bpair_set_apply_raw_phase_locks(set, pairs, 2, 100, 17, &diag), 0);
	failed |= check_i32("phase lock bpair", diag.n_locked_bpair, 1);
	failed |= check_i32("phase lock raw", (int32_t)diag.n_locked_raw, 1);
	failed |= check_i32("phase lock state", set->bpairs[0].locked_state, HK_BLIND_STATE_01);
	hk_blind_bpair_set_update_posterior_from_coords(set, &conf, coords, 1.0f, 1.0f, 2.0f, log_prior, 1.0f);
	failed |= check_close("locked p00", set->bpairs[0].p4[HK_BLIND_STATE_00], 0.0f);
	failed |= check_close("locked p01", set->bpairs[0].p4[HK_BLIND_STATE_01], 1.0f);
	failed |= check_close("locked p10", set->bpairs[0].p4[HK_BLIND_STATE_10], 0.0f);
	failed |= check_close("locked p11", set->bpairs[0].p4[HK_BLIND_STATE_11], 0.0f);
	failed |= check_close("locked pU", set->bpairs[0].pU, 0.0f);
	failed |= check_true("locked real log norm refreshed", fabsf(set->bpairs[0].real_log_norm) > 1e-5f);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_raw_locks_follow_chr_flips(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[2];
	struct hk_pair pairs[2];
	struct hk_blind_bpair_set *set;
	uint8_t chr_flipped[1] = { 1 };
	int failed = 0;

	raw[0] = make_blind_pair(0, 2000000, 0, 5000000);
	raw[1] = make_blind_pair(0, 5000000, 0, 2000000);
	memset(pairs, 0, sizeof(pairs));
	pairs[0].phase[0] = 0;
	pairs[0].phase[1] = 1;
	pairs[1].phase[0] = 1;
	pairs[1].phase[1] = 0;
	set = hk_blind_bpair_set_build(b, 2, raw);
	failed |= check_true("raw lock flip set exists", set != 0);
	if (set == 0)
		return 1;
	failed |= check_i32("raw lock flip apply",
						hk_blind_bpair_set_apply_raw_phase_locks(set, pairs, 2, 100,
																 17, &set->phase_lock_diag), 0);
	failed |= check_i32("raw lock flip forward before", set->raw_locked_state[0],
						HK_BLIND_STATE_01);
	failed |= check_i32("raw lock flip swapped before", set->raw_locked_state[1],
						HK_BLIND_STATE_01);
	failed |= check_i32("raw lock apply chr flips",
						hk_blind_bpair_set_apply_chr_flips(set, b, chr_flipped, 1), 0);
	failed |= check_i32("raw lock flip forward after", set->raw_locked_state[0],
						HK_BLIND_STATE_10);
	failed |= check_i32("raw lock flip swapped after", set->raw_locked_state[1],
						HK_BLIND_STATE_10);
	failed |= check_i32("bpair locked state after flip", set->bpairs[0].locked_state,
						HK_BLIND_STATE_10);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_fixed_p4_locks_follow_chr_flips(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[1];
	struct hk_pair pairs[1];
	struct hk_blind_bpair_set *set;
	uint8_t chr_flipped[1] = { 1 };
	int failed = 0;

	raw[0] = make_blind_pair(0, 2000000, 0, 5000000);
	memset(pairs, 0, sizeof(pairs));
	pairs[0]._.p4[HK_BLIND_STATE_00] = 0.01f;
	pairs[0]._.p4[HK_BLIND_STATE_01] = 0.80f;
	pairs[0]._.p4[HK_BLIND_STATE_10] = 0.10f;
	pairs[0]._.p4[HK_BLIND_STATE_11] = 0.09f;

	set = hk_blind_bpair_set_build(b, 1, raw);
	failed |= check_true("fixed p4 flip set exists", set != 0);
	if (set == 0)
		return 1;
	failed |= check_i32("fixed p4 lock apply",
						hk_blind_bpair_set_apply_imputed_p4_top_locks(set, pairs, 1, 100,
																	 17, 0.75f,
																	 &set->phase_lock_diag), 0);
	failed |= check_i32("fixed p4 raw state before", set->raw_locked_state[0],
						HK_BLIND_STATE_01);
	failed |= check_i32("fixed p4 bpair state before", set->bpairs[0].locked_state,
						HK_BLIND_STATE_01);
	failed |= check_close("fixed p4 p01 before", set->bpairs[0].p4[HK_BLIND_STATE_01], 1.0f);
	failed |= check_i32("fixed p4 apply chr flips",
						hk_blind_bpair_set_apply_chr_flips(set, b, chr_flipped, 1), 0);
	failed |= check_i32("fixed p4 raw state after", set->raw_locked_state[0],
						HK_BLIND_STATE_10);
	failed |= check_i32("fixed p4 bpair state after", set->bpairs[0].locked_state,
						HK_BLIND_STATE_10);
	failed |= check_close("fixed p4 p10 after", set->bpairs[0].p4[HK_BLIND_STATE_10], 1.0f);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_fixed_p4_refreshes_real_log_norm(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[1];
	struct hk_pair pairs[1];
	struct hk_blind_bpair_set *set;
	struct hk_fdg_conf conf;
	fvec3_t coords[12];
	float log_prior[HK_BLIND_N_STATE];
	int failed = 0;

	raw[0] = make_blind_pair(0, 2000000, 0, 5000000);
	memset(pairs, 0, sizeof(pairs));
	pairs[0]._.p4[HK_BLIND_STATE_00] = 0.05f;
	pairs[0]._.p4[HK_BLIND_STATE_01] = 0.80f;
	pairs[0]._.p4[HK_BLIND_STATE_10] = 0.10f;
	pairs[0]._.p4[HK_BLIND_STATE_11] = 0.05f;
	set = hk_blind_bpair_set_build(b, 1, raw);
	failed |= check_true("fixed p4 log norm set exists", set != 0);
	if (set == 0)
		return 1;
	failed |= check_i32("fixed p4 log norm lock apply",
						hk_blind_bpair_set_apply_imputed_p4_top_locks(set, pairs, 1, 100,
																	 17, 0.75f,
																	 &set->phase_lock_diag), 0);
	set->bpairs[0].real_log_norm = 0.0f;
	hk_fdg_conf_init(&conf);
	hk_blind_init_uniform_log_prior(log_prior);
	init_far_coords(coords, 12);
	set_coord(coords[hk_diploid_bid(2, HK_DIPLOID_COPY0)], 100.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(2, HK_DIPLOID_COPY1)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(5, HK_DIPLOID_COPY0)], 1.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(5, HK_DIPLOID_COPY1)], 200.0f, 0.0f, 0.0f);
	hk_blind_bpair_set_update_posterior_from_coords_params(set, &conf, coords, 1.0f, log_prior, 1.0f);
	failed |= check_close("fixed p4 remains p01", set->bpairs[0].p4[HK_BLIND_STATE_01], 1.0f);
	failed |= check_true("fixed p4 real log norm refreshed", fabsf(set->bpairs[0].real_log_norm) > 1e-5f);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_phase_like_columns_still_ignored(const struct hk_map *m, const struct hk_bmap *b)
{
	struct hk_blind_pair raw[4];
	struct hk_blind_bpair_set *set;
	struct hk_fdg_conf conf;
	fvec3_t coords[8];
	float log_prior[HK_BLIND_N_STATE];
	float first_raw_p4[HK_BLIND_N_STATE], raw_p4[HK_BLIND_N_STATE];
	int failed = 0;
	int32_t i;

	failed |= check_i32("fixture n_pairs", m->n_pairs, 4);
	for (i = 0; i < m->n_pairs && i < 4; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);

	set = hk_blind_bpair_set_build(b, 4, raw);
	failed |= check_i32("fixture set n_raw", set->n_raw, 4);
	failed |= check_i32("fixture set n_bpairs", set->n_bpairs, 1);
	failed |= check_key("fixture key", &set->bpairs[0].key, 3, 3);
	failed |= check_i32("fixture n_raw count", set->bpairs[0].n_raw, 4);

	hk_fdg_conf_init(&conf);
	hk_blind_init_uniform_log_prior(log_prior);
	init_far_coords(coords, 8);
	set_coord(coords[hk_diploid_bid(3, HK_DIPLOID_COPY0)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(3, HK_DIPLOID_COPY1)], 1.0f, 0.0f, 0.0f);
	hk_blind_bpair_set_update_posterior_from_coords(set, &conf, coords, 1.0f, 1.0f, 2.0f, log_prior, 1.0f);

	failed |= check_finite_bpair(&set->bpairs[0]);
	failed |= check_sums("fixture prob sum", &set->bpairs[0]);
	for (i = 0; i < 4; ++i) {
		failed |= check_i32("fixture raw bpair", set->raw2binned[i].bpair_id, 0);
		failed |= check_u8("fixture raw swapped", set->raw2binned[i].swapped, 0);
		hk_blind_p4_to_raw_order(set->bpairs[set->raw2binned[i].bpair_id].p4, set->raw2binned[i].swapped, raw_p4);
		if (i == 0) {
			int j;
			for (j = 0; j < HK_BLIND_N_STATE; ++j)
				first_raw_p4[j] = raw_p4[j];
		} else {
			failed |= check_p4("fixture same raw posterior", raw_p4, first_raw_p4);
		}
	}

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_raw_split_filtered_all_keeps_uniform_weak_contacts(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[2];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list wedges;
	int failed = 0, i;

	raw[0] = make_blind_pair(0, 2000000, 0, 5000000);
	raw[1] = make_blind_pair(0, 3000000, 0, 6000000);
	set = hk_blind_bpair_set_build(b, 2, raw);
	failed |= check_true("filtered-all set exists", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		int s;
		for (s = 0; s < HK_BLIND_N_STATE; ++s)
			set->bpairs[i].p4[s] = 1.0f / (float)HK_BLIND_N_STATE;
	}

	hk_blind_wedge_list_init(&wedges);
	failed |= check_i32("filtered-all build",
						hk_blind_wedge_list_build_mstep_graph(&wedges, b, set,
															  1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
															  HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
															  HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f,
															  1.0f, 1.0f,
															  HK_BLIND_STATE_WEIGHT_POSTERIOR,
															  HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL,
															  0.0f, 0.0f, 1.0f,
															  HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
															  0.0f, 1.0f, 1.0f, 0.0f), 0);
	failed |= check_i32("filtered-all candidate split pairs",
						(int32_t)wedges.n_split_candidate_pairs, 8);
	failed |= check_i32("filtered-all selected weak states",
						(int32_t)wedges.n_split_selected_raw, 8);
	failed |= check_i32("filtered-all state00 weak count",
						(int32_t)wedges.n_split_state_raw_count[HK_BLIND_STATE_00], 2);
	failed |= check_i32("filtered-all state01 weak count",
						(int32_t)wedges.n_split_state_raw_count[HK_BLIND_STATE_01], 2);
	failed |= check_i32("filtered-all state10 weak count",
						(int32_t)wedges.n_split_state_raw_count[HK_BLIND_STATE_10], 2);
	failed |= check_i32("filtered-all state11 weak count",
						(int32_t)wedges.n_split_state_raw_count[HK_BLIND_STATE_11], 2);
	failed |= check_i32("filtered-all same-bin skip",
						(int32_t)wedges.n_split_same_bin_skip_raw, 0);
	hk_blind_wedge_list_destroy(&wedges);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

int main(void)
{
	struct hk_map *m;
	struct hk_bmap *b;
	int failed = 0;

	hk_verbose = 0;
	m = hk_map_read("testdata/p9016_blind_fixture.pairs");
	if (m == 0) {
		fprintf(stderr, "failed to read fixture\n");
		return 1;
	}
	b = hk_bmap_gen(m->d, m->n_pairs, m->pairs, 1000000, 1);
	if (b == 0) {
		fprintf(stderr, "failed to build bmap\n");
		hk_map_destroy(m);
		return 1;
	}

	failed |= check_bpair_set_update(b);
	failed |= check_posterior_independent_of_bpair_params(b);
	failed |= check_estep_logdist2_breaks_fdg_flat_ties();
	failed |= check_raw_order_helper();
	failed |= check_raw2binned_inheritance(b);
	failed |= check_raw_phase_locks_record_raw_state(b);
	failed |= check_imputed_p4_top_locks_record_raw_state(b);
	failed |= check_same_bin_posterior_is_unknown(b);
	failed |= check_phase_lock_skips_posterior_update(b);
	failed |= check_raw_locks_follow_chr_flips(b);
	failed |= check_fixed_p4_locks_follow_chr_flips(b);
	failed |= check_fixed_p4_refreshes_real_log_norm(b);
	failed |= check_phase_like_columns_still_ignored(m, b);
	failed |= check_raw_split_filtered_all_keeps_uniform_weak_contacts(b);

	hk_bmap_destroy(b);
	hk_map_destroy(m);
	return failed != 0;
}
