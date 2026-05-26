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

static int check_close_rel(const char *label, float got, float expected, float rel_tol)
{
	float scale = fabsf(expected) > 1.0f? fabsf(expected) : 1.0f;
	if (fabsf(got - expected) > rel_tol * scale) {
		fprintf(stderr, "%s: got %.8g, expected %.8g within %.4g relative\n",
				label, got, expected, rel_tol);
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

static int check_true(const char *label, int ok)
{
	if (!ok) {
		fprintf(stderr, "%s: false\n", label);
		return 1;
	}
	return 0;
}

static void set_bpair(struct hk_blind_bpair *bp, float p00, float p01, float p10, float p11)
{
	memset(bp, 0, sizeof(*bp));
	bp->key.bid[0] = 2;
	bp->key.bid[1] = 5;
	bp->n_raw = 1;
	bp->base_d_scale = 1.0f;
	bp->base_k = 1.0f;
	bp->contact_class = HK_BLIND_CONTACT_CIS;
	bp->p4[HK_BLIND_STATE_00] = p00;
	bp->p4[HK_BLIND_STATE_01] = p01;
	bp->p4[HK_BLIND_STATE_10] = p10;
	bp->p4[HK_BLIND_STATE_11] = p11;
		bp->real_log_norm = 0.0f;
	bp->entropy = 0.0f;
	bp->pmax = 0.0f;
	bp->margin = 0.0f;
	bp->rho_output = 0.0f;
	bp->pU = 0.0f;
	bp->real_log_norm = 0.0f;
}

static int check_edge(const char *label, const struct hk_blind_wedge *edge, int32_t bid0, int32_t bid1,
					  float k, float d_scale, int32_t state)
{
	int failed = 0;
	char buf[64];
	snprintf(buf, sizeof(buf), "%s bid0", label);
	failed |= check_i32(buf, edge->bid[0], bid0);
	snprintf(buf, sizeof(buf), "%s bid1", label);
	failed |= check_i32(buf, edge->bid[1], bid1);
	snprintf(buf, sizeof(buf), "%s k", label);
	failed |= check_close(buf, edge->k, k);
	snprintf(buf, sizeof(buf), "%s d_scale", label);
	failed |= check_close(buf, edge->d_scale, d_scale);
	snprintf(buf, sizeof(buf), "%s state", label);
	failed |= check_i32(buf, edge->state, state);
	snprintf(buf, sizeof(buf), "%s state_mask", label);
	failed |= check_i32(buf, edge->state_mask, HK_BLIND_STATE_MASK(state));
	return failed;
}

static float sum_edge_k(const struct hk_blind_wedge edges[HK_BLIND_N_STATE])
{
	float sum = 0.0f;
	int i;
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		sum += edges[i].k;
	return sum;
}

static float sum_list_edge_k(const struct hk_blind_wedge_list *list)
{
	float sum = 0.0f;
	int32_t i;
	for (i = 0; i < list->n_edges; ++i)
		sum += list->edges[i].k;
	return sum;
}

static float sum_list_edge_d_scale(const struct hk_blind_wedge_list *list)
{
	float sum = 0.0f;
	int32_t i;
	for (i = 0; i < list->n_edges; ++i)
		sum += list->edges[i].d_scale;
	return sum;
}

static int check_uniform_expansion(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge edges[HK_BLIND_N_STATE];
	int failed = 0;

	set_bpair(&bp, 0.25f, 0.25f, 0.25f, 0.25f);
	hk_blind_bpair_expand_weighted_edges(&bp, 2.0f, 0.7f, 1.0f, edges);

	failed |= check_edge("uniform 00", &edges[HK_BLIND_STATE_00], 4, 10, 0.5f, 0.7f, HK_BLIND_STATE_00);
	failed |= check_edge("uniform 01", &edges[HK_BLIND_STATE_01], 4, 11, 0.5f, 0.7f, HK_BLIND_STATE_01);
	failed |= check_edge("uniform 10", &edges[HK_BLIND_STATE_10], 5, 10, 0.5f, 0.7f, HK_BLIND_STATE_10);
	failed |= check_edge("uniform 11", &edges[HK_BLIND_STATE_11], 5, 11, 0.5f, 0.7f, HK_BLIND_STATE_11);
	failed |= check_close("uniform conserved k", sum_edge_k(edges), 2.0f);
	return failed;
}

static int check_nonuniform_expansion(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge edges[HK_BLIND_N_STATE];
	float base_k = 3.0f, rho_train = 0.5f, base_d_scale = 0.7f;
	int failed = 0;

	set_bpair(&bp, 0.7f, 0.1f, 0.15f, 0.05f);
	hk_blind_bpair_expand_weighted_edges(&bp, base_k, base_d_scale, rho_train, edges);

	failed |= check_edge("nonuniform 00", &edges[HK_BLIND_STATE_00], 4, 10, base_k * rho_train * 0.7f, base_d_scale, HK_BLIND_STATE_00);
	failed |= check_edge("nonuniform 01", &edges[HK_BLIND_STATE_01], 4, 11, base_k * rho_train * 0.1f, base_d_scale, HK_BLIND_STATE_01);
	failed |= check_edge("nonuniform 10", &edges[HK_BLIND_STATE_10], 5, 10, base_k * rho_train * 0.15f, base_d_scale, HK_BLIND_STATE_10);
	failed |= check_edge("nonuniform 11", &edges[HK_BLIND_STATE_11], 5, 11, base_k * rho_train * 0.05f, base_d_scale, HK_BLIND_STATE_11);
	failed |= check_close("nonuniform conserved k", sum_edge_k(edges), base_k * rho_train);
	return failed;
}

static int check_zero_rho(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge edges[HK_BLIND_N_STATE];
	int failed = 0, i;

	set_bpair(&bp, 0.7f, 0.1f, 0.15f, 0.05f);
	hk_blind_bpair_expand_weighted_edges(&bp, 3.0f, 0.7f, 0.0f, edges);
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		failed |= check_close("zero rho k", edges[i].k, 0.0f);
		failed |= check_close("zero rho d_scale", edges[i].d_scale, 0.7f);
	}
	return failed;
}

static int check_entropy_rho_training_mode(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge edges[HK_BLIND_N_STATE];
	float base_k = 4.0f, rho_train = 1.5f;
	int failed = 0, i;

	set_bpair(&bp, 0.25f, 0.25f, 0.25f, 0.25f);
	bp.rho_output = 0.0f;
	failed |= check_close("entropy uniform effective rho",
						  hk_blind_bpair_effective_rho_train(&bp, rho_train, HK_BLIND_RHO_TRAIN_ENTROPY),
						  0.0f);
	hk_blind_bpair_expand_weighted_edges_mode(&bp, base_k, 0.7f, rho_train,
											  HK_BLIND_RHO_TRAIN_ENTROPY, HK_BLIND_D_SCALE_RAW_COUNT,
											  1e-6f, edges);
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		failed |= check_close("entropy uniform k", edges[i].k, 0.0f);

	hk_blind_bpair_expand_weighted_edges_mode(&bp, base_k, 0.7f, rho_train,
											  HK_BLIND_RHO_TRAIN_CONSTANT, HK_BLIND_D_SCALE_RAW_COUNT,
											  1e-6f, edges);
	failed |= check_close("constant ignores output uncertainty", sum_edge_k(edges), base_k * rho_train);

	set_bpair(&bp, 1.0f, 0.0f, 0.0f, 0.0f);
	bp.rho_output = 1.0f;
	failed |= check_close("entropy confident effective rho",
						  hk_blind_bpair_effective_rho_train(&bp, rho_train, HK_BLIND_RHO_TRAIN_ENTROPY),
						  rho_train);
	hk_blind_bpair_expand_weighted_edges_mode(&bp, base_k, 0.7f, rho_train,
											  HK_BLIND_RHO_TRAIN_ENTROPY, HK_BLIND_D_SCALE_RAW_COUNT,
											  1e-6f, edges);
	failed |= check_close("entropy confident conserved k", sum_edge_k(edges), base_k * rho_train);
	return failed;
}

static int check_floor_and_cis_trans_rho_training_modes(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge edges[HK_BLIND_N_STATE];
	float base_k = 4.0f, rho_train = 1.5f, rho_floor = 0.5f;
	int failed = 0;

	set_bpair(&bp, 0.25f, 0.25f, 0.25f, 0.25f);
	bp.rho_output = 0.0f;
	bp.contact_class = HK_BLIND_CONTACT_CIS;
	failed |= check_close("entropy floor cis effective rho",
						  hk_blind_bpair_effective_rho_train_floor(&bp, rho_train,
																	HK_BLIND_RHO_TRAIN_ENTROPY_WITH_FLOOR,
																	rho_floor),
						  rho_train * rho_floor);
	hk_blind_bpair_expand_weighted_edges_mode_ex(&bp, base_k, 0.7f, rho_train,
												 HK_BLIND_RHO_TRAIN_ENTROPY_WITH_FLOOR,
												 rho_floor, HK_BLIND_D_SCALE_RAW_COUNT,
												 1e-6f, 1.0f, 1.0f, edges);
	failed |= check_close("entropy floor conserved k", sum_edge_k(edges), base_k * rho_train * rho_floor);

	failed |= check_close("cis entropy still zero",
						  hk_blind_bpair_effective_rho_train_floor(&bp, rho_train,
																	HK_BLIND_RHO_TRAIN_ENTROPY_CIS_CONSTANT_TRANS,
																	rho_floor),
						  0.0f);
	bp.contact_class = HK_BLIND_CONTACT_TRANS;
	failed |= check_close("trans constant effective rho",
						  hk_blind_bpair_effective_rho_train_floor(&bp, rho_train,
																	HK_BLIND_RHO_TRAIN_ENTROPY_CIS_CONSTANT_TRANS,
																	rho_floor),
						  rho_train);
	failed |= check_close("trans floor effective rho",
						  hk_blind_bpair_effective_rho_train_floor(&bp, rho_train,
																	HK_BLIND_RHO_TRAIN_ENTROPY_CIS_FLOOR_TRANS,
																	rho_floor),
						  rho_train * rho_floor);
	hk_blind_bpair_expand_weighted_edges_mode_ex(&bp, base_k, 0.7f, rho_train,
												 HK_BLIND_RHO_TRAIN_ENTROPY_CIS_FLOOR_TRANS,
												 rho_floor, HK_BLIND_D_SCALE_RAW_COUNT,
												 1e-6f, 1.0f, 2.0f, edges);
	failed |= check_close("trans multiplier conserved k", sum_edge_k(edges),
						  2.0f * base_k * rho_train * rho_floor);
	return failed;
}

static int check_expected_count_d_scale_mode(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge raw_edges[HK_BLIND_N_STATE], expected_edges[HK_BLIND_N_STATE];
	float base_d_scale = 0.5f, eps_count = 1e-3f;
	float d_high, d_low, d_zero;
	int failed = 0;

	set_bpair(&bp, 0.75f, 0.25f, 0.0f, 0.0f);
	bp.n_raw = 8;
	hk_blind_bpair_expand_weighted_edges_mode(&bp, 2.0f, base_d_scale, 1.0f,
											  HK_BLIND_RHO_TRAIN_CONSTANT, HK_BLIND_D_SCALE_RAW_COUNT,
											  eps_count, raw_edges);
	hk_blind_bpair_expand_weighted_edges_mode(&bp, 2.0f, base_d_scale, 1.0f,
											  HK_BLIND_RHO_TRAIN_CONSTANT, HK_BLIND_D_SCALE_EXPECTED_COUNT,
											  eps_count, expected_edges);

	d_high = powf(8.0f * 0.75f, -1.0f / 3.0f);
	d_low = powf(8.0f * 0.25f, -1.0f / 3.0f);
	d_zero = powf(eps_count, -1.0f / 3.0f);
	failed |= check_close("raw-count high d_scale", raw_edges[HK_BLIND_STATE_00].d_scale, base_d_scale);
	failed |= check_close("raw-count low d_scale", raw_edges[HK_BLIND_STATE_01].d_scale, base_d_scale);
	failed |= check_close("expected-count high d_scale", expected_edges[HK_BLIND_STATE_00].d_scale, d_high);
	failed |= check_close("expected-count low d_scale", expected_edges[HK_BLIND_STATE_01].d_scale, d_low);
	failed |= check_greater("expected-count low larger than high",
							expected_edges[HK_BLIND_STATE_01].d_scale,
							expected_edges[HK_BLIND_STATE_00].d_scale);
	failed |= check_close("expected-count zero posterior finite",
						  expected_edges[HK_BLIND_STATE_10].d_scale, d_zero);
	failed |= check_close("expected-count conserved k", sum_edge_k(expected_edges), 2.0f);
	return failed;
}

static int check_weighted_scalar_energy(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge_half, edge_full;
	float distance, e_half, e_full;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	edge_half.bid[0] = 4;
	edge_half.bid[1] = 10;
	edge_half.k = 0.5f;
	edge_half.d_scale = 0.7f;
	edge_half.state = HK_BLIND_STATE_00;
	edge_half.state_mask = HK_BLIND_STATE_MASK(HK_BLIND_STATE_00);
	edge_full = edge_half;
	edge_full.k = 1.0f;

	distance = 2.5f * edge_half.d_scale;
	e_half = hk_blind_wedge_contact_energy(&conf, &edge_half, distance, 1.0f);
	e_full = hk_blind_wedge_contact_energy(&conf, &edge_full, distance, 1.0f);
	failed |= check_greater("weighted energy positive", e_full, 0.0f);
	failed |= check_close("weighted energy linear k", 2.0f * e_half, e_full);
	return failed;
}

static int check_edges_match(const char *label, const struct hk_blind_wedge got[HK_BLIND_N_STATE],
							 const struct hk_blind_wedge expected[HK_BLIND_N_STATE])
{
	int failed = 0, s;
	char buf[96];
	for (s = 0; s < HK_BLIND_N_STATE; ++s) {
		snprintf(buf, sizeof(buf), "%s state%d", label, s);
		failed |= check_edge(buf, &got[s], expected[s].bid[0], expected[s].bid[1],
							 expected[s].k, expected[s].d_scale, expected[s].state);
	}
	return failed;
}

static int check_state_weight_posterior_matches_legacy_mode(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge old_edges[HK_BLIND_N_STATE], new_edges[HK_BLIND_N_STATE];
	int failed = 0;

	set_bpair(&bp, 0.55f, 0.15f, 0.25f, 0.05f);
	bp.n_raw = 7;
	bp.rho_output = 0.8f;
	hk_blind_bpair_expand_weighted_edges_mode(&bp, 3.0f, 0.7f, 1.25f,
											  HK_BLIND_RHO_TRAIN_ENTROPY_WITH_FLOOR,
											  HK_BLIND_D_SCALE_EXPECTED_COUNT,
											  1e-3f, old_edges);
	hk_blind_bpair_expand_weighted_edges_mode_state_weight(&bp, 3.0f, 0.7f, 1.25f,
														   HK_BLIND_RHO_TRAIN_ENTROPY_WITH_FLOOR,
														   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
														   HK_BLIND_D_SCALE_EXPECTED_COUNT,
														   1e-3f, 1.0f, 1.0f,
														   HK_BLIND_STATE_WEIGHT_POSTERIOR,
														   0.9f, 0.9f, 1.0f, new_edges);
	failed |= check_edges_match("posterior state-weight matches legacy", new_edges, old_edges);
	return failed;
}

static int check_state_weight_posterior_power_gamma(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge edges[HK_BLIND_N_STATE];
	float norm;
	int failed = 0;

	set_bpair(&bp, 0.50f, 0.25f, 0.25f, 0.0f);
	norm = 0.50f * 0.50f + 0.25f * 0.25f + 0.25f * 0.25f;
	hk_blind_bpair_expand_weighted_edges_mode_state_weight(&bp, 4.0f, 0.5f, 1.0f,
														   HK_BLIND_RHO_TRAIN_CONSTANT,
														   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
														   HK_BLIND_D_SCALE_RAW_COUNT,
														   1e-3f, 1.0f, 1.0f,
														   HK_BLIND_STATE_WEIGHT_POSTERIOR,
														   0.0f, 0.0f, 2.0f, edges);
	failed |= check_close("posterior gamma state00", edges[HK_BLIND_STATE_00].k,
						  4.0f * 0.50f * 0.50f / norm);
	failed |= check_close("posterior gamma state01", edges[HK_BLIND_STATE_01].k,
						  4.0f * 0.25f * 0.25f / norm);
	failed |= check_close("posterior gamma state10", edges[HK_BLIND_STATE_10].k,
						  4.0f * 0.25f * 0.25f / norm);
	failed |= check_close("posterior gamma keeps total k", sum_edge_k(edges), 4.0f);
	return failed;
}

static int check_state_weight_binary_support_uses_full_edges(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge edges[HK_BLIND_N_STATE];
	float eps_count = 1e-3f;
	float d_high, d_low, d_zero;
	int failed = 0;

	set_bpair(&bp, 0.75f, 0.25f, 0.0f, 0.0f);
	bp.n_raw = 8;
	hk_blind_bpair_expand_weighted_edges_mode_state_weight(&bp, 2.0f, 0.5f, 1.0f,
														   HK_BLIND_RHO_TRAIN_CONSTANT,
														   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
														   HK_BLIND_D_SCALE_EXPECTED_COUNT,
														   eps_count, 1.0f, 1.0f,
														   HK_BLIND_STATE_WEIGHT_BINARY_SUPPORT,
														   0.0f, 0.20f, 1.0f, edges);

	d_high = powf(8.0f * 0.75f, -1.0f / 3.0f);
	d_low = powf(8.0f * 0.25f, -1.0f / 3.0f);
	d_zero = powf(eps_count, -1.0f / 3.0f);
	failed |= check_edge("binary-support state00 full edge", &edges[HK_BLIND_STATE_00],
						 hk_diploid_bid(bp.key.bid[0], HK_DIPLOID_COPY0),
						 hk_diploid_bid(bp.key.bid[1], HK_DIPLOID_COPY0),
						 2.0f, d_high, HK_BLIND_STATE_00);
	failed |= check_edge("binary-support state01 full edge", &edges[HK_BLIND_STATE_01],
						 hk_diploid_bid(bp.key.bid[0], HK_DIPLOID_COPY0),
						 hk_diploid_bid(bp.key.bid[1], HK_DIPLOID_COPY1),
						 2.0f, d_low, HK_BLIND_STATE_01);
	failed |= check_edge("binary-support state10 below threshold", &edges[HK_BLIND_STATE_10],
						 hk_diploid_bid(bp.key.bid[0], HK_DIPLOID_COPY1),
						 hk_diploid_bid(bp.key.bid[1], HK_DIPLOID_COPY0),
						 0.0f, d_zero, HK_BLIND_STATE_10);
	failed |= check_edge("binary-support state11 below threshold", &edges[HK_BLIND_STATE_11],
						 hk_diploid_bid(bp.key.bid[0], HK_DIPLOID_COPY1),
						 hk_diploid_bid(bp.key.bid[1], HK_DIPLOID_COPY1),
						 0.0f, d_zero, HK_BLIND_STATE_11);
	failed |= check_close("binary-support does not conserve posterior k", sum_edge_k(edges), 4.0f);
	return failed;
}

static int check_state_weight_top_only_thresholds(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge edges[HK_BLIND_N_STATE];
	int failed = 0;

	set_bpair(&bp, 0.60f, 0.30f, 0.10f, 0.0f);
	bp.n_raw = 10;
	hk_blind_bpair_expand_weighted_edges_mode_state_weight(&bp, 2.0f, 0.5f, 1.0f,
														   HK_BLIND_RHO_TRAIN_CONSTANT,
														   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
														   HK_BLIND_D_SCALE_RAW_COUNT,
														   1e-3f, 1.0f, 1.0f,
														   HK_BLIND_STATE_WEIGHT_TOP_ONLY,
														   0.20f, 0.50f, 1.0f, edges);
	failed |= check_close("top-only keeps top state", edges[HK_BLIND_STATE_00].k, 2.0f);
	failed |= check_close("top-only drops second state", edges[HK_BLIND_STATE_01].k, 0.0f);
	failed |= check_close("top-only drops third state", edges[HK_BLIND_STATE_10].k, 0.0f);
	failed |= check_close("top-only drops fourth state", edges[HK_BLIND_STATE_11].k, 0.0f);

	hk_blind_bpair_expand_weighted_edges_mode_state_weight(&bp, 2.0f, 0.5f, 1.0f,
														   HK_BLIND_RHO_TRAIN_CONSTANT,
														   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
														   HK_BLIND_D_SCALE_RAW_COUNT,
														   1e-3f, 1.0f, 1.0f,
														   HK_BLIND_STATE_WEIGHT_TOP_ONLY,
														   0.40f, 0.50f, 1.0f, edges);
	failed |= check_close("top-only margin gate zeros all states", sum_edge_k(edges), 0.0f);
	return failed;
}

static int check_fixed_p4_lock_is_not_updated(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf conf;
	fvec3_t coords[12];
	float log_prior[HK_BLIND_N_STATE] = {0.0f, 0.0f, 0.0f, 0.0f};
	float old_p4[HK_BLIND_N_STATE];
	int i, a, failed = 0;

	memset(&set, 0, sizeof(set));
	set_bpair(&bp, 0.20f, 0.30f, 0.40f, 0.10f);
	bp.lock_mode = HK_BLIND_LOCK_FIXED_P4;
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		old_p4[i] = bp.p4[i];
	set.bpairs = &bp;
	set.n_bpairs = 1;

	for (i = 0; i < 12; ++i)
		for (a = 0; a < 3; ++a)
			coords[i][a] = (float)(i + a + 1);
	hk_fdg_conf_init(&conf);
	hk_blind_bpair_set_update_posterior_from_coords_params(&set, &conf, coords, 1.0f,
														   log_prior, 1.0f);
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		failed |= check_close("fixed-p4 lock unchanged", bp.p4[i], old_p4[i]);
	return failed;
}

static int check_zero_weight_wedges_are_not_list_edges(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	int failed = 0, ret;

	memset(&set, 0, sizeof(set));
	set_bpair(&bp, 1.0f, 0.0f, 0.0f, 0.0f);
	bp.n_raw = 4;
	bp.base_k = 1.0f;
	bp.base_d_scale = 0.5f;
	set.bpairs = &bp;
	set.n_bpairs = 1;

	hk_blind_wedge_list_init(&list);
	ret = hk_blind_wedge_list_build_from_bpair_set_params_state_weight(&list, &set,
																	   1.0f,
																	   HK_BLIND_RHO_TRAIN_CONSTANT,
																	   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
																	   HK_BLIND_D_SCALE_RAW_COUNT,
																	   1e-6f, 1.0f, 1.0f,
																	   HK_BLIND_STATE_WEIGHT_POSTERIOR,
																	   0.0f, 0.0f, 1.0f);
	failed |= check_i32("zero-weight build ret", ret, 0);
	failed |= check_i32("zero-weight list edges", list.n_edges, 1);
	if (list.n_edges == 1)
		failed |= check_edge("zero-weight retained edge", &list.edges[0],
							 4, 10, 1.0f, 0.5f, HK_BLIND_STATE_00);
	ret = hk_blind_wedge_list_aggregate_exact(&list);
	failed |= check_i32("zero-weight aggregate ret", ret, 0);
	failed |= check_i32("zero-weight aggregate edges", list.n_edges, 1);
	failed |= check_i32("zero-weight aggregate removed", list.n_aggregated_edges_removed, 0);
	if (list.n_edges == 1)
		failed |= check_i32("zero-weight aggregate state mask", list.edges[0].state_mask,
							HK_BLIND_STATE_MASK(HK_BLIND_STATE_00));
	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_aggregate_drops_zero_weight_edges(void)
{
	struct hk_blind_wedge edges[2];
	struct hk_blind_wedge_list list;
	int failed = 0, ret;

	hk_blind_wedge_list_init(&list);
	list.edges = edges;
	list.n_edges = 2;
	list.m_edges = 2;
	edges[0].bid[0] = 4;
	edges[0].bid[1] = 10;
	edges[0].k = 0.0f;
	edges[0].d_scale = 0.5f;
	edges[0].state = HK_BLIND_STATE_00;
	edges[0].state_mask = HK_BLIND_STATE_MASK(HK_BLIND_STATE_00);
	edges[1].bid[0] = 4;
	edges[1].bid[1] = 10;
	edges[1].k = 2.0f;
	edges[1].d_scale = 0.5f;
	edges[1].state = HK_BLIND_STATE_01;
	edges[1].state_mask = HK_BLIND_STATE_MASK(HK_BLIND_STATE_01);

	ret = hk_blind_wedge_list_aggregate_exact(&list);
	failed |= check_i32("zero aggregate ret", ret, 0);
	failed |= check_i32("zero aggregate final edges", list.n_edges, 1);
	failed |= check_i32("zero aggregate removed", list.n_aggregated_edges_removed, 1);
	if (list.n_edges == 1) {
		failed |= check_close("zero aggregate final k", list.edges[0].k, 2.0f);
		failed |= check_i32("zero aggregate final state mask", list.edges[0].state_mask,
							HK_BLIND_STATE_MASK(HK_BLIND_STATE_01));
	}
	// list.edges points to stack storage in this test.
	list.edges = 0;
	list.n_edges = 0;
	list.m_edges = 0;
	return failed;
}

static int check_aggregate_drops_all_zero_weight_edges(void)
{
	struct hk_blind_wedge edges[2];
	struct hk_blind_wedge_list list;
	int failed = 0, ret;

	hk_blind_wedge_list_init(&list);
	list.edges = edges;
	list.n_edges = 2;
	list.m_edges = 2;
	edges[0].bid[0] = 4;
	edges[0].bid[1] = 10;
	edges[0].k = 0.0f;
	edges[0].d_scale = 0.5f;
	edges[0].state = HK_BLIND_STATE_00;
	edges[0].state_mask = HK_BLIND_STATE_MASK(HK_BLIND_STATE_00);
	edges[1].bid[0] = 5;
	edges[1].bid[1] = 11;
	edges[1].k = 0.0f;
	edges[1].d_scale = 0.5f;
	edges[1].state = HK_BLIND_STATE_11;
	edges[1].state_mask = HK_BLIND_STATE_MASK(HK_BLIND_STATE_11);

	ret = hk_blind_wedge_list_aggregate_exact(&list);
	failed |= check_i32("all-zero aggregate ret", ret, 0);
	failed |= check_i32("all-zero aggregate final edges", list.n_edges, 0);
	failed |= check_i32("all-zero aggregate removed", list.n_aggregated_edges_removed, 2);
	// list.edges points to stack storage in this test.
	list.edges = 0;
	list.n_edges = 0;
	list.m_edges = 0;
	return failed;
}

static void init_raw_split_bmap(struct hk_bmap *bmap, struct hk_sdict *dict,
								char **names, int32_t *len,
								struct hk_bead *beads, uint64_t *offcnt)
{
	int32_t i;
	memset(dict, 0, sizeof(*dict));
	dict->n = 1;
	dict->m = 1;
	dict->name = names;
	dict->len = len;
	dict->h = 0;
	memset(bmap, 0, sizeof(*bmap));
	bmap->d = dict;
	bmap->n_beads = 6;
	bmap->unit = 1.0f;
	bmap->beads = beads;
	bmap->offcnt = offcnt;
	for (i = 0; i < bmap->n_beads; ++i) {
		beads[i].chr = 0;
		beads[i].st = i * 1000000;
		beads[i].en = (i + 1) * 1000000;
	}
	offcnt[0] = (uint64_t)0 << 32 | 6u;
	offcnt[1] = (uint64_t)6 << 32 | 0u;
}

static struct hk_blind_pair make_raw_split_pair(int32_t pos0, int32_t pos1)
{
	struct hk_blind_pair p;
	p.chr[0] = 0;
	p.chr[1] = 0;
	p.pos[0] = pos0;
	p.pos[1] = pos1;
	p.strand[0] = 1;
	p.strand[1] = -1;
	return p;
}

static int check_mstep_raw_split_top_uses_raw_contacts(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[8];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list list;
	int failed = 0, ret, i;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	for (i = 0; i < 4; ++i) {
		raw[i] = make_raw_split_pair(i * 1000, 5500000 - i * 1000);
		raw[i + 4] = make_raw_split_pair(5500000 - i * 1000, i * 1000);
	}
	set = hk_blind_bpair_set_build(&bmap, 8, raw);
	failed |= check_true("raw split set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 1.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
	}
	hk_blind_wedge_list_init(&list);
	ret = hk_blind_wedge_list_build_mstep_graph(&list, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_TOP_ONLY,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_TOP,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("raw split build ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&list);
	failed |= check_i32("raw split aggregate ret", ret, 0);
	failed |= check_i32("raw split edge count", list.n_edges, 1);
	if (list.n_edges == 1) {
		failed |= check_i32("raw split bid0", list.edges[0].bid[0], hk_diploid_bid(0, 0));
		failed |= check_i32("raw split bid1", list.edges[0].bid[1], hk_diploid_bid(5, 1));
		failed |= check_close("raw split d_scale", list.edges[0].d_scale, powf(8.0f, -1.0f / 3.0f));
		failed |= check_close("raw split k", list.edges[0].k, 1.0f);
		failed |= check_i32("raw split state mask", list.edges[0].state_mask, 0xff);
	}
	hk_blind_wedge_list_destroy(&list);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int find_edge_by_bids(const struct hk_blind_wedge_list *list, int32_t bid0, int32_t bid1)
{
	int32_t i;
	if (bid1 < bid0) {
		int32_t t = bid0;
		bid0 = bid1;
		bid1 = t;
	}
	for (i = 0; i < list->n_edges; ++i)
		if (list->edges[i].bid[0] == bid0 && list->edges[i].bid[1] == bid1)
			return i;
	return -1;
}

static int check_list_same_edges_for_k_and_d_scale(const char *label,
												   const struct hk_blind_wedge_list *got,
												   const struct hk_blind_wedge_list *k_ref,
												   const struct hk_blind_wedge_list *d_ref)
{
	int failed = 0;
	int32_t i;
	char buf[128];
	assert(label);
	assert(got);
	assert(k_ref);
	assert(d_ref);
	failed |= check_i32("edge compare n_edges got/k_ref", got->n_edges, k_ref->n_edges);
	failed |= check_i32("edge compare n_edges got/d_ref", got->n_edges, d_ref->n_edges);
	for (i = 0; i < got->n_edges; ++i) {
		const struct hk_blind_wedge *e = &got->edges[i];
		int ik = find_edge_by_bids(k_ref, e->bid[0], e->bid[1]);
		int id = find_edge_by_bids(d_ref, e->bid[0], e->bid[1]);
		snprintf(buf, sizeof(buf), "%s k ref edge %d", label, (int)i);
		failed |= check_true(buf, ik >= 0);
		snprintf(buf, sizeof(buf), "%s d ref edge %d", label, (int)i);
		failed |= check_true(buf, id >= 0);
		if (ik >= 0) {
			snprintf(buf, sizeof(buf), "%s edge %d k", label, (int)i);
			failed |= check_close(buf, e->k, k_ref->edges[ik].k);
			snprintf(buf, sizeof(buf), "%s edge %d state mask", label, (int)i);
			failed |= check_i32(buf, e->state_mask, k_ref->edges[ik].state_mask);
		}
		if (id >= 0) {
			snprintf(buf, sizeof(buf), "%s edge %d d_scale", label, (int)i);
			failed |= check_close(buf, e->d_scale, d_ref->edges[id].d_scale);
		}
	}
	return failed;
}

static int check_list_same_edges_for_d_scale(const char *label,
											 const struct hk_blind_wedge_list *got,
											 const struct hk_blind_wedge_list *ref)
{
	int failed = 0;
	int32_t i;
	char buf[128];
	assert(label);
	assert(got);
	assert(ref);
	failed |= check_i32("edge dscale compare n_edges", got->n_edges, ref->n_edges);
	for (i = 0; i < got->n_edges; ++i) {
		const struct hk_blind_wedge *e = &got->edges[i];
		int id = find_edge_by_bids(ref, e->bid[0], e->bid[1]);
		snprintf(buf, sizeof(buf), "%s ref edge %d", label, (int)i);
		failed |= check_true(buf, id >= 0);
		if (id >= 0) {
			snprintf(buf, sizeof(buf), "%s edge %d d_scale", label, (int)i);
			failed |= check_close(buf, e->d_scale, ref->edges[id].d_scale);
		}
	}
	return failed;
}

static int check_raw_split_mstep_posterior_power_preserves_mass(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[8];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list base, gamma2, full;
	int failed = 0, ret, i;
	int e_base_01, e_base_10, e_gamma_01, e_gamma_10;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	for (i = 0; i < 4; ++i) {
		raw[i] = make_raw_split_pair(i * 1000, 5500000 - i * 1000);
		raw[i + 4] = make_raw_split_pair(5500000 - i * 1000, i * 1000);
	}
	set = hk_blind_bpair_set_build(&bmap, 8, raw);
	failed |= check_true("raw split power set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.60f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.40f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
		set->bpairs[i].entropy = 0.67301166f;
		set->bpairs[i].pmax = 0.60f;
	}

	hk_blind_wedge_list_init(&base);
	ret = hk_blind_wedge_list_build_mstep_graph_ex(
		&base, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS, 0.0f, 1.0f);
	failed |= check_i32("raw split power base ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&base);
	failed |= check_i32("raw split power base aggregate ret", ret, 0);

	hk_blind_wedge_list_init(&gamma2);
	ret = hk_blind_wedge_list_build_mstep_graph_ex(
		&gamma2, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		1.0f, 1.0f, 1.0f, 0.0f, 2.0f, 0, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS, 0.0f, 1.0f);
	failed |= check_i32("raw split power gamma2 ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&gamma2);
	failed |= check_i32("raw split power gamma2 aggregate ret", ret, 0);

	hk_blind_wedge_list_init(&full);
	ret = hk_blind_wedge_list_build_mstep_graph_ex(
		&full, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_FULL,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS, 0.0f, 1.0f);
	failed |= check_i32("raw split power full ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&full);
	failed |= check_i32("raw split power full aggregate ret", ret, 0);

	ret = hk_blind_wedge_list_build_mstep_graph_ex(
		&full, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_FULL,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		1.0f, 1.0f, 1.0f, 0.0f, 2.0f, 0, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS, 0.0f, 1.0f);
	failed |= check_i32("raw split power rejects full mode ret", ret, -1);

	failed |= check_i32("raw split power edge count", gamma2.n_edges, base.n_edges);
	failed |= check_i32("raw split power candidate count",
						(int32_t)gamma2.n_split_candidate_pairs,
						(int32_t)base.n_split_candidate_pairs);
	failed |= check_i32("raw split power filter1 count",
						(int32_t)gamma2.n_split_filter1_pairs,
						(int32_t)base.n_split_filter1_pairs);
	failed |= check_i32("raw split power filter2 count",
						(int32_t)gamma2.n_split_filter2_pairs,
						(int32_t)base.n_split_filter2_pairs);
	failed |= check_i32("raw split power bmap count",
						(int32_t)gamma2.n_split_bmap_pairs,
						(int32_t)base.n_split_bmap_pairs);
	failed |= check_close_rel("raw split power final mass nearly preserved",
							  sum_list_edge_k(&gamma2),
							  sum_list_edge_k(&base), 0.01f);
	failed |= check_close("raw split power dscale unchanged",
						  sum_list_edge_d_scale(&gamma2),
						  sum_list_edge_d_scale(&full));
	failed |= check_list_same_edges_for_d_scale("raw split power dscale per-edge",
												&gamma2, &base);
	e_base_01 = find_edge_by_bids(&base, hk_diploid_bid(0, 0), hk_diploid_bid(5, 1));
	e_base_10 = find_edge_by_bids(&base, hk_diploid_bid(0, 1), hk_diploid_bid(5, 0));
	e_gamma_01 = find_edge_by_bids(&gamma2, hk_diploid_bid(0, 0), hk_diploid_bid(5, 1));
	e_gamma_10 = find_edge_by_bids(&gamma2, hk_diploid_bid(0, 1), hk_diploid_bid(5, 0));
	failed |= check_true("raw split power base 01 edge", e_base_01 >= 0);
	failed |= check_true("raw split power base 10 edge", e_base_10 >= 0);
	failed |= check_true("raw split power gamma 01 edge", e_gamma_01 >= 0);
	failed |= check_true("raw split power gamma 10 edge", e_gamma_10 >= 0);
	if (e_base_01 >= 0 && e_base_10 >= 0 && e_gamma_01 >= 0 && e_gamma_10 >= 0) {
		failed |= check_greater("raw split power top state grows",
								gamma2.edges[e_gamma_01].k,
								base.edges[e_base_01].k);
		failed |= check_greater("raw split power weak state shrinks",
								base.edges[e_base_10].k,
								gamma2.edges[e_gamma_10].k);
	}

	hk_blind_wedge_list_destroy(&full);
	hk_blind_wedge_list_destroy(&gamma2);
	hk_blind_wedge_list_destroy(&base);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_raw_split_kweight_final_prob_survives_sort(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[16];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list base, gamma2;
	int failed = 0, ret, i;
	int a01, a10, b01, b10;
	int ga01, ga10, gb01, gb10;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	for (i = 0; i < 8; ++i) {
		raw[i] = make_raw_split_pair(1000000 + i * 1000, 4500000 - i * 1000);
		raw[i + 8] = make_raw_split_pair(i * 1000, 5500000 - i * 1000);
	}
	set = hk_blind_bpair_set_build(&bmap, 16, raw);
	failed |= check_true("raw split kweight sort set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		struct hk_blind_bpair *bp = &set->bpairs[i];
		if (bp->key.bid[0] == 0 && bp->key.bid[1] == 5) {
			bp->p4[HK_BLIND_STATE_00] = 0.0f;
			bp->p4[HK_BLIND_STATE_01] = 0.80f;
			bp->p4[HK_BLIND_STATE_10] = 0.20f;
			bp->p4[HK_BLIND_STATE_11] = 0.0f;
		bp->real_log_norm = 0.0f;
			bp->pmax = 0.80f;
		} else if (bp->key.bid[0] == 1 && bp->key.bid[1] == 4) {
			bp->p4[HK_BLIND_STATE_00] = 0.0f;
			bp->p4[HK_BLIND_STATE_01] = 0.20f;
			bp->p4[HK_BLIND_STATE_10] = 0.80f;
			bp->p4[HK_BLIND_STATE_11] = 0.0f;
		bp->real_log_norm = 0.0f;
			bp->pmax = 0.80f;
		} else {
			failed |= check_true("raw split kweight unexpected bpair", 0);
		}
		bp->entropy = 0.50040245f;
	}

	hk_blind_wedge_list_init(&base);
	ret = hk_blind_wedge_list_build_mstep_graph_ex(
		&base, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS, 0.0f, 1.0f);
	failed |= check_i32("raw split kweight sort base ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&base);
	failed |= check_i32("raw split kweight sort base aggregate ret", ret, 0);

	hk_blind_wedge_list_init(&gamma2);
	ret = hk_blind_wedge_list_build_mstep_graph_ex(
		&gamma2, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		1.0f, 1.0f, 1.0f, 0.0f, 2.0f, 0, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS, HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS, 0.0f, 1.0f);
	failed |= check_i32("raw split kweight sort gamma ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&gamma2);
	failed |= check_i32("raw split kweight sort gamma aggregate ret", ret, 0);

	failed |= check_i32("raw split kweight sort edge count", gamma2.n_edges, base.n_edges);
	a01 = find_edge_by_bids(&base, hk_diploid_bid(0, 0), hk_diploid_bid(5, 1));
	a10 = find_edge_by_bids(&base, hk_diploid_bid(0, 1), hk_diploid_bid(5, 0));
	b01 = find_edge_by_bids(&base, hk_diploid_bid(1, 0), hk_diploid_bid(4, 1));
	b10 = find_edge_by_bids(&base, hk_diploid_bid(1, 1), hk_diploid_bid(4, 0));
	ga01 = find_edge_by_bids(&gamma2, hk_diploid_bid(0, 0), hk_diploid_bid(5, 1));
	ga10 = find_edge_by_bids(&gamma2, hk_diploid_bid(0, 1), hk_diploid_bid(5, 0));
	gb01 = find_edge_by_bids(&gamma2, hk_diploid_bid(1, 0), hk_diploid_bid(4, 1));
	gb10 = find_edge_by_bids(&gamma2, hk_diploid_bid(1, 1), hk_diploid_bid(4, 0));
	failed |= check_true("raw split kweight sort A 01", a01 >= 0 && ga01 >= 0);
	failed |= check_true("raw split kweight sort A 10", a10 >= 0 && ga10 >= 0);
	failed |= check_true("raw split kweight sort B 01", b01 >= 0 && gb01 >= 0);
	failed |= check_true("raw split kweight sort B 10", b10 >= 0 && gb10 >= 0);
	if (a01 >= 0 && a10 >= 0 && b01 >= 0 && b10 >= 0 &&
		ga01 >= 0 && ga10 >= 0 && gb01 >= 0 && gb10 >= 0) {
		failed |= check_greater("raw split kweight sort A top grows",
								gamma2.edges[ga01].k, base.edges[a01].k);
		failed |= check_greater("raw split kweight sort A weak shrinks",
								base.edges[a10].k, gamma2.edges[ga10].k);
		failed |= check_greater("raw split kweight sort B top grows",
								gamma2.edges[gb10].k, base.edges[b10].k);
		failed |= check_greater("raw split kweight sort B weak shrinks",
								base.edges[b01].k, gamma2.edges[gb01].k);
	}

	hk_blind_wedge_list_destroy(&gamma2);
	hk_blind_wedge_list_destroy(&base);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_mstep_raw_split_soft_confidence(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[1];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list list;
	int failed = 0, ret, e01, e10, e11;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	raw[0] = make_raw_split_pair(0, 5500000);
	set = hk_blind_bpair_set_build(&bmap, 1, raw);
	failed |= check_true("raw soft set", set != 0);
	if (set == 0)
		return 1;
	set->bpairs[0].p4[HK_BLIND_STATE_00] = 0.25f;
	set->bpairs[0].p4[HK_BLIND_STATE_01] = 0.25f;
	set->bpairs[0].p4[HK_BLIND_STATE_10] = 0.25f;
	set->bpairs[0].p4[HK_BLIND_STATE_11] = 0.25f;
	set->bpairs[0].entropy = logf(4.0f);
	set->bpairs[0].pmax = 0.25f;

	hk_blind_wedge_list_init(&list);
	ret = hk_blind_wedge_list_build_mstep_graph(&list, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
												0.1f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("raw soft uniform build ret", ret, 0);
	failed |= check_i32("raw soft uniform edges before aggregate", list.n_edges, 4);
	ret = hk_blind_wedge_list_aggregate_exact(&list);
	failed |= check_i32("raw soft uniform aggregate ret", ret, 0);
	failed |= check_i32("raw soft uniform edge count", list.n_edges, 4);
	e01 = find_edge_by_bids(&list, hk_diploid_bid(0, 0), hk_diploid_bid(5, 1));
	failed |= check_true("raw soft uniform 01 edge", e01 >= 0);
	if (e01 >= 0)
		failed |= check_close("raw soft uniform weak k", list.edges[e01].k, 0.025f);
	hk_blind_wedge_list_destroy(&list);

	set->bpairs[0].p4[HK_BLIND_STATE_00] = 0.02f;
	set->bpairs[0].p4[HK_BLIND_STATE_01] = 0.90f;
	set->bpairs[0].p4[HK_BLIND_STATE_10] = 0.04f;
	set->bpairs[0].p4[HK_BLIND_STATE_11] = 0.04f;
	set->bpairs[0].entropy = 0.42870036f;
	set->bpairs[0].pmax = 0.90f;

	hk_blind_wedge_list_init(&list);
	ret = hk_blind_wedge_list_build_mstep_graph(&list, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_PMAX_MARGIN,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("raw soft high build ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&list);
	failed |= check_i32("raw soft high aggregate ret", ret, 0);
	e01 = find_edge_by_bids(&list, hk_diploid_bid(0, 0), hk_diploid_bid(5, 1));
	e11 = find_edge_by_bids(&list, hk_diploid_bid(0, 1), hk_diploid_bid(5, 1));
	failed |= check_true("raw soft high 01 edge", e01 >= 0);
	failed |= check_true("raw soft high 11 edge", e11 >= 0);
	if (e01 >= 0 && e11 >= 0)
		failed |= check_greater("raw soft high top dominates", list.edges[e01].k, list.edges[e11].k);
	hk_blind_wedge_list_destroy(&list);

	hk_blind_wedge_list_init(&list);
	ret = hk_blind_wedge_list_build_mstep_graph(&list, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_TOP1,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("raw soft top1 build ret", ret, 0);
	failed |= check_i32("raw soft top1 edges before aggregate", list.n_edges, 1);
	ret = hk_blind_wedge_list_aggregate_exact(&list);
	failed |= check_i32("raw soft top1 aggregate ret", ret, 0);
	failed |= check_i32("raw soft top1 edge count", list.n_edges, 1);
	e01 = find_edge_by_bids(&list, hk_diploid_bid(0, 0), hk_diploid_bid(5, 1));
	failed |= check_true("raw soft top1 kept top edge", e01 >= 0);
	if (e01 >= 0)
		failed |= check_close("raw soft top1 top k", list.edges[e01].k, 0.62168223f);
	hk_blind_wedge_list_destroy(&list);

	hk_blind_wedge_list_init(&list);
	ret = hk_blind_wedge_list_build_mstep_graph(&list, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_TOP2,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("raw soft top2 build ret", ret, 0);
	failed |= check_i32("raw soft top2 edges before aggregate", list.n_edges, 2);
	ret = hk_blind_wedge_list_aggregate_exact(&list);
	failed |= check_i32("raw soft top2 aggregate ret", ret, 0);
	failed |= check_i32("raw soft top2 edge count", list.n_edges, 2);
	e01 = find_edge_by_bids(&list, hk_diploid_bid(0, 0), hk_diploid_bid(5, 1));
	e10 = find_edge_by_bids(&list, hk_diploid_bid(0, 1), hk_diploid_bid(5, 0));
	e11 = find_edge_by_bids(&list, hk_diploid_bid(0, 1), hk_diploid_bid(5, 1));
	failed |= check_true("raw soft top2 kept top edge", e01 >= 0);
	failed |= check_true("raw soft top2 kept one tied second edge", e10 >= 0 || e11 >= 0);
	if (e01 >= 0)
		failed |= check_close("raw soft top2 top k", list.edges[e01].k, 0.62168223f);
	hk_blind_wedge_list_destroy(&list);

	set->raw[0] = make_raw_split_pair(5500000, 0);
	hk_blind_raw2binned_set(&set->raw2binned[0], 0, 1);
	hk_blind_wedge_list_init(&list);
	ret = hk_blind_wedge_list_build_mstep_graph(&list, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_TOP1,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("raw soft reversed top1 ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&list);
	failed |= check_i32("raw soft reversed top1 aggregate ret", ret, 0);
	e01 = find_edge_by_bids(&list, hk_diploid_bid(0, 0), hk_diploid_bid(5, 1));
	e10 = find_edge_by_bids(&list, hk_diploid_bid(0, 1), hk_diploid_bid(5, 0));
	failed |= check_true("raw soft reversed canonical top edge", e01 >= 0);
	failed |= check_true("raw soft reversed no raw-order swapped edge", e10 < 0);
	if (e01 >= 0)
		failed |= check_close("raw soft reversed canonical top k",
							  list.edges[e01].k, 0.62168223f);
	hk_blind_wedge_list_destroy(&list);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_mstep_filtered_soft_native_support(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[8];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list top, filtered1, filtered2, weighted1, weighted2, sample1, weighted_sample1;
	struct hk_blind_wedge_list bernoulli, bernoulli_uniform, sample_gate;
	int failed = 0, ret, i, e01, e11, ew01, ew11;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	for (i = 0; i < 4; ++i) {
		raw[i] = make_raw_split_pair(i * 1000, 5500000 - i * 1000);
		raw[i + 4] = make_raw_split_pair(5500000 - i * 1000, i * 1000);
	}
	set = hk_blind_bpair_set_build(&bmap, 8, raw);
	failed |= check_true("filtered raw split set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 1.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
	}

	hk_blind_wedge_list_init(&top);
	ret = hk_blind_wedge_list_build_mstep_graph(&top, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_TOP_ONLY,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_TOP,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("filtered top baseline build ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&top);
	failed |= check_i32("filtered top baseline aggregate ret", ret, 0);

	hk_blind_wedge_list_init(&filtered1);
	ret = hk_blind_wedge_list_build_mstep_graph(&filtered1, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_TOP1,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("filtered top1 build ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&filtered1);
	failed |= check_i32("filtered top1 aggregate ret", ret, 0);
	failed |= check_i32("filtered top1 edge count matches hard top",
						filtered1.n_edges, top.n_edges);
	if (top.n_edges == 1 && filtered1.n_edges == 1) {
		failed |= check_i32("filtered top1 bid0", filtered1.edges[0].bid[0], top.edges[0].bid[0]);
		failed |= check_i32("filtered top1 bid1", filtered1.edges[0].bid[1], top.edges[0].bid[1]);
		failed |= check_close("filtered top1 k", filtered1.edges[0].k, top.edges[0].k);
		failed |= check_close("filtered top1 d_scale", filtered1.edges[0].d_scale, top.edges[0].d_scale);
		failed |= check_i32("filtered top1 state mask", filtered1.edges[0].state_mask, 0xff);
	}
		hk_blind_wedge_list_destroy(&top);
		hk_blind_wedge_list_destroy(&filtered1);

		hk_blind_wedge_list_init(&top);
		ret = hk_blind_wedge_list_build_mstep_graph(&top, &bmap, set, 1.0f,
													HK_BLIND_RHO_TRAIN_CONSTANT,
													0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
													1e-6f, 1.0f, 1.0f,
													HK_BLIND_STATE_WEIGHT_TOP_ONLY,
													HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_TOP,
													0.0f, 0.0f, 1.0f,
													HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
													0.0f, 1.0f, 1.0f, 0.0f);
		failed |= check_i32("sample1 hard baseline build ret", ret, 0);
		ret = hk_blind_wedge_list_aggregate_exact(&top);
		failed |= check_i32("sample1 hard baseline aggregate ret", ret, 0);
		hk_blind_wedge_list_init(&sample1);
		ret = hk_blind_wedge_list_build_mstep_graph(&sample1, &bmap, set, 1.0f,
													HK_BLIND_RHO_TRAIN_CONSTANT,
													0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
													1e-6f, 1.0f, 1.0f,
													HK_BLIND_STATE_WEIGHT_POSTERIOR,
													HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1,
													0.0f, 0.0f, 1.0f,
													HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
													0.0f, 1.0f, 1.0f, 0.0f);
		failed |= check_i32("filtered sample1 one-hot build ret", ret, 0);
		ret = hk_blind_wedge_list_aggregate_exact(&sample1);
		failed |= check_i32("filtered sample1 one-hot aggregate ret", ret, 0);
		failed |= check_i32("filtered sample1 one-hot edge count",
							sample1.n_edges, top.n_edges);
		if (top.n_edges == 1 && sample1.n_edges == 1) {
			failed |= check_i32("filtered sample1 bid0", sample1.edges[0].bid[0], top.edges[0].bid[0]);
			failed |= check_i32("filtered sample1 bid1", sample1.edges[0].bid[1], top.edges[0].bid[1]);
			failed |= check_close("filtered sample1 k", sample1.edges[0].k, top.edges[0].k);
			failed |= check_close("filtered sample1 d_scale", sample1.edges[0].d_scale, top.edges[0].d_scale);
			failed |= check_i32("filtered sample1 state mask", sample1.edges[0].state_mask, 0xff);
		}
			hk_blind_wedge_list_destroy(&sample1);
			hk_blind_wedge_list_destroy(&top);

			hk_blind_wedge_list_init(&top);
			ret = hk_blind_wedge_list_build_mstep_graph(&top, &bmap, set, 1.0f,
														HK_BLIND_RHO_TRAIN_CONSTANT,
														0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
														1e-6f, 1.0f, 1.0f,
														HK_BLIND_STATE_WEIGHT_TOP_ONLY,
														HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_TOP,
														0.0f, 0.0f, 1.0f,
														HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
														0.0f, 1.0f, 1.0f, 0.0f);
			failed |= check_i32("bernoulli hard baseline build ret", ret, 0);
			ret = hk_blind_wedge_list_aggregate_exact(&top);
			failed |= check_i32("bernoulli hard baseline aggregate ret", ret, 0);
			hk_blind_wedge_list_init(&bernoulli);
			ret = hk_blind_wedge_list_build_mstep_graph(&bernoulli, &bmap, set, 1.0f,
														HK_BLIND_RHO_TRAIN_CONSTANT,
														0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
														1e-6f, 1.0f, 1.0f,
														HK_BLIND_STATE_WEIGHT_POSTERIOR,
														HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI,
														0.0f, 0.0f, 1.0f,
														HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
														0.0f, 1.0f, 1.0f, 0.0f);
			failed |= check_i32("filtered bernoulli one-hot build ret", ret, 0);
			ret = hk_blind_wedge_list_aggregate_exact(&bernoulli);
			failed |= check_i32("filtered bernoulli one-hot aggregate ret", ret, 0);
			failed |= check_i32("filtered bernoulli one-hot edge count",
								bernoulli.n_edges, top.n_edges);
			if (top.n_edges == 1 && bernoulli.n_edges == 1) {
				failed |= check_i32("filtered bernoulli bid0", bernoulli.edges[0].bid[0], top.edges[0].bid[0]);
				failed |= check_i32("filtered bernoulli bid1", bernoulli.edges[0].bid[1], top.edges[0].bid[1]);
				failed |= check_close("filtered bernoulli k", bernoulli.edges[0].k, top.edges[0].k);
				failed |= check_close("filtered bernoulli d_scale", bernoulli.edges[0].d_scale, top.edges[0].d_scale);
				failed |= check_i32("filtered bernoulli state mask", bernoulli.edges[0].state_mask, 0xff);
			}
			hk_blind_wedge_list_destroy(&bernoulli);
			hk_blind_wedge_list_destroy(&top);

		hk_blind_wedge_list_init(&weighted_sample1);
		ret = hk_blind_wedge_list_build_mstep_graph(&weighted_sample1, &bmap, set, 1.0f,
													HK_BLIND_RHO_TRAIN_CONSTANT,
													0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
													1e-6f, 1.0f, 1.0f,
													HK_BLIND_STATE_WEIGHT_POSTERIOR,
													HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_SAMPLE1,
													0.0f, 0.0f, 1.0f,
													HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
													0.0f, 1.0f, 1.0f, 0.0f);
		failed |= check_i32("filtered weighted sample1 one-hot build ret", ret, 0);
		ret = hk_blind_wedge_list_aggregate_exact(&weighted_sample1);
		failed |= check_i32("filtered weighted sample1 one-hot aggregate ret", ret, 0);
		failed |= check_i32("filtered weighted sample1 one-hot edge count", weighted_sample1.n_edges, 1);
		if (weighted_sample1.n_edges == 1) {
			failed |= check_close("filtered weighted sample1 one-hot k", weighted_sample1.edges[0].k, 1.0f);
			failed |= check_i32("filtered weighted sample1 state mask", weighted_sample1.edges[0].state_mask, 0xff);
		}
			hk_blind_wedge_list_destroy(&weighted_sample1);

			for (i = 0; i < set->n_bpairs; ++i) {
				set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.25f;
			set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.25f;
			set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.25f;
			set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.25f;
		set->bpairs[i].real_log_norm = 0.0f;
		}
		hk_blind_wedge_list_init(&sample_gate);
		ret = hk_blind_wedge_list_build_mstep_graph(&sample_gate, &bmap, set, 1.0f,
													HK_BLIND_RHO_TRAIN_CONSTANT,
													0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
													1e-6f, 1.0f, 1.0f,
													HK_BLIND_STATE_WEIGHT_POSTERIOR,
													HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1,
													0.10f, 0.50f, 1.0f,
													HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
													0.0f, 1.0f, 1.0f, 0.0f);
		failed |= check_i32("filtered sample1 gate build ret", ret, 0);
		ret = hk_blind_wedge_list_aggregate_exact(&sample_gate);
		failed |= check_i32("filtered sample1 gate aggregate ret", ret, 0);
		failed |= check_i32("filtered sample1 gate drops uniform", sample_gate.n_edges, 0);
			hk_blind_wedge_list_destroy(&sample_gate);

			hk_blind_wedge_list_init(&bernoulli_uniform);
			ret = hk_blind_wedge_list_build_mstep_graph(&bernoulli_uniform, &bmap, set, 1.0f,
														HK_BLIND_RHO_TRAIN_CONSTANT,
														0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
														1e-6f, 1.0f, 1.0f,
														HK_BLIND_STATE_WEIGHT_POSTERIOR,
														HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI,
														0.0f, 0.0f, 1.0f,
														HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
														0.0f, 1.0f, 1.0f, 0.0f);
			failed |= check_i32("filtered bernoulli uniform build ret", ret, 0);
			failed |= check_greater("filtered bernoulli uniform selected some states",
									(float)bernoulli_uniform.n_split_selected_raw, 0.0f);
			failed |= check_true("filtered bernoulli uniform selected fewer than all states",
								 bernoulli_uniform.n_split_selected_raw < 32);
			failed |= check_true("filtered bernoulli uniform selected multiple states",
								 bernoulli_uniform.n_split_state_raw_count[HK_BLIND_STATE_00] +
								 bernoulli_uniform.n_split_state_raw_count[HK_BLIND_STATE_01] +
								 bernoulli_uniform.n_split_state_raw_count[HK_BLIND_STATE_10] +
								 bernoulli_uniform.n_split_state_raw_count[HK_BLIND_STATE_11] ==
								 bernoulli_uniform.n_split_selected_raw);
			failed |= check_i32("filtered bernoulli uniform gate skip",
								(int32_t)bernoulli_uniform.n_split_gate_skip_raw, 0);
			hk_blind_wedge_list_destroy(&bernoulli_uniform);

			for (i = 0; i < set->n_bpairs; ++i) {
			set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.0f;
			set->bpairs[i].p4[HK_BLIND_STATE_01] = 1.0f;
			set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.0f;
			set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
		}
		hk_blind_wedge_list_init(&weighted1);
	ret = hk_blind_wedge_list_build_mstep_graph(&weighted1, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP1,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("filtered weighted top1 build ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&weighted1);
	failed |= check_i32("filtered weighted top1 aggregate ret", ret, 0);
	failed |= check_i32("filtered weighted top1 edge count", weighted1.n_edges, 1);
	if (weighted1.n_edges == 1) {
		failed |= check_i32("filtered weighted top1 state mask", weighted1.edges[0].state_mask, 0xff);
		failed |= check_close("filtered weighted top1 one-hot k", weighted1.edges[0].k, 1.0f);
	}
	hk_blind_wedge_list_destroy(&weighted1);

	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.60f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.40f;
		set->bpairs[i].real_log_norm = 0.0f;
	}
	hk_blind_wedge_list_init(&filtered2);
	ret = hk_blind_wedge_list_build_mstep_graph(&filtered2, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_TOP2,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("filtered top2 build ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&filtered2);
	failed |= check_i32("filtered top2 aggregate ret", ret, 0);
	e01 = find_edge_by_bids(&filtered2, hk_diploid_bid(0, 0), hk_diploid_bid(5, 1));
	e11 = find_edge_by_bids(&filtered2, hk_diploid_bid(0, 1), hk_diploid_bid(5, 1));
	failed |= check_true("filtered top2 kept first candidate", e01 >= 0);
	failed |= check_true("filtered top2 kept second candidate", e11 >= 0);
	if (e01 >= 0 && e11 >= 0) {
		failed |= check_greater("filtered top2 first native k positive", filtered2.edges[e01].k, 0.0f);
		failed |= check_greater("filtered top2 second native k positive", filtered2.edges[e11].k, 0.0f);
	}

	hk_blind_wedge_list_init(&weighted2);
	ret = hk_blind_wedge_list_build_mstep_graph(&weighted2, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_WEIGHTED_TOP2,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("filtered weighted top2 build ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&weighted2);
	failed |= check_i32("filtered weighted top2 aggregate ret", ret, 0);
	ew01 = find_edge_by_bids(&weighted2, hk_diploid_bid(0, 0), hk_diploid_bid(5, 1));
	ew11 = find_edge_by_bids(&weighted2, hk_diploid_bid(0, 1), hk_diploid_bid(5, 1));
	failed |= check_true("filtered weighted top2 kept first candidate", ew01 >= 0);
	failed |= check_true("filtered weighted top2 kept second candidate", ew11 >= 0);
	if (e01 >= 0 && e11 >= 0 && ew01 >= 0 && ew11 >= 0) {
		failed |= check_greater("unweighted first exceeds weighted first",
								filtered2.edges[e01].k, weighted2.edges[ew01].k);
		failed |= check_greater("unweighted second exceeds weighted second",
								filtered2.edges[e11].k, weighted2.edges[ew11].k);
	}
	hk_blind_wedge_list_destroy(&weighted2);
	hk_blind_wedge_list_destroy(&filtered2);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_sample1_normalizes_p4_weights(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 7000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[7];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[40];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list a, b, sample1_full, sample1_conf;
	int failed = 0, ret, i;
	int64_t n_a, n_b;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	bmap.n_beads = 7;
	beads[6].chr = 0;
	beads[6].st = 6000000;
	beads[6].en = 7000000;
	offcnt[0] = (uint64_t)0 << 32 | 7u;
	offcnt[1] = (uint64_t)7 << 32 | 0u;
	for (i = 0; i < 20; ++i) {
		raw[i] = make_raw_split_pair(i * 1000, 5500000 - i * 1000);
		raw[i + 20] = make_raw_split_pair(1000000 + i * 1000, 6500000 - i * 1000);
	}
	set = hk_blind_bpair_set_build(&bmap, 40, raw);
	failed |= check_true("sample1 p4 normalization set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 2.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 1.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 1.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
	}
	hk_blind_wedge_list_init(&a);
	ret = hk_blind_wedge_list_build_mstep_graph(&a, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("sample1 unnormalized build ret", ret, 0);
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.50f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
	}
	hk_blind_wedge_list_init(&b);
	ret = hk_blind_wedge_list_build_mstep_graph(&b, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("sample1 normalized build ret", ret, 0);
	failed |= check_i32("sample1 selected raw invariant",
						(int32_t)a.n_split_selected_raw, (int32_t)b.n_split_selected_raw);
	n_a = a.n_split_state_raw_count[HK_BLIND_STATE_00] +
		a.n_split_state_raw_count[HK_BLIND_STATE_01] +
		a.n_split_state_raw_count[HK_BLIND_STATE_10] +
		a.n_split_state_raw_count[HK_BLIND_STATE_11];
	n_b = b.n_split_state_raw_count[HK_BLIND_STATE_00] +
		b.n_split_state_raw_count[HK_BLIND_STATE_01] +
		b.n_split_state_raw_count[HK_BLIND_STATE_10] +
		b.n_split_state_raw_count[HK_BLIND_STATE_11];
	failed |= check_i32("sample1 unnormalized state sum", (int32_t)n_a, (int32_t)a.n_split_selected_raw);
	failed |= check_i32("sample1 normalized state sum", (int32_t)n_b, (int32_t)b.n_split_selected_raw);
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		failed |= check_i32("sample1 proportional state counts",
							(int32_t)a.n_split_state_raw_count[i],
							(int32_t)b.n_split_state_raw_count[i]);
	failed |= check_i32("sample1 zero-prob state11",
						(int32_t)a.n_split_state_raw_count[HK_BLIND_STATE_11], 0);
	hk_blind_wedge_list_destroy(&a);
	hk_blind_wedge_list_destroy(&b);

	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.50f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
		set->bpairs[i].entropy = 0.75f;
		set->bpairs[i].pmax = 0.5f;
	}
	hk_blind_wedge_list_init(&sample1_full);
	ret = hk_blind_wedge_list_build_mstep_graph(&sample1_full, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
												0.25f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("sample1 conf baseline build ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&sample1_full);
	failed |= check_i32("sample1 conf baseline aggregate ret", ret, 0);
	hk_blind_wedge_list_init(&sample1_conf);
	ret = hk_blind_wedge_list_build_mstep_graph(&sample1_conf, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1_CONF_KWEIGHT,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
												0.25f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("sample1 conf kweight build ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&sample1_conf);
	failed |= check_i32("sample1 conf kweight aggregate ret", ret, 0);
	failed |= check_i32("sample1 conf kweight is stochastic",
						(int32_t)sample1_conf.n_split_selected_raw,
						(int32_t)sample1_full.n_split_selected_raw);
	failed |= check_greater("sample1 full k exceeds conf k",
							sum_list_edge_k(&sample1_full),
							sum_list_edge_k(&sample1_conf));
	failed |= check_greater("sample1 conf keeps positive k",
							sum_list_edge_k(&sample1_conf), 0.0f);
	failed |= check_greater("sample1 conf dscale raw-count style",
							sum_list_edge_d_scale(&sample1_conf),
							0.0f);
	hk_blind_wedge_list_destroy(&sample1_conf);
	hk_blind_wedge_list_destroy(&sample1_full);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_sample1_raw_accounting(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 7000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[7];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[4];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list list;
	int failed = 0, ret, i;
	int32_t locked_bpair, gate_bpair;
	int64_t state_sum;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	bmap.n_beads = 7;
	beads[6].chr = 0;
	beads[6].st = 6000000;
	beads[6].en = 7000000;
	offcnt[0] = (uint64_t)0 << 32 | 7u;
	offcnt[1] = (uint64_t)7 << 32 | 0u;
	raw[0] = make_raw_split_pair(0, 5500000);
	raw[1] = make_raw_split_pair(1000000, 6500000);
	raw[2] = make_raw_split_pair(2000000, 2500000);
	raw[3] = make_raw_split_pair(3000000, 6500000);
	set = hk_blind_bpair_set_build(&bmap, 4, raw);
	failed |= check_true("sample1 accounting set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.25f;
		set->bpairs[i].real_log_norm = 0.0f;
	}
	set->raw_locked_state[0] = HK_BLIND_RAW_LOCKED_STATE_SKIP;
	set->raw_locked_state[3] = HK_BLIND_STATE_01;
	locked_bpair = set->raw2binned[3].bpair_id;
	gate_bpair = set->raw2binned[1].bpair_id;
	set->bpairs[locked_bpair].p4[HK_BLIND_STATE_00] = 0.0f;
	set->bpairs[locked_bpair].p4[HK_BLIND_STATE_01] = 0.0f;
	set->bpairs[locked_bpair].p4[HK_BLIND_STATE_10] = 1.0f;
	set->bpairs[locked_bpair].p4[HK_BLIND_STATE_11] = 0.0f;
	set->bpairs[gate_bpair].p4[HK_BLIND_STATE_00] = 0.25f;
	set->bpairs[gate_bpair].p4[HK_BLIND_STATE_01] = 0.25f;
	set->bpairs[gate_bpair].p4[HK_BLIND_STATE_10] = 0.25f;
	set->bpairs[gate_bpair].p4[HK_BLIND_STATE_11] = 0.25f;
	hk_blind_wedge_list_init(&list);
	ret = hk_blind_wedge_list_build_mstep_graph(&list, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_SAMPLE1,
												0.10f, 0.50f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("sample1 accounting build ret", ret, 0);
	failed |= check_i32("sample1 accounting selected raw", (int32_t)list.n_split_selected_raw, 1);
	failed |= check_i32("sample1 accounting gate skip", (int32_t)list.n_split_gate_skip_raw, 1);
	failed |= check_i32("sample1 accounting same-bin skip", (int32_t)list.n_split_same_bin_skip_raw, 1);
	failed |= check_i32("sample1 accounting locked skip", (int32_t)list.n_split_locked_skip_raw, 1);
	failed |= check_i32("sample1 accounting total",
						(int32_t)(list.n_split_selected_raw + list.n_split_gate_skip_raw +
								  list.n_split_same_bin_skip_raw + list.n_split_locked_skip_raw), 4);
	state_sum = list.n_split_state_raw_count[HK_BLIND_STATE_00] +
		list.n_split_state_raw_count[HK_BLIND_STATE_01] +
		list.n_split_state_raw_count[HK_BLIND_STATE_10] +
		list.n_split_state_raw_count[HK_BLIND_STATE_11];
	failed |= check_i32("sample1 accounting state sum", (int32_t)state_sum,
						(int32_t)list.n_split_selected_raw);
	failed |= check_i32("sample1 accounting locked selected state",
						(int32_t)list.n_split_state_raw_count[HK_BLIND_STATE_01], 1);
	hk_blind_wedge_list_destroy(&list);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_pcut_top1_filters_weak_states(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[8];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list list;
	int failed = 0, ret, i;
	int64_t state_sum;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	for (i = 0; i < 4; ++i) {
		raw[i] = make_raw_split_pair(i * 1000, 5500000 - i * 1000);
		raw[i + 4] = make_raw_split_pair(5500000 - i * 1000, i * 1000);
	}
	set = hk_blind_bpair_set_build(&bmap, 8, raw);
	failed |= check_true("pcut top1 set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.10f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.55f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.30f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.05f;
		set->bpairs[i].real_log_norm = 0.0f;
	}
	hk_blind_wedge_list_init(&list);
	ret = hk_blind_wedge_list_build_mstep_graph(&list, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.26f);
	failed |= check_i32("pcut high build ret", ret, 0);
	failed |= check_i32("pcut high selected raw", (int32_t)list.n_split_selected_raw, 16);
	failed |= check_i32("pcut high state00 filtered",
						(int32_t)list.n_split_state_raw_count[HK_BLIND_STATE_00], 0);
	failed |= check_i32("pcut high state01 kept",
						(int32_t)list.n_split_state_raw_count[HK_BLIND_STATE_01], 8);
	failed |= check_i32("pcut high state10 kept",
						(int32_t)list.n_split_state_raw_count[HK_BLIND_STATE_10], 8);
	failed |= check_i32("pcut high state11 filtered",
						(int32_t)list.n_split_state_raw_count[HK_BLIND_STATE_11], 0);
	hk_blind_wedge_list_destroy(&list);

	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.25f;
		set->bpairs[i].real_log_norm = 0.0f;
	}
	hk_blind_wedge_list_init(&list);
	ret = hk_blind_wedge_list_build_mstep_graph(&list, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.26f);
	failed |= check_i32("pcut uniform build ret", ret, 0);
	failed |= check_i32("pcut uniform fallback selected raw",
						(int32_t)list.n_split_selected_raw, 8);
	state_sum = list.n_split_state_raw_count[HK_BLIND_STATE_00] +
		list.n_split_state_raw_count[HK_BLIND_STATE_01] +
		list.n_split_state_raw_count[HK_BLIND_STATE_10] +
		list.n_split_state_raw_count[HK_BLIND_STATE_11];
	failed |= check_i32("pcut uniform state sum", (int32_t)state_sum,
						(int32_t)list.n_split_selected_raw);
	hk_blind_wedge_list_destroy(&list);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_pcut_top1_renorm_preserves_selected_contact_weight(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[8];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list pcut, renorm;
	float pcut_sum, renorm_sum;
	int64_t state_sum;
	int failed = 0, ret, i;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	for (i = 0; i < 4; ++i) {
		raw[i] = make_raw_split_pair(i * 1000, 5500000 - i * 1000);
		raw[i + 4] = make_raw_split_pair(5500000 - i * 1000, i * 1000);
	}
	set = hk_blind_bpair_set_build(&bmap, 8, raw);
	failed |= check_true("pcut renorm set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.10f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.55f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.30f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.05f;
		set->bpairs[i].real_log_norm = 0.0f;
	}

	hk_blind_wedge_list_init(&pcut);
	ret = hk_blind_wedge_list_build_mstep_graph(&pcut, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.26f);
	failed |= check_i32("pcut renorm baseline ret", ret, 0);
	failed |= check_i32("pcut renorm baseline selected raw",
						(int32_t)pcut.n_split_selected_raw, 16);
	ret = hk_blind_wedge_list_aggregate_exact(&pcut);
	failed |= check_i32("pcut renorm baseline aggregate ret", ret, 0);
	pcut_sum = sum_list_edge_k(&pcut);

	hk_blind_wedge_list_init(&renorm);
	ret = hk_blind_wedge_list_build_mstep_graph(&renorm, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1_RENORM,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.26f);
	failed |= check_i32("pcut renorm build ret", ret, 0);
	failed |= check_i32("pcut renorm selected raw",
						(int32_t)renorm.n_split_selected_raw, 16);
	ret = hk_blind_wedge_list_aggregate_exact(&renorm);
	failed |= check_i32("pcut renorm aggregate ret", ret, 0);
	renorm_sum = sum_list_edge_k(&renorm);
	failed |= check_greater("pcut selected graph survives filtering", pcut_sum, 0.0f);
	failed |= check_greater("pcut renorm increases selected contact weight",
							renorm_sum, pcut_sum);
	hk_blind_wedge_list_destroy(&renorm);
	hk_blind_wedge_list_destroy(&pcut);

	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.25f;
		set->bpairs[i].real_log_norm = 0.0f;
	}
	hk_blind_wedge_list_init(&renorm);
	ret = hk_blind_wedge_list_build_mstep_graph(&renorm, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_PCUT_TOP1_RENORM,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.26f);
	failed |= check_i32("pcut renorm uniform fallback ret", ret, 0);
	failed |= check_i32("pcut renorm uniform fallback selected raw",
						(int32_t)renorm.n_split_selected_raw, 8);
	state_sum = renorm.n_split_state_raw_count[HK_BLIND_STATE_00] +
		renorm.n_split_state_raw_count[HK_BLIND_STATE_01] +
		renorm.n_split_state_raw_count[HK_BLIND_STATE_10] +
		renorm.n_split_state_raw_count[HK_BLIND_STATE_11];
	failed |= check_i32("pcut renorm uniform fallback state sum",
						(int32_t)state_sum, (int32_t)renorm.n_split_selected_raw);
	ret = hk_blind_wedge_list_aggregate_exact(&renorm);
	failed |= check_i32("pcut renorm uniform fallback aggregate ret", ret, 0);
	failed |= check_greater("pcut renorm uniform fallback weak positive k",
							sum_list_edge_k(&renorm), 0.0f);
	failed |= check_greater("pcut renorm keeps fallback weak",
							2.1f, sum_list_edge_k(&renorm));
	hk_blind_wedge_list_destroy(&renorm);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_filtered_all_conf_downweights_ambiguous_states(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[8];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list all, conf;
	float all_sum, conf_sum;
	int failed = 0, ret, i;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	for (i = 0; i < 4; ++i) {
		raw[i] = make_raw_split_pair(i * 1000, 5500000 - i * 1000);
		raw[i + 4] = make_raw_split_pair(5500000 - i * 1000, i * 1000);
	}
	set = hk_blind_bpair_set_build(&bmap, 8, raw);
	failed |= check_true("filtered all conf set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.25f;
		set->bpairs[i].real_log_norm = 0.0f;
		set->bpairs[i].entropy = logf(4.0f);
		set->bpairs[i].pmax = 0.25f;
	}

	hk_blind_wedge_list_init(&all);
	ret = hk_blind_wedge_list_build_mstep_graph(&all, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("filtered all conf baseline ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&all);
	failed |= check_i32("filtered all conf baseline aggregate ret", ret, 0);
	all_sum = sum_list_edge_k(&all);

	hk_blind_wedge_list_init(&conf);
	ret = hk_blind_wedge_list_build_mstep_graph(&conf, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
												0.1f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("filtered all conf weak ret", ret, 0);
	failed |= check_i32("filtered all conf selected raw",
						(int32_t)conf.n_split_selected_raw, 32);
	ret = hk_blind_wedge_list_aggregate_exact(&conf);
	failed |= check_i32("filtered all conf weak aggregate ret", ret, 0);
	conf_sum = sum_list_edge_k(&conf);
	failed |= check_greater("filtered all conf baseline positive", all_sum, 0.0f);
	failed |= check_greater("filtered all conf weak remains positive", conf_sum, 0.0f);
	failed |= check_greater("filtered all conf downweights ambiguous", all_sum, conf_sum * 5.0f);
	hk_blind_wedge_list_destroy(&conf);
	hk_blind_wedge_list_destroy(&all);

	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 1.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
		set->bpairs[i].entropy = 0.0f;
		set->bpairs[i].pmax = 1.0f;
	}
	hk_blind_wedge_list_init(&all);
	ret = hk_blind_wedge_list_build_mstep_graph(&all, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("filtered all conf certain baseline ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&all);
	failed |= check_i32("filtered all conf certain baseline aggregate ret", ret, 0);
	all_sum = sum_list_edge_k(&all);
	hk_blind_wedge_list_init(&conf);
	ret = hk_blind_wedge_list_build_mstep_graph(&conf, &bmap, set, 1.0f,
												HK_BLIND_RHO_TRAIN_CONSTANT,
												0.0f, HK_BLIND_D_SCALE_RAW_COUNT,
												1e-6f, 1.0f, 1.0f,
												HK_BLIND_STATE_WEIGHT_POSTERIOR,
												HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF,
												0.0f, 0.0f, 1.0f,
												HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
												0.1f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("filtered all conf certain ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&conf);
	failed |= check_i32("filtered all conf certain aggregate ret", ret, 0);
	conf_sum = sum_list_edge_k(&conf);
	failed |= check_close("filtered all conf keeps confident weight", conf_sum, all_sum);
	hk_blind_wedge_list_destroy(&conf);
	hk_blind_wedge_list_destroy(&all);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_bernoulli_conf_thins_before_filter(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[8];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list hard, thinned;
	int failed = 0, ret, i;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	for (i = 0; i < 4; ++i) {
		raw[i] = make_raw_split_pair(i * 1000, 5500000 - i * 1000);
		raw[i + 4] = make_raw_split_pair(5500000 - i * 1000, i * 1000);
	}
	set = hk_blind_bpair_set_build(&bmap, 8, raw);
	failed |= check_true("bernoulli conf set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.25f;
		set->bpairs[i].real_log_norm = 0.0f;
		set->bpairs[i].entropy = logf(4.0f);
		set->bpairs[i].pmax = 0.25f;
	}
	hk_blind_wedge_list_init(&thinned);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&thinned, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
		0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("bernoulli conf uniform build ret", ret, 0);
	failed |= check_i32("bernoulli conf uniform no selected states",
						(int32_t)thinned.n_split_selected_raw, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&thinned);
	failed |= check_i32("bernoulli conf uniform aggregate ret", ret, 0);
	failed |= check_i32("bernoulli conf uniform no edges", thinned.n_edges, 0);
	hk_blind_wedge_list_destroy(&thinned);

	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 1.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
		set->bpairs[i].entropy = 0.0f;
		set->bpairs[i].pmax = 1.0f;
	}
	hk_blind_wedge_list_init(&hard);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&hard, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_TOP_ONLY,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_TOP,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
		0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("bernoulli conf hard baseline ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&hard);
	failed |= check_i32("bernoulli conf hard baseline aggregate ret", ret, 0);
	hk_blind_wedge_list_init(&thinned);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&thinned, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
		0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("bernoulli conf one-hot ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&thinned);
	failed |= check_i32("bernoulli conf one-hot aggregate ret", ret, 0);
	failed |= check_i32("bernoulli conf one-hot edge count",
						thinned.n_edges, hard.n_edges);
	failed |= check_close("bernoulli conf one-hot sum k",
						  sum_list_edge_k(&thinned), sum_list_edge_k(&hard));
	hk_blind_wedge_list_destroy(&thinned);
	hk_blind_wedge_list_destroy(&hard);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_bernoulli_conf_weight_modes_are_separate(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[8];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list base, full, weighted;
	float base_sum, full_sum, weighted_sum;
	int failed = 0, ret, i;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	for (i = 0; i < 4; ++i) {
		raw[i] = make_raw_split_pair(i * 1000, 5500000 - i * 1000);
		raw[i + 4] = make_raw_split_pair(5500000 - i * 1000, i * 1000);
	}
	set = hk_blind_bpair_set_build(&bmap, 8, raw);
	failed |= check_true("bernoulli conf weight set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.6f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.4f;
		set->bpairs[i].real_log_norm = 0.0f;
		set->bpairs[i].entropy = 0.67301166f;
		set->bpairs[i].pmax = 0.6f;
	}
	hk_blind_wedge_list_init(&base);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&base, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
		0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("bernoulli base ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&base);
	failed |= check_i32("bernoulli base aggregate ret", ret, 0);
	base_sum = sum_list_edge_k(&base);

	hk_blind_wedge_list_init(&full);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&full, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		1.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("bernoulli conf full ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&full);
	failed |= check_i32("bernoulli conf full aggregate ret", ret, 0);
	full_sum = sum_list_edge_k(&full);

	hk_blind_wedge_list_init(&weighted);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&weighted, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_BERNOULLI_CONF_WEIGHTED,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		1.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("bernoulli conf weighted ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&weighted);
	failed |= check_i32("bernoulli conf weighted aggregate ret", ret, 0);
	weighted_sum = sum_list_edge_k(&weighted);

	failed |= check_close("bernoulli full matches base", full_sum, base_sum);
	failed |= check_greater("bernoulli full exceeds weighted", full_sum, weighted_sum);
	hk_blind_wedge_list_destroy(&weighted);
	hk_blind_wedge_list_destroy(&full);
	hk_blind_wedge_list_destroy(&base);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_weighted_filter_uses_confidence_before_bmap(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[8];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list weak, certain, full, kweight, kdweight, kdhalf, native, expected, soft_all, locked, pcut, outlier, kweight_scaled;
	int failed = 0, ret, i;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	for (i = 0; i < 4; ++i) {
		raw[i] = make_raw_split_pair(i * 1000, 5500000 - i * 1000);
		raw[i + 4] = make_raw_split_pair(5500000 - i * 1000, i * 1000);
	}
	set = hk_blind_bpair_set_build(&bmap, 8, raw);
	failed |= check_true("weighted filter set", set != 0);
	if (set == 0)
		return 1;

	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.25f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.25f;
		set->bpairs[i].real_log_norm = 0.0f;
		set->bpairs[i].entropy = logf(4.0f);
		set->bpairs[i].pmax = 0.25f;
	}
	hk_blind_wedge_list_init(&weak);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&weak, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		0.1f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("weighted filter weak ret", ret, 0);
	failed |= check_i32("weighted filter weak selected raw",
						(int32_t)weak.n_split_selected_raw, 32);
	failed |= check_greater("weighted filter weak removes support before bmap",
							32.0f, (float)weak.n_split_filter2_pairs);
	ret = hk_blind_wedge_list_aggregate_exact(&weak);
	failed |= check_i32("weighted filter weak aggregate ret", ret, 0);
	failed |= check_i32("weighted filter weak no final edges", weak.n_edges, 0);

	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 1.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
		set->bpairs[i].entropy = 0.0f;
		set->bpairs[i].pmax = 1.0f;
	}
	hk_blind_wedge_list_init(&certain);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&certain, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
		0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("weighted filter certain ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&certain);
	failed |= check_i32("weighted filter certain aggregate ret", ret, 0);
	failed |= check_greater("weighted filter certain keeps more filter2 support",
							(float)certain.n_split_filter2_pairs,
							(float)weak.n_split_filter2_pairs);
	hk_blind_wedge_list_destroy(&weak);

	hk_blind_wedge_list_init(&native);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&native, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_TOP_ONLY,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_TOP,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
		0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("weighted filter native ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&native);
	failed |= check_i32("weighted filter native aggregate ret", ret, 0);
	failed |= check_i32("weighted filter all-one candidate matches native",
						(int32_t)certain.n_split_candidate_pairs,
						(int32_t)native.n_split_candidate_pairs);
	failed |= check_i32("weighted filter all-one filter1 matches native",
						(int32_t)certain.n_split_filter1_pairs,
						(int32_t)native.n_split_filter1_pairs);
	failed |= check_i32("weighted filter all-one filter2 matches native",
						(int32_t)certain.n_split_filter2_pairs,
						(int32_t)native.n_split_filter2_pairs);
	failed |= check_i32("weighted filter all-one bmap matches native",
						(int32_t)certain.n_split_bmap_pairs,
						(int32_t)native.n_split_bmap_pairs);
	failed |= check_close("weighted filter all-one k matches native",
						  sum_list_edge_k(&certain), sum_list_edge_k(&native));
	hk_blind_wedge_list_init(&expected);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&expected, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_COUNT,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
		0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("raw expected one-hot ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&expected);
	failed |= check_i32("raw expected one-hot aggregate ret", ret, 0);
	failed |= check_i32("raw expected one-hot candidate matches native",
						(int32_t)expected.n_split_candidate_pairs,
						(int32_t)native.n_split_candidate_pairs);
	failed |= check_i32("raw expected one-hot filter1 matches native",
						(int32_t)expected.n_split_filter1_pairs,
						(int32_t)native.n_split_filter1_pairs);
	failed |= check_i32("raw expected one-hot filter2 matches native",
						(int32_t)expected.n_split_filter2_pairs,
						(int32_t)native.n_split_filter2_pairs);
	failed |= check_i32("raw expected one-hot bmap matches native",
						(int32_t)expected.n_split_bmap_pairs,
						(int32_t)native.n_split_bmap_pairs);
	failed |= check_list_same_edges_for_k_and_d_scale("raw expected one-hot per-edge",
													  &expected, &native, &native);
	hk_blind_wedge_list_destroy(&expected);
	hk_blind_wedge_list_destroy(&native);

	hk_blind_wedge_list_init(&full);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&full, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_FULL,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
		0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("weighted filter full ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&full);
	failed |= check_i32("weighted filter full aggregate ret", ret, 0);
	failed |= check_i32("weighted filter full candidate count",
						(int32_t)full.n_split_candidate_pairs,
						(int32_t)certain.n_split_candidate_pairs);
	failed |= check_i32("weighted filter full filter1 count",
						(int32_t)full.n_split_filter1_pairs,
						(int32_t)certain.n_split_filter1_pairs);
	failed |= check_i32("weighted filter full filter2 count",
						(int32_t)full.n_split_filter2_pairs,
						(int32_t)certain.n_split_filter2_pairs);
	failed |= check_i32("weighted filter full bmap count",
						(int32_t)full.n_split_bmap_pairs,
						(int32_t)certain.n_split_bmap_pairs);
	failed |= check_close("weighted filter full matches certain",
				  sum_list_edge_k(&full), sum_list_edge_k(&certain));
	hk_blind_wedge_list_destroy(&full);
	hk_blind_wedge_list_destroy(&certain);

	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 0.5f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.5f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
		set->bpairs[i].entropy = logf(2.0f);
		set->bpairs[i].pmax = 0.5f;
		set->bpairs[i].real_log_norm = 10.0f;
	}
	hk_blind_wedge_list_init(&full);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&full, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_FULL,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		1.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("weighted filter full mixed ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&full);
	failed |= check_i32("weighted filter full mixed aggregate ret", ret, 0);

	hk_blind_wedge_list_init(&certain);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&certain, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		1.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("weighted filter soft mixed ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&certain);
	failed |= check_i32("weighted filter soft mixed aggregate ret", ret, 0);

	hk_blind_wedge_list_init(&kweight);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&kweight, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
		1.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("weighted filter kweight mixed ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&kweight);
	failed |= check_i32("weighted filter kweight mixed aggregate ret", ret, 0);
	failed |= check_i32("weighted filter kweight candidate count",
						(int32_t)kweight.n_split_candidate_pairs,
						(int32_t)certain.n_split_candidate_pairs);
	failed |= check_i32("weighted filter kweight filter1 count",
						(int32_t)kweight.n_split_filter1_pairs,
						(int32_t)certain.n_split_filter1_pairs);
	failed |= check_i32("weighted filter kweight filter2 count",
						(int32_t)kweight.n_split_filter2_pairs,
						(int32_t)certain.n_split_filter2_pairs);
	failed |= check_i32("weighted filter kweight bmap count",
						(int32_t)kweight.n_split_bmap_pairs,
						(int32_t)certain.n_split_bmap_pairs);
	failed |= check_close("weighted filter kweight k matches soft",
				  sum_list_edge_k(&kweight), sum_list_edge_k(&certain));
	failed |= check_greater("weighted filter full exceeds kweight k",
					sum_list_edge_k(&full), sum_list_edge_k(&kweight));
	failed |= check_greater("weighted filter soft uses larger d_scale",
					sum_list_edge_d_scale(&certain),
					sum_list_edge_d_scale(&kweight));
		failed |= check_close("weighted filter kweight d_scale matches full",
				  sum_list_edge_d_scale(&kweight),
				  sum_list_edge_d_scale(&full));
		failed |= check_list_same_edges_for_k_and_d_scale("weighted filter kweight per-edge",
														  &kweight, &certain, &full);
		hk_blind_wedge_list_init(&kdweight);
		ret = hk_blind_wedge_list_build_mstep_graph(
			&kdweight, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
			0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
			HK_BLIND_STATE_WEIGHT_POSTERIOR,
			HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDWEIGHT,
			0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
			1.0f, 1.0f, 1.0f, 0.0f);
		failed |= check_i32("weighted filter kdweight mixed ret", ret, 0);
		ret = hk_blind_wedge_list_aggregate_exact(&kdweight);
		failed |= check_i32("weighted filter kdweight mixed aggregate ret", ret, 0);
		failed |= check_i32("weighted filter kdweight edge count", kdweight.n_edges, kweight.n_edges);
		failed |= check_close("weighted filter kdweight k matches kweight",
							  sum_list_edge_k(&kdweight), sum_list_edge_k(&kweight));
		failed |= check_greater("weighted filter kdweight uses weaker d_scale targets",
								sum_list_edge_d_scale(&kdweight),
								sum_list_edge_d_scale(&kweight));
		hk_blind_wedge_list_init(&expected);
		ret = hk_blind_wedge_list_build_mstep_graph(
			&expected, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
			0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
			HK_BLIND_STATE_WEIGHT_POSTERIOR,
			HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_COUNT,
			0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
			1.0f, 1.0f, 1.0f, 0.0f);
		failed |= check_i32("raw expected mixed ret", ret, 0);
		ret = hk_blind_wedge_list_aggregate_exact(&expected);
		failed |= check_i32("raw expected mixed aggregate ret", ret, 0);
			failed |= check_i32("raw expected mixed edge count", expected.n_edges, kdweight.n_edges);
			failed |= check_list_same_edges_for_k_and_d_scale("raw expected mixed per-edge",
														  &expected, &kdweight, &kdweight);
			hk_blind_wedge_list_init(&soft_all);
			ret = hk_blind_wedge_list_build_mstep_graph(
				&soft_all, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
				0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
				HK_BLIND_STATE_WEIGHT_POSTERIOR,
				HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_ALL,
				0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
				1.0f, 1.0f, 1.0f, 0.0f);
			failed |= check_i32("raw expected soft-all ret", ret, 0);
			ret = hk_blind_wedge_list_aggregate_exact(&soft_all);
			failed |= check_i32("raw expected soft-all aggregate ret", ret, 0);
			failed |= check_i32("raw expected soft-all edge count", soft_all.n_edges, expected.n_edges);
			failed |= check_list_same_edges_for_k_and_d_scale("raw expected soft-all per-edge",
															  &soft_all, &expected, &expected);
			hk_blind_wedge_list_destroy(&soft_all);
			hk_blind_wedge_list_init(&outlier);
			ret = hk_blind_wedge_list_build_mstep_graph_ex(
				&outlier, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
				0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
				HK_BLIND_STATE_WEIGHT_POSTERIOR,
				HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER,
				0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
				1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1,
				HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS,
				HK_BLIND_RAW_SOFT_DEFAULT_MIN_Q, 1.0f);
			failed |= check_i32("raw expected outlier good ret", ret, 0);
			failed |= check_i32("raw expected outlier good seen",
								(int32_t)outlier.raw_outlier_diag.n_seen, 8);
			failed |= check_i32("raw expected outlier good unlocked",
								(int32_t)outlier.raw_outlier_diag.n_unlocked, 8);
			failed |= check_true("raw expected outlier good real mass high",
								 outlier.raw_outlier_diag.real_mass_all > 7.99);
			failed |= check_true("raw expected outlier good qU low",
								 outlier.raw_outlier_diag.max_qU < 1e-4f);
			hk_blind_wedge_list_destroy(&outlier);
			hk_blind_wedge_list_init(&outlier);
			ret = hk_blind_wedge_list_build_mstep_graph(
				&outlier, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
				0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
				HK_BLIND_STATE_WEIGHT_POSTERIOR,
				HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER,
				0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
				1.0f, 1.0f, 1.0f, 0.0f);
			failed |= check_i32("raw expected outlier wrapper ret", ret, 0);
			failed |= check_i32("raw expected outlier wrapper seen",
								(int32_t)outlier.raw_outlier_diag.n_seen, 8);
			failed |= check_true("raw expected outlier wrapper enabled U",
								 outlier.raw_outlier_diag.real_mass_all > 7.99 &&
								 outlier.raw_outlier_diag.real_mass_all < 8.0);
			hk_blind_wedge_list_destroy(&outlier);
			for (i = 0; i < set->n_bpairs; ++i)
				set->bpairs[i].real_log_norm = -10.0f;
			hk_blind_wedge_list_init(&outlier);
			ret = hk_blind_wedge_list_build_mstep_graph_ex(
				&outlier, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
				0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
				HK_BLIND_STATE_WEIGHT_POSTERIOR,
				HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER,
				0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
				1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1,
				HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS,
				HK_BLIND_RAW_SOFT_DEFAULT_MIN_Q, 1.0f);
			failed |= check_i32("raw expected outlier bad ret", ret, 0);
			failed |= check_true("raw expected outlier bad qU high",
								 outlier.raw_outlier_diag.min_qU > 0.999f);
			failed |= check_close("raw expected outlier bad conservation",
								  (float)(outlier.raw_outlier_diag.emitted_real_mass_all +
										  outlier.raw_outlier_diag.truncated_real_mass_all +
										  outlier.raw_outlier_diag.outlier_mass_all),
								  8.0f);
			hk_blind_wedge_list_destroy(&outlier);
			hk_blind_wedge_list_init(&outlier);
			ret = hk_blind_wedge_list_build_mstep_graph_ex(
				&outlier, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
				0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
				HK_BLIND_STATE_WEIGHT_POSTERIOR,
				HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER,
				0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
				1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1,
				HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS,
				HK_BLIND_RAW_SOFT_DEFAULT_MIN_Q, 2.0f);
			failed |= check_i32("raw expected outlier T2 ret", ret, 0);
			failed |= check_greater("raw expected outlier T2 still mostly U",
									outlier.raw_outlier_diag.min_qU, 0.98f);
			hk_blind_wedge_list_destroy(&outlier);
			for (i = 0; i < set->n_bpairs; ++i)
				set->bpairs[i].real_log_norm = 10.0f;
			hk_blind_wedge_list_destroy(&expected);
			set->raw_locked_state[0] = HK_BLIND_STATE_01;
			hk_blind_wedge_list_init(&outlier);
			ret = hk_blind_wedge_list_build_mstep_graph_ex(
				&outlier, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
				0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
				HK_BLIND_STATE_WEIGHT_POSTERIOR,
				HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER,
				0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
				1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1,
				HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS,
				HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS,
				HK_BLIND_RAW_SOFT_DEFAULT_MIN_Q, 1.0f);
			failed |= check_i32("raw expected outlier locked ret", ret, 0);
			failed |= check_i32("raw expected outlier locked count",
								(int32_t)outlier.raw_outlier_diag.n_locked, 1);
			failed |= check_close("raw expected outlier locked real mass",
								  (float)outlier.raw_outlier_diag.real_mass_all,
								  8.0f);
			hk_blind_wedge_list_destroy(&outlier);
			hk_blind_wedge_list_init(&locked);
			ret = hk_blind_wedge_list_build_mstep_graph(
				&locked, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
				0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
				HK_BLIND_STATE_WEIGHT_POSTERIOR,
				HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_LOCKED_ONLY,
				0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
				1.0f, 1.0f, 1.0f, 0.0f);
			failed |= check_i32("raw expected locked-only ret", ret, 0);
			failed |= check_i32("raw expected locked-only selected",
								(int32_t)locked.n_split_selected_raw, 1);
			failed |= check_i32("raw expected locked-only gated",
								(int32_t)locked.n_split_gate_skip_raw, 7);
			ret = hk_blind_wedge_list_aggregate_exact(&locked);
			failed |= check_i32("raw expected locked-only aggregate ret", ret, 0);
			hk_blind_wedge_list_destroy(&locked);
			set->raw_locked_state[0] = HK_BLIND_RAW_LOCKED_STATE_NONE;
			hk_blind_wedge_list_init(&pcut);
			ret = hk_blind_wedge_list_build_mstep_graph(
				&pcut, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
				0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
				HK_BLIND_STATE_WEIGHT_POSTERIOR,
				HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT,
				0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
				1.0f, 1.0f, 1.0f, 0.6f);
			failed |= check_i32("raw expected pcut ret", ret, 0);
			failed |= check_i32("raw expected pcut selected",
								(int32_t)pcut.n_split_selected_raw, 0);
			failed |= check_i32("raw expected pcut gated",
								(int32_t)pcut.n_split_gate_skip_raw, 8);
			ret = hk_blind_wedge_list_aggregate_exact(&pcut);
			failed |= check_i32("raw expected pcut aggregate ret", ret, 0);
			hk_blind_wedge_list_destroy(&pcut);
		hk_blind_wedge_list_init(&kdhalf);
		ret = hk_blind_wedge_list_build_mstep_graph(
			&kdhalf, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
			0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
			HK_BLIND_STATE_WEIGHT_POSTERIOR,
			HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KDHALF,
			0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
			1.0f, 1.0f, 1.0f, 0.0f);
		failed |= check_i32("weighted filter kdhalf mixed ret", ret, 0);
		ret = hk_blind_wedge_list_aggregate_exact(&kdhalf);
		failed |= check_i32("weighted filter kdhalf mixed aggregate ret", ret, 0);
		failed |= check_i32("weighted filter kdhalf edge count", kdhalf.n_edges, kweight.n_edges);
		failed |= check_close("weighted filter kdhalf k matches kweight",
							  sum_list_edge_k(&kdhalf), sum_list_edge_k(&kweight));
		failed |= check_greater("weighted filter kdhalf above kweight d_scale",
								sum_list_edge_d_scale(&kdhalf),
								sum_list_edge_d_scale(&kweight));
		failed |= check_greater("weighted filter kdweight above kdhalf d_scale",
								sum_list_edge_d_scale(&kdweight),
								sum_list_edge_d_scale(&kdhalf));
		hk_blind_wedge_list_destroy(&kdhalf);
		hk_blind_wedge_list_destroy(&kdweight);
			hk_blind_wedge_list_init(&kweight_scaled);
			ret = hk_blind_wedge_list_build_mstep_graph(
			&kweight_scaled, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
			0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 3.0f, 1.0f,
			HK_BLIND_STATE_WEIGHT_POSTERIOR,
			HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_KWEIGHT,
			0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_FLOOR_ENTROPY,
			1.0f, 1.0f, 1.0f, 0.0f);
		failed |= check_i32("weighted filter kweight scaled ret", ret, 0);
		ret = hk_blind_wedge_list_aggregate_exact(&kweight_scaled);
		failed |= check_i32("weighted filter kweight scaled aggregate ret", ret, 0);
		failed |= check_i32("weighted filter kweight scaled edge count",
							kweight_scaled.n_edges, kweight.n_edges);
		failed |= check_close("weighted filter kweight scaled k",
					  sum_list_edge_k(&kweight_scaled),
					  3.0f * sum_list_edge_k(&kweight));
		failed |= check_close("weighted filter kweight scaled d_scale unchanged",
					  sum_list_edge_d_scale(&kweight_scaled),
					  sum_list_edge_d_scale(&kweight));
		hk_blind_wedge_list_destroy(&kweight_scaled);
	hk_blind_wedge_list_destroy(&kweight);
	hk_blind_wedge_list_destroy(&certain);
	hk_blind_wedge_list_destroy(&full);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_weighted_filter_native_radius_boundary(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[2];
	struct hk_blind_bpair_set *set;
	struct hk_blind_wedge_list native, weighted;
	int failed = 0, ret, i;

	init_raw_split_bmap(&bmap, &dict, names, len, beads, offcnt);
	raw[0] = make_raw_split_pair(100000, 2100000);
	raw[1] = make_raw_split_pair(200000, 3100000);
	set = hk_blind_bpair_set_build(&bmap, 2, raw);
	failed |= check_true("weighted filter boundary set", set != 0);
	if (set == 0)
		return 1;
	for (i = 0; i < set->n_bpairs; ++i) {
		set->bpairs[i].p4[HK_BLIND_STATE_00] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_01] = 1.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_10] = 0.0f;
		set->bpairs[i].p4[HK_BLIND_STATE_11] = 0.0f;
		set->bpairs[i].real_log_norm = 0.0f;
		set->bpairs[i].entropy = 0.0f;
		set->bpairs[i].pmax = 1.0f;
	}

	hk_blind_wedge_list_init(&native);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&native, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_TOP_ONLY,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_TOP,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
		0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("weighted filter boundary native ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&native);
	failed |= check_i32("weighted filter boundary native aggregate ret", ret, 0);

	hk_blind_wedge_list_init(&weighted);
	ret = hk_blind_wedge_list_build_mstep_graph(
		&weighted, &bmap, set, 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT,
		0.0f, HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f, 1.0f, 1.0f,
		HK_BLIND_STATE_WEIGHT_POSTERIOR,
		HK_BLIND_MSTEP_GRAPH_RAW_SPLIT_SOFT_FILTERED_ALL_CONF_WFILTER_FULL,
		0.0f, 0.0f, 1.0f, HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
		0.0f, 1.0f, 1.0f, 0.0f);
	failed |= check_i32("weighted filter boundary weighted ret", ret, 0);
	ret = hk_blind_wedge_list_aggregate_exact(&weighted);
	failed |= check_i32("weighted filter boundary weighted aggregate ret", ret, 0);
	failed |= check_i32("weighted filter boundary candidate matches native",
						(int32_t)weighted.n_split_candidate_pairs,
						(int32_t)native.n_split_candidate_pairs);
	failed |= check_i32("weighted filter boundary filter1 matches native",
						(int32_t)weighted.n_split_filter1_pairs,
						(int32_t)native.n_split_filter1_pairs);
	failed |= check_i32("weighted filter boundary filter2 matches native",
						(int32_t)weighted.n_split_filter2_pairs,
						(int32_t)native.n_split_filter2_pairs);
	failed |= check_i32("weighted filter boundary bmap matches native",
						(int32_t)weighted.n_split_bmap_pairs,
						(int32_t)native.n_split_bmap_pairs);
	failed |= check_i32("weighted filter boundary no edges", weighted.n_edges, 0);
	hk_blind_wedge_list_destroy(&weighted);
	hk_blind_wedge_list_destroy(&native);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_uniform_expansion();
	failed |= check_nonuniform_expansion();
	failed |= check_zero_rho();
	failed |= check_entropy_rho_training_mode();
	failed |= check_floor_and_cis_trans_rho_training_modes();
	failed |= check_expected_count_d_scale_mode();
	failed |= check_weighted_scalar_energy();
	failed |= check_state_weight_posterior_matches_legacy_mode();
	failed |= check_state_weight_posterior_power_gamma();
	failed |= check_state_weight_binary_support_uses_full_edges();
	failed |= check_state_weight_top_only_thresholds();
	failed |= check_fixed_p4_lock_is_not_updated();
	failed |= check_zero_weight_wedges_are_not_list_edges();
	failed |= check_aggregate_drops_zero_weight_edges();
	failed |= check_aggregate_drops_all_zero_weight_edges();
	failed |= check_mstep_raw_split_top_uses_raw_contacts();
	failed |= check_raw_split_mstep_posterior_power_preserves_mass();
	failed |= check_raw_split_kweight_final_prob_survives_sort();
	failed |= check_mstep_raw_split_soft_confidence();
	failed |= check_mstep_filtered_soft_native_support();
	failed |= check_sample1_normalizes_p4_weights();
	failed |= check_sample1_raw_accounting();
	failed |= check_pcut_top1_filters_weak_states();
	failed |= check_pcut_top1_renorm_preserves_selected_contact_weight();
	failed |= check_filtered_all_conf_downweights_ambiguous_states();
	failed |= check_bernoulli_conf_thins_before_filter();
	failed |= check_bernoulli_conf_weight_modes_are_separate();
	failed |= check_weighted_filter_uses_confidence_before_bmap();
	failed |= check_weighted_filter_native_radius_boundary();
	return failed != 0;
}
