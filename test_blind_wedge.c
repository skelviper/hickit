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

static int check_greater(const char *label, float a, float b)
{
	if (!(a > b)) {
		fprintf(stderr, "%s: expected %.8g > %.8g\n", label, a, b);
		return 1;
	}
	return 0;
}

static void set_bpair(struct hk_blind_bpair *bp, float p00, float p01, float p10, float p11)
{
	bp->key.bid[0] = 2;
	bp->key.bid[1] = 5;
	bp->n_raw = 1;
	bp->contact_class = HK_BLIND_CONTACT_CIS;
	bp->p4[HK_BLIND_STATE_00] = p00;
	bp->p4[HK_BLIND_STATE_01] = p01;
	bp->p4[HK_BLIND_STATE_10] = p10;
	bp->p4[HK_BLIND_STATE_11] = p11;
	bp->entropy = 0.0f;
	bp->pmax = 0.0f;
	bp->margin = 0.0f;
	bp->rho_output = 0.0f;
	bp->pU = 0.0f;
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
														   0.9f, 0.9f, 3.0f, new_edges);
	failed |= check_edges_match("posterior state-weight matches legacy", new_edges, old_edges);
	return failed;
}

static int check_state_weight_top_if_confident_cis_keeps_posterior(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge posterior_edges[HK_BLIND_N_STATE], top_edges[HK_BLIND_N_STATE];
	int failed = 0;

	set_bpair(&bp, 0.05f, 0.2f, 0.65f, 0.1f);
	bp.contact_class = HK_BLIND_CONTACT_CIS;
	bp.pmax = 0.65f;
	bp.margin = 0.45f;
	hk_blind_bpair_expand_weighted_edges_mode_state_weight(&bp, 2.0f, 0.7f, 1.0f,
														   HK_BLIND_RHO_TRAIN_CONSTANT,
														   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
														   HK_BLIND_D_SCALE_RAW_COUNT,
														   1e-3f, 1.0f, 1.0f,
														   HK_BLIND_STATE_WEIGHT_POSTERIOR,
														   0.4f, 0.6f, 1.0f, posterior_edges);
	hk_blind_bpair_expand_weighted_edges_mode_state_weight(&bp, 2.0f, 0.7f, 1.0f,
														   HK_BLIND_RHO_TRAIN_CONSTANT,
														   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
														   HK_BLIND_D_SCALE_RAW_COUNT,
														   1e-3f, 1.0f, 1.0f,
														   HK_BLIND_STATE_WEIGHT_TOP_IF_CONFIDENT,
														   0.4f, 0.6f, 1.0f, top_edges);
	failed |= check_edges_match("top-if-confident cis keeps posterior", top_edges, posterior_edges);
	return failed;
}

static int check_state_weight_top_if_confident_trans_gate(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge edges[HK_BLIND_N_STATE];
	float base_k = 3.0f, rho_train = 1.5f, trans_multiplier = 2.0f;
	int failed = 0;

	set_bpair(&bp, 0.1f, 0.2f, 0.6f, 0.1f);
	bp.contact_class = HK_BLIND_CONTACT_TRANS;
	bp.pmax = 0.6f;
	bp.margin = 0.4f;
	hk_blind_bpair_expand_weighted_edges_mode_state_weight(&bp, base_k, 0.7f, rho_train,
														   HK_BLIND_RHO_TRAIN_CONSTANT,
														   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
														   HK_BLIND_D_SCALE_RAW_COUNT,
														   1e-3f, 1.0f, trans_multiplier,
														   HK_BLIND_STATE_WEIGHT_TOP_IF_CONFIDENT,
														   0.3f, 0.5f, 1.0f, edges);
	failed |= check_close("top-if-confident trans pass k00", edges[HK_BLIND_STATE_00].k, 0.0f);
	failed |= check_close("top-if-confident trans pass k01", edges[HK_BLIND_STATE_01].k, 0.0f);
	failed |= check_close("top-if-confident trans pass k10", edges[HK_BLIND_STATE_10].k,
						  base_k * rho_train * trans_multiplier);
	failed |= check_close("top-if-confident trans pass k11", edges[HK_BLIND_STATE_11].k, 0.0f);
	failed |= check_close("top-if-confident trans pass sum", sum_edge_k(edges),
						  base_k * rho_train * trans_multiplier);

	bp.margin = 0.2f;
	hk_blind_bpair_expand_weighted_edges_mode_state_weight(&bp, base_k, 0.7f, rho_train,
														   HK_BLIND_RHO_TRAIN_CONSTANT,
														   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
														   HK_BLIND_D_SCALE_RAW_COUNT,
														   1e-3f, 1.0f, trans_multiplier,
														   HK_BLIND_STATE_WEIGHT_TOP_IF_CONFIDENT,
														   0.3f, 0.5f, 1.0f, edges);
	failed |= check_close("top-if-confident trans fail sum", sum_edge_k(edges), 0.0f);
	return failed;
}

static int check_state_weight_power_sharpen_trans(void)
{
	struct hk_blind_bpair bp;
	struct hk_blind_wedge edges[HK_BLIND_N_STATE];
	float p2_sum, total_k;
	int failed = 0;

	set_bpair(&bp, 0.5f, 0.25f, 0.25f, 0.0f);
	bp.contact_class = HK_BLIND_CONTACT_TRANS;
	hk_blind_bpair_expand_weighted_edges_mode_state_weight(&bp, 4.0f, 0.7f, 1.0f,
														   HK_BLIND_RHO_TRAIN_CONSTANT,
														   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
														   HK_BLIND_D_SCALE_RAW_COUNT,
														   1e-3f, 1.0f, 1.0f,
														   HK_BLIND_STATE_WEIGHT_POWER_SHARPEN,
														   0.0f, 0.0f, 2.0f, edges);
	p2_sum = 0.5f * 0.5f + 0.25f * 0.25f + 0.25f * 0.25f;
	total_k = 4.0f;
	failed |= check_close("power-sharpen k00", edges[HK_BLIND_STATE_00].k,
						  total_k * (0.5f * 0.5f) / p2_sum);
	failed |= check_close("power-sharpen k01", edges[HK_BLIND_STATE_01].k,
						  total_k * (0.25f * 0.25f) / p2_sum);
	failed |= check_close("power-sharpen k10", edges[HK_BLIND_STATE_10].k,
						  total_k * (0.25f * 0.25f) / p2_sum);
	failed |= check_close("power-sharpen k11", edges[HK_BLIND_STATE_11].k, 0.0f);
	failed |= check_close("power-sharpen conserved k", sum_edge_k(edges), total_k);
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
	failed |= check_state_weight_top_if_confident_cis_keeps_posterior();
	failed |= check_state_weight_top_if_confident_trans_gate();
	failed |= check_state_weight_power_sharpen_trans();
	return failed != 0;
}
