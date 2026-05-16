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

static int check_i64(const char *label, int64_t got, int64_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %lld, expected %lld\n", label, (long long)got, (long long)expected);
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

static void set_bpair(struct hk_blind_bpair *bp, int32_t bid0, int32_t bid1, float p00, float p01, float p10, float p11)
{
	bp->key.bid[0] = bid0;
	bp->key.bid[1] = bid1;
	bp->n_raw = 1;
	bp->base_d_scale = 1.0f;
	bp->base_k = 1.0f;
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

static void set_bpair_params(struct hk_blind_bpair *bp, float base_k, float base_d_scale)
{
	bp->base_k = base_k;
	bp->base_d_scale = base_d_scale;
}

static void set_bpair_set(struct hk_blind_bpair_set *set, struct hk_blind_bpair *bpairs, int32_t n_bpairs)
{
	set->bpairs = bpairs;
	set->n_bpairs = n_bpairs;
	set->raw2binned = 0;
	set->n_raw = 0;
}

static void set_wedge(struct hk_blind_wedge *edge, int32_t bid0, int32_t bid1, float k, float d_scale, int32_t state)
{
	edge->bid[0] = bid0;
	edge->bid[1] = bid1;
	edge->k = k;
	edge->d_scale = d_scale;
	edge->state = state;
	edge->state_mask = HK_BLIND_STATE_MASK(state);
}

static float sum_edge_k(const struct hk_blind_wedge_list *list)
{
	float sum = 0.0f;
	int32_t i;
	for (i = 0; i < list->n_edges; ++i)
		sum += list->edges[i].k;
	return sum;
}

static int check_edge_ex(const char *label, const struct hk_blind_wedge *edge, int32_t bid0, int32_t bid1,
						 float k, float d_scale, int32_t state, uint8_t state_mask)
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
	failed |= check_i32(buf, edge->state_mask, state_mask);
	return failed;
}

static int check_edge(const char *label, const struct hk_blind_wedge *edge, int32_t bid0, int32_t bid1,
					  float k, float d_scale, int32_t state)
{
	return check_edge_ex(label, edge, bid0, bid1, k, d_scale, state, HK_BLIND_STATE_MASK(state));
}

static int check_stats(const char *label, const struct hk_blind_wedge_list *list, int32_t n_edges,
					   int64_t n_input_bpair, int64_t n_expanded_edges, int64_t n_skipped_self_edges)
{
	int failed = 0;
	char buf[64];
	snprintf(buf, sizeof(buf), "%s n_edges", label);
	failed |= check_i32(buf, list->n_edges, n_edges);
	snprintf(buf, sizeof(buf), "%s n_input_bpair", label);
	failed |= check_i64(buf, list->n_input_bpair, n_input_bpair);
	snprintf(buf, sizeof(buf), "%s n_expanded_edges", label);
	failed |= check_i64(buf, list->n_expanded_edges, n_expanded_edges);
	snprintf(buf, sizeof(buf), "%s n_skipped_self_edges", label);
	failed |= check_i64(buf, list->n_skipped_self_edges, n_skipped_self_edges);
	return failed;
}

static int check_aggregation_stats(const char *label, const struct hk_blind_wedge_list *list,
								   int64_t n_edges_before, int64_t n_removed)
{
	int failed = 0;
	char buf[64];
	snprintf(buf, sizeof(buf), "%s n_edges_before_aggregation", label);
	failed |= check_i64(buf, list->n_edges_before_aggregation, n_edges_before);
	snprintf(buf, sizeof(buf), "%s n_aggregated_edges_removed", label);
	failed |= check_i64(buf, list->n_aggregated_edges_removed, n_removed);
	return failed;
}

static int check_non_diagonal_uniform(void)
{
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	int failed = 0;

	set_bpair(&bpairs[0], 2, 5, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair_set(&set, bpairs, 1);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("uniform build ret", hk_blind_wedge_list_build_from_bpair_set(&list, &set, 2.0f, 0.7f, 1.0f), 0);

	failed |= check_stats("uniform stats", &list, 4, 1, 4, 0);
	failed |= check_aggregation_stats("uniform aggregation stats", &list, 4, 0);
	failed |= check_edge("uniform 00", &list.edges[0], 4, 10, 0.5f, 0.7f, HK_BLIND_STATE_00);
	failed |= check_edge("uniform 01", &list.edges[1], 4, 11, 0.5f, 0.7f, HK_BLIND_STATE_01);
	failed |= check_edge("uniform 10", &list.edges[2], 5, 10, 0.5f, 0.7f, HK_BLIND_STATE_10);
	failed |= check_edge("uniform 11", &list.edges[3], 5, 11, 0.5f, 0.7f, HK_BLIND_STATE_11);
	failed |= check_close("uniform conserved k", sum_edge_k(&list), 2.0f);

	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_non_diagonal_nonuniform(void)
{
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	float base_k = 3.0f, rho_train = 0.5f, base_d_scale = 0.7f;
	int failed = 0;

	set_bpair(&bpairs[0], 2, 5, 0.7f, 0.1f, 0.15f, 0.05f);
	set_bpair_set(&set, bpairs, 1);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("nonuniform build ret", hk_blind_wedge_list_build_from_bpair_set(&list, &set, base_k, base_d_scale, rho_train), 0);

	failed |= check_stats("nonuniform stats", &list, 4, 1, 4, 0);
	failed |= check_aggregation_stats("nonuniform aggregation stats", &list, 4, 0);
	failed |= check_edge("nonuniform 00", &list.edges[0], 4, 10, base_k * rho_train * 0.7f, base_d_scale, HK_BLIND_STATE_00);
	failed |= check_edge("nonuniform 01", &list.edges[1], 4, 11, base_k * rho_train * 0.1f, base_d_scale, HK_BLIND_STATE_01);
	failed |= check_edge("nonuniform 10", &list.edges[2], 5, 10, base_k * rho_train * 0.15f, base_d_scale, HK_BLIND_STATE_10);
	failed |= check_edge("nonuniform 11", &list.edges[3], 5, 11, base_k * rho_train * 0.05f, base_d_scale, HK_BLIND_STATE_11);
	failed |= check_close("nonuniform conserved k", sum_edge_k(&list), base_k * rho_train);

	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_same_bin_self_edges(void)
{
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	int failed = 0;

	set_bpair(&bpairs[0], 3, 3, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair_set(&set, bpairs, 1);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("same-bin build ret", hk_blind_wedge_list_build_from_bpair_set(&list, &set, 2.0f, 0.7f, 1.0f), 0);

	failed |= check_stats("same-bin stats", &list, 0, 1, 0, 0);
	failed |= check_i64("same-bin skipped bpairs", list.n_skipped_same_bin_bpairs, 1);
	failed |= check_aggregation_stats("same-bin aggregation stats", &list, 0, 0);
	failed |= check_close("same-bin filtered k", sum_edge_k(&list), 0.0f);

	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_zero_rho_keeps_edges(void)
{
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	int failed = 0, i;

	set_bpair(&bpairs[0], 2, 5, 0.7f, 0.1f, 0.15f, 0.05f);
	set_bpair_set(&set, bpairs, 1);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("zero-rho build ret", hk_blind_wedge_list_build_from_bpair_set(&list, &set, 3.0f, 0.7f, 0.0f), 0);

	failed |= check_stats("zero-rho stats", &list, 4, 1, 4, 0);
	failed |= check_aggregation_stats("zero-rho aggregation stats", &list, 4, 0);
	for (i = 0; i < list.n_edges; ++i) {
		failed |= check_close("zero-rho k", list.edges[i].k, 0.0f);
		failed |= check_close("zero-rho d_scale", list.edges[i].d_scale, 0.7f);
	}

	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_param_aware_builder_uses_bpair_params(void)
{
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	float rho_train = 0.5f;
	int failed = 0;

	set_bpair(&bpairs[0], 2, 5, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair(&bpairs[1], 3, 7, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair_params(&bpairs[0], 2.0f, 1.0f);
	set_bpair_params(&bpairs[1], 4.0f, 0.5f);
	set_bpair_set(&set, bpairs, 2);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("params build ret", hk_blind_wedge_list_build_from_bpair_set_params(&list, &set, rho_train), 0);

	failed |= check_stats("params stats", &list, 8, 2, 8, 0);
	failed |= check_aggregation_stats("params aggregation stats", &list, 8, 0);
	failed |= check_edge("params bpair0 00", &list.edges[0], 4, 10, 0.25f, 1.0f, HK_BLIND_STATE_00);
	failed |= check_edge("params bpair0 01", &list.edges[1], 4, 11, 0.25f, 1.0f, HK_BLIND_STATE_01);
	failed |= check_edge("params bpair0 10", &list.edges[2], 5, 10, 0.25f, 1.0f, HK_BLIND_STATE_10);
	failed |= check_edge("params bpair0 11", &list.edges[3], 5, 11, 0.25f, 1.0f, HK_BLIND_STATE_11);
	failed |= check_edge("params bpair1 00", &list.edges[4], 6, 14, 0.5f, 0.5f, HK_BLIND_STATE_00);
	failed |= check_edge("params bpair1 01", &list.edges[5], 6, 15, 0.5f, 0.5f, HK_BLIND_STATE_01);
	failed |= check_edge("params bpair1 10", &list.edges[6], 7, 14, 0.5f, 0.5f, HK_BLIND_STATE_10);
	failed |= check_edge("params bpair1 11", &list.edges[7], 7, 15, 0.5f, 0.5f, HK_BLIND_STATE_11);

	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_param_aware_stiffness_conservation(void)
{
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	float rho_train = 2.0f;
	int failed = 0;

	set_bpair(&bpairs[0], 2, 5, 0.1f, 0.2f, 0.3f, 0.4f);
	set_bpair(&bpairs[1], 3, 3, 0.1f, 0.2f, 0.3f, 0.4f);
	set_bpair_params(&bpairs[0], 1.0f, 0.8f);
	set_bpair_params(&bpairs[1], 1.0f, 0.5f);
	set_bpair_set(&set, bpairs, 2);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("params stiffness ret", hk_blind_wedge_list_build_from_bpair_set_params(&list, &set, rho_train), 0);

	failed |= check_stats("params stiffness stats", &list, 4, 2, 4, 0);
	failed |= check_i64("params stiffness same-bin skipped", list.n_skipped_same_bin_bpairs, 1);
	failed |= check_edge("params stiffness 00", &list.edges[0], 4, 10, 0.2f, 0.8f, HK_BLIND_STATE_00);
	failed |= check_edge("params stiffness 01", &list.edges[1], 4, 11, 0.4f, 0.8f, HK_BLIND_STATE_01);
	failed |= check_edge("params stiffness 10", &list.edges[2], 5, 10, 0.6f, 0.8f, HK_BLIND_STATE_10);
	failed |= check_edge("params stiffness 11", &list.edges[3], 5, 11, 0.8f, 0.8f, HK_BLIND_STATE_11);
	failed |= check_close("params stiffness conserved filtered k", sum_edge_k(&list), 2.0f);

	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_param_entropy_training_mode(void)
{
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	float rho_train = 0.5f;
	int failed = 0;

	set_bpair(&bpairs[0], 2, 5, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair(&bpairs[1], 3, 7, 1.0f, 0.0f, 0.0f, 0.0f);
	set_bpair_params(&bpairs[0], 2.0f, 1.0f);
	set_bpair_params(&bpairs[1], 3.0f, 1.0f);
	bpairs[0].rho_output = 0.0f;
	bpairs[1].rho_output = 1.0f;
	set_bpair_set(&set, bpairs, 2);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("entropy params build ret",
						hk_blind_wedge_list_build_from_bpair_set_params_mode(&list, &set, rho_train,
																			  HK_BLIND_RHO_TRAIN_ENTROPY,
																			  HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f),
						0);
	failed |= check_stats("entropy params stats", &list, 8, 2, 8, 0);
	failed |= check_close("entropy params conserved confident k", sum_edge_k(&list), 3.0f * rho_train);
	failed |= check_close("entropy params mean rho", list.mean_rho_train_bpair, 0.25f);
	failed |= check_close("entropy params min rho", list.min_rho_train_bpair, 0.0f);
	failed |= check_close("entropy params max rho", list.max_rho_train_bpair, rho_train);
	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_param_expected_count_d_scale(void)
{
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	float d_high, d_low, d_zero;
	int failed = 0;

	set_bpair(&bpairs[0], 2, 5, 0.75f, 0.25f, 0.0f, 0.0f);
	bpairs[0].n_raw = 8;
	set_bpair_params(&bpairs[0], 2.0f, 0.5f);
	set_bpair_set(&set, bpairs, 1);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("expected d_scale params build ret",
						hk_blind_wedge_list_build_from_bpair_set_params_mode(&list, &set, 1.0f,
																			  HK_BLIND_RHO_TRAIN_CONSTANT,
																			  HK_BLIND_D_SCALE_EXPECTED_COUNT, 1e-3f),
						0);
	d_high = powf(8.0f * 0.75f, -1.0f / 3.0f);
	d_low = powf(8.0f * 0.25f, -1.0f / 3.0f);
	d_zero = powf(1e-3f, -1.0f / 3.0f);
	failed |= check_stats("expected d_scale stats", &list, 4, 1, 4, 0);
	failed |= check_edge("expected d_scale high", &list.edges[0], 4, 10, 1.5f, d_high, HK_BLIND_STATE_00);
	failed |= check_edge("expected d_scale low", &list.edges[1], 4, 11, 0.5f, d_low, HK_BLIND_STATE_01);
	failed |= check_edge("expected d_scale zero", &list.edges[2], 5, 10, 0.0f, d_zero, HK_BLIND_STATE_10);
	failed |= check_close("expected d_scale conserved k", sum_edge_k(&list), 2.0f);
	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_global_builder_ignores_bpair_params(void)
{
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	int failed = 0;

	set_bpair(&bpairs[0], 2, 5, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair_params(&bpairs[0], 99.0f, 0.25f);
	set_bpair_set(&set, bpairs, 1);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("global ignores params ret", hk_blind_wedge_list_build_from_bpair_set(&list, &set, 2.0f, 0.7f, 1.0f), 0);

	failed |= check_stats("global ignores params stats", &list, 4, 1, 4, 0);
	failed |= check_edge("global ignores params 00", &list.edges[0], 4, 10, 0.5f, 0.7f, HK_BLIND_STATE_00);
	failed |= check_edge("global ignores params 01", &list.edges[1], 4, 11, 0.5f, 0.7f, HK_BLIND_STATE_01);
	failed |= check_edge("global ignores params 10", &list.edges[2], 5, 10, 0.5f, 0.7f, HK_BLIND_STATE_10);
	failed |= check_edge("global ignores params 11", &list.edges[3], 5, 11, 0.5f, 0.7f, HK_BLIND_STATE_11);

	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_param_aware_exact_aggregation(void)
{
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	int failed = 0;

	set_bpair(&bpairs[0], 3, 4, 0.0f, 0.2f, 0.3f, 0.0f);
	set_bpair(&bpairs[1], 3, 4, 0.0f, 0.4f, 0.1f, 0.0f);
	set_bpair_params(&bpairs[0], 1.0f, 0.5f);
	set_bpair_params(&bpairs[1], 1.0f, 0.5f);
	set_bpair_set(&set, bpairs, 2);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("params aggregate same d_scale build ret", hk_blind_wedge_list_build_from_bpair_set_params(&list, &set, 1.0f), 0);
	failed |= check_stats("params aggregate same d_scale build stats", &list, 8, 2, 8, 0);
	failed |= check_i32("params aggregate same d_scale ret", hk_blind_wedge_list_aggregate_exact(&list), 0);
	failed |= check_stats("params aggregate same d_scale stats", &list, 4, 2, 8, 0);
	failed |= check_aggregation_stats("params aggregate same d_scale merge stats", &list, 8, 4);
	failed |= check_edge_ex("params aggregate same d_scale edge", &list.edges[1], 6, 9, 0.6f, 0.5f, HK_BLIND_STATE_01, HK_BLIND_STATE_MASK(HK_BLIND_STATE_01));
	hk_blind_wedge_list_destroy(&list);

	set_bpair_params(&bpairs[1], 1.0f, 0.25f);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("params aggregate different d_scale build ret", hk_blind_wedge_list_build_from_bpair_set_params(&list, &set, 1.0f), 0);
	failed |= check_stats("params aggregate different d_scale build stats", &list, 8, 2, 8, 0);
	failed |= check_i32("params aggregate different d_scale ret", hk_blind_wedge_list_aggregate_exact(&list), 0);
	failed |= check_stats("params aggregate different d_scale stats", &list, 8, 2, 8, 0);
	failed |= check_aggregation_stats("params aggregate different d_scale merge stats", &list, 8, 0);

	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_zero_initialized_output(void)
{
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list = {0};
	int failed = 0;

	set_bpair(&bpairs[0], 2, 5, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair_set(&set, bpairs, 1);
	failed |= check_i32("zero-init build ret", hk_blind_wedge_list_build_from_bpair_set(&list, &set, 2.0f, 0.7f, 1.0f), 0);
	failed |= check_stats("zero-init stats", &list, 4, 1, 4, 0);
	failed |= check_aggregation_stats("zero-init aggregation stats", &list, 4, 0);
	failed |= check_edge("zero-init 00", &list.edges[0], 4, 10, 0.5f, 0.7f, HK_BLIND_STATE_00);
	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_deterministic_ordering_and_reuse(void)
{
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	int failed = 0;

	set_bpair(&bpairs[0], 2, 5, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair(&bpairs[1], 3, 3, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair_set(&set, bpairs, 1);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("reuse first build ret", hk_blind_wedge_list_build_from_bpair_set(&list, &set, 2.0f, 0.7f, 1.0f), 0);
	failed |= check_stats("reuse first stats", &list, 4, 1, 4, 0);
	failed |= check_aggregation_stats("reuse first aggregation stats", &list, 4, 0);

	set_bpair_set(&set, bpairs, 2);
	failed |= check_i32("reuse second build ret", hk_blind_wedge_list_build_from_bpair_set(&list, &set, 2.0f, 0.7f, 1.0f), 0);
	failed |= check_stats("reuse second stats", &list, 4, 2, 4, 0);
	failed |= check_i64("reuse same-bin skipped", list.n_skipped_same_bin_bpairs, 1);
	failed |= check_aggregation_stats("reuse second aggregation stats", &list, 4, 0);
	failed |= check_edge("order bpair0 00", &list.edges[0], 4, 10, 0.5f, 0.7f, HK_BLIND_STATE_00);
	failed |= check_edge("order bpair0 01", &list.edges[1], 4, 11, 0.5f, 0.7f, HK_BLIND_STATE_01);
	failed |= check_edge("order bpair0 10", &list.edges[2], 5, 10, 0.5f, 0.7f, HK_BLIND_STATE_10);
	failed |= check_edge("order bpair0 11", &list.edges[3], 5, 11, 0.5f, 0.7f, HK_BLIND_STATE_11);

	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_same_bin_aggregation(void)
{
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	float sum_before, sum_after;
	int failed = 0;

	set_bpair(&bpairs[0], 3, 3, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair_set(&set, bpairs, 1);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("same-bin aggregate build ret", hk_blind_wedge_list_build_from_bpair_set(&list, &set, 2.0f, 0.7f, 1.0f), 0);
	failed |= check_stats("same-bin aggregate build stats", &list, 0, 1, 0, 0);
	failed |= check_i64("same-bin aggregate skipped", list.n_skipped_same_bin_bpairs, 1);
	sum_before = sum_edge_k(&list);

	failed |= check_i32("same-bin aggregate ret", hk_blind_wedge_list_aggregate_exact(&list), 0);
	sum_after = sum_edge_k(&list);
	failed |= check_stats("same-bin aggregate stats", &list, 0, 1, 0, 0);
	failed |= check_aggregation_stats("same-bin aggregate merge stats", &list, 0, 0);
	failed |= check_close("same-bin aggregate conserved k", sum_after, sum_before);

	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_non_diagonal_no_aggregation(void)
{
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list list;
	float sum_before, sum_after;
	int failed = 0;

	set_bpair(&bpairs[0], 2, 5, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair_set(&set, bpairs, 1);
	hk_blind_wedge_list_init(&list);
	failed |= check_i32("non-diagonal aggregate build ret", hk_blind_wedge_list_build_from_bpair_set(&list, &set, 2.0f, 0.7f, 1.0f), 0);
	sum_before = sum_edge_k(&list);

	failed |= check_i32("non-diagonal aggregate ret", hk_blind_wedge_list_aggregate_exact(&list), 0);
	sum_after = sum_edge_k(&list);
	failed |= check_stats("non-diagonal aggregate stats", &list, 4, 1, 4, 0);
	failed |= check_aggregation_stats("non-diagonal aggregate merge stats", &list, 4, 0);
	failed |= check_edge("non-diagonal aggregate 00", &list.edges[0], 4, 10, 0.5f, 0.7f, HK_BLIND_STATE_00);
	failed |= check_edge("non-diagonal aggregate 01", &list.edges[1], 4, 11, 0.5f, 0.7f, HK_BLIND_STATE_01);
	failed |= check_edge("non-diagonal aggregate 10", &list.edges[2], 5, 10, 0.5f, 0.7f, HK_BLIND_STATE_10);
	failed |= check_edge("non-diagonal aggregate 11", &list.edges[3], 5, 11, 0.5f, 0.7f, HK_BLIND_STATE_11);
	failed |= check_close("non-diagonal aggregate conserved k", sum_after, sum_before);

	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_duplicate_different_d_scale(void)
{
	struct hk_blind_wedge edges[2];
	struct hk_blind_wedge_list list = {0};
	float sum_before, sum_after;
	int failed = 0;

	set_wedge(&edges[0], 6, 7, 0.5f, 0.8f, HK_BLIND_STATE_10);
	set_wedge(&edges[1], 6, 7, 0.5f, 0.7f, HK_BLIND_STATE_01);
	list.edges = edges;
	list.n_edges = 2;
	list.m_edges = 2;
	sum_before = sum_edge_k(&list);

	failed |= check_i32("different d_scale aggregate ret", hk_blind_wedge_list_aggregate_exact(&list), 0);
	sum_after = sum_edge_k(&list);
	failed |= check_i32("different d_scale n_edges", list.n_edges, 2);
	failed |= check_aggregation_stats("different d_scale aggregate stats", &list, 2, 0);
	failed |= check_edge("different d_scale first", &list.edges[0], 6, 7, 0.5f, 0.7f, HK_BLIND_STATE_01);
	failed |= check_edge("different d_scale second", &list.edges[1], 6, 7, 0.5f, 0.8f, HK_BLIND_STATE_10);
	failed |= check_close("different d_scale conserved k", sum_after, sum_before);
	return failed;
}

static int check_duplicate_same_d_scale(void)
{
	struct hk_blind_wedge edges[3];
	struct hk_blind_wedge_list list = {0};
	float sum_before, sum_after;
	uint8_t mask = HK_BLIND_STATE_MASK(HK_BLIND_STATE_01) | HK_BLIND_STATE_MASK(HK_BLIND_STATE_10) |
				   HK_BLIND_STATE_MASK(HK_BLIND_STATE_11);
	int failed = 0;

	set_wedge(&edges[0], 6, 7, 0.4f, 0.7f, HK_BLIND_STATE_11);
	set_wedge(&edges[1], 6, 7, 0.2f, 0.7f, HK_BLIND_STATE_01);
	set_wedge(&edges[2], 6, 7, 0.3f, 0.7f, HK_BLIND_STATE_10);
	list.edges = edges;
	list.n_edges = 3;
	list.m_edges = 3;
	sum_before = sum_edge_k(&list);

	failed |= check_i32("same d_scale aggregate ret", hk_blind_wedge_list_aggregate_exact(&list), 0);
	sum_after = sum_edge_k(&list);
	failed |= check_i32("same d_scale n_edges", list.n_edges, 1);
	failed |= check_aggregation_stats("same d_scale aggregate stats", &list, 3, 2);
	failed |= check_edge_ex("same d_scale merged edge", &list.edges[0], 6, 7, 0.9f, 0.7f, HK_BLIND_STATE_01, mask);
	failed |= check_close("same d_scale conserved k", sum_after, sum_before);
	return failed;
}

static int check_zero_k_conservation(void)
{
	struct hk_blind_wedge edges[3];
	struct hk_blind_wedge_list list = {0};
	float sum_before, sum_after;
	uint8_t mask = HK_BLIND_STATE_MASK(HK_BLIND_STATE_01) | HK_BLIND_STATE_MASK(HK_BLIND_STATE_10);
	int failed = 0;

	set_wedge(&edges[0], 8, 9, 0.0f, 0.7f, HK_BLIND_STATE_00);
	set_wedge(&edges[1], 6, 7, 0.2f, 0.7f, HK_BLIND_STATE_01);
	set_wedge(&edges[2], 6, 7, 0.0f, 0.7f, HK_BLIND_STATE_10);
	list.edges = edges;
	list.n_edges = 3;
	list.m_edges = 3;
	sum_before = sum_edge_k(&list);

	failed |= check_i32("zero-k aggregate ret", hk_blind_wedge_list_aggregate_exact(&list), 0);
	sum_after = sum_edge_k(&list);
	failed |= check_i32("zero-k n_edges", list.n_edges, 2);
	failed |= check_aggregation_stats("zero-k aggregate stats", &list, 3, 1);
	failed |= check_edge_ex("zero-k merged edge", &list.edges[0], 6, 7, 0.2f, 0.7f, HK_BLIND_STATE_01, mask);
	failed |= check_edge("zero-k nonduplicate edge", &list.edges[1], 8, 9, 0.0f, 0.7f, HK_BLIND_STATE_00);
	failed |= check_close("zero-k conserved k", sum_after, sum_before);
	return failed;
}

static int check_aggregation_deterministic_ordering(void)
{
	struct hk_blind_wedge edges[4];
	struct hk_blind_wedge_list list = {0};
	int failed = 0;

	set_wedge(&edges[0], 10, 12, 0.1f, 0.9f, HK_BLIND_STATE_00);
	set_wedge(&edges[1], 6, 7, 0.2f, 0.8f, HK_BLIND_STATE_00);
	set_wedge(&edges[2], 6, 7, 0.3f, 0.7f, HK_BLIND_STATE_01);
	set_wedge(&edges[3], 2, 3, 0.4f, 0.7f, HK_BLIND_STATE_11);
	list.edges = edges;
	list.n_edges = 4;
	list.m_edges = 4;

	failed |= check_i32("ordering aggregate ret", hk_blind_wedge_list_aggregate_exact(&list), 0);
	failed |= check_i32("ordering n_edges", list.n_edges, 4);
	failed |= check_aggregation_stats("ordering aggregate stats", &list, 4, 0);
	failed |= check_edge("ordering first", &list.edges[0], 2, 3, 0.4f, 0.7f, HK_BLIND_STATE_11);
	failed |= check_edge("ordering second", &list.edges[1], 6, 7, 0.3f, 0.7f, HK_BLIND_STATE_01);
	failed |= check_edge("ordering third", &list.edges[2], 6, 7, 0.2f, 0.8f, HK_BLIND_STATE_00);
	failed |= check_edge("ordering fourth", &list.edges[3], 10, 12, 0.1f, 0.9f, HK_BLIND_STATE_00);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_non_diagonal_uniform();
	failed |= check_non_diagonal_nonuniform();
	failed |= check_same_bin_self_edges();
	failed |= check_zero_rho_keeps_edges();
	failed |= check_param_aware_builder_uses_bpair_params();
	failed |= check_param_aware_stiffness_conservation();
	failed |= check_param_entropy_training_mode();
	failed |= check_param_expected_count_d_scale();
	failed |= check_global_builder_ignores_bpair_params();
	failed |= check_param_aware_exact_aggregation();
	failed |= check_zero_initialized_output();
	failed |= check_deterministic_ordering_and_reuse();
	failed |= check_same_bin_aggregation();
	failed |= check_non_diagonal_no_aggregation();
	failed |= check_duplicate_different_d_scale();
	failed |= check_duplicate_same_d_scale();
	failed |= check_zero_k_conservation();
	failed |= check_aggregation_deterministic_ordering();
	return failed != 0;
}
