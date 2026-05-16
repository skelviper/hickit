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

static float dist3(const fvec3_t x, const fvec3_t y)
{
	float dx = x[0] - y[0];
	float dy = x[1] - y[1];
	float dz = x[2] - y[2];
	return sqrtf(dx * dx + dy * dy + dz * dz);
}

static void set_test_bmap(struct hk_bmap *bmap, struct hk_bead *beads, int32_t n_haploid)
{
	int32_t i;
	memset(bmap, 0, sizeof(*bmap));
	bmap->n_beads = n_haploid;
	bmap->unit = 1.0f;
	bmap->beads = beads;
	for (i = 0; i < n_haploid; ++i) {
		beads[i].chr = i < 2? 0 : 1;
		beads[i].st = i * 1000000;
		beads[i].en = beads[i].st + 1000000;
	}
}

static void set_bpair(struct hk_blind_bpair *bp, int32_t bid0, int32_t bid1,
					  float p00, float p01, float p10, float p11)
{
	bp->key.bid[0] = bid0;
	bp->key.bid[1] = bid1;
	bp->n_raw = 1;
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

static void set_bpair_set(struct hk_blind_bpair_set *set, struct hk_blind_bpair *bpairs, int32_t n_bpairs)
{
	set->bpairs = bpairs;
	set->n_bpairs = n_bpairs;
	set->raw2binned = 0;
	set->n_raw = 0;
}

static float sum_edge_k(const struct hk_blind_wedge_list *list)
{
	float sum = 0.0f;
	int32_t i;
	for (i = 0; i < list->n_edges; ++i)
		sum += list->edges[i].k;
	return sum;
}

static int check_coords_finite(const fvec3_t *coords, int32_t n_diploid)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n_diploid; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("finite coord", isfinite(coords[i][a]));
	return failed;
}

static int check_wedges_finite(const struct hk_blind_wedge_list *list)
{
	int failed = 0;
	int32_t i;
	for (i = 0; i < list->n_edges; ++i) {
		failed |= check_true("edge finite k", isfinite(list->edges[i].k));
		failed |= check_true("edge finite d_scale", isfinite(list->edges[i].d_scale));
		failed |= check_true("edge nonnegative k", list->edges[i].k >= 0.0f);
		failed |= check_true("edge positive d_scale", list->edges[i].d_scale > 0.0f);
		failed |= check_true("edge not self", list->edges[i].bid[0] != list->edges[i].bid[1]);
	}
	failed |= check_true("edge sum finite", isfinite(sum_edge_k(list)));
	return failed;
}

static int check_relax_ok(const char *label, const struct hk_blind_relax_diag *diag, int32_t n_steps)
{
	int failed = 0;
	char buf[128];
	snprintf(buf, sizeof(buf), "%s n_completed", label);
	failed |= check_i32(buf, diag->n_completed, n_steps);
	snprintf(buf, sizeof(buf), "%s nonfinite step", label);
	failed |= check_i32(buf, diag->n_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s coord nonfinite", label);
	failed |= check_i32(buf, diag->n_coord_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s initial energy finite", label);
	failed |= check_true(buf, isfinite(diag->initial_total_energy));
	snprintf(buf, sizeof(buf), "%s final energy finite", label);
	failed |= check_true(buf, isfinite(diag->final_total_energy));
	snprintf(buf, sizeof(buf), "%s max force finite", label);
	failed |= check_true(buf, isfinite(diag->max_force_l1));
	snprintf(buf, sizeof(buf), "%s final force finite", label);
	failed |= check_true(buf, isfinite(diag->final_force_l1));
	snprintf(buf, sizeof(buf), "%s backbone initial finite", label);
	failed |= check_true(buf, isfinite(diag->initial_backbone_energy));
	snprintf(buf, sizeof(buf), "%s backbone final finite", label);
	failed |= check_true(buf, isfinite(diag->final_backbone_energy));
	snprintf(buf, sizeof(buf), "%s repulsion initial finite", label);
	failed |= check_true(buf, isfinite(diag->initial_repulsion_energy));
	snprintf(buf, sizeof(buf), "%s repulsion final finite", label);
	failed |= check_true(buf, isfinite(diag->final_repulsion_energy));
	snprintf(buf, sizeof(buf), "%s backbone force finite", label);
	failed |= check_true(buf, isfinite(diag->final_backbone_force_l1));
	snprintf(buf, sizeof(buf), "%s repulsion force finite", label);
	failed |= check_true(buf, isfinite(diag->final_repulsion_force_l1));
	snprintf(buf, sizeof(buf), "%s initial total sum", label);
	failed |= check_close(buf, diag->initial_total_energy,
						  diag->initial_contact_energy + diag->initial_backbone_energy +
						  diag->initial_repulsion_energy + diag->initial_sep_energy);
	snprintf(buf, sizeof(buf), "%s final total sum", label);
	failed |= check_close(buf, diag->final_total_energy,
						  diag->final_contact_energy + diag->final_backbone_energy +
						  diag->final_repulsion_energy + diag->final_sep_energy);
	snprintf(buf, sizeof(buf), "%s backbone nonfinite", label);
	failed |= check_i32(buf, diag->n_backbone_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s repulsion nonfinite", label);
	failed |= check_i32(buf, diag->n_repulsion_nonfinite_step, 0);
	return failed;
}

static int build_aggregated_edges(const struct hk_blind_bpair_set *set, float base_k, float base_d_scale, float rho_train,
								  struct hk_blind_wedge_list *list)
{
	int ret;
	hk_blind_wedge_list_init(list);
	ret = hk_blind_wedge_list_build_from_bpair_set(list, set, base_k, base_d_scale, rho_train);
	if (ret != 0) return ret;
	return hk_blind_wedge_list_aggregate_exact(list);
}

static int check_long_distance_weighted_contact(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list edges;
	struct hk_fdg_conf conf;
	struct hk_blind_relax_diag diag;
	fvec3_t coords[4];
	float before, after;
	int failed = 0;

	set_test_bmap(&bmap, beads, 2);
	set_bpair(&bpairs[0], 0, 1, 1.0f, 0.0f, 0.0f, 0.0f);
	set_bpair_set(&set, bpairs, 1);
	failed |= check_i32("long build", build_aggregated_edges(&set, 2.0f, 1.0f, 1.0f, &edges), 0);
	failed |= check_i32("long n_edges", edges.n_edges, 4);
	failed |= check_close("long sum k", sum_edge_k(&edges), 2.0f);
	failed |= check_wedges_finite(&edges);

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 100.0f, 0.0f, 0.0f);
	set_coord(coords[2], 2.5f, 0.0f, 0.0f);
	set_coord(coords[3], 100.0f, 1.0f, 0.0f);
	before = dist3(coords[hk_diploid_bid(0, HK_DIPLOID_COPY0)], coords[hk_diploid_bid(1, HK_DIPLOID_COPY0)]);
	failed |= check_i32("long relax", hk_blind_relax_cpu(&conf, &edges, 0, bmap.n_beads, coords, 1.0f, 0.01f, 6,
														 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, &diag), 0);
	after = dist3(coords[hk_diploid_bid(0, HK_DIPLOID_COPY0)], coords[hk_diploid_bid(1, HK_DIPLOID_COPY0)]);
	failed |= check_true("long 00 distance decreases", after < before);
	failed |= check_relax_ok("long diag", &diag, 6);
	failed |= check_coords_finite(coords, 2 * bmap.n_beads);
	hk_blind_wedge_list_destroy(&edges);
	return failed;
}

static int check_anti_collapse_only_path(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list edges;
	struct hk_fdg_conf conf;
	struct hk_blind_relax_diag diag;
	fvec3_t coords[2];
	float before, after;
	int failed = 0;

	set_test_bmap(&bmap, beads, 1);
	set_bpair_set(&set, 0, 0);
	failed |= check_i32("anti build", build_aggregated_edges(&set, 2.0f, 1.0f, 1.0f, &edges), 0);
	failed |= check_i32("anti n_edges", edges.n_edges, 0);
	failed |= check_close("anti sum k", sum_edge_k(&edges), 0.0f);

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.25f, 0.0f, 0.0f);
	before = dist3(coords[0], coords[1]);
	failed |= check_i32("anti relax", hk_blind_relax_cpu(&conf, &edges, 0, bmap.n_beads, coords, 1.0f, 1.0f, 5,
														 0.5f, 0.05f, 0, HK_BLIND_REPULSION_NONE, &diag), 0);
	after = dist3(coords[0], coords[1]);
	failed |= check_true("anti separation increases", after > before);
	failed |= check_true("anti sep initially active", diag.initial_sep_energy > 0.0f);
	failed |= check_relax_ok("anti diag", &diag, 5);
	failed |= check_coords_finite(coords, 2 * bmap.n_beads);
	hk_blind_wedge_list_destroy(&edges);
	return failed;
}

static float run_weighted_00_distance_decrease(float p00, float p01, int *failed, const char *label)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list edges;
	struct hk_fdg_conf conf;
	struct hk_blind_relax_diag diag;
	fvec3_t coords[4];
	float before, after, sum_k_before, sum_k_after;
	int32_t n_edges_before;

	set_test_bmap(&bmap, beads, 2);
	set_bpair(&bpairs[0], 0, 1, p00, p01, 0.0f, 0.0f);
	set_bpair_set(&set, bpairs, 1);
	if (build_aggregated_edges(&set, 2.0f, 1.0f, 1.0f, &edges) != 0)
		return -1.0f;
	sum_k_before = sum_edge_k(&edges);
	n_edges_before = edges.n_edges;
	*failed |= check_wedges_finite(&edges);

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 100.0f, 0.0f, 0.0f);
	set_coord(coords[2], 2.5f, 0.0f, 0.0f);
	set_coord(coords[3], 1.0f, 0.0f, 0.0f); // The nonzero 01 edge is in Hickit's zero-force plateau.
	before = dist3(coords[0], coords[2]);
	if (hk_blind_relax_cpu(&conf, &edges, 0, bmap.n_beads, coords, 1.0f, 0.01f, 1, 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, &diag) != 0) {
		hk_blind_wedge_list_destroy(&edges);
		return -1.0f;
	}
	after = dist3(coords[0], coords[2]);
	sum_k_after = sum_edge_k(&edges);
	*failed |= check_relax_ok(label, &diag, 1);
	*failed |= check_coords_finite(coords, 2 * bmap.n_beads);
	*failed |= check_i32("stiffness edge count unchanged", edges.n_edges, n_edges_before);
	*failed |= check_close("stiffness sum k unchanged", sum_k_after, sum_k_before);
	hk_blind_wedge_list_destroy(&edges);
	return before - after;
}

static int check_posterior_weight_controls_stiffness(void)
{
	int failed = 0;
	float full_decrease = run_weighted_00_distance_decrease(1.0f, 0.0f, &failed, "full-weight diag");
	float half_decrease = run_weighted_00_distance_decrease(0.5f, 0.5f, &failed, "half-weight diag");

	failed |= check_true("full decrease finite", isfinite(full_decrease));
	failed |= check_true("half decrease finite", isfinite(half_decrease));
	failed |= check_true("full decrease positive", full_decrease > 0.0f);
	failed |= check_true("half decrease positive", half_decrease > 0.0f);
	failed |= check_true("half-weight movement weaker", half_decrease < full_decrease);
	return failed;
}

static int check_same_bin_self_edge_filtering(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	struct hk_blind_wedge_list edges;
	struct hk_fdg_conf conf;
	struct hk_blind_relax_diag diag;
	fvec3_t coords[4];
	int failed = 0;

	set_test_bmap(&bmap, beads, 2);
	set_bpair(&bpairs[0], 1, 1, 0.25f, 0.25f, 0.25f, 0.25f);
	set_bpair_set(&set, bpairs, 1);
	failed |= check_i32("same-bin build", build_aggregated_edges(&set, 2.0f, 1.0f, 1.0f, &edges), 0);
	failed |= check_i64("same-bin skipped self", edges.n_skipped_self_edges, 0);
	failed |= check_i64("same-bin skipped bpair", edges.n_skipped_same_bin_bpairs, 1);
	failed |= check_i32("same-bin aggregate edge count", edges.n_edges, 0);
	failed |= check_close("same-bin sum k", sum_edge_k(&edges), 0.0f);
	failed |= check_wedges_finite(&edges);

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 10.0f, 0.0f, 0.0f);
	set_coord(coords[1], 11.0f, 0.0f, 0.0f);
	set_coord(coords[2], 0.0f, 0.0f, 0.0f);
	set_coord(coords[3], 1.0f, 0.0f, 0.0f);
	failed |= check_i32("same-bin relax", hk_blind_relax_cpu(&conf, &edges, 0, bmap.n_beads, coords, 1.0f, 0.01f, 3,
															 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, &diag), 0);
	failed |= check_relax_ok("same-bin diag", &diag, 3);
	failed |= check_coords_finite(coords, 2 * bmap.n_beads);
	hk_blind_wedge_list_destroy(&edges);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_long_distance_weighted_contact();
	failed |= check_anti_collapse_only_path();
	failed |= check_posterior_weight_controls_stiffness();
	failed |= check_same_bin_self_edge_filtering();
	return failed != 0;
}
