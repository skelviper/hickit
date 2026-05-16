#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "hkpriv.h"

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

static void zero_force(fvec3_t *force, int32_t n)
{
	int32_t i;
	for (i = 0; i < n; ++i)
		set_coord(force[i], 0.0f, 0.0f, 0.0f);
}

static void set_sentinel_force(fvec3_t *force, int32_t n)
{
	int32_t i;
	for (i = 0; i < n; ++i)
		set_coord(force[i], 0.25f * (float)(i + 1), -0.5f * (float)(i + 1), 0.75f * (float)(i + 1));
}

static void copy_force(fvec3_t *dst, const fvec3_t *src, int32_t n)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			dst[i][a] = src[i][a];
}

static int check_force_close(const char *label, const fvec3_t *got, const fvec3_t *expected, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	char buf[96];
	for (i = 0; i < n; ++i) {
		for (a = 0; a < 3; ++a) {
			snprintf(buf, sizeof(buf), "%s bead%d axis%d", label, (int)i, a);
			failed |= check_close(buf, got[i][a], expected[i][a]);
		}
	}
	return failed;
}

static void set_edge(struct hk_blind_wedge *edge, int32_t bid0, int32_t bid1, float k, float d_scale, int32_t state)
{
	edge->bid[0] = bid0;
	edge->bid[1] = bid1;
	edge->k = k;
	edge->d_scale = d_scale;
	edge->state = state;
	edge->state_mask = HK_BLIND_STATE_MASK(state);
}

static void set_list(struct hk_blind_wedge_list *list, struct hk_blind_wedge *edges, int32_t n_edges)
{
	list->edges = edges;
	list->n_edges = n_edges;
	list->m_edges = n_edges;
	list->n_input_bpair = 0;
	list->n_expanded_edges = n_edges;
	list->n_skipped_self_edges = 0;
	list->n_edges_before_aggregation = n_edges;
	list->n_aggregated_edges_removed = 0;
}

static void set_diag_bmap(struct hk_bmap *bmap, struct hk_bead *beads)
{
	int32_t i;
	bmap->d = 0;
	bmap->n_beads = 3;
	bmap->beads = beads;
	for (i = 0; i < 3; ++i) {
		beads[i].chr = i < 2? 0 : 1;
		beads[i].st = i * 1000000;
		beads[i].en = beads[i].st + 1000000;
	}
}

static int check_zero_k_edge(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	fvec3_t coords[2], force[2], expected[2];
	int32_t n_nonfinite = -1;
	float energy;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 2.5f, 0.0f, 0.0f);
	set_sentinel_force(force, 2);
	copy_force(expected, force, 2);
	set_edge(&edge, 0, 1, 0.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	energy = hk_blind_wedge_list_accumulate_contact_force_cpu(&conf, &list, coords, 2, 1.0f, force, &n_nonfinite);
	failed |= check_close("zero-k energy", energy, 0.0f);
	failed |= check_i32("zero-k nonfinite", n_nonfinite, 0);
	failed |= check_force_close("zero-k force unchanged", force, expected, 2);
	return failed;
}

static int check_plateau_region(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	fvec3_t coords[2], force[2], expected[2];
	int32_t n_nonfinite = -1;
	float energy;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 1.0f, 0.0f, 0.0f);
	zero_force(force, 2);
	zero_force(expected, 2);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	energy = hk_blind_wedge_list_accumulate_contact_force_cpu(&conf, &list, coords, 2, 1.0f, force, &n_nonfinite);
	failed |= check_close("plateau energy", energy, 0.0f);
	failed |= check_i32("plateau nonfinite", n_nonfinite, 0);
	failed |= check_force_close("plateau force", force, expected, 2);
	return failed;
}

static int check_short_distance_force_accumulates(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	fvec3_t coords[2], force[2], initial[2], expected[2];
	int32_t n_nonfinite = -1;
	float energy;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.25f, 0.0f, 0.0f);
	set_sentinel_force(force, 2);
	copy_force(initial, force, 2);
	copy_force(expected, initial, 2);
	// Short contact side: delta points negative x; positive Hickit force pushes the beads apart.
	expected[0][0] += -1.0f;
	expected[1][0] +=  1.0f;
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	energy = hk_blind_wedge_list_accumulate_contact_force_cpu(&conf, &list, coords, 2, 1.0f, force, &n_nonfinite);
	failed |= check_close("short energy", energy, 2.0f * 0.25f * 0.25f);
	failed |= check_i32("short nonfinite", n_nonfinite, 0);
	failed |= check_force_close("short accumulated force", force, expected, 2);
	failed |= check_close("short action reaction x", (force[0][0] - initial[0][0]) + (force[1][0] - initial[1][0]), 0.0f);
	return failed;
}

static int check_long_distance_force(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	fvec3_t coords[2], force[2];
	int32_t n_nonfinite = -1;
	float energy;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 2.5f, 0.0f, 0.0f);
	zero_force(force, 2);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	energy = hk_blind_wedge_list_accumulate_contact_force_cpu(&conf, &list, coords, 2, 1.0f, force, &n_nonfinite);
	failed |= check_close("long energy consistency", energy, hk_fdg_contact_energy_dist(&conf, 2.5f, 1.0f, 1.0f, 2.0f));
	failed |= check_true("long energy positive", energy > 0.0f);
	failed |= check_i32("long nonfinite", n_nonfinite, 0);
	// Long contact side: delta points negative x; negative Hickit force pulls bead 0 toward bead 1.
	failed |= check_true("long force pulls i", force[0][0] > 0.0f);
	failed |= check_true("long force pulls j", force[1][0] < 0.0f);
	failed |= check_close("long action reaction x", force[0][0] + force[1][0], 0.0f);
	return failed;
}

static int check_long_distance_force_z_action_reaction(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	fvec3_t coords[2], force[2];
	int32_t n_nonfinite = -1;
	float energy;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 0.0f, 2.5f);
	zero_force(force, 2);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	energy = hk_blind_wedge_list_accumulate_contact_force_cpu(&conf, &list, coords, 2, 1.0f, force, &n_nonfinite);
	failed |= check_close("long z energy consistency", energy,
						  hk_fdg_contact_energy_dist(&conf, 2.5f, 1.0f, 1.0f, 2.0f));
	failed |= check_i32("long z nonfinite", n_nonfinite, 0);
	failed |= check_close("long z action reaction x", force[0][0] + force[1][0], 0.0f);
	failed |= check_close("long z action reaction y", force[0][1] + force[1][1], 0.0f);
	failed |= check_close("long z action reaction z", force[0][2] + force[1][2], 0.0f);
	failed |= check_true("long z pulls i upward", force[0][2] > 0.0f);
	failed |= check_true("long z pulls j downward", force[1][2] < 0.0f);
	failed |= check_true("long z opposite signs", force[0][2] * force[1][2] < 0.0f);
	return failed;
}

static int check_energy_consistency_table(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	float distances[] = {0.0f, 0.25f, 0.5f, 1.0f, 1.75f, 2.5f};
	int failed = 0;
	int i;

	hk_fdg_conf_init(&conf);
	set_edge(&edge, 0, 1, 1.3f, 0.7f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);
	for (i = 0; i < (int)(sizeof(distances) / sizeof(distances[0])); ++i) {
		fvec3_t coords[2], force[2];
		int32_t n_nonfinite = -1;
		float got, expected;
		char buf[96];
		set_coord(coords[0], 0.0f, 0.0f, 0.0f);
		set_coord(coords[1], distances[i], 0.0f, 0.0f);
		zero_force(force, 2);
		got = hk_blind_wedge_list_accumulate_contact_force_cpu(&conf, &list, coords, 2, 1.0f, force, &n_nonfinite);
		expected = hk_fdg_contact_energy_dist(&conf, distances[i], 1.0f, edge.d_scale, edge.k);
		snprintf(buf, sizeof(buf), "energy consistency distance%d", i);
		failed |= check_close(buf, got, expected);
		snprintf(buf, sizeof(buf), "energy consistency nonfinite%d", i);
		failed |= check_i32(buf, n_nonfinite, 0);
	}
	return failed;
}

static int check_multiple_edges_additive(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edges[2];
	struct hk_blind_wedge_list list_both, list_one;
	fvec3_t coords[3], force_both[3], force_sum[3], tmp[3];
	float e_both, e_sum = 0.0f;
	int32_t n_nonfinite = -1;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 1.75f, 0.0f, 0.0f);
	set_coord(coords[2], 2.0f, 0.0f, 0.0f);
	set_edge(&edges[0], 0, 1, 1.0f, 1.0f, HK_BLIND_STATE_00);
	set_edge(&edges[1], 1, 2, 2.0f, 1.0f, HK_BLIND_STATE_01);
	set_list(&list_both, edges, 2);

	zero_force(force_both, 3);
	e_both = hk_blind_wedge_list_accumulate_contact_force_cpu(&conf, &list_both, coords, 3, 1.0f, force_both, &n_nonfinite);
	failed |= check_i32("multi nonfinite", n_nonfinite, 0);

	zero_force(force_sum, 3);
	set_list(&list_one, &edges[0], 1);
	zero_force(tmp, 3);
	e_sum += hk_blind_wedge_list_accumulate_contact_force_cpu(&conf, &list_one, coords, 3, 1.0f, tmp, &n_nonfinite);
	for (int i = 0; i < 3; ++i)
		for (int a = 0; a < 3; ++a)
			force_sum[i][a] += tmp[i][a];
	set_list(&list_one, &edges[1], 1);
	zero_force(tmp, 3);
	e_sum += hk_blind_wedge_list_accumulate_contact_force_cpu(&conf, &list_one, coords, 3, 1.0f, tmp, &n_nonfinite);
	for (int i = 0; i < 3; ++i)
		for (int a = 0; a < 3; ++a)
			force_sum[i][a] += tmp[i][a];

	failed |= check_close("multi energy additive", e_both, e_sum);
	failed |= check_force_close("multi force additive", force_both, force_sum, 3);
	return failed;
}

static int check_exact_overlap_fallback(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	fvec3_t coords[2], force[2];
	int32_t n_nonfinite = -1;
	float energy;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 0.0f, 0.0f);
	zero_force(force, 2);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	energy = hk_blind_wedge_list_accumulate_contact_force_cpu(&conf, &list, coords, 2, 1.0f, force, &n_nonfinite);
	failed |= check_close("overlap energy", energy, 2.0f * 0.5f * 0.5f);
	failed |= check_i32("overlap nonfinite", n_nonfinite, 0);
	failed |= check_close("overlap fallback f0x", force[0][0], 2.0f);
	failed |= check_close("overlap fallback f1x", force[1][0], -2.0f);
	failed |= check_close("overlap fallback f0y", force[0][1], 0.0f);
	failed |= check_close("overlap fallback f0z", force[0][2], 0.0f);
	failed |= check_true("overlap finite energy", isfinite(energy));
	failed |= check_true("overlap finite force", isfinite(force[0][0]) && isfinite(force[1][0]));
	return failed;
}

static int check_contact_class_diag_split(void)
{
	struct hk_fdg_conf conf;
	struct hk_bmap bmap;
	struct hk_bead beads[3];
	struct hk_blind_wedge edges[2];
	struct hk_blind_wedge_list list;
	struct hk_blind_contact_class_diag diag;
	fvec3_t coords[6];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_diag_bmap(&bmap, beads);
	for (int i = 0; i < 6; ++i)
		set_coord(coords[i], 0.0f, 0.0f, 0.0f);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[2], 2.5f, 0.0f, 0.0f);
	set_coord(coords[4], 0.0f, 2.5f, 0.0f);
	set_edge(&edges[0], 0, 2, 1.0f, 1.0f, HK_BLIND_STATE_00);
	set_edge(&edges[1], 0, 4, 3.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, edges, 2);
	failed |= check_i32("class diag ret",
						hk_blind_contact_class_diag_accumulate(&conf, &bmap, &list, coords, bmap.n_beads,
															   1.0f, &diag), 0);
	failed |= check_i32("class diag cis wedges", (int32_t)diag.n_wedges_cis, 1);
	failed |= check_i32("class diag trans wedges", (int32_t)diag.n_wedges_trans, 1);
	failed |= check_i32("class diag nonfinite", diag.n_nonfinite, 0);
	failed |= check_close("class diag cis k", (float)diag.sum_wedge_k_cis, 1.0f);
	failed |= check_close("class diag trans k", (float)diag.sum_wedge_k_trans, 3.0f);
	failed |= check_true("class diag cis energy positive", diag.contact_energy_cis > 0.0f);
	failed |= check_true("class diag trans energy positive", diag.contact_energy_trans > 0.0f);
	failed |= check_true("class diag cis force positive", diag.contact_force_l1_cis > 0.0f);
	failed |= check_true("class diag trans force positive", diag.contact_force_l1_trans > 0.0f);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_zero_k_edge();
	failed |= check_plateau_region();
	failed |= check_short_distance_force_accumulates();
	failed |= check_long_distance_force();
	failed |= check_long_distance_force_z_action_reaction();
	failed |= check_energy_consistency_table();
	failed |= check_multiple_edges_additive();
	failed |= check_exact_overlap_fallback();
	failed |= check_contact_class_diag_split();
	return failed != 0;
}
