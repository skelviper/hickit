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

static void copy_coords(fvec3_t *dst, const fvec3_t *src, int32_t n)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			dst[i][a] = src[i][a];
}

static float dist3(const fvec3_t x, const fvec3_t y)
{
	float dx = x[0] - y[0];
	float dy = x[1] - y[1];
	float dz = x[2] - y[2];
	return sqrtf(dx * dx + dy * dy + dz * dz);
}

static int check_coords_close(const char *label, const fvec3_t *got, const fvec3_t *expected, int32_t n)
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

static int check_coords_finite(const fvec3_t *coords, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("finite coord", isfinite(coords[i][a]));
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

static void set_bmap(struct hk_bmap *bmap, struct hk_bead *beads, int32_t n_beads)
{
	int32_t i;
	memset(bmap, 0, sizeof(*bmap));
	bmap->n_beads = n_beads;
	bmap->beads = beads;
	for (i = 0; i < n_beads; ++i) {
		beads[i].chr = 0;
		beads[i].st = i * 1000000;
		beads[i].en = beads[i].st + 1000000;
	}
}

static int check_clean_relax_diag(const char *label, const struct hk_blind_relax_diag *diag, int32_t expected_steps)
{
	int failed = 0;
	char buf[128];
	snprintf(buf, sizeof(buf), "%s n_steps", label);
	failed |= check_i32(buf, diag->n_steps, expected_steps);
	snprintf(buf, sizeof(buf), "%s n_completed", label);
	failed |= check_i32(buf, diag->n_completed, expected_steps);
	snprintf(buf, sizeof(buf), "%s nonfinite step", label);
	failed |= check_i32(buf, diag->n_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s coord nonfinite", label);
	failed |= check_i32(buf, diag->n_coord_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s initial total finite", label);
	failed |= check_true(buf, isfinite(diag->initial_total_energy));
	snprintf(buf, sizeof(buf), "%s final total finite", label);
	failed |= check_true(buf, isfinite(diag->final_total_energy));
	snprintf(buf, sizeof(buf), "%s initial backbone finite", label);
	failed |= check_true(buf, isfinite(diag->initial_backbone_energy));
	snprintf(buf, sizeof(buf), "%s final backbone finite", label);
	failed |= check_true(buf, isfinite(diag->final_backbone_energy));
	snprintf(buf, sizeof(buf), "%s initial repulsion finite", label);
	failed |= check_true(buf, isfinite(diag->initial_repulsion_energy));
	snprintf(buf, sizeof(buf), "%s final repulsion finite", label);
	failed |= check_true(buf, isfinite(diag->final_repulsion_energy));
	snprintf(buf, sizeof(buf), "%s max force finite", label);
	failed |= check_true(buf, isfinite(diag->max_force_l1));
	snprintf(buf, sizeof(buf), "%s final force finite", label);
	failed |= check_true(buf, isfinite(diag->final_force_l1));
	snprintf(buf, sizeof(buf), "%s max backbone force finite", label);
	failed |= check_true(buf, isfinite(diag->max_backbone_force_l1));
	snprintf(buf, sizeof(buf), "%s final backbone force finite", label);
	failed |= check_true(buf, isfinite(diag->final_backbone_force_l1));
	snprintf(buf, sizeof(buf), "%s max repulsion force finite", label);
	failed |= check_true(buf, isfinite(diag->max_repulsion_force_l1));
	snprintf(buf, sizeof(buf), "%s final repulsion force finite", label);
	failed |= check_true(buf, isfinite(diag->final_repulsion_force_l1));
	snprintf(buf, sizeof(buf), "%s initial total sum", label);
	failed |= check_close(buf, diag->initial_total_energy,
						  diag->initial_contact_energy + diag->initial_backbone_energy +
						  diag->initial_repulsion_energy + diag->initial_sep_energy);
	snprintf(buf, sizeof(buf), "%s final total sum", label);
	failed |= check_close(buf, diag->final_total_energy,
						  diag->final_contact_energy + diag->final_backbone_energy +
						  diag->final_repulsion_energy + diag->final_sep_energy);
	snprintf(buf, sizeof(buf), "%s backbone nonfinite step", label);
	failed |= check_i32(buf, diag->n_backbone_nonfinite_step, 0);
	snprintf(buf, sizeof(buf), "%s repulsion nonfinite step", label);
	failed |= check_i32(buf, diag->n_repulsion_nonfinite_step, 0);
	return failed;
}

static int check_zero_steps(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_relax_diag diag;
	fvec3_t coords[2], before[2];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 2.5f, 0.0f, 0.0f);
	copy_coords(before, coords, 2);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	failed |= check_i32("zero-steps ret", hk_blind_relax_cpu(&conf, &list, 0, 1, coords, 1.0f, 0.01f, 0,
															 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, &diag), 0);
	failed |= check_coords_close("zero-steps coords", coords, before, 2);
	failed |= check_clean_relax_diag("zero-steps diag", &diag, 0);
	failed |= check_close("zero-steps initial energy", diag.initial_total_energy, 0.0f);
	failed |= check_close("zero-steps final energy", diag.final_total_energy, 0.0f);
	return failed;
}

static int check_long_distance_relaxation(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_relax_diag diag;
	fvec3_t coords[2];
	float before, after;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 2.5f, 0.0f, 0.0f);
	before = dist3(coords[0], coords[1]);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	failed |= check_i32("long-relax ret", hk_blind_relax_cpu(&conf, &list, 0, 1, coords, 1.0f, 0.01f, 6,
															0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, &diag), 0);
	after = dist3(coords[0], coords[1]);
	failed |= check_true("long-relax distance decreases", after < before);
	failed |= check_true("long-relax initial contact", diag.initial_contact_energy > 0.0f);
	failed |= check_close("long-relax initial backbone", diag.initial_backbone_energy, 0.0f);
	failed |= check_clean_relax_diag("long-relax diag", &diag, 6);
	failed |= check_coords_finite(coords, 2);
	return failed;
}

static int check_short_distance_relaxation(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_relax_diag diag;
	fvec3_t coords[2];
	float before, after;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.25f, 0.0f, 0.0f);
	before = dist3(coords[0], coords[1]);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	failed |= check_i32("short-relax ret", hk_blind_relax_cpu(&conf, &list, 0, 1, coords, 1.0f, 0.01f, 6,
															 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, &diag), 0);
	after = dist3(coords[0], coords[1]);
	failed |= check_true("short-relax distance increases", after > before);
	failed |= check_true("short-relax initial contact", diag.initial_contact_energy > 0.0f);
	failed |= check_clean_relax_diag("short-relax diag", &diag, 6);
	failed |= check_coords_finite(coords, 2);
	return failed;
}

static int check_anti_collapse_relaxation(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list list;
	struct hk_blind_relax_diag diag;
	fvec3_t coords[2];
	float before, after;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.25f, 0.0f, 0.0f);
	before = dist3(coords[0], coords[1]);
	set_list(&list, 0, 0);

	failed |= check_i32("anti-relax ret", hk_blind_relax_cpu(&conf, &list, 0, 1, coords, 1.0f, 1.0f, 5,
															0.5f, 0.05f, 0, HK_BLIND_REPULSION_NONE, &diag), 0);
	after = dist3(coords[0], coords[1]);
	failed |= check_true("anti-relax separation increases", after > before);
	failed |= check_close("anti-relax initial contact", diag.initial_contact_energy, 0.0f);
	failed |= check_true("anti-relax initial sep", diag.initial_sep_energy > 0.0f);
	failed |= check_clean_relax_diag("anti-relax diag", &diag, 5);
	failed |= check_coords_finite(coords, 2);
	return failed;
}

static int check_contact_and_anti_collapse_relaxation(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_relax_diag diag;
	fvec3_t coords[4];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.25f, 0.0f, 0.0f);
	set_coord(coords[2], 2.0f, 0.0f, 0.0f);
	set_coord(coords[3], 4.5f, 0.0f, 0.0f);
	set_edge(&edge, 2, 3, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	failed |= check_i32("combined-relax ret", hk_blind_relax_cpu(&conf, &list, 0, 2, coords, 1.0f, 0.01f, 5,
																0.5f, 0.05f, 0, HK_BLIND_REPULSION_NONE, &diag), 0);
	failed |= check_true("combined-relax initial contact", diag.initial_contact_energy > 0.0f);
	failed |= check_true("combined-relax initial sep", diag.initial_sep_energy > 0.0f);
	failed |= check_true("combined-relax final force", diag.final_force_l1 > 0.0f);
	failed |= check_clean_relax_diag("combined-relax diag", &diag, 5);
	failed |= check_coords_finite(coords, 4);
	return failed;
}

static int check_step_zero_relaxation(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_relax_diag diag;
	fvec3_t coords[2], before[2];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 2.5f, 0.0f, 0.0f);
	copy_coords(before, coords, 2);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	failed |= check_i32("step-zero-relax ret", hk_blind_relax_cpu(&conf, &list, 0, 1, coords, 1.0f, 0.0f, 4,
																 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, &diag), 0);
	failed |= check_coords_close("step-zero-relax coords", coords, before, 2);
	failed |= check_true("step-zero-relax initial energy", diag.initial_total_energy > 0.0f);
	failed |= check_true("step-zero-relax final energy", diag.final_total_energy > 0.0f);
	failed |= check_true("step-zero-relax force", diag.final_force_l1 > 0.0f);
	failed |= check_clean_relax_diag("step-zero-relax diag", &diag, 4);
	return failed;
}

static int check_backbone_enabled_relaxation(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list list;
	struct hk_blind_relax_diag diag;
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	fvec3_t coords[4];
	float before, after;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_bmap(&bmap, beads, 2);
	set_list(&list, 0, 0);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 2.0f, 0.0f, 0.0f);
	set_coord(coords[3], 2.0f, 1.0f, 0.0f);
	before = dist3(coords[0], coords[2]);
	failed |= check_i32("backbone-relax ret", hk_blind_relax_cpu(&conf, &list, &bmap, bmap.n_beads, coords,
																 1.0f, 0.01f, 5, 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, &diag), 0);
	after = dist3(coords[0], coords[2]);
	failed |= check_true("backbone-relax distance decreases", after < before);
	failed |= check_true("backbone-relax initial energy", diag.initial_backbone_energy > 0.0f);
	failed |= check_true("backbone-relax final energy finite", isfinite(diag.final_backbone_energy));
	failed |= check_true("backbone-relax force finite", isfinite(diag.final_backbone_force_l1));
	failed |= check_true("backbone-relax force active", diag.max_backbone_force_l1 > 0.0f);
	failed |= check_clean_relax_diag("backbone-relax diag", &diag, 5);
	failed |= check_coords_finite(coords, 4);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_zero_steps();
	failed |= check_long_distance_relaxation();
	failed |= check_short_distance_relaxation();
	failed |= check_anti_collapse_relaxation();
	failed |= check_contact_and_anti_collapse_relaxation();
	failed |= check_step_zero_relaxation();
	failed |= check_backbone_enabled_relaxation();
	return failed != 0;
}
