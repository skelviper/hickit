#include <math.h>
#include <stdio.h>
#include <string.h>
#include "hickit.h"

static void set_coord(fvec3_t x, float a, float b, float c)
{
	x[0] = a;
	x[1] = b;
	x[2] = c;
}

static float dist3(const fvec3_t a, const fvec3_t b)
{
	float dx = a[0] - b[0];
	float dy = a[1] - b[1];
	float dz = a[2] - b[2];
	return sqrtf(dx * dx + dy * dy + dz * dz);
}

static int check_true(const char *label, int ok)
{
	if (!ok) {
		fprintf(stderr, "FAIL: %s\n", label);
		return 1;
	}
	return 0;
}

static int check_i64(const char *label, int64_t got, int64_t expected)
{
	if (got != expected) {
		fprintf(stderr, "FAIL: %s got=%lld expected=%lld\n", label, (long long)got, (long long)expected);
		return 1;
	}
	return 0;
}

static int check_i32(const char *label, int32_t got, int32_t expected)
{
	if (got != expected) {
		fprintf(stderr, "FAIL: %s got=%d expected=%d\n", label, got, expected);
		return 1;
	}
	return 0;
}

static int check_close(const char *label, double got, double expected)
{
	double tol = 1e-5 * (fabs(expected) > 1.0? fabs(expected) : 1.0);
	if (!isfinite(got) || fabs(got - expected) > tol) {
		fprintf(stderr, "FAIL: %s got=%.9g expected=%.9g\n", label, got, expected);
		return 1;
	}
	return 0;
}

static void set_wedge(struct hk_blind_wedge *edge, int32_t bid0, int32_t bid1, float k, float d_scale)
{
	edge->bid[0] = bid0;
	edge->bid[1] = bid1;
	edge->k = k;
	edge->d_scale = d_scale;
	edge->state = HK_BLIND_STATE_00;
	edge->state_mask = 1u << HK_BLIND_STATE_00;
}

static void set_list(struct hk_blind_wedge_list *list, struct hk_blind_wedge *edges, int32_t n_edges)
{
	memset(list, 0, sizeof(*list));
	list->edges = edges;
	list->n_edges = n_edges;
	list->m_edges = n_edges;
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

static int test_single_pair_beyond_radius(void)
{
	struct hk_fdg_conf conf;
	fvec3_t x0, x1, f0, f1;
	float e;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(x0, 0.0f, 0.0f, 0.0f);
	set_coord(x1, conf.d_r + 0.5f, 0.0f, 0.0f);
	e = hk_blind_repulsion_energy_force_cpu(&conf, x0, x1, 1.0f, 1.0f, f0, f1);
	failed |= check_close("rep beyond energy", e, 0.0);
	failed |= check_close("rep beyond f0x", f0[0], 0.0);
	failed |= check_close("rep beyond f1x", f1[0], 0.0);
	return failed;
}

static int test_single_pair_inside_radius(void)
{
	struct hk_fdg_conf conf;
	fvec3_t x0, x1, f0, f1;
	float e, t, expected_e, expected_f;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(x0, 0.0f, 0.0f, 0.0f);
	set_coord(x1, 1.0f, 0.0f, 0.0f);
	t = conf.d_r - 1.0f;
	expected_e = conf.k_rel_rep * t * t;
	expected_f = 2.0f * conf.k_rel_rep * t;
	e = hk_blind_repulsion_energy_force_cpu(&conf, x0, x1, 1.0f, 1.0f, f0, f1);
	failed |= check_close("rep inside energy", e, expected_e);
	failed |= check_close("rep inside f0x", f0[0], -expected_f);
	failed |= check_close("rep inside f1x", f1[0], expected_f);
	failed |= check_close("rep inside action reaction", f0[0] + f1[0], 0.0);
	return failed;
}

static int test_single_pair_exact_overlap(void)
{
	struct hk_fdg_conf conf;
	fvec3_t x0, x1, f0, f1;
	float e, expected_f;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(x0, 0.0f, 0.0f, 0.0f);
	set_coord(x1, 0.0f, 0.0f, 0.0f);
	expected_f = 2.0f * conf.k_rel_rep * conf.d_r;
	e = hk_blind_repulsion_energy_force_cpu(&conf, x0, x1, 1.0f, 1.0f, f0, f1);
	failed |= check_true("rep overlap energy finite", isfinite(e));
	failed |= check_close("rep overlap fallback x", f0[0], expected_f);
	failed |= check_close("rep overlap fallback y", f0[1], 0.0);
	failed |= check_close("rep overlap opposite", f0[0] + f1[0], 0.0);
	return failed;
}

static int test_rel_rep_k_scaling(void)
{
	struct hk_fdg_conf conf;
	fvec3_t x0, x1, f0a, f1a, f0b, f1b;
	float e1, e05;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(x0, 0.0f, 0.0f, 0.0f);
	set_coord(x1, 1.0f, 0.0f, 0.0f);
	e1 = hk_blind_repulsion_energy_force_cpu(&conf, x0, x1, 1.0f, 1.0f, f0a, f1a);
	e05 = hk_blind_repulsion_energy_force_cpu(&conf, x0, x1, 1.0f, 0.5f, f0b, f1b);
	failed |= check_close("rep scale energy", e05, 0.5 * e1);
	failed |= check_close("rep scale force", f0b[0], 0.5 * f0a[0]);
	return failed;
}

static int test_schedule(void)
{
	int failed = 0;
	failed |= check_true("rep schedule n0", hk_blind_rel_rep_schedule_at(0, 0) == 1.0f);
	failed |= check_true("rep schedule finite", isfinite(hk_blind_rel_rep_schedule_at(0, 3)));
	failed |= check_true("rep schedule increases", hk_blind_rel_rep_schedule_at(0, 3) < hk_blind_rel_rep_schedule_at(2, 3));
	return failed;
}

static int test_batch_counts_and_accumulation(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_repulsion_diag diag;
	fvec3_t coords[4], force[4];
	float sentinel = 0.25f;
	int failed = 0;
	int i;

	hk_fdg_conf_init(&conf);
	for (i = 0; i < 4; ++i) {
		set_coord(coords[i], 10.0f * i, 0.0f, 0.0f);
		set_coord(force[i], sentinel, sentinel, sentinel);
	}
	failed |= check_close("rep batch energy",
						  hk_blind_repulsion_accumulate_force_cpu(&conf, 4, coords, 1.0f, 1.0f, 0, 0, force, &diag),
						  0.0);
	failed |= check_i64("rep batch considered", diag.n_pairs_considered, 6);
	failed |= check_i64("rep batch blocked", diag.n_pairs_blocked, 0);
	failed |= check_i64("rep batch active", diag.n_pairs_active, 0);
	failed |= check_close("rep batch sentinel", force[0][0], sentinel);
	return failed;
}

static int test_contact_blocking(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_repulsion_diag diag;
	fvec3_t coords[2], force[2];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.5f, 0.0f, 0.0f);
	set_coord(force[0], 0.0f, 0.0f, 0.0f);
	set_coord(force[1], 0.0f, 0.0f, 0.0f);
	set_wedge(&edge, 0, 1, 1.0f, 1.0f);
	set_list(&list, &edge, 1);
	failed |= check_close("rep contact block energy",
						  hk_blind_repulsion_accumulate_force_cpu(&conf, 2, coords, 1.0f, 1.0f, &list, 0, force, &diag),
						  0.0);
	failed |= check_i64("rep contact block considered", diag.n_pairs_considered, 1);
	failed |= check_i64("rep contact block blocked", diag.n_pairs_blocked, 1);
	failed |= check_i64("rep contact block active", diag.n_pairs_active, 0);
	failed |= check_close("rep contact block force", force[0][0], 0.0);
	return failed;
}

static int test_backbone_blocking(void)
{
	struct hk_fdg_conf conf;
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_repulsion_diag diag;
	fvec3_t coords[4], force[4];
	int failed = 0;
	int i;

	hk_fdg_conf_init(&conf);
	set_bmap(&bmap, beads, 2);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[2], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 10.0f, 0.0f, 0.0f);
	set_coord(coords[3], 10.0f, 0.0f, 0.0f);
	for (i = 0; i < 4; ++i)
		set_coord(force[i], 0.0f, 0.0f, 0.0f);
	failed |= check_close("rep backbone block energy",
						  hk_blind_repulsion_accumulate_force_cpu(&conf, 4, coords, 1.0f, 1.0f, 0, &bmap, force, &diag),
						  0.0);
	failed |= check_i64("rep backbone considered", diag.n_pairs_considered, 6);
	failed |= check_i64("rep backbone blocked", diag.n_pairs_blocked, 2);
	failed |= check_i64("rep backbone active", diag.n_pairs_active, 0);
	return failed;
}

static int test_homolog_pair_not_blocked(void)
{
	struct hk_fdg_conf conf;
	struct hk_bmap bmap;
	struct hk_bead bead[1];
	struct hk_blind_repulsion_diag diag;
	fvec3_t coords[2], force[2];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_bmap(&bmap, bead, 1);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.5f, 0.0f, 0.0f);
	set_coord(force[0], 0.0f, 0.0f, 0.0f);
	set_coord(force[1], 0.0f, 0.0f, 0.0f);
	failed |= check_true("rep homolog energy",
						 hk_blind_repulsion_accumulate_force_cpu(&conf, 2, coords, 1.0f, 1.0f,
																 0, &bmap, force, &diag) > 0.0f);
	failed |= check_i64("rep homolog considered", diag.n_pairs_considered, 1);
	failed |= check_i64("rep homolog blocked", diag.n_pairs_blocked, 0);
	failed |= check_i64("rep homolog active", diag.n_pairs_active, 1);
	failed |= check_true("rep homolog force", force[0][0] < 0.0f && force[1][0] > 0.0f);
	return failed;
}

static int test_net_force_zero(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_repulsion_diag diag;
	fvec3_t coords[3], force[3];
	double sx = 0.0, sy = 0.0, sz = 0.0;
	int failed = 0;
	int i;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.5f, 0.0f, 0.0f);
	set_coord(coords[2], 0.0f, 0.5f, 0.0f);
	for (i = 0; i < 3; ++i)
		set_coord(force[i], 0.0f, 0.0f, 0.0f);
	hk_blind_repulsion_accumulate_force_cpu(&conf, 3, coords, 1.0f, 1.0f, 0, 0, force, &diag);
	for (i = 0; i < 3; ++i) {
		sx += force[i][0];
		sy += force[i][1];
		sz += force[i][2];
	}
	failed |= check_true("rep net active", diag.n_pairs_active > 0);
	failed |= check_close("rep net x", sx, 0.0);
	failed |= check_close("rep net y", sy, 0.0);
	failed |= check_close("rep net z", sz, 0.0);
	return failed;
}

static int compare_n2_cell(const char *label, int32_t n, const fvec3_t coords[8],
						   const struct hk_blind_wedge_list *list,
						   const struct hk_bmap *bmap_or_null)
{
	struct hk_fdg_conf conf;
	struct hk_blind_repulsion_diag n2_diag, cell_diag;
	fvec3_t n2_force[8], cell_force[8];
	float e_n2, e_cell;
	char buf[160];
	int failed = 0;
	int32_t i;
	int a;

	hk_fdg_conf_init(&conf);
	for (i = 0; i < n; ++i) {
		set_coord(n2_force[i], 0.0f, 0.0f, 0.0f);
		set_coord(cell_force[i], 0.0f, 0.0f, 0.0f);
	}
	e_n2 = hk_blind_repulsion_accumulate_force_cpu(&conf, n, coords, 1.0f, 1.0f,
												   list, bmap_or_null, n2_force, &n2_diag);
	e_cell = hk_blind_repulsion_accumulate_force_cell_cpu(&conf, n, coords, 1.0f, 1.0f,
														  list, bmap_or_null, cell_force, &cell_diag);
	snprintf(buf, sizeof(buf), "%s energy", label);
	failed |= check_close(buf, e_cell, e_n2);
	snprintf(buf, sizeof(buf), "%s active", label);
	failed |= check_i64(buf, cell_diag.n_pairs_active, n2_diag.n_pairs_active);
	snprintf(buf, sizeof(buf), "%s nonfinite", label);
	failed |= check_i32(buf, cell_diag.n_nonfinite, n2_diag.n_nonfinite);
	for (i = 0; i < n; ++i) {
		for (a = 0; a < 3; ++a) {
			snprintf(buf, sizeof(buf), "%s force %d/%d", label, (int)i, a);
			failed |= check_close(buf, cell_force[i][a], n2_force[i][a]);
		}
	}
	return failed;
}

static int test_cell_matches_n2_simple(void)
{
	fvec3_t coords[8];

	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.5f, 0.0f, 0.0f);
	set_coord(coords[2], 0.0f, 0.5f, 0.0f);
	set_coord(coords[3], 3.0f, 0.0f, 0.0f);
	return compare_n2_cell("rep cell simple", 4, coords, 0, 0);
}

static int test_cell_matches_n2_negative_coords(void)
{
	fvec3_t coords[8];

	set_coord(coords[0], -1.8f, -0.25f, 0.0f);
	set_coord(coords[1], -0.2f, -0.25f, 0.0f);
	set_coord(coords[2], 0.4f, 0.5f, 0.0f);
	set_coord(coords[3], 4.0f, 0.0f, 0.0f);
	return compare_n2_cell("rep cell negative", 4, coords, 0, 0);
}

static int test_cell_boundary_neighbor_pair(void)
{
	fvec3_t coords[8];

	set_coord(coords[0], 1.9f, 0.0f, 0.0f);
	set_coord(coords[1], 2.1f, 0.0f, 0.0f);
	return compare_n2_cell("rep cell boundary", 2, coords, 0, 0);
}

static int test_cell_beyond_cutoff_zero(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_repulsion_diag diag;
	fvec3_t coords[8], force[8];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], conf.d_r + 0.1f, 0.0f, 0.0f);
	set_coord(force[0], 0.0f, 0.0f, 0.0f);
	set_coord(force[1], 0.0f, 0.0f, 0.0f);
	failed |= compare_n2_cell("rep cell beyond", 2, coords, 0, 0);
	hk_blind_repulsion_accumulate_force_cell_cpu(&conf, 2, coords, 1.0f, 1.0f, 0, 0, force, &diag);
	failed |= check_i64("rep cell beyond active", diag.n_pairs_active, 0);
	return failed;
}

static int test_cell_contact_blocking(void)
{
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	fvec3_t coords[8];

	set_wedge(&edge, 0, 1, 1.0f, 1.0f);
	set_list(&list, &edge, 1);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.5f, 0.0f, 0.0f);
	return compare_n2_cell("rep cell contact block", 2, coords, &list, 0);
}

static int test_cell_backbone_blocking(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	fvec3_t coords[8];

	set_bmap(&bmap, beads, 2);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[2], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 10.0f, 0.0f, 0.0f);
	set_coord(coords[3], 10.0f, 0.0f, 0.0f);
	return compare_n2_cell("rep cell backbone block", 4, coords, 0, &bmap);
}

static int test_cell_homolog_not_blocked(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[1];
	fvec3_t coords[8];

	set_bmap(&bmap, beads, 1);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.5f, 0.0f, 0.0f);
	return compare_n2_cell("rep cell homolog active", 2, coords, 0, &bmap);
}

static int test_relax_step_n2_matches_cell(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag n2_diag, cell_diag;
	fvec3_t n2_coords[4], cell_coords[4];
	char buf[128];
	int failed = 0;
	int i, a;

	hk_fdg_conf_init(&conf);
	set_list(&list, 0, 0);
	set_coord(n2_coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(n2_coords[1], 0.5f, 0.0f, 0.0f);
	set_coord(n2_coords[2], 3.0f, 0.0f, 0.0f);
	set_coord(n2_coords[3], 4.0f, 0.0f, 0.0f);
	memcpy(cell_coords, n2_coords, sizeof(n2_coords));
	failed |= check_i32("rep relax n2 step",
						hk_blind_relax_step_cpu(&conf, &list, 0, 2, n2_coords, 1.0f, 0.1f,
												0.0f, 0.0f, 1, HK_BLIND_REPULSION_N2, 1.0f, &n2_diag), 0);
	failed |= check_i32("rep relax cell step",
						hk_blind_relax_step_cpu(&conf, &list, 0, 2, cell_coords, 1.0f, 0.1f,
												0.0f, 0.0f, 1, HK_BLIND_REPULSION_CELL, 1.0f, &cell_diag), 0);
	failed |= check_close("rep relax energy", cell_diag.repulsion_energy, n2_diag.repulsion_energy);
	failed |= check_i64("rep relax active", cell_diag.n_repulsion_pairs_active, n2_diag.n_repulsion_pairs_active);
	for (i = 0; i < 4; ++i) {
		for (a = 0; a < 3; ++a) {
			snprintf(buf, sizeof(buf), "rep relax coord %d/%d", i, a);
			failed |= check_close(buf, cell_coords[i][a], n2_coords[i][a]);
		}
	}
	return failed;
}

static int test_relax_disabled_preserves_old_behavior(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	fvec3_t coords[2], before[2];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_list(&list, 0, 0);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.5f, 0.0f, 0.0f);
	memcpy(before, coords, sizeof(coords));
	failed |= check_i32("rep disabled step",
						hk_blind_relax_step_cpu(&conf, &list, 0, 1, coords, 1.0f, 0.1f,
												0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, 1.0f, &diag), 0);
	failed |= check_close("rep disabled energy", diag.repulsion_energy, 0.0);
	failed |= check_i64("rep disabled considered", diag.n_repulsion_pairs_considered, 0);
	failed |= check_close("rep disabled coord", dist3(coords[0], before[0]), 0.0);
	return failed;
}

static int test_relax_enabled_pushes_non_attractive_apart(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	fvec3_t coords[2];
	float before, after;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_list(&list, 0, 0);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.5f, 0.0f, 0.0f);
	before = dist3(coords[0], coords[1]);
	failed |= check_i32("rep enabled step",
						hk_blind_relax_step_cpu(&conf, &list, 0, 1, coords, 1.0f, 0.1f, 0.0f, 0.0f, 1, HK_BLIND_REPULSION_N2, 1.0f, &diag), 0);
	after = dist3(coords[0], coords[1]);
	failed |= check_true("rep enabled distance", after > before);
	failed |= check_true("rep enabled active", diag.n_repulsion_pairs_active > 0);
	failed |= check_true("rep enabled energy", diag.repulsion_energy > 0.0f);
	return failed;
}

static int test_relax_contact_pair_blocked(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	fvec3_t coords[2], before[2];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 1.0f, 0.0f, 0.0f);
	memcpy(before, coords, sizeof(coords));
	set_wedge(&edge, 0, 1, 1.0f, 1.0f);
	set_list(&list, &edge, 1);
	failed |= check_i32("rep contact blocked step",
						hk_blind_relax_step_cpu(&conf, &list, 0, 1, coords, 1.0f, 0.1f, 0.0f, 0.0f, 1, HK_BLIND_REPULSION_N2, 1.0f, &diag), 0);
	failed |= check_close("rep contact blocked contact", diag.contact_energy, 0.0);
	failed |= check_close("rep contact blocked rep", diag.repulsion_energy, 0.0);
	failed |= check_i64("rep contact blocked count", diag.n_repulsion_pairs_blocked, 1);
	failed |= check_close("rep contact blocked coords", dist3(coords[0], before[0]), 0.0);
	return failed;
}

static int test_relax_backbone_pair_blocked(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	fvec3_t coords[4], before[4];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_bmap(&bmap, beads, 2);
	set_list(&list, 0, 0);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[2], 1.0f, 0.0f, 0.0f);
	set_coord(coords[1], 10.0f, 0.0f, 0.0f);
	set_coord(coords[3], 11.0f, 0.0f, 0.0f);
	memcpy(before, coords, sizeof(coords));
	failed |= check_i32("rep backbone blocked step",
						hk_blind_relax_step_cpu(&conf, &list, &bmap, bmap.n_beads, coords, 1.0f, 0.1f,
												0.0f, 0.0f, 1, HK_BLIND_REPULSION_N2, 1.0f, &diag), 0);
	failed |= check_close("rep backbone blocked backbone", diag.backbone_energy, 0.0);
	failed |= check_close("rep backbone blocked rep", diag.repulsion_energy, 0.0);
	failed |= check_i64("rep backbone blocked count", diag.n_repulsion_pairs_blocked, 2);
	failed |= check_close("rep backbone blocked coords", dist3(coords[0], before[0]), 0.0);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= test_single_pair_beyond_radius();
	failed |= test_single_pair_inside_radius();
	failed |= test_single_pair_exact_overlap();
	failed |= test_rel_rep_k_scaling();
	failed |= test_schedule();
	failed |= test_batch_counts_and_accumulation();
	failed |= test_contact_blocking();
	failed |= test_backbone_blocking();
	failed |= test_homolog_pair_not_blocked();
	failed |= test_net_force_zero();
	failed |= test_cell_matches_n2_simple();
	failed |= test_cell_matches_n2_negative_coords();
	failed |= test_cell_boundary_neighbor_pair();
	failed |= test_cell_beyond_cutoff_zero();
	failed |= test_cell_contact_blocking();
	failed |= test_cell_backbone_blocking();
	failed |= test_cell_homolog_not_blocked();
	failed |= test_relax_step_n2_matches_cell();
	failed |= test_relax_disabled_preserves_old_behavior();
	failed |= test_relax_enabled_pushes_non_attractive_apart();
	failed |= test_relax_contact_pair_blocked();
	failed |= test_relax_backbone_pair_blocked();
	return failed != 0;
}
