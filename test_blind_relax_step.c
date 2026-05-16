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

static void set_bmap(struct hk_bmap *bmap, struct hk_bead *beads, int32_t n_beads, const int32_t *chr)
{
	int32_t i;
	memset(bmap, 0, sizeof(*bmap));
	bmap->n_beads = n_beads;
	bmap->beads = beads;
	for (i = 0; i < n_beads; ++i) {
		beads[i].chr = chr? chr[i] : 0;
		beads[i].st = i * 1000000;
		beads[i].en = beads[i].st + 1000000;
	}
}

static int check_clean_diag(const char *label, const struct hk_blind_step_diag *diag)
{
	int failed = 0;
	char buf[128];
	snprintf(buf, sizeof(buf), "%s contact energy finite", label);
	failed |= check_true(buf, isfinite(diag->contact_energy));
	snprintf(buf, sizeof(buf), "%s backbone energy finite", label);
	failed |= check_true(buf, isfinite(diag->backbone_energy));
	snprintf(buf, sizeof(buf), "%s repulsion energy finite", label);
	failed |= check_true(buf, isfinite(diag->repulsion_energy));
	snprintf(buf, sizeof(buf), "%s sep energy finite", label);
	failed |= check_true(buf, isfinite(diag->sep_energy));
	snprintf(buf, sizeof(buf), "%s total energy finite", label);
	failed |= check_true(buf, isfinite(diag->total_energy));
	snprintf(buf, sizeof(buf), "%s force l1 finite", label);
	failed |= check_true(buf, isfinite(diag->force_l1));
	snprintf(buf, sizeof(buf), "%s backbone force l1 finite", label);
	failed |= check_true(buf, isfinite(diag->backbone_force_l1));
	snprintf(buf, sizeof(buf), "%s repulsion force l1 finite", label);
	failed |= check_true(buf, isfinite(diag->repulsion_force_l1));
	snprintf(buf, sizeof(buf), "%s total energy sum", label);
	failed |= check_close(buf, diag->total_energy,
						  diag->contact_energy + diag->backbone_energy + diag->repulsion_energy + diag->sep_energy);
	snprintf(buf, sizeof(buf), "%s force nonfinite", label);
	failed |= check_i32(buf, diag->n_force_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s contact nonfinite", label);
	failed |= check_i32(buf, diag->n_contact_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s backbone nonfinite", label);
	failed |= check_i32(buf, diag->n_backbone_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s repulsion nonfinite", label);
	failed |= check_i32(buf, diag->n_repulsion_nonfinite, 0);
	return failed;
}

static int check_step_zero(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	fvec3_t coords[2], before[2];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 2.5f, 0.0f, 0.0f);
	copy_coords(before, coords, 2);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	failed |= check_i32("step-zero ret", hk_blind_relax_step_cpu(&conf, &list, 0, 1, coords, 1.0f, 0.0f,
																 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, 0.0f, &diag), 0);
	failed |= check_coords_close("step-zero coords", coords, before, 2);
	failed |= check_true("step-zero contact energy", diag.contact_energy > 0.0f);
	failed |= check_close("step-zero backbone energy", diag.backbone_energy, 0.0f);
	failed |= check_close("step-zero backbone force", diag.backbone_force_l1, 0.0f);
	failed |= check_i64("step-zero backbone edges", diag.n_backbone_edges, 0);
	failed |= check_close("step-zero sep energy", diag.sep_energy, 0.0f);
	failed |= check_true("step-zero force l1", diag.force_l1 > 0.0f);
	failed |= check_clean_diag("step-zero diag", &diag);
	return failed;
}

static int check_long_distance_pulls_closer(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	fvec3_t coords[2];
	float before, after;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 2.5f, 0.0f, 0.0f);
	before = dist3(coords[0], coords[1]);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	failed |= check_i32("long-step ret", hk_blind_relax_step_cpu(&conf, &list, 0, 1, coords, 1.0f, 0.01f,
																0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, 0.0f, &diag), 0);
	after = dist3(coords[0], coords[1]);
	failed |= check_true("long-step pulls closer", after < before);
	failed |= check_true("long-step contact active", diag.contact_energy > 0.0f);
	failed |= check_clean_diag("long-step diag", &diag);
	failed |= check_coords_finite(coords, 2);
	return failed;
}

static int check_short_distance_pushes_apart(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	fvec3_t coords[2];
	float before, after;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.25f, 0.0f, 0.0f);
	before = dist3(coords[0], coords[1]);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	failed |= check_i32("short-step ret", hk_blind_relax_step_cpu(&conf, &list, 0, 1, coords, 1.0f, 0.01f,
																 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, 0.0f, &diag), 0);
	after = dist3(coords[0], coords[1]);
	failed |= check_true("short-step pushes apart", after > before);
	failed |= check_true("short-step contact active", diag.contact_energy > 0.0f);
	failed |= check_clean_diag("short-step diag", &diag);
	failed |= check_coords_finite(coords, 2);
	return failed;
}

static int check_plateau_no_movement(void)
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
	copy_coords(before, coords, 2);
	set_edge(&edge, 0, 1, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	failed |= check_i32("plateau-step ret", hk_blind_relax_step_cpu(&conf, &list, 0, 1, coords, 1.0f, 0.05f,
																	0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, 0.0f, &diag), 0);
	failed |= check_close("plateau-step contact energy", diag.contact_energy, 0.0f);
	failed |= check_close("plateau-step force l1", diag.force_l1, 0.0f);
	failed |= check_coords_close("plateau-step coords", coords, before, 2);
	failed |= check_clean_diag("plateau-step diag", &diag);
	return failed;
}

static int check_anti_collapse_pushes_apart(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	fvec3_t coords[2];
	float before, after;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.25f, 0.0f, 0.0f);
	before = dist3(coords[0], coords[1]);
	set_list(&list, 0, 0);

	failed |= check_i32("anti-step ret", hk_blind_relax_step_cpu(&conf, &list, 0, 1, coords, 1.0f, 1.0f,
																0.5f, 0.05f, 0, HK_BLIND_REPULSION_NONE, 0.0f, &diag), 0);
	after = dist3(coords[0], coords[1]);
	failed |= check_close("anti-step contact energy", diag.contact_energy, 0.0f);
	failed |= check_true("anti-step sep energy", diag.sep_energy > 0.0f);
	failed |= check_true("anti-step pushes apart", after > before);
	failed |= check_clean_diag("anti-step diag", &diag);
	failed |= check_coords_finite(coords, 2);
	return failed;
}

static int check_contact_and_anti_collapse_combine(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	fvec3_t coords[4];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.25f, 0.0f, 0.0f);
	set_coord(coords[2], 2.0f, 0.0f, 0.0f);
	set_coord(coords[3], 4.5f, 0.0f, 0.0f);
	set_edge(&edge, 2, 3, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);

	failed |= check_i32("combined-step ret", hk_blind_relax_step_cpu(&conf, &list, 0, 2, coords, 1.0f, 0.01f,
																	0.5f, 0.05f, 0, HK_BLIND_REPULSION_NONE, 0.0f, &diag), 0);
	failed |= check_true("combined contact energy", diag.contact_energy > 0.0f);
	failed |= check_true("combined sep energy", diag.sep_energy > 0.0f);
	failed |= check_close("combined total", diag.total_energy,
						  diag.contact_energy + diag.backbone_energy + diag.repulsion_energy + diag.sep_energy);
	failed |= check_true("combined force active", diag.force_l1 > 0.0f);
	failed |= check_clean_diag("combined diag", &diag);
	failed |= check_coords_finite(coords, 4);
	return failed;
}

static int check_backbone_long_distance_pulls_closer(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	fvec3_t coords[4];
	float before, after;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_bmap(&bmap, beads, 2, 0);
	set_list(&list, 0, 0);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 2.0f, 0.0f, 0.0f);
	set_coord(coords[3], 2.0f, 1.0f, 0.0f);
	before = dist3(coords[0], coords[2]);
	failed |= check_i32("backbone-long ret", hk_blind_relax_step_cpu(&conf, &list, &bmap, bmap.n_beads, coords,
																	 1.0f, 0.01f, 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, 0.0f, &diag), 0);
	after = dist3(coords[0], coords[2]);
	failed |= check_true("backbone-long pulls closer", after < before);
	failed |= check_close("backbone-long contact", diag.contact_energy, 0.0f);
	failed |= check_true("backbone-long energy", diag.backbone_energy > 0.0f);
	failed |= check_i64("backbone-long edges", diag.n_backbone_edges, 2);
	failed |= check_true("backbone-long force", diag.backbone_force_l1 > 0.0f);
	failed |= check_clean_diag("backbone-long diag", &diag);
	failed |= check_coords_finite(coords, 4);
	return failed;
}

static int check_backbone_short_distance_pushes_apart(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	fvec3_t coords[4];
	float before, after;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_bmap(&bmap, beads, 2, 0);
	set_list(&list, 0, 0);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 0.05f, 0.0f, 0.0f);
	set_coord(coords[3], 0.05f, 1.0f, 0.0f);
	before = dist3(coords[0], coords[2]);
	failed |= check_i32("backbone-short ret", hk_blind_relax_step_cpu(&conf, &list, &bmap, bmap.n_beads, coords,
																	  1.0f, 0.1f, 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, 0.0f, &diag), 0);
	after = dist3(coords[0], coords[2]);
	failed |= check_true("backbone-short pushes apart", after > before);
	failed |= check_true("backbone-short energy", diag.backbone_energy > 0.0f);
	failed |= check_i64("backbone-short edges", diag.n_backbone_edges, 2);
	failed |= check_clean_diag("backbone-short diag", &diag);
	return failed;
}

static int check_backbone_plateau_no_movement(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	fvec3_t coords[4], before[4];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_bmap(&bmap, beads, 2, 0);
	set_list(&list, 0, 0);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 1.0f, 0.0f, 0.0f);
	set_coord(coords[3], 1.0f, 1.0f, 0.0f);
	copy_coords(before, coords, 4);
	failed |= check_i32("backbone-plateau ret", hk_blind_relax_step_cpu(&conf, &list, &bmap, bmap.n_beads, coords,
																		1.0f, 0.1f, 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, 0.0f, &diag), 0);
	failed |= check_close("backbone-plateau energy", diag.backbone_energy, 0.0f);
	failed |= check_close("backbone-plateau force", diag.backbone_force_l1, 0.0f);
	failed |= check_i64("backbone-plateau edges", diag.n_backbone_edges, 2);
	failed |= check_coords_close("backbone-plateau coords", coords, before, 4);
	failed |= check_clean_diag("backbone-plateau diag", &diag);
	return failed;
}

static int check_no_cross_chromosome_backbone(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	fvec3_t coords[8];
	int32_t chr[4] = {0, 0, 1, 1};
	int failed = 0;
	int i;

	hk_fdg_conf_init(&conf);
	set_bmap(&bmap, beads, 4, chr);
	set_list(&list, 0, 0);
	for (i = 0; i < 4; ++i) {
		set_coord(coords[hk_diploid_bid(i, HK_DIPLOID_COPY0)], (float)i, 0.0f, 0.0f);
		set_coord(coords[hk_diploid_bid(i, HK_DIPLOID_COPY1)], (float)i, 1.0f, 0.0f);
	}
	failed |= check_i32("backbone-boundary ret", hk_blind_relax_step_cpu(&conf, &list, &bmap, bmap.n_beads, coords,
																		 1.0f, 0.0f, 0.0f, 0.0f, 0, HK_BLIND_REPULSION_NONE, 0.0f, &diag), 0);
	failed |= check_i64("backbone-boundary edges", diag.n_backbone_edges, 4);
	failed |= check_clean_diag("backbone-boundary diag", &diag);
	return failed;
}

static int check_contact_backbone_and_anti_collapse_combine(void)
{
	struct hk_fdg_conf conf;
	struct hk_blind_wedge edge;
	struct hk_blind_wedge_list list;
	struct hk_blind_step_diag diag;
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	fvec3_t coords[4];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_bmap(&bmap, beads, 2, 0);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.25f, 0.0f, 0.0f);
	set_coord(coords[2], 3.0f, 0.0f, 0.0f);
	set_coord(coords[3], 3.0f, 1.0f, 0.0f);
	set_edge(&edge, 1, 3, 2.0f, 1.0f, HK_BLIND_STATE_00);
	set_list(&list, &edge, 1);
	failed |= check_i32("all-force ret", hk_blind_relax_step_cpu(&conf, &list, &bmap, bmap.n_beads, coords,
																 1.0f, 0.01f, 0.5f, 0.05f, 0, HK_BLIND_REPULSION_NONE, 0.0f, &diag), 0);
	failed |= check_true("all-force contact", diag.contact_energy > 0.0f);
	failed |= check_true("all-force backbone", diag.backbone_energy > 0.0f);
	failed |= check_true("all-force sep", diag.sep_energy > 0.0f);
	failed |= check_close("all-force total", diag.total_energy,
						  diag.contact_energy + diag.backbone_energy + diag.repulsion_energy + diag.sep_energy);
	failed |= check_true("all-force force", diag.force_l1 > 0.0f);
	failed |= check_true("all-force backbone force", diag.backbone_force_l1 > 0.0f);
	failed |= check_clean_diag("all-force diag", &diag);
	failed |= check_coords_finite(coords, 4);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_step_zero();
	failed |= check_long_distance_pulls_closer();
	failed |= check_short_distance_pushes_apart();
	failed |= check_plateau_no_movement();
	failed |= check_anti_collapse_pushes_apart();
	failed |= check_contact_and_anti_collapse_combine();
	failed |= check_backbone_long_distance_pulls_closer();
	failed |= check_backbone_short_distance_pushes_apart();
	failed |= check_backbone_plateau_no_movement();
	failed |= check_no_cross_chromosome_backbone();
	failed |= check_contact_backbone_and_anti_collapse_combine();
	return failed != 0;
}
