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

static void copy_force(fvec3_t *dst, const fvec3_t *src, int32_t n)
{
	int32_t i;
	for (i = 0; i < n; ++i)
		set_coord(dst[i], src[i][0], src[i][1], src[i][2]);
}

static int check_vec_close(const char *label, const fvec3_t got, float x, float y, float z)
{
	char buf[128];
	int failed = 0;
	snprintf(buf, sizeof(buf), "%s x", label);
	failed |= check_close(buf, got[0], x);
	snprintf(buf, sizeof(buf), "%s y", label);
	failed |= check_close(buf, got[1], y);
	snprintf(buf, sizeof(buf), "%s z", label);
	failed |= check_close(buf, got[2], z);
	return failed;
}

static int check_single_edge_action_reaction(void)
{
	struct hk_fdg_conf conf;
	fvec3_t x0, x1, f0, f1;
	float e;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(x0, 0.0f, 0.0f, 0.0f);
	set_coord(x1, 2.0f, 0.0f, 0.0f);
	e = hk_blind_backbone_energy_force_cpu(&conf, x0, x1, 1.0f, 1.0f, f0, f1);
	failed |= check_close("single long energy", e, (2.0f - conf.d_b2) * (2.0f - conf.d_b2));
	failed |= check_vec_close("single long f0", f0, 2.0f * (2.0f - conf.d_b2), 0.0f, 0.0f);
	failed |= check_vec_close("single long f1", f1, -2.0f * (2.0f - conf.d_b2), 0.0f, 0.0f);
	failed |= check_close("single action reaction x", f0[0] + f1[0], 0.0f);
	failed |= check_close("single action reaction y", f0[1] + f1[1], 0.0f);
	failed |= check_close("single action reaction z", f0[2] + f1[2], 0.0f);
	failed |= check_true("single energy finite", isfinite(e));
	failed |= check_true("single force finite", isfinite(f0[0]) && isfinite(f1[0]));
	return failed;
}

static int check_exact_overlap_safety(void)
{
	struct hk_fdg_conf conf;
	fvec3_t x0, x1, f0, f1;
	float e;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(x0, 1.0f, 2.0f, 3.0f);
	set_coord(x1, 1.0f, 2.0f, 3.0f);
	e = hk_blind_backbone_energy_force_cpu(&conf, x0, x1, 1.0f, 1.0f, f0, f1);
	failed |= check_close("overlap energy", e, conf.d_b1 * conf.d_b1);
	failed |= check_vec_close("overlap f0", f0, 2.0f * conf.d_b1, 0.0f, 0.0f);
	failed |= check_vec_close("overlap f1", f1, -2.0f * conf.d_b1, 0.0f, 0.0f);
	failed |= check_true("overlap energy finite", isfinite(e));
	failed |= check_true("overlap force finite", isfinite(f0[0]) && isfinite(f1[0]));
	return failed;
}

static int check_single_edge_branches(void)
{
	struct hk_fdg_conf conf;
	fvec3_t x0, x1, f0, f1;
	float e, t;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	set_coord(x0, 0.0f, 0.0f, 0.0f);

	set_coord(x1, 0.05f, 0.0f, 0.0f);
	t = conf.d_b1 - 0.05f;
	e = hk_blind_backbone_energy_force_cpu(&conf, x0, x1, 1.0f, 1.0f, f0, f1);
	failed |= check_close("short energy", e, t * t);
	failed |= check_vec_close("short f0", f0, -2.0f * t, 0.0f, 0.0f);
	failed |= check_vec_close("short f1", f1, 2.0f * t, 0.0f, 0.0f);

	set_coord(x1, conf.d_b1, 0.0f, 0.0f);
	e = hk_blind_backbone_energy_force_cpu(&conf, x0, x1, 1.0f, 1.0f, f0, f1);
	failed |= check_close("lower boundary energy", e, 0.0f);
	failed |= check_vec_close("lower boundary f0", f0, 0.0f, 0.0f, 0.0f);

	set_coord(x1, 1.0f, 0.0f, 0.0f);
	e = hk_blind_backbone_energy_force_cpu(&conf, x0, x1, 1.0f, 1.0f, f0, f1);
	failed |= check_close("plateau energy", e, 0.0f);
	failed |= check_vec_close("plateau f0", f0, 0.0f, 0.0f, 0.0f);
	failed |= check_vec_close("plateau f1", f1, 0.0f, 0.0f, 0.0f);

	set_coord(x1, conf.d_b2, 0.0f, 0.0f);
	e = hk_blind_backbone_energy_force_cpu(&conf, x0, x1, 1.0f, 1.0f, f0, f1);
	failed |= check_close("upper boundary energy", e, 0.0f);
	failed |= check_vec_close("upper boundary f0", f0, 0.0f, 0.0f, 0.0f);

	set_coord(x1, 3.0f, 0.0f, 0.0f);
	t = 3.0f / 2.0f - conf.d_b2;
	e = hk_blind_backbone_energy_force_cpu(&conf, x0, x1, 1.0f, 2.0f, f0, f1);
	failed |= check_close("dscale energy", e, t * t);
	failed |= check_vec_close("dscale f0", f0, 2.0f * t, 0.0f, 0.0f);
	failed |= check_vec_close("dscale f1", f1, -2.0f * t, 0.0f, 0.0f);
	return failed;
}

static int check_accumulation_semantics(void)
{
	struct hk_fdg_conf conf;
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_backbone_diag diag;
	fvec3_t coords[4], force[4], before[4], f0, f1;
	float e_single, e_batch;
	int failed = 0;
	int i;

	hk_fdg_conf_init(&conf);
	set_bmap(&bmap, beads, 2, 0);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 2.0f, 0.0f, 0.0f);
	set_coord(coords[3], 2.0f, 1.0f, 0.0f);
	for (i = 0; i < 4; ++i)
		set_coord(force[i], 10.0f + (float)i, 20.0f + (float)i, 30.0f + (float)i);
	copy_force(before, force, 4);

	e_batch = hk_blind_backbone_accumulate_force_cpu(&conf, &bmap, coords, 1.0f, force, &diag);
	e_single = hk_blind_backbone_energy_force_cpu(&conf, coords[0], coords[2], 1.0f, 1.0f, f0, f1);
	failed |= check_i64("accum n_edges", diag.n_edges, 2);
	failed |= check_i32("accum nonfinite", diag.n_nonfinite, 0);
	failed |= check_close("accum energy ret", e_batch, 2.0f * e_single);
	failed |= check_close("accum energy diag", diag.energy, e_batch);
	failed |= check_close("accum force l1", diag.force_l1, 4.0f * fabsf(f0[0]));
	failed |= check_vec_close("accum copy0 bead0", force[0], before[0][0] + f0[0], before[0][1] + f0[1], before[0][2] + f0[2]);
	failed |= check_vec_close("accum copy0 bead1", force[2], before[2][0] + f1[0], before[2][1] + f1[1], before[2][2] + f1[2]);
	failed |= check_vec_close("accum copy1 bead0", force[1], before[1][0] + f0[0], before[1][1] + f0[1], before[1][2] + f0[2]);
	failed |= check_vec_close("accum copy1 bead1", force[3], before[3][0] + f1[0], before[3][1] + f1[1], before[3][2] + f1[2]);
	return failed;
}

static int check_chain_edge_counts(void)
{
	struct hk_fdg_conf conf;
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	struct hk_blind_backbone_diag diag;
	fvec3_t coords[8] = {{0}}, force[8] = {{0}};
	int32_t chr_boundary[4] = {0, 0, 1, 1};
	int failed = 0;
	int i;

	hk_fdg_conf_init(&conf);
	for (i = 0; i < 4; ++i) {
		set_coord(coords[hk_diploid_bid(i, HK_DIPLOID_COPY0)], (float)i, 0.0f, 0.0f);
		set_coord(coords[hk_diploid_bid(i, HK_DIPLOID_COPY1)], (float)i, 1.0f, 0.0f);
	}

	set_bmap(&bmap, beads, 4, 0);
	failed |= check_close("same chr energy", hk_blind_backbone_accumulate_force_cpu(&conf, &bmap, coords, 1.0f, force, &diag), 0.0f);
	failed |= check_i64("same chr edges", diag.n_edges, 6);
	failed |= check_i32("same chr nonfinite", diag.n_nonfinite, 0);

	memset(force, 0, sizeof(force));
	set_bmap(&bmap, beads, 4, chr_boundary);
	failed |= check_close("boundary energy", hk_blind_backbone_accumulate_force_cpu(&conf, &bmap, coords, 1.0f, force, &diag), 0.0f);
	failed |= check_i64("boundary edges", diag.n_edges, 4);
	failed |= check_i32("boundary nonfinite", diag.n_nonfinite, 0);
	return failed;
}

static int check_batch_energy_and_net_force(void)
{
	struct hk_fdg_conf conf;
	struct hk_bmap bmap;
	struct hk_bead beads[2];
	struct hk_blind_backbone_diag diag;
	fvec3_t coords[4], force[4] = {{0}}, f0, f1;
	float single, total;
	int failed = 0;
	int i, a;
	float net[3] = {0.0f, 0.0f, 0.0f};

	hk_fdg_conf_init(&conf);
	set_bmap(&bmap, beads, 2, 0);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(coords[2], 2.0f, 0.0f, 0.0f);
	set_coord(coords[3], 2.0f, 1.0f, 0.0f);

	single = hk_blind_backbone_energy_force_cpu(&conf, coords[0], coords[2], 1.0f, 1.0f, f0, f1);
	total = hk_blind_backbone_accumulate_force_cpu(&conf, &bmap, coords, 1.0f, force, &diag);
	failed |= check_close("batch sum energy", total, 2.0f * single);
	failed |= check_close("batch diag energy", diag.energy, total);
	for (i = 0; i < 4; ++i)
		for (a = 0; a < 3; ++a)
			net[a] += force[i][a];
	failed |= check_close("batch net force x", net[0], 0.0f);
	failed |= check_close("batch net force y", net[1], 0.0f);
	failed |= check_close("batch net force z", net[2], 0.0f);
	return failed;
}

static int check_empty_and_single_bead(void)
{
	struct hk_fdg_conf conf;
	struct hk_bmap bmap;
	struct hk_bead bead;
	struct hk_blind_backbone_diag diag;
	fvec3_t coords[2], force[2];
	int failed = 0;

	hk_fdg_conf_init(&conf);
	memset(&bmap, 0, sizeof(bmap));
	failed |= check_close("empty energy", hk_blind_backbone_accumulate_force_cpu(&conf, &bmap, 0, 1.0f, 0, &diag), 0.0f);
	failed |= check_i64("empty edges", diag.n_edges, 0);
	failed |= check_i32("empty nonfinite", diag.n_nonfinite, 0);

	set_bmap(&bmap, &bead, 1, 0);
	set_coord(coords[0], 0.0f, 0.0f, 0.0f);
	set_coord(coords[1], 0.0f, 1.0f, 0.0f);
	set_coord(force[0], 1.0f, 2.0f, 3.0f);
	set_coord(force[1], 4.0f, 5.0f, 6.0f);
	failed |= check_close("single bead energy", hk_blind_backbone_accumulate_force_cpu(&conf, &bmap, coords, 1.0f, force, &diag), 0.0f);
	failed |= check_i64("single bead edges", diag.n_edges, 0);
	failed |= check_i32("single bead nonfinite", diag.n_nonfinite, 0);
	failed |= check_vec_close("single bead force0 unchanged", force[0], 1.0f, 2.0f, 3.0f);
	failed |= check_vec_close("single bead force1 unchanged", force[1], 4.0f, 5.0f, 6.0f);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_single_edge_action_reaction();
	failed |= check_exact_overlap_safety();
	failed |= check_single_edge_branches();
	failed |= check_accumulation_semantics();
	failed |= check_chain_edge_counts();
	failed |= check_batch_energy_and_net_force();
	failed |= check_empty_and_single_bead();
	return failed != 0;
}
