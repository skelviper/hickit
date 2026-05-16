#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hkpriv.h"

static int check_i32(const char *label, int32_t got, int32_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %d, expected %d\n", label, got, expected);
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

static int check_close(const char *label, float got, float expected)
{
	float scale = fabsf(expected) > 1.0f? fabsf(expected) : 1.0f;
	if (fabsf(got - expected) > 1e-5f * scale) {
		fprintf(stderr, "%s: got %.8g, expected %.8g\n", label, got, expected);
		return 1;
	}
	return 0;
}

static struct hk_sdict *make_dict(void)
{
	struct hk_sdict *d = hk_sd_init();
	hk_sd_put(d, "chrSynthetic", 10500000);
	return d;
}

static void destroy_dict(struct hk_sdict *d)
{
	if (d) {
		hk_sd_destroy(d);
		free(d);
	}
}

static int check_map_and_exposure(void)
{
	struct hk_sdict *d = make_dict();
	struct hk_bmap *fine = hk_bmap_gen(d, 0, 0, 1000000, 1);
	struct hk_bmap *coarse = hk_bmap_gen(d, 0, 0, 4000000, 1);
	struct hk_blind_coarse_to_fine_map *map;
	struct hk_blind_bpair_set set;
	struct hk_blind_bpair bpairs[2];
	int failed = 0;

	failed |= check_true("fine bmap exists", fine != 0);
	failed |= check_true("coarse bmap exists", coarse != 0);
	if (fine == 0 || coarse == 0) goto cleanup;
	failed |= check_i32("fine n_beads", fine->n_beads, 10);
	failed |= check_i32("coarse n_beads", coarse->n_beads, 3);
	failed |= check_true("coarse fewer than fine", coarse->n_beads < fine->n_beads);

	map = hk_blind_coarse_to_fine_map_build(coarse, fine);
	failed |= check_true("map exists", map != 0);
	if (map == 0) goto cleanup;
	failed |= check_i32("map n_fine", map->n_fine, fine->n_beads);
	failed |= check_i32("parent0 children", map->parent_child_count[0], 4);
	failed |= check_i32("parent1 children", map->parent_child_count[1], 4);
	failed |= check_i32("terminal parent children", map->parent_child_count[2], 2);
	failed |= check_i32("fine0 parent", map->fine_to_coarse[0], 0);
	failed |= check_i32("fine3 parent", map->fine_to_coarse[3], 0);
	failed |= check_i32("fine4 parent", map->fine_to_coarse[4], 1);
	failed |= check_i32("fine9 parent", map->fine_to_coarse[9], 2);
	failed |= check_i32("fine9 child rank", map->child_rank[9], 1);
	failed |= check_close("cross exposure 4x4", hk_blind_density_exposure_from_child_counts(4, 4, 0), 16.0f);
	failed |= check_close("same parent exposure", hk_blind_density_exposure_from_child_counts(4, 4, 1), 6.0f);
	failed |= check_close("terminal exposure", hk_blind_density_exposure_from_child_counts(2, 4, 0), 8.0f);

	memset(&set, 0, sizeof(set));
	memset(bpairs, 0, sizeof(bpairs));
	set.n_bpairs = 2;
	set.bpairs = bpairs;
	bpairs[0].key.bid[0] = 0;
	bpairs[0].key.bid[1] = 1;
	bpairs[0].n_raw = 32;
	bpairs[1].key.bid[0] = 2;
	bpairs[1].key.bid[1] = 2;
	bpairs[1].n_raw = 30;
	failed |= check_i32("apply exposure", hk_blind_bpair_set_apply_density_exposure_from_map(&set, map), 0);
	failed |= check_i32("bpair0 child0", bpairs[0].density_child_count[0], 4);
	failed |= check_i32("bpair0 child1", bpairs[0].density_child_count[1], 4);
	failed |= check_close("bpair0 exposure", bpairs[0].density_exposure, 16.0f);
	failed |= check_close("density normalized", hk_blind_bpair_density_normalized_raw_count(&bpairs[0], 1e-6f), 2.0f);
	failed |= check_close("raw d_scale old behavior", powf(8.0f, -1.0f / 3.0f), 0.5f);
	failed |= check_close("capped density", hk_blind_bpair_capped_density_raw_count(&bpairs[1], 1e-6f), 30.0f);
	failed |= check_i32("apply density d_scale",
						hk_blind_bpair_set_apply_density_d_scale_from_map(
							&set, map, HK_BLIND_D_SCALE_DENSITY_NORMALIZED_RAW_COUNT, 1e-6f), 0);
	failed |= check_close("density d_scale applied to posterior base",
						  bpairs[0].base_d_scale, powf(2.0f, -1.0f / 3.0f));

	hk_blind_coarse_to_fine_map_destroy(map);
cleanup:
	if (fine) hk_bmap_destroy(fine);
	if (coarse) hk_bmap_destroy(coarse);
	destroy_dict(d);
	return failed;
}

static int check_lift_and_anchor(void)
{
	struct hk_sdict *d = make_dict();
	struct hk_bmap *fine = hk_bmap_gen(d, 0, 0, 1000000, 1);
	struct hk_bmap *coarse = hk_bmap_gen(d, 0, 0, 4000000, 1);
	struct hk_blind_coarse_to_fine_map *map;
	fvec3_t coarse_x[6], fine_x[20], fine_x2[20], force[20];
	int failed = 0, i, a;
	int32_t n_nonfinite = -1;
	float energy;

	if (fine == 0 || coarse == 0) {
		failed |= check_true("bmaps for lift", 0);
		goto cleanup;
	}
	map = hk_blind_coarse_to_fine_map_build(coarse, fine);
	if (map == 0) {
		failed |= check_true("map for lift", 0);
		goto cleanup;
	}
	for (i = 0; i < coarse->n_beads; ++i) {
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY0)][0] = 4.0f * (float)i;
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY0)][1] = 0.0f;
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY0)][2] = 0.0f;
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY1)][0] = 4.0f * (float)i;
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY1)][1] = 2.0f;
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY1)][2] = 0.0f;
	}
	failed |= check_i32("lift ret", hk_blind_lift_from_4mb(coarse, map, coarse_x, fine_x, 0.1f), 0);
	failed |= check_i32("lift ret deterministic", hk_blind_lift_from_4mb(coarse, map, coarse_x, fine_x2, 0.1f), 0);
	for (i = 0; i < fine->n_beads * HK_DIPLOID_N_COPY; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_close("lift deterministic coord", fine_x[i][a], fine_x2[i][a]);
	failed |= check_close("lift first child offset", fine_x[hk_diploid_bid(0, HK_DIPLOID_COPY0)][0], -0.15f);
	failed |= check_close("lift fourth child offset", fine_x[hk_diploid_bid(3, HK_DIPLOID_COPY0)][0], 0.15f);
	failed |= check_close("lift terminal child offset", fine_x[hk_diploid_bid(9, HK_DIPLOID_COPY0)][0], 8.05f);

	for (i = 0; i < fine->n_beads * HK_DIPLOID_N_COPY; ++i)
		for (a = 0; a < 3; ++a)
			force[i][a] = 0.0f;
	for (i = 0; i < fine->n_beads; ++i) {
		if (map->fine_to_coarse[i] == 0)
			fine_x[hk_diploid_bid(i, HK_DIPLOID_COPY0)][0] = 1.0f;
	}
	energy = hk_blind_parent_centroid_anchor_accumulate_force(map, fine_x, coarse_x, 2.0f, force, &n_nonfinite, 0);
	failed |= check_true("anchor energy positive", energy > 0.0f);
	failed |= check_i32("anchor nonfinite", n_nonfinite, 0);
	failed |= check_close("anchor pulls centroid", force[hk_diploid_bid(0, HK_DIPLOID_COPY0)][0], -0.5f);

	hk_blind_coarse_to_fine_map_destroy(map);
cleanup:
	if (fine) hk_bmap_destroy(fine);
	if (coarse) hk_bmap_destroy(coarse);
	destroy_dict(d);
	return failed;
}

static int check_phase_labels_not_copied(void)
{
	struct hk_pair src;
	struct hk_blind_pair dst;
	int failed = 0;
	memset(&src, 0, sizeof(src));
	src.chr = ((uint64_t)0 << 32) | 0u;
	src.pos = ((uint64_t)123 << 32) | 456u;
	src.strand[0] = 1;
	src.strand[1] = -1;
	src.phase[0] = 1;
	src.phase[1] = 0;
	src._.p4[0] = 0.99f;
	hk_blind_pair_from_pair(&dst, &src);
	failed |= check_i32("blind chr0", dst.chr[0], 0);
	failed |= check_i32("blind chr1", dst.chr[1], 0);
	failed |= check_i32("blind pos0", dst.pos[0], 123);
	failed |= check_i32("blind pos1", dst.pos[1], 456);
	failed |= check_i32("blind strand0", dst.strand[0], 1);
	failed |= check_i32("blind strand1", dst.strand[1], -1);
	return failed;
}

int main(void)
{
	int failed = 0;
	hk_verbose = 0;
	failed |= check_map_and_exposure();
	failed |= check_lift_and_anchor();
	failed |= check_phase_labels_not_copied();
	return failed? 1 : 0;
}
