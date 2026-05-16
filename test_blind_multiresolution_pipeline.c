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

static int check_not_close(const char *label, float got, float forbidden)
{
	float scale = fabsf(forbidden) > 1.0f? fabsf(forbidden) : 1.0f;
	if (fabsf(got - forbidden) <= 1e-5f * scale) {
		fprintf(stderr, "%s: got %.8g, unexpectedly matched %.8g\n", label, got, forbidden);
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

static struct hk_pair make_pair(int32_t chr0, int32_t pos0, int32_t chr1, int32_t pos1,
								int8_t phase0, int8_t phase1)
{
	struct hk_pair p;
	memset(&p, 0, sizeof(p));
	p.chr = ((uint64_t)(uint32_t)chr0 << 32) | (uint32_t)chr1;
	p.pos = ((uint64_t)(uint32_t)pos0 << 32) | (uint32_t)pos1;
	p.strand[0] = 1;
	p.strand[1] = -1;
	p.phase[0] = phase0;
	p.phase[1] = phase1;
	p._.p4[0] = 0.97f;
	p._.p4[1] = 0.01f;
	p._.p4[2] = 0.01f;
	p._.p4[3] = 0.01f;
	return p;
}

static void set_anchor_test_schedule_conf(struct hk_blind_iter_schedule_conf *conf)
{
	memset(conf, 0, sizeof(*conf));
	conf->n_iter = 1;
	conf->base_conf.unit = 1.0f;
	conf->base_conf.d_scale = 1.0f;
	conf->base_conf.base_k = 1.0f;
	conf->base_conf.temperature = 1.0f;
	conf->base_conf.rho_train = 1.0f;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.min_sep_unit = 0.0f;
	conf->base_conf.lambda_sep = 0.0f;
	conf->base_conf.relax_step = 0.01f;
	conf->base_conf.relax_steps = 2;
	conf->base_conf.enable_repulsion = 0;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_NONE;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->base_conf.d_scale_eps_count = 1e-6f;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = 1.0f;
	conf->temperature_end = 1.0f;
	conf->rho_train_start = 1.0f;
	conf->rho_train_end = 1.0f;
}

static void copy_diploid_coords(fvec3_t *dst, const fvec3_t *src, int32_t n_diploid)
{
	int32_t i;
	int a;
	for (i = 0; i < n_diploid; ++i)
		for (a = 0; a < 3; ++a)
			dst[i][a] = src[i][a];
}

static int check_phase_labels_not_copied(void)
{
	struct hk_pair src_a, src_b;
	struct hk_blind_pair blind_a, blind_b;
	int failed = 0;

	src_a = make_pair(0, 123, 0, 456, 0, 1);
	src_b = src_a;
	src_b.phase[0] = 1;
	src_b.phase[1] = 0;
	src_b._.p4[0] = 0.05f;
	src_b._.p4[1] = 0.90f;

	memset(&blind_a, 0x7f, sizeof(blind_a));
	memset(&blind_b, 0x3f, sizeof(blind_b));
	hk_blind_pair_from_pair(&blind_a, &src_a);
	hk_blind_pair_from_pair(&blind_b, &src_b);

	// Blind input carries only endpoint coordinates and strands; phase-like labels are intentionally ignored.
	failed |= check_i32("blind chr0", blind_a.chr[0], 0);
	failed |= check_i32("blind chr1", blind_a.chr[1], 0);
	failed |= check_i32("blind pos0", blind_a.pos[0], 123);
	failed |= check_i32("blind pos1", blind_a.pos[1], 456);
	failed |= check_i32("blind strand0", blind_a.strand[0], 1);
	failed |= check_i32("blind strand1", blind_a.strand[1], -1);
	failed |= check_i32("phase-independent chr0", blind_b.chr[0], blind_a.chr[0]);
	failed |= check_i32("phase-independent chr1", blind_b.chr[1], blind_a.chr[1]);
	failed |= check_i32("phase-independent pos0", blind_b.pos[0], blind_a.pos[0]);
	failed |= check_i32("phase-independent pos1", blind_b.pos[1], blind_a.pos[1]);
	failed |= check_i32("phase-independent strand0", blind_b.strand[0], blind_a.strand[0]);
	failed |= check_i32("phase-independent strand1", blind_b.strand[1], blind_a.strand[1]);
	return failed;
}

static int fill_raw_pairs_from_phased_sources(struct hk_blind_pair *raw, int32_t n_raw)
{
	int32_t i;
	for (i = 0; i < n_raw; ++i) {
		struct hk_pair src = make_pair(0, 1000000, 0, 5000000, (int8_t)(i & 1), (int8_t)((i + 1) & 1));
		hk_blind_pair_from_pair(&raw[i], &src);
	}
	return 0;
}

static int find_bpair_id(const struct hk_blind_bpair_set *set, int32_t bid0, int32_t bid1)
{
	int32_t i;
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		if (bp->key.bid[0] == bid0 && bp->key.bid[1] == bid1)
			return i;
	}
	return -1;
}

static int check_density_normalized_edge_expansion(struct hk_blind_bpair *bp)
{
	struct hk_blind_wedge raw_edges[HK_BLIND_N_STATE];
	struct hk_blind_wedge density_edges[HK_BLIND_N_STATE];
	float raw_d_scale, density_d_scale;
	int failed = 0;

	bp->p4[HK_BLIND_STATE_00] = 0.7f;
	bp->p4[HK_BLIND_STATE_01] = 0.1f;
	bp->p4[HK_BLIND_STATE_10] = 0.1f;
	bp->p4[HK_BLIND_STATE_11] = 0.1f;

	hk_blind_bpair_expand_weighted_edges_mode_ex(bp, 2.0f, bp->base_d_scale,
												 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT, 0.0f,
												 HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f,
												 1.0f, 1.0f, raw_edges);
	hk_blind_bpair_expand_weighted_edges_mode_ex(bp, 2.0f, bp->base_d_scale,
												 1.0f, HK_BLIND_RHO_TRAIN_CONSTANT, 0.0f,
												 HK_BLIND_D_SCALE_DENSITY_NORMALIZED_RAW_COUNT, 1e-6f,
												 1.0f, 1.0f, density_edges);

	raw_d_scale = powf((float)bp->n_raw, -1.0f / 3.0f);
	density_d_scale = powf(hk_blind_bpair_density_normalized_raw_count(bp, 1e-6f), -1.0f / 3.0f);
	failed |= check_close("raw-count d_scale", raw_edges[HK_BLIND_STATE_00].d_scale, raw_d_scale);
	failed |= check_close("density-normalized d_scale", density_edges[HK_BLIND_STATE_00].d_scale, density_d_scale);
	failed |= check_not_close("density-normalized differs from raw count",
							  density_edges[HK_BLIND_STATE_00].d_scale,
							  raw_edges[HK_BLIND_STATE_00].d_scale);
	failed |= check_i32("density edge bid0", density_edges[HK_BLIND_STATE_00].bid[0],
						hk_diploid_bid(bp->key.bid[0], HK_DIPLOID_COPY0));
	failed |= check_i32("density edge bid1", density_edges[HK_BLIND_STATE_00].bid[1],
						hk_diploid_bid(bp->key.bid[1], HK_DIPLOID_COPY0));
	return failed;
}

static int check_lift_and_anchor(const struct hk_bmap *coarse, const struct hk_bmap *fine,
								 const struct hk_blind_coarse_to_fine_map *map)
{
	fvec3_t coarse_x[6], fine_x[20], force[20], old_x[20], new_x[20];
	struct hk_fdg_conf conf;
	struct hk_blind_wedge_list empty_edges;
	struct hk_blind_step_diag step_diag;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_iter_loop_diag old_loop_diag, new_loop_diag;
	struct hk_blind_bpair_set old_empty_set, new_empty_set;
	float log_prior[HK_BLIND_N_STATE];
	int32_t n_nonfinite = -1;
	float energy, anchor_force_l1 = 0.0f;
	int failed = 0, i, a, ret;

	for (i = 0; i < coarse->n_beads; ++i) {
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY0)][0] = 4.0f * (float)i;
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY0)][1] = 0.0f;
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY0)][2] = 0.0f;
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY1)][0] = 4.0f * (float)i;
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY1)][1] = 2.0f;
		coarse_x[hk_diploid_bid(i, HK_DIPLOID_COPY1)][2] = 0.0f;
	}

	failed |= check_i32("lift ret", hk_blind_lift_from_4mb(coarse, map, coarse_x, fine_x, 0.1f), 0);
	failed |= check_close("lift child0 x", fine_x[hk_diploid_bid(0, HK_DIPLOID_COPY0)][0], -0.15f);
	failed |= check_close("lift child3 x", fine_x[hk_diploid_bid(3, HK_DIPLOID_COPY0)][0], 0.15f);
	failed |= check_close("lift terminal child x", fine_x[hk_diploid_bid(9, HK_DIPLOID_COPY0)][0], 8.05f);

	for (i = 0; i < fine->n_beads * HK_DIPLOID_N_COPY; ++i)
		for (a = 0; a < 3; ++a)
			force[i][a] = 0.0f;
	for (i = 0; i < fine->n_beads; ++i) {
		if (map->fine_to_coarse[i] == 0)
			fine_x[hk_diploid_bid(i, HK_DIPLOID_COPY0)][0] = 1.0f;
	}
	energy = hk_blind_parent_centroid_anchor_accumulate_force(map, fine_x, coarse_x, 2.0f,
															  force, &n_nonfinite, &anchor_force_l1);
	failed |= check_true("anchor energy positive", energy > 0.0f);
	failed |= check_i32("anchor nonfinite", n_nonfinite, 0);
	failed |= check_true("anchor force l1 positive", anchor_force_l1 > 0.0f);
	failed |= check_close("anchor pulls first child", force[hk_diploid_bid(0, HK_DIPLOID_COPY0)][0], -0.5f);

	hk_fdg_conf_init(&conf);
	hk_blind_wedge_list_init(&empty_edges);
	ret = hk_blind_relax_step_parent_anchor_cpu(&conf, &empty_edges, 0, fine->n_beads, fine_x,
												1.0f, 0.1f, 0.0f, 0.0f, 0,
												HK_BLIND_REPULSION_NONE, 0.0f,
												map, coarse_x, 2.0f, &step_diag);
	failed |= check_i32("anchored relax step ret", ret, 0);
	failed |= check_true("anchored relax step energy", step_diag.anchor_energy > 0.0f);
	failed |= check_true("anchored relax step force", step_diag.anchor_force_l1 > 0.0f);
	failed |= check_i32("anchored relax step nonfinite", step_diag.n_anchor_nonfinite, 0);
	failed |= check_close("anchored relax step moves child",
						  fine_x[hk_diploid_bid(0, HK_DIPLOID_COPY0)][0], 0.95f);

	failed |= check_i32("relift ret", hk_blind_lift_from_4mb(coarse, map, coarse_x, fine_x, 0.1f), 0);
	copy_diploid_coords(old_x, fine_x, fine->n_beads * HK_DIPLOID_N_COPY);
	copy_diploid_coords(new_x, fine_x, fine->n_beads * HK_DIPLOID_N_COPY);
	memset(&old_empty_set, 0, sizeof(old_empty_set));
	memset(&new_empty_set, 0, sizeof(new_empty_set));
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		log_prior[i] = -logf((float)HK_BLIND_N_STATE);
	set_anchor_test_schedule_conf(&schedule_conf);
	ret = hk_blind_run_iter_loop_scheduled_cpu(fine, &old_empty_set, &conf, old_x,
											   log_prior, &schedule_conf, 0, &old_loop_diag);
	failed |= check_i32("scheduled old ret", ret, 0);
	ret = hk_blind_run_iter_loop_scheduled_parent_anchor_cpu(fine, &new_empty_set, &conf, new_x,
															log_prior, &schedule_conf,
															map, coarse_x, 0.0f,
															0, &new_loop_diag);
	failed |= check_i32("scheduled anchor zero ret", ret, 0);
	for (i = 0; i < fine->n_beads * HK_DIPLOID_N_COPY; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_close("scheduled anchor zero coords", new_x[i][a], old_x[i][a]);
	failed |= check_close("scheduled anchor zero energy", new_loop_diag.final_anchor_energy, 0.0f);
	failed |= check_close("scheduled anchor zero force", new_loop_diag.final_anchor_force_l1, 0.0f);
	failed |= check_i32("scheduled anchor zero nonfinite", new_loop_diag.n_anchor_nonfinite_step, 0);
	return failed;
}

static int check_multiresolution_pipeline(void)
{
	struct hk_sdict *dict = make_dict();
	struct hk_bmap *fine = 0, *coarse = 0;
	struct hk_blind_coarse_to_fine_map *map = 0;
	struct hk_blind_pair raw[32];
	struct hk_blind_bpair_set *set = 0;
	struct hk_blind_bpair *coarse_bp;
	int32_t coarse_bpair_id;
	int failed = 0;

	fine = hk_bmap_gen(dict, 0, 0, 1000000, 1);
	coarse = hk_bmap_gen(dict, 0, 0, 4000000, 1);
	failed |= check_true("fine bmap exists", fine != 0);
	failed |= check_true("coarse bmap exists", coarse != 0);
	if (fine == 0 || coarse == 0) goto cleanup;
	failed |= check_i32("fine 1Mb bead count", fine->n_beads, 10);
	failed |= check_i32("coarse 4Mb bead count", coarse->n_beads, 3);

	map = hk_blind_coarse_to_fine_map_build(coarse, fine);
	failed |= check_true("coarse-to-fine map exists", map != 0);
	if (map == 0) goto cleanup;
	failed |= check_i32("parent0 children", map->parent_child_count[0], 4);
	failed |= check_i32("parent1 children", map->parent_child_count[1], 4);
	failed |= check_i32("parent2 children", map->parent_child_count[2], 2);

	fill_raw_pairs_from_phased_sources(raw, (int32_t)(sizeof(raw) / sizeof(raw[0])));
	set = hk_blind_bpair_set_build(coarse, (int32_t)(sizeof(raw) / sizeof(raw[0])), raw);
	failed |= check_true("coarse bpair set exists", set != 0);
	if (set == 0) goto cleanup;
	failed |= check_i32("coarse bpair set n_raw", set->n_raw, 32);
	failed |= check_i32("coarse bpair set n_bpairs", set->n_bpairs, 1);

	coarse_bpair_id = find_bpair_id(set, 0, 1);
	failed |= check_true("coarse bpair 0-1 found", coarse_bpair_id >= 0);
	if (coarse_bpair_id < 0) goto cleanup;
	coarse_bp = &set->bpairs[coarse_bpair_id];
	failed |= check_i32("coarse bpair raw count", coarse_bp->n_raw, 32);
	failed |= check_close("coarse bpair raw base d_scale", coarse_bp->base_d_scale,
						  powf(32.0f, -1.0f / 3.0f));

	failed |= check_i32("apply density exposure",
						hk_blind_bpair_set_apply_density_exposure_from_map(set, map), 0);
	failed |= check_i32("coarse bpair child count0", coarse_bp->density_child_count[0], 4);
	failed |= check_i32("coarse bpair child count1", coarse_bp->density_child_count[1], 4);
	failed |= check_close("coarse bpair exposure", coarse_bp->density_exposure, 16.0f);
	failed |= check_close("density normalized raw count",
						  hk_blind_bpair_density_normalized_raw_count(coarse_bp, 1e-6f), 2.0f);
	failed |= check_density_normalized_edge_expansion(coarse_bp);
	failed |= check_i32("apply density d_scale",
						hk_blind_bpair_set_apply_density_d_scale_from_map(
							set, map, HK_BLIND_D_SCALE_DENSITY_NORMALIZED_RAW_COUNT, 1e-6f), 0);
	failed |= check_close("coarse bpair density base d_scale",
						  coarse_bp->base_d_scale, powf(2.0f, -1.0f / 3.0f));
	failed |= check_lift_and_anchor(coarse, fine, map);

cleanup:
	hk_blind_bpair_set_destroy(set);
	hk_blind_coarse_to_fine_map_destroy(map);
	if (fine) hk_bmap_destroy(fine);
	if (coarse) hk_bmap_destroy(coarse);
	destroy_dict(dict);
	return failed;
}

int main(void)
{
	int failed = 0;
	hk_verbose = 0;
	failed |= check_phase_labels_not_copied();
	failed |= check_multiresolution_pipeline();
	return failed? 1 : 0;
}
