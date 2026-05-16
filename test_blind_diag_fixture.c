#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "hickit.h"

static int check_i32(const char *label, int32_t got, int32_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %d, expected %d\n", label, got, expected);
		return 1;
	}
	return 0;
}

static int check_i64_positive(const char *label, int64_t got)
{
	if (got <= 0) {
		fprintf(stderr, "%s: got %lld, expected positive\n", label, (long long)got);
		return 1;
	}
	return 0;
}

static int check_close(const char *label, double got, double expected)
{
	double tol = 1e-6;
	double scale = fabs(expected) > 1.0? fabs(expected) : 1.0;
	if (fabs(got - expected) > tol * scale) {
		fprintf(stderr, "%s: got %.12g, expected %.12g\n", label, got, expected);
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

static int check_key(const char *label, const struct hk_blind_bpair_key *key, int32_t bid0, int32_t bid1)
{
	int failed = 0;
	char buf[96];
	snprintf(buf, sizeof(buf), "%s bid0", label);
	failed |= check_i32(buf, key->bid[0], bid0);
	snprintf(buf, sizeof(buf), "%s bid1", label);
	failed |= check_i32(buf, key->bid[1], bid1);
	return failed;
}

static int same_blind_pair(const struct hk_blind_pair *a, const struct hk_blind_pair *b)
{
	int failed = 0;
	failed |= check_i32("same chr0", a->chr[0], b->chr[0]);
	failed |= check_i32("same chr1", a->chr[1], b->chr[1]);
	failed |= check_i32("same pos0", a->pos[0], b->pos[0]);
	failed |= check_i32("same pos1", a->pos[1], b->pos[1]);
	failed |= check_i32("same strand0", a->strand[0], b->strand[0]);
	failed |= check_i32("same strand1", a->strand[1], b->strand[1]);
	return failed;
}

static int check_p4_close(const char *label, const float got[HK_BLIND_N_STATE], const float expected[HK_BLIND_N_STATE])
{
	int failed = 0, i;
	char buf[96];
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		snprintf(buf, sizeof(buf), "%s p4[%d]", label, i);
		failed |= check_close(buf, got[i], expected[i]);
	}
	return failed;
}

static void set_haploid_scaffold(const struct hk_bmap *bmap, fvec3_t *haploid)
{
	int32_t i;
	for (i = 0; i < bmap->n_beads; ++i) {
		haploid[i][0] = (float)i;
		haploid[i][1] = (float)(i % 3);
		haploid[i][2] = 0.25f * (float)(i % 5);
	}
}

static int check_blind_adapter_groups(const struct hk_map *m, struct hk_blind_pair *raw)
{
	struct hk_blind_pair forward0, reverse0;
	int have_forward = 0, have_reverse = 0, n_forward = 0, n_reverse = 0, n_same = 0;
	int failed = 0;
	int32_t i;

	for (i = 0; i < m->n_pairs; ++i) {
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);
		if (raw[i].chr[0] == 0 && raw[i].chr[1] == 0 && raw[i].pos[0] == 2000312 && raw[i].pos[1] == 5000312) {
			if (!have_forward) {
				forward0 = raw[i];
				have_forward = 1;
			} else {
				failed |= same_blind_pair(&forward0, &raw[i]);
			}
			++n_forward;
		} else if (raw[i].chr[0] == 0 && raw[i].chr[1] == 0 && raw[i].pos[0] == 5000312 && raw[i].pos[1] == 2000312) {
			if (!have_reverse) {
				reverse0 = raw[i];
				have_reverse = 1;
			} else {
				failed |= same_blind_pair(&reverse0, &raw[i]);
			}
			++n_reverse;
		} else if (raw[i].chr[0] == 0 && raw[i].chr[1] == 0 && raw[i].pos[0] == 3000312 && raw[i].pos[1] == 3008491) {
			++n_same;
		} else {
			fprintf(stderr, "unexpected blind pair at row %d: %d:%d -> %d:%d\n",
					(int)i, raw[i].chr[0], raw[i].pos[0], raw[i].chr[1], raw[i].pos[1]);
			failed = 1;
		}
	}
	failed |= check_i32("forward blind variants", n_forward, 4);
	failed |= check_i32("reverse blind variants", n_reverse, 4);
	failed |= check_i32("same-bin blind rows", n_same, 1);
	return failed;
}

static int check_binned_mapping(const struct hk_blind_bpair_set *set)
{
	int failed = 0;
	int32_t i;
	int n_key25 = 0, n_key33 = 0, n_25_swapped0 = 0, n_25_swapped1 = 0, n_33_swapped0 = 0;

	failed |= check_i32("set n_raw", set->n_raw, 9);
	failed |= check_i32("set n_bpairs", set->n_bpairs, 2);
	failed |= check_key("set key0", &set->bpairs[0].key, 2, 5);
	failed |= check_key("set key1", &set->bpairs[1].key, 3, 3);
	failed |= check_i32("set key0 n_raw", set->bpairs[0].n_raw, 8);
	failed |= check_i32("set key1 n_raw", set->bpairs[1].n_raw, 1);

	for (i = 0; i < set->n_raw; ++i) {
		int32_t bpair_id = set->raw2binned[i].bpair_id;
		if (bpair_id == 0) {
			++n_key25;
			if (set->raw2binned[i].swapped)
				++n_25_swapped1;
			else
				++n_25_swapped0;
		} else if (bpair_id == 1) {
			++n_key33;
			if (!set->raw2binned[i].swapped)
				++n_33_swapped0;
		} else {
			fprintf(stderr, "unexpected bpair id at raw %d: %d\n", (int)i, bpair_id);
			failed = 1;
		}
	}
	failed |= check_i32("key25 raw count", n_key25, 8);
	failed |= check_i32("key25 swapped false", n_25_swapped0, 4);
	failed |= check_i32("key25 swapped true", n_25_swapped1, 4);
	failed |= check_i32("key33 raw count", n_key33, 1);
	failed |= check_i32("key33 swapped false", n_33_swapped0, 1);
	return failed;
}

static int check_diag(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set, const struct hk_blind_iter_diag *diag)
{
	int failed = 0;
	failed |= check_i32("diag n_raw", diag->n_raw, set->n_raw);
	failed |= check_i32("diag n_bpair", diag->n_bpair, set->n_bpairs);
	failed |= check_i32("diag gauge flips", diag->gauge_stats.n_flipped, 0);
	failed |= check_true("diag mean entropy finite", isfinite(diag->mean_entropy));
	failed |= check_true("diag mean entropy range", diag->mean_entropy >= -1e-6f && diag->mean_entropy <= logf((float)HK_BLIND_N_STATE) + 1e-6f);
	failed |= check_true("diag mean pmax finite", isfinite(diag->mean_pmax));
	failed |= check_true("diag mean pmax range", diag->mean_pmax >= 0.25f - 1e-6f && diag->mean_pmax <= 1.0f + 1e-6f);
	failed |= check_true("diag mean pU finite", isfinite(diag->mean_pU));
	failed |= check_true("diag mean pU range", diag->mean_pU >= -1e-6f && diag->mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("diag mean rho finite", isfinite(diag->mean_rho_output));
	failed |= check_true("diag mean rho range", diag->mean_rho_output >= -1e-6f && diag->mean_rho_output <= 1.0f + 1e-6f);
	failed |= check_i32("diag sep finite", diag->sep_stats.n_finite, bmap->n_beads);
	failed |= check_i32("diag sep nonfinite", diag->sep_stats.n_nonfinite, 0);
	failed |= check_true("diag sep energy finite", isfinite(diag->sep_energy));
	failed |= check_i32("diag sep force nonfinite", diag->sep_force_nonfinite, 0);
	failed |= check_true("diag wedge count", diag->n_wedges >= 0);
	failed |= check_true("diag sum wedge finite", isfinite(diag->sum_wedge_k));
	failed |= check_true("diag sum wedge nonnegative", diag->sum_wedge_k >= 0.0);
	failed |= check_i32("diag wedge nonfinite", diag->n_wedge_nonfinite, 0);
	failed |= check_i32("diag posterior nonfinite", diag->n_posterior_nonfinite, 0);
	failed |= check_i32("diag posterior bad sum", diag->n_posterior_bad_sum, 0);
	failed |= check_i32("diag posterior out of range", diag->n_posterior_out_of_range, 0);
	failed |= check_i32("diag uncertainty nonfinite", diag->n_uncertainty_nonfinite, 0);
	failed |= check_i32("diag uncertainty out of range", diag->n_uncertainty_out_of_range, 0);
	failed |= check_i32("diag five-state bad sum", diag->n_five_state_bad_sum, 0);
	failed |= check_i32("diag wedge bad k", diag->n_wedge_bad_k, 0);
	failed |= check_i32("diag wedge bad d_scale", diag->n_wedge_bad_d_scale, 0);
	failed |= check_i32("diag skipped self edges", (int32_t)diag->n_skipped_self_edges, 0);
	failed |= check_i64_positive("diag skipped same-bin bpairs", diag->n_skipped_same_bin_bpairs);
	return failed;
}

static int check_raw_posterior_inheritance(const struct hk_blind_bpair_set *set)
{
	const float *canonical = set->bpairs[0].p4;
	float raw_ordered[HK_BLIND_N_STATE], raw_swapped[HK_BLIND_N_STATE];
	float expected_swapped[HK_BLIND_N_STATE];
	int have_ordered = 0, have_swapped = 0;
	int failed = 0;
	int32_t i;

	for (i = 0; i < set->n_raw; ++i) {
		if (set->raw2binned[i].bpair_id != 0) continue;
		if (set->raw2binned[i].swapped && !have_swapped) {
			hk_blind_p4_to_raw_order(canonical, set->raw2binned[i].swapped, raw_swapped);
			have_swapped = 1;
		} else if (!set->raw2binned[i].swapped && !have_ordered) {
			hk_blind_p4_to_raw_order(canonical, set->raw2binned[i].swapped, raw_ordered);
			have_ordered = 1;
		}
	}
	expected_swapped[HK_BLIND_STATE_00] = canonical[HK_BLIND_STATE_00];
	expected_swapped[HK_BLIND_STATE_01] = canonical[HK_BLIND_STATE_10];
	expected_swapped[HK_BLIND_STATE_10] = canonical[HK_BLIND_STATE_01];
	expected_swapped[HK_BLIND_STATE_11] = canonical[HK_BLIND_STATE_11];
	failed |= check_true("raw inheritance ordered present", have_ordered);
	failed |= check_true("raw inheritance swapped present", have_swapped);
	if (have_ordered)
		failed |= check_p4_close("raw ordered", raw_ordered, canonical);
	if (have_swapped)
		failed |= check_p4_close("raw swapped", raw_swapped, expected_swapped);
	return failed;
}

int main(void)
{
	struct hk_map *m;
	struct hk_bmap *bmap;
	struct hk_blind_pair raw[9];
	struct hk_blind_bpair_set *set;
	struct hk_fdg_conf conf;
	struct hk_blind_iter_diag diag;
	fvec3_t *haploid = 0, *diploid = 0;
	int failed = 0;

	hk_verbose = 0;
	m = hk_map_read("testdata/p9016_blind_diag.pairs");
	if (m == 0) {
		fprintf(stderr, "failed to read diag fixture\n");
		return 1;
	}
	failed |= check_i32("fixture n_pairs", m->n_pairs, 9);
	failed |= check_blind_adapter_groups(m, raw);

	bmap = hk_bmap_gen(m->d, m->n_pairs, m->pairs, 1000000, 1);
	if (bmap == 0) {
		fprintf(stderr, "failed to build bmap\n");
		hk_map_destroy(m);
		return 1;
	}

	set = hk_blind_bpair_set_build(bmap, m->n_pairs, raw);
	failed |= check_binned_mapping(set);

	haploid = (fvec3_t*)calloc(bmap->n_beads, sizeof(fvec3_t));
	diploid = (fvec3_t*)calloc((size_t)bmap->n_beads * HK_DIPLOID_N_COPY, sizeof(fvec3_t));
	if (haploid == 0 || diploid == 0) {
		fprintf(stderr, "failed to allocate coords\n");
		free(haploid);
		free(diploid);
		hk_blind_bpair_set_destroy(set);
		hk_bmap_destroy(bmap);
		hk_map_destroy(m);
		return 1;
	}
	set_haploid_scaffold(bmap, haploid);
	failed |= check_i32("init diploid ret", hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads, diploid,
																					 0.5f, 0.0f, 17), 0);

	hk_fdg_conf_init(&conf);
	failed |= check_i32("diag ret", hk_blind_run_single_iter_diag(bmap, set, &conf, 0, diploid, 1.0f, 1.0f, 2.0f,
																  0, 1.0f, 1.0f, 0.25f, 0.05f, &diag), 0);
	failed |= check_diag(bmap, set, &diag);
	failed |= check_raw_posterior_inheritance(set);

	free(haploid);
	free(diploid);
	hk_blind_bpair_set_destroy(set);
	hk_bmap_destroy(bmap);
	hk_map_destroy(m);
	return failed != 0;
}
