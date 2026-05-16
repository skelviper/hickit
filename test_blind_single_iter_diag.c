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

static int check_contains(const char *label, const char *s, const char *needle)
{
	if (strstr(s, needle) == 0) {
		fprintf(stderr, "%s: missing substring '%s'\n", label, needle);
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

static void set_pair(fvec3_t *coords, int32_t i, float x00, float x01, float x02, float x10, float x11, float x12)
{
	set_coord(coords[hk_diploid_bid(i, HK_DIPLOID_COPY0)], x00, x01, x02);
	set_coord(coords[hk_diploid_bid(i, HK_DIPLOID_COPY1)], x10, x11, x12);
}

static void copy_coords(fvec3_t *dst, const fvec3_t *src, int32_t n_diploid)
{
	int32_t i;
	int a;
	for (i = 0; i < n_diploid; ++i)
		for (a = 0; a < 3; ++a)
			dst[i][a] = src[i][a];
}

static int check_coords_close(const char *label, const fvec3_t *got, const fvec3_t *expected, int32_t n_diploid)
{
	int failed = 0;
	int32_t i;
	int a;
	char buf[96];
	for (i = 0; i < n_diploid; ++i) {
		for (a = 0; a < 3; ++a) {
			snprintf(buf, sizeof(buf), "%s bead%d axis%d", label, (int)i, a);
			failed |= check_close(buf, got[i][a], expected[i][a]);
		}
	}
	return failed;
}

static void set_test_bmap(struct hk_bmap *bmap, struct hk_bead beads[4])
{
	int32_t i;
	for (i = 0; i < 4; ++i) {
		beads[i].chr = i < 2? 0 : 1;
		beads[i].st = i * 1000000;
		beads[i].en = beads[i].st + 1000000;
	}
	bmap->n_beads = 4;
	bmap->n_pairs = 0;
	bmap->unit = 1.0f;
	bmap->d = 0;
	bmap->beads = beads;
	bmap->offcnt = 0;
	bmap->pairs = 0;
	bmap->x = 0;
	bmap->feat = 0;
	bmap->cpg = 0;
	bmap->gc_bias = 0;
	bmap->gc_corrected = 0;
}

static void set_test_coords(fvec3_t coords[8])
{
	set_pair(coords, 0, 0.0f, 0.0f, 0.0f, 5.0f, 0.0f, 0.0f);
	set_pair(coords, 1, 1.0f, 0.0f, 0.0f, 6.0f, 0.0f, 0.0f);
	set_pair(coords, 2, 10.0f, 0.0f, 0.0f, 15.0f, 0.0f, 0.0f);
	set_pair(coords, 3, 11.0f, 0.0f, 0.0f, 16.0f, 0.0f, 0.0f);
}

static void swap_pair_copies(fvec3_t *coords, int32_t i)
{
	int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
	int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
	float t;
	int a;
	for (a = 0; a < 3; ++a) {
		t = coords[b0][a];
		coords[b0][a] = coords[b1][a];
		coords[b1][a] = t;
	}
}

static void set_bpair(struct hk_blind_bpair *bp, int32_t bid0, int32_t bid1)
{
	int i;
	bp->key.bid[0] = bid0;
	bp->key.bid[1] = bid1;
	bp->n_raw = 1;
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		bp->p4[i] = 1.0f / HK_BLIND_N_STATE;
	bp->entropy = logf((float)HK_BLIND_N_STATE);
	bp->pmax = 1.0f / HK_BLIND_N_STATE;
	bp->margin = 0.0f;
	bp->rho_output = 0.0f;
	bp->pU = 1.0f;
}

static void set_test_bpair_set(struct hk_blind_bpair_set *set, struct hk_blind_bpair bpairs[3],
							   struct hk_blind_raw2binned raw2binned[3])
{
	set_bpair(&bpairs[0], 0, 1);
	set_bpair(&bpairs[1], 1, 2);
	set_bpair(&bpairs[2], 2, 2);
	hk_blind_raw2binned_set(&raw2binned[0], 0, 0);
	hk_blind_raw2binned_set(&raw2binned[1], 1, 0);
	hk_blind_raw2binned_set(&raw2binned[2], 2, 0);
	set->bpairs = bpairs;
	set->n_bpairs = 3;
	set->raw2binned = raw2binned;
	set->n_raw = 3;
}

static int check_diag_ranges(const char *label, const struct hk_blind_iter_diag *diag)
{
	int failed = 0;
	char buf[96];
	snprintf(buf, sizeof(buf), "%s mean entropy", label);
	failed |= check_true(buf, diag->mean_entropy >= -1e-6f && diag->mean_entropy <= logf((float)HK_BLIND_N_STATE) + 1e-6f);
	snprintf(buf, sizeof(buf), "%s mean pmax", label);
	failed |= check_true(buf, diag->mean_pmax >= 0.25f - 1e-6f && diag->mean_pmax <= 1.0f + 1e-6f);
	snprintf(buf, sizeof(buf), "%s mean pU", label);
	failed |= check_true(buf, diag->mean_pU >= -1e-6f && diag->mean_pU <= 1.0f + 1e-6f);
	snprintf(buf, sizeof(buf), "%s mean rho", label);
	failed |= check_true(buf, diag->mean_rho_output >= -1e-6f && diag->mean_rho_output <= 1.0f + 1e-6f);
	snprintf(buf, sizeof(buf), "%s mean margin", label);
	failed |= check_true(buf, diag->mean_margin >= -1e-6f && diag->mean_margin <= 1.0f + 1e-6f);
	return failed;
}

static int check_clean_validation_counters(const char *label, const struct hk_blind_iter_diag *diag)
{
	int failed = 0;
	char buf[128];
	snprintf(buf, sizeof(buf), "%s posterior nonfinite", label);
	failed |= check_i32(buf, diag->n_posterior_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s posterior bad sum", label);
	failed |= check_i32(buf, diag->n_posterior_bad_sum, 0);
	snprintf(buf, sizeof(buf), "%s posterior out of range", label);
	failed |= check_i32(buf, diag->n_posterior_out_of_range, 0);
	snprintf(buf, sizeof(buf), "%s uncertainty nonfinite", label);
	failed |= check_i32(buf, diag->n_uncertainty_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s uncertainty out of range", label);
	failed |= check_i32(buf, diag->n_uncertainty_out_of_range, 0);
	snprintf(buf, sizeof(buf), "%s five-state bad sum", label);
	failed |= check_i32(buf, diag->n_five_state_bad_sum, 0);
	snprintf(buf, sizeof(buf), "%s wedge bad k", label);
	failed |= check_i32(buf, diag->n_wedge_bad_k, 0);
	snprintf(buf, sizeof(buf), "%s wedge bad d_scale", label);
	failed |= check_i32(buf, diag->n_wedge_bad_d_scale, 0);
	return failed;
}

static int check_probability_sums(const struct hk_blind_bpair_set *set)
{
	int failed = 0;
	int32_t i;
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		double sum = 0.0;
		int s;
		char buf[96];
		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			snprintf(buf, sizeof(buf), "p4 range bpair%d state%d", (int)i, s);
			failed |= check_true(buf, bp->p4[s] >= -1e-6f && bp->p4[s] <= 1.0f + 1e-6f);
			sum += bp->p4[s];
		}
		snprintf(buf, sizeof(buf), "p4 sum bpair%d", (int)i);
		failed |= check_close(buf, sum, 1.0);
		snprintf(buf, sizeof(buf), "output sum bpair%d", (int)i);
		failed |= check_close(buf, bp->rho_output * sum + bp->pU, 1.0);
	}
	return failed;
}

static int check_basic_no_prev_diag(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	struct hk_blind_bpair bpairs[3];
	struct hk_blind_raw2binned raw2binned[3];
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf conf;
	struct hk_blind_iter_diag diag;
	fvec3_t coords[8], before[8];
	int failed = 0;

	set_test_bmap(&bmap, beads);
	set_test_bpair_set(&set, bpairs, raw2binned);
	set_test_coords(coords);
	copy_coords(before, coords, 8);
	hk_fdg_conf_init(&conf);

	failed |= check_i32("no-prev diag ret", hk_blind_run_single_iter_diag(&bmap, &set, &conf, 0, coords, 1.0f, 1.0f, 2.0f,
																		 0, 1.0f, 1.0f, 6.0f, 0.05f, &diag), 0);
	failed |= check_coords_close("no-prev coords unchanged", coords, before, 8);
	failed |= check_i32("no-prev n_bpair", diag.n_bpair, 3);
	failed |= check_i32("no-prev n_raw", diag.n_raw, 3);
	failed |= check_i32("no-prev gauge flips", diag.gauge_stats.n_flipped, 0);
	failed |= check_i32("no-prev n_chr_flipped", diag.n_chr_flipped, 0);
	failed |= check_i32("no-prev sep n_haploid", diag.sep_stats.n_haploid, 4);
	failed |= check_i32("no-prev sep nonfinite", diag.sep_stats.n_nonfinite, 0);
	failed |= check_i32("no-prev force nonfinite", diag.sep_force_nonfinite, 0);
	failed |= check_i32("no-prev wedge nonfinite", diag.n_wedge_nonfinite, 0);
	failed |= check_i64("no-prev self skipped", diag.n_skipped_self_edges, 0);
	failed |= check_i64("no-prev same-bin skipped", diag.n_skipped_same_bin_bpairs, 1);
	failed |= check_true("no-prev wedges before aggregation", diag.n_wedges_before_aggregation > 0);
	failed |= check_true("no-prev wedges after aggregation", diag.n_wedges > 0);
	failed |= check_true("no-prev sum wedge k", diag.sum_wedge_k >= 0.0);
	failed |= check_true("no-prev sep energy active", diag.sep_energy > 0.0f);
	failed |= check_true("no-prev sep force active", diag.sep_force_l1 > 0.0f);
	failed |= check_diag_ranges("no-prev", &diag);
	failed |= check_clean_validation_counters("no-prev", &diag);
	failed |= check_probability_sums(&set);
	return failed;
}

static int check_prev_gauge_diag(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	struct hk_blind_bpair bpairs[3];
	struct hk_blind_raw2binned raw2binned[3];
	struct hk_blind_bpair_set set;
	struct hk_fdg_conf conf;
	struct hk_blind_iter_diag diag;
	fvec3_t prev[8], cur[8];
	int failed = 0;

	set_test_bmap(&bmap, beads);
	set_test_bpair_set(&set, bpairs, raw2binned);
	set_test_coords(prev);
	copy_coords(cur, prev, 8);
	swap_pair_copies(cur, 0);
	swap_pair_copies(cur, 1);
	hk_fdg_conf_init(&conf);

	failed |= check_i32("prev diag ret", hk_blind_run_single_iter_diag(&bmap, &set, &conf, prev, cur, 1.0f, 1.0f, 2.0f,
																	   0, 1.0f, 1.0f, 6.0f, 0.05f, &diag), 0);
	failed |= check_i32("prev gauge flips", diag.gauge_stats.n_flipped, 1);
	failed |= check_i32("prev n_chr_flipped", diag.n_chr_flipped, 1);
	failed |= check_coords_close("prev coords stabilized", cur, prev, 8);
	failed |= check_diag_ranges("prev", &diag);
	failed |= check_clean_validation_counters("prev", &diag);
	failed |= check_probability_sums(&set);
	return failed;
}

static int check_rho_train_scales_wedges(void)
{
	struct hk_bmap bmap;
	struct hk_bead beads[4];
	struct hk_blind_bpair bpairs_full[3], bpairs_half[3];
	struct hk_blind_raw2binned raw_full[3], raw_half[3];
	struct hk_blind_bpair_set set_full, set_half;
	struct hk_fdg_conf conf;
	struct hk_blind_iter_diag diag_full, diag_half;
	fvec3_t coords_full[8], coords_half[8];
	int failed = 0;

	set_test_bmap(&bmap, beads);
	set_test_bpair_set(&set_full, bpairs_full, raw_full);
	set_test_bpair_set(&set_half, bpairs_half, raw_half);
	set_test_coords(coords_full);
	copy_coords(coords_half, coords_full, 8);
	hk_fdg_conf_init(&conf);

	failed |= check_i32("rho full ret", hk_blind_run_single_iter_diag(&bmap, &set_full, &conf, 0, coords_full, 1.0f, 1.0f, 2.0f,
																	  0, 1.0f, 1.0f, 6.0f, 0.05f, &diag_full), 0);
	failed |= check_i32("rho half ret", hk_blind_run_single_iter_diag(&bmap, &set_half, &conf, 0, coords_half, 1.0f, 1.0f, 2.0f,
																	  0, 1.0f, 0.5f, 6.0f, 0.05f, &diag_half), 0);
	failed |= check_true("rho full positive", diag_full.sum_wedge_k > 0.0);
	failed |= check_close("rho half sum_k", diag_half.sum_wedge_k, 0.5 * diag_full.sum_wedge_k);
	failed |= check_i64("rho self skipped full", diag_full.n_skipped_self_edges, 0);
	failed |= check_i64("rho self skipped half", diag_half.n_skipped_self_edges, 0);
	failed |= check_i64("rho same-bin skipped full", diag_full.n_skipped_same_bin_bpairs, 1);
	failed |= check_i64("rho same-bin skipped half", diag_half.n_skipped_same_bin_bpairs, 1);
	failed |= check_clean_validation_counters("rho full", &diag_full);
	failed |= check_clean_validation_counters("rho half", &diag_half);
	return failed;
}

static int check_manual_bad_posterior_validation(void)
{
	struct hk_blind_bpair bpairs[5];
	struct hk_blind_bpair_set set;
	struct hk_blind_iter_diag diag;
	int failed = 0;
	int i;

	for (i = 0; i < 5; ++i)
		set_bpair(&bpairs[i], 0, 1);

	bpairs[0].p4[HK_BLIND_STATE_00] = 0.5f;
	bpairs[0].p4[HK_BLIND_STATE_01] = 0.5f;
	bpairs[0].p4[HK_BLIND_STATE_10] = 0.5f;
	bpairs[0].p4[HK_BLIND_STATE_11] = 0.5f;
	bpairs[0].pmax = 0.5f;
	bpairs[0].rho_output = 1.0f;
	bpairs[0].pU = 0.0f;

	bpairs[1].p4[HK_BLIND_STATE_00] = 1.2f;
	bpairs[1].p4[HK_BLIND_STATE_01] = -0.1f;
	bpairs[1].p4[HK_BLIND_STATE_10] = 0.0f;
	bpairs[1].p4[HK_BLIND_STATE_11] = -0.1f;
	bpairs[1].entropy = 0.0f;
	bpairs[1].pmax = 1.0f;
	bpairs[1].rho_output = 1.0f;
	bpairs[1].pU = 0.0f;

	bpairs[2].p4[HK_BLIND_STATE_00] = NAN;

	bpairs[3].rho_output = 1.5f;
	bpairs[3].pU = 0.0f;

	bpairs[4].entropy = NAN;

	set.bpairs = bpairs;
	set.n_bpairs = 5;
	set.raw2binned = 0;
	set.n_raw = 0;

	hk_blind_iter_diag_init(&diag);
	hk_blind_iter_diag_validate_bpair_set(&set, &diag);

	failed |= check_i32("manual posterior nonfinite", diag.n_posterior_nonfinite, 1);
	failed |= check_i32("manual posterior bad sum", diag.n_posterior_bad_sum, 1);
	failed |= check_i32("manual posterior out of range", diag.n_posterior_out_of_range, 1);
	failed |= check_i32("manual uncertainty nonfinite", diag.n_uncertainty_nonfinite, 1);
	failed |= check_i32("manual uncertainty out of range", diag.n_uncertainty_out_of_range, 1);
	failed |= check_i32("manual five-state bad sum", diag.n_five_state_bad_sum, 2);
	return failed;
}

static void set_wedge(struct hk_blind_wedge *edge, float k, float d_scale)
{
	edge->bid[0] = 6;
	edge->bid[1] = 7;
	edge->k = k;
	edge->d_scale = d_scale;
	edge->state = HK_BLIND_STATE_01;
	edge->state_mask = (uint8_t)(1u << HK_BLIND_STATE_01);
}

static int check_manual_bad_wedge_validation(void)
{
	struct hk_blind_wedge edges[4];
	struct hk_blind_wedge_list list;
	struct hk_blind_iter_diag diag;
	int failed = 0;

	set_wedge(&edges[0], NAN, 1.0f);
	set_wedge(&edges[1], -1.0f, 1.0f);
	set_wedge(&edges[2], 0.5f, 0.0f);
	set_wedge(&edges[3], 0.5f, NAN);
	list.edges = edges;
	list.n_edges = 4;
	list.m_edges = 4;
	list.n_input_bpair = 0;
	list.n_expanded_edges = 0;
	list.n_skipped_self_edges = 0;
	list.n_edges_before_aggregation = 0;
	list.n_aggregated_edges_removed = 0;

	hk_blind_iter_diag_init(&diag);
	hk_blind_iter_diag_validate_wedge_list(&list, &diag);

	failed |= check_i32("manual wedge nonfinite", diag.n_wedge_nonfinite, 2);
	failed |= check_i32("manual wedge bad k", diag.n_wedge_bad_k, 2);
	failed |= check_i32("manual wedge bad d_scale", diag.n_wedge_bad_d_scale, 2);
	return failed;
}

static void set_summary_diag(struct hk_blind_iter_diag *diag)
{
	hk_blind_iter_diag_init(diag);
	diag->n_raw = 7;
	diag->n_bpair = 3;
	diag->mean_entropy = 1.0f;
	diag->mean_pmax = 0.5f;
	diag->mean_margin = 0.25f;
	diag->mean_pU = 0.75f;
	diag->mean_rho_output = 0.25f;
	diag->sep_stats.n_haploid = 4;
	diag->sep_stats.n_finite = 4;
	diag->sep_stats.n_nonfinite = 0;
	diag->sep_stats.n_collapsed = 1;
	diag->sep_stats.min_sep = 0.1f;
	diag->sep_stats.mean_sep = 0.4f;
	diag->sep_stats.max_sep = 0.7f;
	diag->sep_stats.collapse_threshold = 0.5f;
	diag->sep_energy = 0.01f;
	diag->sep_force_l1 = 0.2f;
	diag->n_wedges = 11;
	diag->n_wedges_before_aggregation = 12;
	diag->n_skipped_self_edges = 2;
	diag->sum_wedge_k = 4.5;
	diag->gauge_stats.n_chr = 2;
	diag->gauge_stats.n_flipped = 1;
	diag->gauge_stats.mismatch_no_flip = 3.0;
	diag->gauge_stats.mismatch_flip = 1.0;
	diag->gauge_stats.mismatch_chosen = 1.0;
	diag->n_chr_flipped = 1;
}

static int check_diag_summary_helper(void)
{
	struct hk_blind_iter_diag diag;
	char buf[2048];
	char tiny[8];
	int n;
	int failed = 0;

	set_summary_diag(&diag);
	n = hk_blind_iter_diag_snprintf(buf, sizeof(buf), &diag);
	failed |= check_true("summary clean length", n > 0 && n < (int)sizeof(buf));
	failed |= check_contains("summary clean status", buf, "status: OK");
	failed |= check_contains("summary n_raw", buf, "n_raw=7");
	failed |= check_contains("summary n_bpair", buf, "n_bpair=3");
	failed |= check_contains("summary mean_entropy", buf, "mean_entropy=");
	failed |= check_contains("summary mean_pmax", buf, "mean_pmax=");
	failed |= check_contains("summary posterior bad sum key", buf, "n_posterior_bad_sum=0");
	failed |= check_contains("summary five-state key", buf, "n_five_state_bad_sum=0");
	failed |= check_contains("summary n_wedges", buf, "n_wedges=11");
	failed |= check_contains("summary self edges", buf, "n_skipped_self_edges=2");
	failed |= check_contains("summary sep energy", buf, "sep_energy=");
	failed |= check_contains("summary chr flipped", buf, "n_chr_flipped=1");

	diag.n_posterior_nonfinite = 1;
	diag.n_five_state_bad_sum = 2;
	diag.n_wedge_bad_k = 3;
	n = hk_blind_iter_diag_snprintf(buf, sizeof(buf), &diag);
	failed |= check_true("summary bad length", n > 0 && n < (int)sizeof(buf));
	failed |= check_contains("summary warn status", buf, "status: WARN");
	failed |= check_contains("summary posterior nonfinite", buf, "posterior_nonfinite=1");
	failed |= check_contains("summary five-state bad", buf, "five_state_bad_sum=2");
	failed |= check_contains("summary wedge bad k", buf, "wedge_bad_k=3");

	memset(tiny, 'X', sizeof(tiny));
	n = hk_blind_iter_diag_snprintf(tiny, sizeof(tiny), &diag);
	failed |= check_true("summary trunc length", n >= (int)sizeof(tiny));
	failed |= check_true("summary trunc nul", tiny[sizeof(tiny) - 1] == '\0');

	n = hk_blind_iter_diag_snprintf(0, 0, &diag);
	failed |= check_true("summary null zero length", n > 0);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_basic_no_prev_diag();
	failed |= check_prev_gauge_diag();
	failed |= check_rho_train_scales_wedges();
	failed |= check_manual_bad_posterior_validation();
	failed |= check_manual_bad_wedge_validation();
	failed |= check_diag_summary_helper();
	return failed != 0;
}
