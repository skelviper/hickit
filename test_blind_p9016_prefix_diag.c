#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hickit.h"

#define HK_BLIND_P9016_PREFIX_N 1000
#define HK_BLIND_P9016_PATH "../pairs/P9016.pairs.gz"

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

static int check_contains(const char *label, const char *s, const char *needle)
{
	if (strstr(s, needle) == 0) {
		fprintf(stderr, "%s: missing substring '%s'\n", label, needle);
		return 1;
	}
	return 0;
}

static int file_exists(const char *path)
{
	FILE *fp = fopen(path, "rb");
	if (fp == 0) return 0;
	fclose(fp);
	return 1;
}

static void set_haploid_scaffold(const struct hk_bmap *bmap, fvec3_t *haploid)
{
	int32_t i;
	for (i = 0; i < bmap->n_beads; ++i) {
		haploid[i][0] = 0.10f * (float)(i % 97);
		haploid[i][1] = 0.07f * (float)((i / 97) % 97);
		haploid[i][2] = 0.03f * (float)(i % 17);
	}
}

static int check_binned_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set, int32_t max_raw)
{
	int failed = 0;
	int32_t i;

	failed |= check_true("prefix set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_true("prefix raw positive", set->n_raw > 0);
	failed |= check_true("prefix raw bounded", set->n_raw <= max_raw);
	failed |= check_true("prefix bpair positive", set->n_bpairs > 0);
	for (i = 0; i < set->n_raw; ++i) {
		failed |= check_true("prefix raw2binned id lower", set->raw2binned[i].bpair_id >= 0);
		failed |= check_true("prefix raw2binned id upper", set->raw2binned[i].bpair_id < set->n_bpairs);
		failed |= check_true("prefix raw2binned swapped", set->raw2binned[i].swapped == 0 || set->raw2binned[i].swapped == 1);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		failed |= check_true("prefix key sorted", bp->key.bid[0] <= bp->key.bid[1]);
		failed |= check_true("prefix key bid0 valid", bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		failed |= check_true("prefix key bid1 valid", bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
	}
	return failed;
}

static int check_diag_health(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
							 const struct hk_blind_iter_diag *diag)
{
	int failed = 0;
	failed |= check_i32("diag n_raw", diag->n_raw, set->n_raw);
	failed |= check_i32("diag n_bpair", diag->n_bpair, set->n_bpairs);
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
	failed |= check_i32("diag posterior nonfinite", diag->n_posterior_nonfinite, 0);
	failed |= check_i32("diag posterior bad sum", diag->n_posterior_bad_sum, 0);
	failed |= check_i32("diag posterior out of range", diag->n_posterior_out_of_range, 0);
	failed |= check_i32("diag uncertainty nonfinite", diag->n_uncertainty_nonfinite, 0);
	failed |= check_i32("diag uncertainty out of range", diag->n_uncertainty_out_of_range, 0);
	failed |= check_i32("diag five-state bad sum", diag->n_five_state_bad_sum, 0);
	failed |= check_i32("diag wedge nonfinite", diag->n_wedge_nonfinite, 0);
	failed |= check_i32("diag wedge bad k", diag->n_wedge_bad_k, 0);
	failed |= check_i32("diag wedge bad d_scale", diag->n_wedge_bad_d_scale, 0);
	return failed;
}

static int check_summary(const struct hk_blind_iter_diag *diag)
{
	char buf[4096];
	int failed = 0;
	int n = hk_blind_iter_diag_snprintf(buf, sizeof(buf), diag);
	failed |= check_true("summary length", n > 0 && n < (int)sizeof(buf));
	failed |= check_contains("summary status", buf, "status: OK");
	failed |= check_contains("summary n_raw", buf, "n_raw");
	failed |= check_contains("summary n_bpair", buf, "n_bpair");
	failed |= check_contains("summary mean_entropy", buf, "mean_entropy");
	failed |= check_contains("summary mean_pmax", buf, "mean_pmax");
	failed |= check_contains("summary posterior bad sum", buf, "n_posterior_bad_sum");
	failed |= check_contains("summary five-state", buf, "n_five_state_bad_sum");
	failed |= check_contains("summary wedges", buf, "n_wedges");
	failed |= check_contains("summary skipped self", buf, "n_skipped_self_edges");
	failed |= check_contains("summary sep energy", buf, "sep_energy");
	failed |= check_contains("summary chr flipped", buf, "n_chr_flipped");
	if (!failed)
		fprintf(stderr, "%s", buf);
	return failed;
}

int main(void)
{
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	struct hk_blind_bpair_set *set = 0;
	struct hk_fdg_conf conf;
	struct hk_blind_iter_diag diag;
	fvec3_t *haploid = 0, *diploid = 0;
	int32_t n_prefix;
	int32_t i;
	int failed = 0;

	if (!file_exists(HK_BLIND_P9016_PATH)) {
		fprintf(stderr, "SKIP: %s not found\n", HK_BLIND_P9016_PATH);
		return 0;
	}

	hk_verbose = 0;
	m = hk_map_read(HK_BLIND_P9016_PATH);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", HK_BLIND_P9016_PATH);
		return 1;
	}
	n_prefix = m->n_pairs < HK_BLIND_P9016_PREFIX_N? m->n_pairs : HK_BLIND_P9016_PREFIX_N;
	if (n_prefix <= 0) {
		fprintf(stderr, "no pairs in %s\n", HK_BLIND_P9016_PATH);
		hk_map_destroy(m);
		return 1;
	}

	bmap = hk_bmap_gen(m->d, n_prefix, m->pairs, 1000000, 1);
	raw = (struct hk_blind_pair*)calloc(n_prefix, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate prefix bmap/raw data\n");
		free(raw);
		if (bmap) hk_bmap_destroy(bmap);
		hk_map_destroy(m);
		return 1;
	}
	for (i = 0; i < n_prefix; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);

	set = hk_blind_bpair_set_build(bmap, n_prefix, raw);
	failed |= check_binned_set(bmap, set, HK_BLIND_P9016_PREFIX_N);

	haploid = (fvec3_t*)calloc(bmap->n_beads, sizeof(*haploid));
	diploid = (fvec3_t*)calloc((size_t)bmap->n_beads * HK_DIPLOID_N_COPY, sizeof(*diploid));
	if (haploid == 0 || diploid == 0) {
		fprintf(stderr, "failed to allocate prefix coords\n");
		failed = 1;
		goto cleanup;
	}
	set_haploid_scaffold(bmap, haploid);
	failed |= check_i32("prefix init diploid", hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
																						 diploid, 0.5f, 0.0f, 17), 0);

	hk_fdg_conf_init(&conf);
	failed |= check_i32("prefix diag ret", hk_blind_run_single_iter_diag(bmap, set, &conf, 0, diploid, 1.0f, 1.0f, 2.0f,
																		 0, 1.0f, 1.0f, 0.25f, 0.05f, &diag), 0);
	failed |= check_diag_health(bmap, set, &diag);
	failed |= check_summary(&diag);

cleanup:
	free(haploid);
	free(diploid);
	hk_blind_bpair_set_destroy(set);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	hk_map_destroy(m);
	return failed != 0;
}
