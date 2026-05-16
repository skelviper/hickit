#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hickit.h"

#define HK_BLIND_P9016_DRY_RELAX_PREFIX_N 1000
#define HK_BLIND_P9016_DRY_RELAX_PATH "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_DRY_RELAX_UNIT 1.0f
#define HK_BLIND_P9016_DRY_RELAX_D_SCALE 1.0f
#define HK_BLIND_P9016_DRY_RELAX_BASE_K 2.0f
#define HK_BLIND_P9016_DRY_RELAX_TEMPERATURE 1.0f
#define HK_BLIND_P9016_DRY_RELAX_RHO_TRAIN 1.0f
#define HK_BLIND_P9016_DRY_RELAX_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_DRY_RELAX_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_DRY_RELAX_STEP 0.001f
#define HK_BLIND_P9016_DRY_RELAX_N_STEPS 5

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
	float tol = 1e-4f;
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

static void copy_coords(fvec3_t *dst, const fvec3_t *src, int32_t n)
{
	int32_t i;
	for (i = 0; i < n; ++i) {
		dst[i][0] = src[i][0];
		dst[i][1] = src[i][1];
		dst[i][2] = src[i][2];
	}
}

static int coords_any_changed(const fvec3_t *before, const fvec3_t *after, int32_t n, float tol)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			if (fabsf(after[i][a] - before[i][a]) > tol)
				return 1;
	return 0;
}

static int check_coords_finite(const fvec3_t *coords, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("dry-relax coord finite", isfinite(coords[i][a]));
	return failed;
}

static int check_binned_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set, int32_t max_raw)
{
	int failed = 0;
	int32_t i;

	failed |= check_true("dry-relax set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_true("dry-relax raw positive", set->n_raw > 0);
	failed |= check_true("dry-relax raw bounded", set->n_raw <= max_raw);
	failed |= check_true("dry-relax bpair positive", set->n_bpairs > 0);
	for (i = 0; i < set->n_raw; ++i) {
		failed |= check_true("dry-relax raw2binned id lower", set->raw2binned[i].bpair_id >= 0);
		failed |= check_true("dry-relax raw2binned id upper", set->raw2binned[i].bpair_id < set->n_bpairs);
		failed |= check_true("dry-relax raw2binned swapped", set->raw2binned[i].swapped == 0 || set->raw2binned[i].swapped == 1);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		failed |= check_true("dry-relax key sorted", bp->key.bid[0] <= bp->key.bid[1]);
		failed |= check_true("dry-relax key bid0 valid", bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		failed |= check_true("dry-relax key bid1 valid", bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
	}
	return failed;
}

static int check_posterior_health(const struct hk_blind_bpair_set *set, float *mean_entropy, float *mean_pmax)
{
	struct hk_blind_iter_diag diag;
	double entropy = 0.0;
	double pmax = 0.0;
	int failed = 0;
	int32_t i;

	hk_blind_iter_diag_init(&diag);
	hk_blind_iter_diag_validate_bpair_set(set, &diag);
	failed |= check_i32("dry-relax posterior nonfinite", diag.n_posterior_nonfinite, 0);
	failed |= check_i32("dry-relax posterior bad sum", diag.n_posterior_bad_sum, 0);
	failed |= check_i32("dry-relax posterior out of range", diag.n_posterior_out_of_range, 0);
	failed |= check_i32("dry-relax uncertainty nonfinite", diag.n_uncertainty_nonfinite, 0);
	failed |= check_i32("dry-relax uncertainty out of range", diag.n_uncertainty_out_of_range, 0);
	failed |= check_i32("dry-relax five-state bad sum", diag.n_five_state_bad_sum, 0);

	for (i = 0; i < set->n_bpairs; ++i) {
		entropy += set->bpairs[i].entropy;
		pmax += set->bpairs[i].pmax;
	}
	if (set->n_bpairs > 0) {
		*mean_entropy = (float)(entropy / set->n_bpairs);
		*mean_pmax = (float)(pmax / set->n_bpairs);
	}
	failed |= check_true("dry-relax mean entropy finite", isfinite(*mean_entropy));
	failed |= check_true("dry-relax mean entropy range", *mean_entropy >= -1e-6f && *mean_entropy <= logf((float)HK_BLIND_N_STATE) + 1e-6f);
	failed |= check_true("dry-relax mean pmax finite", isfinite(*mean_pmax));
	failed |= check_true("dry-relax mean pmax range", *mean_pmax >= 0.25f - 1e-6f && *mean_pmax <= 1.0f + 1e-6f);
	return failed;
}

static double sum_wedge_k(const struct hk_blind_wedge_list *list)
{
	double sum = 0.0;
	int32_t i;
	for (i = 0; i < list->n_edges; ++i)
		sum += list->edges[i].k;
	return sum;
}

static int check_wedge_health(const struct hk_blind_wedge_list *list, int32_t n_diploid, int32_t n_bpair)
{
	struct hk_blind_iter_diag diag;
	double sum_k = sum_wedge_k(list);
	int failed = 0;
	int32_t i;

	hk_blind_iter_diag_init(&diag);
	hk_blind_iter_diag_validate_wedge_list(list, &diag);
	failed |= check_i32("dry-relax wedge nonfinite", diag.n_wedge_nonfinite, 0);
	failed |= check_i32("dry-relax wedge bad k", diag.n_wedge_bad_k, 0);
	failed |= check_i32("dry-relax wedge bad d_scale", diag.n_wedge_bad_d_scale, 0);
	failed |= check_i64("dry-relax wedge input bpair", list->n_input_bpair, n_bpair);
	failed |= check_i64("dry-relax wedge expanded", list->n_expanded_edges, (int64_t)n_bpair * HK_BLIND_N_STATE);
	failed |= check_true("dry-relax wedge aggregate count", list->n_edges <= list->n_edges_before_aggregation);
	failed |= check_true("dry-relax skipped self nonnegative", list->n_skipped_self_edges >= 0);
	failed |= check_true("dry-relax wedge count", list->n_edges >= 0);
	failed |= check_true("dry-relax sum k finite", isfinite(sum_k));
	failed |= check_true("dry-relax sum k nonnegative", sum_k >= 0.0);
	for (i = 0; i < list->n_edges; ++i) {
		const struct hk_blind_wedge *e = &list->edges[i];
		failed |= check_true("dry-relax wedge sorted", e->bid[0] <= e->bid[1]);
		failed |= check_true("dry-relax wedge not self", e->bid[0] != e->bid[1]);
		failed |= check_true("dry-relax wedge bid0 valid", e->bid[0] >= 0 && e->bid[0] < n_diploid);
		failed |= check_true("dry-relax wedge bid1 valid", e->bid[1] >= 0 && e->bid[1] < n_diploid);
	}
	return failed;
}

static int check_sep_stats(const char *label, const struct hk_bmap *bmap, const struct hk_blind_homolog_sep_stats *stats)
{
	int failed = 0;
	char buf[128];
	snprintf(buf, sizeof(buf), "%s finite count", label);
	failed |= check_i32(buf, stats->n_finite, bmap->n_beads);
	snprintf(buf, sizeof(buf), "%s nonfinite count", label);
	failed |= check_i32(buf, stats->n_nonfinite, 0);
	snprintf(buf, sizeof(buf), "%s min finite", label);
	failed |= check_true(buf, isfinite(stats->min_sep));
	snprintf(buf, sizeof(buf), "%s mean finite", label);
	failed |= check_true(buf, isfinite(stats->mean_sep));
	snprintf(buf, sizeof(buf), "%s max finite", label);
	failed |= check_true(buf, isfinite(stats->max_sep));
	return failed;
}

static int check_relax_diag(const struct hk_blind_relax_diag *diag)
{
	int failed = 0;
	failed |= check_i32("dry-relax completed", diag->n_completed, HK_BLIND_P9016_DRY_RELAX_N_STEPS);
	failed |= check_i32("dry-relax nonfinite step", diag->n_nonfinite_step, 0);
	failed |= check_i32("dry-relax coord nonfinite", diag->n_coord_nonfinite, 0);
	failed |= check_true("dry-relax initial contact finite", isfinite(diag->initial_contact_energy));
	failed |= check_true("dry-relax final contact finite", isfinite(diag->final_contact_energy));
	failed |= check_true("dry-relax initial backbone finite", isfinite(diag->initial_backbone_energy));
	failed |= check_true("dry-relax final backbone finite", isfinite(diag->final_backbone_energy));
	failed |= check_true("dry-relax initial sep finite", isfinite(diag->initial_sep_energy));
	failed |= check_true("dry-relax final sep finite", isfinite(diag->final_sep_energy));
	failed |= check_true("dry-relax initial total finite", isfinite(diag->initial_total_energy));
	failed |= check_true("dry-relax final total finite", isfinite(diag->final_total_energy));
	failed |= check_close("dry-relax initial total sum",
						  diag->initial_total_energy,
						  diag->initial_contact_energy + diag->initial_backbone_energy + diag->initial_sep_energy);
	failed |= check_close("dry-relax final total sum",
						  diag->final_total_energy,
						  diag->final_contact_energy + diag->final_backbone_energy + diag->final_sep_energy);
	failed |= check_true("dry-relax max force finite", isfinite(diag->max_force_l1));
	failed |= check_true("dry-relax final force finite", isfinite(diag->final_force_l1));
	failed |= check_true("dry-relax max backbone force finite", isfinite(diag->max_backbone_force_l1));
	failed |= check_true("dry-relax final backbone force finite", isfinite(diag->final_backbone_force_l1));
	failed |= check_i32("dry-relax backbone nonfinite step", diag->n_backbone_nonfinite_step, 0);
	return failed;
}

int main(void)
{
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	struct hk_blind_bpair_set *set = 0;
	struct hk_blind_wedge_list wedges;
	struct hk_fdg_conf conf;
	struct hk_blind_relax_diag relax_diag;
	struct hk_blind_homolog_sep_stats sep_before, sep_after;
	fvec3_t *haploid = 0, *diploid = 0, *diploid_before = 0;
	float log_prior[HK_BLIND_N_STATE];
	float mean_entropy = 0.0f, mean_pmax = 0.0f;
	double wedge_sum_k = 0.0, wedge_sum_k_before_relax = 0.0;
	int32_t n_edges_before_relax = 0;
	int32_t n_prefix, n_diploid;
	int32_t i;
	int moved = 0;
	int failed = 0;

	if (!file_exists(HK_BLIND_P9016_DRY_RELAX_PATH)) {
		fprintf(stderr, "SKIP: %s not found\n", HK_BLIND_P9016_DRY_RELAX_PATH);
		return 0;
	}

	hk_verbose = 0;
	hk_blind_wedge_list_init(&wedges);
	m = hk_map_read(HK_BLIND_P9016_DRY_RELAX_PATH);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", HK_BLIND_P9016_DRY_RELAX_PATH);
		return 1;
	}
	n_prefix = m->n_pairs < HK_BLIND_P9016_DRY_RELAX_PREFIX_N? m->n_pairs : HK_BLIND_P9016_DRY_RELAX_PREFIX_N;
	if (n_prefix <= 0) {
		fprintf(stderr, "no pairs in %s\n", HK_BLIND_P9016_DRY_RELAX_PATH);
		failed = 1;
		goto cleanup;
	}

	bmap = hk_bmap_gen(m->d, n_prefix, m->pairs, 1000000, 1);
	raw = (struct hk_blind_pair*)calloc(n_prefix, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate P9016 prefix bmap/raw data\n");
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n_prefix; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);

	set = hk_blind_bpair_set_build(bmap, n_prefix, raw);
	failed |= check_binned_set(bmap, set, HK_BLIND_P9016_DRY_RELAX_PREFIX_N);
	if (failed) goto cleanup;

	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	haploid = (fvec3_t*)calloc(bmap->n_beads, sizeof(*haploid));
	diploid = (fvec3_t*)calloc(n_diploid, sizeof(*diploid));
	diploid_before = (fvec3_t*)calloc(n_diploid, sizeof(*diploid_before));
	if (haploid == 0 || diploid == 0 || diploid_before == 0) {
		fprintf(stderr, "failed to allocate P9016 prefix coords\n");
		failed = 1;
		goto cleanup;
	}
	set_haploid_scaffold(bmap, haploid);
	failed |= check_i32("dry-relax init diploid",
						hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads, diploid, 0.5f, 0.0f, 17), 0);
	failed |= check_coords_finite(diploid, n_diploid);

	hk_fdg_conf_init(&conf);
	hk_blind_init_uniform_log_prior(log_prior);
	hk_blind_bpair_set_update_posterior_from_coords(set, &conf, diploid,
													 HK_BLIND_P9016_DRY_RELAX_UNIT,
													 HK_BLIND_P9016_DRY_RELAX_D_SCALE,
													 HK_BLIND_P9016_DRY_RELAX_BASE_K,
													 log_prior,
													 HK_BLIND_P9016_DRY_RELAX_TEMPERATURE);
	failed |= check_posterior_health(set, &mean_entropy, &mean_pmax);

	failed |= check_i32("dry-relax sep before",
						hk_blind_homolog_sep_compute_stats(bmap->n_beads, diploid,
														   HK_BLIND_P9016_DRY_RELAX_MIN_SEP_UNIT * HK_BLIND_P9016_DRY_RELAX_UNIT,
														   &sep_before), 0);
	failed |= check_sep_stats("dry-relax sep before", bmap, &sep_before);

	failed |= check_i32("dry-relax wedge build",
						hk_blind_wedge_list_build_from_bpair_set(&wedges, set,
																 HK_BLIND_P9016_DRY_RELAX_BASE_K,
																 HK_BLIND_P9016_DRY_RELAX_D_SCALE,
																 HK_BLIND_P9016_DRY_RELAX_RHO_TRAIN), 0);
	failed |= check_i32("dry-relax wedge aggregate", hk_blind_wedge_list_aggregate_exact(&wedges), 0);
	failed |= check_wedge_health(&wedges, n_diploid, set->n_bpairs);
	wedge_sum_k = sum_wedge_k(&wedges);
	wedge_sum_k_before_relax = wedge_sum_k;
	n_edges_before_relax = wedges.n_edges;

	copy_coords(diploid_before, diploid, n_diploid);
	failed |= check_i32("dry-relax relax ret",
						hk_blind_relax_cpu(&conf, &wedges, bmap, bmap->n_beads, diploid,
										   HK_BLIND_P9016_DRY_RELAX_UNIT,
										   HK_BLIND_P9016_DRY_RELAX_STEP,
										   HK_BLIND_P9016_DRY_RELAX_N_STEPS,
										   HK_BLIND_P9016_DRY_RELAX_MIN_SEP_UNIT,
										   HK_BLIND_P9016_DRY_RELAX_LAMBDA_SEP,
										   0, &relax_diag), 0);
	failed |= check_relax_diag(&relax_diag);
	failed |= check_i32("dry-relax edge count stable", wedges.n_edges, n_edges_before_relax);
	failed |= check_true("dry-relax edge sum stable",
						 fabs(sum_wedge_k(&wedges) - wedge_sum_k_before_relax) <=
						 1e-6 * (fabs(wedge_sum_k_before_relax) > 1.0? fabs(wedge_sum_k_before_relax) : 1.0));
	failed |= check_coords_finite(diploid, n_diploid);
	moved = coords_any_changed(diploid_before, diploid, n_diploid, 1e-8f);
	if (relax_diag.max_force_l1 > 1e-12f)
		failed |= check_true("dry-relax coordinate movement", moved);
	else
		fprintf(stderr, "NOTE: P9016 dry-relax force_l1 is zero; no movement expected\n");

	failed |= check_i32("dry-relax sep after",
						hk_blind_homolog_sep_compute_stats(bmap->n_beads, diploid,
														   HK_BLIND_P9016_DRY_RELAX_MIN_SEP_UNIT * HK_BLIND_P9016_DRY_RELAX_UNIT,
														   &sep_after), 0);
	failed |= check_sep_stats("dry-relax sep after", bmap, &sep_after);

	if (!failed) {
		fprintf(stderr,
				"P9016 dry-relax: n_raw=%d n_bpair=%d n_edges=%d skipped_self=%lld "
				"sum_k=%.8g mean_entropy=%.8g mean_pmax=%.8g "
				"force_l1=%.8g backbone_energy=%.8g moved=%d sep_mean_before=%.8g sep_mean_after=%.8g\n",
				set->n_raw, set->n_bpairs, wedges.n_edges, (long long)wedges.n_skipped_self_edges,
				wedge_sum_k, mean_entropy, mean_pmax, relax_diag.max_force_l1, relax_diag.initial_backbone_energy, moved,
				sep_before.mean_sep, sep_after.mean_sep);
	}

cleanup:
	free(haploid);
	free(diploid);
	free(diploid_before);
	hk_blind_wedge_list_destroy(&wedges);
	hk_blind_bpair_set_destroy(set);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
