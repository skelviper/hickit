#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "hickit.h"

#define DEFAULT_PAIRS "test_res/20260522_p9016_positive_control_mechanism/native_hickit/P9016.impute.pairs.gz"
#define BIN_SIZE_BP 1000000
#define NEIGHBOR_RADIUS 10000000
#define HARD_PC_SEED 17ULL

static const char *env_or_default(const char *name, const char *fallback)
{
	const char *s = getenv(name);
	return s && s[0]? s : fallback;
}

static float env_float_or_default(const char *name, float fallback)
{
	const char *s = getenv(name);
	char *end = 0;
	double v;
	if (s == 0 || s[0] == 0) return fallback;
	v = strtod(s, &end);
	if (end == s || *end != 0 || !isfinite(v)) return fallback;
	return (float)v;
}

static int mkdir_if_missing(const char *path)
{
	if (mkdir(path, 0775) == 0 || errno == EEXIST) return 0;
	fprintf(stderr, "failed to create %s: %s\n", path, strerror(errno));
	return -1;
}

static int mkdir_p(const char *path)
{
	char buf[4096];
	size_t i, n;
	if (path == 0 || path[0] == 0) return -1;
	n = strlen(path);
	if (n >= sizeof(buf)) return -1;
	memcpy(buf, path, n + 1);
	for (i = 1; i < n; ++i) {
		if (buf[i] != '/') continue;
		buf[i] = 0;
		if (buf[0] != 0 && mkdir_if_missing(buf) != 0) return -1;
		buf[i] = '/';
	}
	return mkdir_if_missing(buf);
}

static void path_join(char *out, size_t out_size, const char *a, const char *b)
{
	int n = snprintf(out, out_size, "%s/%s", a, b);
	assert(n >= 0 && (size_t)n < out_size);
}

static void rstrip(char *s)
{
	size_t n = strlen(s);
	while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r' ||
					 isspace((unsigned char)s[n - 1])))
		s[--n] = 0;
}

static int split_tsv(char *line, char **fields, int max_fields)
{
	int n = 0;
	char *p = line;
	while (n < max_fields) {
		fields[n++] = p;
		p = strchr(p, '\t');
		if (p == 0) break;
		*p++ = 0;
	}
	return n;
}

static int find_bpair_by_bids(const struct hk_blind_bpair_set *set, int32_t bid0, int32_t bid1)
{
	int32_t lo = 0, hi = set->n_bpairs - 1;
	if (bid1 < bid0) {
		int32_t t = bid0; bid0 = bid1; bid1 = t;
	}
	while (lo <= hi) {
		int32_t mid = lo + (hi - lo) / 2;
		const struct hk_blind_bpair *bp = &set->bpairs[mid];
		if (bp->key.bid[0] == bid0 && bp->key.bid[1] == bid1) return mid;
		if (bp->key.bid[0] < bid0 ||
			(bp->key.bid[0] == bid0 && bp->key.bid[1] < bid1))
			lo = mid + 1;
		else hi = mid - 1;
	}
	return -1;
}

static void refresh_diag_from_p4(struct hk_blind_bpair *bp)
{
	const float log4 = logf((float)HK_BLIND_N_STATE);
	float top1 = -1.0f, top2 = -1.0f;
	float entropy = 0.0f;
	int s;
	for (s = 0; s < HK_BLIND_N_STATE; ++s) {
		float p = bp->p4[s];
		if (p > 0.0f) entropy -= p * logf(p);
		if (p > top1) {
			top2 = top1;
			top1 = p;
		} else if (p > top2) {
			top2 = p;
		}
	}
	if (top2 < 0.0f) top2 = 0.0f;
	bp->entropy = entropy;
	bp->pmax = top1;
	bp->margin = top1 - top2;
	bp->rho_output = 1.0f - entropy / log4;
	if (bp->rho_output < 0.0f) bp->rho_output = 0.0f;
	if (bp->rho_output > 1.0f) bp->rho_output = 1.0f;
	bp->pU = 1.0f - bp->rho_output;
}

static int load_posterior_into_set(const char *path, struct hk_blind_bpair_set *set)
{
	FILE *fp = fopen(path, "r");
	char line[8192];
	int64_t n_loaded = 0;
	if (fp == 0) {
		fprintf(stderr, "failed to open posterior %s: %s\n", path, strerror(errno));
		return -1;
	}
	if (fgets(line, sizeof(line), fp) == 0) {
		fclose(fp);
		return -1;
	}
	while (fgets(line, sizeof(line), fp)) {
		char *fields[32];
		int n, id, s;
		double sum = 0.0;
		int32_t bid0, bid1;
		rstrip(line);
		n = split_tsv(line, fields, 32);
		if (n < 16) {
			fclose(fp);
			return -1;
		}
		bid0 = (int32_t)strtol(fields[6], 0, 10);
		bid1 = (int32_t)strtol(fields[7], 0, 10);
		id = find_bpair_by_bids(set, bid0, bid1);
		if (id < 0) {
			fprintf(stderr, "posterior bpair not found: %d %d\n", bid0, bid1);
			fclose(fp);
			return -1;
		}
		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			float v = (float)strtod(fields[11 + s], 0);
			set->bpairs[id].p4[s] = v;
			sum += v;
		}
		if (sum <= 0.0 || !isfinite(sum)) {
			fclose(fp);
			return -1;
		}
		for (s = 0; s < HK_BLIND_N_STATE; ++s)
			set->bpairs[id].p4[s] = (float)((double)set->bpairs[id].p4[s] / sum);
		refresh_diag_from_p4(&set->bpairs[id]);
		++n_loaded;
	}
	fclose(fp);
	if (n_loaded != set->n_bpairs) {
		fprintf(stderr, "posterior rows loaded %lld but set has %d bpairs\n",
				(long long)n_loaded, set->n_bpairs);
		return -1;
	}
	return 0;
}

static int hard_pc_kind_from_string(const char *s)
{
	if (s == 0 || s[0] == 0 || strcmp(s, "state_freq") == 0 ||
		strcmp(s, "statefreq") == 0)
		return 1;
	if (strcmp(s, "imputed_p4_top") == 0 || strcmp(s, "impp4top") == 0)
		return 2;
	if (strcmp(s, "none") == 0 || strcmp(s, "no_pc") == 0)
		return 0;
	fprintf(stderr, "unknown hard pc kind: %s\n", s);
	return -1;
}

static int mstep_mode_from_string(const char *s)
{
	if (strcmp(s, "raw_expected_pcut") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT;
	if (strcmp(s, "raw_expected_pcut_oracle_cis") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_CIS;
	if (strcmp(s, "raw_expected_pcut_oracle_all") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_PCUT_ORACLE_ALL;
	if (strcmp(s, "raw_expected_oracle_all") == 0 ||
		strcmp(s, "raw_oracle_all") == 0 ||
		strcmp(s, "oracle_all") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_ORACLE_ALL;
	if (strcmp(s, "raw_expected_count") == 0 || strcmp(s, "raw_expected") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_COUNT;
	if (strcmp(s, "raw_expected_soft_all") == 0 ||
		strcmp(s, "soft_all") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_ALL;
	if (strcmp(s, "raw_expected_soft_outlier") == 0 ||
		strcmp(s, "soft_outlier") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER;
	if (strcmp(s, "raw_expected_locked_only") == 0)
		return HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_LOCKED_ONLY;
	fprintf(stderr, "unknown mstep mode: %s\n", s);
	return -1;
}

static int write_edges_tsv(const char *path, const struct hk_bmap *bmap,
						   const struct hk_blind_wedge_list *list)
{
	FILE *fp = fopen(path, "w");
	int32_t i;
	if (fp == 0) {
		fprintf(stderr, "failed to open %s: %s\n", path, strerror(errno));
		return -1;
	}
	fprintf(fp, "diploid_bid0\tdiploid_bid1\thaploid_bid0\tcopy0\thaploid_bid1\tcopy1\tclass\tk\td_scale\teffective_count\n");
	for (i = 0; i < list->n_edges; ++i) {
		const struct hk_blind_wedge *e = &list->edges[i];
		int32_t hb0 = hk_diploid_haploid_bid(e->bid[0]);
		int32_t hb1 = hk_diploid_haploid_bid(e->bid[1]);
		int cls = bmap->beads[hb0].chr == bmap->beads[hb1].chr?
			HK_BLIND_CONTACT_CIS : HK_BLIND_CONTACT_TRANS;
		double eff = e->d_scale > 0.0f? pow((double)e->d_scale, -3.0) : NAN;
		fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%s\t%.9g\t%.9g\t%.9g\n",
				e->bid[0], e->bid[1], hb0, hk_diploid_copy(e->bid[0]),
				hb1, hk_diploid_copy(e->bid[1]),
				hk_blind_contact_class_name(cls), e->k, e->d_scale, eff);
	}
	if (fclose(fp) != 0)
		return -1;
	return 0;
}

static int write_stats_tsv(const char *path, const struct hk_bmap *bmap,
						   const struct hk_blind_wedge_list *list)
{
	FILE *fp = fopen(path, "w");
	int32_t i;
	int64_t n_cis = 0, n_trans = 0;
	double sum_k = 0.0, sum_k_cis = 0.0, sum_k_trans = 0.0;
	double sum_eff = 0.0, sum_eff_cis = 0.0, sum_eff_trans = 0.0;
	if (fp == 0) {
		fprintf(stderr, "failed to open %s: %s\n", path, strerror(errno));
		return -1;
	}
	for (i = 0; i < list->n_edges; ++i) {
		const struct hk_blind_wedge *e = &list->edges[i];
		int32_t hb0 = hk_diploid_haploid_bid(e->bid[0]);
		int32_t hb1 = hk_diploid_haploid_bid(e->bid[1]);
		double eff = e->d_scale > 0.0f? pow((double)e->d_scale, -3.0) : 0.0;
		int cis = bmap->beads[hb0].chr == bmap->beads[hb1].chr;
		sum_k += e->k;
		sum_eff += eff;
		if (cis) {
			++n_cis;
			sum_k_cis += e->k;
			sum_eff_cis += eff;
		} else {
			++n_trans;
			sum_k_trans += e->k;
			sum_eff_trans += eff;
		}
	}
	fprintf(fp, "n_edges\tn_cis\tn_trans\tsum_k\tsum_k_cis\tsum_k_trans\tsum_k_cis_trans_ratio\tsum_effective_count\tsum_effective_count_cis\tsum_effective_count_trans\teffective_cis_trans_ratio\n");
	fprintf(fp, "%d\t%lld\t%lld\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
			list->n_edges, (long long)n_cis, (long long)n_trans,
			sum_k, sum_k_cis, sum_k_trans,
			sum_k_trans > 0.0? sum_k_cis / sum_k_trans : NAN,
			sum_eff, sum_eff_cis, sum_eff_trans,
			sum_eff_trans > 0.0? sum_eff_cis / sum_eff_trans : NAN);
	if (fclose(fp) != 0)
		return -1;
	return 0;
}

int main(int argc, char **argv)
{
	const char *pairs_path = argc > 1? argv[1] : env_or_default("HK_P9016_IMPUTED_PAIRS", DEFAULT_PAIRS);
	const char *posterior_path = argc > 2? argv[2] : env_or_default("HK_P9016_BPAIR_POSTERIOR", "");
	const char *out_dir = argc > 3? argv[3] : env_or_default("HK_P9016_GRAPH_OUT", "test_res/tmp_mstep_graph_dump");
	const char *mstep_name = argc > 4? argv[4] : env_or_default("HK_P9016_MSTEP_MODE", "raw_expected_pcut");
	const char *hard_kind_name = argc > 5? argv[5] : env_or_default("HK_P9016_HARD_PC_KIND", "state_freq");
	float state_p_min = env_float_or_default("HK_P9016_RAW_SPLIT_STATE_P_MIN", 0.90f);
	float contact_k_multiplier_cis = env_float_or_default("HK_P9016_CONTACT_K_MULTIPLIER_CIS", 1.0f);
	float contact_k_multiplier_trans = env_float_or_default("HK_P9016_CONTACT_K_MULTIPLIER_TRANS", 1.0f);
	float trans_d_scale_multiplier = env_float_or_default("HK_P9016_TRANS_D_SCALE_MULTIPLIER", 1.0f);
	char edge_path[4096], stats_path[4096];
	struct hk_map *map = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	struct hk_blind_bpair_set *set = 0;
	struct hk_blind_wedge_list wedges;
	struct hk_blind_phase_lock_diag lock_diag;
	int mstep_mode = mstep_mode_from_string(mstep_name);
	int hard_kind = hard_pc_kind_from_string(hard_kind_name);
	int32_t i;
	int ret = 1;

	hk_blind_wedge_list_init(&wedges);
	hk_blind_phase_lock_diag_init(&lock_diag);
	if (mstep_mode < 0 || hard_kind < 0)
		return 1;
	if (trans_d_scale_multiplier <= 0.0f)
		trans_d_scale_multiplier = 1.0f;
	if (posterior_path == 0 || posterior_path[0] == 0) {
		fprintf(stderr, "posterior path is required\n");
		return 1;
	}
	if (mkdir_p(out_dir) != 0)
		return 1;
	path_join(edge_path, sizeof(edge_path), out_dir, "mstep_edges.tsv");
	path_join(stats_path, sizeof(stats_path), out_dir, "mstep_graph_stats.tsv");

	map = hk_map_read(pairs_path);
	if (map == 0 || map->pairs == 0 || map->n_pairs <= 0)
		goto cleanup;
	hk_pair_count_nei(map->n_pairs, map->pairs, NEIGHBOR_RADIUS, NEIGHBOR_RADIUS);
	bmap = hk_bmap_gen(map->d, map->n_pairs, map->pairs, BIN_SIZE_BP, 1);
	if (bmap == 0)
		goto cleanup;
	raw = (struct hk_blind_pair*)calloc((size_t)map->n_pairs, sizeof(*raw));
	if (raw == 0)
		goto cleanup;
	for (i = 0; i < map->n_pairs; ++i)
		hk_blind_pair_from_pair(&raw[i], &map->pairs[i]);
	set = hk_blind_bpair_set_build(bmap, map->n_pairs, raw);
	if (set == 0)
		goto cleanup;
	if (hk_blind_bpair_set_apply_base_k_mode(bmap, set, HK_BLIND_BASE_K_NEIGHBOR_MEDIAN) != 0)
		goto cleanup;
	hk_blind_bpair_set_init_uniform_prior(set);
	if (hard_kind == 1) {
		if (hk_blind_bpair_set_apply_raw_phase_state_freq_locks(set, map->pairs, map->n_pairs,
																100, HARD_PC_SEED, &lock_diag) != 0)
			goto cleanup;
	} else if (hard_kind == 2) {
		if (hk_blind_bpair_set_apply_imputed_p4_top_locks(set, map->pairs, map->n_pairs,
														  100, HARD_PC_SEED, 0.75f, &lock_diag) != 0)
			goto cleanup;
	}
	if (load_posterior_into_set(posterior_path, set) != 0)
		goto cleanup;
	if (mstep_mode == HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_OUTLIER) {
		fprintf(stderr, "raw_expected_soft_outlier dump requires geometry log-normalizers from a live E-step; falling back to soft-all semantics for loaded posterior TSV\n");
		mstep_mode = HK_BLIND_MSTEP_GRAPH_RAW_EXPECTED_SOFT_ALL;
	}
	if (hk_blind_wedge_list_build_mstep_graph_ex(&wedges, bmap, set, 1.0f,
												 HK_BLIND_RHO_TRAIN_CONSTANT, 0.0f,
												 HK_BLIND_D_SCALE_RAW_COUNT, 1e-6f,
												 contact_k_multiplier_cis,
												 contact_k_multiplier_trans,
												 HK_BLIND_STATE_WEIGHT_POSTERIOR,
												 mstep_mode,
												 0.0f, 0.0f, 1.0f,
												 HK_BLIND_RAW_SPLIT_CONF_ENTROPY_LINEAR,
												 0.0f, 1.0f, 1.0f,
												 state_p_min, 1.0f, 0,
												 HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_CIS,
												 HK_BLIND_RAW_OUTLIER_DEFAULT_BETA_TRANS,
												 HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_CIS,
												 HK_BLIND_RAW_OUTLIER_DEFAULT_PRIOR_TRANS,
												 0.0f, 1.0f) != 0)
		goto cleanup;
	hk_blind_wedge_list_apply_trans_d_scale_multiplier(bmap, &wedges,
													   trans_d_scale_multiplier);
	if (hk_blind_wedge_list_aggregate_exact(&wedges) != 0)
		goto cleanup;
	if (write_edges_tsv(edge_path, bmap, &wedges) != 0)
		goto cleanup;
	if (write_stats_tsv(stats_path, bmap, &wedges) != 0)
		goto cleanup;
	fprintf(stderr, "wrote %s and %s\n", edge_path, stats_path);
	ret = 0;

cleanup:
	hk_blind_wedge_list_destroy(&wedges);
	if (set) hk_blind_bpair_set_destroy(set);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (map) hk_map_destroy(map);
	return ret;
}
