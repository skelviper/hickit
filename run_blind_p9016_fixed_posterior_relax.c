#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "hickit.h"
#include "hkpriv.h"

#define HK_FIXED_UNIT 1.0f
#define HK_FIXED_MAX_LINE 8192
#define HK_FIXED_REPULSION_MULTIPLIER 1.0f

enum fixed_graph_mode {
	FIXED_GRAPH_SOFTALL = 0,
	FIXED_GRAPH_TOP1_ONLY = 1,
	FIXED_GRAPH_TOP2_ONLY = 2,
	FIXED_GRAPH_SAME_CROSS_GATED = 3
};

struct fixed_params {
	const char *source_dir;
	const char *output_dir;
	const char *config_name;
	const char *sample;
	const char *input_contact_source;
	int graph_mode;
	int d_scale_mode;
	float d_scale_eps_count;
	float d_scale_posterior_gamma;
	float trans_dscale_multiplier;
	float min_sep_unit;
	float lambda_sep;
	float relax_step;
	float repulsion_multiplier;
	int32_t relax_steps;
};

struct fixed_source {
	struct hk_bmap bmap;
	fvec3_t *coords;
	struct hk_blind_bpair_set set;
};

struct sep_dist_stats {
	int32_t n_beads;
	float min_sep_unit;
	float lambda_sep;
	float sep_min;
	float sep_p01;
	float sep_p05;
	float sep_p10;
	float sep_p25;
	float sep_median;
	float sep_mean;
	float sep_p75;
	float sep_p90;
	float sep_p95;
	float sep_p99;
	float sep_max;
	int32_t n_sep_below_min;
	float frac_sep_below_min;
};

static const char *env_or_default(const char *key, const char *default_value)
{
	const char *v = getenv(key);
	return (v && v[0])? v : default_value;
}

static int env_int_or_default(const char *key, int default_value)
{
	const char *v = getenv(key);
	char *end = 0;
	long x;
	if (v == 0 || v[0] == 0) return default_value;
	errno = 0;
	x = strtol(v, &end, 10);
	if (errno != 0 || end == v || *end != 0 || x < INT32_MIN || x > INT32_MAX) {
		fprintf(stderr, "invalid integer env %s=%s\n", key, v);
		exit(2);
	}
	return (int)x;
}

static float env_float_or_default(const char *key, float default_value, float min_value)
{
	const char *v = getenv(key);
	char *end = 0;
	float x;
	if (v == 0 || v[0] == 0) return default_value;
	errno = 0;
	x = strtof(v, &end);
	if (errno != 0 || end == v || *end != 0 || !isfinite(x) || x < min_value) {
		fprintf(stderr, "invalid float env %s=%s\n", key, v);
		exit(2);
	}
	return x;
}

static void make_parent_dirs(const char *path)
{
	char tmp[4096];
	size_t n, i;
	if (path == 0 || path[0] == 0) return;
	n = strlen(path);
	if (n >= sizeof(tmp)) {
		fprintf(stderr, "path too long: %s\n", path);
		exit(2);
	}
	memcpy(tmp, path, n + 1);
	for (i = 1; i < n; ++i) {
		if (tmp[i] == '/') {
			tmp[i] = 0;
			if (mkdir(tmp, 0775) != 0 && errno != EEXIST) {
				perror(tmp);
				exit(2);
			}
			tmp[i] = '/';
		}
	}
	if (mkdir(tmp, 0775) != 0 && errno != EEXIST) {
		perror(tmp);
		exit(2);
	}
}

static void path_join(char *dst, size_t dst_size, const char *dir, const char *name)
{
	int n = snprintf(dst, dst_size, "%s/%s", dir, name);
	if (n < 0 || (size_t)n >= dst_size) {
		fprintf(stderr, "path too long: %s/%s\n", dir, name);
		exit(2);
	}
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
	if (n > 0) {
		char *q = fields[n - 1] + strlen(fields[n - 1]);
		while (q > fields[n - 1] && (q[-1] == '\n' || q[-1] == '\r'))
			*--q = 0;
	}
	return n;
}

static int parse_contact_class(const char *s)
{
	if (strcmp(s, "cis") == 0) return HK_BLIND_CONTACT_CIS;
	if (strcmp(s, "trans") == 0) return HK_BLIND_CONTACT_TRANS;
	return HK_BLIND_CONTACT_CIS;
}

static int parse_graph_mode(const char *s)
{
	if (strcmp(s, "softall") == 0) return FIXED_GRAPH_SOFTALL;
	if (strcmp(s, "top1_only") == 0) return FIXED_GRAPH_TOP1_ONLY;
	if (strcmp(s, "top2_only") == 0) return FIXED_GRAPH_TOP2_ONLY;
	if (strcmp(s, "same_cross_gated") == 0) return FIXED_GRAPH_SAME_CROSS_GATED;
	fprintf(stderr, "invalid HK_FIXED_POSTERIOR_GRAPH_MODE=%s\n", s);
	exit(2);
}

static const char *graph_mode_name(int mode)
{
	switch (mode) {
	case FIXED_GRAPH_SOFTALL: return "softall";
	case FIXED_GRAPH_TOP1_ONLY: return "top1_only";
	case FIXED_GRAPH_TOP2_ONLY: return "top2_only";
	case FIXED_GRAPH_SAME_CROSS_GATED: return "same_cross_gated";
	default: return "unknown";
	}
}

static int parse_d_scale_mode(const char *s)
{
	if (strcmp(s, "raw_count") == 0) return HK_BLIND_D_SCALE_RAW_COUNT;
	if (strcmp(s, "posterior_count") == 0) return HK_BLIND_D_SCALE_POSTERIOR_COUNT;
	if (strcmp(s, "expected_count") == 0) return HK_BLIND_D_SCALE_POSTERIOR_COUNT;
	if (strcmp(s, "tempered_posterior_count") == 0) return HK_BLIND_D_SCALE_TEMPERED_POSTERIOR_COUNT;
	fprintf(stderr, "invalid HK_FIXED_POSTERIOR_D_SCALE_MODE=%s\n", s);
	exit(2);
}

static int p4_top_state(const float p4[HK_BLIND_N_STATE])
{
	int best = 0, s;
	for (s = 1; s < HK_BLIND_N_STATE; ++s)
		if (p4[s] > p4[best])
			best = s;
	return best;
}

static void p4_top2_states(const float p4[HK_BLIND_N_STATE], int *top0, int *top1)
{
	int a = 0, b = 1, s;
	if (p4[b] > p4[a]) {
		int t = a;
		a = b;
		b = t;
	}
	for (s = 2; s < HK_BLIND_N_STATE; ++s) {
		if (p4[s] > p4[a]) {
			b = a;
			a = s;
		} else if (p4[s] > p4[b]) {
			b = s;
		}
	}
	*top0 = a;
	*top1 = b;
}

static int state_is_same(int state)
{
	return state == HK_BLIND_STATE_00 || state == HK_BLIND_STATE_11;
}

static int keep_state_for_mode(const struct hk_blind_bpair *bp, int mode, int state)
{
	if (mode == FIXED_GRAPH_SOFTALL)
		return 1;
	if (mode == FIXED_GRAPH_TOP1_ONLY)
		return state == p4_top_state(bp->p4);
	if (mode == FIXED_GRAPH_TOP2_ONLY) {
		int a, b;
		p4_top2_states(bp->p4, &a, &b);
		return state == a || state == b;
	}
	if (mode == FIXED_GRAPH_SAME_CROSS_GATED) {
		float psame = bp->p4[HK_BLIND_STATE_00] + bp->p4[HK_BLIND_STATE_11];
		float pcross = bp->p4[HK_BLIND_STATE_01] + bp->p4[HK_BLIND_STATE_10];
		int keep_same = psame >= pcross;
		return state_is_same(state) == keep_same;
	}
	return 0;
}

static void fixed_wedge_set(struct hk_blind_wedge *edge, int32_t bid0, int32_t bid1,
							float k, float d_scale, int state)
{
	if (bid1 < bid0) {
		int32_t t = bid0;
		bid0 = bid1;
		bid1 = t;
	}
	edge->bid[0] = bid0;
	edge->bid[1] = bid1;
	edge->k = k;
	edge->d_scale = d_scale;
	edge->state = (int8_t)state;
	edge->state_mask = HK_BLIND_STATE_MASK(state);
}

static void fixed_wedge_push(struct hk_blind_wedge_list *list, const struct hk_blind_wedge *edge)
{
	if (list->n_edges == list->m_edges)
		EXPAND(list->edges, list->m_edges);
	if (list->edges == 0) {
		fprintf(stderr, "out of memory while appending wedges\n");
		exit(2);
	}
	list->edges[list->n_edges++] = *edge;
}

static float state_d_scale(const struct hk_blind_bpair *bp, int state,
						   int d_scale_mode, float eps_count, float gamma,
						   float trans_multiplier)
{
	float d_scale = bp->base_d_scale;
	if (d_scale_mode == HK_BLIND_D_SCALE_POSTERIOR_COUNT ||
		d_scale_mode == HK_BLIND_D_SCALE_TEMPERED_POSTERIOR_COUNT) {
		float p = bp->p4[state];
		float n_eff;
		if (p < 0.0f) p = 0.0f;
		if (p > 1.0f) p = 1.0f;
		n_eff = (float)bp->n_raw * powf(p, gamma);
		if (!isfinite(n_eff) || n_eff < eps_count)
			n_eff = eps_count;
		d_scale = powf(n_eff, -1.0f / 3.0f);
	}
	if (bp->contact_class == HK_BLIND_CONTACT_TRANS)
		d_scale *= trans_multiplier;
	return d_scale;
}

static int build_fixed_wedges(struct hk_blind_wedge_list *out,
							  const struct hk_blind_bpair_set *set,
							  const struct fixed_params *p)
{
	int32_t i;
	hk_blind_wedge_list_init(out);
	out->n_input_bpair = set->n_bpairs;
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		int s;
		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			struct hk_blind_wedge edge;
			int32_t b0, b1;
			float k, d_scale;
			if (!keep_state_for_mode(bp, p->graph_mode, s))
				continue;
			k = bp->base_k * bp->p4[s];
			if (!isfinite(k) || k <= 0.0f)
				continue;
			switch (s) {
			case HK_BLIND_STATE_00:
				b0 = hk_diploid_bid(bp->key.bid[0], HK_DIPLOID_COPY0);
				b1 = hk_diploid_bid(bp->key.bid[1], HK_DIPLOID_COPY0);
				break;
			case HK_BLIND_STATE_01:
				b0 = hk_diploid_bid(bp->key.bid[0], HK_DIPLOID_COPY0);
				b1 = hk_diploid_bid(bp->key.bid[1], HK_DIPLOID_COPY1);
				break;
			case HK_BLIND_STATE_10:
				b0 = hk_diploid_bid(bp->key.bid[0], HK_DIPLOID_COPY1);
				b1 = hk_diploid_bid(bp->key.bid[1], HK_DIPLOID_COPY0);
				break;
			default:
				b0 = hk_diploid_bid(bp->key.bid[0], HK_DIPLOID_COPY1);
				b1 = hk_diploid_bid(bp->key.bid[1], HK_DIPLOID_COPY1);
				break;
			}
			if (b0 == b1) {
				++out->n_skipped_self_edges;
				continue;
			}
			d_scale = state_d_scale(bp, s, p->d_scale_mode, p->d_scale_eps_count,
									p->d_scale_posterior_gamma, p->trans_dscale_multiplier);
			if (!isfinite(d_scale) || d_scale <= 0.0f)
				continue;
			fixed_wedge_set(&edge, b0, b1, k, d_scale, s);
			fixed_wedge_push(out, &edge);
		}
	}
	out->n_expanded_edges = out->n_edges;
	out->n_edges_before_aggregation = out->n_edges;
	if (hk_blind_wedge_list_aggregate_exact(out) != 0)
		return -1;
	out->mean_rho_train_bpair = 1.0f;
	out->min_rho_train_bpair = 1.0f;
	out->max_rho_train_bpair = 1.0f;
	out->sum_rho_train_bpair = (double)out->n_edges;
	return 0;
}

static int chr_id_for_name(struct hk_sdict *d, const char *name)
{
	int32_t id = hk_sd_get(d, name);
	if (id >= 0) return id;
	id = hk_sd_put(d, name, -1);
	if (id < 0) {
		fprintf(stderr, "failed to add chromosome %s\n", name);
		exit(2);
	}
	return id;
}

static void ensure_bead_capacity(struct hk_bmap *bmap, int32_t n)
{
	int32_t old = bmap->n_beads;
	if (n <= bmap->n_beads) return;
	REALLOC(bmap->beads, n);
	if (bmap->beads == 0) {
		fprintf(stderr, "out of memory for beads\n");
		exit(2);
	}
	memset(bmap->beads + old, 0, (size_t)(n - old) * sizeof(*bmap->beads));
	bmap->n_beads = n;
}

static void ensure_coord_capacity(fvec3_t **coords, int32_t *m_coords, int32_t n)
{
	int32_t old = *m_coords;
	if (n <= *m_coords) return;
	*m_coords = n + n / 2 + 16;
	REALLOC(*coords, *m_coords);
	if (*coords == 0) {
		fprintf(stderr, "out of memory for coords\n");
		exit(2);
	}
	memset((*coords) + old, 0, (size_t)(*m_coords - old) * sizeof(**coords));
}

static int load_coords(const char *path, struct hk_bmap *bmap, fvec3_t **coords_out)
{
	FILE *fp = fopen(path, "r");
	char line[HK_FIXED_MAX_LINE];
	char *f[16];
	int m_coords = 0;
	fvec3_t *coords = 0;
	int32_t n_rows = 0;
	if (fp == 0) {
		perror(path);
		return -1;
	}
	memset(bmap, 0, sizeof(*bmap));
	bmap->d = hk_sd_init();
	bmap->unit = HK_FIXED_UNIT;
	if (bmap->d == 0)
		return -1;
	if (fgets(line, sizeof(line), fp) == 0) {
		fclose(fp);
		return -1;
	}
	while (fgets(line, sizeof(line), fp)) {
		int nf = split_tsv(line, f, 16);
		int chr, bid, diploid_bid;
		if (nf < 9) {
			fprintf(stderr, "bad coords row in %s\n", path);
			fclose(fp);
			return -1;
		}
		chr = chr_id_for_name(bmap->d, f[0]);
		bid = atoi(f[3]);
		diploid_bid = atoi(f[5]);
		if (bid < 0 || diploid_bid < 0) {
			fprintf(stderr, "negative bid in %s\n", path);
			fclose(fp);
			return -1;
		}
		ensure_bead_capacity(bmap, bid + 1);
		ensure_coord_capacity(&coords, &m_coords, diploid_bid + 1);
		bmap->beads[bid].chr = chr;
		bmap->beads[bid].st = atoi(f[1]);
		bmap->beads[bid].en = atoi(f[2]);
		coords[diploid_bid][0] = strtof(f[6], 0);
		coords[diploid_bid][1] = strtof(f[7], 0);
		coords[diploid_bid][2] = strtof(f[8], 0);
		++n_rows;
	}
	if (fclose(fp) != 0)
		return -1;
	if (n_rows != 2 * bmap->n_beads) {
		fprintf(stderr, "coords row count mismatch: rows=%d n_beads=%d\n",
				n_rows, bmap->n_beads);
		return -1;
	}
	hk_bmap_set_offcnt(bmap);
	*coords_out = coords;
	return 0;
}

static int load_posterior(const char *path, struct hk_blind_bpair_set *set)
{
	FILE *fp = fopen(path, "r");
	char line[HK_FIXED_MAX_LINE];
	char *f[32];
	int32_t m = 0;
	if (fp == 0) {
		perror(path);
		return -1;
	}
	memset(set, 0, sizeof(*set));
	if (fgets(line, sizeof(line), fp) == 0) {
		fclose(fp);
		return -1;
	}
	while (fgets(line, sizeof(line), fp)) {
		int nf = split_tsv(line, f, 32);
		struct hk_blind_bpair *bp;
		int s;
		if (nf < 23) {
			fprintf(stderr, "bad posterior row in %s\n", path);
			fclose(fp);
			return -1;
		}
		if (set->n_bpairs == m) {
			m = m? m + (m >> 1) : 1024;
			REALLOC(set->bpairs, m);
			if (set->bpairs == 0) {
				fprintf(stderr, "out of memory for bpairs\n");
				fclose(fp);
				return -1;
			}
		}
		bp = &set->bpairs[set->n_bpairs++];
		memset(bp, 0, sizeof(*bp));
		bp->key.bid[0] = atoi(f[6]);
		bp->key.bid[1] = atoi(f[7]);
		bp->n_raw = atoi(f[8]);
		bp->base_d_scale = strtof(f[9], 0);
		bp->base_k = strtof(f[10], 0);
		bp->p4[0] = strtof(f[11], 0);
		bp->p4[1] = strtof(f[12], 0);
		bp->p4[2] = strtof(f[13], 0);
		bp->p4[3] = strtof(f[14], 0);
		bp->pU = strtof(f[15], 0);
		bp->entropy = strtof(f[18], 0);
		bp->margin = strtof(f[19], 0);
		bp->pmax = strtof(f[20], 0);
		bp->rho_output = strtof(f[21], 0);
		bp->contact_class = (int8_t)parse_contact_class(f[22]);
		for (s = 0; s < HK_BLIND_N_STATE; ++s)
			bp->log_prior[s] = logf(0.25f);
		set->n_raw += bp->n_raw;
		if (bp->contact_class == HK_BLIND_CONTACT_CIS) {
			set->n_raw_cis += bp->n_raw;
			++set->n_bpair_cis;
		} else {
			set->n_raw_trans += bp->n_raw;
			++set->n_bpair_trans;
		}
	}
	if (fclose(fp) != 0)
		return -1;
	return 0;
}

static int fixed_source_load(const char *dir, struct fixed_source *src)
{
	char coords_path[4096], posterior_path[4096];
	memset(src, 0, sizeof(*src));
	path_join(coords_path, sizeof(coords_path), dir, "p9016_full.coords.tsv");
	path_join(posterior_path, sizeof(posterior_path), dir, "p9016_full.bpair_posterior.tsv");
	if (load_coords(coords_path, &src->bmap, &src->coords) != 0)
		return -1;
	if (load_posterior(posterior_path, &src->set) != 0)
		return -1;
	return 0;
}

static void fixed_source_destroy(struct fixed_source *src)
{
	if (src->bmap.d)
		hk_sd_destroy(src->bmap.d);
	free(src->bmap.beads);
	free(src->bmap.offcnt);
	free(src->bmap.pairs);
	free(src->bmap.x);
	free(src->bmap.feat);
	free(src->bmap.cpg);
	free(src->bmap.gc_bias);
	free(src->coords);
	free(src->set.bpairs);
	memset(src, 0, sizeof(*src));
}

static int float_cmp(const void *a_, const void *b_)
{
	const float a = *(const float*)a_;
	const float b = *(const float*)b_;
	return (a > b) - (a < b);
}

static float percentile_sorted(const float *x, int32_t n, float p)
{
	double idx, lo_f, frac;
	int32_t lo, hi;
	if (n <= 0) return NAN;
	idx = (double)(n - 1) * (double)p;
	lo_f = floor(idx);
	lo = (int32_t)lo_f;
	hi = lo + 1;
	if (hi >= n) return x[n - 1];
	frac = idx - lo_f;
	return (float)((double)x[lo] * (1.0 - frac) + (double)x[hi] * frac);
}

static void compute_sep_dist_stats(const struct hk_bmap *bmap, const fvec3_t *coords,
								   float min_sep_unit, float lambda_sep,
								   struct sep_dist_stats *out)
{
	int32_t i;
	float *sep = MALLOC(float, bmap->n_beads);
	double sum = 0.0;
	memset(out, 0, sizeof(*out));
	out->n_beads = bmap->n_beads;
	out->min_sep_unit = min_sep_unit;
	out->lambda_sep = lambda_sep;
	for (i = 0; i < bmap->n_beads; ++i) {
		int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		double dx = (double)coords[b0][0] - coords[b1][0];
		double dy = (double)coords[b0][1] - coords[b1][1];
		double dz = (double)coords[b0][2] - coords[b1][2];
		sep[i] = (float)sqrt(dx * dx + dy * dy + dz * dz);
		sum += sep[i];
		if (min_sep_unit > 0.0f && sep[i] < min_sep_unit)
			++out->n_sep_below_min;
	}
	qsort(sep, bmap->n_beads, sizeof(*sep), float_cmp);
	out->sep_min = percentile_sorted(sep, bmap->n_beads, 0.0f);
	out->sep_p01 = percentile_sorted(sep, bmap->n_beads, 0.01f);
	out->sep_p05 = percentile_sorted(sep, bmap->n_beads, 0.05f);
	out->sep_p10 = percentile_sorted(sep, bmap->n_beads, 0.10f);
	out->sep_p25 = percentile_sorted(sep, bmap->n_beads, 0.25f);
	out->sep_median = percentile_sorted(sep, bmap->n_beads, 0.50f);
	out->sep_mean = bmap->n_beads > 0? (float)(sum / bmap->n_beads) : NAN;
	out->sep_p75 = percentile_sorted(sep, bmap->n_beads, 0.75f);
	out->sep_p90 = percentile_sorted(sep, bmap->n_beads, 0.90f);
	out->sep_p95 = percentile_sorted(sep, bmap->n_beads, 0.95f);
	out->sep_p99 = percentile_sorted(sep, bmap->n_beads, 0.99f);
	out->sep_max = percentile_sorted(sep, bmap->n_beads, 1.0f);
	out->frac_sep_below_min = bmap->n_beads > 0?
		(float)out->n_sep_below_min / (float)bmap->n_beads : NAN;
	free(sep);
}

static int write_sep_diag(const char *path, const struct sep_dist_stats *s)
{
	FILE *fp = fopen(path, "w");
	if (fp == 0) return -1;
	if (fprintf(fp,
				"n_beads\tmin_sep_unit\tlambda_sep\tsep_min\tsep_p01\tsep_p05\tsep_p10\t"
				"sep_p25\tsep_median\tsep_mean\tsep_p75\tsep_p90\tsep_p95\tsep_p99\t"
				"sep_max\tn_sep_below_min\tfrac_sep_below_min\n"
				"%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
				"%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\n",
				s->n_beads, s->min_sep_unit, s->lambda_sep, s->sep_min, s->sep_p01,
				s->sep_p05, s->sep_p10, s->sep_p25, s->sep_median, s->sep_mean,
				s->sep_p75, s->sep_p90, s->sep_p95, s->sep_p99, s->sep_max,
				s->n_sep_below_min, s->frac_sep_below_min) < 0) {
		fclose(fp);
		return -1;
	}
	return fclose(fp) == 0? 0 : -1;
}

static int write_coords_outputs(const char *coords_path, const char *coords_gz_path,
								const struct hk_bmap *bmap, const fvec3_t *coords)
{
	FILE *fp = fopen(coords_path, "w");
	int ret;
	if (fp == 0) return -1;
	ret = hk_blind_write_diploid_coords_tsv(fp, bmap, coords);
	if (fclose(fp) != 0) ret = -1;
	if (ret != 0) return ret;
	return hk_blind_write_diploid_coords_tsv_gz(coords_gz_path, bmap, coords);
}

static int write_posterior_output(const char *path, const struct hk_bmap *bmap,
								  const struct hk_blind_bpair_set *set)
{
	FILE *fp = fopen(path, "w");
	int ret;
	if (fp == 0) return -1;
	ret = hk_blind_write_bpair_posterior_tsv(fp, bmap, set);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static double mean_entropy(const struct hk_blind_bpair_set *set)
{
	int32_t i;
	double sum = 0.0, wsum = 0.0;
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		sum += (double)bp->entropy * (double)bp->n_raw;
		wsum += bp->n_raw;
	}
	return wsum > 0.0? sum / wsum : NAN;
}

static double mean_pU(const struct hk_blind_bpair_set *set)
{
	int32_t i;
	double sum = 0.0, wsum = 0.0;
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		sum += (double)bp->pU * (double)bp->n_raw;
		wsum += bp->n_raw;
	}
	return wsum > 0.0? sum / wsum : NAN;
}

static int write_loop_diag(const char *path, const struct hk_blind_relax_diag *diag,
						   const struct sep_dist_stats *sep,
						   const struct hk_blind_wedge_list *wedges,
						   const struct hk_blind_bpair_set *set)
{
	FILE *fp = fopen(path, "w");
	double sum_wedge_k = 0.0;
	int32_t i;
	if (fp == 0) return -1;
	for (i = 0; i < wedges->n_edges; ++i)
		sum_wedge_k += wedges->edges[i].k;
	if (fprintf(fp,
				"n_iter\tn_completed\tinitial_mean_entropy\tfinal_mean_entropy\t"
				"initial_mean_pU\tfinal_mean_pU\tinitial_temperature\tfinal_temperature\t"
				"initial_rho_train\tfinal_rho_train\ttotal_chr_flipped\tn_bad_iter\t"
				"n_relax_nonfinite_iter\tn_coord_nonfinite\tfinal_mean_sep\tfinal_min_sep\t"
				"final_max_sep\tfinal_sum_wedge_k\tfinal_mean_rho_train_bpair\t"
				"final_min_rho_train_bpair\tfinal_max_rho_train_bpair\tfinal_n_expanded_edges\t"
				"final_n_softall_candidate_pairs\tfinal_n_softall_filter1_pairs\t"
				"final_n_softall_filter2_pairs\tfinal_n_softall_bmap_pairs\t"
				"final_n_softall_selected_raw\tfinal_n_softall_gate_skip_raw\t"
				"final_n_softall_same_bin_skip_raw\tfinal_n_skipped_same_bin_bpairs\t"
				"final_contact_energy\tfinal_repulsion_energy\tfinal_backbone_energy\t"
				"final_sep_energy\tfinal_sep_force_l1\tfinal_repulsion_force_l1\t"
				"final_backbone_force_l1\trepulsion_mode\tposterior_refreshed_after_final_relax\n") < 0) {
		fclose(fp);
		return -1;
	}
	if (fprintf(fp,
				"%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t1\t1\t1\t1\t0\t0\t%d\t%d\t"
				"%.9g\t%.9g\t%.9g\t%.17g\t1\t1\t1\t%lld\t0\t0\t0\t0\t0\t0\t0\t0\t"
				"%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t0\n",
				diag->n_steps, diag->n_completed,
				mean_entropy(set), mean_entropy(set),
				mean_pU(set), mean_pU(set),
				diag->n_nonfinite_step, diag->n_coord_nonfinite,
				sep->sep_mean, sep->sep_min, sep->sep_max,
				sum_wedge_k,
				(long long)wedges->n_expanded_edges,
				diag->final_contact_energy, diag->final_repulsion_energy,
				diag->final_backbone_energy, diag->final_sep_energy,
				diag->final_sep_force_l1, diag->final_repulsion_force_l1,
				diag->final_backbone_force_l1, diag->repulsion_mode) < 0) {
		fclose(fp);
		return -1;
	}
	return fclose(fp) == 0? 0 : -1;
}

static int write_force_class_diag(const char *path, const struct hk_fdg_conf *conf,
								  const struct hk_bmap *bmap,
								  const struct hk_blind_wedge_list *wedges,
								  const fvec3_t *coords,
								  const struct hk_blind_relax_diag *relax_diag,
								  struct hk_blind_contact_class_diag *diag_out)
{
	struct hk_blind_contact_class_diag diag;
	FILE *fp;
	if (hk_blind_contact_class_diag_accumulate(conf, bmap, wedges, coords,
											   bmap->n_beads, HK_FIXED_UNIT, &diag) != 0)
		return -1;
	if (diag_out)
		*diag_out = diag;
	fp = fopen(path, "w");
	if (fp == 0)
		return -1;
	if (fprintf(fp,
				"class\tn_wedges\tsum_wedge_k\tcontact_energy\tcontact_force_l1\t"
				"backbone_force_l1\trepulsion_force_l1\thomolog_sep_force_l1\tn_nonfinite\n"
				"cis\t%lld\t%.17g\t%.9g\t%.9g\t0\t0\t0\t%d\n"
				"trans\t%lld\t%.17g\t%.9g\t%.9g\t0\t0\t0\t%d\n"
				"total\t%lld\t%.17g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\n",
				(long long)diag.n_wedges_cis, diag.sum_wedge_k_cis,
				diag.contact_energy_cis, diag.contact_force_l1_cis, diag.n_nonfinite,
				(long long)diag.n_wedges_trans, diag.sum_wedge_k_trans,
				diag.contact_energy_trans, diag.contact_force_l1_trans, diag.n_nonfinite,
				(long long)(diag.n_wedges_cis + diag.n_wedges_trans),
				diag.sum_wedge_k_cis + diag.sum_wedge_k_trans,
				diag.contact_energy_cis + diag.contact_energy_trans,
				diag.contact_force_l1_cis + diag.contact_force_l1_trans,
				relax_diag? relax_diag->final_backbone_force_l1 : 0.0f,
				relax_diag? relax_diag->final_repulsion_force_l1 : 0.0f,
				relax_diag? relax_diag->final_sep_force_l1 : 0.0f,
				diag.n_nonfinite) < 0) {
		fclose(fp);
		return -1;
	}
	return fclose(fp) == 0? 0 : -1;
}

static int write_manifest(const char *path, const struct fixed_params *p,
						  const struct fixed_source *src,
						  const struct hk_blind_wedge_list *wedges,
						  const struct hk_blind_relax_diag *relax_diag,
						  const struct sep_dist_stats *sep,
						  const struct hk_blind_contact_class_diag *force_diag,
						  const char *posterior_path, const char *coords_path,
						  const char *coords_gz_path, const char *loop_diag_path,
						  const char *force_diag_path, const char *sep_diag_path)
{
	FILE *fp = fopen(path, "w");
	if (fp == 0) return -1;
	if (fprintf(fp,
				"key\tvalue\n"
				"sample\t%s\n"
				"runner_family\tp9016_fixed_posterior_diagnostic\n"
				"runner_version\tfixed_posterior_diagnostic_2026-06-18\n"
				"default_profile\tfixed_posterior_graph_semantics\n"
				"config_name\t%s\n"
				"git_commit\t%s\n"
				"git_dirty_count\t%s\n"
				"binary_hash\t%s\n"
				"input_path\t%s\n"
				"input_contact_source\t%s\n"
				"output_dir\t%s\n"
				"n_raw\t%d\n"
				"n_bpair\t%d\n"
				"n_beads\t%d\n"
				"resolution\t1000000\n"
				"bin_size_bp\t1000000\n"
				"resolution_label\t1Mb\n"
				"relax_backend\tcpu\n"
				"unit\t1\n"
				"d_scale\t1\n"
				"base_k_mode\tfixed_source_neighbor_median\n"
				"init_mode\tfixed_source_coords\n"
				"init_seed\t0\n"
				"init_eps\t0\n"
				"init_noise_scale\t0\n"
				"init_scale\t0\n"
				"init_eps_effective\t0\n"
				"init_noise_scale_effective\t0\n"
				"init_scale_effective\t0\n"
				"prior_mode\tuniform\n"
				"rho_train_mode\tconstant\n"
				"rho_train\t1\n"
				"rho_train_floor\t0\n"
				"baseline\tfixed_posterior_diagnostic\n"
				"mstep_graph_mode\tfixed_bpair_state_edges\n"
				"training_graph_mode\tfixed_bpair_state_edges_%s\n"
				"native_softall_replay\t0\n"
				"fixed_bpair_state_edges\t1\n"
				"fixed_posterior\t1\n"
				"fixed_posterior_source_dir\t%s\n"
				"fixed_posterior_graph_mode\t%s\n"
				"fixed_posterior_topology_note\tposterior_not_refreshed_after_relax\n"
				"edge_k_probability_weighted\t1\n"
				"dscale_mode\t%s\n"
				"d_scale_mode\t%s\n"
				"d_scale_mode_input_string\t%s\n"
				"d_scale_posterior_gamma\t%.9g\n"
				"dscale_posterior_gamma\t%.9g\n"
				"dscale_effective_count_formula\t%s\n"
				"dscale_probability_weighted\t%d\n"
				"legacy_expected_count_alias_used\t0\n"
				"d_scale_eps_count\t%.9g\n"
				"trans_dscale_multiplier\t%.9g\n"
				"estep_score_mode\tfdg_flat\n"
				"copy_labels_are_gauge_only\t1\n"
				"uses_phase_labels\t0\n"
				"uses_charm_or_reference\t1\n"
				"uses_charm_for_training\t1\n"
				"reference_derived_positive_control\t1\n"
				"reference_training_source_3dg\t%s\n"
				"min_sep_unit\t%.9g\n"
				"lambda_sep\t%.9g\n"
				"chr_sep_unit\t0\n"
				"lambda_chr_sep\t0\n"
				"relax_step\t%.9g\n"
				"relax_steps\t%d\n"
				"n_iter\t0\n"
				"temperature_start\t1\n"
				"temperature_end\t1\n"
				"rho_train_start\t1\n"
				"rho_train_end\t1\n"
				"enable_repulsion\t1\n"
				"repulsion_mode\t%d\n"
				"repulsion_multiplier\t%.9g\n"
				"k_rel_rep_effective\t%.9g\n"
				"n_fixed_edges\t%d\n"
				"n_fixed_edges_before_aggregation\t%lld\n"
				"n_fixed_edges_removed_by_aggregation\t%lld\n"
				"final_mean_entropy\t%.9g\n"
				"final_mean_pU\t%.9g\n"
				"final_mean_sep\t%.9g\n"
				"final_min_sep\t%.9g\n"
				"final_max_sep\t%.9g\n"
				"final_contact_energy\t%.9g\n"
				"final_repulsion_energy\t%.9g\n"
				"final_backbone_energy\t%.9g\n"
				"final_sep_energy\t%.9g\n"
				"final_sep_force_l1\t%.9g\n"
				"force_diag_total_contact_energy\t%.9g\n"
				"force_diag_total_contact_force_l1\t%.9g\n"
				"output_bpair_posterior\t%s\n"
				"output_coords\t%s\n"
				"output_coords_gz\t%s\n"
				"output_loop_diag\t%s\n"
				"output_force_class_diag\t%s\n"
				"output_sep_diag\t%s\n"
				"status\tOK\n",
				p->sample, p->config_name, env_or_default("HK_BLIND_GIT_COMMIT", "NA"),
				env_or_default("HK_BLIND_GIT_DIRTY_COUNT", "NA"),
				env_or_default("HK_BLIND_BINARY_HASH", "NA"),
				p->source_dir, p->input_contact_source, p->output_dir,
				src->set.n_raw, src->set.n_bpairs, src->bmap.n_beads,
				graph_mode_name(p->graph_mode), p->source_dir, graph_mode_name(p->graph_mode),
				hk_blind_d_scale_mode_name(p->d_scale_mode),
				hk_blind_d_scale_mode_name(p->d_scale_mode),
				hk_blind_d_scale_mode_name(p->d_scale_mode),
				p->d_scale_posterior_gamma, p->d_scale_posterior_gamma,
				hk_blind_d_scale_effective_count_formula(p->d_scale_mode, p->d_scale_posterior_gamma),
				((p->d_scale_mode == HK_BLIND_D_SCALE_POSTERIOR_COUNT ||
				  p->d_scale_mode == HK_BLIND_D_SCALE_TEMPERED_POSTERIOR_COUNT) &&
				 p->d_scale_posterior_gamma > 0.0f)? 1 : 0,
				p->d_scale_eps_count, p->trans_dscale_multiplier,
				env_or_default("HK_BLIND_P9016_REFERENCE_SOURCE_TDG", "CHARM_3DG_derived_positive_control"),
				p->min_sep_unit, p->lambda_sep, p->relax_step, p->relax_steps,
				HK_BLIND_REPULSION_CELL, p->repulsion_multiplier,
				0.05f * p->repulsion_multiplier,
				wedges->n_edges, (long long)wedges->n_edges_before_aggregation,
				(long long)wedges->n_aggregated_edges_removed,
				mean_entropy(&src->set), mean_pU(&src->set),
				sep->sep_mean, sep->sep_min, sep->sep_max,
				relax_diag->final_contact_energy, relax_diag->final_repulsion_energy,
				relax_diag->final_backbone_energy, relax_diag->final_sep_energy,
				relax_diag->final_sep_force_l1,
				force_diag->contact_energy_cis + force_diag->contact_energy_trans,
				force_diag->contact_force_l1_cis + force_diag->contact_force_l1_trans,
				posterior_path, coords_path, coords_gz_path, loop_diag_path,
				force_diag_path, sep_diag_path) < 0) {
		fclose(fp);
		return -1;
	}
	return fclose(fp) == 0? 0 : -1;
}

static void load_params(struct fixed_params *p)
{
	const char *dmode = env_or_default("HK_FIXED_POSTERIOR_D_SCALE_MODE", "posterior_count");
	memset(p, 0, sizeof(*p));
	p->source_dir = env_or_default("HK_FIXED_POSTERIOR_SOURCE_DIR", "");
	p->output_dir = env_or_default("HK_FIXED_POSTERIOR_OUTPUT_DIR", "");
	p->config_name = env_or_default("HK_FIXED_POSTERIOR_CONFIG_NAME", "fixed_posterior_softall");
	p->sample = env_or_default("HK_BLIND_SAMPLE", env_or_default("HK_BLIND_P9016_SAMPLE", "P9016"));
	p->input_contact_source = env_or_default("HK_FIXED_POSTERIOR_INPUT_CONTACT_SOURCE", "charm3dg_derived_pairs");
	p->graph_mode = parse_graph_mode(env_or_default("HK_FIXED_POSTERIOR_GRAPH_MODE", "softall"));
	p->d_scale_mode = parse_d_scale_mode(dmode);
	p->d_scale_eps_count = env_float_or_default("HK_FIXED_POSTERIOR_D_SCALE_EPS_COUNT", 1e-6f, 0.0f);
	p->d_scale_posterior_gamma = env_float_or_default("HK_FIXED_POSTERIOR_D_SCALE_POSTERIOR_GAMMA", 1.0f, 0.0f);
	p->trans_dscale_multiplier = env_float_or_default("HK_FIXED_POSTERIOR_TRANS_DSCALE_MULTIPLIER", 1.0f, 0.0f);
	p->min_sep_unit = env_float_or_default("HK_FIXED_POSTERIOR_MIN_SEP_UNIT", 1.5f, 0.0f);
	p->lambda_sep = env_float_or_default("HK_FIXED_POSTERIOR_LAMBDA_SEP", 1.0f, 0.0f);
	p->relax_step = env_float_or_default("HK_FIXED_POSTERIOR_RELAX_STEP", 0.012f, 0.0f);
	p->repulsion_multiplier = env_float_or_default("HK_FIXED_POSTERIOR_REPULSION_MULTIPLIER",
												   HK_FIXED_REPULSION_MULTIPLIER, 0.0f);
	p->relax_steps = env_int_or_default("HK_FIXED_POSTERIOR_RELAX_STEPS", 100);
	if (p->source_dir[0] == 0 || p->output_dir[0] == 0) {
		fprintf(stderr, "HK_FIXED_POSTERIOR_SOURCE_DIR and HK_FIXED_POSTERIOR_OUTPUT_DIR are required\n");
		exit(2);
	}
	if (p->d_scale_eps_count <= 0.0f || p->relax_steps < 0 || p->relax_step <= 0.0f) {
		fprintf(stderr, "invalid fixed-posterior numeric settings\n");
		exit(2);
	}
}

int main(void)
{
	struct fixed_params p;
	struct fixed_source src;
	struct hk_blind_wedge_list wedges;
	struct hk_fdg_conf fdg_conf;
	struct hk_blind_relax_diag relax_diag;
	struct hk_blind_contact_class_diag force_diag;
	struct sep_dist_stats sep_stats;
	char posterior_path[4096], coords_path[4096], coords_gz_path[4096];
	char loop_diag_path[4096], force_diag_path[4096], sep_diag_path[4096], manifest_path[4096];
	int ret = 1;

	load_params(&p);
	make_parent_dirs(p.output_dir);
	if (fixed_source_load(p.source_dir, &src) != 0)
		return 1;
	if (build_fixed_wedges(&wedges, &src.set, &p) != 0)
		goto cleanup_src;

	hk_fdg_conf_init(&fdg_conf);
	fdg_conf.backend = HK_FDG_BACKEND_CPU;
	fdg_conf.k_rel_rep *= p.repulsion_multiplier;
	if (hk_blind_relax_cpu(&fdg_conf, &wedges, &src.bmap, src.bmap.n_beads, src.coords,
						   HK_FIXED_UNIT, p.relax_step, p.relax_steps,
						   p.min_sep_unit, p.lambda_sep, 1, HK_BLIND_REPULSION_CELL,
						   &relax_diag) != 0) {
		fprintf(stderr, "fixed-posterior relax failed\n");
		goto cleanup_wedges;
	}

	compute_sep_dist_stats(&src.bmap, src.coords, p.min_sep_unit, p.lambda_sep, &sep_stats);
	path_join(posterior_path, sizeof(posterior_path), p.output_dir, "p9016_full.bpair_posterior.tsv");
	path_join(coords_path, sizeof(coords_path), p.output_dir, "p9016_full.coords.tsv");
	path_join(coords_gz_path, sizeof(coords_gz_path), p.output_dir, "p9016_full.coords.tsv.gz");
	path_join(loop_diag_path, sizeof(loop_diag_path), p.output_dir, "p9016_full.loop_diag.tsv");
	path_join(force_diag_path, sizeof(force_diag_path), p.output_dir, "p9016_full.force_class_diag.tsv");
	path_join(sep_diag_path, sizeof(sep_diag_path), p.output_dir, "p9016_full.sep_diag.tsv");
	path_join(manifest_path, sizeof(manifest_path), p.output_dir, "p9016_full.manifest.tsv");
	if (write_posterior_output(posterior_path, &src.bmap, &src.set) != 0 ||
		write_coords_outputs(coords_path, coords_gz_path, &src.bmap, src.coords) != 0 ||
		write_sep_diag(sep_diag_path, &sep_stats) != 0 ||
		write_loop_diag(loop_diag_path, &relax_diag, &sep_stats, &wedges, &src.set) != 0 ||
		write_force_class_diag(force_diag_path, &fdg_conf, &src.bmap, &wedges, src.coords,
							   &relax_diag, &force_diag) != 0 ||
		write_manifest(manifest_path, &p, &src, &wedges, &relax_diag, &sep_stats,
					   &force_diag, posterior_path, coords_path, coords_gz_path,
					   loop_diag_path, force_diag_path, sep_diag_path) != 0) {
		fprintf(stderr, "failed to write fixed-posterior outputs\n");
		goto cleanup_wedges;
	}
	fprintf(stderr, "[M::fixed-posterior] wrote %s graph=%s edges=%d\n",
			p.output_dir, graph_mode_name(p.graph_mode), wedges.n_edges);
	ret = 0;

cleanup_wedges:
	hk_blind_wedge_list_destroy(&wedges);
cleanup_src:
	fixed_source_destroy(&src);
	return ret;
}
