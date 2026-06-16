#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "hkpriv.h"
#include "fdg_gpu.h"

struct hk_blind_bpair_aux {
	struct hk_blind_bpair_key key;
	int32_t raw_id;
	uint8_t swapped;
};

struct hk_blind_base_k_aux {
	struct hk_blind_bpair_key key;
	int32_t max_nei;
};

struct hk_blind_softall_aux {
	struct hk_map *map;
	float *final_phased_prob;
	int32_t n_final_phased_prob;
	int32_t m_final_phased_prob;
	int64_t n_selected_raw;
	int64_t n_gate_skip_raw;
	int64_t n_same_bin_skip_raw;
	int64_t state_count[HK_BLIND_N_STATE];
};

static int32_t hk_blind_bmap_n_chr(const struct hk_bmap *bmap);

static int hk_blind_bpair_aux_cmp(const void *a_, const void *b_)
{
	const struct hk_blind_bpair_aux *a = (const struct hk_blind_bpair_aux*)a_;
	const struct hk_blind_bpair_aux *b = (const struct hk_blind_bpair_aux*)b_;
	if (a->key.bid[0] != b->key.bid[0])
		return a->key.bid[0] < b->key.bid[0]? -1 : 1;
	if (a->key.bid[1] != b->key.bid[1])
		return a->key.bid[1] < b->key.bid[1]? -1 : 1;
	if (a->raw_id != b->raw_id)
		return a->raw_id < b->raw_id? -1 : 1;
	return 0;
}

static int hk_blind_base_k_aux_cmp(const void *a_, const void *b_)
{
	const struct hk_blind_base_k_aux *a = (const struct hk_blind_base_k_aux*)a_;
	const struct hk_blind_base_k_aux *b = (const struct hk_blind_base_k_aux*)b_;
	if (a->key.bid[0] != b->key.bid[0])
		return a->key.bid[0] < b->key.bid[0]? -1 : 1;
	if (a->key.bid[1] != b->key.bid[1])
		return a->key.bid[1] < b->key.bid[1]? -1 : 1;
	return 0;
}

static void hk_blind_bpair_posterior_from_coords_score_mode_with_log_norm(const struct hk_fdg_conf *conf,
																		  int32_t bid0, int32_t bid1,
																		  const fvec3_t *coords,
																		  float unit, float d_scale, float k,
																		  const float log_prior[HK_BLIND_N_STATE],
																		  float temperature, int estep_score_mode,
																		  float energy_out[HK_BLIND_N_STATE],
																		  float p4[HK_BLIND_N_STATE],
																		  float *entropy, float *pmax, float *margin,
																		  float *rho_output, float *pU,
																		  float *real_log_norm);

static int hk_blind_bpair_key_eq(const struct hk_blind_bpair_key *a, const struct hk_blind_bpair_key *b)
{
	return a->bid[0] == b->bid[0] && a->bid[1] == b->bid[1];
}

static float hk_blind_coord_dist(const fvec3_t x, const fvec3_t y)
{
	float dx = x[0] - y[0];
	float dy = x[1] - y[1];
	float dz = x[2] - y[2];
	return sqrtf(dx * dx + dy * dy + dz * dz);
}

static double hk_blind_coord_dist2(const fvec3_t x, const fvec3_t y)
{
	double dx = (double)x[0] - y[0];
	double dy = (double)x[1] - y[1];
	double dz = (double)x[2] - y[2];
	return dx * dx + dy * dy + dz * dz;
}

static float hk_blind_bpair_base_d_scale_from_n_raw(int32_t n_raw)
{
	assert(n_raw > 0);
	return powf((float)n_raw, -1.0f / 3.0f);
}

const char *hk_blind_init_mode_name(int mode)
{
	switch (mode) {
	case HK_BLIND_INIT_TOY_SPLIT: return "toy_split";
	case HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT: return "unphased_scaffold_split";
	case HK_BLIND_INIT_RANDOM_DIPLOID: return "random_diploid";
	case HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT: return "random_haploid_split";
	default: return "unknown";
	}
}

const char *hk_blind_prior_mode_name(int mode)
{
	switch (mode) {
	case HK_BLIND_PRIOR_UNIFORM: return "uniform";
	case HK_BLIND_PRIOR_CIS_INTER_RATIO: return "cis_inter_ratio";
	default: return "unknown";
	}
}

const char *hk_blind_rho_train_mode_name(int mode)
{
	switch (mode) {
	case HK_BLIND_RHO_TRAIN_CONSTANT: return "constant";
	case HK_BLIND_RHO_TRAIN_ENTROPY: return "entropy";
	case HK_BLIND_RHO_TRAIN_ENTROPY_WITH_FLOOR: return "entropy_with_floor";
	case HK_BLIND_RHO_TRAIN_ENTROPY_CIS_CONSTANT_TRANS: return "entropy_cis_constant_trans";
	case HK_BLIND_RHO_TRAIN_ENTROPY_CIS_FLOOR_TRANS: return "entropy_cis_floor_trans";
	default: return "unknown";
	}
}

const char *hk_blind_d_scale_mode_name(int mode)
{
	switch (mode) {
	case HK_BLIND_D_SCALE_RAW_COUNT: return "raw_count";
	case HK_BLIND_D_SCALE_EXPECTED_COUNT: return "expected_count";
	case HK_BLIND_D_SCALE_DENSITY_NORMALIZED_RAW_COUNT: return "density_normalized_raw_count";
	case HK_BLIND_D_SCALE_CAPPED_DENSITY_RAW_COUNT: return "capped_density_raw_count";
	default: return "unknown";
	}
}

const char *hk_blind_estep_score_mode_name(int mode)
{
	switch (mode) {
	case HK_BLIND_ESTEP_SCORE_FDG_FLAT: return "fdg_flat";
	case HK_BLIND_ESTEP_SCORE_LOGDIST2: return "logdist2";
	case HK_BLIND_ESTEP_SCORE_DIST2: return "dist2";
	default: return "unknown";
	}
}

const char *hk_blind_base_k_mode_name(int mode)
{
	switch (mode) {
	case HK_BLIND_BASE_K_UNIFORM: return "uniform";
	case HK_BLIND_BASE_K_NEIGHBOR_MEDIAN: return "neighbor_median";
	default: return "unknown";
	}
}

const char *hk_blind_contact_class_name(int contact_class)
{
	switch (contact_class) {
	case HK_BLIND_CONTACT_CIS: return "cis";
	case HK_BLIND_CONTACT_TRANS: return "trans";
	default: return "unknown";
	}
}

int hk_blind_init_mode_valid(int mode)
{
	return mode == HK_BLIND_INIT_TOY_SPLIT ||
		   mode == HK_BLIND_INIT_UNPHASED_SCAFFOLD_SPLIT ||
		   mode == HK_BLIND_INIT_RANDOM_DIPLOID ||
		   mode == HK_BLIND_INIT_RANDOM_HAPLOID_SPLIT;
}

int hk_blind_prior_mode_valid(int mode)
{
	return mode == HK_BLIND_PRIOR_UNIFORM ||
		   mode == HK_BLIND_PRIOR_CIS_INTER_RATIO;
}

int hk_blind_rho_train_mode_valid(int mode)
{
	return mode == HK_BLIND_RHO_TRAIN_CONSTANT ||
		   mode == HK_BLIND_RHO_TRAIN_ENTROPY ||
		   mode == HK_BLIND_RHO_TRAIN_ENTROPY_WITH_FLOOR ||
		   mode == HK_BLIND_RHO_TRAIN_ENTROPY_CIS_CONSTANT_TRANS ||
		   mode == HK_BLIND_RHO_TRAIN_ENTROPY_CIS_FLOOR_TRANS;
}

int hk_blind_d_scale_mode_valid(int mode)
{
	return mode == HK_BLIND_D_SCALE_RAW_COUNT ||
		mode == HK_BLIND_D_SCALE_EXPECTED_COUNT ||
		mode == HK_BLIND_D_SCALE_DENSITY_NORMALIZED_RAW_COUNT ||
		mode == HK_BLIND_D_SCALE_CAPPED_DENSITY_RAW_COUNT;
}

int hk_blind_estep_score_mode_valid(int mode)
{
	return mode == HK_BLIND_ESTEP_SCORE_FDG_FLAT ||
		   mode == HK_BLIND_ESTEP_SCORE_LOGDIST2 ||
		   mode == HK_BLIND_ESTEP_SCORE_DIST2;
}

int hk_blind_base_k_mode_valid(int mode)
{
	return mode == HK_BLIND_BASE_K_UNIFORM ||
		   mode == HK_BLIND_BASE_K_NEIGHBOR_MEDIAN;
}

int hk_blind_contact_class_valid(int contact_class)
{
	return contact_class == HK_BLIND_CONTACT_CIS ||
		   contact_class == HK_BLIND_CONTACT_TRANS;
}

static int hk_blind_bpair_contact_class_or_default(const struct hk_blind_bpair *bp)
{
	assert(bp);
	return hk_blind_contact_class_valid(bp->contact_class)?
		bp->contact_class : HK_BLIND_CONTACT_CIS;
}

static int hk_blind_float_isfinite(float x)
{
	union { float f; uint32_t u; } v;
	v.f = x;
	return (v.u & 0x7f800000u) != 0x7f800000u;
}

static float hk_blind_clip01(float x)
{
	if (x < 0.0f) return 0.0f;
	if (x > 1.0f) return 1.0f;
	return x;
}

static float hk_blind_logsumexp4(const float score[HK_BLIND_N_STATE])
{
	float max_score, sum_exp = 0.0f;
	int i;
	assert(score);
	max_score = score[0];
	for (i = 1; i < HK_BLIND_N_STATE; ++i)
		if (score[i] > max_score)
			max_score = score[i];
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		sum_exp += expf(score[i] - max_score);
	return max_score + logf(sum_exp);
}

static void hk_blind_p4_uncertainty(const float p4[HK_BLIND_N_STATE], float *entropy, float *pmax, float *margin,
									float *rho_output, float *pU)
{
	float top1 = 0.0f, top2 = 0.0f;
	int i;

	assert(p4);
	assert(entropy);
	assert(pmax);
	assert(margin);
	assert(rho_output);
	assert(pU);

	*entropy = 0.0f;
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		assert(isfinite(p4[i]));
		assert(p4[i] >= 0.0f);
		if (p4[i] > 0.0f)
			*entropy -= p4[i] * logf(p4[i]);
		if (p4[i] > top1) {
			top2 = top1;
			top1 = p4[i];
		} else if (p4[i] > top2) {
			top2 = p4[i];
		}
	}
	*pmax = top1;
	*margin = top1 - top2;
	*rho_output = hk_blind_clip01(1.0f - *entropy / logf((float)HK_BLIND_N_STATE));
	*pU = 1.0f - *rho_output;
}

static void hk_blind_bpair_set_unknown_posterior(struct hk_blind_bpair *p)
{
	int i;
	assert(p);
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		p->p4[i] = 1.0f / HK_BLIND_N_STATE;
	p->entropy = logf((float)HK_BLIND_N_STATE);
	p->pmax = 1.0f / HK_BLIND_N_STATE;
	p->margin = 0.0f;
	p->rho_output = 0.0f;
	p->pU = 1.0f;
	p->real_log_norm = 0.0f;
}

static uint64_t hk_blind_chr_seed(uint64_t seed, int32_t chr)
{
	return kr_splitmix64(seed ^ 0x6a09e667f3bcc909ULL ^ ((uint64_t)(uint32_t)chr + 1) * 0x9e3779b97f4a7c15ULL);
}

static void hk_blind_chr_unit_vector(uint64_t seed, int32_t chr, float v[3])
{
	krng_t rng;
	float z, theta, xy;

	assert(chr >= 0);
	assert(v);

	kr_srand_r(&rng, hk_blind_chr_seed(seed, chr));
	z = (float)(2.0 * kr_drand_r(&rng) - 1.0);
	theta = (float)(6.2831853071795864769 * kr_drand_r(&rng));
	xy = sqrtf(fmaxf(0.0f, 1.0f - z * z));
	v[0] = xy * cosf(theta);
	v[1] = xy * sinf(theta);
	v[2] = z;
}

static float hk_blind_noise_coord(krng_t *rng, float noise_scale)
{
	return noise_scale * (float)(2.0 * kr_drand_r(rng) - 1.0);
}

static void hk_blind_bpair_init(struct hk_blind_bpair *p, const struct hk_blind_bpair_key *key)
{
	int i;
	assert(p);
	assert(key);
	p->key = *key;
	p->n_raw = 0;
	p->base_d_scale = 1.0f;
	p->base_k = 1.0f;
	p->contact_class = HK_BLIND_CONTACT_CIS;
	p->density_child_count[0] = 1;
	p->density_child_count[1] = 1;
	p->density_exposure = 1.0f;
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		p->log_prior[i] = -logf((float)HK_BLIND_N_STATE);
	hk_blind_bpair_set_unknown_posterior(p);
}

static int hk_blind_bpair_key_contact_class(const struct hk_bmap *bmap, const struct hk_blind_bpair_key *key)
{
	const struct hk_bead *b0, *b1;
	assert(bmap);
	assert(key);
	assert(key->bid[0] >= 0 && key->bid[0] < bmap->n_beads);
	assert(key->bid[1] >= 0 && key->bid[1] < bmap->n_beads);
	assert(bmap->beads);
	b0 = &bmap->beads[key->bid[0]];
	b1 = &bmap->beads[key->bid[1]];
	return b0->chr == b1->chr? HK_BLIND_CONTACT_CIS : HK_BLIND_CONTACT_TRANS;
}

void hk_blind_init_uniform_log_prior(float log_prior[HK_BLIND_N_STATE])
{
	int i;
	assert(log_prior);
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		log_prior[i] = -logf((float)HK_BLIND_N_STATE);
}

void hk_blind_bpair_set_init_uniform_prior(struct hk_blind_bpair_set *set)
{
	int32_t i;
	assert(set);
	assert(set->n_bpairs >= 0);
	assert(set->n_bpairs == 0 || set->bpairs);
	for (i = 0; i < set->n_bpairs; ++i)
		hk_blind_init_uniform_log_prior(set->bpairs[i].log_prior);
}

static int hk_blind_float_cmp(const void *a_, const void *b_)
{
	float a = *(const float*)a_;
	float b = *(const float*)b_;
	return (a > b) - (a < b);
}

static void hk_blind_set_same_cross_prior(struct hk_blind_bpair *bp, float alpha_cross)
{
	float same, cross;
	assert(bp);
	assert(isfinite(alpha_cross));
	assert(alpha_cross > 0.0f && alpha_cross <= 0.5f);
	same = 0.5f * (1.0f - alpha_cross);
	cross = 0.5f * alpha_cross;
	bp->log_prior[HK_BLIND_STATE_00] = logf(same);
	bp->log_prior[HK_BLIND_STATE_11] = logf(same);
	bp->log_prior[HK_BLIND_STATE_01] = logf(cross);
	bp->log_prior[HK_BLIND_STATE_10] = logf(cross);
}

static double hk_blind_density_smooth3(const double *obs, const double *possible, int32_t max_delta, int32_t delta)
{
	double obs_sum = 0.0, possible_sum = 0.0;
	int32_t d;
	assert(obs);
	assert(possible);
	assert(delta >= 0 && delta <= max_delta);
	for (d = delta - 1; d <= delta + 1; ++d) {
		if (d < 1 || d > max_delta) continue;
		obs_sum += obs[d];
		possible_sum += possible[d];
	}
	return possible_sum > 0.0? obs_sum / possible_sum : 0.0;
}

int hk_blind_bpair_set_init_cis_inter_ratio_prior(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
												  float eps_prior, struct hk_blind_prior_diag *diag)
{
	double *obs_cis = 0, *possible_cis = 0;
	float *alpha_values = 0;
	double observed_inter = 0.0, possible_inter = 0.0, inter_density;
	int32_t max_delta = 0, n_alpha = 0;
	int32_t i, c0, c1;
	int ret = -1;

	assert(bmap);
	assert(set);
	assert(set->n_bpairs >= 0);
	assert(set->n_bpairs == 0 || set->bpairs);
	assert(bmap->n_beads >= 0);
	assert(bmap->n_beads == 0 || bmap->beads);
	assert(bmap->d);
	assert(bmap->offcnt);
	if (!isfinite(eps_prior) || eps_prior <= 0.0f)
		eps_prior = 1e-6f;
	if (eps_prior > 0.5f)
		eps_prior = 0.5f;
	for (c0 = 0; c0 < bmap->d->n; ++c0) {
		int32_t cnt0 = (int32_t)bmap->offcnt[c0];
		if (cnt0 - 1 > max_delta)
			max_delta = cnt0 - 1;
		for (c1 = c0 + 1; c1 < bmap->d->n; ++c1) {
			int32_t cnt1 = (int32_t)bmap->offcnt[c1];
			possible_inter += (double)cnt0 * (double)cnt1;
		}
	}
	obs_cis = CALLOC(double, (size_t)max_delta + 1);
	possible_cis = CALLOC(double, (size_t)max_delta + 1);
	alpha_values = set->n_bpairs > 0? MALLOC(float, set->n_bpairs) : 0;
	if ((max_delta > 0 && (obs_cis == 0 || possible_cis == 0)) ||
		(set->n_bpairs > 0 && alpha_values == 0))
		goto cleanup;
	for (c0 = 0; c0 < bmap->d->n; ++c0) {
		int32_t cnt = (int32_t)bmap->offcnt[c0];
		int32_t d;
		for (d = 1; d < cnt; ++d)
			possible_cis[d] += (double)(cnt - d);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		struct hk_blind_bpair *bp = &set->bpairs[i];
		int32_t bid0 = bp->key.bid[0], bid1 = bp->key.bid[1];
		const struct hk_bead *b0, *b1;
		assert(bid0 >= 0 && bid0 < bmap->n_beads);
		assert(bid1 >= 0 && bid1 < bmap->n_beads);
		b0 = &bmap->beads[bid0];
		b1 = &bmap->beads[bid1];
		if (b0->chr != b1->chr) {
			observed_inter += bp->n_raw;
		} else if (bid0 != bid1) {
			int32_t delta = bid1 > bid0? bid1 - bid0 : bid0 - bid1;
			if (delta >= 1 && delta <= max_delta)
				obs_cis[delta] += bp->n_raw;
		}
	}
	inter_density = possible_inter > 0.0? observed_inter / possible_inter : 0.0;
	hk_blind_bpair_set_init_uniform_prior(set);
	for (i = 0; i < set->n_bpairs; ++i) {
		struct hk_blind_bpair *bp = &set->bpairs[i];
		int32_t bid0 = bp->key.bid[0], bid1 = bp->key.bid[1];
		const struct hk_bead *b0 = &bmap->beads[bid0];
		const struct hk_bead *b1 = &bmap->beads[bid1];
		if (b0->chr == b1->chr && bid0 != bid1) {
			int32_t delta = bid1 > bid0? bid1 - bid0 : bid0 - bid1;
			double cis_density = hk_blind_density_smooth3(obs_cis, possible_cis, max_delta, delta);
			double denom = cis_density > eps_prior? cis_density : eps_prior;
			float alpha = (float)(inter_density / denom);
			if (!isfinite(alpha)) alpha = 0.5f;
			if (alpha < eps_prior) alpha = eps_prior;
			if (alpha > 0.5f) alpha = 0.5f;
			hk_blind_set_same_cross_prior(bp, alpha);
			alpha_values[n_alpha++] = alpha;
		}
	}
	if (diag) {
		memset(diag, 0, sizeof(*diag));
		diag->prior_mode = HK_BLIND_PRIOR_CIS_INTER_RATIO;
		diag->n_distance_bins = max_delta;
		diag->n_alpha = n_alpha;
		diag->observed_inter = observed_inter;
		diag->possible_inter = possible_inter;
		diag->inter_density = inter_density;
		diag->eps_prior = eps_prior;
		if (n_alpha > 0) {
			qsort(alpha_values, n_alpha, sizeof(*alpha_values), hk_blind_float_cmp);
			diag->alpha_min = alpha_values[0];
			diag->alpha_median = alpha_values[n_alpha / 2];
			diag->alpha_max = alpha_values[n_alpha - 1];
		}
	}
	ret = 0;

cleanup:
	free(obs_cis);
	free(possible_cis);
	free(alpha_values);
	return ret;
}

void hk_blind_bpair_set_base_k_stats(const struct hk_blind_bpair_set *set, struct hk_blind_base_k_stats *stats)
{
	double sum = 0.0;
	int32_t i;
	assert(set);
	assert(stats);
	memset(stats, 0, sizeof(*stats));
	stats->min = INFINITY;
	stats->max = 0.0f;
	for (i = 0; i < set->n_bpairs; ++i) {
		float k = set->bpairs[i].base_k;
		if (!isfinite(k)) {
			++stats->n_nonfinite;
			continue;
		}
		if (stats->n == 0 || k < stats->min) stats->min = k;
		if (stats->n == 0 || k > stats->max) stats->max = k;
		sum += k;
		++stats->n;
	}
	if (stats->n > 0)
		stats->mean = (float)(sum / stats->n);
	else
		stats->min = stats->max = stats->mean = 0.0f;
}

static void hk_blind_posterior_from_energy_with_log_norm(const float energy[HK_BLIND_N_STATE],
														 const float log_prior[HK_BLIND_N_STATE],
														 float temperature,
														 float p4[HK_BLIND_N_STATE],
														 float *entropy, float *pmax,
														 float *margin, float *rho_output,
														 float *pU, float *real_log_norm)
{
	float score[HK_BLIND_N_STATE], max_score, sum_exp = 0.0f;
	int i;

	assert(energy);
	assert(log_prior);
	assert(p4);
	assert(entropy);
	assert(pmax);
	assert(margin);
	assert(rho_output);
	assert(pU);
	assert(temperature > 0.0f);

	max_score = (-energy[0] + log_prior[0]) / temperature;
	score[0] = max_score;
	for (i = 1; i < HK_BLIND_N_STATE; ++i) {
		score[i] = (-energy[i] + log_prior[i]) / temperature;
		if (score[i] > max_score) max_score = score[i];
	}
	if (real_log_norm)
		*real_log_norm = hk_blind_logsumexp4(score);
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		p4[i] = expf(score[i] - max_score);
		sum_exp += p4[i];
	}
	assert(sum_exp > 0.0f);

	*entropy = 0.0f;
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		p4[i] /= sum_exp;
	}
	hk_blind_p4_uncertainty(p4, entropy, pmax, margin, rho_output, pU);
}

void hk_blind_posterior_from_energy(const float energy[HK_BLIND_N_STATE], const float log_prior[HK_BLIND_N_STATE], float temperature,
									float p4[HK_BLIND_N_STATE], float *entropy, float *pmax, float *margin, float *rho_output, float *pU)
{
	hk_blind_posterior_from_energy_with_log_norm(energy, log_prior, temperature,
												 p4, entropy, pmax, margin,
												 rho_output, pU, 0);
}

void hk_blind_bpair_posterior_from_coords_score_mode(const struct hk_fdg_conf *conf,
													 int32_t bid0, int32_t bid1,
													 const fvec3_t *coords,
													 float unit, float d_scale, float k,
													 const float log_prior[HK_BLIND_N_STATE],
													 float temperature, int estep_score_mode,
													 float energy_out[HK_BLIND_N_STATE],
													 float p4[HK_BLIND_N_STATE],
													 float *entropy, float *pmax, float *margin,
													 float *rho_output, float *pU)
{
	hk_blind_bpair_posterior_from_coords_score_mode_with_log_norm(
		conf, bid0, bid1, coords, unit, d_scale, k, log_prior, temperature,
		estep_score_mode, energy_out, p4, entropy, pmax, margin, rho_output,
		pU, 0);
}

static void hk_blind_bpair_posterior_from_coords_score_mode_with_log_norm(const struct hk_fdg_conf *conf,
																		  int32_t bid0, int32_t bid1,
																		  const fvec3_t *coords,
																		  float unit, float d_scale, float k,
																		  const float log_prior[HK_BLIND_N_STATE],
																		  float temperature, int estep_score_mode,
																		  float energy_out[HK_BLIND_N_STATE],
																		  float p4[HK_BLIND_N_STATE],
																		  float *entropy, float *pmax, float *margin,
																		  float *rho_output, float *pU,
																		  float *real_log_norm)
{
	float energy[HK_BLIND_N_STATE];
	int32_t i0, i1, j0, j1;
	float dist[HK_BLIND_N_STATE];
	int i;

	assert(conf);
	assert(bid0 >= 0);
	assert(bid1 >= 0);
	assert(coords);
	assert(unit > 0.0f);
	assert(d_scale > 0.0f);
	assert(k >= 0.0f);
	assert(log_prior);
	assert(p4);
	assert(hk_blind_estep_score_mode_valid(estep_score_mode));

	i0 = hk_diploid_bid(bid0, HK_DIPLOID_COPY0);
	i1 = hk_diploid_bid(bid0, HK_DIPLOID_COPY1);
	j0 = hk_diploid_bid(bid1, HK_DIPLOID_COPY0);
	j1 = hk_diploid_bid(bid1, HK_DIPLOID_COPY1);

	dist[HK_BLIND_STATE_00] = hk_blind_coord_dist(coords[i0], coords[j0]);
	dist[HK_BLIND_STATE_01] = hk_blind_coord_dist(coords[i0], coords[j1]);
	dist[HK_BLIND_STATE_10] = hk_blind_coord_dist(coords[i1], coords[j0]);
	dist[HK_BLIND_STATE_11] = hk_blind_coord_dist(coords[i1], coords[j1]);
	if (estep_score_mode == HK_BLIND_ESTEP_SCORE_FDG_FLAT) {
		for (i = 0; i < HK_BLIND_N_STATE; ++i)
			energy[i] = hk_fdg_contact_energy_dist(conf, dist[i], unit, d_scale, k);
	} else if (estep_score_mode == HK_BLIND_ESTEP_SCORE_LOGDIST2) {
		const float eps = 1.0e-6f;
		for (i = 0; i < HK_BLIND_N_STATE; ++i) {
			float r = dist[i] / (unit * d_scale);
			energy[i] = k * log1pf((r * r) / eps);
		}
	} else {
		for (i = 0; i < HK_BLIND_N_STATE; ++i) {
			float r = dist[i] / (unit * d_scale);
			energy[i] = k * r * r;
		}
	}

	if (energy_out)
		for (i = 0; i < HK_BLIND_N_STATE; ++i)
			energy_out[i] = energy[i];
	hk_blind_posterior_from_energy_with_log_norm(energy, log_prior, temperature,
												 p4, entropy, pmax, margin,
												 rho_output, pU, real_log_norm);
}

void hk_blind_bpair_posterior_from_coords(const struct hk_fdg_conf *conf, int32_t bid0, int32_t bid1, const fvec3_t *coords,
										  float unit, float d_scale, float k, const float log_prior[HK_BLIND_N_STATE], float temperature,
										  float energy_out[HK_BLIND_N_STATE], float p4[HK_BLIND_N_STATE],
										  float *entropy, float *pmax, float *margin, float *rho_output, float *pU)
{
	hk_blind_bpair_posterior_from_coords_score_mode(conf, bid0, bid1, coords, unit, d_scale, k,
													log_prior, temperature,
													HK_BLIND_ESTEP_SCORE_FDG_FLAT,
													energy_out, p4, entropy, pmax, margin,
													rho_output, pU);
}

void hk_blind_bpair_set_update_posterior_from_coords(struct hk_blind_bpair_set *set, const struct hk_fdg_conf *conf,
													 const fvec3_t *coords, float unit, float d_scale, float k,
													 const float log_prior[HK_BLIND_N_STATE], float temperature)
{
	int32_t i;

	assert(set);
	assert(set->n_bpairs >= 0);
	assert(set->n_bpairs == 0 || set->bpairs);
	assert(conf);
	assert(coords);
	assert(unit > 0.0f);
	assert(d_scale > 0.0f);
	assert(k >= 0.0f);
	assert(log_prior);
	assert(temperature > 0.0f);

	for (i = 0; i < set->n_bpairs; ++i) {
		struct hk_blind_bpair *p = &set->bpairs[i];
		if (hk_blind_bpair_is_same_bin(p)) {
			hk_blind_bpair_set_unknown_posterior(p);
			continue;
		}
		hk_blind_bpair_posterior_from_coords_score_mode_with_log_norm(
			conf, p->key.bid[0], p->key.bid[1], coords,
			unit, d_scale, k, log_prior, temperature,
			HK_BLIND_ESTEP_SCORE_FDG_FLAT, 0, p->p4, &p->entropy,
			&p->pmax, &p->margin, &p->rho_output, &p->pU,
			&p->real_log_norm);
	}
}

void hk_blind_bpair_set_update_posterior_from_coords_params(struct hk_blind_bpair_set *set, const struct hk_fdg_conf *conf,
															const fvec3_t *coords, float unit,
															const float log_prior[HK_BLIND_N_STATE], float temperature)
{
	hk_blind_bpair_set_update_posterior_from_coords_params_score_mode(set, conf, coords, unit,
																	  log_prior, temperature,
																	  HK_BLIND_ESTEP_SCORE_FDG_FLAT);
}

void hk_blind_bpair_set_update_posterior_from_coords_params_score_mode(struct hk_blind_bpair_set *set,
																	   const struct hk_fdg_conf *conf,
																	   const fvec3_t *coords, float unit,
																	   const float log_prior[HK_BLIND_N_STATE],
																	   float temperature,
																	   int estep_score_mode)
{
	int32_t i;

	assert(set);
	assert(set->n_bpairs >= 0);
	assert(set->n_bpairs == 0 || set->bpairs);
	assert(conf);
	assert(set->n_bpairs == 0 || coords);
	assert(unit > 0.0f);
	assert(temperature > 0.0f);
	assert(hk_blind_estep_score_mode_valid(estep_score_mode));

	for (i = 0; i < set->n_bpairs; ++i) {
		struct hk_blind_bpair *p = &set->bpairs[i];
		int s;
		assert(isfinite(p->base_d_scale));
		assert(p->base_d_scale > 0.0f);
		assert(isfinite(p->base_k));
		assert(p->base_k >= 0.0f);
		if (hk_blind_bpair_is_same_bin(p)) {
			hk_blind_bpair_set_unknown_posterior(p);
			continue;
		}
		if (log_prior == 0) {
			for (s = 0; s < HK_BLIND_N_STATE; ++s)
				assert(isfinite(p->log_prior[s]));
		}
		hk_blind_bpair_posterior_from_coords_score_mode_with_log_norm(
			conf, p->key.bid[0], p->key.bid[1], coords,
			unit, p->base_d_scale, p->base_k,
			log_prior? log_prior : p->log_prior, temperature,
			estep_score_mode, 0, p->p4, &p->entropy,
			&p->pmax, &p->margin, &p->rho_output, &p->pU,
			&p->real_log_norm);
		for (s = 0; s < HK_BLIND_N_STATE; ++s)
			assert(isfinite(p->p4[s]));
		assert(isfinite(p->entropy));
		assert(isfinite(p->pmax));
		assert(isfinite(p->margin));
		assert(isfinite(p->rho_output));
		assert(isfinite(p->pU));
		assert(isfinite(p->real_log_norm));
	}
}

int hk_blind_init_diploid_coords_from_haploid(const struct hk_bmap *bmap, const fvec3_t *haploid_x, int32_t n_haploid,
											 fvec3_t *diploid_x, float eps, float noise_scale, uint64_t seed)
{
	krng_t noise_rng;
	int32_t i;

	assert(bmap);
	assert(n_haploid >= 0);
	assert(bmap->n_beads == n_haploid);
	assert(n_haploid == 0 || bmap->beads);
	assert(n_haploid == 0 || haploid_x);
	assert(n_haploid == 0 || diploid_x);
	assert(isfinite(eps));
	assert(eps >= 0.0f);
	assert(isfinite(noise_scale));
	assert(noise_scale >= 0.0f);

	kr_srand_r(&noise_rng, kr_splitmix64(seed ^ 0xbb67ae8584caa73bULL));
	for (i = 0; i < n_haploid; ++i) {
		float v[3], n0[3] = {0.0f, 0.0f, 0.0f}, n1[3] = {0.0f, 0.0f, 0.0f};
		int32_t d0, d1;
		int a;

		assert(bmap->beads[i].chr >= 0);
		for (a = 0; a < 3; ++a)
			assert(isfinite(haploid_x[i][a]));

		hk_blind_chr_unit_vector(seed, bmap->beads[i].chr, v);
		if (noise_scale > 0.0f) {
			for (a = 0; a < 3; ++a) n0[a] = hk_blind_noise_coord(&noise_rng, noise_scale);
			for (a = 0; a < 3; ++a) n1[a] = hk_blind_noise_coord(&noise_rng, noise_scale);
		}

		d0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		d1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		for (a = 0; a < 3; ++a) {
			diploid_x[d0][a] = haploid_x[i][a] + eps * v[a] + n0[a];
			diploid_x[d1][a] = haploid_x[i][a] - eps * v[a] + n1[a];
			assert(isfinite(diploid_x[d0][a]));
			assert(isfinite(diploid_x[d1][a]));
		}
	}
	return 0;
}

int hk_blind_init_toy_haploid_scaffold(const struct hk_bmap *bmap, fvec3_t *haploid_x)
{
	int32_t i;
	assert(bmap);
	assert(bmap->n_beads >= 0);
	assert(bmap->n_beads == 0 || haploid_x);
	for (i = 0; i < bmap->n_beads; ++i) {
		haploid_x[i][0] = 0.10f * (float)(i % 97);
		haploid_x[i][1] = 0.07f * (float)((i / 97) % 97);
		haploid_x[i][2] = 0.03f * (float)(i % 17);
	}
	return 0;
}

int hk_blind_init_random_haploid_scaffold(const struct hk_bmap *bmap, fvec3_t *haploid_x, float scale, uint64_t seed)
{
	krng_t rng;
	int32_t i;
	int a;
	assert(bmap);
	assert(bmap->n_beads >= 0);
	assert(bmap->n_beads == 0 || haploid_x);
	assert(isfinite(scale));
	assert(scale >= 0.0f);
	kr_srand_r(&rng, kr_splitmix64(seed ^ 0x510e527fade682d1ULL));
	for (i = 0; i < bmap->n_beads; ++i) {
		for (a = 0; a < 3; ++a) {
			haploid_x[i][a] = scale * (float)(2.0 * kr_drand_r(&rng) - 1.0);
			assert(isfinite(haploid_x[i][a]));
		}
	}
	return 0;
}

int hk_blind_init_random_diploid_coords(const struct hk_bmap *bmap, fvec3_t *diploid_x, float scale, uint64_t seed)
{
	krng_t rng;
	int32_t n_diploid, i;
	int a;
	assert(bmap);
	assert(bmap->n_beads >= 0);
	assert(bmap->n_beads <= INT32_MAX / HK_DIPLOID_N_COPY);
	assert(bmap->n_beads == 0 || diploid_x);
	assert(isfinite(scale));
	assert(scale >= 0.0f);
	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	kr_srand_r(&rng, kr_splitmix64(seed ^ 0x1f83d9abfb41bd6bULL));
	for (i = 0; i < n_diploid; ++i) {
		for (a = 0; a < 3; ++a) {
			diploid_x[i][a] = scale * (float)(2.0 * kr_drand_r(&rng) - 1.0);
			assert(isfinite(diploid_x[i][a]));
		}
	}
	return 0;
}

int hk_blind_init_haploid_scaffold_from_bmap_fdg(struct hk_bmap *bmap, const struct hk_fdg_conf *fdg_conf,
												 fvec3_t *haploid_x, uint64_t seed)
{
	struct hk_fdg_conf local_conf;
	krng_t rng;
	int32_t i;
	int a;

	assert(bmap);
	assert(bmap->n_beads >= 0);
	assert(bmap->n_beads == 0 || haploid_x);
	if (fdg_conf == 0) {
		hk_fdg_conf_init(&local_conf);
		fdg_conf = &local_conf;
	}
	kr_srand_r(&rng, kr_splitmix64(seed ^ 0x5be0cd19137e2179ULL));
	if (bmap->x) {
		free(bmap->x);
		bmap->x = 0;
	}
	hk_fdg(fdg_conf, bmap, 0, &rng);
	if (bmap->n_beads > 0 && bmap->x == 0)
		return -1;
	for (i = 0; i < bmap->n_beads; ++i) {
		for (a = 0; a < 3; ++a) {
			haploid_x[i][a] = bmap->x[i][a];
			if (!isfinite(haploid_x[i][a]))
				return -1;
		}
	}
	return 0;
}

float hk_blind_homolog_sep_energy(const fvec3_t x0, const fvec3_t x1, float unit, float min_sep_unit, float lambda_sep)
{
	float d, d_hat, delta;
	int a;

	assert(x0);
	assert(x1);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(isfinite(min_sep_unit));
	assert(min_sep_unit >= 0.0f);
	assert(isfinite(lambda_sep));
	assert(lambda_sep >= 0.0f);
	for (a = 0; a < 3; ++a) {
		assert(isfinite(x0[a]));
		assert(isfinite(x1[a]));
	}

	d = hk_blind_coord_dist(x0, x1);
	d_hat = d / unit;
	if (d_hat >= min_sep_unit || lambda_sep == 0.0f)
		return 0.0f;
	delta = min_sep_unit - d_hat;
	return lambda_sep * delta * delta;
}

float hk_blind_homolog_sep_energy_force(const fvec3_t x0, const fvec3_t x1, float unit, float min_sep_unit, float lambda_sep,
										fvec3_t f0, fvec3_t f1)
{
	const float eps = 1e-20f;
	float dx, dy, dz, d, d_hat, delta, force_mag, dir[3];
	int a;

	assert(x0);
	assert(x1);
	assert(f0);
	assert(f1);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(isfinite(min_sep_unit));
	assert(min_sep_unit >= 0.0f);
	assert(isfinite(lambda_sep));
	assert(lambda_sep >= 0.0f);

	for (a = 0; a < 3; ++a) {
		assert(isfinite(x0[a]));
		assert(isfinite(x1[a]));
		f0[a] = 0.0f;
		f1[a] = 0.0f;
	}

	dx = x0[0] - x1[0];
	dy = x0[1] - x1[1];
	dz = x0[2] - x1[2];
	d = sqrtf(dx * dx + dy * dy + dz * dz);
	d_hat = d / unit;
	if (d_hat >= min_sep_unit || lambda_sep == 0.0f)
		return 0.0f;

	delta = min_sep_unit - d_hat;
	force_mag = 2.0f * lambda_sep * delta / unit;
	if (d > eps) {
		dir[0] = dx / d;
		dir[1] = dy / d;
		dir[2] = dz / d;
	} else {
		// Deterministic fallback for overlapping or numerically indistinguishable homolog copies.
		dir[0] = 1.0f;
		dir[1] = 0.0f;
		dir[2] = 0.0f;
	}
	for (a = 0; a < 3; ++a) {
		f0[a] = force_mag * dir[a];
		f1[a] = -f0[a];
		assert(isfinite(f0[a]));
		assert(isfinite(f1[a]));
	}
	return lambda_sep * delta * delta;
}

float hk_blind_homolog_sep_accumulate_force_ex(int32_t n_haploid, const fvec3_t *coords, fvec3_t *force,
											  float unit, float min_sep_unit, float lambda_sep,
											  float *force_l1, int32_t *n_nonfinite)
{
	double total = 0.0;
	double l1 = 0.0;
	int32_t n_bad = 0;
	int32_t i;

	assert(n_haploid >= 0);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(isfinite(min_sep_unit));
	assert(min_sep_unit >= 0.0f);
	assert(isfinite(lambda_sep));
	assert(lambda_sep >= 0.0f);
	if (force_l1) *force_l1 = 0.0f;
	if (n_nonfinite) *n_nonfinite = 0;
	if (n_haploid == 0)
		return 0.0f;
	assert(coords);
	assert(force);

	for (i = 0; i < n_haploid; ++i) {
		int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		fvec3_t f0, f1;
		float e;
		int a;

		e = hk_blind_homolog_sep_energy_force(coords[b0], coords[b1], unit, min_sep_unit, lambda_sep, f0, f1);
		total += e;
		for (a = 0; a < 3; ++a) {
			force[b0][a] += f0[a];
			force[b1][a] += f1[a];
			if (isfinite(f0[a])) l1 += fabs((double)f0[a]);
			else ++n_bad;
			if (isfinite(f1[a])) l1 += fabs((double)f1[a]);
			else ++n_bad;
			assert(isfinite(force[b0][a]));
			assert(isfinite(force[b1][a]));
		}
	}
	if (force_l1) *force_l1 = (float)l1;
	if (n_nonfinite) *n_nonfinite = n_bad;
	return (float)total;
}

float hk_blind_homolog_sep_accumulate_force(int32_t n_haploid, const fvec3_t *coords, fvec3_t *force,
											float unit, float min_sep_unit, float lambda_sep)
{
	return hk_blind_homolog_sep_accumulate_force_ex(n_haploid, coords, force, unit,
													min_sep_unit, lambda_sep, 0, 0);
}

static float hk_blind_chr_centroid_sep_accumulate_force(const struct hk_bmap *bmap,
														const fvec3_t *coords,
														fvec3_t *force,
														float unit,
														float min_sep_unit,
														float lambda_sep,
														float *force_l1,
														int32_t *n_nonfinite)
{
	int32_t n_chr, c, i;
	int32_t *counts = 0;
	fvec3_t *centroid0 = 0, *centroid1 = 0;
	double total = 0.0, l1 = 0.0;
	int32_t n_bad = 0;

	assert(bmap);
	assert(bmap->n_beads >= 0);
	assert(bmap->n_beads == 0 || bmap->beads);
	assert(coords);
	assert(force);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(isfinite(min_sep_unit));
	assert(min_sep_unit >= 0.0f);
	assert(isfinite(lambda_sep));
	assert(lambda_sep >= 0.0f);
	if (force_l1) *force_l1 = 0.0f;
	if (n_nonfinite) *n_nonfinite = 0;
	if (bmap->n_beads == 0 || min_sep_unit <= 0.0f || lambda_sep <= 0.0f)
		return 0.0f;
	n_chr = hk_blind_bmap_n_chr(bmap);
	if (n_chr <= 0)
		return 0.0f;
	counts = CALLOC(int32_t, n_chr);
	centroid0 = CALLOC(fvec3_t, n_chr);
	centroid1 = CALLOC(fvec3_t, n_chr);
	if (counts == 0 || centroid0 == 0 || centroid1 == 0)
		goto cleanup;
	for (i = 0; i < bmap->n_beads; ++i) {
		int chr = bmap->beads[i].chr;
		int a;
		assert(chr >= 0 && chr < n_chr);
		++counts[chr];
		for (a = 0; a < 3; ++a) {
			centroid0[chr][a] += coords[hk_diploid_bid(i, HK_DIPLOID_COPY0)][a];
			centroid1[chr][a] += coords[hk_diploid_bid(i, HK_DIPLOID_COPY1)][a];
		}
	}
	for (c = 0; c < n_chr; ++c) {
		fvec3_t f0, f1;
		float e;
		int a;
		if (counts[c] <= 0)
			continue;
		for (a = 0; a < 3; ++a) {
			centroid0[c][a] /= (float)counts[c];
			centroid1[c][a] /= (float)counts[c];
		}
		e = hk_blind_homolog_sep_energy_force(centroid0[c], centroid1[c], unit,
											  min_sep_unit, lambda_sep, f0, f1);
		total += e;
		if (e <= 0.0f)
			continue;
		for (i = 0; i < bmap->n_beads; ++i) {
			int a2;
			if (bmap->beads[i].chr != c)
				continue;
			for (a2 = 0; a2 < 3; ++a2) {
				float g0 = f0[a2] / (float)counts[c];
				float g1 = f1[a2] / (float)counts[c];
				int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
				int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
				force[b0][a2] += g0;
				force[b1][a2] += g1;
				if (isfinite(g0)) l1 += fabs((double)g0);
				else ++n_bad;
				if (isfinite(g1)) l1 += fabs((double)g1);
				else ++n_bad;
			}
		}
	}
cleanup:
	free(counts);
	free(centroid0);
	free(centroid1);
	if (force_l1) *force_l1 = (float)l1;
	if (n_nonfinite) *n_nonfinite = n_bad;
	return (float)total;
}

void hk_blind_homolog_sep_stats_init(struct hk_blind_homolog_sep_stats *stats)
{
	assert(stats);
	stats->n_haploid = 0;
	stats->n_finite = 0;
	stats->n_nonfinite = 0;
	stats->n_collapsed = 0;
	stats->min_sep = 0.0f;
	stats->max_sep = 0.0f;
	stats->mean_sep = 0.0f;
	stats->collapse_threshold = 0.0f;
}

static int hk_blind_homolog_pair_finite(const fvec3_t x0, const fvec3_t x1)
{
	int a;
	for (a = 0; a < 3; ++a)
		if (!hk_blind_float_isfinite(x0[a]) || !hk_blind_float_isfinite(x1[a]))
			return 0;
	return 1;
}

int hk_blind_homolog_sep_compute_stats(int32_t n_haploid, const fvec3_t *coords, float collapse_threshold,
									   struct hk_blind_homolog_sep_stats *stats)
{
	double sum = 0.0;
	int32_t i;

	assert(n_haploid >= 0);
	assert(stats);
	assert(isfinite(collapse_threshold));
	assert(collapse_threshold >= 0.0f);
	if (n_haploid > 0)
		assert(coords);

	hk_blind_homolog_sep_stats_init(stats);
	stats->n_haploid = n_haploid;
	stats->collapse_threshold = collapse_threshold;
	for (i = 0; i < n_haploid; ++i) {
		int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
		int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
		float sep;

		if (!hk_blind_homolog_pair_finite(coords[b0], coords[b1])) {
			++stats->n_nonfinite;
			continue;
		}
		sep = hk_blind_coord_dist(coords[b0], coords[b1]);
		if (stats->n_finite == 0) {
			stats->min_sep = sep;
			stats->max_sep = sep;
		} else {
			if (sep < stats->min_sep) stats->min_sep = sep;
			if (sep > stats->max_sep) stats->max_sep = sep;
		}
		sum += sep;
		++stats->n_finite;
		if (sep < collapse_threshold)
			++stats->n_collapsed;
	}
	if (stats->n_finite > 0)
		stats->mean_sep = (float)(sum / stats->n_finite);
	return 0;
}

void hk_blind_gauge_stats_init(struct hk_blind_gauge_stats *stats)
{
	assert(stats);
	stats->n_chr = 0;
	stats->n_flipped = 0;
	stats->mismatch_no_flip = 0.0;
	stats->mismatch_flip = 0.0;
	stats->mismatch_chosen = 0.0;
}

static void hk_blind_swap_coord(fvec3_t x, fvec3_t y)
{
	float t;
	int a;
	for (a = 0; a < 3; ++a) {
		t = x[a];
		x[a] = y[a];
		y[a] = t;
	}
}

int hk_blind_temporal_gauge_stabilize_coords(int32_t n_haploid, const int32_t *chr_by_haploid, int32_t n_chr,
											 const fvec3_t *prev_coords, fvec3_t *cur_coords,
											 uint8_t *chr_flipped, struct hk_blind_gauge_stats *stats)
{
	double mismatch_no_total = 0.0, mismatch_flip_total = 0.0, mismatch_chosen_total = 0.0;
	int32_t c, i;
	int32_t n_flipped = 0;

	assert(n_haploid >= 0);
	assert(n_chr >= 0);
	if (n_haploid > 0) {
		assert(chr_by_haploid);
		assert(prev_coords);
		assert(cur_coords);
	}
	if (stats)
		hk_blind_gauge_stats_init(stats);
	for (c = 0; c < n_chr; ++c)
		if (chr_flipped) chr_flipped[c] = 0;
	for (i = 0; i < n_haploid; ++i)
		assert(chr_by_haploid[i] >= 0 && chr_by_haploid[i] < n_chr);

	for (c = 0; c < n_chr; ++c) {
		double mismatch_no = 0.0, mismatch_flip = 0.0;
		for (i = 0; i < n_haploid; ++i) {
			int32_t chr = chr_by_haploid[i];
			int32_t b0, b1;
			int a;

			if (chr != c) continue;
			b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
			b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
			for (a = 0; a < 3; ++a) {
				assert(hk_blind_float_isfinite(prev_coords[b0][a]));
				assert(hk_blind_float_isfinite(prev_coords[b1][a]));
				assert(hk_blind_float_isfinite(cur_coords[b0][a]));
				assert(hk_blind_float_isfinite(cur_coords[b1][a]));
			}
			mismatch_no += hk_blind_coord_dist2(cur_coords[b0], prev_coords[b0]);
			mismatch_no += hk_blind_coord_dist2(cur_coords[b1], prev_coords[b1]);
			mismatch_flip += hk_blind_coord_dist2(cur_coords[b1], prev_coords[b0]);
			mismatch_flip += hk_blind_coord_dist2(cur_coords[b0], prev_coords[b1]);
		}
		mismatch_no_total += mismatch_no;
		mismatch_flip_total += mismatch_flip;
		if (mismatch_flip < mismatch_no) {
			for (i = 0; i < n_haploid; ++i) {
				if (chr_by_haploid[i] == c) {
					int32_t b0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
					int32_t b1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
					hk_blind_swap_coord(cur_coords[b0], cur_coords[b1]);
				}
			}
			if (chr_flipped) chr_flipped[c] = 1;
			++n_flipped;
			mismatch_chosen_total += mismatch_flip;
		} else {
			mismatch_chosen_total += mismatch_no;
		}
	}
	if (stats) {
		stats->n_chr = n_chr;
		stats->n_flipped = n_flipped;
		stats->mismatch_no_flip = mismatch_no_total;
		stats->mismatch_flip = mismatch_flip_total;
		stats->mismatch_chosen = mismatch_chosen_total;
	}
	return 0;
}

static int32_t hk_blind_bmap_n_chr(const struct hk_bmap *bmap)
{
	int32_t n_chr = 0, i;

	assert(bmap);
	assert(bmap->n_beads >= 0);
	// Prefer the sequence dictionary count when present; otherwise use max bead chr + 1.
	if (bmap->d) {
		assert(bmap->d->n >= 0);
		n_chr = bmap->d->n;
	}
	if (bmap->n_beads > 0)
		assert(bmap->beads);
	for (i = 0; i < bmap->n_beads; ++i) {
		int32_t chr = bmap->beads[i].chr;
		assert(chr >= 0);
		if (bmap->d)
			assert(chr < bmap->d->n);
		if (chr + 1 > n_chr)
			n_chr = chr + 1;
	}
	return n_chr;
}

int hk_blind_temporal_gauge_stabilize_bmap(const struct hk_bmap *bmap, const fvec3_t *prev_coords, fvec3_t *cur_coords,
										   uint8_t *chr_flipped, struct hk_blind_gauge_stats *stats)
{
	int32_t *chr_by_haploid = 0;
	int32_t n_chr, i;
	int ret;

	assert(bmap);
	assert(bmap->n_beads >= 0);
	n_chr = hk_blind_bmap_n_chr(bmap);
	if (bmap->n_beads > 0) {
		assert(bmap->beads);
		chr_by_haploid = MALLOC(int32_t, bmap->n_beads);
		for (i = 0; i < bmap->n_beads; ++i)
			chr_by_haploid[i] = bmap->beads[i].chr;
	}
	ret = hk_blind_temporal_gauge_stabilize_coords(bmap->n_beads, chr_by_haploid, n_chr, prev_coords, cur_coords, chr_flipped, stats);
	free(chr_by_haploid);
	return ret;
}

void hk_blind_p4_apply_endpoint_flips(const float in_p4[HK_BLIND_N_STATE], int flip_endpoint0, int flip_endpoint1,
									  float out_p4[HK_BLIND_N_STATE])
{
	float in[HK_BLIND_N_STATE];
	int s;

	assert(in_p4);
	assert(out_p4);
	for (s = 0; s < HK_BLIND_N_STATE; ++s) {
		assert(isfinite(in_p4[s]));
		in[s] = in_p4[s];
	}

	if (flip_endpoint0) {
		if (flip_endpoint1) {
			out_p4[HK_BLIND_STATE_00] = in[HK_BLIND_STATE_11];
			out_p4[HK_BLIND_STATE_01] = in[HK_BLIND_STATE_10];
			out_p4[HK_BLIND_STATE_10] = in[HK_BLIND_STATE_01];
			out_p4[HK_BLIND_STATE_11] = in[HK_BLIND_STATE_00];
		} else {
			out_p4[HK_BLIND_STATE_00] = in[HK_BLIND_STATE_10];
			out_p4[HK_BLIND_STATE_01] = in[HK_BLIND_STATE_11];
			out_p4[HK_BLIND_STATE_10] = in[HK_BLIND_STATE_00];
			out_p4[HK_BLIND_STATE_11] = in[HK_BLIND_STATE_01];
		}
	} else {
		if (flip_endpoint1) {
			out_p4[HK_BLIND_STATE_00] = in[HK_BLIND_STATE_01];
			out_p4[HK_BLIND_STATE_01] = in[HK_BLIND_STATE_00];
			out_p4[HK_BLIND_STATE_10] = in[HK_BLIND_STATE_11];
			out_p4[HK_BLIND_STATE_11] = in[HK_BLIND_STATE_10];
		} else {
			for (s = 0; s < HK_BLIND_N_STATE; ++s)
				out_p4[s] = in[s];
		}
	}
}

int hk_blind_bpair_set_apply_chr_flips(struct hk_blind_bpair_set *set, const struct hk_bmap *bmap,
									   const uint8_t *chr_flipped, int32_t n_chr)
{
	int32_t i;

	assert(set);
	assert(bmap);
	assert(n_chr >= 0);
	assert(set->n_bpairs >= 0);
	assert(set->n_bpairs == 0 || set->bpairs);
	assert(bmap->n_beads >= 0);
	if (bmap->n_beads > 0)
		assert(bmap->beads);
	if (n_chr > 0)
		assert(chr_flipped);

	for (i = 0; i < set->n_bpairs; ++i) {
		struct hk_blind_bpair *bp = &set->bpairs[i];
		int32_t bid0 = bp->key.bid[0], bid1 = bp->key.bid[1];
		int32_t chr0, chr1;
		int flip0, flip1;

		assert(bid0 >= 0 && bid0 < bmap->n_beads);
		assert(bid1 >= 0 && bid1 < bmap->n_beads);
		chr0 = bmap->beads[bid0].chr;
		chr1 = bmap->beads[bid1].chr;
		assert(chr0 >= 0 && chr0 < n_chr);
		assert(chr1 >= 0 && chr1 < n_chr);
		flip0 = chr_flipped[chr0] != 0;
		flip1 = chr_flipped[chr1] != 0;
		hk_blind_p4_apply_endpoint_flips(bp->p4, flip0, flip1, bp->p4);
		hk_blind_p4_uncertainty(bp->p4, &bp->entropy, &bp->pmax, &bp->margin, &bp->rho_output, &bp->pU);
	}
	return 0;
}

int hk_blind_raw_state_to_canonical_state(int raw_state, uint8_t swapped)
{
	int a, b;
	assert(raw_state >= 0 && raw_state < HK_BLIND_N_STATE);
	if (!swapped)
		return raw_state;
	a = (raw_state >> 1) & 1;
	b = raw_state & 1;
	return (b << 1) | a;
}

static int hk_blind_canonical_state_to_raw_state(int canonical_state, uint8_t swapped)
{
	return hk_blind_raw_state_to_canonical_state(canonical_state, swapped);
}

void hk_blind_p4_to_raw_order(const float canonical_p4[HK_BLIND_N_STATE], uint8_t swapped, float raw_p4[HK_BLIND_N_STATE])
{
	float p00, p01, p10, p11;
	assert(canonical_p4);
	assert(raw_p4);
	if (swapped) {
		p00 = canonical_p4[HK_BLIND_STATE_00];
		p01 = canonical_p4[HK_BLIND_STATE_01];
		p10 = canonical_p4[HK_BLIND_STATE_10];
		p11 = canonical_p4[HK_BLIND_STATE_11];
		raw_p4[HK_BLIND_STATE_00] = p00;
		raw_p4[HK_BLIND_STATE_01] = p10;
		raw_p4[HK_BLIND_STATE_10] = p01;
		raw_p4[HK_BLIND_STATE_11] = p11;
	} else {
		int i;
		for (i = 0; i < HK_BLIND_N_STATE; ++i)
			raw_p4[i] = canonical_p4[i];
	}
}

static void hk_blind_wedge_set(struct hk_blind_wedge *edge, int32_t bid0, int32_t bid1, float k, float d_scale, int8_t state)
{
	assert(edge);
	assert(bid0 >= 0);
	assert(bid1 >= 0);
	assert(k >= 0.0f);
	assert(d_scale > 0.0f);
	assert(state >= 0 && state < HK_BLIND_N_STATE);
	edge->bid[0] = bid0;
	edge->bid[1] = bid1;
	edge->k = k;
	edge->d_scale = d_scale;
	edge->state = state;
	edge->state_mask = HK_BLIND_STATE_MASK(state);
}

void hk_blind_bpair_expand_weighted_edges(const struct hk_blind_bpair *bp, float base_k, float base_d_scale, float rho_train,
										  struct hk_blind_wedge out_edges[HK_BLIND_N_STATE])
{
	hk_blind_bpair_expand_weighted_edges_mode(bp, base_k, base_d_scale, rho_train,
											 HK_BLIND_RHO_TRAIN_CONSTANT, HK_BLIND_D_SCALE_RAW_COUNT,
											 1e-6f, out_edges);
}

float hk_blind_bpair_effective_rho_train(const struct hk_blind_bpair *bp, float rho_train, int rho_train_mode)
{
	return hk_blind_bpair_effective_rho_train_floor(bp, rho_train, rho_train_mode,
												   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR);
}

static float hk_blind_rho_train_floor_or_default(float rho_train_floor)
{
	if (rho_train_floor == 0.0f)
		return HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	assert(isfinite(rho_train_floor));
	assert(rho_train_floor >= 0.0f);
	assert(rho_train_floor <= 1.0f);
	return rho_train_floor;
}

float hk_blind_bpair_effective_rho_train_floor(const struct hk_blind_bpair *bp, float rho_train,
											   int rho_train_mode, float rho_train_floor)
{
	assert(bp);
	assert(isfinite(rho_train));
	assert(rho_train >= 0.0f);
	assert(hk_blind_rho_train_mode_valid(rho_train_mode));
	if (rho_train_mode == HK_BLIND_RHO_TRAIN_ENTROPY ||
		rho_train_mode == HK_BLIND_RHO_TRAIN_ENTROPY_WITH_FLOOR ||
		rho_train_mode == HK_BLIND_RHO_TRAIN_ENTROPY_CIS_CONSTANT_TRANS ||
		rho_train_mode == HK_BLIND_RHO_TRAIN_ENTROPY_CIS_FLOOR_TRANS) {
		float rho_output;
		float floor_factor;
		assert(isfinite(bp->rho_output));
		rho_output = hk_blind_clip01(bp->rho_output);
		floor_factor = hk_blind_rho_train_floor_or_default(rho_train_floor);
		if (rho_train_mode == HK_BLIND_RHO_TRAIN_ENTROPY)
			return rho_train * rho_output;
		if (rho_train_mode == HK_BLIND_RHO_TRAIN_ENTROPY_WITH_FLOOR)
			return rho_train * (rho_output > floor_factor? rho_output : floor_factor);
		if (hk_blind_bpair_contact_class_or_default(bp) == HK_BLIND_CONTACT_TRANS) {
			if (rho_train_mode == HK_BLIND_RHO_TRAIN_ENTROPY_CIS_CONSTANT_TRANS)
				return rho_train;
			return rho_train * (rho_output > floor_factor? rho_output : floor_factor);
		}
		return rho_train * rho_output;
	}
	return rho_train;
}

void hk_blind_bpair_expand_weighted_edges_mode(const struct hk_blind_bpair *bp, float base_k, float base_d_scale,
											   float rho_train, int rho_train_mode, int d_scale_mode,
											   float d_scale_eps_count,
											   struct hk_blind_wedge out_edges[HK_BLIND_N_STATE])
{
	hk_blind_bpair_expand_weighted_edges_mode_ex(bp, base_k, base_d_scale, rho_train,
												 rho_train_mode, HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
												 d_scale_mode, d_scale_eps_count, out_edges);
}

void hk_blind_bpair_expand_weighted_edges_mode_ex(const struct hk_blind_bpair *bp, float base_k, float base_d_scale,
												  float rho_train, int rho_train_mode, float rho_train_floor,
												  int d_scale_mode, float d_scale_eps_count,
												  struct hk_blind_wedge out_edges[HK_BLIND_N_STATE])
{
	int32_t i, j, i0, i1, j0, j1;
	float rho_eff;
	int s;

	assert(bp);
	assert(bp->key.bid[0] >= 0);
	assert(bp->key.bid[1] >= 0);
	assert(out_edges);
	assert(base_k >= 0.0f);
	assert(base_d_scale > 0.0f);
	assert(rho_train >= 0.0f);
	assert(hk_blind_rho_train_mode_valid(rho_train_mode));
	assert(hk_blind_d_scale_mode_valid(d_scale_mode));
	assert(isfinite(d_scale_eps_count));
	assert(d_scale_eps_count > 0.0f);
	assert(bp->n_raw > 0);
	for (s = 0; s < HK_BLIND_N_STATE; ++s) {
		assert(isfinite(bp->p4[s]));
		assert(bp->p4[s] >= 0.0f);
	}

	i = bp->key.bid[0];
	j = bp->key.bid[1];
	i0 = hk_diploid_bid(i, HK_DIPLOID_COPY0);
	i1 = hk_diploid_bid(i, HK_DIPLOID_COPY1);
	j0 = hk_diploid_bid(j, HK_DIPLOID_COPY0);
	j1 = hk_diploid_bid(j, HK_DIPLOID_COPY1);

	rho_eff = hk_blind_bpair_effective_rho_train_floor(bp, rho_train, rho_train_mode, rho_train_floor);
	for (s = 0; s < HK_BLIND_N_STATE; ++s) {
		float d_scale = base_d_scale;
		float k = base_k * rho_eff * bp->p4[s];
		if (d_scale_mode == HK_BLIND_D_SCALE_EXPECTED_COUNT) {
			float n_eff = (float)bp->n_raw * bp->p4[s];
			if (n_eff < d_scale_eps_count) n_eff = d_scale_eps_count;
			d_scale = powf(n_eff, -1.0f / 3.0f);
		} else if (d_scale_mode == HK_BLIND_D_SCALE_DENSITY_NORMALIZED_RAW_COUNT) {
			float n_eff = hk_blind_bpair_density_normalized_raw_count(bp, d_scale_eps_count);
			d_scale = powf(n_eff, -1.0f / 3.0f);
		} else if (d_scale_mode == HK_BLIND_D_SCALE_CAPPED_DENSITY_RAW_COUNT) {
			float n_eff = hk_blind_bpair_capped_density_raw_count(bp, d_scale_eps_count);
			d_scale = powf(n_eff, -1.0f / 3.0f);
		}
		switch (s) {
		case HK_BLIND_STATE_00:
			hk_blind_wedge_set(&out_edges[s], i0, j0, k, d_scale, s);
			break;
		case HK_BLIND_STATE_01:
			hk_blind_wedge_set(&out_edges[s], i0, j1, k, d_scale, s);
			break;
		case HK_BLIND_STATE_10:
			hk_blind_wedge_set(&out_edges[s], i1, j0, k, d_scale, s);
			break;
		case HK_BLIND_STATE_11:
			hk_blind_wedge_set(&out_edges[s], i1, j1, k, d_scale, s);
			break;
		}
	}
}

float hk_blind_wedge_contact_energy(const struct hk_fdg_conf *conf, const struct hk_blind_wedge *edge, float distance, float unit)
{
	assert(conf);
	assert(edge);
	assert(edge->k >= 0.0f);
	assert(edge->d_scale > 0.0f);
	assert(distance >= 0.0f);
	assert(unit > 0.0f);
	return hk_fdg_contact_energy_dist(conf, distance, unit, edge->d_scale, edge->k);
}

static float hk_blind_contact_force_r(const struct hk_fdg_conf *conf, float r, float k)
{
	float t;

	assert(conf);
	assert(r >= 0.0f);
	assert(k >= 0.0f);
	if (r < conf->d_c1) {
		t = conf->d_c1 - r;
		return 2.0f * k * t;
	} else if (r <= conf->d_c2) {
		return 0.0f;
	} else if (r <= conf->d_c3) {
		t = r - conf->d_c2;
		return -2.0f * k * t;
	} else {
		t = r - conf->d_c2;
		return -k * (conf->c_c1 - conf->c_c2 / (t * t));
	}
}

void hk_blind_contact_class_diag_init(struct hk_blind_contact_class_diag *diag)
{
	assert(diag);
	diag->n_wedges_cis = 0;
	diag->n_wedges_trans = 0;
	diag->n_nonfinite = 0;
	diag->sum_wedge_k_cis = 0.0;
	diag->sum_wedge_k_trans = 0.0;
	diag->contact_energy_cis = 0.0f;
	diag->contact_energy_trans = 0.0f;
	diag->contact_force_l1_cis = 0.0f;
	diag->contact_force_l1_trans = 0.0f;
}

static int hk_blind_wedge_contact_class_from_bmap(const struct hk_bmap *bmap, const struct hk_blind_wedge *edge)
{
	int32_t h0, h1;
	assert(bmap);
	assert(edge);
	h0 = hk_diploid_haploid_bid(edge->bid[0]);
	h1 = hk_diploid_haploid_bid(edge->bid[1]);
	assert(h0 >= 0 && h0 < bmap->n_beads);
	assert(h1 >= 0 && h1 < bmap->n_beads);
	return bmap->beads[h0].chr == bmap->beads[h1].chr? HK_BLIND_CONTACT_CIS : HK_BLIND_CONTACT_TRANS;
}

int hk_blind_contact_class_diag_accumulate(const struct hk_fdg_conf *conf, const struct hk_bmap *bmap,
										   const struct hk_blind_wedge_list *edges, const fvec3_t *coords,
										   int32_t n_haploid, float unit, struct hk_blind_contact_class_diag *diag)
{
	int32_t n_diploid, eidx;

	assert(conf);
	assert(bmap);
	assert(edges);
	assert(coords || n_haploid == 0);
	assert(diag);
	assert(n_haploid >= 0);
	assert(bmap->n_beads == n_haploid);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(n_haploid <= INT32_MAX / HK_DIPLOID_N_COPY);
	n_diploid = n_haploid * HK_DIPLOID_N_COPY;
	hk_blind_contact_class_diag_init(diag);
	for (eidx = 0; eidx < edges->n_edges; ++eidx) {
		const struct hk_blind_wedge *edge = &edges->edges[eidx];
		int32_t i = edge->bid[0], j = edge->bid[1];
		int contact_class;
		float dx, dy, dz, d, inv_d, dir[3], r, energy, fmag, force_l1;
		assert(i >= 0 && i < n_diploid);
		assert(j >= 0 && j < n_diploid);
		assert(i != j);
		assert(edge->k >= 0.0f);
		assert(edge->d_scale > 0.0f);
		contact_class = hk_blind_wedge_contact_class_from_bmap(bmap, edge);
		if (contact_class == HK_BLIND_CONTACT_TRANS) {
			++diag->n_wedges_trans;
			diag->sum_wedge_k_trans += edge->k;
		} else {
			++diag->n_wedges_cis;
			diag->sum_wedge_k_cis += edge->k;
		}
		dx = coords[i][0] - coords[j][0];
		dy = coords[i][1] - coords[j][1];
		dz = coords[i][2] - coords[j][2];
		d = sqrtf(dx * dx + dy * dy + dz * dz);
		if (d > 0.0f) {
			inv_d = 1.0f / d;
			dir[0] = dx * inv_d;
			dir[1] = dy * inv_d;
			dir[2] = dz * inv_d;
		} else {
			dir[0] = 1.0f;
			dir[1] = 0.0f;
			dir[2] = 0.0f;
		}
		r = (d / unit) / edge->d_scale;
		energy = hk_fdg_contact_energy_dist(conf, d, unit, edge->d_scale, edge->k);
		fmag = hk_blind_contact_force_r(conf, r, edge->k);
		force_l1 = 2.0f * fabsf(fmag) * (fabsf(dir[0]) + fabsf(dir[1]) + fabsf(dir[2]));
		if (!hk_blind_float_isfinite(energy) || !hk_blind_float_isfinite(force_l1)) {
			++diag->n_nonfinite;
			continue;
		}
		if (contact_class == HK_BLIND_CONTACT_TRANS) {
			diag->contact_energy_trans += energy;
			diag->contact_force_l1_trans += force_l1;
		} else {
			diag->contact_energy_cis += energy;
			diag->contact_force_l1_cis += force_l1;
		}
	}
	return 0;
}

float hk_blind_wedge_list_accumulate_contact_force_cpu(const struct hk_fdg_conf *conf, const struct hk_blind_wedge_list *edges,
													   const fvec3_t *coords, int32_t n_diploid, float unit,
													   fvec3_t *force, int32_t *n_nonfinite)
{
	double total_energy = 0.0;
	int32_t eidx;

	assert(conf);
	assert(edges);
	assert(edges->n_edges >= 0);
	assert(edges->n_edges == 0 || edges->edges);
	assert(n_diploid >= 0);
	assert(unit > 0.0f);
	if (n_diploid > 0) {
		assert(coords);
		assert(force);
	}
	if (n_nonfinite)
		*n_nonfinite = 0;

	for (eidx = 0; eidx < edges->n_edges; ++eidx) {
		const struct hk_blind_wedge *edge = &edges->edges[eidx];
		int32_t i = edge->bid[0], j = edge->bid[1];
		float dx, dy, dz, d, inv_d, dir[3], r, energy, fmag, f[3];
		int bad = 0;

		assert(i >= 0 && i < n_diploid);
		assert(j >= 0 && j < n_diploid);
		assert(i != j);
		assert(edge->k >= 0.0f);
		assert(edge->d_scale > 0.0f);

		dx = coords[i][0] - coords[j][0];
		dy = coords[i][1] - coords[j][1];
		dz = coords[i][2] - coords[j][2];
		d = sqrtf(dx * dx + dy * dy + dz * dz);
		if (d > 0.0f) {
			inv_d = 1.0f / d;
			dir[0] = dx * inv_d;
			dir[1] = dy * inv_d;
			dir[2] = dz * inv_d;
		} else {
			// Existing FDG asserts nonzero contact distance; blind reference uses a deterministic safe fallback.
			dir[0] = 1.0f;
			dir[1] = 0.0f;
			dir[2] = 0.0f;
		}
		r = (d / unit) / edge->d_scale;
		energy = hk_fdg_contact_energy_dist(conf, d, unit, edge->d_scale, edge->k);
		fmag = hk_blind_contact_force_r(conf, r, edge->k);
		f[0] = fmag * dir[0];
		f[1] = fmag * dir[1];
		f[2] = fmag * dir[2];

		if (!hk_blind_float_isfinite(energy) ||
			!hk_blind_float_isfinite(f[0]) || !hk_blind_float_isfinite(f[1]) || !hk_blind_float_isfinite(f[2]))
			bad = 1;
		if (bad) {
			if (n_nonfinite)
				++(*n_nonfinite);
			continue;
		}
		total_energy += energy;
		force[i][0] += f[0];
		force[i][1] += f[1];
		force[i][2] += f[2];
		force[j][0] -= f[0];
		force[j][1] -= f[1];
		force[j][2] -= f[2];
	}
	return (float)total_energy;
}

void hk_blind_backbone_diag_init(struct hk_blind_backbone_diag *diag)
{
	assert(diag);
	diag->n_edges = 0;
	diag->n_nonfinite = 0;
	diag->energy = 0.0f;
	diag->force_l1 = 0.0f;
}

static float hk_blind_backbone_energy_r(const struct hk_fdg_conf *conf, float r)
{
	float t;

	assert(conf);
	assert(r >= 0.0f);
	if (r < conf->d_b1) {
		t = conf->d_b1 - r;
		return t * t;
	} else if (r <= conf->d_b2) {
		return 0.0f;
	} else {
		t = r - conf->d_b2;
		return t * t;
	}
}

static float hk_blind_backbone_force_r(const struct hk_fdg_conf *conf, float r)
{
	float t;

	assert(conf);
	assert(r >= 0.0f);
	if (r < conf->d_b1) {
		t = conf->d_b1 - r;
		return 2.0f * t;
	} else if (r <= conf->d_b2) {
		return 0.0f;
	} else {
		t = r - conf->d_b2;
		return -2.0f * t;
	}
}

float hk_blind_backbone_energy_force_cpu(const struct hk_fdg_conf *conf, const fvec3_t x0, const fvec3_t x1,
										 float unit, float d_scale, fvec3_t f0, fvec3_t f1)
{
	float dx, dy, dz, d, inv_d, r, energy, fmag, dir[3];
	int a;

	assert(conf);
	assert(x0);
	assert(x1);
	assert(f0);
	assert(f1);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(isfinite(d_scale));
	assert(d_scale > 0.0f);
	assert(conf->d_b1 >= 0.0f);
	assert(conf->d_b2 >= conf->d_b1);
	for (a = 0; a < 3; ++a) {
		assert(isfinite(x0[a]));
		assert(isfinite(x1[a]));
	}

	dx = x0[0] - x1[0];
	dy = x0[1] - x1[1];
	dz = x0[2] - x1[2];
	d = sqrtf(dx * dx + dy * dy + dz * dz);
	if (d > 0.0f) {
		inv_d = 1.0f / d;
		dir[0] = dx * inv_d;
		dir[1] = dy * inv_d;
		dir[2] = dz * inv_d;
	} else {
		// Existing FDG asserts nonzero backbone distance; blind reference uses a deterministic safe fallback.
		dir[0] = 1.0f;
		dir[1] = 0.0f;
		dir[2] = 0.0f;
	}
	r = (d / unit) / d_scale;
	energy = hk_blind_backbone_energy_r(conf, r);
	fmag = hk_blind_backbone_force_r(conf, r);
	for (a = 0; a < 3; ++a) {
		f0[a] = fmag * dir[a];
		f1[a] = -f0[a];
		assert(isfinite(f0[a]));
		assert(isfinite(f1[a]));
	}
	assert(isfinite(energy));
	return energy;
}

static int32_t hk_blind_bmap_mid_bead_size(const struct hk_bmap *bmap)
{
	int32_t *sizes;
	int32_t mid;
	int32_t i;

	assert(bmap);
	assert(bmap->n_beads > 0);
	assert(bmap->beads);
	sizes = MALLOC(int32_t, bmap->n_beads);
	for (i = 0; i < bmap->n_beads; ++i) {
		sizes[i] = bmap->beads[i].en - bmap->beads[i].st;
		assert(sizes[i] > 0);
	}
	mid = ks_ksmall_int32_t(bmap->n_beads, sizes, (int)(bmap->n_beads * 0.5));
	free(sizes);
	assert(mid > 0);
	return mid;
}

static float hk_blind_backbone_d_scale(const struct hk_bead *a, const struct hk_bead *b, int32_t mid_dist)
{
	int32_t d;

	assert(a);
	assert(b);
	assert(mid_dist > 0);
	d = ((a->en - a->st) + (b->en - b->st)) / 2;
	assert(d > 0);
	return pow((double)d / mid_dist, 1.0 / 3.0);
}

float hk_blind_backbone_accumulate_force_cpu(const struct hk_fdg_conf *conf, const struct hk_bmap *bmap,
											 const fvec3_t *coords, float unit, fvec3_t *force,
											 struct hk_blind_backbone_diag *diag)
{
	struct hk_blind_backbone_diag local_diag;
	double total_energy = 0.0;
	double force_l1 = 0.0;
	int32_t mid_dist;
	int32_t i;

	assert(conf);
	assert(bmap);
	assert(bmap->n_beads >= 0);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	if (diag == 0)
		diag = &local_diag;
	hk_blind_backbone_diag_init(diag);
	if (bmap->n_beads == 0)
		return 0.0f;
	assert(bmap->beads);
	assert(coords);
	assert(force);

	mid_dist = hk_blind_bmap_mid_bead_size(bmap);
	for (i = 1; i < bmap->n_beads; ++i) {
		const struct hk_bead *prev = &bmap->beads[i - 1];
		const struct hk_bead *cur = &bmap->beads[i];
		float d_scale;
		int copy;

		if (prev->chr != cur->chr)
			continue;
		d_scale = hk_blind_backbone_d_scale(prev, cur, mid_dist);
		for (copy = 0; copy < HK_DIPLOID_N_COPY; ++copy) {
			int32_t b0 = hk_diploid_bid(i - 1, copy);
			int32_t b1 = hk_diploid_bid(i, copy);
			fvec3_t f0, f1;
			float e;
			int bad = 0;
			int a;

			e = hk_blind_backbone_energy_force_cpu(conf, coords[b0], coords[b1], unit, d_scale, f0, f1);
			if (!hk_blind_float_isfinite(e))
				bad = 1;
			for (a = 0; a < 3; ++a)
				if (!hk_blind_float_isfinite(f0[a]) || !hk_blind_float_isfinite(f1[a]))
					bad = 1;
			if (bad) {
				++diag->n_nonfinite;
				continue;
			}
			total_energy += e;
			++diag->n_edges;
			for (a = 0; a < 3; ++a) {
				force[b0][a] += f0[a];
				force[b1][a] += f1[a];
				force_l1 += fabs((double)f0[a]) + fabs((double)f1[a]);
				if (!hk_blind_float_isfinite(force[b0][a]) || !hk_blind_float_isfinite(force[b1][a]))
					++diag->n_nonfinite;
			}
		}
	}
	diag->energy = (float)total_energy;
	diag->force_l1 = (float)force_l1;
	return diag->energy;
}

float hk_blind_rel_rep_schedule_at(int32_t iter, int32_t n_iter)
{
	float progress;

	assert(iter >= 0);
	if (n_iter <= 0)
		return 1.0f;
	assert(iter < n_iter);
	progress = (float)(iter + 1) / (float)n_iter;
	return 1.0f / (1.0f + expf(-10.0f * (progress - 0.333f)));
}

void hk_blind_repulsion_diag_init(struct hk_blind_repulsion_diag *diag)
{
	assert(diag);
	diag->n_pairs_considered = 0;
	diag->n_pairs_blocked = 0;
	diag->n_pairs_active = 0;
	diag->n_nonfinite = 0;
	diag->energy = 0.0f;
	diag->force_l1 = 0.0f;
}

static int hk_blind_repulsion_mode_valid(int mode)
{
	return mode == HK_BLIND_REPULSION_NONE ||
		   mode == HK_BLIND_REPULSION_N2 ||
		   mode == HK_BLIND_REPULSION_CELL;
}

float hk_blind_repulsion_energy_force_cpu(const struct hk_fdg_conf *conf, const fvec3_t x0, const fvec3_t x1,
										  float unit, float rel_rep_k, fvec3_t f0, fvec3_t f1)
{
	float dx, dy, dz, d2, d, inv_d, r, t, k, energy, fmag, dir[3];
	int a;

	assert(conf);
	assert(x0);
	assert(x1);
	assert(f0);
	assert(f1);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(isfinite(rel_rep_k));
	assert(rel_rep_k >= 0.0f);
	assert(conf->d_r >= 0.0f);
	assert(conf->k_rel_rep >= 0.0f);
	for (a = 0; a < 3; ++a) {
		assert(isfinite(x0[a]));
		assert(isfinite(x1[a]));
	}

	dx = x0[0] - x1[0];
	dy = x0[1] - x1[1];
	dz = x0[2] - x1[2];
	d2 = dx * dx + dy * dy + dz * dz;
	d = sqrtf(d2);
	if (d > 0.0f) {
		inv_d = 1.0f / d;
		dir[0] = dx * inv_d;
		dir[1] = dy * inv_d;
		dir[2] = dz * inv_d;
	} else {
		// Existing FDG asserts nonzero repulsion distance; blind reference uses a deterministic safe fallback.
		dir[0] = 1.0f;
		dir[1] = 0.0f;
		dir[2] = 0.0f;
	}
	r = d / unit;
	k = conf->k_rel_rep * rel_rep_k;
	if (r >= conf->d_r) {
		energy = 0.0f;
		fmag = 0.0f;
	} else {
		t = conf->d_r - r;
		energy = k * t * t;
		fmag = 2.0f * k * t;
	}
	for (a = 0; a < 3; ++a) {
		f0[a] = fmag * dir[a];
		f1[a] = -f0[a];
		assert(isfinite(f0[a]));
		assert(isfinite(f1[a]));
	}
	assert(isfinite(energy));
	return energy;
}

static uint64_t hk_blind_block_pair_key(int32_t i, int32_t j)
{
	uint32_t a, b;

	assert(i >= 0);
	assert(j >= 0);
	assert(i != j);
	if (i < j) {
		a = (uint32_t)i;
		b = (uint32_t)j;
	} else {
		a = (uint32_t)j;
		b = (uint32_t)i;
	}
	return ((uint64_t)a << 32) | (uint64_t)b;
}

static int hk_blind_u64_cmp(const void *a_, const void *b_)
{
	uint64_t a = *(const uint64_t*)a_;
	uint64_t b = *(const uint64_t*)b_;
	return a < b? -1 : a > b? 1 : 0;
}

static int64_t hk_blind_u64_unique(uint64_t *a, int64_t n)
{
	int64_t i, m;

	assert(n >= 0);
	if (n == 0)
		return 0;
	assert(a);
	qsort(a, (size_t)n, sizeof(*a), hk_blind_u64_cmp);
	for (i = m = 1; i < n; ++i)
		if (a[i] != a[m - 1])
			a[m++] = a[i];
	return m;
}

static int hk_blind_u64_contains(const uint64_t *a, int64_t n, uint64_t key)
{
	int64_t lo = 0, hi = n;

	assert(n >= 0);
	while (lo < hi) {
		int64_t mid = lo + (hi - lo) / 2;
		if (a[mid] < key) lo = mid + 1;
		else hi = mid;
	}
	return lo < n && a[lo] == key;
}

static void hk_blind_add_block_pair(uint64_t *blocked, int64_t *n_blocked, int64_t m_blocked,
									int32_t i, int32_t j, int32_t n_diploid)
{
	assert(blocked);
	assert(n_blocked);
	assert(*n_blocked >= 0 && *n_blocked < m_blocked);
	assert(i >= 0 && i < n_diploid);
	assert(j >= 0 && j < n_diploid);
	if (i == j)
		return;
	blocked[(*n_blocked)++] = hk_blind_block_pair_key(i, j);
}

static int64_t hk_blind_build_repulsion_blocked_pairs(const struct hk_blind_wedge_list *contact_edges_or_null,
													  const struct hk_bmap *bmap_or_null, int32_t n_diploid,
													  float contact_k_min,
													  uint64_t **blocked_out)
{
	uint64_t *blocked = 0;
	int64_t n_blocked = 0, m_blocked = 0;

	assert(blocked_out);
	assert(n_diploid >= 0);
	assert(isfinite(contact_k_min));
	assert(contact_k_min >= 0.0f);
	*blocked_out = 0;
	if (contact_edges_or_null) {
		assert(contact_edges_or_null->n_edges >= 0);
		assert(contact_edges_or_null->n_edges == 0 || contact_edges_or_null->edges);
		m_blocked += contact_edges_or_null->n_edges;
	}
	if (bmap_or_null) {
		assert(bmap_or_null->n_beads >= 0);
		assert(bmap_or_null->n_beads <= INT32_MAX / HK_DIPLOID_N_COPY);
		assert(bmap_or_null->n_beads * HK_DIPLOID_N_COPY == n_diploid);
		if (bmap_or_null->n_beads > 0)
			assert(bmap_or_null->beads);
		m_blocked += (int64_t)HK_DIPLOID_N_COPY * (bmap_or_null->n_beads > 0? bmap_or_null->n_beads - 1 : 0);
	}
	if (m_blocked == 0)
		return 0;
	blocked = MALLOC(uint64_t, m_blocked);
	if (contact_edges_or_null) {
		int32_t e;
		for (e = 0; e < contact_edges_or_null->n_edges; ++e) {
				const struct hk_blind_wedge *edge = &contact_edges_or_null->edges[e];
				assert(edge->bid[0] >= 0 && edge->bid[0] < n_diploid);
				assert(edge->bid[1] >= 0 && edge->bid[1] < n_diploid);
				assert(isfinite(edge->k));
				if (edge->k < contact_k_min)
					continue;
				hk_blind_add_block_pair(blocked, &n_blocked, m_blocked, edge->bid[0], edge->bid[1], n_diploid);
			}
		}
	if (bmap_or_null) {
		int32_t i;
		for (i = 1; i < bmap_or_null->n_beads; ++i) {
			const struct hk_bead *prev = &bmap_or_null->beads[i - 1];
			const struct hk_bead *cur = &bmap_or_null->beads[i];
			int copy;
			if (prev->chr != cur->chr)
				continue;
			for (copy = 0; copy < HK_DIPLOID_N_COPY; ++copy)
				hk_blind_add_block_pair(blocked, &n_blocked, m_blocked,
										hk_diploid_bid(i - 1, copy), hk_diploid_bid(i, copy), n_diploid);
		}
	}
	n_blocked = hk_blind_u64_unique(blocked, n_blocked);
	*blocked_out = blocked;
	return n_blocked;
}

static void hk_blind_repulsion_accumulate_candidate(const struct hk_fdg_conf *conf, int32_t i, int32_t j,
													const fvec3_t *coords, float unit, float rel_rep_k,
													float rep_radius2, const uint64_t *blocked, int64_t n_blocked,
													fvec3_t *force, struct hk_blind_repulsion_diag *diag,
													double *total_energy, double *force_l1)
{
	uint64_t key;
	float dx, dy, dz, d2;

	assert(conf);
	assert(i >= 0);
	assert(j >= 0);
	assert(i != j);
	assert(coords);
	assert(force);
	assert(diag);
	assert(total_energy);
	assert(force_l1);
	key = hk_blind_block_pair_key(i, j);
	++diag->n_pairs_considered;
	if (hk_blind_u64_contains(blocked, n_blocked, key)) {
		++diag->n_pairs_blocked;
		return;
	}
	dx = coords[i][0] - coords[j][0];
	dy = coords[i][1] - coords[j][1];
	dz = coords[i][2] - coords[j][2];
	d2 = dx * dx + dy * dy + dz * dz;
	if (d2 < rep_radius2) {
		fvec3_t f0, f1;
		float e;
		int a, bad = 0;
		e = hk_blind_repulsion_energy_force_cpu(conf, coords[i], coords[j], unit, rel_rep_k, f0, f1);
		if (!hk_blind_float_isfinite(e))
			bad = 1;
		for (a = 0; a < 3; ++a)
			if (!hk_blind_float_isfinite(f0[a]) || !hk_blind_float_isfinite(f1[a]))
				bad = 1;
		if (bad) {
			++diag->n_nonfinite;
			return;
		}
		*total_energy += e;
		if (e > 0.0f)
			++diag->n_pairs_active;
		for (a = 0; a < 3; ++a) {
			force[i][a] += f0[a];
			force[j][a] += f1[a];
			*force_l1 += fabs((double)f0[a]) + fabs((double)f1[a]);
			if (!hk_blind_float_isfinite(force[i][a]) || !hk_blind_float_isfinite(force[j][a]))
				++diag->n_nonfinite;
		}
	}
}

static float hk_blind_repulsion_accumulate_force_cpu_impl(const struct hk_fdg_conf *conf, int32_t n_diploid,
														  const fvec3_t *coords, float unit, float rel_rep_k,
														  float contact_block_k_min,
														  const struct hk_blind_wedge_list *contact_edges_or_null,
														  const struct hk_bmap *bmap_or_null, fvec3_t *force,
														  struct hk_blind_repulsion_diag *diag)
{
	struct hk_blind_repulsion_diag local_diag;
	uint64_t *blocked = 0;
	int64_t n_blocked = 0;
	double total_energy = 0.0, force_l1 = 0.0;
	float rep_radius2;
	int32_t i, j;

	assert(conf);
	assert(n_diploid >= 0);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(isfinite(rel_rep_k));
	assert(rel_rep_k >= 0.0f);
	assert(isfinite(contact_block_k_min));
	assert(contact_block_k_min >= 0.0f);
	assert(conf->d_r >= 0.0f);
	if (n_diploid > 0) {
		assert(coords);
		assert(force);
	}
	if (diag == 0)
		diag = &local_diag;
	hk_blind_repulsion_diag_init(diag);
	if (n_diploid <= 1)
		return 0.0f;

	n_blocked = hk_blind_build_repulsion_blocked_pairs(contact_edges_or_null, bmap_or_null, n_diploid,
													   contact_block_k_min, &blocked);
	rep_radius2 = unit * conf->d_r;
	rep_radius2 *= rep_radius2;
	for (i = 0; i < n_diploid; ++i) {
		for (j = i + 1; j < n_diploid; ++j) {
			hk_blind_repulsion_accumulate_candidate(conf, i, j, coords, unit, rel_rep_k, rep_radius2,
													blocked, n_blocked, force, diag, &total_energy, &force_l1);
		}
	}
	free(blocked);
	diag->energy = (float)total_energy;
	diag->force_l1 = (float)force_l1;
	return diag->energy;
}

float hk_blind_repulsion_accumulate_force_cpu(const struct hk_fdg_conf *conf, int32_t n_diploid,
											  const fvec3_t *coords, float unit, float rel_rep_k,
											  const struct hk_blind_wedge_list *contact_edges_or_null,
											  const struct hk_bmap *bmap_or_null, fvec3_t *force,
											  struct hk_blind_repulsion_diag *diag)
{
	return hk_blind_repulsion_accumulate_force_cpu_impl(conf, n_diploid, coords, unit, rel_rep_k, 0.0f,
													   contact_edges_or_null, bmap_or_null, force, diag);
}

struct hk_blind_repulsion_cell_entry {
	int64_t c[3];
	int32_t bid;
};

struct hk_blind_repulsion_cell_group {
	int64_t c[3];
	int32_t start;
	int32_t end;
};

static int hk_blind_repulsion_cell_entry_cmp(const void *a_, const void *b_)
{
	const struct hk_blind_repulsion_cell_entry *a = (const struct hk_blind_repulsion_cell_entry*)a_;
	const struct hk_blind_repulsion_cell_entry *b = (const struct hk_blind_repulsion_cell_entry*)b_;
	int axis;
	for (axis = 0; axis < 3; ++axis) {
		if (a->c[axis] != b->c[axis])
			return a->c[axis] < b->c[axis]? -1 : 1;
	}
	if (a->bid != b->bid)
		return a->bid < b->bid? -1 : 1;
	return 0;
}

static int64_t hk_blind_repulsion_cell_coord(float x, float cell_size)
{
	assert(isfinite(x));
	assert(isfinite(cell_size));
	assert(cell_size > 0.0f);
	return (int64_t)floor((double)x / (double)cell_size);
}

static int hk_blind_repulsion_cell_group_find(const struct hk_blind_repulsion_cell_group *groups, int32_t n_groups,
											  int64_t cx, int64_t cy, int64_t cz)
{
	int32_t lo = 0, hi = n_groups;

	assert(n_groups >= 0);
	while (lo < hi) {
		int32_t mid = lo + (hi - lo) / 2;
		const struct hk_blind_repulsion_cell_group *g = &groups[mid];
		int cmp;
		if (g->c[0] != cx) cmp = g->c[0] < cx? -1 : 1;
		else if (g->c[1] != cy) cmp = g->c[1] < cy? -1 : 1;
		else if (g->c[2] != cz) cmp = g->c[2] < cz? -1 : 1;
		else cmp = 0;
		if (cmp < 0) lo = mid + 1;
		else hi = mid;
	}
	if (lo < n_groups &&
		groups[lo].c[0] == cx && groups[lo].c[1] == cy && groups[lo].c[2] == cz)
		return lo;
	return -1;
}

static float hk_blind_repulsion_accumulate_force_cell_cpu_impl(const struct hk_fdg_conf *conf, int32_t n_diploid,
															   const fvec3_t *coords, float unit, float rel_rep_k,
															   float contact_block_k_min,
															   const struct hk_blind_wedge_list *contact_edges_or_null,
															   const struct hk_bmap *bmap_or_null, fvec3_t *force,
															   struct hk_blind_repulsion_diag *diag)
{
	struct hk_blind_repulsion_diag local_diag;
	struct hk_blind_repulsion_cell_entry *entries = 0;
	struct hk_blind_repulsion_cell_group *groups = 0;
	uint64_t *blocked = 0;
	int64_t n_blocked = 0;
	double total_energy = 0.0, force_l1 = 0.0;
	float cutoff, rep_radius2;
	int32_t n_groups = 0;
	int32_t i, g;

	assert(conf);
	assert(n_diploid >= 0);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(isfinite(rel_rep_k));
	assert(rel_rep_k >= 0.0f);
	assert(isfinite(contact_block_k_min));
	assert(contact_block_k_min >= 0.0f);
	assert(conf->d_r >= 0.0f);
	if (n_diploid > 0) {
		assert(coords);
		assert(force);
	}
	if (diag == 0)
		diag = &local_diag;
	hk_blind_repulsion_diag_init(diag);
	if (n_diploid <= 1)
		return 0.0f;

	cutoff = unit * conf->d_r;
	if (!(cutoff > 0.0f))
		return 0.0f;
	rep_radius2 = cutoff * cutoff;
	n_blocked = hk_blind_build_repulsion_blocked_pairs(contact_edges_or_null, bmap_or_null, n_diploid,
													   contact_block_k_min, &blocked);
	entries = MALLOC(struct hk_blind_repulsion_cell_entry, n_diploid);
	groups = MALLOC(struct hk_blind_repulsion_cell_group, n_diploid);
	for (i = 0; i < n_diploid; ++i) {
		int a;
		for (a = 0; a < 3; ++a)
			assert(isfinite(coords[i][a]));
		entries[i].c[0] = hk_blind_repulsion_cell_coord(coords[i][0], cutoff);
		entries[i].c[1] = hk_blind_repulsion_cell_coord(coords[i][1], cutoff);
		entries[i].c[2] = hk_blind_repulsion_cell_coord(coords[i][2], cutoff);
		entries[i].bid = i;
	}
	qsort(entries, (size_t)n_diploid, sizeof(*entries), hk_blind_repulsion_cell_entry_cmp);
	for (i = 0; i < n_diploid; ) {
		int32_t j = i + 1;
		groups[n_groups].c[0] = entries[i].c[0];
		groups[n_groups].c[1] = entries[i].c[1];
		groups[n_groups].c[2] = entries[i].c[2];
		groups[n_groups].start = i;
		while (j < n_diploid &&
			   entries[j].c[0] == entries[i].c[0] &&
			   entries[j].c[1] == entries[i].c[1] &&
			   entries[j].c[2] == entries[i].c[2])
			++j;
		groups[n_groups].end = j;
		++n_groups;
		i = j;
	}
	for (g = 0; g < n_groups; ++g) {
		int dx, dy, dz;
		const struct hk_blind_repulsion_cell_group *cg = &groups[g];
		for (dx = -1; dx <= 1; ++dx) {
			for (dy = -1; dy <= 1; ++dy) {
				for (dz = -1; dz <= 1; ++dz) {
					int h = hk_blind_repulsion_cell_group_find(groups, n_groups,
															   cg->c[0] + dx, cg->c[1] + dy, cg->c[2] + dz);
					int32_t a, b;
					if (h < 0 || h < g)
						continue;
					if (h == g) {
						for (a = cg->start; a < cg->end; ++a)
							for (b = a + 1; b < cg->end; ++b)
								hk_blind_repulsion_accumulate_candidate(conf, entries[a].bid, entries[b].bid,
																		coords, unit, rel_rep_k, rep_radius2,
																		blocked, n_blocked, force, diag,
																		&total_energy, &force_l1);
					} else {
						const struct hk_blind_repulsion_cell_group *ng = &groups[h];
						for (a = cg->start; a < cg->end; ++a)
							for (b = ng->start; b < ng->end; ++b)
								hk_blind_repulsion_accumulate_candidate(conf, entries[a].bid, entries[b].bid,
																		coords, unit, rel_rep_k, rep_radius2,
																		blocked, n_blocked, force, diag,
																		&total_energy, &force_l1);
					}
				}
			}
		}
	}
	free(groups);
	free(entries);
	free(blocked);
	diag->energy = (float)total_energy;
	diag->force_l1 = (float)force_l1;
	return diag->energy;
}

float hk_blind_repulsion_accumulate_force_cell_cpu(const struct hk_fdg_conf *conf, int32_t n_diploid,
												   const fvec3_t *coords, float unit, float rel_rep_k,
												   const struct hk_blind_wedge_list *contact_edges_or_null,
												   const struct hk_bmap *bmap_or_null, fvec3_t *force,
												   struct hk_blind_repulsion_diag *diag)
{
	return hk_blind_repulsion_accumulate_force_cell_cpu_impl(conf, n_diploid, coords, unit, rel_rep_k, 0.0f,
															contact_edges_or_null, bmap_or_null, force, diag);
}

void hk_blind_step_diag_init(struct hk_blind_step_diag *diag)
{
	assert(diag);
	diag->contact_energy = 0.0f;
	diag->backbone_energy = 0.0f;
	diag->repulsion_energy = 0.0f;
	diag->sep_energy = 0.0f;
	diag->anchor_energy = 0.0f;
	diag->total_energy = 0.0f;
	diag->force_l1 = 0.0f;
	diag->backbone_force_l1 = 0.0f;
	diag->repulsion_force_l1 = 0.0f;
	diag->sep_force_l1 = 0.0f;
	diag->anchor_force_l1 = 0.0f;
	diag->n_force_nonfinite = 0;
	diag->n_contact_nonfinite = 0;
	diag->n_backbone_nonfinite = 0;
	diag->n_repulsion_nonfinite = 0;
	diag->n_sep_nonfinite = 0;
	diag->n_anchor_nonfinite = 0;
	diag->n_backbone_edges = 0;
	diag->n_repulsion_pairs_considered = 0;
	diag->n_repulsion_pairs_blocked = 0;
	diag->n_repulsion_pairs_active = 0;
	diag->repulsion_mode = HK_BLIND_REPULSION_NONE;
}

static int hk_blind_relax_step_cpu_impl(const struct hk_fdg_conf *conf, const struct hk_blind_wedge_list *edges,
										const struct hk_bmap *bmap_or_null, int32_t n_haploid, fvec3_t *coords, float unit, float step,
										float min_sep_unit, float lambda_sep, int enable_repulsion, int repulsion_mode, float rel_rep_k,
										float repulsion_block_k_min,
										float chr_sep_unit, float lambda_chr_sep,
										const struct hk_blind_coarse_to_fine_map *anchor_map,
										const fvec3_t *coarse_diploid_coords, float anchor_k,
										struct hk_blind_step_diag *diag)
{
	fvec3_t *force = 0;
	struct hk_blind_backbone_diag backbone_diag;
	struct hk_blind_repulsion_diag repulsion_diag;
	float contact_energy, backbone_energy = 0.0f, repulsion_energy = 0.0f, sep_energy;
	float sep_force_l1 = 0.0f, anchor_energy = 0.0f, anchor_force_l1 = 0.0f;
	double force_l1 = 0.0;
	int32_t n_diploid;
	int32_t n_contact_nonfinite = 0;
	int32_t n_sep_nonfinite = 0;
	int32_t n_anchor_nonfinite = 0;
	int32_t n_force_nonfinite = 0;
	int32_t i;
	int a;

	assert(conf);
	assert(edges);
	assert(edges->n_edges >= 0);
	assert(edges->n_edges == 0 || edges->edges);
	assert(n_haploid >= 0);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(isfinite(step));
	assert(step >= 0.0f);
	assert(isfinite(min_sep_unit));
	assert(min_sep_unit >= 0.0f);
	assert(isfinite(lambda_sep));
	assert(lambda_sep >= 0.0f);
	assert(enable_repulsion == 0 || enable_repulsion == 1);
	assert(hk_blind_repulsion_mode_valid(repulsion_mode));
	assert(!enable_repulsion || repulsion_mode == HK_BLIND_REPULSION_N2 || repulsion_mode == HK_BLIND_REPULSION_CELL);
	assert(isfinite(rel_rep_k));
	assert(rel_rep_k >= 0.0f);
	assert(isfinite(repulsion_block_k_min));
	assert(repulsion_block_k_min >= 0.0f);
	assert(isfinite(chr_sep_unit));
	assert(chr_sep_unit >= 0.0f);
	assert(isfinite(lambda_chr_sep));
	assert(lambda_chr_sep >= 0.0f);
	assert(isfinite(anchor_k));
	assert(anchor_k >= 0.0f);
	if (bmap_or_null) {
		assert(bmap_or_null->n_beads == n_haploid);
		assert(n_haploid == 0 || bmap_or_null->beads);
	}
	if (anchor_map) {
		assert(anchor_map->n_fine == n_haploid);
		assert(anchor_map->n_coarse >= 0);
		assert(coarse_diploid_coords || anchor_map->n_coarse == 0);
	} else {
		assert(anchor_k == 0.0f);
	}

	assert(n_haploid <= INT32_MAX / HK_DIPLOID_N_COPY);
	n_diploid = n_haploid * HK_DIPLOID_N_COPY;
	if (n_diploid > 0)
		assert(coords);
	else
		assert(edges->n_edges == 0);

	if (diag)
		hk_blind_step_diag_init(diag);

	if (n_diploid > 0)
		force = CALLOC(fvec3_t, n_diploid);

	contact_energy = hk_blind_wedge_list_accumulate_contact_force_cpu(conf, edges, coords, n_diploid, unit,
																	  force, &n_contact_nonfinite);
	hk_blind_backbone_diag_init(&backbone_diag);
	if (bmap_or_null)
		backbone_energy = hk_blind_backbone_accumulate_force_cpu(conf, bmap_or_null, coords, unit, force, &backbone_diag);
	hk_blind_repulsion_diag_init(&repulsion_diag);
	if (enable_repulsion) {
		if (repulsion_mode == HK_BLIND_REPULSION_CELL)
			repulsion_energy = hk_blind_repulsion_accumulate_force_cell_cpu_impl(conf, n_diploid, coords, unit, rel_rep_k,
																				 repulsion_block_k_min,
																				 edges, bmap_or_null, force, &repulsion_diag);
		else
			repulsion_energy = hk_blind_repulsion_accumulate_force_cpu_impl(conf, n_diploid, coords, unit, rel_rep_k,
																			repulsion_block_k_min,
																			edges, bmap_or_null, force, &repulsion_diag);
	}
	sep_energy = hk_blind_homolog_sep_accumulate_force_ex(n_haploid, coords, force,
														  unit, min_sep_unit, lambda_sep,
														  &sep_force_l1, &n_sep_nonfinite);
	if (bmap_or_null && chr_sep_unit > 0.0f && lambda_chr_sep > 0.0f) {
		float chr_sep_force_l1 = 0.0f;
		int32_t chr_sep_nonfinite = 0;
		sep_energy += hk_blind_chr_centroid_sep_accumulate_force(bmap_or_null, coords, force,
																 unit, chr_sep_unit, lambda_chr_sep,
																 &chr_sep_force_l1, &chr_sep_nonfinite);
		sep_force_l1 += chr_sep_force_l1;
		n_sep_nonfinite += chr_sep_nonfinite;
	}
	if (anchor_map && anchor_k > 0.0f)
		anchor_energy = hk_blind_parent_centroid_anchor_accumulate_force(anchor_map, coords, coarse_diploid_coords,
																		 anchor_k, force, &n_anchor_nonfinite,
																		 &anchor_force_l1);

	for (i = 0; i < n_diploid; ++i) {
		for (a = 0; a < 3; ++a) {
			float f = force[i][a];
			if (hk_blind_float_isfinite(f)) {
				force_l1 += fabs((double)f);
				if (step > 0.0f)
					coords[i][a] += step * f;
			} else {
				++n_force_nonfinite;
			}
		}
	}

	if (diag) {
		diag->contact_energy = contact_energy;
		diag->backbone_energy = backbone_energy;
		diag->repulsion_energy = repulsion_energy;
		diag->sep_energy = sep_energy;
		diag->anchor_energy = anchor_energy;
		diag->total_energy = contact_energy + backbone_energy + repulsion_energy + sep_energy + anchor_energy;
			diag->force_l1 = (float)force_l1;
			diag->backbone_force_l1 = backbone_diag.force_l1;
			diag->repulsion_force_l1 = repulsion_diag.force_l1;
			diag->sep_force_l1 = sep_force_l1;
			diag->anchor_force_l1 = anchor_force_l1;
			diag->n_force_nonfinite = n_force_nonfinite;
			diag->n_contact_nonfinite = n_contact_nonfinite;
			diag->n_backbone_nonfinite = backbone_diag.n_nonfinite;
			diag->n_repulsion_nonfinite = repulsion_diag.n_nonfinite;
			diag->n_sep_nonfinite = n_sep_nonfinite;
			diag->n_anchor_nonfinite = n_anchor_nonfinite;
		diag->n_backbone_edges = backbone_diag.n_edges;
		diag->n_repulsion_pairs_considered = repulsion_diag.n_pairs_considered;
		diag->n_repulsion_pairs_blocked = repulsion_diag.n_pairs_blocked;
		diag->n_repulsion_pairs_active = repulsion_diag.n_pairs_active;
		diag->repulsion_mode = enable_repulsion? repulsion_mode : HK_BLIND_REPULSION_NONE;
	}
	free(force);
	return 0;
}

int hk_blind_relax_step_cpu(const struct hk_fdg_conf *conf, const struct hk_blind_wedge_list *edges,
							const struct hk_bmap *bmap_or_null, int32_t n_haploid, fvec3_t *coords, float unit, float step,
							float min_sep_unit, float lambda_sep, int enable_repulsion, int repulsion_mode, float rel_rep_k,
							struct hk_blind_step_diag *diag)
{
	return hk_blind_relax_step_cpu_impl(conf, edges, bmap_or_null, n_haploid, coords, unit, step,
										min_sep_unit, lambda_sep, enable_repulsion, repulsion_mode, rel_rep_k,
										0.0f, 0.0f, 0.0f, 0, 0, 0.0f, diag);
}

int hk_blind_relax_step_parent_anchor_cpu(const struct hk_fdg_conf *conf, const struct hk_blind_wedge_list *edges,
										  const struct hk_bmap *bmap_or_null, int32_t n_haploid, fvec3_t *coords,
										  float unit, float step, float min_sep_unit, float lambda_sep,
										  int enable_repulsion, int repulsion_mode, float rel_rep_k,
										  const struct hk_blind_coarse_to_fine_map *anchor_map,
										  const fvec3_t *coarse_diploid_coords, float anchor_k,
										  struct hk_blind_step_diag *diag)
{
	return hk_blind_relax_step_cpu_impl(conf, edges, bmap_or_null, n_haploid, coords, unit, step,
										min_sep_unit, lambda_sep, enable_repulsion, repulsion_mode, rel_rep_k,
										0.0f, 0.0f, 0.0f, anchor_map, coarse_diploid_coords, anchor_k, diag);
}

void hk_blind_relax_diag_init(struct hk_blind_relax_diag *diag)
{
	assert(diag);
	diag->n_steps = 0;
	diag->n_completed = 0;
	diag->initial_total_energy = 0.0f;
	diag->final_total_energy = 0.0f;
	diag->initial_contact_energy = 0.0f;
	diag->final_contact_energy = 0.0f;
	diag->initial_backbone_energy = 0.0f;
	diag->final_backbone_energy = 0.0f;
	diag->initial_repulsion_energy = 0.0f;
	diag->final_repulsion_energy = 0.0f;
	diag->initial_sep_energy = 0.0f;
	diag->final_sep_energy = 0.0f;
	diag->initial_anchor_energy = 0.0f;
	diag->final_anchor_energy = 0.0f;
	diag->max_force_l1 = 0.0f;
	diag->final_force_l1 = 0.0f;
	diag->max_backbone_force_l1 = 0.0f;
	diag->final_backbone_force_l1 = 0.0f;
	diag->max_repulsion_force_l1 = 0.0f;
	diag->final_repulsion_force_l1 = 0.0f;
	diag->max_sep_force_l1 = 0.0f;
	diag->final_sep_force_l1 = 0.0f;
	diag->max_anchor_force_l1 = 0.0f;
	diag->final_anchor_force_l1 = 0.0f;
	diag->final_n_repulsion_pairs_considered = 0;
	diag->final_n_repulsion_pairs_blocked = 0;
	diag->final_n_repulsion_pairs_active = 0;
	diag->n_nonfinite_step = 0;
	diag->n_backbone_nonfinite_step = 0;
	diag->n_repulsion_nonfinite_step = 0;
	diag->n_sep_nonfinite_step = 0;
	diag->n_anchor_nonfinite_step = 0;
	diag->n_coord_nonfinite = 0;
	diag->repulsion_mode = HK_BLIND_REPULSION_NONE;
}

static int32_t hk_blind_count_nonfinite_coords(const fvec3_t *coords, int32_t n_diploid)
{
	int32_t n = 0;
	int32_t i;
	int a;

	assert(n_diploid >= 0);
	if (n_diploid == 0)
		return 0;
	assert(coords);
	for (i = 0; i < n_diploid; ++i)
		for (a = 0; a < 3; ++a)
			if (!hk_blind_float_isfinite(coords[i][a]))
				++n;
	return n;
}

static int hk_blind_step_diag_has_nonfinite(const struct hk_blind_step_diag *diag)
{
	assert(diag);
	return diag->n_contact_nonfinite != 0 ||
		   diag->n_backbone_nonfinite != 0 ||
		   diag->n_repulsion_nonfinite != 0 ||
		   diag->n_sep_nonfinite != 0 ||
		   diag->n_anchor_nonfinite != 0 ||
		   diag->n_force_nonfinite != 0 ||
		   !hk_blind_float_isfinite(diag->contact_energy) ||
		   !hk_blind_float_isfinite(diag->backbone_energy) ||
		   !hk_blind_float_isfinite(diag->repulsion_energy) ||
		   !hk_blind_float_isfinite(diag->sep_energy) ||
		   !hk_blind_float_isfinite(diag->anchor_energy) ||
		   !hk_blind_float_isfinite(diag->total_energy) ||
		   !hk_blind_float_isfinite(diag->backbone_force_l1) ||
		   !hk_blind_float_isfinite(diag->repulsion_force_l1) ||
		   !hk_blind_float_isfinite(diag->sep_force_l1) ||
		   !hk_blind_float_isfinite(diag->anchor_force_l1) ||
		   !hk_blind_float_isfinite(diag->force_l1);
}

static int hk_blind_relax_cpu_impl(const struct hk_fdg_conf *conf, const struct hk_blind_wedge_list *edges,
								   const struct hk_bmap *bmap_or_null, int32_t n_haploid, fvec3_t *coords,
								   float unit, float step, int32_t n_steps, float min_sep_unit, float lambda_sep,
								   int enable_repulsion, int repulsion_mode, float repulsion_block_k_min,
								   float chr_sep_unit, float lambda_chr_sep,
								   const struct hk_blind_coarse_to_fine_map *anchor_map,
								   const fvec3_t *coarse_diploid_coords, float anchor_k,
								   struct hk_blind_relax_diag *diag)
{
	int32_t n_diploid;
	int32_t t;

	assert(conf);
	assert(edges);
	assert(edges->n_edges >= 0);
	assert(edges->n_edges == 0 || edges->edges);
	assert(n_haploid >= 0);
	assert(n_haploid <= INT32_MAX / HK_DIPLOID_N_COPY);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(isfinite(step));
	assert(step >= 0.0f);
	assert(n_steps >= 0);
	assert(isfinite(min_sep_unit));
	assert(min_sep_unit >= 0.0f);
	assert(isfinite(lambda_sep));
	assert(lambda_sep >= 0.0f);
	assert(enable_repulsion == 0 || enable_repulsion == 1);
	assert(hk_blind_repulsion_mode_valid(repulsion_mode));
	assert(!enable_repulsion || repulsion_mode == HK_BLIND_REPULSION_N2 || repulsion_mode == HK_BLIND_REPULSION_CELL);
	assert(isfinite(repulsion_block_k_min));
	assert(repulsion_block_k_min >= 0.0f);
	assert(isfinite(chr_sep_unit));
	assert(chr_sep_unit >= 0.0f);
	assert(isfinite(lambda_chr_sep));
	assert(lambda_chr_sep >= 0.0f);
	assert(isfinite(anchor_k));
	assert(anchor_k >= 0.0f);
	assert(diag);
	if (bmap_or_null) {
		assert(bmap_or_null->n_beads == n_haploid);
		assert(n_haploid == 0 || bmap_or_null->beads);
	}
	if (anchor_map) {
		assert(anchor_map->n_fine == n_haploid);
		assert(coarse_diploid_coords || anchor_map->n_coarse == 0);
	} else {
		assert(anchor_k == 0.0f);
	}

	n_diploid = n_haploid * HK_DIPLOID_N_COPY;
	if (n_diploid > 0)
		assert(coords);
	else
		assert(edges->n_edges == 0);

	hk_blind_relax_diag_init(diag);
	diag->n_steps = n_steps;
	diag->repulsion_mode = enable_repulsion? repulsion_mode : HK_BLIND_REPULSION_NONE;
	diag->n_coord_nonfinite = hk_blind_count_nonfinite_coords(coords, n_diploid);
	if (diag->n_coord_nonfinite != 0) {
		++diag->n_nonfinite_step;
		return -1;
	}
	if (n_steps == 0)
		return 0;

	for (t = 0; t < n_steps; ++t) {
		struct hk_blind_step_diag step_diag;
		int ret;

		ret = hk_blind_relax_step_cpu_impl(conf, edges, bmap_or_null, n_haploid, coords, unit, step, min_sep_unit, lambda_sep,
											   enable_repulsion, repulsion_mode,
											   enable_repulsion? hk_blind_rel_rep_schedule_at(t, n_steps) : 0.0f,
											   repulsion_block_k_min,
											   chr_sep_unit, lambda_chr_sep,
											   anchor_map, coarse_diploid_coords, anchor_k,
											   &step_diag);
		if (t == 0) {
			diag->initial_total_energy = step_diag.total_energy;
			diag->initial_contact_energy = step_diag.contact_energy;
			diag->initial_backbone_energy = step_diag.backbone_energy;
			diag->initial_repulsion_energy = step_diag.repulsion_energy;
			diag->initial_sep_energy = step_diag.sep_energy;
			diag->initial_anchor_energy = step_diag.anchor_energy;
		}
		diag->final_total_energy = step_diag.total_energy;
		diag->final_contact_energy = step_diag.contact_energy;
		diag->final_backbone_energy = step_diag.backbone_energy;
		diag->final_repulsion_energy = step_diag.repulsion_energy;
		diag->final_sep_energy = step_diag.sep_energy;
		diag->final_anchor_energy = step_diag.anchor_energy;
		diag->final_force_l1 = step_diag.force_l1;
		diag->final_backbone_force_l1 = step_diag.backbone_force_l1;
		diag->final_repulsion_force_l1 = step_diag.repulsion_force_l1;
		diag->final_sep_force_l1 = step_diag.sep_force_l1;
		diag->final_anchor_force_l1 = step_diag.anchor_force_l1;
		diag->final_n_repulsion_pairs_considered = step_diag.n_repulsion_pairs_considered;
		diag->final_n_repulsion_pairs_blocked = step_diag.n_repulsion_pairs_blocked;
		diag->final_n_repulsion_pairs_active = step_diag.n_repulsion_pairs_active;
		if (hk_blind_float_isfinite(step_diag.force_l1) && step_diag.force_l1 > diag->max_force_l1)
			diag->max_force_l1 = step_diag.force_l1;
		if (hk_blind_float_isfinite(step_diag.backbone_force_l1) && step_diag.backbone_force_l1 > diag->max_backbone_force_l1)
			diag->max_backbone_force_l1 = step_diag.backbone_force_l1;
		if (hk_blind_float_isfinite(step_diag.repulsion_force_l1) && step_diag.repulsion_force_l1 > diag->max_repulsion_force_l1)
			diag->max_repulsion_force_l1 = step_diag.repulsion_force_l1;
		if (hk_blind_float_isfinite(step_diag.sep_force_l1) && step_diag.sep_force_l1 > diag->max_sep_force_l1)
			diag->max_sep_force_l1 = step_diag.sep_force_l1;
		if (hk_blind_float_isfinite(step_diag.anchor_force_l1) && step_diag.anchor_force_l1 > diag->max_anchor_force_l1)
			diag->max_anchor_force_l1 = step_diag.anchor_force_l1;
		if (step_diag.n_backbone_nonfinite != 0)
			++diag->n_backbone_nonfinite_step;
		if (step_diag.n_repulsion_nonfinite != 0)
			++diag->n_repulsion_nonfinite_step;
		if (step_diag.n_sep_nonfinite != 0)
			++diag->n_sep_nonfinite_step;
		if (step_diag.n_anchor_nonfinite != 0)
			++diag->n_anchor_nonfinite_step;

		diag->n_coord_nonfinite = hk_blind_count_nonfinite_coords(coords, n_diploid);
		if (ret != 0 || hk_blind_step_diag_has_nonfinite(&step_diag) || diag->n_coord_nonfinite != 0) {
			++diag->n_nonfinite_step;
			return -1;
		}
		++diag->n_completed;
	}
	return 0;
}

static int hk_blind_gpu_pair_push(struct hk_fdg_gpu_pair **pairs, size_t *n_pairs, size_t *m_pairs,
								  int32_t i, int32_t j, float k, float d_scale, uint8_t type)
{
	struct hk_fdg_gpu_pair *p;
	size_t new_m;
	assert(pairs);
	assert(n_pairs);
	assert(m_pairs);
	if (i == j || k <= 0.0f)
		return 0;
	if (!isfinite(k) || !isfinite(d_scale) || d_scale <= 0.0f) {
		fprintf(stderr, "[E::blind-gpu] bad GPU pair i=%d j=%d type=%u k=%.9g d_scale=%.9g\n",
				i, j, (unsigned)type, k, d_scale);
		return -1;
	}
	if (*n_pairs == *m_pairs) {
		new_m = *m_pairs? (*m_pairs << 1) : 1024;
		p = (struct hk_fdg_gpu_pair*)realloc(*pairs, new_m * sizeof(**pairs));
		if (p == 0) {
			fprintf(stderr, "[E::blind-gpu] failed to grow GPU pair buffer to %zu pairs\n", new_m);
			return -1;
		}
		*pairs = p;
		*m_pairs = new_m;
	}
	p = &(*pairs)[(*n_pairs)++];
	p->i = i;
	p->j = j;
	p->k = k;
	p->d_scale = d_scale;
	p->type = type;
	p->pad[0] = p->pad[1] = p->pad[2] = 0;
	return 0;
}

static int hk_blind_gpu_build_pairs(const struct hk_bmap *bmap_or_null, const struct hk_blind_wedge_list *edges,
									int32_t n_diploid, struct hk_fdg_gpu_pair **pairs_out,
									size_t *n_pairs_out, uint64_t **block_keys_out, size_t *n_block_keys_out)
{
	struct hk_fdg_gpu_pair *pairs = 0;
	uint64_t *block_keys = 0;
	size_t n_pairs = 0, m_pairs = 0;
	int64_t n_blocked = 0;
	int32_t i, eidx;

	assert(edges);
	assert(pairs_out);
	assert(n_pairs_out);
	assert(block_keys_out);
	assert(n_block_keys_out);
	*pairs_out = 0;
	*n_pairs_out = 0;
	*block_keys_out = 0;
	*n_block_keys_out = 0;

	if (bmap_or_null && bmap_or_null->n_beads > 0) {
		int32_t mid_dist = hk_blind_bmap_mid_bead_size(bmap_or_null);
		for (i = 1; i < bmap_or_null->n_beads; ++i) {
			const struct hk_bead *prev = &bmap_or_null->beads[i - 1];
			const struct hk_bead *cur = &bmap_or_null->beads[i];
			float d_scale;
			int copy;
			if (prev->chr != cur->chr)
				continue;
			d_scale = hk_blind_backbone_d_scale(prev, cur, mid_dist);
			for (copy = 0; copy < HK_DIPLOID_N_COPY; ++copy) {
				if (hk_blind_gpu_pair_push(&pairs, &n_pairs, &m_pairs,
										   hk_diploid_bid(i - 1, copy), hk_diploid_bid(i, copy),
										   1.0f, d_scale, HK_FDG_PAIR_TYPE_BACKBONE) != 0)
					goto fail;
			}
		}
	}
	for (eidx = 0; eidx < edges->n_edges; ++eidx) {
		const struct hk_blind_wedge *edge = &edges->edges[eidx];
		assert(edge->bid[0] >= 0 && edge->bid[0] < n_diploid);
		assert(edge->bid[1] >= 0 && edge->bid[1] < n_diploid);
		if (hk_blind_gpu_pair_push(&pairs, &n_pairs, &m_pairs,
								   edge->bid[0], edge->bid[1], edge->k, edge->d_scale,
								   HK_FDG_PAIR_TYPE_CONTACT) != 0)
			goto fail;
	}
	n_blocked = hk_blind_build_repulsion_blocked_pairs(edges, bmap_or_null, n_diploid, 0.0f, &block_keys);
	if (n_blocked < 0)
		goto fail;
	if (n_blocked > 0) {
		uint64_t *expanded = MALLOC(uint64_t, (size_t)n_blocked * 2);
		int64_t k;
		if (expanded == 0)
			goto fail;
		for (k = 0; k < n_blocked; ++k) {
			uint32_t i0 = (uint32_t)(block_keys[k] >> 32);
			uint32_t i1 = (uint32_t)block_keys[k];
			expanded[2 * k] = ((uint64_t)i0 << 32) | (uint64_t)i1;
			expanded[2 * k + 1] = ((uint64_t)i1 << 32) | (uint64_t)i0;
		}
		free(block_keys);
		block_keys = expanded;
		n_blocked *= 2;
	}
	*pairs_out = pairs;
	*n_pairs_out = n_pairs;
	*block_keys_out = block_keys;
	*n_block_keys_out = (size_t)n_blocked;
	return 0;

fail:
	free(block_keys);
	free(pairs);
	return -1;
}

static int hk_blind_apply_extra_cpu_forces(const struct hk_blind_coarse_to_fine_map *anchor_map,
										   const fvec3_t *coarse_diploid_coords, float anchor_k,
										   int32_t n_haploid, fvec3_t *coords, float unit, float step,
										   float min_sep_unit, float lambda_sep,
										   fvec3_t *extra_force, struct hk_blind_step_diag *step_diag)
{
	float sep_force_l1 = 0.0f, anchor_force_l1 = 0.0f;
	int32_t n_sep_nonfinite = 0, n_anchor_nonfinite = 0;
	int32_t n_diploid = n_haploid * HK_DIPLOID_N_COPY;
	int32_t i;
	int a;

	assert(coords || n_diploid == 0);
	assert(extra_force || n_diploid == 0);
	assert(step_diag);
	if (n_diploid == 0)
		return 0;
	step_diag->sep_energy = hk_blind_homolog_sep_accumulate_force_ex(n_haploid, coords, extra_force,
																	 unit, min_sep_unit, lambda_sep,
																	 &sep_force_l1, &n_sep_nonfinite);
	step_diag->sep_force_l1 = sep_force_l1;
	step_diag->n_sep_nonfinite = n_sep_nonfinite;
	if (anchor_map && anchor_k > 0.0f) {
		step_diag->anchor_energy = hk_blind_parent_centroid_anchor_accumulate_force(anchor_map, coords,
																					coarse_diploid_coords,
																					anchor_k, extra_force,
																					&n_anchor_nonfinite,
																					&anchor_force_l1);
		step_diag->anchor_force_l1 = anchor_force_l1;
		step_diag->n_anchor_nonfinite = n_anchor_nonfinite;
	}
	for (i = 0; i < n_diploid; ++i) {
		for (a = 0; a < 3; ++a) {
			float f = extra_force[i][a];
			if (hk_blind_float_isfinite(f)) {
				step_diag->force_l1 += fabsf(f);
				if (step > 0.0f)
					coords[i][a] += step * f;
			} else {
				++step_diag->n_force_nonfinite;
			}
		}
	}
	return 0;
}

static int hk_blind_relax_gpu_impl(const struct hk_fdg_conf *conf, const struct hk_blind_wedge_list *edges,
								   const struct hk_bmap *bmap_or_null, int32_t n_haploid, fvec3_t *coords,
								   float unit, float step, int32_t n_steps, float min_sep_unit, float lambda_sep,
								   int enable_repulsion, int repulsion_mode, float repulsion_block_k_min,
								   const struct hk_blind_coarse_to_fine_map *anchor_map,
								   const fvec3_t *coarse_diploid_coords, float anchor_k,
								   struct hk_blind_relax_diag *diag)
{
	struct hk_fdg_gpu_pair *gpu_pairs = 0;
	uint64_t *block_keys = 0;
	struct hk_fdg_gpu_ctx *gpu_ctx = 0;
	size_t n_gpu_pairs = 0, n_block_keys = 0;
	int32_t n_diploid;
	int32_t t;
	int ret = -1;

	assert(conf);
	assert(edges);
	assert(n_haploid >= 0);
	assert(n_haploid <= INT32_MAX / HK_DIPLOID_N_COPY);
	assert(isfinite(unit));
	assert(unit > 0.0f);
	assert(isfinite(step));
	assert(step >= 0.0f);
	assert(n_steps >= 0);
	assert(enable_repulsion == 0 || enable_repulsion == 1);
	assert(repulsion_mode == HK_BLIND_REPULSION_CELL);
	assert(repulsion_block_k_min == 0.0f);
	assert(diag);
	if (bmap_or_null) {
		assert(bmap_or_null->n_beads == n_haploid);
		assert(n_haploid == 0 || bmap_or_null->beads);
	}
	if (anchor_map) {
		assert(anchor_map->n_fine == n_haploid);
		assert(coarse_diploid_coords || anchor_map->n_coarse == 0);
	} else {
		assert(anchor_k == 0.0f);
	}
	n_diploid = n_haploid * HK_DIPLOID_N_COPY;
	if (n_diploid > 0)
		assert(coords);
	else
		assert(edges->n_edges == 0);

	hk_blind_relax_diag_init(diag);
	diag->n_steps = n_steps;
	diag->repulsion_mode = enable_repulsion? repulsion_mode : HK_BLIND_REPULSION_NONE;
	diag->n_coord_nonfinite = hk_blind_count_nonfinite_coords(coords, n_diploid);
	if (diag->n_coord_nonfinite != 0) {
		++diag->n_nonfinite_step;
		return -1;
	}
	if (n_steps == 0)
		return 0;
	if (!hk_fdg_gpu_is_available()) {
		fprintf(stderr, "[E::blind-gpu] GPU backend requested but CUDA backend/device is unavailable; rebuild with make gpu=1 and check CUDA device visibility\n");
		return -1;
	}
	if (hk_blind_gpu_build_pairs(bmap_or_null, edges, n_diploid, &gpu_pairs, &n_gpu_pairs,
								 &block_keys, &n_block_keys) != 0) {
		fprintf(stderr, "[E::blind-gpu] build GPU pairs failed: n_diploid=%d n_edges=%d\n",
				n_diploid, edges->n_edges);
		goto cleanup;
	}
	fprintf(stderr, "[M::blind-gpu] prepared pair buffers: n_diploid=%d n_pairs=%zu n_block_keys=%zu n_steps=%d\n",
			n_diploid, n_gpu_pairs, n_block_keys, n_steps);
	gpu_ctx = hk_fdg_gpu_create(n_diploid);
	if (gpu_ctx == 0) {
		fprintf(stderr, "[E::blind-gpu] create GPU context failed: n_diploid=%d\n", n_diploid);
		goto cleanup;
	}
	if (hk_fdg_gpu_prepare(gpu_ctx, n_diploid, n_block_keys) != 0) {
		fprintf(stderr, "[E::blind-gpu] prepare GPU context failed: n_diploid=%d n_block_keys=%zu\n",
				n_diploid, n_block_keys);
		goto cleanup;
	}
	if (hk_fdg_gpu_set_blocklist(gpu_ctx, block_keys, n_block_keys) != 0) {
		fprintf(stderr, "[E::blind-gpu] set GPU blocklist failed: n_block_keys=%zu\n", n_block_keys);
		goto cleanup;
	}
	if (hk_fdg_gpu_upload_positions(gpu_ctx, (const fvec3_t*)coords, n_diploid) != 0) {
		fprintf(stderr, "[E::blind-gpu] upload GPU positions failed: n_diploid=%d\n", n_diploid);
		goto cleanup;
	}

	for (t = 0; t < n_steps; ++t) {
		struct hk_blind_step_diag step_diag;
		struct hk_fdg_conf gpu_conf = *conf;
		struct hk_fdg_gpu_stats stats;
		fvec3_t *extra_force = 0;
		double rms_force = 0.0;
		int need_sync = 1;

		hk_blind_step_diag_init(&step_diag);
		gpu_conf.step = step / unit;
		gpu_conf.coef_moment = 0.0f;
		gpu_conf.max_f = 0.0f;
		if (!enable_repulsion)
			gpu_conf.k_rel_rep = 0.0f;
		if (hk_fdg_gpu_compute(gpu_ctx, &gpu_conf, t == 0? gpu_pairs : 0, n_gpu_pairs,
							   unit,
							   enable_repulsion? hk_blind_rel_rep_schedule_at(t, n_steps) : 0.0f,
							   enable_repulsion? gpu_conf.d_r : 0.0f,
							   &stats, &rms_force, need_sync) != 0) {
			fprintf(stderr, "[E::blind-gpu] GPU compute failed at relax_step=%d/%d n_pairs=%zu n_diploid=%d rel_rep=%.9g rep_radius=%.9g\n",
					t + 1, n_steps, n_gpu_pairs, n_diploid,
					enable_repulsion? hk_blind_rel_rep_schedule_at(t, n_steps) : 0.0f,
					enable_repulsion? gpu_conf.d_r : 0.0f);
			goto cleanup;
		}
		if (!isfinite(rms_force) ||
			!isfinite(stats.energy[HK_FDG_PAIR_TYPE_CONTACT]) ||
			!isfinite(stats.energy[HK_FDG_PAIR_TYPE_BACKBONE]) ||
			!isfinite(stats.energy[HK_FDG_PAIR_TYPE_REPEL])) {
			fprintf(stderr, "[E::blind-gpu] non-finite GPU stats at relax_step=%d/%d rms=%.9g energy=%.9g,%.9g,%.9g\n",
					t + 1, n_steps, rms_force,
					stats.energy[HK_FDG_PAIR_TYPE_CONTACT],
					stats.energy[HK_FDG_PAIR_TYPE_BACKBONE],
					stats.energy[HK_FDG_PAIR_TYPE_REPEL]);
			goto cleanup;
		}
		if (hk_fdg_gpu_download_positions(gpu_ctx, coords, n_diploid) != 0) {
			fprintf(stderr, "[E::blind-gpu] download GPU positions failed at relax_step=%d/%d n_diploid=%d\n",
					t + 1, n_steps, n_diploid);
			goto cleanup;
		}
		step_diag.contact_energy = stats.energy[HK_FDG_PAIR_TYPE_CONTACT];
		step_diag.backbone_energy = stats.energy[HK_FDG_PAIR_TYPE_BACKBONE];
		step_diag.repulsion_energy = stats.energy[HK_FDG_PAIR_TYPE_REPEL];
		step_diag.force_l1 = (float)(rms_force * sqrt((double)n_diploid));
		step_diag.n_backbone_edges = stats.active[HK_FDG_PAIR_TYPE_BACKBONE];
		step_diag.n_repulsion_pairs_active = stats.active[HK_FDG_PAIR_TYPE_REPEL];
		step_diag.n_repulsion_pairs_considered = stats.active[HK_FDG_PAIR_TYPE_REPEL];
		step_diag.n_repulsion_pairs_blocked = 0;
		step_diag.repulsion_mode = enable_repulsion? repulsion_mode : HK_BLIND_REPULSION_NONE;
		if (min_sep_unit > 0.0f || lambda_sep > 0.0f || (anchor_map && anchor_k > 0.0f)) {
			extra_force = n_diploid > 0? CALLOC(fvec3_t, n_diploid) : 0;
			if (n_diploid > 0 && extra_force == 0)
				goto cleanup;
			if (hk_blind_apply_extra_cpu_forces(anchor_map, coarse_diploid_coords, anchor_k,
												n_haploid, coords, unit, step,
												min_sep_unit, lambda_sep,
												extra_force, &step_diag) != 0) {
				free(extra_force);
				goto cleanup;
			}
			free(extra_force);
				if (hk_fdg_gpu_upload_positions(gpu_ctx, (const fvec3_t*)coords, n_diploid) != 0) {
					fprintf(stderr, "[E::blind-gpu] upload GPU positions after extra forces failed at relax_step=%d/%d\n",
							t + 1, n_steps);
					goto cleanup;
				}
		}
		step_diag.total_energy = step_diag.contact_energy + step_diag.backbone_energy +
			step_diag.repulsion_energy + step_diag.sep_energy + step_diag.anchor_energy;
		if (t == 0) {
			diag->initial_total_energy = step_diag.total_energy;
			diag->initial_contact_energy = step_diag.contact_energy;
			diag->initial_backbone_energy = step_diag.backbone_energy;
			diag->initial_repulsion_energy = step_diag.repulsion_energy;
			diag->initial_sep_energy = step_diag.sep_energy;
			diag->initial_anchor_energy = step_diag.anchor_energy;
		}
		diag->final_total_energy = step_diag.total_energy;
		diag->final_contact_energy = step_diag.contact_energy;
		diag->final_backbone_energy = step_diag.backbone_energy;
		diag->final_repulsion_energy = step_diag.repulsion_energy;
		diag->final_sep_energy = step_diag.sep_energy;
		diag->final_anchor_energy = step_diag.anchor_energy;
		diag->final_force_l1 = step_diag.force_l1;
		diag->final_backbone_force_l1 = step_diag.backbone_force_l1;
		diag->final_repulsion_force_l1 = step_diag.repulsion_force_l1;
		diag->final_sep_force_l1 = step_diag.sep_force_l1;
		diag->final_anchor_force_l1 = step_diag.anchor_force_l1;
		diag->final_n_repulsion_pairs_considered = step_diag.n_repulsion_pairs_considered;
		diag->final_n_repulsion_pairs_blocked = step_diag.n_repulsion_pairs_blocked;
		diag->final_n_repulsion_pairs_active = step_diag.n_repulsion_pairs_active;
		if (hk_blind_float_isfinite(step_diag.force_l1) && step_diag.force_l1 > diag->max_force_l1)
			diag->max_force_l1 = step_diag.force_l1;
		if (hk_blind_float_isfinite(step_diag.anchor_force_l1) && step_diag.anchor_force_l1 > diag->max_anchor_force_l1)
			diag->max_anchor_force_l1 = step_diag.anchor_force_l1;
		if (step_diag.n_backbone_nonfinite != 0)
			++diag->n_backbone_nonfinite_step;
		if (step_diag.n_repulsion_nonfinite != 0)
			++diag->n_repulsion_nonfinite_step;
		if (step_diag.n_sep_nonfinite != 0)
			++diag->n_sep_nonfinite_step;
		if (step_diag.n_anchor_nonfinite != 0)
			++diag->n_anchor_nonfinite_step;
		diag->n_coord_nonfinite = hk_blind_count_nonfinite_coords(coords, n_diploid);
		if (hk_blind_step_diag_has_nonfinite(&step_diag) || diag->n_coord_nonfinite != 0) {
			++diag->n_nonfinite_step;
			goto cleanup;
		}
		++diag->n_completed;
	}
	ret = 0;

cleanup:
	if (gpu_ctx) hk_fdg_gpu_destroy(gpu_ctx);
	free(block_keys);
	free(gpu_pairs);
	return ret;
}

static int hk_blind_relax_impl(const struct hk_fdg_conf *conf, const struct hk_blind_wedge_list *edges,
							   const struct hk_bmap *bmap_or_null, int32_t n_haploid, fvec3_t *coords,
							   float unit, float step, int32_t n_steps, float min_sep_unit, float lambda_sep,
							   int enable_repulsion, int repulsion_mode, float repulsion_block_k_min,
							   float chr_sep_unit, float lambda_chr_sep,
							   const struct hk_blind_coarse_to_fine_map *anchor_map,
							   const fvec3_t *coarse_diploid_coords, float anchor_k,
							   struct hk_blind_relax_diag *diag)
{
	if (conf->backend == HK_FDG_BACKEND_GPU) {
		if (chr_sep_unit > 0.0f || lambda_chr_sep > 0.0f ||
			repulsion_mode != HK_BLIND_REPULSION_CELL || repulsion_block_k_min != 0.0f)
			return -1;
		return hk_blind_relax_gpu_impl(conf, edges, bmap_or_null, n_haploid, coords, unit, step,
									   n_steps, min_sep_unit, lambda_sep, enable_repulsion,
									   repulsion_mode, repulsion_block_k_min,
									   anchor_map, coarse_diploid_coords, anchor_k, diag);
	}
	return hk_blind_relax_cpu_impl(conf, edges, bmap_or_null, n_haploid, coords, unit, step,
								   n_steps, min_sep_unit, lambda_sep, enable_repulsion,
								   repulsion_mode, repulsion_block_k_min, chr_sep_unit,
								   lambda_chr_sep, anchor_map, coarse_diploid_coords,
								   anchor_k, diag);
}

int hk_blind_relax_cpu(const struct hk_fdg_conf *conf, const struct hk_blind_wedge_list *edges,
					   const struct hk_bmap *bmap_or_null, int32_t n_haploid, fvec3_t *coords, float unit, float step, int32_t n_steps,
					   float min_sep_unit, float lambda_sep, int enable_repulsion, int repulsion_mode, struct hk_blind_relax_diag *diag)
{
	return hk_blind_relax_cpu_impl(conf, edges, bmap_or_null, n_haploid, coords, unit, step, n_steps,
								   min_sep_unit, lambda_sep, enable_repulsion, repulsion_mode,
								   0.0f, 0.0f, 0.0f, 0, 0, 0.0f, diag);
}

int hk_blind_relax_parent_anchor_cpu(const struct hk_fdg_conf *conf, const struct hk_blind_wedge_list *edges,
									 const struct hk_bmap *bmap_or_null, int32_t n_haploid, fvec3_t *coords,
									 float unit, float step, int32_t n_steps, float min_sep_unit, float lambda_sep,
									 int enable_repulsion, int repulsion_mode,
									 const struct hk_blind_coarse_to_fine_map *anchor_map,
									 const fvec3_t *coarse_diploid_coords, float anchor_k,
									 struct hk_blind_relax_diag *diag)
{
	return hk_blind_relax_cpu_impl(conf, edges, bmap_or_null, n_haploid, coords, unit, step, n_steps,
								   min_sep_unit, lambda_sep, enable_repulsion, repulsion_mode,
								   0.0f, 0.0f, 0.0f, anchor_map, coarse_diploid_coords, anchor_k, diag);
}

void hk_blind_wedge_list_init(struct hk_blind_wedge_list *list)
{
	assert(list);
	list->edges = 0;
	list->n_edges = 0;
	list->m_edges = 0;
	list->n_input_bpair = 0;
	list->n_expanded_edges = 0;
	list->n_softall_candidate_pairs = 0;
	list->n_softall_filter1_pairs = 0;
	list->n_softall_filter2_pairs = 0;
	list->n_softall_bmap_pairs = 0;
	list->n_softall_selected_raw = 0;
	list->n_softall_gate_skip_raw = 0;
	list->n_softall_same_bin_skip_raw = 0;
	memset(list->n_softall_state_raw_count, 0, sizeof(list->n_softall_state_raw_count));
	list->n_skipped_self_edges = 0;
	list->n_skipped_same_bin_bpairs = 0;
	list->n_edges_before_aggregation = 0;
	list->n_aggregated_edges_removed = 0;
	list->sum_rho_train_bpair = 0.0;
	list->mean_rho_train_bpair = 0.0f;
	list->min_rho_train_bpair = 0.0f;
	list->max_rho_train_bpair = 0.0f;
}

void hk_blind_wedge_list_destroy(struct hk_blind_wedge_list *list)
{
	if (list == 0) return;
	free(list->edges);
	hk_blind_wedge_list_init(list);
}

static void hk_blind_wedge_list_push(struct hk_blind_wedge_list *list, const struct hk_blind_wedge *edge)
{
	assert(list);
	assert(edge);
	if (list->n_edges == list->m_edges)
		EXPAND(list->edges, list->m_edges);
	list->edges[list->n_edges++] = *edge;
}

static void hk_blind_wedge_canonicalize_bid(struct hk_blind_wedge *edge)
{
	int32_t t;
	assert(edge);
	if (edge->bid[1] < edge->bid[0]) {
		t = edge->bid[0];
		edge->bid[0] = edge->bid[1];
		edge->bid[1] = t;
	}
}

int hk_blind_wedge_list_build_from_bpair_set(struct hk_blind_wedge_list *out, const struct hk_blind_bpair_set *set,
											 float base_k, float base_d_scale, float rho_train)
{
	int32_t i;
	int32_t n_rho = 0;

	assert(out);
	assert(set);
	assert(set->n_bpairs >= 0);
	assert(set->n_bpairs == 0 || set->bpairs);
	assert(base_k >= 0.0f);
	assert(base_d_scale > 0.0f);
	assert(rho_train >= 0.0f);

	out->n_edges = 0;
	out->n_input_bpair = set->n_bpairs;
	out->n_expanded_edges = 0;
	out->n_skipped_self_edges = 0;
	out->n_skipped_same_bin_bpairs = 0;
	out->n_edges_before_aggregation = 0;
	out->n_aggregated_edges_removed = 0;
	out->sum_rho_train_bpair = 0.0;
	out->mean_rho_train_bpair = 0.0f;
	out->min_rho_train_bpair = 0.0f;
	out->max_rho_train_bpair = 0.0f;
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		struct hk_blind_wedge edges[HK_BLIND_N_STATE];
		float rho_eff;
		int s;
		if (hk_blind_bpair_is_same_bin(bp)) {
			++out->n_skipped_same_bin_bpairs;
			continue;
		}
		rho_eff = hk_blind_bpair_effective_rho_train(bp, rho_train, HK_BLIND_RHO_TRAIN_CONSTANT);
		if (n_rho == 0 || rho_eff < out->min_rho_train_bpair) out->min_rho_train_bpair = rho_eff;
		if (n_rho == 0 || rho_eff > out->max_rho_train_bpair) out->max_rho_train_bpair = rho_eff;
		out->sum_rho_train_bpair += rho_eff;
		++n_rho;
		hk_blind_bpair_expand_weighted_edges(bp, base_k, base_d_scale, rho_train, edges);
		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			++out->n_expanded_edges;
			if (edges[s].bid[0] == edges[s].bid[1]) {
				++out->n_skipped_self_edges;
				continue;
			}
			// Zero-weight contacts must not block repulsion between non-contacting copies.
			if (edges[s].k <= 0.0f)
				continue;
			hk_blind_wedge_canonicalize_bid(&edges[s]);
			hk_blind_wedge_list_push(out, &edges[s]);
		}
	}
	if (n_rho > 0)
		out->mean_rho_train_bpair = (float)(out->sum_rho_train_bpair / n_rho);
	out->n_edges_before_aggregation = out->n_edges;
	return 0;
}

int hk_blind_wedge_list_build_from_bpair_set_params(struct hk_blind_wedge_list *out, const struct hk_blind_bpair_set *set,
													float rho_train)
{
	return hk_blind_wedge_list_build_from_bpair_set_params_mode(out, set, rho_train,
																HK_BLIND_RHO_TRAIN_CONSTANT,
																HK_BLIND_D_SCALE_RAW_COUNT,
																1e-6f);
}

int hk_blind_wedge_list_build_from_bpair_set_params_mode(struct hk_blind_wedge_list *out, const struct hk_blind_bpair_set *set,
														 float rho_train, int rho_train_mode,
														 int d_scale_mode, float d_scale_eps_count)
{
	return hk_blind_wedge_list_build_from_bpair_set_params_mode_ex(out, set, rho_train, rho_train_mode,
																   HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR,
																   d_scale_mode, d_scale_eps_count);
}

int hk_blind_wedge_list_build_from_bpair_set_params_mode_ex(struct hk_blind_wedge_list *out,
															const struct hk_blind_bpair_set *set,
															float rho_train, int rho_train_mode,
															float rho_train_floor, int d_scale_mode,
															float d_scale_eps_count)
{
	int32_t i;
	int32_t n_rho = 0;

	assert(out);
	assert(set);
	assert(set->n_bpairs >= 0);
	assert(set->n_bpairs == 0 || set->bpairs);
	assert(rho_train >= 0.0f);
	assert(hk_blind_rho_train_mode_valid(rho_train_mode));
	assert(hk_blind_d_scale_mode_valid(d_scale_mode));
	assert(isfinite(d_scale_eps_count));
	assert(d_scale_eps_count > 0.0f);

	out->n_edges = 0;
	out->n_input_bpair = set->n_bpairs;
	out->n_expanded_edges = 0;
	out->n_skipped_self_edges = 0;
	out->n_skipped_same_bin_bpairs = 0;
	out->n_edges_before_aggregation = 0;
	out->n_aggregated_edges_removed = 0;
	out->sum_rho_train_bpair = 0.0;
	out->mean_rho_train_bpair = 0.0f;
	out->min_rho_train_bpair = 0.0f;
	out->max_rho_train_bpair = 0.0f;
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		struct hk_blind_wedge edges[HK_BLIND_N_STATE];
		float rho_eff;
		int s;
		assert(isfinite(bp->base_k));
		assert(bp->base_k >= 0.0f);
		assert(isfinite(bp->base_d_scale));
		assert(bp->base_d_scale > 0.0f);
		if (hk_blind_bpair_is_same_bin(bp)) {
			++out->n_skipped_same_bin_bpairs;
			continue;
		}
		rho_eff = hk_blind_bpair_effective_rho_train_floor(bp, rho_train, rho_train_mode,
														   rho_train_floor);
		if (n_rho == 0 || rho_eff < out->min_rho_train_bpair) out->min_rho_train_bpair = rho_eff;
		if (n_rho == 0 || rho_eff > out->max_rho_train_bpair) out->max_rho_train_bpair = rho_eff;
		out->sum_rho_train_bpair += rho_eff;
		++n_rho;
		hk_blind_bpair_expand_weighted_edges_mode_ex(bp, bp->base_k, bp->base_d_scale,
													 rho_train, rho_train_mode, rho_train_floor,
													 d_scale_mode, d_scale_eps_count, edges);
		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			++out->n_expanded_edges;
			if (edges[s].bid[0] == edges[s].bid[1]) {
				++out->n_skipped_self_edges;
				continue;
			}
			// Zero-weight contacts must not block repulsion between non-contacting copies.
			if (edges[s].k <= 0.0f)
				continue;
			hk_blind_wedge_canonicalize_bid(&edges[s]);
			hk_blind_wedge_list_push(out, &edges[s]);
		}
	}
	if (n_rho > 0)
		out->mean_rho_train_bpair = (float)(out->sum_rho_train_bpair / n_rho);
	out->n_edges_before_aggregation = out->n_edges;
	return 0;
}

static int hk_blind_cmp_i32(const void *a_, const void *b_)
{
	int32_t a = *(const int32_t*)a_;
	int32_t b = *(const int32_t*)b_;
	return (a > b) - (a < b);
}

static int hk_blind_split_graph_median_nei(const struct hk_bmap *bmap)
{
	int32_t *tmp;
	int32_t i, med;
	if (bmap == 0 || bmap->n_pairs <= 0)
		return 0;
	tmp = CALLOC(int32_t, bmap->n_pairs);
	if (tmp == 0)
		return 0;
	for (i = 0; i < bmap->n_pairs; ++i)
		tmp[i] = bmap->pairs[i].max_nei;
	qsort(tmp, (size_t)bmap->n_pairs, sizeof(*tmp), hk_blind_cmp_i32);
	med = tmp[bmap->n_pairs / 2];
	free(tmp);
	return med;
}

static float hk_blind_native_fdg_k_from_nei(const struct hk_bpair *p, int32_t median_nei)
{
	const double a_third = 1.0 / 3.0;
	assert(p);
	if (median_nei <= 0)
		return 1.0f;
	if (p->max_nei >= median_nei)
		return 1.0f;
	if (p->max_nei <= 0)
		return 0.0f;
	return powf((float)((double)p->max_nei / (double)median_nei), (float)a_third);
}

static int32_t hk_blind_base_chr_from_split_name(const char *name, int *copy,
												 const struct hk_sdict *base_d)
{
	size_t len;
	char buf[256];
	int32_t i;
	assert(name);
	assert(copy);
	assert(base_d);
	len = strlen(name);
	if (len >= 2 && (name[len - 1] == 'a' || name[len - 1] == 'b')) {
		if (len >= sizeof(buf))
			return -1;
		memcpy(buf, name, len - 1);
		buf[len - 1] = 0;
		*copy = name[len - 1] == 'b';
		for (i = 0; i < base_d->n; ++i)
			if (strcmp(base_d->name[i], buf) == 0)
				return i;
		return -1;
	}
	*copy = 0;
	for (i = 0; i < base_d->n; ++i)
		if (strcmp(base_d->name[i], name) == 0)
			return i;
	return -1;
}

static int hk_blind_build_split_lookup(const struct hk_bmap *split_bmap,
									   const struct hk_bmap *base_bmap,
									   int32_t *split_to_diploid)
{
	int32_t i;
	assert(split_bmap);
	assert(base_bmap);
	assert(split_to_diploid || split_bmap->n_beads == 0);
	for (i = 0; i < split_bmap->n_beads; ++i) {
		const struct hk_bead *b = &split_bmap->beads[i];
		const char *split_name = split_bmap->d->name[b->chr];
		int copy = 0;
		int32_t base_chr = hk_blind_base_chr_from_split_name(split_name, &copy, base_bmap->d);
		int32_t base_bid;
		if (base_chr < 0)
			return -1;
		base_bid = hk_bmap_pos2bid(base_bmap, base_chr, b->st);
		split_to_diploid[i] = hk_diploid_bid(base_bid, copy);
	}
	return 0;
}


struct hk_blind_pair_weight_aux {
	int32_t i;
	double weight;
};

static void hk_blind_softall_map_push_pair(struct hk_map *m,
									   int32_t *m_pairs,
									   const struct hk_blind_pair *q,
									   const int32_t *ploidy_XY,
									   const int32_t *old2new,
									   int canonical_state,
									   uint8_t swapped,
									   float phased_prob,
									   float final_phased_prob,
									   struct hk_blind_softall_aux *aux)
{
	struct hk_pair *r;
	int raw_state;
	int raw_copy[2];
	int32_t chr[2];
	assert(m);
	assert(m_pairs);
	assert(q);
	assert(ploidy_XY);
	assert(old2new);
	assert(canonical_state >= 0 && canonical_state < HK_BLIND_N_STATE);
	assert(isfinite(phased_prob));
	assert(isfinite(final_phased_prob));
	assert(aux);
	raw_state = hk_blind_canonical_state_to_raw_state(canonical_state, swapped);
	raw_copy[0] = (raw_state >> 1) & 1;
	raw_copy[1] = raw_state & 1;
	chr[0] = old2new[q->chr[0]] +
		((ploidy_XY[q->chr[0]] >> 8) == 1? 0 : raw_copy[0]);
	chr[1] = old2new[q->chr[1]] +
		((ploidy_XY[q->chr[1]] >> 8) == 1? 0 : raw_copy[1]);
	if (m->n_pairs == *m_pairs)
		EXPAND(m->pairs, *m_pairs);
	if (aux->n_final_phased_prob == aux->m_final_phased_prob)
		EXPAND(aux->final_phased_prob, aux->m_final_phased_prob);
	if (aux->n_final_phased_prob > INT32_MAX)
		abort();
	r = &m->pairs[m->n_pairs++];
	memset(r, 0, sizeof(*r));
	r->n_ctn = (uint32_t)aux->n_final_phased_prob;
	if (chr[0] > chr[1]) {
		r->chr = (uint64_t)chr[1] << 32 | chr[0];
		r->pos = (uint64_t)(uint32_t)q->pos[1] << 32 | (uint32_t)q->pos[0];
		r->strand[0] = q->strand[1];
		r->strand[1] = q->strand[0];
		r->phase[0] = raw_copy[1];
		r->phase[1] = raw_copy[0];
	} else {
		r->chr = (uint64_t)chr[0] << 32 | chr[1];
		r->pos = (uint64_t)(uint32_t)q->pos[0] << 32 | (uint32_t)q->pos[1];
		r->strand[0] = q->strand[0];
		r->strand[1] = q->strand[1];
		r->phase[0] = raw_copy[0];
		r->phase[1] = raw_copy[1];
	}
	r->_.phased_prob = phased_prob;
	aux->final_phased_prob[aux->n_final_phased_prob++] = final_phased_prob;
}

static int hk_blind_build_softall_map(const struct hk_bmap *bmap,
									 const struct hk_blind_bpair_set *set,
									 struct hk_blind_softall_aux *aux)
{
	struct hk_map *m = 0;
	int32_t *ploidy_XY = 0, *old2new = 0;
	int32_t i, m_pairs = 0;

	assert(bmap);
	assert(bmap->d);
	assert(set);
	assert(aux);
	assert(set->n_raw == 0 || set->raw);
	assert(set->n_raw == 0 || set->raw2binned);
	memset(aux, 0, sizeof(*aux));

	ploidy_XY = hk_sd_ploidy_XY(bmap->d, 0);
	old2new = CALLOC(int32_t, bmap->d->n);
	m = CALLOC(struct hk_map, 1);
	if (ploidy_XY == 0 || old2new == 0 || m == 0)
		goto fail;
	for (i = 1; i < bmap->d->n; ++i)
		old2new[i] = old2new[i - 1] + (ploidy_XY[i - 1] >> 8);
	m->d = hk_sd_split_phase(bmap->d, ploidy_XY);
	if (m->d == 0)
		goto fail;

	for (i = 0; i < set->n_raw; ++i) {
		const struct hk_blind_raw2binned *r2b = &set->raw2binned[i];
		const struct hk_blind_bpair *bp;
		const struct hk_blind_pair *q = &set->raw[i];
		int s;
		assert(r2b->bpair_id >= 0 && r2b->bpair_id < set->n_bpairs);
		bp = &set->bpairs[r2b->bpair_id];
		if (hk_blind_bpair_is_same_bin(bp)) {
			++aux->n_same_bin_skip_raw;
			continue;
		}
		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			float p = bp->p4[s];
			assert(isfinite(p));
			assert(p >= 0.0f);
			if (p <= 0.0f)
				continue;
			hk_blind_softall_map_push_pair(m, &m_pairs, q, ploidy_XY,
									 old2new, s, r2b->swapped, p, p, aux);
			++aux->n_selected_raw;
			++aux->state_count[s];
		}
	}
	hk_pair_sort(m->n_pairs, m->pairs);
	m->cols = 1 << 8;
	free(old2new);
	free(ploidy_XY);
	aux->map = m;
	return 0;

fail:
	free(old2new);
	free(ploidy_XY);
	free(aux->final_phased_prob);
	aux->final_phased_prob = 0;
	aux->n_final_phased_prob = aux->m_final_phased_prob = 0;
	if (m) {
		if (m->d)
			hk_map_destroy(m);
		else {
			free(m->pairs);
			free(m);
		}
	}
	return -1;
}

static int hk_blind_i32_cmp(const void *a_, const void *b_)
{
	const int32_t a = *(const int32_t*)a_;
	const int32_t b = *(const int32_t*)b_;
	return a < b? -1 : a > b? 1 : 0;
}

static int32_t hk_blind_lower_bound_i32(const int32_t *a, int32_t n, int32_t x)
{
	int32_t lo = 0, hi = n;
	while (lo < hi) {
		int32_t mid = lo + ((hi - lo) >> 1);
		if (a[mid] < x) lo = mid + 1;
		else hi = mid;
	}
	return lo;
}

static int32_t hk_blind_upper_bound_i32(const int32_t *a, int32_t n, int32_t x)
{
	int32_t lo = 0, hi = n;
	while (lo < hi) {
		int32_t mid = lo + ((hi - lo) >> 1);
		if (a[mid] <= x) lo = mid + 1;
		else hi = mid;
	}
	return lo;
}

static void hk_blind_fenwick_add(double *bit, int32_t n, int32_t idx, double delta)
{
	for (; idx <= n; idx += idx & -idx)
		bit[idx] += delta;
}

static double hk_blind_fenwick_sum(const double *bit, int32_t idx)
{
	double sum = 0.0;
	for (; idx > 0; idx -= idx & -idx)
		sum += bit[idx];
	return sum;
}

static double hk_blind_fenwick_range_sum(const double *bit, int32_t lo, int32_t hi)
{
	if (hi < lo)
		return 0.0;
	return hk_blind_fenwick_sum(bit, hi) - hk_blind_fenwick_sum(bit, lo - 1);
}

static void hk_blind_fenwick_range_add(double *bit, int32_t n,
								   int32_t lo, int32_t hi, double delta)
{
	if (hi < lo)
		return;
	hk_blind_fenwick_add(bit, n, lo, delta);
	if (hi + 1 <= n)
		hk_blind_fenwick_add(bit, n, hi + 1, -delta);
}

static int hk_blind_pair_count_nei_weighted(int32_t n_pairs,
										const struct hk_pair *pairs,
										int32_t r1,
										int32_t r2,
										double *weights)
{
	int32_t k, st;
	assert(weights);
	for (k = 1; k < n_pairs; ++k)
		if (pairs[k - 1].chr > pairs[k].chr ||
			(pairs[k - 1].chr == pairs[k].chr &&
			 pairs[k - 1].pos > pairs[k].pos))
			break;
	assert(k == n_pairs);
	for (k = 0; k < n_pairs; ++k)
		weights[k] = 0.0;
	for (k = 1, st = 0; k <= n_pairs; ++k) {
		if (k == n_pairs || pairs[k - 1].chr != pairs[k].chr) {
			int32_t n = k - st;
			int32_t *pos2 = 0;
			double *active_bit = 0, *add_bit = 0, *insert_add = 0;
			int32_t i, j, n_pos = 0, left = 0;
			if (n <= 0) {
				st = k;
				continue;
			}
			pos2 = MALLOC(int32_t, n);
			active_bit = CALLOC(double, n + 2);
			add_bit = CALLOC(double, n + 3);
			insert_add = CALLOC(double, n);
			if (pos2 == 0 || active_bit == 0 || add_bit == 0 || insert_add == 0) {
				free(pos2);
				free(active_bit);
				free(add_bit);
				free(insert_add);
				return -1;
			}
			for (i = 0; i < n; ++i)
				pos2[i] = hk_ppos2(&pairs[st + i]);
			qsort(pos2, (size_t)n, sizeof(*pos2), hk_blind_i32_cmp);
			for (i = 0; i < n; ++i)
				if (i == 0 || pos2[i] != pos2[n_pos - 1])
					pos2[n_pos++] = pos2[i];
			for (i = 0; i < n; ++i) {
				const struct hk_pair *p = &pairs[st + i];
				int32_t p1 = hk_ppos1(p);
				int32_t p2 = hk_ppos2(p);
				int32_t idx = hk_blind_lower_bound_i32(pos2, n_pos, p2) + 1;
				int32_t lo, hi;
				double w = hk_blind_clip01(p->_.phased_prob);
				while (left < i &&
					   p1 - hk_ppos1(&pairs[st + left]) >= r1) {
					const struct hk_pair *q = &pairs[st + left];
					int32_t qidx = hk_blind_lower_bound_i32(pos2, n_pos,
												hk_ppos2(q)) + 1;
					double qw = hk_blind_clip01(q->_.phased_prob);
					weights[st + left] +=
						hk_blind_fenwick_sum(add_bit, qidx) - insert_add[left];
					if (qw > 0.0)
						hk_blind_fenwick_add(active_bit, n_pos, qidx, -qw);
					++left;
				}
				lo = hk_blind_upper_bound_i32(pos2, n_pos,
										p2 > r2? p2 - r2 : 0) + 1;
				hi = hk_blind_lower_bound_i32(pos2, n_pos, p2 + r2);
				weights[st + i] += hk_blind_fenwick_range_sum(active_bit, lo, hi);
				if (w > 0.0)
					hk_blind_fenwick_range_add(add_bit, n_pos, lo, hi, w);
				hk_blind_fenwick_add(active_bit, n_pos, idx, w);
				insert_add[i] = hk_blind_fenwick_sum(add_bit, idx);
			}
			for (j = left; j < n; ++j) {
				const struct hk_pair *q = &pairs[st + j];
				int32_t qidx = hk_blind_lower_bound_i32(pos2, n_pos,
											hk_ppos2(q)) + 1;
				weights[st + j] +=
					hk_blind_fenwick_sum(add_bit, qidx) - insert_add[j];
			}
			free(pos2);
			free(active_bit);
			free(add_bit);
			free(insert_add);
			st = k;
		}
	}
	return 0;
}

static int32_t hk_blind_pair_filter_isolated_weighted(int32_t n_pairs,
										  struct hk_pair *pairs,
										  int32_t max_radius,
										  double min_weight,
										  int *failed)
{
	struct hk_blind_pair_weight_aux *weights = 0;
	double *raw_weights = 0;
	int32_t i, k;
	assert(failed);
	*failed = 0;
	if (n_pairs <= 0)
		return 0;
	weights = CALLOC(struct hk_blind_pair_weight_aux, n_pairs);
	raw_weights = CALLOC(double, n_pairs);
	if (weights == 0 || raw_weights == 0) {
		free(weights);
		free(raw_weights);
		*failed = 1;
		return 0;
	}
	if (hk_blind_pair_count_nei_weighted(n_pairs, pairs, max_radius,
									   max_radius, raw_weights) != 0) {
		free(weights);
		free(raw_weights);
		*failed = 1;
		return 0;
	}
	for (i = 0; i < n_pairs; ++i) {
		weights[i].i = i;
		weights[i].weight = raw_weights[i];
	}
	free(raw_weights);
	for (i = k = 0; i < n_pairs; ++i) {
		if (weights[i].weight >= min_weight)
			pairs[k++] = pairs[i];
	}
	free(weights);
	return k;
}

static int hk_blind_pair_count_nei_weighted_to_int(int32_t n_pairs,
										   struct hk_pair *pairs,
										   int32_t r1,
										   int32_t r2)
{
	double *weights = 0;
	int32_t i;
	if (n_pairs <= 0)
		return 0;
	weights = CALLOC(double, n_pairs);
	if (weights == 0)
		return -1;
	if (hk_blind_pair_count_nei_weighted(n_pairs, pairs, r1, r2, weights) != 0) {
		free(weights);
		return -1;
	}
	for (i = 0; i < n_pairs; ++i) {
		double w = weights[i];
		if (w < 0.0)
			w = 0.0;
		if (w > (double)INT32_MAX)
			w = (double)INT32_MAX;
		pairs[i].n_nei = (uint32_t)(w + 0.5);
		pairs[i].n_nei_corner = pairs[i].n_nei;
	}
	free(weights);
	return 0;
}

static int hk_blind_softall_apply_final_prob_to_filtered_pairs(struct hk_map *m,
											   const float *final_prob,
											   int32_t n_final_prob)
{
	int32_t i;
	assert(m);
	assert(final_prob || n_final_prob == 0);
	assert(n_final_prob >= 0);
	for (i = 0; i < m->n_pairs; ++i) {
		uint32_t id = m->pairs[i].n_ctn;
		if (id >= (uint32_t)n_final_prob)
			return -1;
		m->pairs[i]._.phased_prob = final_prob[id];
	}
	return 0;
}

static int hk_blind_wedge_list_build_softall_from_map(struct hk_blind_wedge_list *out,
											  const struct hk_bmap *bmap,
											  struct hk_blind_softall_aux *aux,
											  int d_scale_mode, float d_scale_eps_count)
{
	struct hk_bmap *split_bmap = 0;
	int32_t *split_to_diploid = 0;
	int32_t i;
	int32_t median_nei;
	int ret = -1;

	assert(out);
	assert(bmap);
	assert(aux);
	assert(aux->map);
	assert(aux->n_final_phased_prob == aux->map->n_pairs);
	assert(hk_blind_d_scale_mode_valid(d_scale_mode));
	assert(isfinite(d_scale_eps_count));
	assert(d_scale_eps_count > 0.0f);

	out->n_softall_selected_raw = aux->n_selected_raw;
	out->n_softall_gate_skip_raw = aux->n_gate_skip_raw;
	out->n_softall_same_bin_skip_raw = aux->n_same_bin_skip_raw;
	memcpy(out->n_softall_state_raw_count, aux->state_count,
		   sizeof(out->n_softall_state_raw_count));
	if (aux->map->n_pairs == 0)
		return 0;
	out->n_softall_candidate_pairs = aux->map->n_pairs;
	{
		int filter_failed = 0;
		aux->map->n_pairs = hk_blind_pair_filter_isolated_weighted(
			aux->map->n_pairs, aux->map->pairs, 1000000, 1.0, &filter_failed);
		if (filter_failed)
			goto cleanup;
	}
	out->n_softall_filter1_pairs = aux->map->n_pairs;
	aux->map->cols |= 1 << 6 | 1 << 9;
	{
		int filter_failed = 0;
		aux->map->n_pairs = hk_blind_pair_filter_isolated_weighted(
			aux->map->n_pairs, aux->map->pairs, 10000000, 5.0, &filter_failed);
		if (filter_failed)
			goto cleanup;
	}
	out->n_softall_filter2_pairs = aux->map->n_pairs;
	aux->map->cols |= 1 << 6 | 1 << 9;
	if (aux->map->n_pairs == 0)
		return 0;
	hk_pair_sort(aux->map->n_pairs, aux->map->pairs);
	if (hk_blind_pair_count_nei_weighted_to_int(aux->map->n_pairs,
										  aux->map->pairs,
										  10000000, 10000000) != 0)
		goto cleanup;
	if (hk_blind_softall_apply_final_prob_to_filtered_pairs(
			aux->map, aux->final_phased_prob, aux->n_final_phased_prob) != 0)
		goto cleanup;
	split_bmap = hk_bmap_gen(aux->map->d, aux->map->n_pairs,
						 aux->map->pairs, 1000000, 1);
	if (split_bmap == 0)
		goto cleanup;
	out->n_softall_bmap_pairs = split_bmap->n_pairs;
	split_to_diploid = CALLOC(int32_t, split_bmap->n_beads);
	if (split_bmap->n_beads > 0 && split_to_diploid == 0)
		goto cleanup;
	if (hk_blind_build_split_lookup(split_bmap, bmap, split_to_diploid) != 0)
		goto cleanup;
	median_nei = hk_blind_split_graph_median_nei(split_bmap);
	out->n_input_bpair = split_bmap->n_pairs;
	for (i = 0; i < split_bmap->n_pairs; ++i) {
		const struct hk_bpair *p = &split_bmap->pairs[i];
		struct hk_blind_wedge edge;
		int32_t d0, d1;
		float k, d_scale, prob, effective_n;
		if (p->bid[0] == p->bid[1]) {
			++out->n_skipped_self_edges;
			continue;
		}
		d0 = split_to_diploid[p->bid[0]];
		d1 = split_to_diploid[p->bid[1]];
		if (d0 == d1) {
			++out->n_skipped_self_edges;
			continue;
		}
		prob = hk_blind_clip01(p->p);
		k = hk_blind_native_fdg_k_from_nei(p, median_nei) * prob;
		if (k <= 0.0f)
			continue;
		effective_n = (float)p->n;
		if (d_scale_mode == HK_BLIND_D_SCALE_EXPECTED_COUNT)
			effective_n *= prob;
		if (!isfinite(effective_n) || effective_n < d_scale_eps_count)
			effective_n = d_scale_eps_count;
		d_scale = powf(effective_n, -1.0f / 3.0f);
		if (d1 < d0) {
			int32_t t = d0;
			d0 = d1;
			d1 = t;
		}
		hk_blind_wedge_set(&edge, d0, d1, k, d_scale, HK_BLIND_STATE_00);
		edge.state_mask = 0xffu;
		hk_blind_wedge_list_push(out, &edge);
	}
	out->n_expanded_edges = out->n_edges;
	out->n_edges_before_aggregation = out->n_edges;
	if (out->n_edges > 0) {
		out->mean_rho_train_bpair = 1.0f;
		out->min_rho_train_bpair = 1.0f;
		out->max_rho_train_bpair = 1.0f;
		out->sum_rho_train_bpair = (double)out->n_edges;
	}
	ret = 0;

cleanup:
	free(split_to_diploid);
	if (split_bmap)
		hk_bmap_destroy(split_bmap);
	return ret;
}

int hk_blind_wedge_list_build_softall(struct hk_blind_wedge_list *out,
									  const struct hk_bmap *bmap,
									  const struct hk_blind_bpair_set *set)
{
	return hk_blind_wedge_list_build_softall_mode(out, bmap, set,
												 HK_BLIND_D_SCALE_EXPECTED_COUNT,
												 1e-6f);
}

int hk_blind_wedge_list_build_softall_mode(struct hk_blind_wedge_list *out,
										   const struct hk_bmap *bmap,
										   const struct hk_blind_bpair_set *set,
										   int d_scale_mode,
										   float d_scale_eps_count)
{
	struct hk_blind_softall_aux aux;
	int ret = -1;
	assert(out);
	assert(bmap);
	assert(set);
	assert(hk_blind_d_scale_mode_valid(d_scale_mode));
	assert(isfinite(d_scale_eps_count));
	assert(d_scale_eps_count > 0.0f);
	memset(&aux, 0, sizeof(aux));

	out->n_edges = 0;
	out->n_input_bpair = set->n_raw;
	out->n_expanded_edges = 0;
	out->n_softall_candidate_pairs = 0;
	out->n_softall_filter1_pairs = 0;
	out->n_softall_filter2_pairs = 0;
	out->n_softall_bmap_pairs = 0;
	out->n_softall_selected_raw = 0;
	out->n_softall_gate_skip_raw = 0;
	out->n_softall_same_bin_skip_raw = 0;
	memset(out->n_softall_state_raw_count, 0, sizeof(out->n_softall_state_raw_count));
	out->n_skipped_self_edges = 0;
	out->n_skipped_same_bin_bpairs = 0;
	out->n_edges_before_aggregation = 0;
	out->n_aggregated_edges_removed = 0;
	out->sum_rho_train_bpair = 0.0;
	out->mean_rho_train_bpair = 0.0f;
	out->min_rho_train_bpair = 0.0f;
	out->max_rho_train_bpair = 0.0f;

	if (hk_blind_build_softall_map(bmap, set, &aux) != 0)
		goto cleanup;
	ret = hk_blind_wedge_list_build_softall_from_map(out, bmap, &aux,
													 d_scale_mode, d_scale_eps_count);

cleanup:
	free(aux.final_phased_prob);
	if (aux.map)
		hk_map_destroy(aux.map);
	return ret;
}

static int hk_blind_wedge_cmp(const void *a_, const void *b_)
{
	const struct hk_blind_wedge *a = (const struct hk_blind_wedge*)a_;
	const struct hk_blind_wedge *b = (const struct hk_blind_wedge*)b_;
	if (a->bid[0] != b->bid[0])
		return a->bid[0] < b->bid[0]? -1 : 1;
	if (a->bid[1] != b->bid[1])
		return a->bid[1] < b->bid[1]? -1 : 1;
	if (a->d_scale != b->d_scale)
		return a->d_scale < b->d_scale? -1 : 1;
	if (a->state != b->state)
		return a->state < b->state? -1 : 1;
	if (a->state_mask != b->state_mask)
		return a->state_mask < b->state_mask? -1 : 1;
	if (a->k != b->k)
		return a->k < b->k? -1 : 1;
	return 0;
}

static int hk_blind_wedge_same_key(const struct hk_blind_wedge *a, const struct hk_blind_wedge *b)
{
	return a->bid[0] == b->bid[0] && a->bid[1] == b->bid[1] && a->d_scale == b->d_scale;
}

int hk_blind_wedge_list_aggregate_exact(struct hk_blind_wedge_list *list)
{
	int32_t i, dst, n0, n_pos;

	assert(list);
	assert(list->n_edges >= 0);
	assert(list->n_edges == 0 || list->edges);

	n0 = list->n_edges;
	list->n_edges_before_aggregation = n0;
	list->n_aggregated_edges_removed = 0;
	if (n0 == 0) return 0;

	for (i = 0; i < n0; ++i) {
		assert(list->edges[i].bid[0] >= 0);
		assert(list->edges[i].bid[1] >= 0);
		assert(list->edges[i].bid[0] <= list->edges[i].bid[1]);
		assert(isfinite(list->edges[i].k));
		assert(list->edges[i].k >= 0.0f);
		assert(isfinite(list->edges[i].d_scale));
		assert(list->edges[i].d_scale > 0.0f);
		assert(list->edges[i].state >= 0 && list->edges[i].state < HK_BLIND_N_STATE);
		assert(list->edges[i].state_mask != 0);
	}
	n_pos = 0;
	for (i = 0; i < n0; ++i) {
		if (list->edges[i].k <= 0.0f)
			continue;
		if (n_pos != i)
			list->edges[n_pos] = list->edges[i];
		++n_pos;
	}
	if (n_pos == 0) {
		list->n_edges = 0;
		list->n_aggregated_edges_removed = n0;
		return 0;
	}
	n0 = n_pos;

	qsort(list->edges, n0, sizeof(*list->edges), hk_blind_wedge_cmp);
	dst = 0;
	for (i = 1; i < n0; ++i) {
		if (hk_blind_wedge_same_key(&list->edges[dst], &list->edges[i])) {
			list->edges[dst].k += list->edges[i].k;
			list->edges[dst].state_mask |= list->edges[i].state_mask;
		} else {
			list->edges[++dst] = list->edges[i];
		}
	}
	list->n_edges = dst + 1;
	list->n_aggregated_edges_removed = list->n_edges_before_aggregation - list->n_edges;
	return 0;
}

void hk_blind_iter_diag_init(struct hk_blind_iter_diag *diag)
{
	assert(diag);
	diag->n_bpair = 0;
	diag->n_raw = 0;
	diag->mean_entropy = 0.0f;
	diag->mean_pmax = 0.0f;
	diag->mean_margin = 0.0f;
	diag->mean_pU = 0.0f;
	diag->mean_rho_output = 0.0f;
	diag->n_posterior_nonfinite = 0;
	diag->n_posterior_bad_sum = 0;
	diag->n_posterior_out_of_range = 0;
	diag->n_uncertainty_nonfinite = 0;
	diag->n_uncertainty_out_of_range = 0;
	diag->n_five_state_bad_sum = 0;
	hk_blind_homolog_sep_stats_init(&diag->sep_stats);
	diag->sep_energy = 0.0f;
	diag->sep_force_l1 = 0.0f;
	diag->sep_force_nonfinite = 0;
	diag->n_wedges = 0;
	diag->n_wedges_before_aggregation = 0;
	diag->n_expanded_edges = 0;
	diag->n_softall_candidate_pairs = 0;
	diag->n_softall_filter1_pairs = 0;
	diag->n_softall_filter2_pairs = 0;
	diag->n_softall_bmap_pairs = 0;
	diag->n_softall_selected_raw = 0;
	diag->n_softall_gate_skip_raw = 0;
	diag->n_softall_same_bin_skip_raw = 0;
	memset(diag->n_softall_state_raw_count, 0, sizeof(diag->n_softall_state_raw_count));
	diag->n_skipped_self_edges = 0;
	diag->n_skipped_same_bin_bpairs = 0;
	diag->sum_wedge_k = 0.0;
	diag->mean_rho_train_bpair = 0.0f;
	diag->min_rho_train_bpair = 0.0f;
	diag->max_rho_train_bpair = 0.0f;
	diag->n_wedge_nonfinite = 0;
	diag->n_wedge_bad_k = 0;
	diag->n_wedge_bad_d_scale = 0;
	hk_blind_gauge_stats_init(&diag->gauge_stats);
	diag->n_chr_flipped = 0;
}

static void hk_blind_iter_diag_check_bpair_set_shape(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set)
{
	int32_t i;

	assert(bmap);
	assert(set);
	assert(bmap->n_beads >= 0);
	assert(set->n_bpairs >= 0);
	assert(set->n_bpairs == 0 || set->bpairs);
	assert(set->n_raw >= 0);
	for (i = 0; i < set->n_bpairs; ++i) {
		assert(set->bpairs[i].key.bid[0] >= 0 && set->bpairs[i].key.bid[0] < bmap->n_beads);
		assert(set->bpairs[i].key.bid[1] >= 0 && set->bpairs[i].key.bid[1] < bmap->n_beads);
	}
}

void hk_blind_iter_diag_validate_bpair_set(const struct hk_blind_bpair_set *set, struct hk_blind_iter_diag *diag)
{
	const float tol = 1e-4f;
	const float log_n_state = logf((float)HK_BLIND_N_STATE);
	const float min_pmax = 1.0f / (float)HK_BLIND_N_STATE;
	int32_t i;
	int s;

	assert(set);
	assert(diag);
	assert(set->n_bpairs >= 0);
	assert(set->n_bpairs == 0 || set->bpairs);

	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *p = &set->bpairs[i];
		double sum_p4 = 0.0;
		int posterior_nonfinite = 0;
		int posterior_out_of_range = 0;
		int uncertainty_nonfinite = 0;
		int uncertainty_out_of_range = 0;

		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			float ps = p->p4[s];
			if (!hk_blind_float_isfinite(ps)) {
				posterior_nonfinite = 1;
			} else {
				sum_p4 += ps;
				if (ps < -tol || ps > 1.0f + tol)
					posterior_out_of_range = 1;
			}
		}
		if (posterior_nonfinite) {
			++diag->n_posterior_nonfinite;
		} else if (fabs(sum_p4 - 1.0) > tol) {
			++diag->n_posterior_bad_sum;
		}
		if (posterior_out_of_range)
			++diag->n_posterior_out_of_range;

		if (!hk_blind_float_isfinite(p->entropy) || !hk_blind_float_isfinite(p->pmax) ||
			!hk_blind_float_isfinite(p->margin) || !hk_blind_float_isfinite(p->rho_output) ||
			!hk_blind_float_isfinite(p->pU)) {
			uncertainty_nonfinite = 1;
		} else {
			if (p->entropy < -tol || p->entropy > log_n_state + tol)
				uncertainty_out_of_range = 1;
			if (p->pmax < min_pmax - tol || p->pmax > 1.0f + tol)
				uncertainty_out_of_range = 1;
			if (p->margin < -tol || p->margin > 1.0f + tol)
				uncertainty_out_of_range = 1;
			if (p->rho_output < -tol || p->rho_output > 1.0f + tol)
				uncertainty_out_of_range = 1;
			if (p->pU < -tol || p->pU > 1.0f + tol)
				uncertainty_out_of_range = 1;
		}
		if (uncertainty_nonfinite)
			++diag->n_uncertainty_nonfinite;
		else if (uncertainty_out_of_range)
			++diag->n_uncertainty_out_of_range;

		if (!posterior_nonfinite && !uncertainty_nonfinite) {
			double five_sum = (double)p->rho_output * sum_p4 + p->pU;
			if (fabs(five_sum - 1.0) > tol)
				++diag->n_five_state_bad_sum;
		}
	}
}

static void hk_blind_iter_diag_update_posterior_stats(const struct hk_blind_bpair_set *set, struct hk_blind_iter_diag *diag)
{
	double entropy = 0.0, pmax = 0.0, margin = 0.0, pU = 0.0, rho = 0.0;
	int32_t i;

	assert(set);
	assert(diag);
	if (set->n_bpairs == 0)
		return;
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *p = &set->bpairs[i];
		assert(hk_blind_float_isfinite(p->entropy));
		assert(hk_blind_float_isfinite(p->pmax));
		assert(hk_blind_float_isfinite(p->margin));
		assert(hk_blind_float_isfinite(p->pU));
		assert(hk_blind_float_isfinite(p->rho_output));
		entropy += p->entropy;
		pmax += p->pmax;
		margin += p->margin;
		pU += p->pU;
		rho += p->rho_output;
	}
	diag->mean_entropy = (float)(entropy / set->n_bpairs);
	diag->mean_pmax = (float)(pmax / set->n_bpairs);
	diag->mean_margin = (float)(margin / set->n_bpairs);
	diag->mean_pU = (float)(pU / set->n_bpairs);
	diag->mean_rho_output = (float)(rho / set->n_bpairs);
}

static void hk_blind_iter_diag_update_force_stats(const fvec3_t *force, int32_t n_diploid, struct hk_blind_iter_diag *diag)
{
	double l1 = 0.0;
	int32_t i;
	int a;

	assert(diag);
	if (n_diploid == 0)
		return;
	assert(force);
	for (i = 0; i < n_diploid; ++i) {
		for (a = 0; a < 3; ++a) {
			if (hk_blind_float_isfinite(force[i][a]))
				l1 += fabs((double)force[i][a]);
			else
				++diag->sep_force_nonfinite;
		}
	}
	diag->sep_force_l1 = (float)l1;
}

void hk_blind_iter_diag_validate_wedge_list(const struct hk_blind_wedge_list *list, struct hk_blind_iter_diag *diag)
{
	double sum_k = 0.0;
	int32_t i;

	assert(list);
	assert(diag);
	for (i = 0; i < list->n_edges; ++i) {
		const struct hk_blind_wedge *e = &list->edges[i];
		int k_finite = hk_blind_float_isfinite(e->k);
		int d_scale_finite = hk_blind_float_isfinite(e->d_scale);
		if (!k_finite || !d_scale_finite)
			++diag->n_wedge_nonfinite;
		if (!k_finite || e->k < 0.0f)
			++diag->n_wedge_bad_k;
		if (!d_scale_finite || e->d_scale <= 0.0f)
			++diag->n_wedge_bad_d_scale;
		if (k_finite)
			sum_k += e->k;
	}
	diag->sum_wedge_k = sum_k;
}

static int hk_blind_iter_diag_status_ok(const struct hk_blind_iter_diag *diag)
{
	assert(diag);
	return diag->n_posterior_nonfinite == 0 &&
		   diag->n_posterior_bad_sum == 0 &&
		   diag->n_posterior_out_of_range == 0 &&
		   diag->n_uncertainty_nonfinite == 0 &&
		   diag->n_uncertainty_out_of_range == 0 &&
		   diag->n_five_state_bad_sum == 0 &&
		   diag->n_wedge_nonfinite == 0 &&
		   diag->n_wedge_bad_k == 0 &&
		   diag->n_wedge_bad_d_scale == 0 &&
		   diag->sep_force_nonfinite == 0;
}

int hk_blind_iter_diag_snprintf(char *buf, size_t buf_size, const struct hk_blind_iter_diag *diag)
{
	const char *status;

	assert(diag);
	assert(buf_size == 0 || buf);

	status = hk_blind_iter_diag_status_ok(diag)? "OK" : "WARN";
	return snprintf(buf, buf_size,
					"status: %s\n"
					"counts: n_raw=%d n_bpair=%d\n"
					"posterior: mean_entropy=%.9g mean_pmax=%.9g mean_margin=%.9g mean_pU=%.9g mean_rho_output=%.9g\n"
					"bad: n_posterior_nonfinite=%d n_posterior_bad_sum=%d n_posterior_out_of_range=%d "
					"n_uncertainty_nonfinite=%d n_uncertainty_out_of_range=%d n_five_state_bad_sum=%d "
					"sep_force_nonfinite=%d n_wedge_nonfinite=%d n_wedge_bad_k=%d n_wedge_bad_d_scale=%d\n"
					"sep: sep_stats.n_haploid=%d sep_stats.n_finite=%d sep_stats.n_nonfinite=%d "
					"sep_stats.n_collapsed=%d sep_stats.min_sep=%.9g sep_stats.mean_sep=%.9g "
					"sep_stats.max_sep=%.9g sep_stats.collapse_threshold=%.9g\n"
					"anti_collapse: sep_energy=%.9g sep_force_l1=%.9g sep_force_nonfinite=%d\n"
					"wedges: n_wedges=%d n_wedges_before_aggregation=%d n_skipped_self_edges=%lld "
					"n_skipped_same_bin_bpairs=%lld sum_wedge_k=%.17g "
					"mean_rho_train_bpair=%.9g min_rho_train_bpair=%.9g max_rho_train_bpair=%.9g "
					"n_wedge_nonfinite=%d n_wedge_bad_k=%d n_wedge_bad_d_scale=%d\n"
					"gauge: gauge_stats.n_chr=%d gauge_stats.n_flipped=%d "
					"gauge_stats.mismatch_no_flip=%.17g gauge_stats.mismatch_flip=%.17g "
					"gauge_stats.mismatch_chosen=%.17g n_chr_flipped=%d\n",
					status,
					diag->n_raw, diag->n_bpair,
					diag->mean_entropy, diag->mean_pmax, diag->mean_margin, diag->mean_pU, diag->mean_rho_output,
					diag->n_posterior_nonfinite, diag->n_posterior_bad_sum, diag->n_posterior_out_of_range,
					diag->n_uncertainty_nonfinite, diag->n_uncertainty_out_of_range, diag->n_five_state_bad_sum,
					diag->sep_force_nonfinite, diag->n_wedge_nonfinite, diag->n_wedge_bad_k, diag->n_wedge_bad_d_scale,
					diag->sep_stats.n_haploid, diag->sep_stats.n_finite, diag->sep_stats.n_nonfinite,
					diag->sep_stats.n_collapsed, diag->sep_stats.min_sep, diag->sep_stats.mean_sep,
					diag->sep_stats.max_sep, diag->sep_stats.collapse_threshold,
					diag->sep_energy, diag->sep_force_l1, diag->sep_force_nonfinite,
					diag->n_wedges, diag->n_wedges_before_aggregation, (long long)diag->n_skipped_self_edges,
					(long long)diag->n_skipped_same_bin_bpairs, diag->sum_wedge_k,
					diag->mean_rho_train_bpair, diag->min_rho_train_bpair, diag->max_rho_train_bpair,
					diag->n_wedge_nonfinite, diag->n_wedge_bad_k, diag->n_wedge_bad_d_scale,
					diag->gauge_stats.n_chr, diag->gauge_stats.n_flipped,
					diag->gauge_stats.mismatch_no_flip, diag->gauge_stats.mismatch_flip,
					diag->gauge_stats.mismatch_chosen, diag->n_chr_flipped);
}

int hk_blind_run_single_iter_diag(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
								  const struct hk_fdg_conf *conf, const fvec3_t *prev_coords,
								  fvec3_t *cur_coords, float unit, float d_scale, float base_k,
								  const float log_prior[HK_BLIND_N_STATE], float temperature,
								  float rho_train, float min_sep_unit, float lambda_sep,
								  struct hk_blind_iter_diag *diag)
{
	struct hk_blind_wedge_list wedges;
	float local_log_prior[HK_BLIND_N_STATE];
	fvec3_t *sep_force = 0;
	int32_t n_diploid;
	int ret;

	assert(bmap);
	assert(set);
	assert(conf);
	assert(cur_coords);
	assert(unit > 0.0f);
	assert(d_scale > 0.0f);
	assert(base_k >= 0.0f);
	assert(temperature > 0.0f);
	assert(rho_train >= 0.0f);
	assert(min_sep_unit >= 0.0f);
	assert(lambda_sep >= 0.0f);
	assert(diag);
	hk_blind_iter_diag_check_bpair_set_shape(bmap, set);

	hk_blind_iter_diag_init(diag);
	diag->n_bpair = set->n_bpairs;
	diag->n_raw = set->n_raw;
	diag->gauge_stats.n_chr = hk_blind_bmap_n_chr(bmap);

	if (log_prior == 0) {
		hk_blind_init_uniform_log_prior(local_log_prior);
		log_prior = local_log_prior;
	}

	if (prev_coords) {
		ret = hk_blind_temporal_gauge_stabilize_bmap(bmap, prev_coords, cur_coords, 0, &diag->gauge_stats);
		if (ret != 0) return ret;
		diag->n_chr_flipped = diag->gauge_stats.n_flipped;
	}

	// Step 10A computes a fresh coordinate posterior after coordinate gauge stabilization.
	// This no-FDG diagnostic driver intentionally reports caller-scalar posterior params.
	// hk_blind_run_single_iter_cpu() uses per-bpair params for core E-step scoring.
	hk_blind_bpair_set_update_posterior_from_coords(set, conf, cur_coords, unit, d_scale, base_k, log_prior, temperature);
	hk_blind_iter_diag_update_posterior_stats(set, diag);
	hk_blind_iter_diag_validate_bpair_set(set, diag);

	ret = hk_blind_homolog_sep_compute_stats(bmap->n_beads, cur_coords, min_sep_unit * unit, &diag->sep_stats);
	if (ret != 0) return ret;

	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	if (n_diploid > 0)
		sep_force = CALLOC(fvec3_t, n_diploid);
	diag->sep_energy = hk_blind_homolog_sep_accumulate_force(bmap->n_beads, cur_coords, sep_force, unit, min_sep_unit, lambda_sep);
	hk_blind_iter_diag_update_force_stats(sep_force, n_diploid, diag);
	free(sep_force);

	hk_blind_wedge_list_init(&wedges);
	// This no-FDG diagnostic driver intentionally reports caller-scalar wedge params.
	// hk_blind_run_single_iter_cpu() uses per-bpair params for M-step relaxation.
	ret = hk_blind_wedge_list_build_from_bpair_set(&wedges, set, base_k, d_scale, rho_train);
	if (ret != 0) {
		hk_blind_wedge_list_destroy(&wedges);
		return ret;
	}
	diag->n_wedges_before_aggregation = wedges.n_edges_before_aggregation;
	diag->n_skipped_self_edges = wedges.n_skipped_self_edges;
	diag->n_skipped_same_bin_bpairs = wedges.n_skipped_same_bin_bpairs;
	diag->mean_rho_train_bpair = wedges.mean_rho_train_bpair;
	diag->min_rho_train_bpair = wedges.min_rho_train_bpair;
	diag->max_rho_train_bpair = wedges.max_rho_train_bpair;
	ret = hk_blind_wedge_list_aggregate_exact(&wedges);
	if (ret != 0) {
		hk_blind_wedge_list_destroy(&wedges);
		return ret;
	}
	diag->n_wedges = wedges.n_edges;
	hk_blind_iter_diag_validate_wedge_list(&wedges, diag);
	hk_blind_wedge_list_destroy(&wedges);
	return 0;
}

void hk_blind_single_iter_diag_init(struct hk_blind_single_iter_diag *diag)
{
	assert(diag);
	hk_blind_iter_diag_init(&diag->pre_relax_diag);
	hk_blind_relax_diag_init(&diag->relax_diag);
	hk_blind_gauge_stats_init(&diag->gauge_stats);
	diag->n_chr_flipped = 0;
	diag->n_wedges = 0;
	diag->n_wedges_before_aggregation = 0;
	diag->n_expanded_edges = 0;
	diag->n_softall_candidate_pairs = 0;
	diag->n_softall_filter1_pairs = 0;
	diag->n_softall_filter2_pairs = 0;
	diag->n_softall_bmap_pairs = 0;
	diag->n_softall_selected_raw = 0;
	diag->n_softall_gate_skip_raw = 0;
	diag->n_softall_same_bin_skip_raw = 0;
	memset(diag->n_softall_state_raw_count, 0, sizeof(diag->n_softall_state_raw_count));
	diag->n_skipped_self_edges = 0;
	diag->n_skipped_same_bin_bpairs = 0;
	diag->sum_wedge_k = 0.0;
	diag->mean_rho_train_bpair = 0.0f;
	diag->min_rho_train_bpair = 0.0f;
	diag->max_rho_train_bpair = 0.0f;
}

static int hk_blind_run_single_iter_cpu_impl(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
											 const struct hk_fdg_conf *fdg_conf, fvec3_t *coords,
											 const float log_prior[HK_BLIND_N_STATE],
											 const struct hk_blind_single_iter_conf *iter_conf,
											 const struct hk_blind_coarse_to_fine_map *anchor_map,
											 const fvec3_t *coarse_diploid_coords, float anchor_k,
											 struct hk_blind_single_iter_diag *diag)
{
	struct hk_blind_wedge_list wedges;
	struct hk_blind_iter_diag *pre_diag;
	fvec3_t *prev_coords = 0, *sep_force = 0;
	int32_t n_diploid;
	int ret;

	assert(bmap);
	assert(set);
	assert(fdg_conf);
	assert(coords);
	assert(iter_conf);
	assert(diag);
	assert(bmap->n_beads >= 0);
	assert(bmap->n_beads <= INT32_MAX / HK_DIPLOID_N_COPY);
	assert(isfinite(iter_conf->unit));
	assert(iter_conf->unit > 0.0f);
	assert(isfinite(iter_conf->d_scale));
	assert(iter_conf->d_scale > 0.0f);
	assert(isfinite(iter_conf->base_k));
	assert(iter_conf->base_k >= 0.0f);
	assert(isfinite(iter_conf->temperature));
	assert(iter_conf->temperature > 0.0f);
	assert(isfinite(iter_conf->rho_train));
	assert(iter_conf->rho_train >= 0.0f);
	assert(iter_conf->rho_train_floor == 0.0f ||
		   (isfinite(iter_conf->rho_train_floor) &&
			iter_conf->rho_train_floor >= 0.0f &&
			iter_conf->rho_train_floor <= 1.0f));
	assert(isfinite(iter_conf->min_sep_unit));
	assert(iter_conf->min_sep_unit >= 0.0f);
	assert(isfinite(iter_conf->lambda_sep));
	assert(iter_conf->lambda_sep >= 0.0f);
	assert(isfinite(iter_conf->chr_sep_unit));
	assert(iter_conf->chr_sep_unit >= 0.0f);
	assert(isfinite(iter_conf->lambda_chr_sep));
	assert(iter_conf->lambda_chr_sep >= 0.0f);
	assert(isfinite(iter_conf->relax_step));
	assert(iter_conf->relax_step >= 0.0f);
	assert(iter_conf->relax_steps >= 0);
	assert(iter_conf->enable_repulsion == 0 || iter_conf->enable_repulsion == 1);
	assert(hk_blind_repulsion_mode_valid(iter_conf->repulsion_mode));
	assert(!iter_conf->enable_repulsion ||
		   iter_conf->repulsion_mode == HK_BLIND_REPULSION_N2 ||
		   iter_conf->repulsion_mode == HK_BLIND_REPULSION_CELL);
	assert(hk_blind_rho_train_mode_valid(iter_conf->rho_train_mode));
	assert(hk_blind_d_scale_mode_valid(iter_conf->d_scale_mode));
	assert(isfinite(iter_conf->d_scale_eps_count));
	assert(iter_conf->d_scale_eps_count > 0.0f);
	assert(hk_blind_estep_score_mode_valid(iter_conf->estep_score_mode));
	assert(isfinite(anchor_k));
	assert(anchor_k >= 0.0f);
	if (anchor_map) {
		assert(anchor_map->n_fine == bmap->n_beads);
		assert(coarse_diploid_coords || anchor_map->n_coarse == 0);
	} else {
		assert(anchor_k == 0.0f);
	}
	hk_blind_iter_diag_check_bpair_set_shape(bmap, set);

	hk_blind_single_iter_diag_init(diag);
	pre_diag = &diag->pre_relax_diag;
	pre_diag->n_bpair = set->n_bpairs;
	pre_diag->n_raw = set->n_raw;
	pre_diag->gauge_stats.n_chr = hk_blind_bmap_n_chr(bmap);

	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	if (n_diploid > 0) {
		prev_coords = MALLOC(fvec3_t, n_diploid);
		memcpy(prev_coords, coords, (size_t)n_diploid * sizeof(*coords));
	}

	// Core CPU blind iteration uses per-bpair Hickit-style contact params for E-step scoring.
	// iter_conf->d_scale/base_k are retained for legacy/global diagnostics and ablations.
	hk_blind_bpair_set_update_posterior_from_coords_params_score_mode(set, fdg_conf, coords,
																	  iter_conf->unit,
																	  log_prior,
																	  iter_conf->temperature,
																	  iter_conf->estep_score_mode);
	hk_blind_iter_diag_update_posterior_stats(set, pre_diag);
	hk_blind_iter_diag_validate_bpair_set(set, pre_diag);

	ret = hk_blind_homolog_sep_compute_stats(bmap->n_beads, coords,
											 iter_conf->min_sep_unit * iter_conf->unit,
											 &pre_diag->sep_stats);
	if (ret != 0) {
		free(prev_coords);
		return ret;
	}
	if (n_diploid > 0)
		sep_force = CALLOC(fvec3_t, n_diploid);
	pre_diag->sep_energy = hk_blind_homolog_sep_accumulate_force(bmap->n_beads, coords, sep_force,
																 iter_conf->unit, iter_conf->min_sep_unit,
																 iter_conf->lambda_sep);
	hk_blind_iter_diag_update_force_stats(sep_force, n_diploid, pre_diag);
	free(sep_force);

	hk_blind_wedge_list_init(&wedges);
	ret = hk_blind_wedge_list_build_softall_mode(&wedges, bmap, set,
												 iter_conf->d_scale_mode,
												 iter_conf->d_scale_eps_count);
	if (ret != 0) {
		hk_blind_wedge_list_destroy(&wedges);
		free(prev_coords);
		return ret;
	}
	pre_diag->n_wedges_before_aggregation = wedges.n_edges_before_aggregation;
	pre_diag->n_expanded_edges = wedges.n_expanded_edges;
	pre_diag->n_softall_candidate_pairs = wedges.n_softall_candidate_pairs;
	pre_diag->n_softall_filter1_pairs = wedges.n_softall_filter1_pairs;
	pre_diag->n_softall_filter2_pairs = wedges.n_softall_filter2_pairs;
	pre_diag->n_softall_bmap_pairs = wedges.n_softall_bmap_pairs;
	pre_diag->n_softall_selected_raw = wedges.n_softall_selected_raw;
	pre_diag->n_softall_gate_skip_raw = wedges.n_softall_gate_skip_raw;
	pre_diag->n_softall_same_bin_skip_raw = wedges.n_softall_same_bin_skip_raw;
	memcpy(pre_diag->n_softall_state_raw_count, wedges.n_softall_state_raw_count,
		   sizeof(pre_diag->n_softall_state_raw_count));
	pre_diag->n_skipped_self_edges = wedges.n_skipped_self_edges;
	pre_diag->n_skipped_same_bin_bpairs = wedges.n_skipped_same_bin_bpairs;
	pre_diag->mean_rho_train_bpair = wedges.mean_rho_train_bpair;
	pre_diag->min_rho_train_bpair = wedges.min_rho_train_bpair;
	pre_diag->max_rho_train_bpair = wedges.max_rho_train_bpair;
	diag->n_wedges_before_aggregation = wedges.n_edges_before_aggregation;
	diag->n_expanded_edges = wedges.n_expanded_edges;
	diag->n_softall_candidate_pairs = wedges.n_softall_candidate_pairs;
	diag->n_softall_filter1_pairs = wedges.n_softall_filter1_pairs;
	diag->n_softall_filter2_pairs = wedges.n_softall_filter2_pairs;
	diag->n_softall_bmap_pairs = wedges.n_softall_bmap_pairs;
	diag->n_softall_selected_raw = wedges.n_softall_selected_raw;
	diag->n_softall_gate_skip_raw = wedges.n_softall_gate_skip_raw;
	diag->n_softall_same_bin_skip_raw = wedges.n_softall_same_bin_skip_raw;
	memcpy(diag->n_softall_state_raw_count, wedges.n_softall_state_raw_count,
		   sizeof(diag->n_softall_state_raw_count));
	diag->n_skipped_self_edges = wedges.n_skipped_self_edges;
	diag->n_skipped_same_bin_bpairs = wedges.n_skipped_same_bin_bpairs;
	diag->mean_rho_train_bpair = wedges.mean_rho_train_bpair;
	diag->min_rho_train_bpair = wedges.min_rho_train_bpair;
	diag->max_rho_train_bpair = wedges.max_rho_train_bpair;
	ret = hk_blind_wedge_list_aggregate_exact(&wedges);
	if (ret != 0) {
		hk_blind_wedge_list_destroy(&wedges);
		free(prev_coords);
		return ret;
	}
	pre_diag->n_wedges = wedges.n_edges;
	diag->n_wedges = wedges.n_edges;
	hk_blind_iter_diag_validate_wedge_list(&wedges, pre_diag);
	diag->sum_wedge_k = pre_diag->sum_wedge_k;

	ret = hk_blind_relax_impl(fdg_conf, &wedges, bmap, bmap->n_beads, coords, iter_conf->unit,
							  iter_conf->relax_step, iter_conf->relax_steps, iter_conf->min_sep_unit,
							  iter_conf->lambda_sep, iter_conf->enable_repulsion,
							  iter_conf->repulsion_mode, iter_conf->repulsion_block_k_min,
							  iter_conf->chr_sep_unit, iter_conf->lambda_chr_sep,
							  anchor_map, coarse_diploid_coords, anchor_k,
							  &diag->relax_diag);
	hk_blind_wedge_list_destroy(&wedges);
	if (ret != 0) {
		free(prev_coords);
		return ret;
	}

	{
		int32_t n_chr = hk_blind_bmap_n_chr(bmap);
		uint8_t *chr_flipped = n_chr > 0? CALLOC(uint8_t, n_chr) : 0;
		if (n_chr > 0 && chr_flipped == 0) {
			free(prev_coords);
			return -1;
		}
		ret = hk_blind_temporal_gauge_stabilize_bmap(bmap, prev_coords, coords, chr_flipped, &diag->gauge_stats);
		if (ret == 0 && diag->gauge_stats.n_flipped > 0)
			ret = hk_blind_bpair_set_apply_chr_flips(set, bmap, chr_flipped, n_chr);
		free(chr_flipped);
	}
	free(prev_coords);
	if (ret != 0)
		return ret;
	diag->n_chr_flipped = diag->gauge_stats.n_flipped;

	// Remap stored p4 to the same chromosome-copy gauge as the stabilized coordinates.
	// The probabilities still come from pre-relax geometry; the next iteration, or the
	// final posterior refresh after the last iteration, recomputes p4 from stabilized coordinates.
	return 0;
}

int hk_blind_run_single_iter_cpu(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
								 const struct hk_fdg_conf *fdg_conf, fvec3_t *coords,
								 const float log_prior[HK_BLIND_N_STATE],
								 const struct hk_blind_single_iter_conf *iter_conf,
								 struct hk_blind_single_iter_diag *diag)
{
	return hk_blind_run_single_iter_cpu_impl(bmap, set, fdg_conf, coords, log_prior,
											 iter_conf, 0, 0, 0.0f, diag);
}

int hk_blind_run_single_iter_parent_anchor_cpu(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
											   const struct hk_fdg_conf *fdg_conf, fvec3_t *coords,
											   const float log_prior[HK_BLIND_N_STATE],
											   const struct hk_blind_single_iter_conf *iter_conf,
											   const struct hk_blind_coarse_to_fine_map *anchor_map,
											   const fvec3_t *coarse_diploid_coords, float anchor_k,
											   struct hk_blind_single_iter_diag *diag)
{
	return hk_blind_run_single_iter_cpu_impl(bmap, set, fdg_conf, coords, log_prior,
											 iter_conf, anchor_map, coarse_diploid_coords, anchor_k, diag);
}

void hk_blind_iter_loop_diag_init(struct hk_blind_iter_loop_diag *diag)
{
	assert(diag);
	diag->n_iter = 0;
	diag->n_completed = 0;
	diag->initial_temperature = 0.0f;
	diag->final_temperature = 0.0f;
	diag->initial_rho_train = 0.0f;
	diag->final_rho_train = 0.0f;
	diag->initial_mean_entropy = 0.0f;
	diag->final_mean_entropy = 0.0f;
	diag->initial_mean_pU = 0.0f;
	diag->final_mean_pU = 0.0f;
	diag->total_chr_flipped = 0;
	diag->n_bad_iter = 0;
	diag->n_relax_nonfinite_iter = 0;
	diag->n_coord_nonfinite = 0;
	diag->final_mean_sep = 0.0f;
	diag->final_min_sep = 0.0f;
	diag->final_max_sep = 0.0f;
	diag->final_sum_wedge_k = 0.0;
	diag->final_mean_rho_train_bpair = 0.0f;
	diag->final_min_rho_train_bpair = 0.0f;
	diag->final_max_rho_train_bpair = 0.0f;
	diag->final_n_expanded_edges = 0;
	diag->final_n_softall_candidate_pairs = 0;
	diag->final_n_softall_filter1_pairs = 0;
	diag->final_n_softall_filter2_pairs = 0;
	diag->final_n_softall_bmap_pairs = 0;
	diag->final_n_softall_selected_raw = 0;
	diag->final_n_softall_gate_skip_raw = 0;
	diag->final_n_softall_same_bin_skip_raw = 0;
	memset(diag->final_n_softall_state_raw_count, 0, sizeof(diag->final_n_softall_state_raw_count));
	diag->final_n_skipped_same_bin_bpairs = 0;
	diag->final_repulsion_energy = 0.0f;
	diag->final_repulsion_force_l1 = 0.0f;
	diag->final_n_repulsion_pairs_considered = 0;
	diag->final_n_repulsion_pairs_blocked = 0;
	diag->final_n_repulsion_pairs_active = 0;
	diag->final_anchor_energy = 0.0f;
	diag->final_anchor_force_l1 = 0.0f;
	diag->n_repulsion_nonfinite_step = 0;
	diag->n_anchor_nonfinite_step = 0;
	diag->repulsion_mode = HK_BLIND_REPULSION_NONE;
	diag->posterior_refreshed_after_final_relax = 0;
	diag->posterior_refresh_temperature = 0.0f;
	diag->posterior_refresh_mean_kl = 0.0;
	diag->posterior_refresh_top_state_switch_frac = 0.0f;
	diag->posterior_refresh_mean_pU_before = 0.0f;
	diag->posterior_refresh_mean_pU_after = 0.0f;
}

static int hk_blind_single_iter_diag_has_bad_numeric(const struct hk_blind_single_iter_diag *diag)
{
	const struct hk_blind_iter_diag *pre;

	assert(diag);
	pre = &diag->pre_relax_diag;
	return pre->n_posterior_nonfinite != 0 ||
		   pre->n_posterior_bad_sum != 0 ||
		   pre->n_posterior_out_of_range != 0 ||
		   pre->n_uncertainty_nonfinite != 0 ||
		   pre->n_uncertainty_out_of_range != 0 ||
		   pre->n_five_state_bad_sum != 0 ||
		   pre->n_wedge_nonfinite != 0 ||
		   pre->n_wedge_bad_k != 0 ||
		   pre->n_wedge_bad_d_scale != 0 ||
		   pre->sep_force_nonfinite != 0;
}

static int hk_blind_single_iter_diag_has_relax_nonfinite(const struct hk_blind_single_iter_diag *diag)
{
	assert(diag);
	return diag->relax_diag.n_nonfinite_step != 0 ||
		   diag->relax_diag.n_coord_nonfinite != 0 ||
		   diag->relax_diag.n_backbone_nonfinite_step != 0 ||
		   diag->relax_diag.n_repulsion_nonfinite_step != 0;
}

static void hk_blind_assert_single_iter_conf(const struct hk_blind_single_iter_conf *conf)
{
	assert(conf);
	assert(isfinite(conf->unit));
	assert(conf->unit > 0.0f);
	assert(isfinite(conf->d_scale));
	assert(conf->d_scale > 0.0f);
	assert(isfinite(conf->base_k));
	assert(conf->base_k >= 0.0f);
	assert(isfinite(conf->temperature));
	assert(conf->temperature > 0.0f);
	assert(isfinite(conf->rho_train));
	assert(conf->rho_train >= 0.0f);
	assert(conf->rho_train_floor == 0.0f ||
		   (isfinite(conf->rho_train_floor) &&
			conf->rho_train_floor >= 0.0f &&
			conf->rho_train_floor <= 1.0f));
	assert(isfinite(conf->min_sep_unit));
	assert(conf->min_sep_unit >= 0.0f);
	assert(isfinite(conf->lambda_sep));
	assert(conf->lambda_sep >= 0.0f);
	assert(isfinite(conf->chr_sep_unit));
	assert(conf->chr_sep_unit >= 0.0f);
	assert(isfinite(conf->lambda_chr_sep));
	assert(conf->lambda_chr_sep >= 0.0f);
	assert(isfinite(conf->relax_step));
	assert(conf->relax_step >= 0.0f);
	assert(conf->relax_steps >= 0);
	assert(conf->enable_repulsion == 0 || conf->enable_repulsion == 1);
	assert(hk_blind_repulsion_mode_valid(conf->repulsion_mode));
	assert(!conf->enable_repulsion ||
		   conf->repulsion_mode == HK_BLIND_REPULSION_N2 ||
		   conf->repulsion_mode == HK_BLIND_REPULSION_CELL);
	assert(isfinite(conf->repulsion_block_k_min));
	assert(conf->repulsion_block_k_min >= 0.0f);
	assert(hk_blind_rho_train_mode_valid(conf->rho_train_mode));
	assert(hk_blind_d_scale_mode_valid(conf->d_scale_mode));
	assert(isfinite(conf->d_scale_eps_count));
	assert(conf->d_scale_eps_count > 0.0f);
	assert(hk_blind_estep_score_mode_valid(conf->estep_score_mode));
}

static void hk_blind_assert_iter_schedule_conf(const struct hk_blind_iter_schedule_conf *conf)
{
	assert(conf);
	assert(conf->n_iter >= 0);
	hk_blind_assert_single_iter_conf(&conf->base_conf);
	assert(isfinite(conf->temperature_start));
	assert(conf->temperature_start > 0.0f);
	assert(isfinite(conf->temperature_end));
	assert(conf->temperature_end > 0.0f);
	assert(isfinite(conf->rho_train_start));
	assert(conf->rho_train_start >= 0.0f);
	assert(isfinite(conf->rho_train_end));
	assert(conf->rho_train_end >= 0.0f);
}

float hk_blind_iter_schedule_temperature_at(const struct hk_blind_iter_schedule_conf *schedule_conf, int32_t t)
{
	float a;

	hk_blind_assert_iter_schedule_conf(schedule_conf);
	assert(t >= 0);
	assert(schedule_conf->n_iter <= 0 || t < schedule_conf->n_iter);
	if (schedule_conf->n_iter <= 1)
		return schedule_conf->temperature_start;
	a = (float)t / (float)(schedule_conf->n_iter - 1);
	return schedule_conf->temperature_start * powf(schedule_conf->temperature_end / schedule_conf->temperature_start, a);
}

float hk_blind_iter_schedule_rho_train_at(const struct hk_blind_iter_schedule_conf *schedule_conf, int32_t t)
{
	float a;

	hk_blind_assert_iter_schedule_conf(schedule_conf);
	assert(t >= 0);
	assert(schedule_conf->n_iter <= 0 || t < schedule_conf->n_iter);
	if (schedule_conf->n_iter <= 1)
		return schedule_conf->rho_train_start;
	a = (float)t / (float)(schedule_conf->n_iter - 1);
	return schedule_conf->rho_train_start + (schedule_conf->rho_train_end - schedule_conf->rho_train_start) * a;
}

static int hk_blind_p4_top_state(const float p4[HK_BLIND_N_STATE])
{
	int best = 0;
	int s;
	assert(p4);
	for (s = 1; s < HK_BLIND_N_STATE; ++s)
		if (p4[s] > p4[best])
			best = s;
	return best;
}

static int hk_blind_refresh_final_posterior(struct hk_blind_bpair_set *set,
											const struct hk_fdg_conf *fdg_conf,
											const fvec3_t *coords, float unit,
											const float log_prior[HK_BLIND_N_STATE],
											float temperature,
											int estep_score_mode,
											struct hk_blind_iter_loop_diag *loop_diag)
{
	float (*old_p4)[HK_BLIND_N_STATE] = 0;
	double kl_sum = 0.0, pU_before = 0.0, pU_after = 0.0, entropy_after = 0.0;
	int32_t i, n_switch = 0;
	const float eps = 1e-12f;

	assert(set);
	assert(fdg_conf);
	assert(loop_diag);
	assert(temperature > 0.0f);
	assert(hk_blind_estep_score_mode_valid(estep_score_mode));
	if (set->n_bpairs > 0) {
		old_p4 = (float (*)[HK_BLIND_N_STATE])malloc((size_t)set->n_bpairs * sizeof(*old_p4));
		if (old_p4 == 0)
			return -1;
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		int s;
		for (s = 0; s < HK_BLIND_N_STATE; ++s)
			old_p4[i][s] = set->bpairs[i].p4[s];
		pU_before += set->bpairs[i].pU;
	}
	hk_blind_bpair_set_update_posterior_from_coords_params_score_mode(set, fdg_conf, coords,
																	  unit, log_prior, temperature,
																	  estep_score_mode);
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		int old_top = hk_blind_p4_top_state(old_p4[i]);
		int new_top = hk_blind_p4_top_state(bp->p4);
		int s;
		if (old_top != new_top) ++n_switch;
		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			float oldp = old_p4[i][s] > eps? old_p4[i][s] : eps;
			float newp = bp->p4[s] > eps? bp->p4[s] : eps;
			kl_sum += oldp * log((double)oldp / (double)newp);
		}
		pU_after += bp->pU;
		entropy_after += bp->entropy;
	}
	loop_diag->posterior_refreshed_after_final_relax = 1;
	loop_diag->posterior_refresh_temperature = temperature;
	if (set->n_bpairs > 0) {
		loop_diag->posterior_refresh_mean_kl = kl_sum / set->n_bpairs;
		loop_diag->posterior_refresh_top_state_switch_frac = (float)n_switch / (float)set->n_bpairs;
		loop_diag->posterior_refresh_mean_pU_before = (float)(pU_before / set->n_bpairs);
		loop_diag->posterior_refresh_mean_pU_after = (float)(pU_after / set->n_bpairs);
		loop_diag->final_mean_pU = (float)(pU_after / set->n_bpairs);
		loop_diag->final_mean_entropy = (float)(entropy_after / set->n_bpairs);
	}
	free(old_p4);
	return 0;
}

int hk_blind_run_iter_loop_cpu(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
							   const struct hk_fdg_conf *fdg_conf, fvec3_t *coords,
							   const float log_prior[HK_BLIND_N_STATE],
							   const struct hk_blind_iter_loop_conf *loop_conf,
							   struct hk_blind_single_iter_diag *per_iter_diag_or_null,
							   struct hk_blind_iter_loop_diag *loop_diag)
{
	struct hk_blind_single_iter_diag local_iter_diag;
	struct hk_blind_homolog_sep_stats final_sep;
	int32_t n_diploid;
	int32_t t;

	assert(bmap);
	assert(set);
	assert(fdg_conf);
	assert(loop_conf);
	assert(loop_diag);
	assert(loop_conf->n_iter >= 0);
	assert(bmap->n_beads >= 0);
	assert(bmap->n_beads <= INT32_MAX / HK_DIPLOID_N_COPY);
	hk_blind_assert_single_iter_conf(&loop_conf->single_iter_conf);
	if (bmap->n_beads > 0)
		assert(coords);
	hk_blind_iter_diag_check_bpair_set_shape(bmap, set);

	hk_blind_iter_loop_diag_init(loop_diag);
	loop_diag->n_iter = loop_conf->n_iter;
	if (loop_conf->n_iter > 0) {
		loop_diag->initial_temperature = loop_conf->single_iter_conf.temperature;
		loop_diag->final_temperature = loop_conf->single_iter_conf.temperature;
		loop_diag->initial_rho_train = loop_conf->single_iter_conf.rho_train;
		loop_diag->final_rho_train = loop_conf->single_iter_conf.rho_train;
	}
	loop_diag->repulsion_mode = loop_conf->single_iter_conf.enable_repulsion?
		loop_conf->single_iter_conf.repulsion_mode : HK_BLIND_REPULSION_NONE;
	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	if (loop_conf->n_iter == 0) {
		loop_diag->n_coord_nonfinite = hk_blind_count_nonfinite_coords(coords, n_diploid);
		if (loop_diag->n_coord_nonfinite == 0) {
			hk_blind_homolog_sep_compute_stats(bmap->n_beads, coords,
											   loop_conf->single_iter_conf.min_sep_unit * loop_conf->single_iter_conf.unit,
											   &final_sep);
			loop_diag->final_mean_sep = final_sep.mean_sep;
			loop_diag->final_min_sep = final_sep.min_sep;
			loop_diag->final_max_sep = final_sep.max_sep;
		}
		return 0;
	}

	for (t = 0; t < loop_conf->n_iter; ++t) {
		struct hk_blind_single_iter_diag *iter_diag;
		int ret;

		iter_diag = per_iter_diag_or_null? &per_iter_diag_or_null[t] : &local_iter_diag;
		ret = hk_blind_run_single_iter_cpu(bmap, set, fdg_conf, coords, log_prior,
										   &loop_conf->single_iter_conf, iter_diag);
		if (ret != 0)
			return ret;

		if (t == 0) {
			loop_diag->initial_mean_entropy = iter_diag->pre_relax_diag.mean_entropy;
			loop_diag->initial_mean_pU = iter_diag->pre_relax_diag.mean_pU;
		}
		loop_diag->final_mean_entropy = iter_diag->pre_relax_diag.mean_entropy;
		loop_diag->final_mean_pU = iter_diag->pre_relax_diag.mean_pU;
		loop_diag->final_sum_wedge_k = iter_diag->sum_wedge_k;
		loop_diag->final_mean_rho_train_bpair = iter_diag->mean_rho_train_bpair;
		loop_diag->final_min_rho_train_bpair = iter_diag->min_rho_train_bpair;
		loop_diag->final_max_rho_train_bpair = iter_diag->max_rho_train_bpair;
		loop_diag->final_n_expanded_edges = iter_diag->n_expanded_edges;
		loop_diag->final_n_softall_candidate_pairs = iter_diag->n_softall_candidate_pairs;
		loop_diag->final_n_softall_filter1_pairs = iter_diag->n_softall_filter1_pairs;
		loop_diag->final_n_softall_filter2_pairs = iter_diag->n_softall_filter2_pairs;
		loop_diag->final_n_softall_bmap_pairs = iter_diag->n_softall_bmap_pairs;
		loop_diag->final_n_softall_selected_raw = iter_diag->n_softall_selected_raw;
		loop_diag->final_n_softall_gate_skip_raw = iter_diag->n_softall_gate_skip_raw;
		loop_diag->final_n_softall_same_bin_skip_raw = iter_diag->n_softall_same_bin_skip_raw;
		memcpy(loop_diag->final_n_softall_state_raw_count, iter_diag->n_softall_state_raw_count,
			   sizeof(loop_diag->final_n_softall_state_raw_count));
		loop_diag->final_n_skipped_same_bin_bpairs = iter_diag->n_skipped_same_bin_bpairs;
		loop_diag->final_repulsion_energy = iter_diag->relax_diag.final_repulsion_energy;
		loop_diag->final_repulsion_force_l1 = iter_diag->relax_diag.final_repulsion_force_l1;
		loop_diag->final_n_repulsion_pairs_considered = iter_diag->relax_diag.final_n_repulsion_pairs_considered;
		loop_diag->final_n_repulsion_pairs_blocked = iter_diag->relax_diag.final_n_repulsion_pairs_blocked;
		loop_diag->final_n_repulsion_pairs_active = iter_diag->relax_diag.final_n_repulsion_pairs_active;
		loop_diag->final_anchor_energy = iter_diag->relax_diag.final_anchor_energy;
		loop_diag->final_anchor_force_l1 = iter_diag->relax_diag.final_anchor_force_l1;
		loop_diag->repulsion_mode = iter_diag->relax_diag.repulsion_mode;
		loop_diag->n_repulsion_nonfinite_step += iter_diag->relax_diag.n_repulsion_nonfinite_step;
		loop_diag->n_anchor_nonfinite_step += iter_diag->relax_diag.n_anchor_nonfinite_step;
		loop_diag->total_chr_flipped += iter_diag->n_chr_flipped;
		if (hk_blind_single_iter_diag_has_bad_numeric(iter_diag))
			++loop_diag->n_bad_iter;
		if (hk_blind_single_iter_diag_has_relax_nonfinite(iter_diag))
			++loop_diag->n_relax_nonfinite_iter;
		loop_diag->n_coord_nonfinite = hk_blind_count_nonfinite_coords(coords, n_diploid);
		if (loop_diag->n_coord_nonfinite != 0) {
			++loop_diag->n_relax_nonfinite_iter;
			return -1;
		}
		++loop_diag->n_completed;
	}

	if (hk_blind_refresh_final_posterior(set, fdg_conf, coords,
										 loop_conf->single_iter_conf.unit, log_prior,
										 loop_conf->single_iter_conf.temperature,
										 loop_conf->single_iter_conf.estep_score_mode,
										 loop_diag) != 0)
		return -1;
	if (hk_blind_homolog_sep_compute_stats(bmap->n_beads, coords,
										   loop_conf->single_iter_conf.min_sep_unit * loop_conf->single_iter_conf.unit,
										   &final_sep) != 0)
		return -1;
	loop_diag->final_mean_sep = final_sep.mean_sep;
	loop_diag->final_min_sep = final_sep.min_sep;
	loop_diag->final_max_sep = final_sep.max_sep;
	return 0;
}

static int hk_blind_run_iter_loop_scheduled_cpu_impl(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
													 const struct hk_fdg_conf *fdg_conf, fvec3_t *coords,
													 const float log_prior[HK_BLIND_N_STATE],
													 const struct hk_blind_iter_schedule_conf *schedule_conf,
													 const struct hk_blind_coarse_to_fine_map *anchor_map,
													 const fvec3_t *coarse_diploid_coords, float anchor_k,
													 struct hk_blind_single_iter_diag *per_iter_diag_or_null,
													 struct hk_blind_iter_loop_diag *loop_diag)
{
	struct hk_blind_single_iter_diag local_iter_diag;
	struct hk_blind_homolog_sep_stats final_sep;
	int32_t n_diploid;
	int32_t t;

	assert(bmap);
	assert(set);
	assert(fdg_conf);
	assert(loop_diag);
	assert(bmap->n_beads >= 0);
	assert(bmap->n_beads <= INT32_MAX / HK_DIPLOID_N_COPY);
	hk_blind_assert_iter_schedule_conf(schedule_conf);
	assert(isfinite(anchor_k));
	assert(anchor_k >= 0.0f);
	if (anchor_map) {
		assert(anchor_map->n_fine == bmap->n_beads);
		assert(coarse_diploid_coords || anchor_map->n_coarse == 0);
	} else {
		assert(anchor_k == 0.0f);
	}
	if (bmap->n_beads > 0)
		assert(coords);
	hk_blind_iter_diag_check_bpair_set_shape(bmap, set);

	hk_blind_iter_loop_diag_init(loop_diag);
	loop_diag->n_iter = schedule_conf->n_iter;
	if (schedule_conf->n_iter > 0) {
		loop_diag->initial_temperature = hk_blind_iter_schedule_temperature_at(schedule_conf, 0);
		loop_diag->final_temperature = hk_blind_iter_schedule_temperature_at(schedule_conf, schedule_conf->n_iter - 1);
		loop_diag->initial_rho_train = hk_blind_iter_schedule_rho_train_at(schedule_conf, 0);
		loop_diag->final_rho_train = hk_blind_iter_schedule_rho_train_at(schedule_conf, schedule_conf->n_iter - 1);
	}
	loop_diag->repulsion_mode = schedule_conf->base_conf.enable_repulsion?
		schedule_conf->base_conf.repulsion_mode : HK_BLIND_REPULSION_NONE;
	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	if (schedule_conf->n_iter == 0) {
		loop_diag->n_coord_nonfinite = hk_blind_count_nonfinite_coords(coords, n_diploid);
		if (loop_diag->n_coord_nonfinite == 0) {
			hk_blind_homolog_sep_compute_stats(bmap->n_beads, coords,
											   schedule_conf->base_conf.min_sep_unit * schedule_conf->base_conf.unit,
											   &final_sep);
			loop_diag->final_mean_sep = final_sep.mean_sep;
			loop_diag->final_min_sep = final_sep.min_sep;
			loop_diag->final_max_sep = final_sep.max_sep;
		}
		return 0;
	}

	for (t = 0; t < schedule_conf->n_iter; ++t) {
		struct hk_blind_single_iter_diag *iter_diag;
		struct hk_blind_single_iter_conf iter_conf;
		int ret;

		iter_conf = schedule_conf->base_conf;
		iter_conf.temperature = hk_blind_iter_schedule_temperature_at(schedule_conf, t);
		iter_conf.rho_train = hk_blind_iter_schedule_rho_train_at(schedule_conf, t);
		iter_diag = per_iter_diag_or_null? &per_iter_diag_or_null[t] : &local_iter_diag;
		ret = hk_blind_run_single_iter_cpu_impl(bmap, set, fdg_conf, coords, log_prior, &iter_conf,
												anchor_map, coarse_diploid_coords, anchor_k, iter_diag);
		if (ret != 0)
			return ret;

		if (t == 0) {
			loop_diag->initial_mean_entropy = iter_diag->pre_relax_diag.mean_entropy;
			loop_diag->initial_mean_pU = iter_diag->pre_relax_diag.mean_pU;
		}
		loop_diag->final_mean_entropy = iter_diag->pre_relax_diag.mean_entropy;
		loop_diag->final_mean_pU = iter_diag->pre_relax_diag.mean_pU;
		loop_diag->final_sum_wedge_k = iter_diag->sum_wedge_k;
		loop_diag->final_mean_rho_train_bpair = iter_diag->mean_rho_train_bpair;
		loop_diag->final_min_rho_train_bpair = iter_diag->min_rho_train_bpair;
		loop_diag->final_max_rho_train_bpair = iter_diag->max_rho_train_bpair;
		loop_diag->final_n_expanded_edges = iter_diag->n_expanded_edges;
		loop_diag->final_n_softall_candidate_pairs = iter_diag->n_softall_candidate_pairs;
		loop_diag->final_n_softall_filter1_pairs = iter_diag->n_softall_filter1_pairs;
		loop_diag->final_n_softall_filter2_pairs = iter_diag->n_softall_filter2_pairs;
		loop_diag->final_n_softall_bmap_pairs = iter_diag->n_softall_bmap_pairs;
		loop_diag->final_n_softall_selected_raw = iter_diag->n_softall_selected_raw;
		loop_diag->final_n_softall_gate_skip_raw = iter_diag->n_softall_gate_skip_raw;
		loop_diag->final_n_softall_same_bin_skip_raw = iter_diag->n_softall_same_bin_skip_raw;
		memcpy(loop_diag->final_n_softall_state_raw_count, iter_diag->n_softall_state_raw_count,
			   sizeof(loop_diag->final_n_softall_state_raw_count));
		loop_diag->final_n_skipped_same_bin_bpairs = iter_diag->n_skipped_same_bin_bpairs;
		loop_diag->final_repulsion_energy = iter_diag->relax_diag.final_repulsion_energy;
		loop_diag->final_repulsion_force_l1 = iter_diag->relax_diag.final_repulsion_force_l1;
		loop_diag->final_n_repulsion_pairs_considered = iter_diag->relax_diag.final_n_repulsion_pairs_considered;
		loop_diag->final_n_repulsion_pairs_blocked = iter_diag->relax_diag.final_n_repulsion_pairs_blocked;
		loop_diag->final_n_repulsion_pairs_active = iter_diag->relax_diag.final_n_repulsion_pairs_active;
		loop_diag->final_anchor_energy = iter_diag->relax_diag.final_anchor_energy;
		loop_diag->final_anchor_force_l1 = iter_diag->relax_diag.final_anchor_force_l1;
		loop_diag->repulsion_mode = iter_diag->relax_diag.repulsion_mode;
		loop_diag->n_repulsion_nonfinite_step += iter_diag->relax_diag.n_repulsion_nonfinite_step;
		loop_diag->n_anchor_nonfinite_step += iter_diag->relax_diag.n_anchor_nonfinite_step;
		loop_diag->total_chr_flipped += iter_diag->n_chr_flipped;
		if (hk_blind_single_iter_diag_has_bad_numeric(iter_diag))
			++loop_diag->n_bad_iter;
		if (hk_blind_single_iter_diag_has_relax_nonfinite(iter_diag))
			++loop_diag->n_relax_nonfinite_iter;
		loop_diag->n_coord_nonfinite = hk_blind_count_nonfinite_coords(coords, n_diploid);
		if (loop_diag->n_coord_nonfinite != 0) {
			++loop_diag->n_relax_nonfinite_iter;
			return -1;
		}
		++loop_diag->n_completed;
	}

	if (hk_blind_refresh_final_posterior(set, fdg_conf, coords,
										 schedule_conf->base_conf.unit, log_prior,
										 hk_blind_iter_schedule_temperature_at(schedule_conf, schedule_conf->n_iter - 1),
										 schedule_conf->base_conf.estep_score_mode,
										 loop_diag) != 0)
		return -1;
	if (hk_blind_homolog_sep_compute_stats(bmap->n_beads, coords,
										   schedule_conf->base_conf.min_sep_unit * schedule_conf->base_conf.unit,
										   &final_sep) != 0)
		return -1;
	loop_diag->final_mean_sep = final_sep.mean_sep;
	loop_diag->final_min_sep = final_sep.min_sep;
	loop_diag->final_max_sep = final_sep.max_sep;
	return 0;
}

int hk_blind_run_iter_loop_scheduled_cpu(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
										 const struct hk_fdg_conf *fdg_conf, fvec3_t *coords,
										 const float log_prior[HK_BLIND_N_STATE],
										 const struct hk_blind_iter_schedule_conf *schedule_conf,
										 struct hk_blind_single_iter_diag *per_iter_diag_or_null,
										 struct hk_blind_iter_loop_diag *loop_diag)
{
	return hk_blind_run_iter_loop_scheduled_cpu_impl(bmap, set, fdg_conf, coords, log_prior,
													 schedule_conf, 0, 0, 0.0f,
													 per_iter_diag_or_null, loop_diag);
}

int hk_blind_run_iter_loop_scheduled_parent_anchor_cpu(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
													   const struct hk_fdg_conf *fdg_conf, fvec3_t *coords,
													   const float log_prior[HK_BLIND_N_STATE],
													   const struct hk_blind_iter_schedule_conf *schedule_conf,
													   const struct hk_blind_coarse_to_fine_map *anchor_map,
													   const fvec3_t *coarse_diploid_coords, float anchor_k,
													   struct hk_blind_single_iter_diag *per_iter_diag_or_null,
													   struct hk_blind_iter_loop_diag *loop_diag)
{
	return hk_blind_run_iter_loop_scheduled_cpu_impl(bmap, set, fdg_conf, coords, log_prior,
													 schedule_conf, anchor_map, coarse_diploid_coords, anchor_k,
													 per_iter_diag_or_null, loop_diag);
}

static const char *hk_blind_bead_chr_label(const struct hk_bmap *bmap, int32_t bid, char fallback[32])
{
	const struct hk_bead *bead;

	assert(bmap);
	assert(bid >= 0 && bid < bmap->n_beads);
	assert(fallback);
	bead = &bmap->beads[bid];
	if (bmap->d && bmap->d->name && bead->chr >= 0 && bead->chr < bmap->d->n && bmap->d->name[bead->chr])
		return bmap->d->name[bead->chr];
	snprintf(fallback, 32, "%d", bead->chr);
	return fallback;
}

static const char *hk_blind_chr_label(const struct hk_bmap *bmap, int32_t chr, char fallback[32])
{
	assert(bmap);
	assert(chr >= 0);
	assert(fallback);
	if (bmap->d && bmap->d->name && chr < bmap->d->n && bmap->d->name[chr])
		return bmap->d->name[chr];
	snprintf(fallback, 32, "%d", chr);
	return fallback;
}

int hk_blind_write_bpair_posterior_tsv(FILE *fp, const struct hk_bmap *bmap,
									   const struct hk_blind_bpair_set *set)
{
	int32_t i;

	assert(fp);
	assert(bmap);
	assert(set);
	assert(bmap->n_beads >= 0);
	if (bmap->n_beads > 0)
		assert(bmap->beads);

	if (fprintf(fp, "chr1\tstart1\tend1\tchr2\tstart2\tend2\tbid1\tbid2\tn_raw\t"
					"base_d_scale\tbase_k\tp00\tp01\tp10\tp11\tpU\tpsame_raw\tpcross_raw\t"
					"entropy\tmargin\tpmax\trho_output\tcontact_class\n") < 0)
		return -1;
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		int32_t bid0 = bp->key.bid[0], bid1 = bp->key.bid[1];
		const struct hk_bead *b0, *b1;
		char chr0_buf[32], chr1_buf[32];
		const char *chr0, *chr1;
		float psame_raw, pcross_raw;

		assert(bid0 >= 0 && bid0 < bmap->n_beads);
		assert(bid1 >= 0 && bid1 < bmap->n_beads);
		b0 = &bmap->beads[bid0];
		b1 = &bmap->beads[bid1];
		chr0 = hk_blind_bead_chr_label(bmap, bid0, chr0_buf);
		chr1 = hk_blind_bead_chr_label(bmap, bid1, chr1_buf);
		psame_raw = bp->p4[HK_BLIND_STATE_00] + bp->p4[HK_BLIND_STATE_11];
		pcross_raw = bp->p4[HK_BLIND_STATE_01] + bp->p4[HK_BLIND_STATE_10];
		if (fprintf(fp, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t"
						"%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%s\n",
					chr0, b0->st, b0->en, chr1, b1->st, b1->en, bid0, bid1, bp->n_raw,
					bp->base_d_scale, bp->base_k,
					bp->p4[HK_BLIND_STATE_00], bp->p4[HK_BLIND_STATE_01],
					bp->p4[HK_BLIND_STATE_10], bp->p4[HK_BLIND_STATE_11],
					bp->pU, psame_raw, pcross_raw, bp->entropy, bp->margin,
					bp->pmax, bp->rho_output,
					hk_blind_contact_class_name(hk_blind_bpair_contact_class_or_default(bp))) < 0)
			return -1;
	}
	return ferror(fp)? -1 : 0;
}

int hk_blind_write_raw_contact_posterior_tsv(FILE *fp, const struct hk_blind_pair *raw, int32_t n_raw,
											 const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set)
{
	int32_t i;

	assert(fp);
	assert(n_raw >= 0);
	assert(bmap);
	assert(set);
	assert(set->n_raw == n_raw);
	assert(bmap->n_beads >= 0);
	if (n_raw > 0) {
		assert(raw);
		assert(set->raw2binned);
		assert(set->bpairs);
		assert(bmap->beads);
		assert(bmap->d);
		assert(bmap->offcnt);
	}

	if (fprintf(fp, "raw_id\tchr1\tpos1\tchr2\tpos2\tbid1_raw\tbid2_raw\tbpair_id\t"
					"bid1_canonical\tbid2_canonical\tswapped\tp00\tp01\tp10\tp11\tpU\t"
					"psame_raw\tpcross_raw\tentropy\tmargin\tpmax\trho_output\tcontact_class\n") < 0)
		return -1;
	for (i = 0; i < n_raw; ++i) {
		const struct hk_blind_raw2binned *r2b = &set->raw2binned[i];
		const struct hk_blind_bpair *bp;
		float raw_p4[HK_BLIND_N_STATE];
		float psame_raw, pcross_raw;
		int32_t raw_bid[2];
		char chr0_buf[32], chr1_buf[32];
		const char *chr0, *chr1;

		assert(r2b->bpair_id >= 0 && r2b->bpair_id < set->n_bpairs);
		bp = &set->bpairs[r2b->bpair_id];
		hk_blind_pair_to_bids(bmap, &raw[i], raw_bid);
		hk_blind_p4_to_raw_order(bp->p4, r2b->swapped, raw_p4);
		psame_raw = raw_p4[HK_BLIND_STATE_00] + raw_p4[HK_BLIND_STATE_11];
		pcross_raw = raw_p4[HK_BLIND_STATE_01] + raw_p4[HK_BLIND_STATE_10];
		chr0 = hk_blind_chr_label(bmap, raw[i].chr[0], chr0_buf);
		chr1 = hk_blind_chr_label(bmap, raw[i].chr[1], chr1_buf);
		if (fprintf(fp, "%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%u\t"
						"%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%s\n",
					i, chr0, raw[i].pos[0], chr1, raw[i].pos[1], raw_bid[0], raw_bid[1],
					r2b->bpair_id, bp->key.bid[0], bp->key.bid[1], (unsigned)r2b->swapped,
					raw_p4[HK_BLIND_STATE_00], raw_p4[HK_BLIND_STATE_01],
					raw_p4[HK_BLIND_STATE_10], raw_p4[HK_BLIND_STATE_11],
					bp->pU, psame_raw, pcross_raw, bp->entropy, bp->margin,
					bp->pmax, bp->rho_output,
					hk_blind_contact_class_name(hk_blind_bpair_contact_class_or_default(bp))) < 0)
			return -1;
	}
	return ferror(fp)? -1 : 0;
}

int hk_blind_write_diploid_coords_tsv(FILE *fp, const struct hk_bmap *bmap, const fvec3_t *coords)
{
	int32_t i;
	int copy;

	assert(fp);
	assert(bmap);
	assert(bmap->n_beads >= 0);
	if (bmap->n_beads > 0) {
		assert(bmap->beads);
		assert(coords);
	}

	if (fprintf(fp, "chr\tstart\tend\tbid\tcopy\tdiploid_bid\tx\ty\tz\n") < 0)
		return -1;
	for (i = 0; i < bmap->n_beads; ++i) {
		const struct hk_bead *bead = &bmap->beads[i];
		char chr_buf[32];
		const char *chr = hk_blind_bead_chr_label(bmap, i, chr_buf);
		for (copy = 0; copy < HK_DIPLOID_N_COPY; ++copy) {
			int32_t diploid_bid = hk_diploid_bid(i, copy);
			if (fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\n",
						chr, bead->st, bead->en, i, copy, diploid_bid,
						coords[diploid_bid][0], coords[diploid_bid][1], coords[diploid_bid][2]) < 0)
				return -1;
		}
	}
	return ferror(fp)? -1 : 0;
}

int hk_blind_write_diploid_coords_tsv_gz(const char *path, const struct hk_bmap *bmap, const fvec3_t *coords)
{
	gzFile fp;
	int32_t i;
	int copy;
	int failed = 0;
	int close_ret;

	assert(path);
	assert(bmap);
	assert(bmap->n_beads >= 0);
	if (bmap->n_beads > 0) {
		assert(bmap->beads);
		assert(coords);
	}

	fp = gzopen(path, "wb");
	if (fp == 0)
		return -1;
	if (gzprintf(fp, "chr\tstart\tend\tbid\tcopy\tdiploid_bid\tx\ty\tz\n") < 0)
		failed = 1;
	for (i = 0; !failed && i < bmap->n_beads; ++i) {
		const struct hk_bead *bead = &bmap->beads[i];
		char chr_buf[32];
		const char *chr = hk_blind_bead_chr_label(bmap, i, chr_buf);
		for (copy = 0; copy < HK_DIPLOID_N_COPY; ++copy) {
			int32_t diploid_bid = hk_diploid_bid(i, copy);
			if (gzprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\n",
						 chr, bead->st, bead->en, i, copy, diploid_bid,
						 coords[diploid_bid][0], coords[diploid_bid][1],
						 coords[diploid_bid][2]) < 0) {
				failed = 1;
				break;
			}
		}
	}
	close_ret = gzclose(fp);
	if (close_ret != Z_OK)
		failed = 1;
	return failed? -1 : 0;
}

int hk_blind_write_iter_loop_diag_tsv(FILE *fp, const struct hk_blind_iter_loop_diag *diag)
{
	assert(fp);
	assert(diag);
	if (fprintf(fp, "n_iter\tn_completed\tinitial_mean_entropy\tfinal_mean_entropy\t"
					"initial_mean_pU\tfinal_mean_pU\tinitial_temperature\tfinal_temperature\t"
					"initial_rho_train\tfinal_rho_train\ttotal_chr_flipped\tn_bad_iter\t"
					"n_relax_nonfinite_iter\tn_coord_nonfinite\tfinal_mean_sep\tfinal_min_sep\t"
					"final_max_sep\tfinal_sum_wedge_k\tfinal_mean_rho_train_bpair\t"
					"final_min_rho_train_bpair\tfinal_max_rho_train_bpair\t"
					"final_n_expanded_edges\tfinal_n_softall_candidate_pairs\t"
					"final_n_softall_filter1_pairs\tfinal_n_softall_filter2_pairs\t"
					"final_n_softall_bmap_pairs\t"
					"final_n_softall_selected_raw\tfinal_n_softall_gate_skip_raw\t"
					"final_n_softall_same_bin_skip_raw\t"
					"final_n_softall_state00_raw\tfinal_n_softall_state01_raw\t"
					"final_n_softall_state10_raw\tfinal_n_softall_state11_raw\t"
					"final_n_skipped_same_bin_bpairs\tfinal_repulsion_energy\t"
					"final_repulsion_force_l1\tfinal_n_repulsion_pairs_considered\t"
					"final_n_repulsion_pairs_blocked\tfinal_n_repulsion_pairs_active\t"
					"n_repulsion_nonfinite_step\trepulsion_mode\t"
					"posterior_refreshed_after_final_relax\tposterior_refresh_temperature\t"
					"posterior_refresh_mean_kl\tposterior_refresh_top_state_switch_frac\t"
					"posterior_refresh_mean_pU_before\tposterior_refresh_mean_pU_after\n") < 0)
		return -1;
	if (fprintf(fp, "%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t"
					"%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.17g\t%.9g\t%.9g\t%.9g\t"
					"%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t"
					"%lld\t%lld\t%lld\t%lld\t%.9g\t%.9g\t%lld\t%lld\t%lld\t%d\t%d\t"
					"%d\t%.9g\t%.17g\t%.9g\t%.9g\t%.9g\n",
				diag->n_iter, diag->n_completed,
				diag->initial_mean_entropy, diag->final_mean_entropy,
				diag->initial_mean_pU, diag->final_mean_pU,
				diag->initial_temperature, diag->final_temperature,
				diag->initial_rho_train, diag->final_rho_train,
				diag->total_chr_flipped, diag->n_bad_iter,
				diag->n_relax_nonfinite_iter, diag->n_coord_nonfinite,
				diag->final_mean_sep, diag->final_min_sep, diag->final_max_sep,
				diag->final_sum_wedge_k, diag->final_mean_rho_train_bpair,
				diag->final_min_rho_train_bpair, diag->final_max_rho_train_bpair,
				(long long)diag->final_n_expanded_edges,
				(long long)diag->final_n_softall_candidate_pairs,
				(long long)diag->final_n_softall_filter1_pairs,
				(long long)diag->final_n_softall_filter2_pairs,
				(long long)diag->final_n_softall_bmap_pairs,
				(long long)diag->final_n_softall_selected_raw,
				(long long)diag->final_n_softall_gate_skip_raw,
				(long long)diag->final_n_softall_same_bin_skip_raw,
				(long long)diag->final_n_softall_state_raw_count[HK_BLIND_STATE_00],
				(long long)diag->final_n_softall_state_raw_count[HK_BLIND_STATE_01],
				(long long)diag->final_n_softall_state_raw_count[HK_BLIND_STATE_10],
				(long long)diag->final_n_softall_state_raw_count[HK_BLIND_STATE_11],
				(long long)diag->final_n_skipped_same_bin_bpairs,
				diag->final_repulsion_energy,
				diag->final_repulsion_force_l1,
				(long long)diag->final_n_repulsion_pairs_considered,
				(long long)diag->final_n_repulsion_pairs_blocked,
				(long long)diag->final_n_repulsion_pairs_active,
				diag->n_repulsion_nonfinite_step,
				diag->repulsion_mode, diag->posterior_refreshed_after_final_relax,
				diag->posterior_refresh_temperature, diag->posterior_refresh_mean_kl,
				diag->posterior_refresh_top_state_switch_frac,
				diag->posterior_refresh_mean_pU_before,
				diag->posterior_refresh_mean_pU_after) < 0)
		return -1;
	return ferror(fp)? -1 : 0;
}

struct hk_blind_bpair_set *hk_blind_bpair_set_build(const struct hk_bmap *m, int32_t n_raw, const struct hk_blind_pair *raw)
{
	struct hk_blind_bpair_set *set;
	struct hk_blind_bpair_aux *aux = 0;
	int32_t i, bpair_id = -1;

	assert(m);
	assert(n_raw >= 0);
	assert(n_raw == 0 || raw);

	set = CALLOC(struct hk_blind_bpair_set, 1);
	set->n_raw = n_raw;
	set->same_bin_filter_enabled = 1;
	if (n_raw == 0) return set;

	set->raw2binned = CALLOC(struct hk_blind_raw2binned, n_raw);
	set->raw = CALLOC(struct hk_blind_pair, n_raw);
	set->bpairs = CALLOC(struct hk_blind_bpair, n_raw);
	aux = MALLOC(struct hk_blind_bpair_aux, n_raw);
	if (set->raw == 0 || set->raw2binned == 0 ||
		set->bpairs == 0 || aux == 0) {
		free(aux);
		hk_blind_bpair_set_destroy(set);
		return 0;
	}
	memcpy(set->raw, raw, (size_t)n_raw * sizeof(*raw));

	// Invalid coordinates follow hk_bmap_pos2bid(): assert/abort instead of silent skipping.
	for (i = 0; i < n_raw; ++i) {
		hk_blind_pair_to_bpair_key(m, &raw[i], &aux[i].key, &aux[i].swapped);
		aux[i].raw_id = i;
		if (hk_blind_bpair_key_is_same_bin(&aux[i].key))
			++set->n_raw_same_bin_excluded;
	}

	qsort(aux, n_raw, sizeof(*aux), hk_blind_bpair_aux_cmp);
	for (i = 0; i < n_raw; ++i) {
		if (i == 0 || !hk_blind_bpair_key_eq(&aux[i].key, &aux[i - 1].key)) {
			bpair_id = set->n_bpairs++;
			hk_blind_bpair_init(&set->bpairs[bpair_id], &aux[i].key);
			set->bpairs[bpair_id].contact_class = hk_blind_bpair_key_contact_class(m, &aux[i].key);
			if (set->bpairs[bpair_id].contact_class == HK_BLIND_CONTACT_TRANS)
				++set->n_bpair_trans;
			else
				++set->n_bpair_cis;
			if (hk_blind_bpair_key_is_same_bin(&aux[i].key))
				++set->n_bpair_same_bin_excluded;
		}
		++set->bpairs[bpair_id].n_raw;
		if (set->bpairs[bpair_id].contact_class == HK_BLIND_CONTACT_TRANS)
			++set->n_raw_trans;
		else
			++set->n_raw_cis;
		hk_blind_raw2binned_set(&set->raw2binned[aux[i].raw_id], bpair_id, aux[i].swapped);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		assert(set->bpairs[i].n_raw > 0);
		set->bpairs[i].base_d_scale = hk_blind_bpair_base_d_scale_from_n_raw(set->bpairs[i].n_raw);
		set->bpairs[i].base_k = 1.0f;
	}

	free(aux);
	return set;
}

static int hk_blind_bpair_key_find_aux(const struct hk_blind_base_k_aux *aux, int32_t n_aux,
									   const struct hk_blind_bpair_key *key)
{
	int32_t lo = 0, hi = n_aux - 1;
	assert(key);
	while (lo <= hi) {
		int32_t mid = lo + (hi - lo) / 2;
		struct hk_blind_base_k_aux t;
		int cmp;
		t.key = *key;
		t.max_nei = 0;
		cmp = hk_blind_base_k_aux_cmp(&aux[mid], &t);
		if (cmp == 0) return mid;
		if (cmp < 0) lo = mid + 1;
		else hi = mid - 1;
	}
	return -1;
}

int hk_blind_bpair_set_apply_neighbor_median_base_k(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set)
{
	const double a_third = 1.0 / 3.0;
	struct hk_blind_base_k_aux *aux = 0;
	int32_t *nei = 0;
	int32_t i, n_aux = 0, n_nei = 0, median_nei;

	assert(set);
	for (i = 0; i < set->n_bpairs; ++i)
		set->bpairs[i].base_k = 1.0f;
	if (set->n_bpairs == 0)
		return 0;
	if (bmap == 0 || bmap->n_pairs <= 0 || bmap->pairs == 0)
		return 0;

	aux = CALLOC(struct hk_blind_base_k_aux, bmap->n_pairs);
	nei = CALLOC(int32_t, bmap->n_pairs);
	if (aux == 0 || nei == 0) {
		free(aux);
		free(nei);
		return -1;
	}

	for (i = 0; i < bmap->n_pairs; ++i) {
		const struct hk_bpair *p = &bmap->pairs[i];
		int32_t bid0 = p->bid[0], bid1 = p->bid[1];
		if (bid0 <= bid1) {
			aux[n_aux].key.bid[0] = bid0;
			aux[n_aux].key.bid[1] = bid1;
		} else {
			aux[n_aux].key.bid[0] = bid1;
			aux[n_aux].key.bid[1] = bid0;
		}
		aux[n_aux].max_nei = p->max_nei;
		++n_aux;
	}
	qsort(aux, n_aux, sizeof(*aux), hk_blind_base_k_aux_cmp);
	for (i = 0; i < n_aux; ++i) {
		if (i == 0 || !hk_blind_bpair_key_eq(&aux[i].key, &aux[n_nei - 1].key)) {
			aux[n_nei++] = aux[i];
		} else if (aux[i].max_nei > aux[n_nei - 1].max_nei) {
			aux[n_nei - 1].max_nei = aux[i].max_nei;
		}
	}
	for (i = 0; i < n_nei; ++i)
		nei[i] = aux[i].max_nei;
	median_nei = n_nei > 0? ks_ksmall_int32_t(n_nei, nei, (size_t)(n_nei / 2)) : 0;
	if (median_nei <= 0) {
		free(aux);
		free(nei);
		return 0;
	}

	for (i = 0; i < set->n_bpairs; ++i) {
		int32_t j = hk_blind_bpair_key_find_aux(aux, n_nei, &set->bpairs[i].key);
		int32_t max_nei = j >= 0? aux[j].max_nei : median_nei;
		if (max_nei >= median_nei)
			set->bpairs[i].base_k = 1.0f;
		else if (max_nei <= 0)
			set->bpairs[i].base_k = 0.0f;
		else
			set->bpairs[i].base_k = powf((float)((double)max_nei / (double)median_nei), (float)a_third);
	}

	free(aux);
	free(nei);
	return 0;
}

int hk_blind_bpair_set_apply_base_k_mode(const struct hk_bmap *bmap, struct hk_blind_bpair_set *set,
										 int base_k_mode)
{
	int32_t i;
	assert(set);
	assert(hk_blind_base_k_mode_valid(base_k_mode));
	if (base_k_mode == HK_BLIND_BASE_K_UNIFORM) {
		for (i = 0; i < set->n_bpairs; ++i)
			set->bpairs[i].base_k = 1.0f;
		return 0;
	}
	return hk_blind_bpair_set_apply_neighbor_median_base_k(bmap, set);
}

void hk_blind_heldout_diag_init(struct hk_blind_heldout_diag *diag)
{
	assert(diag);
	memset(diag, 0, sizeof(*diag));
}

int hk_blind_eval_heldout_bpair_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
									const struct hk_fdg_conf *conf, const fvec3_t *coords,
									float unit, float temperature,
									struct hk_blind_heldout_diag *diag)
{
	double sum_expected_energy = 0.0;
	double sum_min_energy = 0.0;
	double sum_entropy = 0.0;
	double sum_pU = 0.0;
	double sum_best_norm_dist = 0.0;
	double sum_short = 0.0;
	int64_t weight_sum = 0;
	int32_t i;

	assert(bmap);
	assert(set);
	assert(conf);
	assert(coords);
	assert(unit > 0.0f);
	assert(temperature > 0.0f);
	assert(diag);

	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		float energy[HK_BLIND_N_STATE];
		float p4[HK_BLIND_N_STATE];
		float entropy, pmax, margin, rho_output, pU;
		float dist[HK_BLIND_N_STATE];
		double expected_energy = 0.0;
		float min_energy;
		float best_norm_dist;
		int32_t i0, i1, j0, j1;
		int s;
		int64_t w;

		assert(bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		assert(bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
		if (hk_blind_bpair_is_same_bin(bp)) {
			++diag->n_bpair_same_bin_skipped;
			diag->n_raw_same_bin_skipped += bp->n_raw;
			continue;
		}

		hk_blind_bpair_posterior_from_coords(conf, bp->key.bid[0], bp->key.bid[1],
											 coords, unit, bp->base_d_scale, bp->base_k,
											 bp->log_prior, temperature,
											 energy, p4, &entropy, &pmax, &margin,
											 &rho_output, &pU);
		i0 = hk_diploid_bid(bp->key.bid[0], HK_DIPLOID_COPY0);
		i1 = hk_diploid_bid(bp->key.bid[0], HK_DIPLOID_COPY1);
		j0 = hk_diploid_bid(bp->key.bid[1], HK_DIPLOID_COPY0);
		j1 = hk_diploid_bid(bp->key.bid[1], HK_DIPLOID_COPY1);
		dist[HK_BLIND_STATE_00] = hk_blind_coord_dist(coords[i0], coords[j0]);
		dist[HK_BLIND_STATE_01] = hk_blind_coord_dist(coords[i0], coords[j1]);
		dist[HK_BLIND_STATE_10] = hk_blind_coord_dist(coords[i1], coords[j0]);
		dist[HK_BLIND_STATE_11] = hk_blind_coord_dist(coords[i1], coords[j1]);

		min_energy = energy[0];
		best_norm_dist = dist[0] / (unit * bp->base_d_scale);
		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			float norm_dist = dist[s] / (unit * bp->base_d_scale);
			expected_energy += (double)p4[s] * energy[s];
			if (s == 0 || energy[s] < min_energy) min_energy = energy[s];
			if (s == 0 || norm_dist < best_norm_dist) best_norm_dist = norm_dist;
		}

		w = bp->n_raw > 0? bp->n_raw : 1;
		sum_expected_energy += (double)w * expected_energy;
		sum_min_energy += (double)w * min_energy;
		sum_entropy += (double)w * entropy;
		sum_pU += (double)w * pU;
		sum_best_norm_dist += (double)w * best_norm_dist;
		if (best_norm_dist <= conf->d_c2)
			sum_short += (double)w;
		weight_sum += w;
		++diag->n_bpair_eval;
	}

	diag->n_raw_eval = weight_sum;
	if (weight_sum > 0) {
		diag->mean_expected_energy = sum_expected_energy / (double)weight_sum;
		diag->mean_min_energy = sum_min_energy / (double)weight_sum;
		diag->mean_entropy = sum_entropy / (double)weight_sum;
		diag->mean_pU = sum_pU / (double)weight_sum;
		diag->mean_best_normalized_distance = sum_best_norm_dist / (double)weight_sum;
		diag->short_distance_frac = sum_short / (double)weight_sum;
	}
	return 0;
}

void hk_blind_bpair_set_destroy(struct hk_blind_bpair_set *set)
{
	if (set == 0) return;
	free(set->bpairs);
	free(set->raw2binned);
	free(set->raw);
	free(set);
}

static const char *hk_blind_bmap_chr_name(const struct hk_bmap *bmap, int32_t chr, char fallback[32])
{
	assert(fallback);
	if (bmap && bmap->d && chr >= 0 && chr < bmap->d->n && bmap->d->name && bmap->d->name[chr])
		return bmap->d->name[chr];
	snprintf(fallback, 32, "%d", chr);
	return fallback;
}

int hk_blind_write_bmap_summary_tsv(FILE *fp, const struct hk_bmap *bmap,
									const int32_t *n_child_1mb_or_null)
{
	int32_t i;
	assert(fp);
	assert(bmap);
	assert(bmap->n_beads >= 0);
	assert(bmap->n_beads == 0 || bmap->beads);
	if (fprintf(fp, "chr\tbid\tstart\tend\tbead_size_bp\tn_child_1mb\n") < 0)
		return -1;
	for (i = 0; i < bmap->n_beads; ++i) {
		const struct hk_bead *b = &bmap->beads[i];
		char chr_buf[32];
		const char *chr = hk_blind_bmap_chr_name(bmap, b->chr, chr_buf);
		int32_t n_child = n_child_1mb_or_null? n_child_1mb_or_null[i] : -1;
		if (fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\n",
					chr, i, b->st, b->en, b->en - b->st, n_child) < 0)
			return -1;
	}
	return ferror(fp)? -1 : 0;
}

struct hk_blind_coarse_to_fine_map *hk_blind_coarse_to_fine_map_build(const struct hk_bmap *coarse,
																	  const struct hk_bmap *fine)
{
	struct hk_blind_coarse_to_fine_map *map = 0;
	int32_t i;
	assert(coarse);
	assert(fine);
	assert(coarse->d && fine->d);
	assert(coarse->d->n == fine->d->n);
	assert(coarse->n_beads >= 0);
	assert(fine->n_beads >= 0);
	assert(coarse->n_beads == 0 || coarse->beads);
	assert(fine->n_beads == 0 || fine->beads);

	map = CALLOC(struct hk_blind_coarse_to_fine_map, 1);
	map->n_coarse = coarse->n_beads;
	map->n_fine = fine->n_beads;
		if (fine->n_beads > 0) {
			map->fine_to_coarse = MALLOC(int32_t, fine->n_beads);
			map->child_rank = MALLOC(int32_t, fine->n_beads);
			if (map->fine_to_coarse == 0 || map->child_rank == 0)
				goto fail;
		}
		if (coarse->n_beads > 0) {
			map->parent_first_fine = MALLOC(int32_t, coarse->n_beads);
			map->parent_child_count = CALLOC(int32_t, coarse->n_beads);
			if (map->parent_first_fine == 0 || map->parent_child_count == 0)
				goto fail;
			for (i = 0; i < coarse->n_beads; ++i)
				map->parent_first_fine[i] = -1;
		}
	for (i = 0; i < fine->n_beads; ++i) {
		const struct hk_bead *fb = &fine->beads[i];
		const struct hk_bead *cb;
		int64_t mid64;
		int32_t mid, parent;
		assert(fb->en > fb->st);
		mid64 = (int64_t)fb->st + ((int64_t)fb->en - fb->st) / 2;
		if (mid64 >= fb->en) mid64 = fb->en - 1;
		mid = (int32_t)mid64;
			parent = hk_bmap_pos2bid(coarse, fb->chr, mid);
			if (parent < 0 || parent >= coarse->n_beads)
				goto fail;
			cb = &coarse->beads[parent];
			if (cb->chr != fb->chr || fb->st < cb->st || fb->en > cb->en)
				goto fail;
		map->fine_to_coarse[i] = parent;
		map->child_rank[i] = map->parent_child_count[parent]++;
		if (map->parent_first_fine[parent] < 0)
			map->parent_first_fine[parent] = i;
	}
	return map;

fail:
	hk_blind_coarse_to_fine_map_destroy(map);
	return 0;
}

void hk_blind_coarse_to_fine_map_destroy(struct hk_blind_coarse_to_fine_map *map)
{
	if (map == 0) return;
	free(map->fine_to_coarse);
	free(map->child_rank);
	free(map->parent_first_fine);
	free(map->parent_child_count);
	free(map);
}

int hk_blind_write_coarse_to_fine_map_tsv(FILE *fp, const struct hk_bmap *coarse,
										  const struct hk_bmap *fine,
										  const struct hk_blind_coarse_to_fine_map *map)
{
	int32_t i;
	assert(fp);
	assert(coarse);
	assert(fine);
	assert(map);
	assert(map->n_coarse == coarse->n_beads);
	assert(map->n_fine == fine->n_beads);
	if (fprintf(fp, "coarse_chr\tcoarse_bid\tcoarse_start\tcoarse_end\tfine_bid\tfine_start\tfine_end\tchild_index\tn_children_for_coarse\tchild_offset_bp\n") < 0)
		return -1;
	for (i = 0; i < fine->n_beads; ++i) {
		int32_t parent = map->fine_to_coarse[i];
		const struct hk_bead *cb, *fb;
		char chr_buf[32];
		const char *chr;
		assert(parent >= 0 && parent < coarse->n_beads);
		cb = &coarse->beads[parent];
		fb = &fine->beads[i];
		chr = hk_blind_bmap_chr_name(coarse, cb->chr, chr_buf);
		if (fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					chr, parent, cb->st, cb->en, i, fb->st, fb->en,
					map->child_rank[i], map->parent_child_count[parent],
					fb->st - cb->st) < 0)
			return -1;
	}
	return ferror(fp)? -1 : 0;
}

float hk_blind_density_exposure_from_child_counts(int32_t child_count0, int32_t child_count1,
												  int same_parent)
{
	float exposure;
	if (child_count0 < 1) child_count0 = 1;
	if (child_count1 < 1) child_count1 = 1;
	if (same_parent) {
		if (child_count0 < 2) return 1.0f;
		exposure = 0.5f * (float)child_count0 * (float)(child_count0 - 1);
	} else {
		exposure = (float)child_count0 * (float)child_count1;
	}
	if (!isfinite(exposure) || exposure < 1.0f)
		exposure = 1.0f;
	return exposure;
}

float hk_blind_bpair_density_exposure(const struct hk_blind_bpair *bp)
{
	assert(bp);
	if (isfinite(bp->density_exposure) && bp->density_exposure >= 1.0f)
		return bp->density_exposure;
	return 1.0f;
}

float hk_blind_bpair_density_normalized_raw_count(const struct hk_blind_bpair *bp,
												  float eps_count)
{
	float n_eff;
	assert(bp);
	assert(bp->n_raw > 0);
	assert(isfinite(eps_count));
	assert(eps_count > 0.0f);
	n_eff = (float)bp->n_raw / hk_blind_bpair_density_exposure(bp);
	if (!isfinite(n_eff) || n_eff < eps_count)
		n_eff = eps_count;
	return n_eff;
}

float hk_blind_bpair_capped_density_raw_count(const struct hk_blind_bpair *bp,
											  float eps_count)
{
	float n_eff = hk_blind_bpair_density_normalized_raw_count(bp, eps_count);
	if (n_eff > HK_BLIND_D_SCALE_DEFAULT_COUNT_CAP)
		n_eff = HK_BLIND_D_SCALE_DEFAULT_COUNT_CAP;
	if (n_eff < eps_count)
		n_eff = eps_count;
	return n_eff;
}

int hk_blind_bpair_set_apply_density_exposure_from_map(struct hk_blind_bpair_set *set,
													   const struct hk_blind_coarse_to_fine_map *map)
{
	int32_t i;
	assert(set);
	assert(map);
	assert(map->n_coarse >= 0);
	assert(map->n_coarse == 0 || map->parent_child_count);
	for (i = 0; i < set->n_bpairs; ++i) {
		struct hk_blind_bpair *bp = &set->bpairs[i];
		int32_t b0 = bp->key.bid[0], b1 = bp->key.bid[1];
		int32_t c0, c1;
		if (b0 < 0 || b0 >= map->n_coarse || b1 < 0 || b1 >= map->n_coarse)
			return -1;
		c0 = map->parent_child_count[b0];
		c1 = map->parent_child_count[b1];
		bp->density_child_count[0] = c0;
		bp->density_child_count[1] = c1;
		bp->density_exposure = hk_blind_density_exposure_from_child_counts(c0, c1, b0 == b1);
	}
	return 0;
}

int hk_blind_bpair_set_apply_density_d_scale_from_map(struct hk_blind_bpair_set *set,
													  const struct hk_blind_coarse_to_fine_map *map,
													  int d_scale_mode, float eps_count)
{
	int32_t i;
	assert(set);
	assert(map);
	assert(hk_blind_d_scale_mode_valid(d_scale_mode));
	assert(isfinite(eps_count));
	assert(eps_count > 0.0f);
	if (hk_blind_bpair_set_apply_density_exposure_from_map(set, map) != 0)
		return -1;
	for (i = 0; i < set->n_bpairs; ++i) {
		struct hk_blind_bpair *bp = &set->bpairs[i];
		float n_eff = (float)bp->n_raw;
		if (d_scale_mode == HK_BLIND_D_SCALE_DENSITY_NORMALIZED_RAW_COUNT)
			n_eff = hk_blind_bpair_density_normalized_raw_count(bp, eps_count);
		else if (d_scale_mode == HK_BLIND_D_SCALE_CAPPED_DENSITY_RAW_COUNT)
			n_eff = hk_blind_bpair_capped_density_raw_count(bp, eps_count);
		if (!isfinite(n_eff) || n_eff < eps_count)
			n_eff = eps_count;
		bp->base_d_scale = powf(n_eff, -1.0f / 3.0f);
		if (!isfinite(bp->base_d_scale) || bp->base_d_scale <= 0.0f)
			return -1;
	}
	return 0;
}

static void hk_blind_lift_axis_for_parent(const struct hk_bmap *coarse, const fvec3_t *coarse_diploid_coords,
										  int32_t parent, int32_t copy, float axis[3])
{
	int32_t prev = -1, next = -1, bid;
	float norm;
	int a;
	assert(coarse);
	assert(coarse_diploid_coords);
	assert(axis);
	assert(parent >= 0 && parent < coarse->n_beads);
	for (bid = parent - 1; bid >= 0; --bid) {
		if (coarse->beads[bid].chr != coarse->beads[parent].chr) break;
		prev = bid;
		break;
	}
	for (bid = parent + 1; bid < coarse->n_beads; ++bid) {
		if (coarse->beads[bid].chr != coarse->beads[parent].chr) break;
		next = bid;
		break;
	}
	if (prev >= 0 && next >= 0) {
		int32_t dprev = hk_diploid_bid(prev, copy);
		int32_t dnext = hk_diploid_bid(next, copy);
		for (a = 0; a < 3; ++a)
			axis[a] = coarse_diploid_coords[dnext][a] - coarse_diploid_coords[dprev][a];
	} else if (next >= 0) {
		int32_t dparent = hk_diploid_bid(parent, copy);
		int32_t dnext = hk_diploid_bid(next, copy);
		for (a = 0; a < 3; ++a)
			axis[a] = coarse_diploid_coords[dnext][a] - coarse_diploid_coords[dparent][a];
	} else if (prev >= 0) {
		int32_t dprev = hk_diploid_bid(prev, copy);
		int32_t dparent = hk_diploid_bid(parent, copy);
		for (a = 0; a < 3; ++a)
			axis[a] = coarse_diploid_coords[dparent][a] - coarse_diploid_coords[dprev][a];
	} else {
		axis[0] = 1.0f; axis[1] = 0.0f; axis[2] = 0.0f;
	}
	norm = sqrtf(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
	if (!isfinite(norm) || norm <= 0.0f) {
		axis[0] = 1.0f; axis[1] = 0.0f; axis[2] = 0.0f;
	} else {
		for (a = 0; a < 3; ++a)
			axis[a] /= norm;
	}
}

int hk_blind_lift_from_4mb(const struct hk_bmap *coarse,
						   const struct hk_blind_coarse_to_fine_map *map,
						   const fvec3_t *coarse_diploid_coords,
						   fvec3_t *fine_diploid_coords,
						   float child_offset_step)
{
	int32_t fine_bid;
	assert(coarse);
	assert(map);
	assert(coarse_diploid_coords);
	assert(fine_diploid_coords);
	assert(map->n_coarse == coarse->n_beads);
	assert(isfinite(child_offset_step));
	assert(child_offset_step >= 0.0f);
	for (fine_bid = 0; fine_bid < map->n_fine; ++fine_bid) {
		int32_t parent = map->fine_to_coarse[fine_bid];
		int32_t n_child, rank, copy;
		float centered;
		assert(parent >= 0 && parent < map->n_coarse);
		n_child = map->parent_child_count[parent];
		rank = map->child_rank[fine_bid];
		assert(n_child > 0);
		assert(rank >= 0 && rank < n_child);
		centered = (float)rank - 0.5f * (float)(n_child - 1);
		for (copy = 0; copy < HK_DIPLOID_N_COPY; ++copy) {
			float axis[3];
			int32_t src = hk_diploid_bid(parent, copy);
			int32_t dst = hk_diploid_bid(fine_bid, copy);
			int a;
			hk_blind_lift_axis_for_parent(coarse, coarse_diploid_coords, parent, copy, axis);
			for (a = 0; a < 3; ++a) {
				fine_diploid_coords[dst][a] = coarse_diploid_coords[src][a] +
					centered * child_offset_step * axis[a];
				if (!isfinite(fine_diploid_coords[dst][a]))
					return -1;
			}
		}
	}
	return 0;
}

float hk_blind_parent_centroid_anchor_accumulate_force(const struct hk_blind_coarse_to_fine_map *map,
													   const fvec3_t *fine_diploid_coords,
													   const fvec3_t *coarse_diploid_coords,
													   float k_anchor,
													   fvec3_t *fine_force,
													   int32_t *n_nonfinite,
													   float *force_l1_or_null)
{
	double energy = 0.0;
	double force_l1 = 0.0;
	int32_t parent, copy;
	assert(map);
	assert(fine_diploid_coords || map->n_fine == 0);
	assert(coarse_diploid_coords || map->n_coarse == 0);
	assert(fine_force || map->n_fine == 0);
	assert(isfinite(k_anchor));
	assert(k_anchor >= 0.0f);
	if (n_nonfinite) *n_nonfinite = 0;
	if (force_l1_or_null) *force_l1_or_null = 0.0f;
	for (parent = 0; parent < map->n_coarse; ++parent) {
		int32_t n_child = map->parent_child_count[parent];
		int32_t first_fine = map->parent_first_fine[parent];
		if (n_child <= 0) continue;
		if (first_fine < 0 || first_fine + n_child > map->n_fine) {
			if (n_nonfinite) ++(*n_nonfinite);
			continue;
		}
		for (copy = 0; copy < HK_DIPLOID_N_COPY; ++copy) {
			double centroid[3] = {0.0, 0.0, 0.0};
			float delta[3], f[3];
			int32_t child;
			int a;
			for (child = 0; child < n_child; ++child) {
				int32_t fine_bid = first_fine + child;
				if (map->fine_to_coarse[fine_bid] != parent) {
					if (n_nonfinite) ++(*n_nonfinite);
					goto next_parent_copy;
				}
				for (a = 0; a < 3; ++a)
					centroid[a] += fine_diploid_coords[hk_diploid_bid(fine_bid, copy)][a];
			}
			for (a = 0; a < 3; ++a) {
				centroid[a] /= n_child;
				delta[a] = coarse_diploid_coords[hk_diploid_bid(parent, copy)][a] - (float)centroid[a];
				f[a] = k_anchor * delta[a] / (float)n_child;
				if (!isfinite(delta[a]) || !isfinite(f[a])) {
					if (n_nonfinite) ++(*n_nonfinite);
					goto next_parent_copy;
				}
			}
			energy += 0.5 * (double)k_anchor *
				((double)delta[0] * delta[0] + (double)delta[1] * delta[1] + (double)delta[2] * delta[2]);
			for (child = 0; child < n_child; ++child) {
				int32_t fine_bid = first_fine + child;
				for (a = 0; a < 3; ++a) {
					fine_force[hk_diploid_bid(fine_bid, copy)][a] += f[a];
					force_l1 += fabs((double)f[a]);
				}
			}
next_parent_copy:
			;
		}
	}
	if (force_l1_or_null) *force_l1_or_null = (float)force_l1;
	return (float)energy;
}
