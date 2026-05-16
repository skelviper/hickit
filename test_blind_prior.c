#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "hickit.h"

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

static void set_test_bmap(struct hk_bmap *bmap, struct hk_sdict *d, int32_t len[2],
						  struct hk_bead beads[5], uint64_t offcnt[3])
{
	memset(d, 0, sizeof(*d));
	d->n = 2;
	d->len = len;
	len[0] = 3000000;
	len[1] = 2000000;
	beads[0].chr = 0; beads[0].st = 0;       beads[0].en = 1000000;
	beads[1].chr = 0; beads[1].st = 1000000; beads[1].en = 2000000;
	beads[2].chr = 0; beads[2].st = 2000000; beads[2].en = 3000000;
	beads[3].chr = 1; beads[3].st = 0;       beads[3].en = 1000000;
	beads[4].chr = 1; beads[4].st = 1000000; beads[4].en = 2000000;
	offcnt[0] = ((uint64_t)0 << 32) | 3u;
	offcnt[1] = ((uint64_t)3 << 32) | 2u;
	offcnt[2] = ((uint64_t)5 << 32) | 0u;
	memset(bmap, 0, sizeof(*bmap));
	bmap->d = d;
	bmap->n_beads = 5;
	bmap->beads = beads;
	bmap->offcnt = offcnt;
}

static void set_bpair(struct hk_blind_bpair *bp, int32_t bid0, int32_t bid1, int32_t n_raw)
{
	memset(bp, 0, sizeof(*bp));
	bp->key.bid[0] = bid0;
	bp->key.bid[1] = bid1;
	bp->n_raw = n_raw;
	bp->base_d_scale = powf((float)n_raw, -1.0f / 3.0f);
	bp->base_k = 1.0f;
	hk_blind_init_uniform_log_prior(bp->log_prior);
}

static void set_bpair_set(struct hk_blind_bpair_set *set, struct hk_blind_bpair *bpairs, int32_t n_bpairs)
{
	memset(set, 0, sizeof(*set));
	set->bpairs = bpairs;
	set->n_bpairs = n_bpairs;
}

static float prior_p(const struct hk_blind_bpair *bp, int state)
{
	return expf(bp->log_prior[state]);
}

static int check_strong_cis_same_copy_prior(void)
{
	struct hk_bmap bmap;
	struct hk_sdict d;
	int32_t len[2];
	struct hk_bead beads[5];
	uint64_t offcnt[3];
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_bpair_set set;
	struct hk_blind_prior_diag diag;
	float same, cross;
	int failed = 0;

	set_test_bmap(&bmap, &d, len, beads, offcnt);
	set_bpair(&bpairs[0], 0, 1, 100);
	set_bpair(&bpairs[1], 0, 3, 1);
	set_bpair_set(&set, bpairs, 2);
	failed |= check_true("cis-inter prior ret",
						 hk_blind_bpair_set_init_cis_inter_ratio_prior(&bmap, &set, 1e-6f, &diag) == 0);
	same = prior_p(&bpairs[0], HK_BLIND_STATE_00) + prior_p(&bpairs[0], HK_BLIND_STATE_11);
	cross = prior_p(&bpairs[0], HK_BLIND_STATE_01) + prior_p(&bpairs[0], HK_BLIND_STATE_10);
	failed |= check_true("strong cis same-copy mass", same > cross);
	failed |= check_true("strong cis alpha below uniform", diag.alpha_min < 0.5f);
	failed |= check_close("trans uniform p00", prior_p(&bpairs[1], HK_BLIND_STATE_00), 0.25f);
	failed |= check_close("trans uniform p01", prior_p(&bpairs[1], HK_BLIND_STATE_01), 0.25f);
	failed |= check_true("inter density uses possible pairs", diag.possible_inter == 6.0);
	return failed;
}

static int check_uniform_limit_and_finiteness(void)
{
	struct hk_bmap bmap;
	struct hk_sdict d;
	int32_t len[2];
	struct hk_bead beads[5];
	uint64_t offcnt[3];
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_bpair_set set;
	struct hk_blind_prior_diag diag;
	float same, cross, sum;
	int failed = 0, i, s;

	set_test_bmap(&bmap, &d, len, beads, offcnt);
	set_bpair(&bpairs[0], 0, 1, 1);
	set_bpair(&bpairs[1], 0, 3, 4);
	set_bpair_set(&set, bpairs, 2);
	failed |= check_true("uniform-limit prior ret",
						 hk_blind_bpair_set_init_cis_inter_ratio_prior(&bmap, &set, 1e-6f, &diag) == 0);
	same = prior_p(&bpairs[0], HK_BLIND_STATE_00) + prior_p(&bpairs[0], HK_BLIND_STATE_11);
	cross = prior_p(&bpairs[0], HK_BLIND_STATE_01) + prior_p(&bpairs[0], HK_BLIND_STATE_10);
	failed |= check_close("uniform-limit same mass", same, 0.5f);
	failed |= check_close("uniform-limit cross mass", cross, 0.5f);
	for (i = 0; i < 2; ++i) {
		sum = 0.0f;
		for (s = 0; s < HK_BLIND_N_STATE; ++s) {
			failed |= check_true("prior finite", isfinite(bpairs[i].log_prior[s]));
			sum += prior_p(&bpairs[i], s);
		}
		failed |= check_close("prior sums to one", sum, 1.0f);
	}
	return failed;
}

static int check_uniform_mode_is_old_behavior(void)
{
	struct hk_blind_bpair bpairs[1];
	struct hk_blind_bpair_set set;
	int failed = 0, s;

	set_bpair(&bpairs[0], 0, 1, 10);
	set_bpair_set(&set, bpairs, 1);
	hk_blind_bpair_set_init_uniform_prior(&set);
	for (s = 0; s < HK_BLIND_N_STATE; ++s)
		failed |= check_close("uniform mode prior", prior_p(&bpairs[0], s), 0.25f);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_strong_cis_same_copy_prior();
	failed |= check_uniform_limit_and_finiteness();
	failed |= check_uniform_mode_is_old_behavior();
	return failed != 0;
}
