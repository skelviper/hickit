#include <stdint.h>
#include <stdio.h>
#include <math.h>
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

static int check_u8(const char *label, uint8_t got, uint8_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %u, expected %u\n", label, got, expected);
		return 1;
	}
	return 0;
}

static int check_float(const char *label, float got, float expected)
{
	if (fabsf(got - expected) > 1e-6f) {
		fprintf(stderr, "%s: got %.8g, expected %.8g\n", label, got, expected);
		return 1;
	}
	return 0;
}

static struct hk_blind_pair make_blind_pair(int32_t chr0, int32_t pos0, int32_t chr1, int32_t pos1)
{
	struct hk_blind_pair p;
	p.chr[0] = chr0;
	p.chr[1] = chr1;
	p.pos[0] = pos0;
	p.pos[1] = pos1;
	p.strand[0] = 1;
	p.strand[1] = -1;
	return p;
}

static int check_key(const char *label, const struct hk_blind_bpair_key *key, int32_t bid0, int32_t bid1)
{
	int failed = 0;
	char buf[64];
	snprintf(buf, sizeof(buf), "%s bid0", label);
	failed |= check_i32(buf, key->bid[0], bid0);
	snprintf(buf, sizeof(buf), "%s bid1", label);
	failed |= check_i32(buf, key->bid[1], bid1);
	return failed;
}

static int check_pos_to_bid(const struct hk_bmap *b, int32_t pos, int32_t expected_bid)
{
	struct hk_blind_pair p = make_blind_pair(0, pos, 0, pos);
	int32_t bid[2];
	int failed = 0;
	hk_blind_pair_to_bids(b, &p, bid);
	failed |= check_i32("pos endpoint0 bid", bid[0], expected_bid);
	failed |= check_i32("pos endpoint1 bid", bid[1], expected_bid);
	return failed;
}

static int check_uniform_posterior(const char *label, const struct hk_blind_bpair *p)
{
	int failed = 0, i;
	char buf[64];
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		snprintf(buf, sizeof(buf), "%s p4[%d]", label, i);
		failed |= check_float(buf, p->p4[i], 0.25f);
	}
	failed |= check_float("entropy", p->entropy, logf(4.0f));
	failed |= check_float("pmax", p->pmax, 0.25f);
	failed |= check_float("margin", p->margin, 0.0f);
	failed |= check_float("rho_output", p->rho_output, 0.0f);
	failed |= check_float("pU", p->pU, 1.0f);
	return failed;
}

static int check_bpair_params(const char *label, const struct hk_blind_bpair *p, int32_t n_raw, float base_d_scale)
{
	int failed = 0;
	char buf[64];
	snprintf(buf, sizeof(buf), "%s n_raw", label);
	failed |= check_i32(buf, p->n_raw, n_raw);
	snprintf(buf, sizeof(buf), "%s base_d_scale", label);
	failed |= check_float(buf, p->base_d_scale, base_d_scale);
	snprintf(buf, sizeof(buf), "%s base_k", label);
	failed |= check_float(buf, p->base_k, 1.0f);
	return failed;
}

static int check_canonical_keys(void)
{
	struct hk_blind_bpair_key key;
	struct hk_blind_raw2binned raw2binned;
	uint8_t swapped = 99;
	int failed = 0;

	hk_blind_bpair_key_from_bids(&key, 2, 5, &swapped);
	failed |= check_key("2,5 key", &key, 2, 5);
	failed |= check_u8("2,5 swapped", swapped, 0);

	hk_blind_bpair_key_from_bids(&key, 5, 2, &swapped);
	failed |= check_key("5,2 key", &key, 2, 5);
	failed |= check_u8("5,2 swapped", swapped, 1);

	hk_blind_bpair_key_from_bids(&key, 3, 3, &swapped);
	failed |= check_key("3,3 key", &key, 3, 3);
	failed |= check_u8("3,3 swapped", swapped, 0);

	hk_blind_raw2binned_set(&raw2binned, 17, swapped);
	failed |= check_i32("raw2binned bpair_id", raw2binned.bpair_id, 17);
	failed |= check_u8("raw2binned swapped", raw2binned.swapped, 0);

	hk_blind_raw2binned_set(&raw2binned, 17, 1);
	failed |= check_i32("raw2binned reversed bpair_id", raw2binned.bpair_id, 17);
	failed |= check_u8("raw2binned reversed swapped", raw2binned.swapped, 1);
	return failed;
}

static int check_bpair_set_counting(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[3];
	struct hk_blind_bpair_set *set;
	int failed = 0;

	raw[0] = make_blind_pair(0, 2000000, 0, 5000000);
	raw[1] = make_blind_pair(0, 5000000, 0, 2000000);
	raw[2] = make_blind_pair(0, 2000000, 0, 5000000);
	set = hk_blind_bpair_set_build(b, 3, raw);

	failed |= check_i32("counting n_raw", set->n_raw, 3);
	failed |= check_i32("counting n_bpairs", set->n_bpairs, 1);
	failed |= check_key("counting key", &set->bpairs[0].key, 2, 5);
	failed |= check_i32("counting bpair n_raw", set->bpairs[0].n_raw, 3);
	failed |= check_i32("counting raw0 bpair_id", set->raw2binned[0].bpair_id, 0);
	failed |= check_i32("counting raw1 bpair_id", set->raw2binned[1].bpair_id, 0);
	failed |= check_i32("counting raw2 bpair_id", set->raw2binned[2].bpair_id, 0);
	failed |= check_u8("counting raw0 swapped", set->raw2binned[0].swapped, 0);
	failed |= check_u8("counting raw1 swapped", set->raw2binned[1].swapped, 1);
	failed |= check_u8("counting raw2 swapped", set->raw2binned[2].swapped, 0);
	failed |= check_bpair_params("counting bpair params", &set->bpairs[0], 3, powf(3.0f, -1.0f / 3.0f));
	failed |= check_uniform_posterior("counting posterior", &set->bpairs[0]);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_bpair_set_count_params(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[36];
	struct hk_blind_bpair_set *set;
	int failed = 0;
	int32_t i, n = 0;
	float before_d_scale, before_base_k;

	raw[n++] = make_blind_pair(0, 2000000, 0, 5000000);
	for (i = 0; i < 8; ++i)
		raw[n++] = make_blind_pair(0, 3000000, 0, 7000000);
	for (i = 0; i < 27; ++i)
		raw[n++] = make_blind_pair(0, 4000000, 0, 6000000);
	set = hk_blind_bpair_set_build(b, n, raw);

	failed |= check_i32("count params n_raw", set->n_raw, 36);
	failed |= check_i32("count params n_bpairs", set->n_bpairs, 3);
	failed |= check_key("count params key0", &set->bpairs[0].key, 2, 5);
	failed |= check_key("count params key1", &set->bpairs[1].key, 3, 7);
	failed |= check_key("count params key2", &set->bpairs[2].key, 4, 6);
	failed |= check_bpair_params("count params n1", &set->bpairs[0], 1, 1.0f);
	failed |= check_bpair_params("count params n8", &set->bpairs[1], 8, 0.5f);
	failed |= check_bpair_params("count params n27", &set->bpairs[2], 27, 1.0f / 3.0f);

	before_d_scale = set->bpairs[1].base_d_scale;
	before_base_k = set->bpairs[1].base_k;
	set->bpairs[1].p4[HK_BLIND_STATE_00] = 0.7f;
	set->bpairs[1].p4[HK_BLIND_STATE_01] = 0.1f;
	set->bpairs[1].p4[HK_BLIND_STATE_10] = 0.15f;
	set->bpairs[1].p4[HK_BLIND_STATE_11] = 0.05f;
	failed |= check_float("count params posterior-independent d_scale", set->bpairs[1].base_d_scale, before_d_scale);
	failed |= check_float("count params posterior-independent base_k", set->bpairs[1].base_k, before_base_k);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_bpair_set_sorted_keys(const struct hk_bmap *b)
{
	struct hk_blind_pair raw[3];
	struct hk_blind_bpair_set *set;
	int failed = 0;

	raw[0] = make_blind_pair(0, 3000000, 0, 7000000);
	raw[1] = make_blind_pair(0, 2000000, 0, 5000000);
	raw[2] = make_blind_pair(0, 3000000, 0, 3000000);
	set = hk_blind_bpair_set_build(b, 3, raw);

	failed |= check_i32("sorted n_raw", set->n_raw, 3);
	failed |= check_i32("sorted n_bpairs", set->n_bpairs, 3);
	failed |= check_key("sorted key0", &set->bpairs[0].key, 2, 5);
	failed |= check_key("sorted key1", &set->bpairs[1].key, 3, 3);
	failed |= check_key("sorted key2", &set->bpairs[2].key, 3, 7);
	failed |= check_i32("sorted key0 n_raw", set->bpairs[0].n_raw, 1);
	failed |= check_i32("sorted key1 n_raw", set->bpairs[1].n_raw, 1);
	failed |= check_i32("sorted key2 n_raw", set->bpairs[2].n_raw, 1);
	failed |= check_i32("sorted raw0 bpair_id", set->raw2binned[0].bpair_id, 2);
	failed |= check_i32("sorted raw1 bpair_id", set->raw2binned[1].bpair_id, 0);
	failed |= check_i32("sorted raw2 bpair_id", set->raw2binned[2].bpair_id, 1);
	failed |= check_u8("sorted raw0 swapped", set->raw2binned[0].swapped, 0);
	failed |= check_u8("sorted raw1 swapped", set->raw2binned[1].swapped, 0);
	failed |= check_u8("sorted raw2 swapped", set->raw2binned[2].swapped, 0);
	failed |= check_uniform_posterior("sorted posterior0", &set->bpairs[0]);
	failed |= check_uniform_posterior("sorted posterior1", &set->bpairs[1]);
	failed |= check_uniform_posterior("sorted posterior2", &set->bpairs[2]);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_bpair_set_contact_classes(void)
{
	struct hk_sdict dict;
	int32_t len[2] = {2000000, 2000000};
	struct hk_bead beads[4];
	uint64_t offcnt[3];
	struct hk_bmap bmap;
	struct hk_blind_pair raw[2];
	struct hk_blind_bpair_set *set;
	int failed = 0;

	dict.n = 2;
	dict.m = 2;
	dict.name = 0;
	dict.len = len;
	dict.h = 0;
	beads[0].chr = 0; beads[0].st = 0;       beads[0].en = 1000000;
	beads[1].chr = 0; beads[1].st = 1000000; beads[1].en = 2000000;
	beads[2].chr = 1; beads[2].st = 0;       beads[2].en = 1000000;
	beads[3].chr = 1; beads[3].st = 1000000; beads[3].en = 2000000;
	offcnt[0] = ((uint64_t)0 << 32) | 2u;
	offcnt[1] = ((uint64_t)2 << 32) | 2u;
	offcnt[2] = ((uint64_t)4 << 32) | 0u;
	bmap.n_beads = 4;
	bmap.n_pairs = 0;
	bmap.unit = 1.0f;
	bmap.d = &dict;
	bmap.beads = beads;
	bmap.offcnt = offcnt;
	bmap.pairs = 0;
	bmap.x = 0;
	bmap.feat = 0;
	bmap.cpg = 0;
	bmap.gc_bias = 0;
	bmap.gc_corrected = 0;

	raw[0] = make_blind_pair(0, 1000, 0, 1001000);
	raw[1] = make_blind_pair(0, 1000, 1, 1000);
	set = hk_blind_bpair_set_build(&bmap, 2, raw);

	failed |= check_i32("class n_raw_cis", (int32_t)set->n_raw_cis, 1);
	failed |= check_i32("class n_raw_trans", (int32_t)set->n_raw_trans, 1);
	failed |= check_i32("class n_bpair_cis", set->n_bpair_cis, 1);
	failed |= check_i32("class n_bpair_trans", set->n_bpair_trans, 1);
	failed |= check_i32("class first cis", set->bpairs[0].contact_class, HK_BLIND_CONTACT_CIS);
	failed |= check_i32("class second trans", set->bpairs[1].contact_class, HK_BLIND_CONTACT_TRANS);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_neighbor_median_base_k(const struct hk_bmap *b)
{
	struct hk_bmap bmap = *b;
	struct hk_bpair pairs[3];
	struct hk_blind_pair raw[3];
	struct hk_blind_bpair_set *set;
	struct hk_blind_base_k_stats stats;
	int failed = 0;

	memset(pairs, 0, sizeof(pairs));
	pairs[0].bid[0] = 5; pairs[0].bid[1] = 2; pairs[0].max_nei = 8;
	pairs[1].bid[0] = 3; pairs[1].bid[1] = 7; pairs[1].max_nei = 1;
	pairs[2].bid[0] = 4; pairs[2].bid[1] = 6; pairs[2].max_nei = 27;
	bmap.n_pairs = 3;
	bmap.pairs = pairs;

	raw[0] = make_blind_pair(0, 2000000, 0, 5000000);
	raw[1] = make_blind_pair(0, 3000000, 0, 7000000);
	raw[2] = make_blind_pair(0, 4000000, 0, 6000000);
	set = hk_blind_bpair_set_build(&bmap, 3, raw);
	failed |= check_i32("neighbor base_k apply", hk_blind_bpair_set_apply_neighbor_median_base_k(&bmap, set), 0);
	failed |= check_float("neighbor base_k reversed canonical", set->bpairs[0].base_k, 1.0f);
	failed |= check_float("neighbor base_k weak", set->bpairs[1].base_k, 0.5f);
	failed |= check_float("neighbor base_k saturated", set->bpairs[2].base_k, 1.0f);
	hk_blind_bpair_set_base_k_stats(set, &stats);
	failed |= check_i32("neighbor stats n", stats.n, 3);
	failed |= check_i32("neighbor stats finite", stats.n_nonfinite, 0);
	failed |= check_float("neighbor stats min", stats.min, 0.5f);
	failed |= check_float("neighbor stats max", stats.max, 1.0f);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int check_phase_like_columns_are_ignored(const struct hk_map *m, const struct hk_bmap *b)
{
	struct hk_blind_bpair_key first_key;
	uint8_t first_swapped = 0;
	int have_first = 0;
	int failed = 0;
	int32_t i;

	for (i = 0; i < m->n_pairs; ++i) {
		struct hk_blind_pair blind;
		struct hk_blind_bpair_key key;
		uint8_t swapped = 99;

		hk_blind_pair_from_pair(&blind, &m->pairs[i]);
		hk_blind_pair_to_bpair_key(b, &blind, &key, &swapped);

		failed |= check_key("fixture key", &key, 3, 3);
		failed |= check_u8("fixture swapped", swapped, 0);
		if (!have_first) {
			first_key = key;
			first_swapped = swapped;
			have_first = 1;
		} else {
			failed |= check_key("fixture same key", &key, first_key.bid[0], first_key.bid[1]);
			failed |= check_u8("fixture same swapped", swapped, first_swapped);
		}
	}
	return failed;
}

static int check_bpair_set_phase_like_columns_are_ignored(const struct hk_map *m, const struct hk_bmap *b)
{
	struct hk_blind_pair raw[4];
	struct hk_blind_bpair_set *set;
	int failed = 0;
	int32_t i;

	failed |= check_i32("fixture n_pairs for set", m->n_pairs, 4);
	for (i = 0; i < m->n_pairs && i < 4; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);
	set = hk_blind_bpair_set_build(b, 4, raw);

	failed |= check_i32("fixture set n_raw", set->n_raw, 4);
	failed |= check_i32("fixture set n_bpairs", set->n_bpairs, 1);
	failed |= check_key("fixture set key", &set->bpairs[0].key, 3, 3);
	failed |= check_i32("fixture set bpair n_raw", set->bpairs[0].n_raw, 4);
	for (i = 0; i < 4; ++i) {
		failed |= check_i32("fixture set raw bpair_id", set->raw2binned[i].bpair_id, 0);
		failed |= check_u8("fixture set raw swapped", set->raw2binned[i].swapped, 0);
	}
	failed |= check_uniform_posterior("fixture set posterior", &set->bpairs[0]);

	hk_blind_bpair_set_destroy(set);
	return failed;
}

int main(void)
{
	struct hk_map *m;
	struct hk_bmap *b;
	struct hk_blind_pair p;
	struct hk_blind_bpair_key key;
	int32_t bid[2];
	uint8_t swapped;
	int failed = 0;

	hk_verbose = 0;
	m = hk_map_read("testdata/p9016_blind_fixture.pairs");
	if (m == 0) {
		fprintf(stderr, "failed to read fixture\n");
		return 1;
	}
	b = hk_bmap_gen(m->d, m->n_pairs, m->pairs, 1000000, 1);
	if (b == 0) {
		fprintf(stderr, "failed to build bmap\n");
		hk_map_destroy(m);
		return 1;
	}

	failed |= check_pos_to_bid(b, 0, 0);
	failed |= check_pos_to_bid(b, 999999, 0);
	failed |= check_pos_to_bid(b, 1000000, 1);
	failed |= check_pos_to_bid(b, 1999999, 1);

	p = make_blind_pair(0, 2000000, 0, 5000000);
	hk_blind_pair_to_bids(b, &p, bid);
	failed |= check_i32("pair bid0", bid[0], 2);
	failed |= check_i32("pair bid1", bid[1], 5);
	hk_blind_pair_to_bpair_key(b, &p, &key, &swapped);
	failed |= check_key("pair key", &key, 2, 5);
	failed |= check_u8("pair swapped", swapped, 0);

	p = make_blind_pair(0, 5000000, 0, 2000000);
	hk_blind_pair_to_bpair_key(b, &p, &key, &swapped);
	failed |= check_key("reversed pair key", &key, 2, 5);
	failed |= check_u8("reversed pair swapped", swapped, 1);

	failed |= check_canonical_keys();
	failed |= check_phase_like_columns_are_ignored(m, b);
	failed |= check_bpair_set_counting(b);
	failed |= check_bpair_set_count_params(b);
	failed |= check_neighbor_median_base_k(b);
	failed |= check_bpair_set_sorted_keys(b);
	failed |= check_bpair_set_contact_classes();
	failed |= check_bpair_set_phase_like_columns_are_ignored(m, b);

	hk_bmap_destroy(b);
	hk_map_destroy(m);
	return failed != 0;
}
