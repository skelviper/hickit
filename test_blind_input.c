#include <stdint.h>
#include <stdio.h>
#include "hickit.h"

static int check_i32(const char *label, int32_t got, int32_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %d, expected %d\n", label, got, expected);
		return 1;
	}
	return 0;
}

static int check_blind_pair(const struct hk_blind_pair *p)
{
	int failed = 0;
	failed |= check_i32("blind chr0", p->chr[0], 0);
	failed |= check_i32("blind chr1", p->chr[1], 0);
	failed |= check_i32("blind pos0", p->pos[0], 3000312);
	failed |= check_i32("blind pos1", p->pos[1], 3013909);
	failed |= check_i32("blind strand0", p->strand[0], 1);
	failed |= check_i32("blind strand1", p->strand[1], 1);
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

int main(void)
{
	struct hk_map *m;
	struct hk_blind_pair first;
	int failed = 0;
	int count_dot_dot = 0, count_0_dot = 0, count_dot_1 = 0, count_0_1 = 0;
	int32_t i;

	m = hk_map_read("testdata/p9016_blind_fixture.pairs");
	if (m == 0) {
		fprintf(stderr, "failed to read fixture\n");
		return 1;
	}
	failed |= check_i32("n_pairs", m->n_pairs, 4);

	for (i = 0; i < m->n_pairs; ++i) {
		const struct hk_pair *src = &m->pairs[i];
		struct hk_blind_pair blind;

		if (src->phase[0] < 0 && src->phase[1] < 0) ++count_dot_dot;
		else if (src->phase[0] == 0 && src->phase[1] < 0) ++count_0_dot;
		else if (src->phase[0] < 0 && src->phase[1] == 1) ++count_dot_1;
		else if (src->phase[0] == 0 && src->phase[1] == 1) ++count_0_1;

		hk_blind_pair_from_pair(&blind, src);
		failed |= check_blind_pair(&blind);
		if (i == 0) first = blind;
		else failed |= same_blind_pair(&first, &blind);
	}

	failed |= check_i32("parser phase ./. count", count_dot_dot, 1);
	failed |= check_i32("parser phase 0/. count", count_0_dot, 1);
	failed |= check_i32("parser phase ./1 count", count_dot_1, 1);
	failed |= check_i32("parser phase 0/1 count", count_0_1, 1);

	hk_map_destroy(m);
	return failed != 0;
}
