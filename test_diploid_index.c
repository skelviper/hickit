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

int main(int argc, char *argv[])
{
	int failed = 0;
	int32_t haploid_ids[] = {0, 1, 2, 7, 31, 1024};
	int32_t i, copy;

	if (argc == 2 && strcmp(argv[1], "--invalid-copy") == 0) {
		(void)hk_diploid_bid(0, 2);
		return 0;
	}

	failed |= check_i32("haploid=0 copy=0", hk_diploid_bid(0, 0), 0);
	failed |= check_i32("haploid=0 copy=1", hk_diploid_bid(0, 1), 1);
	failed |= check_i32("haploid=7 copy=0", hk_diploid_bid(7, 0), 14);
	failed |= check_i32("haploid=7 copy=1", hk_diploid_bid(7, 1), 15);
	failed |= check_i32("diploid=14 haploid", hk_diploid_haploid_bid(14), 7);
	failed |= check_i32("diploid=14 copy", hk_diploid_copy(14), 0);
	failed |= check_i32("diploid=15 haploid", hk_diploid_haploid_bid(15), 7);
	failed |= check_i32("diploid=15 copy", hk_diploid_copy(15), 1);

	for (i = 0; i < (int32_t)(sizeof(haploid_ids) / sizeof(haploid_ids[0])); ++i) {
		for (copy = 0; copy < HK_DIPLOID_N_COPY; ++copy) {
			int32_t haploid = haploid_ids[i];
			int32_t diploid = hk_diploid_bid(haploid, copy);
			failed |= check_i32("roundtrip haploid", hk_diploid_haploid_bid(diploid), haploid);
			failed |= check_i32("roundtrip copy", hk_diploid_copy(diploid), copy);
		}
	}

	return failed != 0;
}
