#include <math.h>
#include <stdio.h>
#include "hkpriv.h"

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

static float square(float x)
{
	return x * x;
}

int main(void)
{
	struct hk_fdg_conf conf;
	const float k = 2.0f;
	const float unit = 0.75f;
	const float d_scale = 1.25f;
	int failed = 0;
	int i;
	struct {
		const char *label;
		float r;
		float expected;
	} cases[7];

	hk_fdg_conf_init(&conf);
	cases[0].label = "below d_c1";
	cases[0].r = 0.0f;
	cases[0].expected = k * square(conf.d_c1 - cases[0].r);
	cases[1].label = "at d_c1";
	cases[1].r = conf.d_c1;
	cases[1].expected = 0.0f;
	cases[2].label = "between d_c1 and d_c2";
	cases[2].r = 0.5f * (conf.d_c1 + conf.d_c2);
	cases[2].expected = 0.0f;
	cases[3].label = "at d_c2";
	cases[3].r = conf.d_c2;
	cases[3].expected = 0.0f;
	cases[4].label = "between d_c2 and d_c3";
	cases[4].r = 0.5f * (conf.d_c2 + conf.d_c3);
	cases[4].expected = k * square(cases[4].r - conf.d_c2);
	cases[5].label = "at d_c3";
	cases[5].r = conf.d_c3;
	cases[5].expected = k * square(conf.d_c3 - conf.d_c2);
	cases[6].label = "beyond d_c3";
	cases[6].r = conf.d_c3 + 0.75f;
	cases[6].expected = k * (conf.c_c1 * (cases[6].r - conf.d_c3) + conf.c_c2 / (cases[6].r - conf.d_c2));

	for (i = 0; i < 7; ++i) {
		float got_r = hk_fdg_contact_energy_r(&conf, cases[i].r, k);
		float distance = cases[i].r * unit * d_scale;
		float got_dist = hk_fdg_contact_energy_dist(&conf, distance, unit, d_scale, k);
		failed |= check_close(cases[i].label, got_r, cases[i].expected);
		failed |= check_close("distance wrapper", got_dist, got_r);
	}
	return failed != 0;
}
