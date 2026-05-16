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

static int check_greater(const char *label, float a, float b)
{
	if (!(a > b)) {
		fprintf(stderr, "%s: expected %.8g > %.8g\n", label, a, b);
		return 1;
	}
	return 0;
}

static int check_finite_outputs(const float energy[HK_BLIND_N_STATE], const float p4[HK_BLIND_N_STATE],
								float entropy, float pmax, float margin, float rho_output, float pU)
{
	int failed = 0, i;
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		if (!isfinite(energy[i])) {
			fprintf(stderr, "energy[%d] is not finite: %.8g\n", i, energy[i]);
			failed = 1;
		}
		if (!isfinite(p4[i])) {
			fprintf(stderr, "p4[%d] is not finite: %.8g\n", i, p4[i]);
			failed = 1;
		}
	}
	if (!isfinite(entropy)) {
		fprintf(stderr, "entropy is not finite: %.8g\n", entropy);
		failed = 1;
	}
	if (!isfinite(pmax)) {
		fprintf(stderr, "pmax is not finite: %.8g\n", pmax);
		failed = 1;
	}
	if (!isfinite(margin)) {
		fprintf(stderr, "margin is not finite: %.8g\n", margin);
		failed = 1;
	}
	if (!isfinite(rho_output)) {
		fprintf(stderr, "rho_output is not finite: %.8g\n", rho_output);
		failed = 1;
	}
	if (!isfinite(pU)) {
		fprintf(stderr, "pU is not finite: %.8g\n", pU);
		failed = 1;
	}
	return failed;
}

static int check_sums(const char *label, const float p4[HK_BLIND_N_STATE], float rho_output, float pU)
{
	float p_sum = 0.0f, out_sum;
	int i, failed = 0;
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		p_sum += p4[i];
	out_sum = rho_output * p_sum + pU;
	failed |= check_close(label, p_sum, 1.0f);
	failed |= check_close("output-with-U sum", out_sum, 1.0f);
	return failed;
}

static void set_coord(fvec3_t x, float x0, float x1, float x2)
{
	x[0] = x0;
	x[1] = x1;
	x[2] = x2;
}

static void run_coord_posterior(const struct hk_fdg_conf *conf, const fvec3_t *coords,
								float energy[HK_BLIND_N_STATE], float p4[HK_BLIND_N_STATE],
								float *entropy, float *pmax, float *margin, float *rho_output, float *pU)
{
	float log_prior[HK_BLIND_N_STATE];
	hk_blind_init_uniform_log_prior(log_prior);
	hk_blind_bpair_posterior_from_coords(conf, 0, 1, coords, 1.0f, 1.0f, 2.0f, log_prior, 1.0f,
										 energy, p4, entropy, pmax, margin, rho_output, pU);
}

static int check_equal_distances(const struct hk_fdg_conf *conf)
{
	fvec3_t coords[4];
	float energy[HK_BLIND_N_STATE], p4[HK_BLIND_N_STATE];
	float entropy, pmax, margin, rho_output, pU;
	int failed = 0, i;

	set_coord(coords[hk_diploid_bid(0, HK_DIPLOID_COPY0)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(0, HK_DIPLOID_COPY1)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(1, HK_DIPLOID_COPY0)], 1.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(1, HK_DIPLOID_COPY1)], 1.0f, 0.0f, 0.0f);

	run_coord_posterior(conf, coords, energy, p4, &entropy, &pmax, &margin, &rho_output, &pU);
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
		failed |= check_close("equal energy", energy[i], hk_fdg_contact_energy_dist(conf, 1.0f, 1.0f, 1.0f, 2.0f));
		failed |= check_close("equal p4", p4[i], 0.25f);
	}
	failed |= check_close("equal entropy", entropy, logf(4.0f));
	failed |= check_close("equal pmax", pmax, 0.25f);
	failed |= check_close("equal margin", margin, 0.0f);
	failed |= check_close("equal rho_output", rho_output, 0.0f);
	failed |= check_close("equal pU", pU, 1.0f);
	failed |= check_sums("equal prob sum", p4, rho_output, pU);
	return failed;
}

static int check_state_00_best(const struct hk_fdg_conf *conf)
{
	fvec3_t coords[4];
	float energy[HK_BLIND_N_STATE], p4[HK_BLIND_N_STATE];
	float entropy, pmax, margin, rho_output, pU;
	int failed = 0;

	set_coord(coords[hk_diploid_bid(0, HK_DIPLOID_COPY0)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(0, HK_DIPLOID_COPY1)], 100.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(1, HK_DIPLOID_COPY0)], 1.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(1, HK_DIPLOID_COPY1)], 200.0f, 0.0f, 0.0f);

	run_coord_posterior(conf, coords, energy, p4, &entropy, &pmax, &margin, &rho_output, &pU);
	failed |= check_close("00 energy", energy[HK_BLIND_STATE_00], hk_fdg_contact_energy_dist(conf, 1.0f, 1.0f, 1.0f, 2.0f));
	failed |= check_greater("00 p > 01", p4[HK_BLIND_STATE_00], p4[HK_BLIND_STATE_01]);
	failed |= check_greater("00 p > 10", p4[HK_BLIND_STATE_00], p4[HK_BLIND_STATE_10]);
	failed |= check_greater("00 p > 11", p4[HK_BLIND_STATE_00], p4[HK_BLIND_STATE_11]);
	failed |= check_close("00 pmax", pmax, p4[HK_BLIND_STATE_00]);
	failed |= check_finite_outputs(energy, p4, entropy, pmax, margin, rho_output, pU);
	failed |= check_sums("00 prob sum", p4, rho_output, pU);
	return failed;
}

static int check_state_11_best(const struct hk_fdg_conf *conf)
{
	fvec3_t coords[4];
	float energy[HK_BLIND_N_STATE], p4[HK_BLIND_N_STATE];
	float entropy, pmax, margin, rho_output, pU;
	int failed = 0;

	set_coord(coords[hk_diploid_bid(0, HK_DIPLOID_COPY0)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(0, HK_DIPLOID_COPY1)], 100.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(1, HK_DIPLOID_COPY0)], 200.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(1, HK_DIPLOID_COPY1)], 101.0f, 0.0f, 0.0f);

	run_coord_posterior(conf, coords, energy, p4, &entropy, &pmax, &margin, &rho_output, &pU);
	failed |= check_close("11 energy", energy[HK_BLIND_STATE_11], hk_fdg_contact_energy_dist(conf, 1.0f, 1.0f, 1.0f, 2.0f));
	failed |= check_greater("11 p > 00", p4[HK_BLIND_STATE_11], p4[HK_BLIND_STATE_00]);
	failed |= check_greater("11 p > 01", p4[HK_BLIND_STATE_11], p4[HK_BLIND_STATE_01]);
	failed |= check_greater("11 p > 10", p4[HK_BLIND_STATE_11], p4[HK_BLIND_STATE_10]);
	failed |= check_close("11 pmax", pmax, p4[HK_BLIND_STATE_11]);
	failed |= check_finite_outputs(energy, p4, entropy, pmax, margin, rho_output, pU);
	failed |= check_sums("11 prob sum", p4, rho_output, pU);
	return failed;
}

static int check_state_order(const struct hk_fdg_conf *conf)
{
	fvec3_t coords[4];
	float energy[HK_BLIND_N_STATE], p4[HK_BLIND_N_STATE];
	float entropy, pmax, margin, rho_output, pU;
	int failed = 0;

	set_coord(coords[hk_diploid_bid(0, HK_DIPLOID_COPY0)], 0.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(0, HK_DIPLOID_COPY1)], 1.0f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(1, HK_DIPLOID_COPY0)], -2.5f, 0.0f, 0.0f);
	set_coord(coords[hk_diploid_bid(1, HK_DIPLOID_COPY1)], -3.0f, 0.0f, 0.0f);

	run_coord_posterior(conf, coords, energy, p4, &entropy, &pmax, &margin, &rho_output, &pU);
	failed |= check_close("ordered E00", energy[HK_BLIND_STATE_00], hk_fdg_contact_energy_dist(conf, 2.5f, 1.0f, 1.0f, 2.0f));
	failed |= check_close("ordered E01", energy[HK_BLIND_STATE_01], hk_fdg_contact_energy_dist(conf, 3.0f, 1.0f, 1.0f, 2.0f));
	failed |= check_close("ordered E10", energy[HK_BLIND_STATE_10], hk_fdg_contact_energy_dist(conf, 3.5f, 1.0f, 1.0f, 2.0f));
	failed |= check_close("ordered E11", energy[HK_BLIND_STATE_11], hk_fdg_contact_energy_dist(conf, 4.0f, 1.0f, 1.0f, 2.0f));
	failed |= check_greater("ordered p00 > p01", p4[HK_BLIND_STATE_00], p4[HK_BLIND_STATE_01]);
	failed |= check_greater("ordered p01 > p10", p4[HK_BLIND_STATE_01], p4[HK_BLIND_STATE_10]);
	failed |= check_greater("ordered p10 > p11", p4[HK_BLIND_STATE_10], p4[HK_BLIND_STATE_11]);
	failed |= check_finite_outputs(energy, p4, entropy, pmax, margin, rho_output, pU);
	failed |= check_sums("ordered prob sum", p4, rho_output, pU);
	return failed;
}

int main(void)
{
	struct hk_fdg_conf conf;
	int failed = 0;

	hk_fdg_conf_init(&conf);
	failed |= check_equal_distances(&conf);
	failed |= check_state_00_best(&conf);
	failed |= check_state_11_best(&conf);
	failed |= check_state_order(&conf);
	return failed != 0;
}
