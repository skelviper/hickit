#include <math.h>
#include <stdio.h>
#include "hickit.h"

static int check_i32(const char *label, int got, int expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %d, expected %d\n", label, got, expected);
		return 1;
	}
	return 0;
}

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

static int check_less(const char *label, float a, float b)
{
	if (!(a < b)) {
		fprintf(stderr, "%s: expected %.8g < %.8g\n", label, a, b);
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

static int check_finite_posterior(const float p4[HK_BLIND_N_STATE], float entropy, float pmax, float margin, float rho_output, float pU)
{
	int failed = 0, i;
	for (i = 0; i < HK_BLIND_N_STATE; ++i) {
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

static void run_posterior(const float energy[HK_BLIND_N_STATE], const float log_prior[HK_BLIND_N_STATE], float temperature,
						  float p4[HK_BLIND_N_STATE], float *entropy, float *pmax, float *margin, float *rho_output, float *pU)
{
	hk_blind_posterior_from_energy(energy, log_prior, temperature, p4, entropy, pmax, margin, rho_output, pU);
}

static int check_uniform_energies(void)
{
	float energy[HK_BLIND_N_STATE] = {0.0f, 0.0f, 0.0f, 0.0f};
	float log_prior[HK_BLIND_N_STATE], p4[HK_BLIND_N_STATE];
	float entropy, pmax, margin, rho_output, pU;
	int failed = 0, i;

	hk_blind_init_uniform_log_prior(log_prior);
	run_posterior(energy, log_prior, 1.0f, p4, &entropy, &pmax, &margin, &rho_output, &pU);
	for (i = 0; i < HK_BLIND_N_STATE; ++i)
		failed |= check_close("uniform p4", p4[i], 0.25f);
	failed |= check_close("uniform entropy", entropy, logf(4.0f));
	failed |= check_close("uniform pmax", pmax, 0.25f);
	failed |= check_close("uniform margin", margin, 0.0f);
	failed |= check_close("uniform rho_output", rho_output, 0.0f);
	failed |= check_close("uniform pU", pU, 1.0f);
	failed |= check_sums("uniform prob sum", p4, rho_output, pU);
	return failed;
}

static int check_best_state(void)
{
	float energy[HK_BLIND_N_STATE] = {0.0f, 10.0f, 10.0f, 10.0f};
	float log_prior[HK_BLIND_N_STATE], p4[HK_BLIND_N_STATE];
	float entropy, pmax, margin, rho_output, pU;
	int failed = 0;

	hk_blind_init_uniform_log_prior(log_prior);
	run_posterior(energy, log_prior, 1.0f, p4, &entropy, &pmax, &margin, &rho_output, &pU);
	failed |= check_greater("best p00 > p01", p4[HK_BLIND_STATE_00], p4[HK_BLIND_STATE_01]);
	failed |= check_greater("best p00 > p10", p4[HK_BLIND_STATE_00], p4[HK_BLIND_STATE_10]);
	failed |= check_greater("best p00 > p11", p4[HK_BLIND_STATE_00], p4[HK_BLIND_STATE_11]);
	failed |= check_close("best pmax", pmax, p4[HK_BLIND_STATE_00]);
	failed |= check_less("best entropy", entropy, logf(4.0f));
	failed |= check_less("best pU", pU, 1.0f);
	failed |= check_greater("best rho_output", rho_output, 0.0f);
	failed |= check_greater("best margin", margin, 0.0f);
	failed |= check_sums("best prob sum", p4, rho_output, pU);
	return failed;
}

static int check_temperature_effect(void)
{
	float energy[HK_BLIND_N_STATE] = {0.0f, 2.0f, 2.0f, 2.0f};
	float log_prior[HK_BLIND_N_STATE], p4_t1[HK_BLIND_N_STATE], p4_t2[HK_BLIND_N_STATE];
	float entropy_t1, pmax_t1, margin_t1, rho_t1, pU_t1;
	float entropy_t2, pmax_t2, margin_t2, rho_t2, pU_t2;
	int failed = 0;

	hk_blind_init_uniform_log_prior(log_prior);
	run_posterior(energy, log_prior, 1.0f, p4_t1, &entropy_t1, &pmax_t1, &margin_t1, &rho_t1, &pU_t1);
	run_posterior(energy, log_prior, 2.0f, p4_t2, &entropy_t2, &pmax_t2, &margin_t2, &rho_t2, &pU_t2);
	failed |= check_less("T=2 p00 less sharp", p4_t2[HK_BLIND_STATE_00], p4_t1[HK_BLIND_STATE_00]);
	failed |= check_greater("T=2 entropy higher", entropy_t2, entropy_t1);
	failed |= check_greater("T=2 pU higher", pU_t2, pU_t1);
	failed |= check_sums("T=1 prob sum", p4_t1, rho_t1, pU_t1);
	failed |= check_sums("T=2 prob sum", p4_t2, rho_t2, pU_t2);
	return failed;
}

static int check_stable_softmax(void)
{
	float energy[HK_BLIND_N_STATE] = {10000.0f, 10001.0f, 10002.0f, 10003.0f};
	float log_prior[HK_BLIND_N_STATE], p4[HK_BLIND_N_STATE];
	float entropy, pmax, margin, rho_output, pU;
	int failed = 0;

	hk_blind_init_uniform_log_prior(log_prior);
	run_posterior(energy, log_prior, 1.0f, p4, &entropy, &pmax, &margin, &rho_output, &pU);
	failed |= check_finite_posterior(p4, entropy, pmax, margin, rho_output, pU);
	failed |= check_greater("stable p00 > p01", p4[HK_BLIND_STATE_00], p4[HK_BLIND_STATE_01]);
	failed |= check_greater("stable p01 > p10", p4[HK_BLIND_STATE_01], p4[HK_BLIND_STATE_10]);
	failed |= check_greater("stable p10 > p11", p4[HK_BLIND_STATE_10], p4[HK_BLIND_STATE_11]);
	failed |= check_sums("stable prob sum", p4, rho_output, pU);
	return failed;
}

static int check_nonuniform_prior(void)
{
	float energy[HK_BLIND_N_STATE] = {0.0f, 0.0f, 0.0f, 0.0f};
	float log_prior[HK_BLIND_N_STATE], p4[HK_BLIND_N_STATE];
	float entropy, pmax, margin, rho_output, pU;
	int failed = 0;

	log_prior[HK_BLIND_STATE_00] = logf(0.7f);
	log_prior[HK_BLIND_STATE_01] = logf(0.1f);
	log_prior[HK_BLIND_STATE_10] = logf(0.1f);
	log_prior[HK_BLIND_STATE_11] = logf(0.1f);
	run_posterior(energy, log_prior, 1.0f, p4, &entropy, &pmax, &margin, &rho_output, &pU);
	failed |= check_close("prior p00", p4[HK_BLIND_STATE_00], 0.7f);
	failed |= check_close("prior p01", p4[HK_BLIND_STATE_01], 0.1f);
	failed |= check_close("prior p10", p4[HK_BLIND_STATE_10], 0.1f);
	failed |= check_close("prior p11", p4[HK_BLIND_STATE_11], 0.1f);
	failed |= check_close("prior pmax", pmax, 0.7f);
	failed |= check_close("prior margin", margin, 0.6f);
	failed |= check_sums("prior prob sum", p4, rho_output, pU);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_i32("state 00", HK_BLIND_STATE_00, 0);
	failed |= check_i32("state 01", HK_BLIND_STATE_01, 1);
	failed |= check_i32("state 10", HK_BLIND_STATE_10, 2);
	failed |= check_i32("state 11", HK_BLIND_STATE_11, 3);
	failed |= check_i32("n state", HK_BLIND_N_STATE, 4);
	failed |= check_uniform_energies();
	failed |= check_best_state();
	failed |= check_temperature_effect();
	failed |= check_stable_softmax();
	failed |= check_nonuniform_prior();
	return failed != 0;
}
