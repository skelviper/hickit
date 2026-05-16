#include <math.h>
#include <stdio.h>
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

static void set_coord(fvec3_t x, float x0, float x1, float x2)
{
	x[0] = x0;
	x[1] = x1;
	x[2] = x2;
}

static int check_vec_zero(const char *label, const fvec3_t x)
{
	int failed = 0;
	char buf[64];
	snprintf(buf, sizeof(buf), "%s x", label);
	failed |= check_close(buf, x[0], 0.0f);
	snprintf(buf, sizeof(buf), "%s y", label);
	failed |= check_close(buf, x[1], 0.0f);
	snprintf(buf, sizeof(buf), "%s z", label);
	failed |= check_close(buf, x[2], 0.0f);
	return failed;
}

static int check_force_finite_pair(const fvec3_t f0, const fvec3_t f1)
{
	int failed = 0, a;
	for (a = 0; a < 3; ++a) {
		failed |= check_true("finite f0", isfinite(f0[a]));
		failed |= check_true("finite f1", isfinite(f1[a]));
		failed |= check_close("opposite forces", f0[a] + f1[a], 0.0f);
	}
	return failed;
}

static int check_beyond_min_sep(void)
{
	fvec3_t x0, x1, f0, f1;
	float e, ef;
	int failed = 0;

	set_coord(x0, 0.0f, 0.0f, 0.0f);
	set_coord(x1, 1.0f, 0.0f, 0.0f);
	e = hk_blind_homolog_sep_energy(x0, x1, 1.0f, 0.5f, 0.05f);
	ef = hk_blind_homolog_sep_energy_force(x0, x1, 1.0f, 0.5f, 0.05f, f0, f1);
	failed |= check_true("beyond finite energy", isfinite(e));
	failed |= check_true("beyond finite energy force", isfinite(ef));
	failed |= check_close("beyond energy", e, 0.0f);
	failed |= check_close("beyond force energy", ef, e);
	failed |= check_vec_zero("beyond f0", f0);
	failed |= check_vec_zero("beyond f1", f1);
	return failed;
}

static int check_below_min_sep(void)
{
	fvec3_t x0, x1, f0, f1;
	float lambda_sep = 0.05f, min_sep_unit = 0.5f;
	float expected_e = lambda_sep * 0.25f * 0.25f;
	float expected_f = 2.0f * lambda_sep * 0.25f;
	float e, ef;
	int failed = 0;

	set_coord(x0, 0.0f, 0.0f, 0.0f);
	set_coord(x1, 0.25f, 0.0f, 0.0f);
	e = hk_blind_homolog_sep_energy(x0, x1, 1.0f, min_sep_unit, lambda_sep);
	ef = hk_blind_homolog_sep_energy_force(x0, x1, 1.0f, min_sep_unit, lambda_sep, f0, f1);
	failed |= check_true("below finite energy", isfinite(e));
	failed |= check_true("below finite energy force", isfinite(ef));
	failed |= check_close("below energy", e, expected_e);
	failed |= check_close("below force energy", ef, expected_e);
	failed |= check_force_finite_pair(f0, f1);
	failed |= check_close("below f0 x", f0[0], -expected_f);
	failed |= check_close("below f1 x", f1[0], expected_f);
	failed |= check_close("below f0 y", f0[1], 0.0f);
	failed |= check_close("below f0 z", f0[2], 0.0f);
	failed |= check_true("below f0 pushes negative x", f0[0] < 0.0f);
	failed |= check_true("below f1 pushes positive x", f1[0] > 0.0f);
	return failed;
}

static int check_at_boundary(void)
{
	fvec3_t x0, x1, f0, f1;
	float e, ef;
	int failed = 0;

	set_coord(x0, 0.0f, 0.0f, 0.0f);
	set_coord(x1, 0.5f, 0.0f, 0.0f);
	e = hk_blind_homolog_sep_energy(x0, x1, 1.0f, 0.5f, 0.05f);
	ef = hk_blind_homolog_sep_energy_force(x0, x1, 1.0f, 0.5f, 0.05f, f0, f1);
	failed |= check_true("boundary finite energy", isfinite(e));
	failed |= check_true("boundary finite force energy", isfinite(ef));
	failed |= check_close("boundary energy", e, 0.0f);
	failed |= check_close("boundary force energy", ef, e);
	failed |= check_vec_zero("boundary f0", f0);
	failed |= check_vec_zero("boundary f1", f1);
	return failed;
}

static int check_zero_distance(void)
{
	fvec3_t x0, x1, f0, f1;
	float lambda_sep = 0.05f, min_sep_unit = 0.5f;
	float expected_e = lambda_sep * min_sep_unit * min_sep_unit;
	float expected_f = 2.0f * lambda_sep * min_sep_unit;
	float e, ef;
	int failed = 0;

	set_coord(x0, 0.0f, 0.0f, 0.0f);
	set_coord(x1, 0.0f, 0.0f, 0.0f);
	e = hk_blind_homolog_sep_energy(x0, x1, 1.0f, min_sep_unit, lambda_sep);
	ef = hk_blind_homolog_sep_energy_force(x0, x1, 1.0f, min_sep_unit, lambda_sep, f0, f1);
	failed |= check_true("zero finite energy", isfinite(e));
	failed |= check_true("zero finite force energy", isfinite(ef));
	failed |= check_close("zero energy", e, expected_e);
	failed |= check_close("zero force energy", ef, expected_e);
	failed |= check_force_finite_pair(f0, f1);
	failed |= check_close("zero f0 x", f0[0], expected_f);
	failed |= check_close("zero f1 x", f1[0], -expected_f);
	failed |= check_close("zero f0 y", f0[1], 0.0f);
	failed |= check_close("zero f0 z", f0[2], 0.0f);
	failed |= check_close("zero f1 y", f1[1], 0.0f);
	failed |= check_close("zero f1 z", f1[2], 0.0f);
	return failed;
}

static int check_lambda_zero(void)
{
	fvec3_t x0, x1, f0, f1;
	float e, ef;
	int failed = 0;

	set_coord(x0, 0.0f, 0.0f, 0.0f);
	set_coord(x1, 0.25f, 0.0f, 0.0f);
	e = hk_blind_homolog_sep_energy(x0, x1, 1.0f, 0.5f, 0.0f);
	ef = hk_blind_homolog_sep_energy_force(x0, x1, 1.0f, 0.5f, 0.0f, f0, f1);
	failed |= check_true("lambda0 finite energy", isfinite(e));
	failed |= check_true("lambda0 finite force energy", isfinite(ef));
	failed |= check_close("lambda0 energy", e, 0.0f);
	failed |= check_close("lambda0 force energy", ef, 0.0f);
	failed |= check_vec_zero("lambda0 f0", f0);
	failed |= check_vec_zero("lambda0 f1", f1);
	return failed;
}

static int check_unit_scaling(void)
{
	fvec3_t x0, x1, f0, f1;
	float lambda_sep = 0.05f;
	float expected_e = lambda_sep * 0.25f * 0.25f;
	float e, ef;
	int failed = 0;

	set_coord(x0, 0.0f, 0.0f, 0.0f);
	set_coord(x1, 0.5f, 0.0f, 0.0f);
	e = hk_blind_homolog_sep_energy(x0, x1, 2.0f, 0.5f, lambda_sep);
	ef = hk_blind_homolog_sep_energy_force(x0, x1, 2.0f, 0.5f, lambda_sep, f0, f1);
	failed |= check_true("unit finite energy", isfinite(e));
	failed |= check_true("unit finite force energy", isfinite(ef));
	failed |= check_close("unit-scaled energy", e, expected_e);
	failed |= check_close("unit-scaled force energy", ef, expected_e);
	failed |= check_force_finite_pair(f0, f1);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= check_beyond_min_sep();
	failed |= check_below_min_sep();
	failed |= check_at_boundary();
	failed |= check_zero_distance();
	failed |= check_lambda_zero();
	failed |= check_unit_scaling();
	return failed != 0;
}
