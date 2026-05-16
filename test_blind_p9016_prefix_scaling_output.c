#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "hickit.h"

#define HK_BLIND_P9016_SCALING_PATH "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_SCALING_N_PREFIX 2
#define HK_BLIND_P9016_SCALING_N_ITER 3
#define HK_BLIND_P9016_SCALING_UNIT 1.0f
#define HK_BLIND_P9016_SCALING_D_SCALE 1.0f
#define HK_BLIND_P9016_SCALING_BASE_K 2.0f
#define HK_BLIND_P9016_SCALING_TEMPERATURE_START 2.0f
#define HK_BLIND_P9016_SCALING_TEMPERATURE_END 1.0f
#define HK_BLIND_P9016_SCALING_RHO_TRAIN_START 1.0f
#define HK_BLIND_P9016_SCALING_RHO_TRAIN_END 1.0f
#define HK_BLIND_P9016_SCALING_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_SCALING_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_SCALING_STEP 0.001f
#define HK_BLIND_P9016_SCALING_RELAX_STEPS 5

struct scaling_result {
	int32_t requested_n;
	int32_t n_raw;
	int32_t n_bpair;
	int32_t n_beads;
	float final_mean_entropy;
	float final_mean_pU;
	double final_sum_wedge_k;
	int moved;
};

static const int32_t hk_blind_scaling_prefix_sizes[HK_BLIND_P9016_SCALING_N_PREFIX] = { 1000, 10000 };

static int check_i32(const char *label, int32_t got, int32_t expected)
{
	if (got != expected) {
		fprintf(stderr, "%s: got %d, expected %d\n", label, got, expected);
		return 1;
	}
	return 0;
}

static int check_close(const char *label, double got, double expected)
{
	double tol = 1e-5;
	double scale = fabs(expected) > 1.0? fabs(expected) : 1.0;
	if (fabs(got - expected) > tol * scale) {
		fprintf(stderr, "%s: got %.12g, expected %.12g\n", label, got, expected);
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

static int file_exists(const char *path)
{
	FILE *fp = fopen(path, "rb");
	if (fp == 0) return 0;
	fclose(fp);
	return 1;
}

static char *read_file(const char *path)
{
	FILE *fp = fopen(path, "rb");
	long len;
	char *buf;
	size_t n_read;

	if (fp == 0) return 0;
	if (fseek(fp, 0, SEEK_END) != 0) {
		fclose(fp);
		return 0;
	}
	len = ftell(fp);
	if (len < 0) {
		fclose(fp);
		return 0;
	}
	if (fseek(fp, 0, SEEK_SET) != 0) {
		fclose(fp);
		return 0;
	}
	buf = (char*)malloc((size_t)len + 1);
	if (buf == 0) {
		fclose(fp);
		return 0;
	}
	n_read = fread(buf, 1, (size_t)len, fp);
	buf[n_read] = 0;
	fclose(fp);
	return buf;
}

static int count_lines(const char *s)
{
	int n = 0;
	for (; *s; ++s)
		if (*s == '\n') ++n;
	return n;
}

static int check_no_forbidden_strings(const char *label, const char *s)
{
	int failed = 0;
	char buf[128];
	const char *forbidden[] = { "phase0", "phase1", "truth", "oracle" };
	int i;
	for (i = 0; i < 4; ++i) {
		snprintf(buf, sizeof(buf), "%s no %s", label, forbidden[i]);
		failed |= check_true(buf, strstr(s, forbidden[i]) == 0);
	}
	return failed;
}

static int make_temp_dir(char *dir, size_t dir_size)
{
	int i;
	for (i = 0; i < 100; ++i) {
		snprintf(dir, dir_size, "/tmp/hk_blind_p9016_scaling_%ld_%d", (long)getpid(), i);
		if (mkdir(dir, 0700) == 0)
			return 0;
		if (errno != EEXIST)
			return -1;
	}
	return -1;
}

static void remove_run_outputs(const char *run_dir)
{
	char path[512];
	if (run_dir[0] == 0) return;
	snprintf(path, sizeof(path), "%s/blind_prefix.bpair_posterior.tsv", run_dir);
	remove(path);
	snprintf(path, sizeof(path), "%s/blind_prefix.raw_posterior.tsv", run_dir);
	remove(path);
	snprintf(path, sizeof(path), "%s/blind_prefix.coords.tsv", run_dir);
	remove(path);
	snprintf(path, sizeof(path), "%s/blind_prefix.loop_diag.tsv", run_dir);
	remove(path);
	rmdir(run_dir);
}

static void set_haploid_scaffold(const struct hk_bmap *bmap, fvec3_t *haploid)
{
	int32_t i;
	for (i = 0; i < bmap->n_beads; ++i) {
		haploid[i][0] = 0.10f * (float)(i % 97);
		haploid[i][1] = 0.07f * (float)((i / 97) % 97);
		haploid[i][2] = 0.03f * (float)(i % 17);
	}
}

static void copy_coords(fvec3_t *dst, const fvec3_t *src, int32_t n)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			dst[i][a] = src[i][a];
}

static int coords_any_changed(const fvec3_t *before, const fvec3_t *after, int32_t n, float tol)
{
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			if (fabsf(after[i][a] - before[i][a]) > tol)
				return 1;
	return 0;
}

static int check_coords_finite(const fvec3_t *coords, int32_t n)
{
	int failed = 0;
	int32_t i;
	int a;
	for (i = 0; i < n; ++i)
		for (a = 0; a < 3; ++a)
			failed |= check_true("scaling coord finite", isfinite(coords[i][a]));
	return failed;
}

static void set_schedule_conf(struct hk_blind_iter_schedule_conf *conf)
{
	conf->n_iter = HK_BLIND_P9016_SCALING_N_ITER;
	conf->base_conf.unit = HK_BLIND_P9016_SCALING_UNIT;
	conf->base_conf.d_scale = HK_BLIND_P9016_SCALING_D_SCALE;
	conf->base_conf.base_k = HK_BLIND_P9016_SCALING_BASE_K;
	conf->base_conf.temperature = HK_BLIND_P9016_SCALING_TEMPERATURE_START;
	conf->base_conf.rho_train = HK_BLIND_P9016_SCALING_RHO_TRAIN_START;
	conf->base_conf.min_sep_unit = HK_BLIND_P9016_SCALING_MIN_SEP_UNIT;
	conf->base_conf.lambda_sep = HK_BLIND_P9016_SCALING_LAMBDA_SEP;
	conf->base_conf.relax_step = HK_BLIND_P9016_SCALING_STEP;
	conf->base_conf.relax_steps = HK_BLIND_P9016_SCALING_RELAX_STEPS;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->base_conf.d_scale_eps_count = 1e-6f;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = HK_BLIND_P9016_SCALING_TEMPERATURE_START;
	conf->temperature_end = HK_BLIND_P9016_SCALING_TEMPERATURE_END;
	conf->rho_train_start = HK_BLIND_P9016_SCALING_RHO_TRAIN_START;
	conf->rho_train_end = HK_BLIND_P9016_SCALING_RHO_TRAIN_END;
}

static int check_binned_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set, int32_t max_raw)
{
	int failed = 0;
	int32_t i;

	failed |= check_true("scaling set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_true("scaling raw positive", set->n_raw > 0);
	failed |= check_true("scaling raw bounded", set->n_raw <= max_raw);
	failed |= check_true("scaling bpair positive", set->n_bpairs > 0);
	for (i = 0; i < set->n_raw; ++i) {
		failed |= check_true("scaling raw2binned id lower", set->raw2binned[i].bpair_id >= 0);
		failed |= check_true("scaling raw2binned id upper", set->raw2binned[i].bpair_id < set->n_bpairs);
		failed |= check_true("scaling raw2binned swapped", set->raw2binned[i].swapped == 0 || set->raw2binned[i].swapped == 1);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		failed |= check_true("scaling key sorted", bp->key.bid[0] <= bp->key.bid[1]);
		failed |= check_true("scaling key bid0 valid", bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		failed |= check_true("scaling key bid1 valid", bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
	}
	return failed;
}

static int check_loop_health(const struct hk_blind_iter_loop_diag *diag,
							 const struct hk_blind_single_iter_diag *per_iter)
{
	int failed = 0;
	int32_t i;

	failed |= check_i32("scaling loop n_iter", diag->n_iter, HK_BLIND_P9016_SCALING_N_ITER);
	failed |= check_i32("scaling loop completed", diag->n_completed, HK_BLIND_P9016_SCALING_N_ITER);
	failed |= check_i32("scaling loop bad iter", diag->n_bad_iter, 0);
	failed |= check_i32("scaling loop relax bad", diag->n_relax_nonfinite_iter, 0);
	failed |= check_i32("scaling loop coord bad", diag->n_coord_nonfinite, 0);
	failed |= check_close("scaling initial temperature", diag->initial_temperature, HK_BLIND_P9016_SCALING_TEMPERATURE_START);
	failed |= check_close("scaling final temperature", diag->final_temperature, HK_BLIND_P9016_SCALING_TEMPERATURE_END);
	failed |= check_close("scaling initial rho", diag->initial_rho_train, HK_BLIND_P9016_SCALING_RHO_TRAIN_START);
	failed |= check_close("scaling final rho", diag->final_rho_train, HK_BLIND_P9016_SCALING_RHO_TRAIN_END);
	failed |= check_true("scaling final entropy finite", isfinite(diag->final_mean_entropy));
	failed |= check_true("scaling final entropy range", diag->final_mean_entropy >= -1e-6f &&
						 diag->final_mean_entropy <= logf((float)HK_BLIND_N_STATE) + 1e-6f);
	failed |= check_true("scaling final pU finite", isfinite(diag->final_mean_pU));
	failed |= check_true("scaling final pU range", diag->final_mean_pU >= -1e-6f && diag->final_mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("scaling final sep finite", isfinite(diag->final_mean_sep));
	failed |= check_true("scaling final sum k finite", isfinite(diag->final_sum_wedge_k));
	failed |= check_true("scaling final sum k nonnegative", diag->final_sum_wedge_k >= 0.0);
	for (i = 0; i < HK_BLIND_P9016_SCALING_N_ITER; ++i) {
		const struct hk_blind_iter_diag *pre = &per_iter[i].pre_relax_diag;
		failed |= check_i32("scaling posterior nonfinite", pre->n_posterior_nonfinite, 0);
		failed |= check_i32("scaling posterior bad sum", pre->n_posterior_bad_sum, 0);
		failed |= check_i32("scaling posterior out of range", pre->n_posterior_out_of_range, 0);
		failed |= check_i32("scaling uncertainty nonfinite", pre->n_uncertainty_nonfinite, 0);
		failed |= check_i32("scaling uncertainty out of range", pre->n_uncertainty_out_of_range, 0);
		failed |= check_i32("scaling five-state bad sum", pre->n_five_state_bad_sum, 0);
		failed |= check_i32("scaling wedge nonfinite", pre->n_wedge_nonfinite, 0);
		failed |= check_i32("scaling wedge bad k", pre->n_wedge_bad_k, 0);
		failed |= check_i32("scaling wedge bad d_scale", pre->n_wedge_bad_d_scale, 0);
		failed |= check_i32("scaling sep force nonfinite", pre->sep_force_nonfinite, 0);
		failed |= check_i32("scaling relax completed", per_iter[i].relax_diag.n_completed,
							HK_BLIND_P9016_SCALING_RELAX_STEPS);
		failed |= check_i32("scaling relax nonfinite step", per_iter[i].relax_diag.n_nonfinite_step, 0);
		failed |= check_i32("scaling relax coord bad", per_iter[i].relax_diag.n_coord_nonfinite, 0);
		failed |= check_i32("scaling backbone bad", per_iter[i].relax_diag.n_backbone_nonfinite_step, 0);
		failed |= check_i32("scaling repulsion bad", per_iter[i].relax_diag.n_repulsion_nonfinite_step, 0);
		failed |= check_true("scaling relax total finite", isfinite(per_iter[i].relax_diag.final_total_energy));
		failed |= check_true("scaling repulsion finite", isfinite(per_iter[i].relax_diag.final_repulsion_energy));
	}
	return failed;
}

static int write_outputs(const char *posterior_path, const char *raw_path,
						 const char *coords_path, const char *diag_path,
						 const struct hk_blind_pair *raw, int32_t n_raw,
						 const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set,
						 const fvec3_t *coords, const struct hk_blind_iter_loop_diag *loop_diag)
{
	FILE *fp;
	int ret;

	fp = fopen(posterior_path, "w");
	if (fp == 0) return -1;
	ret = hk_blind_write_bpair_posterior_tsv(fp, bmap, set);
	if (fclose(fp) != 0) ret = -1;
	if (ret != 0) return -1;

	fp = fopen(raw_path, "w");
	if (fp == 0) return -1;
	ret = hk_blind_write_raw_contact_posterior_tsv(fp, raw, n_raw, bmap, set);
	if (fclose(fp) != 0) ret = -1;
	if (ret != 0) return -1;

	fp = fopen(coords_path, "w");
	if (fp == 0) return -1;
	ret = hk_blind_write_diploid_coords_tsv(fp, bmap, coords);
	if (fclose(fp) != 0) ret = -1;
	if (ret != 0) return -1;

	fp = fopen(diag_path, "w");
	if (fp == 0) return -1;
	ret = hk_blind_write_iter_loop_diag_tsv(fp, loop_diag);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}

static int next_line(char **p, char *line, size_t line_size)
{
	char *start = *p;
	char *end;
	size_t len;

	if (start == 0 || *start == 0) return 0;
	end = strchr(start, '\n');
	len = end? (size_t)(end - start) : strlen(start);
	if (len >= line_size) len = line_size - 1;
	memcpy(line, start, len);
	line[len] = 0;
	*p = end? end + 1 : start + strlen(start);
	return 1;
}

static int validate_binned_posterior_file(const char *path, const struct hk_blind_bpair_set *set)
{
	char *text = read_file(path);
	char *p;
	char line[4096];
	int failed = 0;
	int row = 0;

	failed |= check_true("scaling binned posterior exists", file_exists(path));
	failed |= check_true("scaling binned posterior text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("scaling binned header p00", strstr(text, "p00") != 0);
	failed |= check_true("scaling binned header pU", strstr(text, "pU") != 0);
	failed |= check_true("scaling binned header base_d_scale", strstr(text, "base_d_scale") != 0);
	failed |= check_true("scaling binned header base_k", strstr(text, "base_k") != 0);
	failed |= check_true("scaling binned header rho", strstr(text, "rho_output") != 0);
	failed |= check_no_forbidden_strings("scaling binned posterior", text);
	failed |= check_i32("scaling binned row count", count_lines(text) - 1, set->n_bpairs);
	p = text;
	next_line(&p, line, sizeof(line)); // header
	while (next_line(&p, line, sizeof(line))) {
		char chr1[128], chr2[128];
		int st1, en1, st2, en2, bid1, bid2, n_raw;
		double base_d_scale, base_k;
		double p00, p01, p10, p11, pU, psame, pcross, entropy, margin, pmax, rho;
		failed |= check_i32("scaling binned parsed",
							sscanf(line, "%127s\t%d\t%d\t%127s\t%d\t%d\t%d\t%d\t%d\t"
								   "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
								   chr1, &st1, &en1, chr2, &st2, &en2, &bid1, &bid2, &n_raw,
								   &base_d_scale, &base_k, &p00, &p01, &p10, &p11, &pU, &psame, &pcross,
								   &entropy, &margin, &pmax, &rho), 22);
		failed |= check_true("scaling binned chr nonempty", chr1[0] != 0 && chr2[0] != 0);
		failed |= check_true("scaling binned interval", st1 <= en1 && st2 <= en2);
		failed |= check_true("scaling binned bids", bid1 >= 0 && bid2 >= bid1);
		failed |= check_true("scaling binned n_raw", n_raw > 0);
		failed |= check_close("scaling binned base_d_scale", base_d_scale, pow((double)n_raw, -1.0 / 3.0));
		failed |= check_close("scaling binned base_k", base_k, 1.0);
		failed |= check_true("scaling binned finite", isfinite(p00) && isfinite(p01) &&
							 isfinite(p10) && isfinite(p11) && isfinite(pU) &&
							 isfinite(entropy) && isfinite(margin) && isfinite(pmax) &&
							 isfinite(rho));
		failed |= check_close("scaling binned p4 sum", p00 + p01 + p10 + p11, 1.0);
		failed |= check_close("scaling binned five-state sum", rho * (p00 + p01 + p10 + p11) + pU, 1.0);
		failed |= check_close("scaling binned psame", psame, p00 + p11);
		failed |= check_close("scaling binned pcross", pcross, p01 + p10);
		++row;
	}
	failed |= check_i32("scaling binned row loop", row, set->n_bpairs);
	free(text);
	return failed;
}

static int validate_raw_posterior_file(const char *path, const struct hk_blind_bpair_set *set)
{
	char *text = read_file(path);
	char *p;
	char line[4096];
	int failed = 0;
	int row = 0;

	failed |= check_true("scaling raw posterior exists", file_exists(path));
	failed |= check_true("scaling raw posterior text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("scaling raw header raw_id", strstr(text, "raw_id") != 0);
	failed |= check_true("scaling raw header swapped", strstr(text, "swapped") != 0);
	failed |= check_true("scaling raw header p00", strstr(text, "p00") != 0);
	failed |= check_true("scaling raw header pU", strstr(text, "pU") != 0);
	failed |= check_no_forbidden_strings("scaling raw posterior", text);
	failed |= check_i32("scaling raw row count", count_lines(text) - 1, set->n_raw);
	p = text;
	next_line(&p, line, sizeof(line)); // header
	while (next_line(&p, line, sizeof(line))) {
		char chr1[128], chr2[128];
		int raw_id, pos1, pos2, bid1_raw, bid2_raw, bpair_id, bid1_can, bid2_can, swapped;
		double p00, p01, p10, p11, pU, psame, pcross, entropy, margin, pmax, rho;
		failed |= check_i32("scaling raw parsed",
							sscanf(line, "%d\t%127s\t%d\t%127s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t"
								   "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
								   &raw_id, chr1, &pos1, chr2, &pos2, &bid1_raw, &bid2_raw,
								   &bpair_id, &bid1_can, &bid2_can, &swapped,
								   &p00, &p01, &p10, &p11, &pU, &psame, &pcross,
								   &entropy, &margin, &pmax, &rho), 22);
		failed |= check_i32("scaling raw id", raw_id, row);
		failed |= check_true("scaling raw chr nonempty", chr1[0] != 0 && chr2[0] != 0);
		failed |= check_true("scaling raw positions", pos1 >= 0 && pos2 >= 0);
		failed |= check_true("scaling raw bids", bid1_raw >= 0 && bid2_raw >= 0);
		failed |= check_true("scaling raw bpair id", bpair_id >= 0 && bpair_id < set->n_bpairs);
		failed |= check_true("scaling raw canonical", bid1_can <= bid2_can);
		failed |= check_true("scaling raw swapped", swapped == 0 || swapped == 1);
		failed |= check_true("scaling raw finite", isfinite(p00) && isfinite(p01) &&
							 isfinite(p10) && isfinite(p11) && isfinite(pU) &&
							 isfinite(entropy) && isfinite(margin) && isfinite(pmax) &&
							 isfinite(rho));
		failed |= check_close("scaling raw p4 sum", p00 + p01 + p10 + p11, 1.0);
		failed |= check_close("scaling raw five-state sum", rho * (p00 + p01 + p10 + p11) + pU, 1.0);
		failed |= check_close("scaling raw psame", psame, p00 + p11);
		failed |= check_close("scaling raw pcross", pcross, p01 + p10);
		++row;
	}
	failed |= check_i32("scaling raw row loop", row, set->n_raw);
	free(text);
	return failed;
}

static int validate_coords_file(const char *path, const struct hk_bmap *bmap)
{
	char *text = read_file(path);
	char *p;
	char line[2048];
	int failed = 0;
	int row = 0;

	failed |= check_true("scaling coords exists", file_exists(path));
	failed |= check_true("scaling coords text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("scaling coords header", strstr(text, "diploid_bid") != 0);
	failed |= check_no_forbidden_strings("scaling coords", text);
	failed |= check_i32("scaling coords row count", count_lines(text) - 1, 2 * bmap->n_beads);
	p = text;
	next_line(&p, line, sizeof(line)); // header
	while (next_line(&p, line, sizeof(line))) {
		char chr[128];
		int st, en, bid, copy, diploid_bid;
		double x, y, z;
		failed |= check_i32("scaling coords parsed",
							sscanf(line, "%127s\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf",
								   chr, &st, &en, &bid, &copy, &diploid_bid, &x, &y, &z), 9);
		failed |= check_true("scaling coords chr", chr[0] != 0);
		failed |= check_true("scaling coords interval", st <= en);
		failed |= check_i32("scaling coords bid", bid, row / 2);
		failed |= check_i32("scaling coords copy", copy, row % 2);
		failed |= check_i32("scaling coords diploid", diploid_bid, hk_diploid_bid(bid, copy));
		failed |= check_true("scaling coords finite", isfinite(x) && isfinite(y) && isfinite(z));
		++row;
	}
	failed |= check_i32("scaling coords row loop", row, 2 * bmap->n_beads);
	free(text);
	return failed;
}

static int validate_diag_file(const char *path, const struct hk_blind_iter_loop_diag *diag)
{
	char *text = read_file(path);
	char *p;
	char line[2048];
	int n_iter, n_completed, total_chr_flipped, n_bad_iter, n_relax_bad, n_coord_bad;
	double initial_entropy, final_entropy, initial_pU, final_pU;
	double initial_temperature, final_temperature, initial_rho, final_rho;
	double mean_sep, min_sep, max_sep, sum_k;
	int failed = 0;

	failed |= check_true("scaling diag exists", file_exists(path));
	failed |= check_true("scaling diag text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("scaling diag header n_iter", strstr(text, "n_iter") != 0);
	failed |= check_true("scaling diag header final entropy", strstr(text, "final_mean_entropy") != 0);
	failed |= check_true("scaling diag header final pU", strstr(text, "final_mean_pU") != 0);
	failed |= check_true("scaling diag header sum k", strstr(text, "final_sum_wedge_k") != 0);
	failed |= check_no_forbidden_strings("scaling diag", text);
	failed |= check_i32("scaling diag row count", count_lines(text) - 1, 1);
	p = text;
	next_line(&p, line, sizeof(line)); // header
	failed |= check_true("scaling diag data", next_line(&p, line, sizeof(line)));
	failed |= check_i32("scaling diag parsed",
						sscanf(line, "%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t"
							   "%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf",
							   &n_iter, &n_completed, &initial_entropy, &final_entropy,
							   &initial_pU, &final_pU, &initial_temperature, &final_temperature,
							   &initial_rho, &final_rho, &total_chr_flipped, &n_bad_iter,
							   &n_relax_bad, &n_coord_bad, &mean_sep, &min_sep, &max_sep, &sum_k), 18);
	failed |= check_i32("scaling diag n_iter", n_iter, diag->n_iter);
	failed |= check_i32("scaling diag completed", n_completed, diag->n_completed);
	failed |= check_i32("scaling diag bad", n_bad_iter, 0);
	failed |= check_i32("scaling diag relax bad", n_relax_bad, 0);
	failed |= check_i32("scaling diag coord bad", n_coord_bad, 0);
	failed |= check_true("scaling diag entropy finite", isfinite(initial_entropy) && isfinite(final_entropy));
	failed |= check_true("scaling diag pU finite", isfinite(initial_pU) && isfinite(final_pU));
	failed |= check_close("scaling diag temperature start", initial_temperature, HK_BLIND_P9016_SCALING_TEMPERATURE_START);
	failed |= check_close("scaling diag temperature end", final_temperature, HK_BLIND_P9016_SCALING_TEMPERATURE_END);
	failed |= check_close("scaling diag rho start", initial_rho, HK_BLIND_P9016_SCALING_RHO_TRAIN_START);
	failed |= check_close("scaling diag rho end", final_rho, HK_BLIND_P9016_SCALING_RHO_TRAIN_END);
	failed |= check_true("scaling diag sep finite", isfinite(mean_sep) && isfinite(min_sep) && isfinite(max_sep));
	failed |= check_true("scaling diag sum k finite", isfinite(sum_k) && sum_k >= 0.0);
	failed |= check_true("scaling diag flips", total_chr_flipped >= 0);
	free(text);
	return failed;
}

static int run_one_prefix(const struct hk_map *m, int32_t requested_n, const char *base_dir,
						  struct scaling_result *result)
{
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	struct hk_blind_bpair_set *set = 0;
	struct hk_fdg_conf conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	struct hk_blind_single_iter_diag *per_iter = 0;
	fvec3_t *haploid = 0, *diploid = 0, *diploid_before = 0;
	char run_dir[512] = {0};
	char posterior_path[512], raw_path[512], coords_path[512], diag_path[512];
	int32_t n_prefix, n_diploid;
	int32_t i;
	float max_force_l1 = 0.0f;
	int failed = 0;

	memset(result, 0, sizeof(*result));
	result->requested_n = requested_n;
	n_prefix = m->n_pairs < requested_n? m->n_pairs : requested_n;
	failed |= check_true("scaling prefix positive", n_prefix > 0);
	failed |= check_true("scaling prefix bounded", n_prefix <= requested_n);
	if (failed) return 1;

	snprintf(run_dir, sizeof(run_dir), "%s/N_%d", base_dir, requested_n);
	failed |= check_i32("scaling run dir", mkdir(run_dir, 0700), 0);
	if (failed) goto cleanup;
	snprintf(posterior_path, sizeof(posterior_path), "%s/blind_prefix.bpair_posterior.tsv", run_dir);
	snprintf(raw_path, sizeof(raw_path), "%s/blind_prefix.raw_posterior.tsv", run_dir);
	snprintf(coords_path, sizeof(coords_path), "%s/blind_prefix.coords.tsv", run_dir);
	snprintf(diag_path, sizeof(diag_path), "%s/blind_prefix.loop_diag.tsv", run_dir);

	bmap = hk_bmap_gen(m->d, n_prefix, m->pairs, 1000000, 1);
	raw = (struct hk_blind_pair*)calloc(n_prefix, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate scaling bmap/raw data for N=%d\n", requested_n);
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n_prefix; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);

	set = hk_blind_bpair_set_build(bmap, n_prefix, raw);
	failed |= check_binned_set(bmap, set, requested_n);
	if (failed) goto cleanup;

	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	haploid = (fvec3_t*)calloc(bmap->n_beads, sizeof(*haploid));
	diploid = (fvec3_t*)calloc(n_diploid, sizeof(*diploid));
	diploid_before = (fvec3_t*)calloc(n_diploid, sizeof(*diploid_before));
	per_iter = (struct hk_blind_single_iter_diag*)calloc(HK_BLIND_P9016_SCALING_N_ITER, sizeof(*per_iter));
	if (haploid == 0 || diploid == 0 || diploid_before == 0 || per_iter == 0) {
		fprintf(stderr, "failed to allocate scaling data for N=%d\n", requested_n);
		failed = 1;
		goto cleanup;
	}

	set_haploid_scaffold(bmap, haploid);
	failed |= check_i32("scaling init diploid",
						hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
																  diploid, 0.5f, 0.0f, 17), 0);
	failed |= check_coords_finite(diploid, n_diploid);
	if (failed) goto cleanup;
	copy_coords(diploid_before, diploid, n_diploid);

	hk_fdg_conf_init(&conf);
	set_schedule_conf(&schedule_conf);
	failed |= check_i32("scaling scheduled ret",
						hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &conf, diploid, 0,
															 &schedule_conf, per_iter, &loop_diag), 0);
	failed |= check_loop_health(&loop_diag, per_iter);
	for (i = 0; i < HK_BLIND_P9016_SCALING_N_ITER; ++i)
		if (per_iter[i].relax_diag.max_force_l1 > max_force_l1)
			max_force_l1 = per_iter[i].relax_diag.max_force_l1;
	failed |= check_coords_finite(diploid, n_diploid);
	result->moved = coords_any_changed(diploid_before, diploid, n_diploid, 1e-8f);
	if (max_force_l1 > 1e-12f)
		failed |= check_true("scaling coordinate movement", result->moved);
	else
		fprintf(stderr, "NOTE: P9016 scaling force_l1 is zero for N=%d; no movement expected\n", requested_n);
	if (failed) goto cleanup;

	failed |= check_i32("scaling write outputs",
						write_outputs(posterior_path, raw_path, coords_path, diag_path,
									  raw, n_prefix, bmap, set, diploid, &loop_diag), 0);
	failed |= validate_binned_posterior_file(posterior_path, set);
	failed |= validate_raw_posterior_file(raw_path, set);
	failed |= validate_coords_file(coords_path, bmap);
	failed |= validate_diag_file(diag_path, &loop_diag);
	if (failed) goto cleanup;

	result->n_raw = set->n_raw;
	result->n_bpair = set->n_bpairs;
	result->n_beads = bmap->n_beads;
	result->final_mean_entropy = loop_diag.final_mean_entropy;
	result->final_mean_pU = loop_diag.final_mean_pU;
	result->final_sum_wedge_k = loop_diag.final_sum_wedge_k;

	fprintf(stderr,
			"P9016 scaling output: N=%d n_raw=%d n_bpair=%d n_beads=%d "
			"final_entropy=%.8g final_pU=%.8g final_sum_wedge_k=%.8g moved=%d\n",
			requested_n, result->n_raw, result->n_bpair, result->n_beads,
			result->final_mean_entropy, result->final_mean_pU,
			result->final_sum_wedge_k, result->moved);

cleanup:
	remove_run_outputs(run_dir);
	free(haploid);
	free(diploid);
	free(diploid_before);
	free(per_iter);
	hk_blind_bpair_set_destroy(set);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	return failed != 0;
}

int main(void)
{
	struct hk_map *m = 0;
	struct scaling_result results[HK_BLIND_P9016_SCALING_N_PREFIX];
	char base_dir[256] = {0};
	int failed = 0;
	int i;

	if (!file_exists(HK_BLIND_P9016_SCALING_PATH)) {
		fprintf(stderr, "SKIP: %s not found\n", HK_BLIND_P9016_SCALING_PATH);
		return 0;
	}

	hk_verbose = 0;
	m = hk_map_read(HK_BLIND_P9016_SCALING_PATH);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", HK_BLIND_P9016_SCALING_PATH);
		return 1;
	}
	failed |= check_i32("scaling temp dir", make_temp_dir(base_dir, sizeof(base_dir)), 0);
	if (failed) goto cleanup;

	for (i = 0; i < HK_BLIND_P9016_SCALING_N_PREFIX; ++i) {
		failed |= run_one_prefix(m, hk_blind_scaling_prefix_sizes[i], base_dir, &results[i]);
		if (failed) goto cleanup;
		if (i > 0) {
			failed |= check_true("scaling requested increasing",
								 results[i].requested_n > results[i - 1].requested_n);
			failed |= check_true("scaling n_raw increasing", results[i].n_raw > results[i - 1].n_raw);
			failed |= check_true("scaling n_bpair nondecreasing", results[i].n_bpair >= results[i - 1].n_bpair);
		}
	}

cleanup:
	for (i = 0; i < HK_BLIND_P9016_SCALING_N_PREFIX; ++i) {
		char run_dir[512];
		snprintf(run_dir, sizeof(run_dir), "%s/N_%d", base_dir, hk_blind_scaling_prefix_sizes[i]);
		remove_run_outputs(run_dir);
	}
	if (base_dir[0]) rmdir(base_dir);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
