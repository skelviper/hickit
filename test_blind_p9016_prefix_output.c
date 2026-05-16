#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "hickit.h"

#define HK_BLIND_P9016_OUTPUT_PREFIX_N 1000
#define HK_BLIND_P9016_OUTPUT_PATH "../pairs/P9016.pairs.gz"
#define HK_BLIND_P9016_OUTPUT_N_ITER 3
#define HK_BLIND_P9016_OUTPUT_UNIT 1.0f
#define HK_BLIND_P9016_OUTPUT_D_SCALE 1.0f
#define HK_BLIND_P9016_OUTPUT_BASE_K 2.0f
#define HK_BLIND_P9016_OUTPUT_TEMPERATURE_START 2.0f
#define HK_BLIND_P9016_OUTPUT_TEMPERATURE_END 1.0f
#define HK_BLIND_P9016_OUTPUT_RHO_TRAIN_START 1.0f
#define HK_BLIND_P9016_OUTPUT_RHO_TRAIN_END 1.0f
#define HK_BLIND_P9016_OUTPUT_MIN_SEP_UNIT 0.25f
#define HK_BLIND_P9016_OUTPUT_LAMBDA_SEP 0.05f
#define HK_BLIND_P9016_OUTPUT_STEP 0.001f
#define HK_BLIND_P9016_OUTPUT_RELAX_STEPS 5

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
		snprintf(dir, dir_size, "/tmp/hk_blind_p9016_output_%ld_%d", (long)getpid(), i);
		if (mkdir(dir, 0700) == 0)
			return 0;
		if (errno != EEXIST)
			return -1;
	}
	return -1;
}

static void cleanup_outputs(const char *dir, const char *posterior_path, const char *raw_path,
							const char *coords_path, const char *diag_path)
{
	if (posterior_path[0]) remove(posterior_path);
	if (raw_path[0]) remove(raw_path);
	if (coords_path[0]) remove(coords_path);
	if (diag_path[0]) remove(diag_path);
	if (dir[0]) rmdir(dir);
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
			failed |= check_true("output coord finite", isfinite(coords[i][a]));
	return failed;
}

static void set_schedule_conf(struct hk_blind_iter_schedule_conf *conf)
{
	conf->n_iter = HK_BLIND_P9016_OUTPUT_N_ITER;
	conf->base_conf.unit = HK_BLIND_P9016_OUTPUT_UNIT;
	conf->base_conf.d_scale = HK_BLIND_P9016_OUTPUT_D_SCALE;
	conf->base_conf.base_k = HK_BLIND_P9016_OUTPUT_BASE_K;
	conf->base_conf.temperature = HK_BLIND_P9016_OUTPUT_TEMPERATURE_START;
	conf->base_conf.rho_train = HK_BLIND_P9016_OUTPUT_RHO_TRAIN_START;
	conf->base_conf.min_sep_unit = HK_BLIND_P9016_OUTPUT_MIN_SEP_UNIT;
	conf->base_conf.lambda_sep = HK_BLIND_P9016_OUTPUT_LAMBDA_SEP;
	conf->base_conf.relax_step = HK_BLIND_P9016_OUTPUT_STEP;
	conf->base_conf.relax_steps = HK_BLIND_P9016_OUTPUT_RELAX_STEPS;
	conf->base_conf.enable_repulsion = 1;
	conf->base_conf.repulsion_mode = HK_BLIND_REPULSION_CELL;
	conf->base_conf.rho_train_mode = HK_BLIND_RHO_TRAIN_CONSTANT;
	conf->base_conf.d_scale_mode = HK_BLIND_D_SCALE_RAW_COUNT;
	conf->base_conf.d_scale_eps_count = 1e-6f;
	conf->base_conf.rho_train_floor = HK_BLIND_RHO_TRAIN_DEFAULT_FLOOR;
	conf->base_conf.contact_k_multiplier_cis = 1.0f;
	conf->base_conf.contact_k_multiplier_trans = 1.0f;
	conf->temperature_start = HK_BLIND_P9016_OUTPUT_TEMPERATURE_START;
	conf->temperature_end = HK_BLIND_P9016_OUTPUT_TEMPERATURE_END;
	conf->rho_train_start = HK_BLIND_P9016_OUTPUT_RHO_TRAIN_START;
	conf->rho_train_end = HK_BLIND_P9016_OUTPUT_RHO_TRAIN_END;
}

static int check_binned_set(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set, int32_t max_raw)
{
	int failed = 0;
	int32_t i;
	int have_multi_count_bpair = 0;

	failed |= check_true("output set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_true("output raw positive", set->n_raw > 0);
	failed |= check_true("output raw bounded", set->n_raw <= max_raw);
	failed |= check_true("output bpair positive", set->n_bpairs > 0);
	for (i = 0; i < set->n_raw; ++i) {
		failed |= check_true("output raw2binned id lower", set->raw2binned[i].bpair_id >= 0);
		failed |= check_true("output raw2binned id upper", set->raw2binned[i].bpair_id < set->n_bpairs);
		failed |= check_true("output raw2binned swapped", set->raw2binned[i].swapped == 0 || set->raw2binned[i].swapped == 1);
	}
	for (i = 0; i < set->n_bpairs; ++i) {
		const struct hk_blind_bpair *bp = &set->bpairs[i];
		failed |= check_true("output key sorted", bp->key.bid[0] <= bp->key.bid[1]);
		failed |= check_true("output key bid0 valid", bp->key.bid[0] >= 0 && bp->key.bid[0] < bmap->n_beads);
		failed |= check_true("output key bid1 valid", bp->key.bid[1] >= 0 && bp->key.bid[1] < bmap->n_beads);
		failed |= check_true("output bpair n_raw positive", bp->n_raw > 0);
		failed |= check_true("output bpair base_d_scale finite", isfinite(bp->base_d_scale));
		failed |= check_true("output bpair base_d_scale positive", bp->base_d_scale > 0.0f);
		failed |= check_close("output bpair base_d_scale count", bp->base_d_scale, powf((float)bp->n_raw, -1.0f / 3.0f));
		failed |= check_close("output bpair base_k", bp->base_k, 1.0f);
		if (bp->n_raw > 1) {
			have_multi_count_bpair = 1;
			failed |= check_true("output multi-count d_scale below one", bp->base_d_scale < 1.0f);
		}
	}
	failed |= check_true("output has multi-count bpair", have_multi_count_bpair);
	return failed;
}

static int check_param_aware_wedges(const struct hk_bmap *bmap, const struct hk_blind_bpair_set *set)
{
	struct hk_blind_wedge_list list;
	double sum_k = 0.0;
	int failed = 0;
	int32_t i, before_aggregation;

	hk_blind_wedge_list_init(&list);
	failed |= check_i32("output param wedge build", hk_blind_wedge_list_build_from_bpair_set_params(&list, set, 1.0f), 0);
	failed |= check_true("output param wedge input stats", list.n_input_bpair == set->n_bpairs);
	failed |= check_true("output param wedge expanded stats",
						 list.n_expanded_edges == (int64_t)HK_BLIND_N_STATE *
						 (set->n_bpairs - list.n_skipped_same_bin_bpairs));
	failed |= check_true("output param wedge skipped stats", list.n_skipped_self_edges >= 0);
	failed |= check_i32("output param wedge before aggregation", list.n_edges_before_aggregation, list.n_edges);
	for (i = 0; i < list.n_edges; ++i) {
		const struct hk_blind_wedge *edge = &list.edges[i];
		failed |= check_true("output param wedge bid0 valid", edge->bid[0] >= 0 && edge->bid[0] < 2 * bmap->n_beads);
		failed |= check_true("output param wedge bid1 valid", edge->bid[1] >= 0 && edge->bid[1] < 2 * bmap->n_beads);
		failed |= check_true("output param wedge not self", edge->bid[0] < edge->bid[1]);
		failed |= check_true("output param wedge k finite", isfinite(edge->k) && edge->k >= 0.0f);
		failed |= check_true("output param wedge d_scale finite", isfinite(edge->d_scale) && edge->d_scale > 0.0f);
		sum_k += edge->k;
	}
	failed |= check_true("output param wedge sum finite", isfinite(sum_k) && sum_k >= 0.0);
	before_aggregation = list.n_edges;
	failed |= check_i32("output param wedge aggregate", hk_blind_wedge_list_aggregate_exact(&list), 0);
	failed |= check_true("output param wedge aggregate count", list.n_edges <= before_aggregation);
	for (i = 0; i < list.n_edges; ++i) {
		const struct hk_blind_wedge *edge = &list.edges[i];
		failed |= check_true("output param aggregate k finite", isfinite(edge->k) && edge->k >= 0.0f);
		failed |= check_true("output param aggregate d_scale finite", isfinite(edge->d_scale) && edge->d_scale > 0.0f);
	}
	hk_blind_wedge_list_destroy(&list);
	return failed;
}

static int check_loop_health(const struct hk_blind_iter_loop_diag *diag,
							 const struct hk_blind_single_iter_diag *per_iter)
{
	int failed = 0;
	int32_t i;

	failed |= check_i32("output loop n_iter", diag->n_iter, HK_BLIND_P9016_OUTPUT_N_ITER);
	failed |= check_i32("output loop completed", diag->n_completed, HK_BLIND_P9016_OUTPUT_N_ITER);
	failed |= check_i32("output loop bad iter", diag->n_bad_iter, 0);
	failed |= check_i32("output loop relax bad", diag->n_relax_nonfinite_iter, 0);
	failed |= check_i32("output loop coord bad", diag->n_coord_nonfinite, 0);
	failed |= check_close("output initial temperature", diag->initial_temperature, HK_BLIND_P9016_OUTPUT_TEMPERATURE_START);
	failed |= check_close("output final temperature", diag->final_temperature, HK_BLIND_P9016_OUTPUT_TEMPERATURE_END);
	failed |= check_close("output initial rho", diag->initial_rho_train, HK_BLIND_P9016_OUTPUT_RHO_TRAIN_START);
	failed |= check_close("output final rho", diag->final_rho_train, HK_BLIND_P9016_OUTPUT_RHO_TRAIN_END);
	failed |= check_true("output final entropy finite", isfinite(diag->final_mean_entropy));
	failed |= check_true("output final entropy range", diag->final_mean_entropy >= -1e-6f &&
						 diag->final_mean_entropy <= logf((float)HK_BLIND_N_STATE) + 1e-6f);
	failed |= check_true("output final pU finite", isfinite(diag->final_mean_pU));
	failed |= check_true("output final pU range", diag->final_mean_pU >= -1e-6f && diag->final_mean_pU <= 1.0f + 1e-6f);
	failed |= check_true("output final sep finite", isfinite(diag->final_mean_sep));
	failed |= check_true("output final sum k finite", isfinite(diag->final_sum_wedge_k));
	failed |= check_true("output final sum k nonnegative", diag->final_sum_wedge_k >= 0.0);
	for (i = 0; i < HK_BLIND_P9016_OUTPUT_N_ITER; ++i) {
		const struct hk_blind_iter_diag *pre = &per_iter[i].pre_relax_diag;
		failed |= check_i32("output posterior nonfinite", pre->n_posterior_nonfinite, 0);
		failed |= check_i32("output posterior bad sum", pre->n_posterior_bad_sum, 0);
		failed |= check_i32("output posterior out of range", pre->n_posterior_out_of_range, 0);
		failed |= check_i32("output uncertainty nonfinite", pre->n_uncertainty_nonfinite, 0);
		failed |= check_i32("output uncertainty out of range", pre->n_uncertainty_out_of_range, 0);
		failed |= check_i32("output five-state bad sum", pre->n_five_state_bad_sum, 0);
		failed |= check_i32("output wedge nonfinite", pre->n_wedge_nonfinite, 0);
		failed |= check_i32("output wedge bad k", pre->n_wedge_bad_k, 0);
		failed |= check_i32("output wedge bad d_scale", pre->n_wedge_bad_d_scale, 0);
		failed |= check_i32("output sep force nonfinite", pre->sep_force_nonfinite, 0);
		failed |= check_i32("output relax completed", per_iter[i].relax_diag.n_completed, HK_BLIND_P9016_OUTPUT_RELAX_STEPS);
		failed |= check_i32("output relax nonfinite step", per_iter[i].relax_diag.n_nonfinite_step, 0);
		failed |= check_i32("output relax coord bad", per_iter[i].relax_diag.n_coord_nonfinite, 0);
		failed |= check_i32("output backbone bad", per_iter[i].relax_diag.n_backbone_nonfinite_step, 0);
		failed |= check_i32("output repulsion bad", per_iter[i].relax_diag.n_repulsion_nonfinite_step, 0);
		failed |= check_true("output relax total finite", isfinite(per_iter[i].relax_diag.final_total_energy));
		failed |= check_true("output repulsion finite", isfinite(per_iter[i].relax_diag.final_repulsion_energy));
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

static int validate_posterior_file(const char *path, const struct hk_blind_bpair_set *set)
{
	char *text = read_file(path);
	char *p;
	char line[4096];
	int failed = 0;
	int row = 0, checked = 0;

	failed |= check_true("posterior file text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("posterior header p00", strstr(text, "p00") != 0);
	failed |= check_true("posterior header p11", strstr(text, "p11") != 0);
	failed |= check_true("posterior header pU", strstr(text, "pU") != 0);
	failed |= check_true("posterior header base_d_scale", strstr(text, "base_d_scale") != 0);
	failed |= check_true("posterior header base_k", strstr(text, "base_k") != 0);
	failed |= check_true("posterior header entropy", strstr(text, "entropy") != 0);
	failed |= check_true("posterior header margin", strstr(text, "margin") != 0);
	failed |= check_true("posterior header pmax", strstr(text, "pmax") != 0);
	failed |= check_true("posterior header rho", strstr(text, "rho_output") != 0);
	failed |= check_no_forbidden_strings("posterior output", text);
	failed |= check_i32("posterior output rows", count_lines(text) - 1, set->n_bpairs);
	p = text;
	next_line(&p, line, sizeof(line)); // header
	while (next_line(&p, line, sizeof(line))) {
		char chr1[128], chr2[128];
		int st1, en1, st2, en2, bid1, bid2, n_raw;
		double base_d_scale, base_k;
		double p00, p01, p10, p11, pU, psame, pcross, entropy, margin, pmax, rho;
		if (checked < 8) {
			failed |= check_i32("posterior output parsed",
								sscanf(line, "%127s\t%d\t%d\t%127s\t%d\t%d\t%d\t%d\t%d\t"
									   "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
									   chr1, &st1, &en1, chr2, &st2, &en2, &bid1, &bid2, &n_raw,
									   &base_d_scale, &base_k, &p00, &p01, &p10, &p11, &pU, &psame, &pcross,
									   &entropy, &margin, &pmax, &rho), 22);
			failed |= check_true("posterior output chr nonempty", chr1[0] != 0 && chr2[0] != 0);
			failed |= check_true("posterior output intervals", st1 <= en1 && st2 <= en2);
			failed |= check_true("posterior output bids valid", bid1 >= 0 && bid2 >= bid1);
			failed |= check_true("posterior output n_raw positive", n_raw > 0);
			failed |= check_close("posterior output base_d_scale", base_d_scale, pow((double)n_raw, -1.0 / 3.0));
			failed |= check_close("posterior output base_k", base_k, 1.0);
			failed |= check_true("posterior output finite", isfinite(p00) && isfinite(p01) &&
								 isfinite(p10) && isfinite(p11) && isfinite(pU) &&
								 isfinite(entropy) && isfinite(margin) && isfinite(pmax) &&
								 isfinite(rho));
			failed |= check_close("posterior output p4 sum", p00 + p01 + p10 + p11, 1.0);
			failed |= check_close("posterior output five-state sum", rho * (p00 + p01 + p10 + p11) + pU, 1.0);
			failed |= check_close("posterior output psame", psame, p00 + p11);
			failed |= check_close("posterior output pcross", pcross, p01 + p10);
			++checked;
		}
		++row;
	}
	failed |= check_i32("posterior output row loop", row, set->n_bpairs);
	free(text);
	return failed;
}

static int validate_raw_posterior_file(const char *path, const struct hk_blind_bpair_set *set)
{
	char *text = read_file(path);
	char *p;
	char line[4096];
	int failed = 0;
	int row = 0, checked = 0;

	failed |= check_true("raw posterior file text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("raw posterior header raw_id", strstr(text, "raw_id") != 0);
	failed |= check_true("raw posterior header chr1", strstr(text, "chr1") != 0);
	failed |= check_true("raw posterior header pos1", strstr(text, "pos1") != 0);
	failed |= check_true("raw posterior header swapped", strstr(text, "swapped") != 0);
	failed |= check_true("raw posterior header p00", strstr(text, "p00") != 0);
	failed |= check_true("raw posterior header p11", strstr(text, "p11") != 0);
	failed |= check_true("raw posterior header pU", strstr(text, "pU") != 0);
	failed |= check_true("raw posterior header entropy", strstr(text, "entropy") != 0);
	failed |= check_true("raw posterior header margin", strstr(text, "margin") != 0);
	failed |= check_true("raw posterior header pmax", strstr(text, "pmax") != 0);
	failed |= check_true("raw posterior header rho", strstr(text, "rho_output") != 0);
	failed |= check_no_forbidden_strings("raw posterior output", text);
	failed |= check_i32("raw posterior output rows", count_lines(text) - 1, set->n_raw);
	p = text;
	next_line(&p, line, sizeof(line)); // header
	while (next_line(&p, line, sizeof(line))) {
		char chr1[128], chr2[128];
		int raw_id, pos1, pos2, bid1_raw, bid2_raw, bpair_id, bid1_can, bid2_can, swapped;
		double p00, p01, p10, p11, pU, psame, pcross, entropy, margin, pmax, rho;
		if (checked < 12) {
			failed |= check_i32("raw posterior parsed",
								sscanf(line, "%d\t%127s\t%d\t%127s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t"
									   "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
									   &raw_id, chr1, &pos1, chr2, &pos2, &bid1_raw, &bid2_raw,
									   &bpair_id, &bid1_can, &bid2_can, &swapped,
									   &p00, &p01, &p10, &p11, &pU, &psame, &pcross,
									   &entropy, &margin, &pmax, &rho), 22);
			failed |= check_i32("raw posterior raw id", raw_id, row);
			failed |= check_true("raw posterior chr nonempty", chr1[0] != 0 && chr2[0] != 0);
			failed |= check_true("raw posterior positions nonnegative", pos1 >= 0 && pos2 >= 0);
			failed |= check_true("raw posterior raw bids valid", bid1_raw >= 0 && bid2_raw >= 0);
			failed |= check_true("raw posterior bpair id valid", bpair_id >= 0 && bpair_id < set->n_bpairs);
			failed |= check_true("raw posterior canonical sorted", bid1_can <= bid2_can);
			failed |= check_true("raw posterior swapped boolean", swapped == 0 || swapped == 1);
			failed |= check_true("raw posterior finite", isfinite(p00) && isfinite(p01) &&
								 isfinite(p10) && isfinite(p11) && isfinite(pU) &&
								 isfinite(entropy) && isfinite(margin) && isfinite(pmax) &&
								 isfinite(rho));
			failed |= check_close("raw posterior p4 sum", p00 + p01 + p10 + p11, 1.0);
			failed |= check_close("raw posterior five-state sum", rho * (p00 + p01 + p10 + p11) + pU, 1.0);
			failed |= check_close("raw posterior psame", psame, p00 + p11);
			failed |= check_close("raw posterior pcross", pcross, p01 + p10);
			++checked;
		}
		++row;
	}
	failed |= check_i32("raw posterior row loop", row, set->n_raw);
	free(text);
	return failed;
}

static int validate_coords_file(const char *path, const struct hk_bmap *bmap)
{
	char *text = read_file(path);
	char *p;
	char line[2048];
	int failed = 0;
	int row = 0, checked = 0;

	failed |= check_true("coords file text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("coords header diploid", strstr(text, "diploid_bid") != 0);
	failed |= check_no_forbidden_strings("coords output", text);
	failed |= check_i32("coords output rows", count_lines(text) - 1, 2 * bmap->n_beads);
	p = text;
	next_line(&p, line, sizeof(line)); // header
	while (next_line(&p, line, sizeof(line))) {
		char chr[128];
		int st, en, bid, copy, diploid_bid;
		double x, y, z;
		if (checked < 12) {
			failed |= check_i32("coords output parsed",
								sscanf(line, "%127s\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf",
									   chr, &st, &en, &bid, &copy, &diploid_bid, &x, &y, &z), 9);
			failed |= check_true("coords output chr nonempty", chr[0] != 0);
			failed |= check_true("coords output interval", st <= en);
			failed |= check_i32("coords output bid order", bid, row / 2);
			failed |= check_i32("coords output copy order", copy, row % 2);
			failed |= check_i32("coords output diploid", diploid_bid, hk_diploid_bid(bid, copy));
			failed |= check_true("coords output finite", isfinite(x) && isfinite(y) && isfinite(z));
			++checked;
		}
		++row;
	}
	failed |= check_i32("coords output row loop", row, 2 * bmap->n_beads);
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

	failed |= check_true("diag file text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("diag header n_iter", strstr(text, "n_iter") != 0);
	failed |= check_true("diag header n_completed", strstr(text, "n_completed") != 0);
	failed |= check_true("diag header final entropy", strstr(text, "final_mean_entropy") != 0);
	failed |= check_true("diag header final pU", strstr(text, "final_mean_pU") != 0);
	failed |= check_true("diag header final sum k", strstr(text, "final_sum_wedge_k") != 0);
	failed |= check_no_forbidden_strings("diag output", text);
	failed |= check_i32("diag output rows", count_lines(text) - 1, 1);
	p = text;
	next_line(&p, line, sizeof(line)); // header
	failed |= check_true("diag data line", next_line(&p, line, sizeof(line)));
	failed |= check_i32("diag output parsed",
						sscanf(line, "%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t"
							   "%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf",
							   &n_iter, &n_completed, &initial_entropy, &final_entropy,
							   &initial_pU, &final_pU, &initial_temperature, &final_temperature,
							   &initial_rho, &final_rho, &total_chr_flipped, &n_bad_iter,
							   &n_relax_bad, &n_coord_bad, &mean_sep, &min_sep, &max_sep, &sum_k), 18);
	failed |= check_i32("diag output n_iter", n_iter, diag->n_iter);
	failed |= check_i32("diag output completed", n_completed, diag->n_completed);
	failed |= check_i32("diag output bad iter", n_bad_iter, 0);
	failed |= check_i32("diag output relax bad", n_relax_bad, 0);
	failed |= check_i32("diag output coord bad", n_coord_bad, 0);
	failed |= check_true("diag output initial entropy finite", isfinite(initial_entropy));
	failed |= check_close("diag output final entropy", final_entropy, diag->final_mean_entropy);
	failed |= check_true("diag output initial pU finite", isfinite(initial_pU));
	failed |= check_close("diag output final pU", final_pU, diag->final_mean_pU);
	failed |= check_close("diag output temperature start", initial_temperature, HK_BLIND_P9016_OUTPUT_TEMPERATURE_START);
	failed |= check_close("diag output temperature end", final_temperature, HK_BLIND_P9016_OUTPUT_TEMPERATURE_END);
	failed |= check_close("diag output rho start", initial_rho, HK_BLIND_P9016_OUTPUT_RHO_TRAIN_START);
	failed |= check_close("diag output rho end", final_rho, HK_BLIND_P9016_OUTPUT_RHO_TRAIN_END);
	failed |= check_true("diag output mean sep finite", isfinite(mean_sep) && isfinite(min_sep) && isfinite(max_sep));
	failed |= check_true("diag output sum k finite", isfinite(sum_k) && sum_k >= 0.0);
	failed |= check_true("diag output flips nonnegative", total_chr_flipped >= 0);
	free(text);
	return failed;
}

int main(void)
{
	struct hk_map *m = 0;
	struct hk_bmap *bmap = 0;
	struct hk_blind_pair *raw = 0;
	struct hk_blind_bpair_set *set = 0;
	struct hk_fdg_conf conf;
	struct hk_blind_iter_schedule_conf schedule_conf;
	struct hk_blind_iter_loop_diag loop_diag;
	struct hk_blind_single_iter_diag *per_iter = 0;
	fvec3_t *haploid = 0, *diploid = 0, *diploid_before = 0;
	char temp_dir[256] = {0};
	char posterior_path[512] = {0};
	char raw_path[512] = {0};
	char coords_path[512] = {0};
	char diag_path[512] = {0};
	int32_t n_prefix, n_diploid;
	int32_t i;
	float max_force_l1 = 0.0f;
	int moved = 0;
	int failed = 0;

	if (!file_exists(HK_BLIND_P9016_OUTPUT_PATH)) {
		fprintf(stderr, "SKIP: %s not found\n", HK_BLIND_P9016_OUTPUT_PATH);
		return 0;
	}

	hk_verbose = 0;
	m = hk_map_read(HK_BLIND_P9016_OUTPUT_PATH);
	if (m == 0) {
		fprintf(stderr, "failed to read %s\n", HK_BLIND_P9016_OUTPUT_PATH);
		return 1;
	}
	n_prefix = m->n_pairs < HK_BLIND_P9016_OUTPUT_PREFIX_N? m->n_pairs : HK_BLIND_P9016_OUTPUT_PREFIX_N;
	if (n_prefix <= 0) {
		fprintf(stderr, "no pairs in %s\n", HK_BLIND_P9016_OUTPUT_PATH);
		failed = 1;
		goto cleanup;
	}

	bmap = hk_bmap_gen(m->d, n_prefix, m->pairs, 1000000, 1);
	raw = (struct hk_blind_pair*)calloc(n_prefix, sizeof(*raw));
	if (bmap == 0 || raw == 0) {
		fprintf(stderr, "failed to allocate P9016 output bmap/raw data\n");
		failed = 1;
		goto cleanup;
	}
	for (i = 0; i < n_prefix; ++i)
		hk_blind_pair_from_pair(&raw[i], &m->pairs[i]);

	set = hk_blind_bpair_set_build(bmap, n_prefix, raw);
	failed |= check_binned_set(bmap, set, HK_BLIND_P9016_OUTPUT_PREFIX_N);
	failed |= check_param_aware_wedges(bmap, set);
	if (failed) goto cleanup;

	n_diploid = bmap->n_beads * HK_DIPLOID_N_COPY;
	haploid = (fvec3_t*)calloc(bmap->n_beads, sizeof(*haploid));
	diploid = (fvec3_t*)calloc(n_diploid, sizeof(*diploid));
	diploid_before = (fvec3_t*)calloc(n_diploid, sizeof(*diploid_before));
	per_iter = (struct hk_blind_single_iter_diag*)calloc(HK_BLIND_P9016_OUTPUT_N_ITER, sizeof(*per_iter));
	if (haploid == 0 || diploid == 0 || diploid_before == 0 || per_iter == 0) {
		fprintf(stderr, "failed to allocate P9016 output data\n");
		failed = 1;
		goto cleanup;
	}
	set_haploid_scaffold(bmap, haploid);
	failed |= check_i32("output init diploid",
						hk_blind_init_diploid_coords_from_haploid(bmap, haploid, bmap->n_beads,
																  diploid, 0.5f, 0.0f, 17), 0);
	failed |= check_coords_finite(diploid, n_diploid);
	if (failed) goto cleanup;
	copy_coords(diploid_before, diploid, n_diploid);

	hk_fdg_conf_init(&conf);
	set_schedule_conf(&schedule_conf);
	failed |= check_i32("output scheduled ret",
						hk_blind_run_iter_loop_scheduled_cpu(bmap, set, &conf, diploid, 0,
															 &schedule_conf, per_iter, &loop_diag), 0);
	failed |= check_loop_health(&loop_diag, per_iter);
	for (i = 0; i < HK_BLIND_P9016_OUTPUT_N_ITER; ++i)
		if (per_iter[i].relax_diag.max_force_l1 > max_force_l1)
			max_force_l1 = per_iter[i].relax_diag.max_force_l1;
	failed |= check_coords_finite(diploid, n_diploid);
	moved = coords_any_changed(diploid_before, diploid, n_diploid, 1e-8f);
	if (max_force_l1 > 1e-12f)
		failed |= check_true("output coordinate movement", moved);
	else
		fprintf(stderr, "NOTE: P9016 output smoke force_l1 is zero; no movement expected\n");
	if (failed) goto cleanup;

	failed |= check_i32("output temp dir", make_temp_dir(temp_dir, sizeof(temp_dir)), 0);
	if (failed) goto cleanup;
	snprintf(posterior_path, sizeof(posterior_path), "%s/blind_prefix.bpair_posterior.tsv", temp_dir);
	snprintf(raw_path, sizeof(raw_path), "%s/blind_prefix.raw_posterior.tsv", temp_dir);
	snprintf(coords_path, sizeof(coords_path), "%s/blind_prefix.coords.tsv", temp_dir);
	snprintf(diag_path, sizeof(diag_path), "%s/blind_prefix.loop_diag.tsv", temp_dir);
	failed |= check_i32("output write files",
						write_outputs(posterior_path, raw_path, coords_path, diag_path,
									  raw, n_prefix, bmap, set, diploid, &loop_diag), 0);
	failed |= validate_posterior_file(posterior_path, set);
	failed |= validate_raw_posterior_file(raw_path, set);
	failed |= validate_coords_file(coords_path, bmap);
	failed |= validate_diag_file(diag_path, &loop_diag);

	if (!failed) {
		fprintf(stderr,
				"P9016 output smoke: n_raw=%d n_bpair=%d n_iter=%d output_dir=%s "
				"final_entropy=%.8g final_pU=%.8g final_sum_wedge_k=%.8g moved=%d\n",
				set->n_raw, set->n_bpairs, loop_diag.n_completed, temp_dir,
				loop_diag.final_mean_entropy, loop_diag.final_mean_pU,
				loop_diag.final_sum_wedge_k, moved);
	}

cleanup:
	cleanup_outputs(temp_dir, posterior_path, raw_path, coords_path, diag_path);
	free(haploid);
	free(diploid);
	free(diploid_before);
	free(per_iter);
	hk_blind_bpair_set_destroy(set);
	free(raw);
	if (bmap) hk_bmap_destroy(bmap);
	if (m) hk_map_destroy(m);
	return failed != 0;
}
