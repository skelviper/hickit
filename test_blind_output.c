#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "hickit.h"

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

static char *read_stream(FILE *fp)
{
	long len;
	char *buf;
	size_t n_read;

	fflush(fp);
	if (fseek(fp, 0, SEEK_END) != 0) return 0;
	len = ftell(fp);
	if (len < 0) return 0;
	if (fseek(fp, 0, SEEK_SET) != 0) return 0;
	buf = (char*)malloc((size_t)len + 1);
	if (buf == 0) return 0;
	n_read = fread(buf, 1, (size_t)len, fp);
	buf[n_read] = 0;
	return buf;
}

static char *read_gz_path(const char *path)
{
	gzFile fp;
	char tmp[4096];
	char *buf = 0, *next;
	size_t cap = 0, len = 0;
	int n;

	fp = gzopen(path, "rb");
	if (fp == 0)
		return 0;
	for (;;) {
		n = gzread(fp, tmp, sizeof(tmp));
		if (n < 0) {
			free(buf);
			gzclose(fp);
			return 0;
		}
		if (n == 0)
			break;
		if (len + (size_t)n + 1 > cap) {
			size_t new_cap = cap? cap * 2 : 8192;
			while (new_cap < len + (size_t)n + 1)
				new_cap *= 2;
			next = (char*)realloc(buf, new_cap);
			if (next == 0) {
				free(buf);
				gzclose(fp);
				return 0;
			}
			buf = next;
			cap = new_cap;
		}
		memcpy(buf + len, tmp, (size_t)n);
		len += (size_t)n;
	}
	if (gzclose(fp) != Z_OK) {
		free(buf);
		return 0;
	}
	if (buf == 0) {
		buf = (char*)malloc(1);
		if (buf == 0)
			return 0;
	}
	buf[len] = 0;
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
	char buf[96];
	const char *forbidden[] = { "phase0", "phase1", "truth", "oracle" };
	int i;
	for (i = 0; i < 4; ++i) {
		snprintf(buf, sizeof(buf), "%s no %s", label, forbidden[i]);
		failed |= check_true(buf, strstr(s, forbidden[i]) == 0);
	}
	return failed;
}

static int tsv_col_index(const char *header, const char *name)
{
	char *copy, *tok;
	int idx = 0;
	copy = strdup(header);
	if (copy == 0)
		return -1;
	for (tok = strtok(copy, "\t\n"); tok; tok = strtok(0, "\t\n"), ++idx) {
		if (strcmp(tok, name) == 0) {
			free(copy);
			return idx;
		}
	}
	free(copy);
	return -1;
}

static int tsv_field_at(const char *line, int target_col, char *out, size_t out_size)
{
	const char *p, *start;
	int col = 0;
	size_t len;
	assert(out);
	assert(out_size > 0);
	out[0] = 0;
	if (line == 0 || target_col < 0)
		return -1;
	p = line;
	while (col < target_col) {
		p = strchr(p, '\t');
		if (p == 0)
			return -1;
		++p;
		++col;
	}
	start = p;
	while (*p && *p != '\t' && *p != '\n')
		++p;
	len = (size_t)(p - start);
	if (len >= out_size)
		len = out_size - 1;
	memcpy(out, start, len);
	out[len] = 0;
	return 0;
}

static int tsv_get_ll(const char *header, const char *line, const char *name, long long *out)
{
	char buf[128];
	int col = tsv_col_index(header, name);
	if (col < 0 || tsv_field_at(line, col, buf, sizeof(buf)) != 0)
		return -1;
	*out = strtoll(buf, 0, 10);
	return 0;
}

static int tsv_get_double(const char *header, const char *line, const char *name, double *out)
{
	char buf[128];
	int col = tsv_col_index(header, name);
	if (col < 0 || tsv_field_at(line, col, buf, sizeof(buf)) != 0)
		return -1;
	*out = strtod(buf, 0);
	return 0;
}

static int tsv_get_int(const char *header, const char *line, const char *name, int *out)
{
	long long v;
	if (tsv_get_ll(header, line, name, &v) != 0)
		return -1;
	*out = (int)v;
	return 0;
}

static void set_bmap(struct hk_bmap *bmap, struct hk_sdict *dict, char **names, int32_t *len,
					 struct hk_bead *beads)
{
	memset(dict, 0, sizeof(*dict));
	dict->n = 2;
	dict->m = 2;
	dict->name = names;
	dict->len = len;
	memset(bmap, 0, sizeof(*bmap));
	bmap->d = dict;
	bmap->n_beads = 3;
	bmap->unit = 1.0f;
	bmap->beads = beads;
	beads[0].chr = 0; beads[0].st = 0;       beads[0].en = 1000000;
	beads[1].chr = 0; beads[1].st = 1000000; beads[1].en = 2000000;
	beads[2].chr = 1; beads[2].st = 0;       beads[2].en = 1000000;
}

static void set_bpair(struct hk_blind_bpair *bp, int32_t bid0, int32_t bid1, int32_t n_raw,
					  float p00, float p01, float p10, float p11,
					  float pU, float entropy, float margin, float pmax, float rho_output)
{
	bp->key.bid[0] = bid0;
	bp->key.bid[1] = bid1;
	bp->n_raw = n_raw;
	bp->n_locked_raw = 0;
	bp->base_d_scale = powf((float)n_raw, -1.0f / 3.0f);
	bp->base_k = 1.0f;
	bp->contact_class = bid0 == 0 && bid1 == 2? HK_BLIND_CONTACT_TRANS : HK_BLIND_CONTACT_CIS;
	bp->lock_mode = HK_BLIND_LOCK_NONE;
	bp->locked_state = -1;
	bp->lock_conflict = 0;
	bp->p4[HK_BLIND_STATE_00] = p00;
	bp->p4[HK_BLIND_STATE_01] = p01;
	bp->p4[HK_BLIND_STATE_10] = p10;
	bp->p4[HK_BLIND_STATE_11] = p11;
	bp->pU = pU;
	bp->entropy = entropy;
	bp->margin = margin;
	bp->pmax = pmax;
	bp->rho_output = rho_output;
}

static void set_bpair_set(struct hk_blind_bpair_set *set, struct hk_blind_bpair *bpairs,
						  struct hk_blind_raw2binned *raw2binned)
{
	set_bpair(&bpairs[0], 0, 2, 3, 0.70f, 0.10f, 0.15f, 0.05f, 0.20f, 0.85f, 0.55f, 0.70f, 0.80f);
	set_bpair(&bpairs[1], 1, 1, 1, 0.25f, 0.25f, 0.25f, 0.25f, 1.00f, logf(4.0f), 0.0f, 0.25f, 0.0f);
	hk_blind_raw2binned_set(&raw2binned[0], 0, 0);
	hk_blind_raw2binned_set(&raw2binned[1], 0, 1);
	hk_blind_raw2binned_set(&raw2binned[2], 0, 0);
	hk_blind_raw2binned_set(&raw2binned[3], 1, 0);
	set->bpairs = bpairs;
	set->n_bpairs = 2;
	set->raw2binned = raw2binned;
	set->n_raw = 4;
}

static void set_coords(fvec3_t coords[6])
{
	int32_t i;
	for (i = 0; i < 6; ++i) {
		coords[i][0] = 0.5f * (float)i;
		coords[i][1] = 1.0f + 0.25f * (float)i;
		coords[i][2] = -0.125f * (float)i;
	}
}

static struct hk_blind_pair make_raw_pair(int32_t bid0, int32_t bid1)
{
	struct hk_blind_pair p;
	p.chr[0] = 0;
	p.chr[1] = 0;
	p.pos[0] = bid0 * 1000000;
	p.pos[1] = bid1 * 1000000;
	p.strand[0] = 1;
	p.strand[1] = -1;
	return p;
}

static void set_raw_writer_bmap(struct hk_bmap *bmap, struct hk_sdict *dict, char **names,
								int32_t *len, struct hk_bead *beads, uint64_t *offcnt)
{
	int32_t i;

	memset(dict, 0, sizeof(*dict));
	dict->n = 1;
	dict->m = 1;
	dict->name = names;
	dict->len = len;
	memset(bmap, 0, sizeof(*bmap));
	bmap->d = dict;
	bmap->n_beads = 6;
	bmap->unit = 1.0f;
	bmap->beads = beads;
	bmap->offcnt = offcnt;
	for (i = 0; i < bmap->n_beads; ++i) {
		beads[i].chr = 0;
		beads[i].st = i * 1000000;
		beads[i].en = (i + 1) * 1000000;
	}
	offcnt[0] = (uint64_t)0 << 32 | 6u;
	offcnt[1] = (uint64_t)6 << 32 | 0u;
}

static int test_posterior_writer(void)
{
	struct hk_sdict dict;
	char *names[2] = { "chrA", "chrB" };
	int32_t len[2] = { 2000000, 1000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[3];
	struct hk_blind_bpair bpairs[2];
	struct hk_blind_raw2binned raw2binned[4];
	struct hk_blind_bpair_set set;
	FILE *fp;
	char *text = 0, *line;
	char chr1[64], chr2[64];
	int st1, en1, st2, en2, bid1, bid2, n_raw;
	double base_d_scale, base_k;
	double p00, p01, p10, p11, pU, psame, pcross, entropy, margin, pmax, rho;
	int failed = 0;

	set_bmap(&bmap, &dict, names, len, beads);
	set_bpair_set(&set, bpairs, raw2binned);
	fp = tmpfile();
	failed |= check_true("posterior tmpfile", fp != 0);
	if (fp == 0) return 1;
	failed |= check_i32("posterior writer ret", hk_blind_write_bpair_posterior_tsv(fp, &bmap, &set), 0);
	text = read_stream(fp);
	fclose(fp);
	failed |= check_true("posterior text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("posterior header p00", strstr(text, "p00") != 0);
	failed |= check_true("posterior header pU", strstr(text, "pU") != 0);
	failed |= check_true("posterior header base_d_scale", strstr(text, "base_d_scale") != 0);
	failed |= check_true("posterior header base_k", strstr(text, "base_k") != 0);
	failed |= check_true("posterior header rho", strstr(text, "rho_output") != 0);
	failed |= check_true("posterior header contact class", strstr(text, "contact_class") != 0);
	failed |= check_true("posterior trans contact class", strstr(text, "\ttrans\t") != 0);
	failed |= check_i32("posterior row count", count_lines(text) - 1, set.n_bpairs);
	failed |= check_no_forbidden_strings("posterior", text);
	line = strchr(text, '\n');
	failed |= check_true("posterior data line", line != 0);
	if (line) {
		++line;
		failed |= check_i32("posterior parsed fields",
							sscanf(line, "%63s\t%d\t%d\t%63s\t%d\t%d\t%d\t%d\t%d\t"
								   "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
								   chr1, &st1, &en1, chr2, &st2, &en2, &bid1, &bid2, &n_raw,
								   &base_d_scale, &base_k, &p00, &p01, &p10, &p11, &pU, &psame, &pcross,
								   &entropy, &margin, &pmax, &rho), 22);
		failed |= check_true("posterior chr1", strcmp(chr1, "chrA") == 0);
		failed |= check_true("posterior chr2", strcmp(chr2, "chrB") == 0);
		failed |= check_i32("posterior start1", st1, 0);
		failed |= check_i32("posterior end1", en1, 1000000);
		failed |= check_i32("posterior start2", st2, 0);
		failed |= check_i32("posterior end2", en2, 1000000);
		failed |= check_i32("posterior bid1", bid1, 0);
		failed |= check_i32("posterior bid2", bid2, 2);
		failed |= check_i32("posterior n_raw", n_raw, 3);
		failed |= check_close("posterior base_d_scale", base_d_scale, pow(3.0, -1.0 / 3.0));
		failed |= check_close("posterior base_k", base_k, 1.0);
		failed |= check_close("posterior p4 sum", p00 + p01 + p10 + p11, 1.0);
		failed |= check_close("posterior five-state sum", rho * (p00 + p01 + p10 + p11) + pU, 1.0);
		failed |= check_close("posterior psame", psame, p00 + p11);
		failed |= check_close("posterior pcross", pcross, p01 + p10);
		failed |= check_true("posterior entropy finite", isfinite(entropy));
		failed |= check_true("posterior margin finite", isfinite(margin));
		failed |= check_true("posterior pmax finite", isfinite(pmax));
	}
	free(text);
	return failed;
}

static int test_posterior_writer_contact_params(void)
{
	struct hk_sdict dict;
	char *names[2] = { "chrA", "chrB" };
	int32_t len[2] = { 2000000, 1000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[3];
	struct hk_blind_bpair bpairs[3];
	struct hk_blind_bpair_set set;
	FILE *fp;
	char *text = 0, *p, *line;
	double expected[3] = { 1.0, 0.5, 1.0 / 3.0 };
	int failed = 0;
	int i;

	set_bmap(&bmap, &dict, names, len, beads);
	set_bpair(&bpairs[0], 0, 1, 1, 0.25f, 0.25f, 0.25f, 0.25f, 1.0f, logf(4.0f), 0.0f, 0.25f, 0.0f);
	set_bpair(&bpairs[1], 0, 2, 8, 0.25f, 0.25f, 0.25f, 0.25f, 1.0f, logf(4.0f), 0.0f, 0.25f, 0.0f);
	set_bpair(&bpairs[2], 1, 2, 27, 0.25f, 0.25f, 0.25f, 0.25f, 1.0f, logf(4.0f), 0.0f, 0.25f, 0.0f);
	set.bpairs = bpairs;
	set.n_bpairs = 3;
	set.raw2binned = 0;
	set.n_raw = 36;
	fp = tmpfile();
	failed |= check_true("posterior params tmpfile", fp != 0);
	if (fp == 0) return 1;
	failed |= check_i32("posterior params writer ret", hk_blind_write_bpair_posterior_tsv(fp, &bmap, &set), 0);
	text = read_stream(fp);
	fclose(fp);
	failed |= check_true("posterior params text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("posterior params header base_d_scale", strstr(text, "base_d_scale") != 0);
	p = text;
	line = strchr(p, '\n');
	failed |= check_true("posterior params header line", line != 0);
	if (line) p = line + 1;
	for (i = 0; i < 3 && p && *p; ++i) {
		char chr1[64], chr2[64];
		int st1, en1, st2, en2, bid1, bid2, n_raw;
		double base_d_scale, base_k, p00, p01, p10, p11, pU, psame, pcross, entropy, margin, pmax, rho;
		failed |= check_i32("posterior params parsed",
							sscanf(p, "%63s\t%d\t%d\t%63s\t%d\t%d\t%d\t%d\t%d\t"
								   "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
								   chr1, &st1, &en1, chr2, &st2, &en2, &bid1, &bid2, &n_raw,
								   &base_d_scale, &base_k, &p00, &p01, &p10, &p11, &pU, &psame, &pcross,
								   &entropy, &margin, &pmax, &rho), 22);
		failed |= check_close("posterior params base_d_scale", base_d_scale, expected[i]);
		failed |= check_close("posterior params base_k", base_k, 1.0);
		p = strchr(p, '\n');
		if (p) ++p;
	}
	free(text);
	return failed;
}

static int test_raw_contact_writer(void)
{
	struct hk_sdict dict;
	char *names[1] = { "chrA" };
	int32_t len[1] = { 6000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[6];
	uint64_t offcnt[2];
	struct hk_blind_pair raw[2];
	struct hk_blind_bpair_set *set;
	FILE *fp;
	char *text = 0, *line0, *line1;
	char chr1[64], chr2[64];
	int raw_id, pos1, pos2, bid1_raw, bid2_raw, bpair_id, bid1_can, bid2_can, swapped;
	double p00, p01, p10, p11, pU, psame, pcross, entropy, margin, pmax, rho;
	int failed = 0;

	set_raw_writer_bmap(&bmap, &dict, names, len, beads, offcnt);
	raw[0] = make_raw_pair(2, 5);
	raw[1] = make_raw_pair(5, 2);
	set = hk_blind_bpair_set_build(&bmap, 2, raw);
	failed |= check_true("raw set exists", set != 0);
	if (set == 0) return 1;
	failed |= check_i32("raw set n_bpairs", set->n_bpairs, 1);
	failed |= check_i32("raw set key0", set->bpairs[0].key.bid[0], 2);
	failed |= check_i32("raw set key1", set->bpairs[0].key.bid[1], 5);
	set_bpair(&set->bpairs[0], 2, 5, 2, 0.10f, 0.20f, 0.30f, 0.40f,
			  0.25f, 1.20f, 0.10f, 0.40f, 0.75f);

	fp = tmpfile();
	failed |= check_true("raw writer tmpfile", fp != 0);
	if (fp == 0) {
		hk_blind_bpair_set_destroy(set);
		return 1;
	}
	failed |= check_i32("raw writer ret", hk_blind_write_raw_contact_posterior_tsv(fp, raw, 2, &bmap, set), 0);
	text = read_stream(fp);
	fclose(fp);
	failed |= check_true("raw writer text", text != 0);
	if (text == 0) {
		hk_blind_bpair_set_destroy(set);
		return 1;
	}
	failed |= check_true("raw header raw_id", strstr(text, "raw_id") != 0);
	failed |= check_true("raw header chr1", strstr(text, "chr1") != 0);
	failed |= check_true("raw header pos1", strstr(text, "pos1") != 0);
	failed |= check_true("raw header swapped", strstr(text, "swapped") != 0);
	failed |= check_true("raw header p00", strstr(text, "p00") != 0);
	failed |= check_true("raw header pU", strstr(text, "pU") != 0);
	failed |= check_true("raw header entropy", strstr(text, "entropy") != 0);
	failed |= check_true("raw header margin", strstr(text, "margin") != 0);
	failed |= check_true("raw header pmax", strstr(text, "pmax") != 0);
	failed |= check_true("raw header rho", strstr(text, "rho_output") != 0);
	failed |= check_true("raw header contact class", strstr(text, "contact_class") != 0);
	failed |= check_i32("raw writer row count", count_lines(text) - 1, 2);
	failed |= check_no_forbidden_strings("raw writer", text);

	line0 = strchr(text, '\n');
	failed |= check_true("raw first line", line0 != 0);
	if (line0) {
		++line0;
		failed |= check_i32("raw first parsed",
							sscanf(line0, "%d\t%63s\t%d\t%63s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t"
								   "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
								   &raw_id, chr1, &pos1, chr2, &pos2, &bid1_raw, &bid2_raw,
								   &bpair_id, &bid1_can, &bid2_can, &swapped,
								   &p00, &p01, &p10, &p11, &pU, &psame, &pcross,
								   &entropy, &margin, &pmax, &rho), 22);
		failed |= check_i32("raw first id", raw_id, 0);
		failed |= check_true("raw first chr1", strcmp(chr1, "chrA") == 0);
		failed |= check_true("raw first chr2", strcmp(chr2, "chrA") == 0);
		failed |= check_i32("raw first pos1", pos1, 2000000);
		failed |= check_i32("raw first pos2", pos2, 5000000);
		failed |= check_i32("raw first bid1 raw", bid1_raw, 2);
		failed |= check_i32("raw first bid2 raw", bid2_raw, 5);
		failed |= check_i32("raw first bpair id", bpair_id, 0);
		failed |= check_i32("raw first canonical0", bid1_can, 2);
		failed |= check_i32("raw first canonical1", bid2_can, 5);
		failed |= check_i32("raw first swapped", swapped, 0);
		failed |= check_close("raw first p00", p00, 0.10);
		failed |= check_close("raw first p01", p01, 0.20);
		failed |= check_close("raw first p10", p10, 0.30);
		failed |= check_close("raw first p11", p11, 0.40);
		failed |= check_close("raw first pU", pU, 0.25);
		failed |= check_close("raw first p4 sum", p00 + p01 + p10 + p11, 1.0);
		failed |= check_close("raw first five-state", rho * (p00 + p01 + p10 + p11) + pU, 1.0);
		failed |= check_close("raw first psame", psame, p00 + p11);
		failed |= check_close("raw first pcross", pcross, p01 + p10);
		failed |= check_close("raw first entropy", entropy, 1.20);
		failed |= check_close("raw first margin", margin, 0.10);
		failed |= check_close("raw first pmax", pmax, 0.40);
		line1 = strchr(line0, '\n');
		failed |= check_true("raw second line", line1 != 0);
		if (line1) {
			++line1;
			failed |= check_i32("raw second parsed",
								sscanf(line1, "%d\t%63s\t%d\t%63s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t"
									   "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
									   &raw_id, chr1, &pos1, chr2, &pos2, &bid1_raw, &bid2_raw,
									   &bpair_id, &bid1_can, &bid2_can, &swapped,
									   &p00, &p01, &p10, &p11, &pU, &psame, &pcross,
									   &entropy, &margin, &pmax, &rho), 22);
			failed |= check_i32("raw second id", raw_id, 1);
			failed |= check_i32("raw second bid1 raw", bid1_raw, 5);
			failed |= check_i32("raw second bid2 raw", bid2_raw, 2);
			failed |= check_i32("raw second swapped", swapped, 1);
			failed |= check_close("raw second p00", p00, 0.10);
			failed |= check_close("raw second p01 swapped", p01, 0.30);
			failed |= check_close("raw second p10 swapped", p10, 0.20);
			failed |= check_close("raw second p11", p11, 0.40);
			failed |= check_close("raw second p4 sum", p00 + p01 + p10 + p11, 1.0);
			failed |= check_close("raw second five-state", rho * (p00 + p01 + p10 + p11) + pU, 1.0);
		}
	}
	free(text);
	hk_blind_bpair_set_destroy(set);
	return failed;
}

static int test_coords_writer(void)
{
	struct hk_sdict dict;
	char *names[2] = { "chrA", "chrB" };
	int32_t len[2] = { 2000000, 1000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[3];
	fvec3_t coords[6];
	FILE *fp;
	char *text = 0, *line0, *line1;
	char chr[64];
	int st, en, bid, copy, diploid_bid;
	double x, y, z;
	int failed = 0;

	set_bmap(&bmap, &dict, names, len, beads);
	set_coords(coords);
	fp = tmpfile();
	failed |= check_true("coords tmpfile", fp != 0);
	if (fp == 0) return 1;
	failed |= check_i32("coords writer ret", hk_blind_write_diploid_coords_tsv(fp, &bmap, coords), 0);
	text = read_stream(fp);
	fclose(fp);
	failed |= check_true("coords text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("coords header copy", strstr(text, "copy") != 0);
	failed |= check_true("coords header diploid", strstr(text, "diploid_bid") != 0);
	failed |= check_i32("coords row count", count_lines(text) - 1, 2 * bmap.n_beads);
	failed |= check_no_forbidden_strings("coords", text);
	line0 = strchr(text, '\n');
	failed |= check_true("coords first data", line0 != 0);
	if (line0) {
		++line0;
		failed |= check_i32("coords first parsed",
							sscanf(line0, "%63s\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf",
								   chr, &st, &en, &bid, &copy, &diploid_bid, &x, &y, &z), 9);
		failed |= check_true("coords first chr", strcmp(chr, "chrA") == 0);
		failed |= check_i32("coords first start", st, 0);
		failed |= check_i32("coords first end", en, 1000000);
		failed |= check_i32("coords first bid", bid, 0);
		failed |= check_i32("coords first copy", copy, 0);
		failed |= check_i32("coords first diploid", diploid_bid, hk_diploid_bid(0, 0));
		failed |= check_true("coords first finite", isfinite(x) && isfinite(y) && isfinite(z));
		line1 = strchr(line0, '\n');
		failed |= check_true("coords second data", line1 != 0);
		if (line1) {
			++line1;
			failed |= check_i32("coords second parsed",
								sscanf(line1, "%63s\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf",
									   chr, &st, &en, &bid, &copy, &diploid_bid, &x, &y, &z), 9);
			failed |= check_i32("coords second start", st, 0);
			failed |= check_i32("coords second end", en, 1000000);
			failed |= check_i32("coords second bid", bid, 0);
			failed |= check_i32("coords second copy", copy, 1);
			failed |= check_i32("coords second diploid", diploid_bid, hk_diploid_bid(0, 1));
		}
	}
	free(text);
	return failed;
}

static int test_coords_writer_gz(void)
{
	struct hk_sdict dict;
	char *names[2] = { "chrA", "chrB" };
	int32_t len[2] = { 2000000, 1000000 };
	struct hk_bmap bmap;
	struct hk_bead beads[3];
	fvec3_t coords[6];
	char path[256];
	char *text = 0, *line;
	char chr[64];
	int st, en, bid, copy, diploid_bid;
	double x, y, z;
	int row_count = 0;
	int copy_mask[3] = {0, 0, 0};
	int failed = 0;

	set_bmap(&bmap, &dict, names, len, beads);
	bmap.n_beads = 2;
	set_coords(coords);
	snprintf(path, sizeof(path), "/tmp/hk_blind_coords_writer_gz_%ld.tsv.gz", (long)getpid());
	unlink(path);
	failed |= check_i32("coords gz writer ret",
						hk_blind_write_diploid_coords_tsv_gz(path, &bmap, coords), 0);
	text = read_gz_path(path);
	failed |= check_true("coords gz text", text != 0);
	if (text == 0) {
		unlink(path);
		return 1;
	}
	failed |= check_true("coords gz header",
						 strstr(text, "chr\tstart\tend\tbid\tcopy\tdiploid_bid\tx\ty\tz\n") == text);
	failed |= check_i32("coords gz row count", count_lines(text) - 1, 2 * bmap.n_beads);
	failed |= check_no_forbidden_strings("coords gz", text);
	line = strchr(text, '\n');
	failed |= check_true("coords gz first data", line != 0);
	if (line) {
		++line;
		while (*line) {
			char *next = strchr(line, '\n');
			int n_parsed;
			bid = copy = diploid_bid = -1;
			x = y = z = 0.0;
			n_parsed = sscanf(line, "%63s\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf",
							  chr, &st, &en, &bid, &copy, &diploid_bid, &x, &y, &z);
			failed |= check_i32("coords gz parsed", n_parsed, 9);
			failed |= check_true("coords gz bid valid", bid >= 0 && bid < bmap.n_beads);
			failed |= check_true("coords gz copy valid", copy == 0 || copy == 1);
			if (bid >= 0 && bid < bmap.n_beads && (copy == 0 || copy == 1)) {
				copy_mask[bid] |= 1 << copy;
				failed |= check_i32("coords gz diploid", diploid_bid, hk_diploid_bid(bid, copy));
			}
			failed |= check_true("coords gz finite", isfinite(x) && isfinite(y) && isfinite(z));
			++row_count;
			if (next == 0)
				break;
			line = next + 1;
		}
	}
	failed |= check_i32("coords gz parsed row count", row_count, 2 * bmap.n_beads);
	failed |= check_i32("coords gz copy mask 0", copy_mask[0], 3);
	failed |= check_i32("coords gz copy mask 1", copy_mask[1], 3);
	free(text);
	unlink(path);
	return failed;
}

static int test_loop_diag_writer(void)
{
	struct hk_blind_iter_loop_diag diag;
	FILE *fp;
	char *text = 0, *line, *header;
	int n_iter, n_completed, total_chr_flipped, n_bad_iter, n_relax_bad, n_coord_bad;
	double initial_entropy, final_entropy, initial_pU, final_pU;
	double initial_temperature, final_temperature, initial_rho, final_rho;
	double mean_sep, min_sep, max_sep, sum_k;
	double mean_rho_train, min_rho_train, max_rho_train;
	double rep_energy, rep_force;
	double refresh_temp, refresh_kl, refresh_switch, refresh_pU_before, refresh_pU_after;
	long long skipped_same_bin;
	int n_rep_bad, rep_mode, refreshed;
	int failed = 0;

	hk_blind_iter_loop_diag_init(&diag);
	diag.n_iter = 3;
	diag.n_completed = 3;
	diag.initial_mean_entropy = 1.3f;
	diag.final_mean_entropy = 1.1f;
	diag.initial_mean_pU = 0.9f;
	diag.final_mean_pU = 0.8f;
	diag.initial_temperature = 2.0f;
	diag.final_temperature = 1.0f;
	diag.initial_rho_train = 1.0f;
	diag.final_rho_train = 1.0f;
	diag.total_chr_flipped = 2;
	diag.n_bad_iter = 0;
	diag.n_relax_nonfinite_iter = 0;
	diag.n_coord_nonfinite = 0;
	diag.final_mean_sep = 1.5f;
	diag.final_min_sep = 1.0f;
	diag.final_max_sep = 2.0f;
	diag.final_sum_wedge_k = 17.25;
	diag.final_mean_rho_train_bpair = 0.75f;
	diag.final_min_rho_train_bpair = 0.5f;
	diag.final_max_rho_train_bpair = 1.0f;
	diag.final_n_skipped_same_bin_bpairs = 4;
	diag.final_repulsion_energy = 3.5f;
	diag.final_repulsion_force_l1 = 4.5f;
	diag.n_repulsion_nonfinite_step = 0;
	diag.repulsion_mode = HK_BLIND_REPULSION_CELL;
	diag.posterior_refreshed_after_final_relax = 1;
	diag.posterior_refresh_temperature = 1.0f;
	diag.posterior_refresh_mean_kl = 0.01;
	diag.posterior_refresh_top_state_switch_frac = 0.25f;
	diag.posterior_refresh_mean_pU_before = 0.8f;
	diag.posterior_refresh_mean_pU_after = 0.7f;

	fp = tmpfile();
	failed |= check_true("diag tmpfile", fp != 0);
	if (fp == 0) return 1;
	failed |= check_i32("diag writer ret", hk_blind_write_iter_loop_diag_tsv(fp, &diag), 0);
	text = read_stream(fp);
	fclose(fp);
	failed |= check_true("diag text", text != 0);
	if (text == 0) return 1;
	failed |= check_true("diag header n_iter", strstr(text, "n_iter") != 0);
	failed |= check_true("diag header final pU", strstr(text, "final_mean_pU") != 0);
	failed |= check_true("diag header sum k", strstr(text, "final_sum_wedge_k") != 0);
	failed |= check_true("diag header repulsion", strstr(text, "final_repulsion_energy") != 0);
	failed |= check_true("diag header repulsion mode", strstr(text, "repulsion_mode") != 0);
	failed |= check_i32("diag row count", count_lines(text) - 1, 1);
	failed |= check_no_forbidden_strings("diag", text);
	line = strchr(text, '\n');
	failed |= check_true("diag data line", line != 0);
	if (line) {
		*line = 0;
		header = text;
		++line;
		failed |= check_i32("diag n_iter parse", tsv_get_int(header, line, "n_iter", &n_iter), 0);
		failed |= check_i32("diag n_completed parse", tsv_get_int(header, line, "n_completed", &n_completed), 0);
		failed |= check_i32("diag initial entropy parse", tsv_get_double(header, line, "initial_mean_entropy", &initial_entropy), 0);
		failed |= check_i32("diag final entropy parse", tsv_get_double(header, line, "final_mean_entropy", &final_entropy), 0);
		failed |= check_i32("diag initial pU parse", tsv_get_double(header, line, "initial_mean_pU", &initial_pU), 0);
		failed |= check_i32("diag final pU parse", tsv_get_double(header, line, "final_mean_pU", &final_pU), 0);
		failed |= check_i32("diag initial temp parse", tsv_get_double(header, line, "initial_temperature", &initial_temperature), 0);
		failed |= check_i32("diag final temp parse", tsv_get_double(header, line, "final_temperature", &final_temperature), 0);
		failed |= check_i32("diag initial rho parse", tsv_get_double(header, line, "initial_rho_train", &initial_rho), 0);
		failed |= check_i32("diag final rho parse", tsv_get_double(header, line, "final_rho_train", &final_rho), 0);
		failed |= check_i32("diag flips parse", tsv_get_int(header, line, "total_chr_flipped", &total_chr_flipped), 0);
		failed |= check_i32("diag bad iter parse", tsv_get_int(header, line, "n_bad_iter", &n_bad_iter), 0);
		failed |= check_i32("diag relax bad parse", tsv_get_int(header, line, "n_relax_nonfinite_iter", &n_relax_bad), 0);
		failed |= check_i32("diag coord bad parse", tsv_get_int(header, line, "n_coord_nonfinite", &n_coord_bad), 0);
		failed |= check_i32("diag mean sep parse", tsv_get_double(header, line, "final_mean_sep", &mean_sep), 0);
		failed |= check_i32("diag min sep parse", tsv_get_double(header, line, "final_min_sep", &min_sep), 0);
		failed |= check_i32("diag max sep parse", tsv_get_double(header, line, "final_max_sep", &max_sep), 0);
		failed |= check_i32("diag sum k parse", tsv_get_double(header, line, "final_sum_wedge_k", &sum_k), 0);
		failed |= check_i32("diag mean rho train parse", tsv_get_double(header, line, "final_mean_rho_train_bpair", &mean_rho_train), 0);
		failed |= check_i32("diag min rho train parse", tsv_get_double(header, line, "final_min_rho_train_bpair", &min_rho_train), 0);
		failed |= check_i32("diag max rho train parse", tsv_get_double(header, line, "final_max_rho_train_bpair", &max_rho_train), 0);
		failed |= check_i32("diag skipped same-bin parse", tsv_get_ll(header, line, "final_n_skipped_same_bin_bpairs", &skipped_same_bin), 0);
		failed |= check_i32("diag rep energy parse", tsv_get_double(header, line, "final_repulsion_energy", &rep_energy), 0);
		failed |= check_i32("diag rep force parse", tsv_get_double(header, line, "final_repulsion_force_l1", &rep_force), 0);
		failed |= check_i32("diag rep bad parse", tsv_get_int(header, line, "n_repulsion_nonfinite_step", &n_rep_bad), 0);
		failed |= check_i32("diag rep mode parse", tsv_get_int(header, line, "repulsion_mode", &rep_mode), 0);
		failed |= check_i32("diag refreshed parse", tsv_get_int(header, line, "posterior_refreshed_after_final_relax", &refreshed), 0);
		failed |= check_i32("diag refresh temp parse", tsv_get_double(header, line, "posterior_refresh_temperature", &refresh_temp), 0);
		failed |= check_i32("diag refresh kl parse", tsv_get_double(header, line, "posterior_refresh_mean_kl", &refresh_kl), 0);
		failed |= check_i32("diag refresh switch parse", tsv_get_double(header, line, "posterior_refresh_top_state_switch_frac", &refresh_switch), 0);
		failed |= check_i32("diag refresh before parse", tsv_get_double(header, line, "posterior_refresh_mean_pU_before", &refresh_pU_before), 0);
		failed |= check_i32("diag refresh after parse", tsv_get_double(header, line, "posterior_refresh_mean_pU_after", &refresh_pU_after), 0);
		failed |= check_i32("diag n_iter", n_iter, 3);
		failed |= check_i32("diag n_completed", n_completed, 3);
		failed |= check_i32("diag n_bad_iter", n_bad_iter, 0);
		failed |= check_close("diag initial entropy", initial_entropy, 1.3);
		failed |= check_close("diag final entropy", final_entropy, 1.1);
		failed |= check_close("diag initial pU", initial_pU, 0.9);
		failed |= check_close("diag final pU", final_pU, 0.8);
		failed |= check_close("diag initial temperature", initial_temperature, 2.0);
		failed |= check_close("diag final temperature", final_temperature, 1.0);
		failed |= check_close("diag initial rho", initial_rho, 1.0);
		failed |= check_close("diag final rho", final_rho, 1.0);
		failed |= check_i32("diag total flips", total_chr_flipped, 2);
		failed |= check_i32("diag relax bad", n_relax_bad, 0);
		failed |= check_i32("diag coord bad", n_coord_bad, 0);
		failed |= check_close("diag mean sep", mean_sep, 1.5);
		failed |= check_close("diag min sep", min_sep, 1.0);
		failed |= check_close("diag max sep", max_sep, 2.0);
		failed |= check_close("diag sum k", sum_k, 17.25);
		failed |= check_close("diag mean rho train", mean_rho_train, 0.75);
		failed |= check_close("diag min rho train", min_rho_train, 0.5);
		failed |= check_close("diag max rho train", max_rho_train, 1.0);
		failed |= check_i32("diag skipped same-bin", (int)skipped_same_bin, 4);
		failed |= check_close("diag repulsion energy", rep_energy, 3.5);
		failed |= check_close("diag repulsion force", rep_force, 4.5);
		failed |= check_i32("diag repulsion bad", n_rep_bad, 0);
		failed |= check_i32("diag repulsion mode", rep_mode, HK_BLIND_REPULSION_CELL);
		failed |= check_i32("diag refreshed", refreshed, 1);
		failed |= check_close("diag refresh temp", refresh_temp, 1.0);
		failed |= check_close("diag refresh kl", refresh_kl, 0.01);
		failed |= check_close("diag refresh switch", refresh_switch, 0.25);
		failed |= check_close("diag refresh pU before", refresh_pU_before, 0.8);
		failed |= check_close("diag refresh pU after", refresh_pU_after, 0.7);
	}
	free(text);
	return failed;
}

int main(void)
{
	int failed = 0;
	failed |= test_posterior_writer();
	failed |= test_posterior_writer_contact_params();
	failed |= test_raw_contact_writer();
	failed |= test_coords_writer();
	failed |= test_coords_writer_gz();
	failed |= test_loop_diag_writer();
	return failed != 0;
}
