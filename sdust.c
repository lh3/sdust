#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "kalloc.h"
#include "khashl.h"
#include "kdq.h"
#include "sdust.h"

KHASHL_MAP_INIT(KH_LOCAL, sd_knt_t, sd_knt, uint32_t, int32_t, kh_hash_uint32, kh_eq_generic)

typedef struct {
	int64_t st, en;
	int64_t r, l;
} perf_intv_t;

KDQ_INIT(uint32_t)

static unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

struct sdust_s {
	int32_t kmer, ws;
	double thres;
	kdq_t(uint32_t) *w;
	int64_t n_perf, m_perf;
	perf_intv_t *perf; // the list of perfect intervals for the current window, sorted by descending start and then by ascending finish
	int64_t n_out, m_out;
	sdust_intv_t *out;
	sd_knt_t *cw, *cv;
	void *km;
};

sdust_t *sdust_init(void *km, int32_t kmer, int32_t ws, double thres)
{
	sdust_t *sd;
	sd = Kcalloc(km, sdust_t, 1);
	sd->km = km;
	sd->kmer = kmer, sd->ws = ws, sd->thres = thres;
	sd->w = kdq_init(uint32_t, sd->km);
	kdq_resize(uint32_t, sd->w, 8);
	sd->cw = sd_knt_init2(km);
	sd->cv = sd_knt_init2(km);
	return sd;
}

void sdust_destroy(sdust_t *sd)
{
	if (sd == 0) return;
	kdq_destroy(uint32_t, sd->w);
	sd_knt_destroy(sd->cv);
	sd_knt_destroy(sd->cw);
	kfree(sd->km, sd->perf);
	kfree(sd->km, sd->out);
	kfree(sd->km, sd);
}

static inline int32_t sd_knt_dec(sd_knt_t *h, uint32_t t)
{
	int32_t ret;
	khint_t k;
	k = sd_knt_get(h, t);
	if (k == kh_end(h)) return -1;
	ret = kh_val(h, k) - 1;
	if (ret == 1)
		sd_knt_del(h, k);
	return ret;
}

static inline int32_t sd_knt_inc(sd_knt_t *h, uint32_t t)
{
	int32_t ret, absent;
	khint_t k;
	k = sd_knt_put(h, t, &absent);
	if (absent)
		kh_val(h, k) = 0;
	ret = kh_val(h, k);
	++kh_val(h, k);
	return ret;
}

static inline int32_t sd_knt_cnt(sd_knt_t *h, uint32_t t)
{
	khint_t k;
	k = sd_knt_get(h, t);
	return k == kh_end(h)? 0 : kh_val(h, k);
}

static inline void shift_window(sdust_t *sd, uint32_t t, int *L, int *rw, int *rv)
{
	uint32_t s;
	if (kdq_size(sd->w) >= sd->ws - sd->kmer + 1) {
		s = *kdq_shift(uint32_t, sd->w);
		*rw -= sd_knt_dec(sd->cw, s);
		if (*L > kdq_size(sd->w))
			--*L, *rv -= sd_knt_dec(sd->cv, s);
	}
	kdq_push(uint32_t, sd->w, t);
	++*L;
	*rw += sd_knt_inc(sd->cw, t);
	*rv += sd_knt_inc(sd->cv, t);
	if (sd_knt_cnt(sd->cv, t) > 2.0 * sd->thres) {
		do {
			s = kdq_at(sd->w, kdq_size(sd->w) - *L);
			*rv -= sd_knt_dec(sd->cv, s);
			--*L;
		} while (s != t);
	}
}

static inline void save_masked_regions(sdust_t *sd, int64_t st)
{
	int saved = 0;
	int64_t i;
	const perf_intv_t *p;

	if (sd->n_perf == 0 || sd->perf[sd->n_perf - 1].st >= st) return;
	p = &sd->perf[sd->n_perf - 1];
	if (sd->n_out > 0) {
		sdust_intv_t *q = &sd->out[sd->n_out - 1];
		if (p->st <= q->en) // if overlapping with or adjacent to the previous interval
			saved = 1, q->en = q->en > p->en? q->en : p->en;
	}
	if (!saved) {
		Kgrow(sd->km, sdust_intv_t, sd->out, sd->n_out, sd->m_out);
		sd->out[sd->n_out].st = p->st;
		sd->out[sd->n_out++].en = p->en;
	}
	for (i = sd->n_perf - 1; i >= 0 && sd->perf[i].st < st; --i) {} // remove perfect intervals that have falled out of the window
	sd->n_perf = i + 1;
}

static void find_perfect(sdust_t *sd, int64_t st, int L, int rv)
{
	sd_knt_t *c;
	int r = rv, i, max_r = 0, max_l = 0;
	c = sd_knt_init2(sd->km);
	for (i = kdq_size(sd->w) - L - 1; i >= 0; --i) {
		uint32_t t = kdq_at(sd->w, i);
		int j, new_r, new_l;
		r += sd_knt_inc(c, t) + sd_knt_cnt(sd->cv, t);
		new_r = r, new_l = kdq_size(sd->w) - i - 1;
		if (new_r > sd->thres * new_l) {
			for (j = 0; j < sd->n_perf && sd->perf[j].st >= i + st; ++j) { // find insertion position
				const perf_intv_t *p = &sd->perf[j];
				if (max_r == 0 || p->r * max_l > max_r * p->l)
					max_r = p->r, max_l = p->l;
			}
			if (max_r == 0 || new_r * max_l >= max_r * new_l) { // then insert
				max_r = new_r, max_l = new_l;
				Kgrow(sd->km, perf_intv_t, sd->perf, sd->n_perf, sd->m_perf);
				memmove(&sd->perf[j+1], &sd->perf[j], (sd->n_perf - j) * sizeof(perf_intv_t)); // make room
				++sd->n_perf;
				sd->perf[j].st = i + st, sd->perf[j].en = kdq_size(sd->w) + (sd->kmer - 1) + st;
				sd->perf[j].r = new_r, sd->perf[j].l = new_l;
			}
		}
	}
	sd_knt_destroy(c);
}

const sdust_intv_t *sdust_seq(sdust_t *sd, const uint8_t *seq, int64_t l_seq, int64_t *n)
{
	int rv = 0, rw = 0, L = 0;
	int64_t i, st, l; // _st_: start of the current window; _l_: length of a contiguous A/C/G/T (sub)sequence
	uint32_t t, mask = (1U<<2*sd->kmer) - 1; // current word

	sd->n_perf = sd->n_out = 0;
	sd->w->front = sd->w->count = 0;
	sd_knt_clear(sd->cv);
	sd_knt_clear(sd->cw);
	if (l_seq < 0) l_seq = strlen((const char*)seq);
	for (i = l = t = 0; i <= l_seq; ++i) {
		int b = i < l_seq? seq_nt4_table[seq[i]] : 4;
		if (b < 4) { // an A/C/G/T base
			++l, t = (t<<2 | b) & mask;
			if (l >= sd->kmer) { // we have seen a word
				st = (l - sd->ws > 0? l - sd->ws : 0) + (i + 1 - l); // set the start of the current window
				save_masked_regions(sd, st); // save intervals falling out of the current window?
				shift_window(sd, t, &L, &rw, &rv);
				if (rw > L * sd->thres)
					find_perfect(sd, st, L, rv);
			}
		} else { // N or the end of sequence; N effectively breaks input into pieces of independent sequences
			st = (l - sd->ws + 1 > 0? l - sd->ws + 1 : 0) + (i + 1 - l);
			while (sd->n_perf > 0) save_masked_regions(sd, st++); // clear up unsaved perfect intervals
			l = t = 0;
		}
	}
	*n = sd->n_out;
	return sd->out;
}

#include <zlib.h>
#include <stdio.h>
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *ks;
	int ws = 64, T = 20, kmer = 3, print_len = 0, c;
	ketopt_t o = KETOPT_INIT;
	sdust_t *sd;

	while ((c = ketopt(&o, argc, argv, 1, "w:t:k:l", 0)) >= 0) {
		if (c == 'w') ws = atoi(o.arg);
		else if (c == 't') T = atoi(o.arg);
		else if (c == 'k') kmer = atoi(o.arg);
		else if (c == 'l') print_len = 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: sdust [options] <in.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT      k-mer size [%d]\n", kmer);
		fprintf(stderr, "  -w INT      window size [%d]\n", ws);
		fprintf(stderr, "  -t INT      threshold*10 [%d]\n", T);
		fprintf(stderr, "  -l          print contig length\n");
		return 1;
	}
	sd = sdust_init(0, kmer, ws, (double)T * 0.1);
	fp = strcmp(argv[o.ind], "-")? gzopen(argv[o.ind], "r") : gzdopen(fileno(stdin), "r");
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		const sdust_intv_t *r;
		int64_t i, n;
		r = sdust_seq(sd, (uint8_t*)ks->seq.s, ks->seq.l, &n);
		for (i = 0; i < n; ++i) {
			if (print_len)
				printf("%s\t%ld\t%ld\t%ld\n", ks->name.s, (long)r[i].st, (long)r[i].en, (long)ks->seq.l);
			else
				printf("%s\t%ld\t%ld\n", ks->name.s, (long)r[i].st, (long)r[i].en);
		}
	}
	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}
