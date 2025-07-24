#ifndef SDUST_H
#define SDUST_H

#define SDUST_VERSION "0.1-r2"

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct sdust_s;
typedef struct sdust_s sdust_t;

typedef struct {
	int64_t st, en;
} sdust_intv_t;

sdust_t *sdust_init(void *km, int32_t kmer, int32_t ws, double thres);
void sdust_destroy(sdust_t *buf);
const sdust_intv_t *sdust_seq(sdust_t *sd, const uint8_t *seq, int64_t l_seq, int64_t *n);

#ifdef __cplusplus
}
#endif

#endif
