#ifndef GARBLED_CIRCUIT_H_
#define GARBLED_CIRCUIT_H_

#include <emmintrin.h>
#include <vector>

#include "Env.h"
#include "Prng.h"
#include "Hash.h"
#include "Aes.h"

#ifdef __CPLUSPLUS
extern "C" {
#endif

#include "../pcflib.h"

void *copy_key(void *);
void delete_key(void *);

void *gen_next_gate(struct PCFState *st, struct PCFGate *gate);
void *evl_next_gate(struct PCFState *st, struct PCFGate *gate);

#ifdef __CPLUSPLUS
}
#endif

#define  _mm_extract_epi8(x, imm) \
	((((imm) & 0x1) == 0) ?   \
	_mm_extract_epi16((x), (imm) >> 1) & 0xff : \
	_mm_extract_epi16( _mm_srli_epi16((x), 8), (imm) >> 1))


#endif
