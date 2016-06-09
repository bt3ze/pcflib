/*
  garbling.h 
  garbling functions written for speed
 */
#ifndef _garbling_h
#define _garbling_h

#include <time.h>
#include <emmintrin.h>
#include <vector>
#include "Aes.h"
#include "Bytes.h"
#include "../pcflib.h"



#define  _mm_extract_epi8(x, imm) \
        ((((imm) & 0x1) == 0) ?   \
	_mm_extract_epi16((x), (imm) >> 1) & 0xff : \
	_mm_extract_epi16( _mm_srli_epi16((x), 8), (imm) >> 1))





/**
   These functions serve as intermediaries
   to get back to the Garbled Circuit object
   since we go through such crazy calling loops with 
   the pcf interpreter struct
 */
void * gen_next_gate(PCFState *st, PCFGate *current_gate);
void * evl_next_gate(PCFState *st, PCFGate *current_gate);


// garbling accessories

void print128_num(__m128i var);
void Double(__m128i & key,const __m128i & clear_mask);


/**
   These functions compute the permutation on an input key
   it is destructive of key so must make copies of the inputs first
   it returns H(K) = pi(L) xor L where L = 2key ^ tweak
 */
void H_Pi(__m128i & destination, __m128i &key, __m128i & tweak,
            const __m128i & clear_mask,
            const AES_KEY_J & fixed_key);

void H_Pi256(__m128i & destination, __m128i &key1, __m128i &key2,
               __m128i & tweak, const __m128i & clear_mask,
               const AES_KEY_J & fixed_key);

void xor_Gate(__m128i & key1, __m128i & key2, __m128i &current_key);


// half gate declarations
void genHalfGatePair(__m128i& out_key, __m128i & key1, __m128i & key2, Bytes & out_bufr,
                       const byte a1, const byte a2,const byte a3,const size_t keysize,
                       const __m128i & clear_mask, const AES_KEY_J & fixed_key ,
                       const __m128i & R, const uint32_t j1, const uint32_t j2);

void evlHalfGatePair(__m128i& cur_key, __m128i & key1, __m128i & key2, Bytes & in_bufr,
                       const size_t keysize, const __m128i & clear_mask,
                       const AES_KEY_J & fixed_key, const uint32_t j1,const uint32_t j2);

// standard gate declarations
void genStandardGate(__m128i& current_key, __m128i & key1, __m128i & key2,
                       Bytes & out_bufr, uint8_t truth_table, const size_t keysize,
                       const __m128i & clear_mask, const AES_KEY_J & fixed_key,
                       const __m128i & R, const uint32_t j1);

void evlStandardGate(__m128i& current_key, __m128i & key1, __m128i & key2,
                       Bytes & in_bufr, const size_t keysize,
                       const __m128i & clear_mask, const AES_KEY_J & fixed_key,
                       const uint32_t j1);



#endif
