#ifndef GARBLED_MAL_CPP
#define GARBLED_MAL_CPP

#include "GarbledBase.h"
#include "GarbledMal.h"

GarbledMal::GarbledMal() : GarbledBase() {
    
  // clear out bufr
  m_out_bufr.clear();
  
  // input hash set to zeros
  m_gen_inp_hash.assign(Env::key_size_in_bytes(), 0);
  
  Bytes tmp(16);
  for (size_t ix = 0; ix < Env::k(); ix++)
    { tmp.set_ith_bit(ix, 1); }
  // load the clear mask with ones for as many digits as K
  m_clear_mask = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));

}


void * GarbledMal::gen_next_gate(){
  return 0;
}

void * GarbledMal::evl_next_gate(){
  return 0;
}

#endif
