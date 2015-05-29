#ifndef GARBLED_HBC_CPP_
#define GARBLED_HBC_CPP_

#include "GarbledBase.h"
#include "GarbledHBC.h"

GarbledHBC::GarbledHBC(): GarbledBase() {
  
}


// virtual functions of GarbledBase to implement:


void GarbledHBC::gen_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed){
  m_ot_keys = &ot_keys;
  m_gen_inp_mask = gen_inp_mask;
  m_prng.seed_rand(seed);

  Bytes tmp;

  tmp = m_prng.rand_bits(Env::k());
  tmp.set_ith_bit(0, 1);
  tmp.resize(16, 0);
  m_R = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));
  
  // pick zero-keys for constant wires
  // filling up extra bits with 0
  tmp = m_prng.rand_bits(Env::k());
  tmp.resize(16, 0);
  m_const_wire[0] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));
  
  tmp = m_prng.rand_bits(Env::k());
  tmp.resize(16, 0);
  m_const_wire[1] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));
  
  
}


void GarbledHBC::evl_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp){
  m_ot_keys = &ot_keys;
  m_gen_inp_mask = masked_gen_inp;
  m_evl_inp = evl_inp;

  m_evl_out.reserve(MAX_OUTPUT_SIZE);
  m_evl_out.clear(); // will grow dynamically
  m_gen_out.reserve(MAX_OUTPUT_SIZE);
  m_gen_out.clear(); // will grow dynamically

  //m_bufr.reserve(CIRCUIT_HASH_BUFFER_SIZE);
  //m_bufr.clear();
  
}

#endif
