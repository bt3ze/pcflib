#ifndef GARBLED_CIRCUIT_HBC_H_
#define GARBLED_CIRCUIT_HBC_H_

#include "garbled_circuit.h"

typedef struct
{
  Bytes               m_bufr;
  Hash                m_hash;

  __m128i             m_R; // constant for free-XOR
  
  const std::vector<Bytes>  *m_ot_keys; // pointer to keys
  
  Prng                m_prng; // generator
  
  uint64_t            m_gate_ix;

  // input indices
  uint32_t            m_gen_inp_ix;
  uint32_t            m_evl_inp_ix;
  uint32_t            m_gen_out_ix;
  uint32_t            m_evl_out_ix;
  
  __m128i             m_clear_mask;

  Bytes               m_gen_inp_mask;
  Bytes               m_gen_inp;
  Bytes               m_evl_inp;
  Bytes               m_gen_out;
  Bytes               m_evl_out;
  
  Bytes               m_out_bufr; // in buffer
  Bytes               m_in_bufr; // out buffer
  Bytes::iterator     m_in_bufr_ix;
  
  struct PCFState    *m_st; // pointer to the PCF state
  __m128i             m_const_wire[2]; // keys for constant 0 and 1

}
garbled_circuit_t;


// prototypes with analogous versions in malicious garbled circuit

void gen_init_circuit(garbled_circuit_t &cct, const std::vector<Bytes> &keys, const Bytes &gen_inp_mask, const Bytes &seed);
void evl_init_circuit(garbled_circuit_t &cct, const std::vector<Bytes> &keys, const Bytes &masked_gen_inp, const Bytes &seed);

inline void trim_output(garbled_circuit_t &cct)
{
	cct.m_gen_out.resize((cct.m_gen_out_ix+7)/8);
	cct.m_evl_out.resize((cct.m_evl_out_ix+7)/8);
}

inline void clear_and_replace_in_bufr(garbled_circuit_t &cct, const Bytes &i_data)
{
	cct.m_in_bufr.clear();
	cct.m_in_bufr += i_data;
        //assert(cct.m_in_bufr.size() > 0);
	cct.m_in_bufr_ix = cct.m_in_bufr.begin();
}

inline const Bytes get_and_clear_out_bufr(garbled_circuit_t &cct)
{
	static Bytes o_data;
	o_data.swap(cct.m_out_bufr);
	cct.m_out_bufr.clear();
	return o_data;
}


void set_const_key(garbled_circuit_t &cct, byte c, const Bytes &key);
const Bytes &get_const_key(garbled_circuit_t &cct, byte c, byte b);


#endif
