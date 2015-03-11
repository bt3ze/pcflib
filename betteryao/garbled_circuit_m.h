#ifndef GARBLED_CIRCUIT_M_H_
#define GARBLED_CIRCUIT_M_H_

#include "garbled_circuit.h"

typedef struct
{
  Bytes               m_bufr;
  Hash                m_hash;

  __m128i             m_R; // constant used for free-XOR

  const std::vector<Bytes>  *m_ot_keys; // pointer to constant ot keys
  
  Prng                m_prng; // pseudorandom number generator
  
  uint64_t            m_gate_ix; // gate index?

  // input indices
  uint32_t            m_gen_inp_ix;
  uint32_t            m_evl_inp_ix;
  uint32_t            m_gen_out_ix;
  uint32_t            m_evl_out_ix;

  __m128i             m_clear_mask; // what is the purpose of this mask?

  // gen and eval inputs
  Bytes               m_gen_inp_mask;
  Bytes               m_gen_inp;
  Bytes               m_evl_inp;
  Bytes               m_gen_out;
  Bytes               m_evl_out;

  Bytes               m_o_bufr; // out buffer
  Bytes               m_i_bufr; // in buffer
  Bytes::iterator     m_i_bufr_ix;

  struct PCFState    *m_st; // pointer to the PCF state
  __m128i             m_const_wire[2]; // keys for constant 0 and 1


  // specific to malicious circuit
  std::vector<Bytes>  m_gen_inp_com; // commitments?
  std::vector<Bytes>  m_gen_inp_decom; // decommitments?
  Bytes               m_gen_inp_hash;
    
} garbled_circuit_m_t;


// prototypes with analogous versions in hbc garbled circuit


void gen_init_circuit(garbled_circuit_m_t &cct, const std::vector<Bytes> &keys, const Bytes &gen_inp_mask, const Bytes &seed);
void evl_init_circuit(garbled_circuit_m_t &cct, const std::vector<Bytes> &keys, const Bytes &masked_gen_inp, const Bytes &seed);


inline void trim_output(garbled_circuit_m_t &cct)
{
	cct.m_gen_out.resize((cct.m_gen_out_ix+7)/8);
	cct.m_evl_out.resize((cct.m_evl_out_ix+7)/8);
}


inline void recv(garbled_circuit_m_t &cct, const Bytes &i_data)
{
	cct.m_i_bufr.clear();
	cct.m_i_bufr += i_data;
        //assert(cct.m_i_bufr.size() > 0);
	cct.m_i_bufr_ix = cct.m_i_bufr.begin();
}


inline const Bytes get_and_clear_out_bufr(garbled_circuit_m_t &cct){
  static Bytes o_data;
  o_data.swap(cct.m_o_bufr);
  cct.m_o_bufr.clear();
  return o_data;
}


void set_const_key(garbled_circuit_m_t &cct, byte c, const Bytes &key);
const Bytes get_const_key(garbled_circuit_m_t &cct, byte c, byte b);
// flag the above &


// now, prototypes and functions specific to malicious

void * evl_next_gate_m(struct PCFState *st, struct PCFGate *current_gate);
void * gen_next_gate_m(struct PCFState *st, struct PCFGate *current_gate);
void evl_next_gen_inp_com(garbled_circuit_m_t &cct, const Bytes &row, size_t kx);
void gen_next_gen_inp_com(garbled_circuit_m_t &cct, const Bytes &row, size_t kx);

inline bool pass_check(const garbled_circuit_m_t &cct)
{
	assert(cct.m_gen_inp_decom.size() == cct.m_gen_inp_com.size());

	bool pass_chk = true;
	for (size_t ix = 0; ix < cct.m_gen_inp_decom.size(); ix++)
	{
		pass_chk &= (cct.m_gen_inp_decom[ix].hash(Env::k()) == cct.m_gen_inp_com[ix]);
	}
	return pass_chk;
}


/**
   set all gate and input indices to 0
   clear the out buffer
   fill the inp hash with as many 0s as the security param
   set <security param> bits in the clear mask to 1 
 */
inline void initialize_circuit_mal(garbled_circuit_m_t &cct)
{
	cct.m_gate_ix = 0;

	cct.m_gen_inp_ix = 0;
	cct.m_evl_inp_ix = 0;
	cct.m_gen_out_ix = 0;
	cct.m_evl_out_ix = 0;

	cct.m_o_bufr.clear();

	cct.m_gen_inp_hash.assign(Env::key_size_in_bytes(), 0);

	Bytes tmp(16);
	for (size_t ix = 0; ix < Env::k(); ix++)
          { tmp.set_ith_bit(ix, 1); }
	cct.m_clear_mask = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
}


#endif
