#ifndef GARBLED_BASE_CPP_
#define GARBLED_BASE_CPP_

#include "Bytes.h"
#include "GarbledBase.h"


GarbledBase::GarbledBase(){
  // m_gate_ix(0), m_gen_inp_ix(0), m_evl_inp_ix(0), m_gen_out_ix(0), m_evl_out_ix(0) {
  /*
  assert(m_gate_ix == 0);
  assert(m_gen_inp_ix == 0);
  assert(m_evl_inp_ix == 0);
  assert(m_gen_out_ix == 0);
  assert(m_evl_out_ix == 0);
  */
  // set clear mask to as many 1s as k
  Bytes tmp(16);

  for (size_t ix = 0; ix < Env::k(); ix++) {
    tmp.set_ith_bit(ix, 1);
  }
  m_clear_mask = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));

  // clear in and out buffers
  m_in_bufr.clear();
  m_out_bufr.clear();
  

}


inline void GarbledBase::trim_output()
{
  assert(m_gen_out.size() > 0);
  m_gen_out.resize((m_gen_out_ix+7)/8);
  m_evl_out.resize((m_evl_out_ix+7)/8);
}

inline void GarbledBase::clear_and_replace_in_bufr(const Bytes &i_data)
{
  m_in_bufr.clear();
  m_in_bufr += i_data;
  //assert(cct.m_i_bufr.size() > 0);
  m_in_bufr_ix = m_in_bufr.begin();
}

inline const Bytes GarbledBase::get_and_clear_out_bufr()
{
  static Bytes o_data;
  o_data.swap(m_out_bufr);
  m_out_bufr.clear();
  return o_data;
}


void GarbledBase::set_const_key(byte c, const Bytes &key)
{
	assert(c == 0 || c == 1); // wire for constant 0 or 1
	Bytes tmp = key;
	tmp.resize(16);
	m_const_wire[c] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
}


const Bytes GarbledBase::get_const_key(byte c, byte b)
{
	assert(c == 0 || c == 1); // wire for constant 0 or 1
	assert(b == 0 || b == 1); // with bit value 0 or 1
	Bytes tmp(16);

	tmp.resize(16);
	if (b)
	{
		_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), _mm_xor_si128(m_R, m_const_wire[c]));
	}
	else
	{
		_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), m_const_wire[c]);
	}

	tmp.resize(Env::key_size_in_bytes());
	return tmp;
}





#endif
