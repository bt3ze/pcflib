#ifndef GARBLED_CIRCUIT_CPP
#define GARBLED_CIRCUIT_CPP

#include "GarbledCircuit.h"


GarbledCircuit::GarbledCircuit(): m_gate_index(0) {

  // initialize the key Mask
  Bytes tmp(16);
  for(size_t ix = 0; ix< Env::k(); ix++){
    tmp.set_ith_bit(ix,1);
  }
  m_clear_mask = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));

}

void GarbledCircuit::init_Generation_Circuit(const std::vector<Bytes>* gen_keys, const std::vector<Bytes> * evl_keys){
  Prng prng = Prng();
  Bytes tmp;
  
  // initialize the value of R
  tmp = prng.rand_bits(Env::k());
  tmp.resize(16,0);
  tmp.set_ith_bit(0,1);
  m_R = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  // initialize the constant wires
  tmp = prng.rand_bits(Env::k());
  tmp.resize(16,0);
  m_const_wire[0] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  tmp = prng.rand_bits(Env::k());
  tmp.resize(16,0);
  m_const_wire[1] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  // set the input keys
  gen_inputs = gen_keys;
  evl_inputs = evl_keys;

}


void GarbledCircuit::init_Evaluation_Circuit(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys){
  
  gen_inputs = gen_keys;
  evl_inputs = evl_keys;

}

bool GarbledCircuit::gen_Next_Gate(){

}

bool GarbledCircuit::evl_Next_Gate(){

}


void GarbledCircuit::set_const_key(byte c, const Bytes &key)
{
  assert(c == 0 || c == 1); // wire for constant 0 or 1
  Bytes tmp = key;
  tmp.resize(16);
  m_const_wire[c] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
}


const Bytes GarbledCircuit::get_const_key(byte c, byte b)
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
