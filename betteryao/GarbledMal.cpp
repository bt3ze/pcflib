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


void GarbledMal::evl_next_gen_inp_com(const Bytes &row, size_t kx){
  Bytes out(this->m_gen_inp_decom[0].size());
  
  for (size_t jx = 0; jx < this->m_gen_inp_decom.size(); jx++)
    {
      if (row.get_ith_bit(jx)) { out ^= this->m_gen_inp_decom[jx]; }
    }
  
  byte bit = out.get_ith_bit(0);
  
  static Bytes tmp;
  
  Bytes::iterator it = this->m_in_bufr_ix + bit*Env::key_size_in_bytes();
  
  __m128i aes_key, aes_plaintext, aes_ciphertext, out_key;
  
  tmp.assign(out.begin(), out.begin()+Env::key_size_in_bytes());
  tmp.resize(16, 0);
  aes_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
  
  aes_plaintext = _mm_set1_epi64x((uint64_t)kx+10);
  
  KDF128((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)&aes_key);
  aes_ciphertext = _mm_and_si128(aes_ciphertext, this->m_clear_mask);
  
  tmp.assign(it, it+Env::key_size_in_bytes());
  tmp.resize(16, 0);
  out_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
  out_key = _mm_xor_si128(out_key, aes_ciphertext);
  
  bit = _mm_extract_epi8(out_key, 0) & 0x01;
  this->m_gen_inp_hash.set_ith_bit(kx, bit);
  
  this->m_in_bufr_ix += 2*Env::key_size_in_bytes();

}

void GarbledMal::gen_next_gen_inp_com(const Bytes &row, size_t kx){
  Bytes tmp;
  
  __m128i out_key[2];
  tmp = this->m_prng.rand_bits(Env::k());
  tmp.set_ith_bit(0, 0);
  tmp.resize(16, 0);
  out_key[0] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
  out_key[1] = _mm_xor_si128(out_key[0], this->m_R);
  
  assert(this->m_gen_inp_decom.size() % 2 == 0);
  
  Bytes msg(this->m_gen_inp_decom[0].size());
  for (size_t jx = 0; jx < this->m_gen_inp_decom.size()/2; jx++)
    {
      if (row.get_ith_bit(jx))
        {
          byte bit = this->m_gen_inp_mask.get_ith_bit(jx);
          msg ^= this->m_gen_inp_decom[2*jx+bit];
        }
    }
  
  __m128i in_key[2], aes_plaintext, aes_ciphertext;
  
  aes_plaintext = _mm_set1_epi64x((uint64_t)kx+10);
  
  tmp.assign(msg.begin(), msg.begin()+Env::key_size_in_bytes());
  tmp.resize(16, 0);
  in_key[0] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
  in_key[1] = _mm_xor_si128(in_key[0], this->m_R);
  
  KDF128((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)&in_key[0]);
  aes_ciphertext = _mm_and_si128(aes_ciphertext, this->m_clear_mask);
  out_key[0] = _mm_xor_si128(out_key[0], aes_ciphertext);

  KDF128((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)&in_key[1]);
  aes_ciphertext = _mm_and_si128(aes_ciphertext, this->m_clear_mask);
  out_key[1] = _mm_xor_si128(out_key[1], aes_ciphertext);

  const byte bit = msg.get_ith_bit(0);

  tmp.resize(16);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), out_key[  bit]);
  this->m_out_bufr.insert(this->m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
        
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), out_key[1-bit]);
  this->m_out_bufr.insert(this->m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());

}


inline bool GarbledMal::pass_check(){
  	assert(m_gen_inp_decom.size() == m_gen_inp_com.size());

	bool pass_chk = true;
	for (size_t ix = 0; ix < m_gen_inp_decom.size(); ix++)
	{
		pass_chk &= (m_gen_inp_decom[ix].hash(Env::k()) == m_gen_inp_com[ix]);
	}
	return pass_chk;
}


// virtual functions of GarbledBase to implement:


void GarbledMal::gen_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed) {

}    

void GarbledMal::evl_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp) {

}

void GarbledMal::initialize_circuit(){

}

void * GarbledMal::gen_next_gate(PCFState *st, PCFGate *current_gate){
  return 0;
}

void * GarbledMal::evl_next_gate(PCFState *st, PCFGate *current_gate){
  return 0;
}


#endif
