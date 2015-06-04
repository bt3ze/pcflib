#ifndef GARBLED_MAL_CPP
#define GARBLED_MAL_CPP

#include "GarbledBase.h"
#include "GarbledMal.h"
#include <stdio.h>
#include <iostream>

GarbledMal::GarbledMal() : GarbledBase() {
  // input, output, and gate indices set to 0 by parent constructor
  // input and output buffers cleared
  // clear mask set to length of security parameter
  
  // set the input hash to zeros
  // (do we actually use the input hash?)
  m_gen_inp_hash.assign(Env::key_size_in_bytes(), 0);
  
  // clear the input commitment and decommitment vectors
  m_gen_inp_com.clear();
  m_gen_inp_decom.clear();
  
  
  
}

// virtual functions of GarbledBase to implement:

// to be replaced with the function it calls
void GarbledMal::gen_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed) {
  this->initialize_gen_circuit(ot_keys, gen_inp_mask, seed);
}    


/**
   To initialize an evaluation/garbling circuit, we must
   - set its OT keys,
   - generate the free-XOR offset R,
   - and seed its PRNG.
   - and set Gen's masked input -- why?
   Gen calls this function to initialize all of his circuits.
   Eval calls this function to initialize her check circuits.
 */
void GarbledMal::initialize_gen_circuit(const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed){
  
  // set the circuit's ot keys
  m_ot_keys = &ot_keys;

  // set the circuit's input mask (is this necessary for new version (BY5)?)
  m_gen_inp_mask = gen_inp_mask;

  // must seed the circuit's random number generator
  m_prng.seed_rand(seed);

  // at this point (now that we have a prng),
  //we can determine the free-XOR value R for the circuit
  Bytes tmp;
  tmp = m_prng.rand_bits(Env::k());
  tmp.set_ith_bit(0,1);
  tmp.resize(16,0);
  m_R = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  // must also pick zero-keys for the constant wires
  tmp = m_prng.rand_bits(Env::k());
  tmp.resize(16, 0);
  m_const_wire[0] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  tmp = m_prng.rand_bits(Env::k());
  tmp.resize(16, 0);
  m_const_wire[1] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

}

// to be replaced with the function it calls
void GarbledMal::evl_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp) {
  this->initialize_eval_circuit(ot_keys, masked_gen_inp, evl_inp);
}

/**
   To initialize an evaluation circuit, we need to
   - set the OT keys
   - set Eval's input
   - set Gen's masked input (why?)
   - and clear the circuit's output buffers
   Eval calls this function for her evaluation circuits.
 */
void GarbledMal::initialize_eval_circuit(const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp){
  m_ot_keys = &ot_keys;
  m_gen_inp_mask = masked_gen_inp; // do we really need to mask Gen's input?
  m_evl_inp = evl_inp;


  // Eval has to figure out the outputs for both Gen and Eval,
  // so the output buffers are reserved and cleared here in preparation
  m_evl_out.reserve(MAX_OUTPUT_SIZE); // MAX_OUTPUT_SIZE in GarbledBase.h
  m_evl_out.clear(); // will grow dynamically
  m_gen_out.reserve(MAX_OUTPUT_SIZE);
  m_gen_out.clear(); // will grow dynamically

}


/**
   Evaluate and Generate Input Commitment Functions
 */

void GarbledMal::evl_next_gen_inp_com(const Bytes &row, size_t kx){
  // get Bytes array the length of a decommitment
  // it defaults to 0
  Bytes out(this->m_gen_inp_decom[0].size());
  
  std::cout << "m_gen_inp_decom.size(): " << this->m_gen_inp_decom.size() << std::endl;
  std::cout << "Out (bytes): " << out.to_hex() << std::endl;
  // this takes the XOR of all of Gen's input decommitments for which
  // the ith value in the 2-UHF row is set
  // this is equivalent to evaluating a circuit of only
  // XOR gates that evaluates the 2-UHF
  // since every time a bit is set in the UHF,
  // we XOR the current key with all of the previous
  // and we end up with the output of the 2-UHF
  // below, we will extract the semantics of the final bit
  for (size_t jx = 0; jx < this->m_gen_inp_decom.size(); jx++)
    {
      if (row.get_ith_bit(jx)) { out ^= this->m_gen_inp_decom[jx];
        std::cout << "Out (bytes): " << out.to_hex() << std::endl;
      }
    }
  
  std::cout << "Out (bytes): " << out.to_hex() << std::endl;
  
  // now get the 0th bit
  // this should be the permutation bit
  // and since it is an output bit,
  // it carries the semantics of the wire
  byte bit = out.get_ith_bit(0);
  
  //  std::cout << "bit: " << (bit == 0 ? "0" : "1") << std::endl;
  //  assert(bit == 0 || bit == 1);

  //static Bytes tmp;
  Bytes tmp;

  // this piece is a little troubling: where is m_in_bufr_ix?
  //  std::cout << "m_in_bufr_ix: " << (this->m_in_bufr_ix) << std::endl;
  Bytes::iterator it = this->m_in_bufr_ix + bit*Env::key_size_in_bytes();
  
  __m128i aes_key, aes_plaintext, aes_ciphertext, out_key;
  
  tmp.assign(out.begin(), out.begin()+Env::key_size_in_bytes());
  tmp.resize(16, 0);
  aes_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));

  std::cout << "aes key: " << tmp.to_hex() << std::endl;
  aes_plaintext = _mm_set1_epi64x((uint64_t)kx+10); // why +10??
  
  KDF128((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)&aes_key);
  aes_ciphertext = _mm_and_si128(aes_ciphertext, this->m_clear_mask);
  
  tmp.assign(it, it+Env::key_size_in_bytes());
  tmp.resize(16, 0);
  out_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
  out_key = _mm_xor_si128(out_key, aes_ciphertext);
  
  std::cout << "out key: " << tmp.to_hex() << std::endl;

  bit = _mm_extract_epi8(out_key, 0) & 0x01;

  std::cout << "hash bit: " << (bit == 0? "0" : "1") << std::endl;
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
  std::cout << "msg (1): " << msg.to_hex() << std::endl;
  for (size_t jx = 0; jx < this->m_gen_inp_decom.size()/2; jx++)
    {
      if (row.get_ith_bit(jx))
        {
          byte bit = this->m_gen_inp_mask.get_ith_bit(jx);
          msg ^= this->m_gen_inp_decom[2*jx+bit];
        }
    }
  std::cout << "msg (2): " << msg.to_hex() << std::endl;
  
  

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



/**
   These are the gen_next_gate functions
   polynomial interpolation is just out of scope of these functions
   they simply generate the next gate for the individual circuit
   that they represent
 */

void * gen_next_malicious_gate(PCFState *st, PCFGate *current_gate){
  //  garbled_circuit_m_t &cct =
  //  *reinterpret_cast<garbled_circuit_m_t*>(get_external_circuit(st));
  GarbledMal &cct =
    *reinterpret_cast<GarbledMal*>(get_external_circuit(st));
  
  
  static __m128i current_zero_key;
  
  if (current_gate->tag == TAG_INPUT_A) // this is a Gen input
    {
      __m128i a[2];
      
      Bytes tmp = cct.m_prng.rand_bits(Env::k());
      tmp.resize(16, 0);
      current_zero_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
      
      uint32_t gen_inp_ix = current_gate->wire1;
      
      a[0] = current_zero_key;
      a[1] = _mm_xor_si128(current_zero_key, cct.m_R);
      
      if(!(gen_inp_ix < cct.m_gen_inp_mask.size()*8)){
        fprintf(stderr,"error: gen inp index not less than m_gen_inp_mask size * 8: %u !< %lu*8 (%lu)\n",gen_inp_ix, cct.m_gen_inp_mask.size(), cct.m_gen_inp_mask.size()*8);
      }
      uint8_t bit = cct.m_gen_inp_mask.get_ith_bit(gen_inp_ix);
      
      // std::cout <<"dcomsize: "<<cct.m_gen_inp_decom.size()<<"  input ix: "<< gen_inp_ix<<"\n";
      assert(cct.m_gen_inp_decom.size() == 2*cct.m_gen_inp_ix);
      
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), a[bit]);
      cct.m_gen_inp_decom
        .push_back(Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes())
                   +cct.m_prng.rand_bits(Env::k()));
      
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), a[1-bit]);
      cct.m_gen_inp_decom
        .push_back(Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes())
                   +cct.m_prng.rand_bits(Env::k()));
      
      cct.m_out_bufr += cct.m_gen_inp_decom[2*cct.m_gen_inp_ix+0].hash(Env::k());
      cct.m_out_bufr += cct.m_gen_inp_decom[2*cct.m_gen_inp_ix+1].hash(Env::k());
      
      //cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
      
      
      cct.m_gen_inp_ix++; // after PCF compiler, this isn't really necessary
      
      
      __m128i onev = _mm_xor_si128(cct.m_R, current_zero_key);
      
      /*
        if(Env::is_root())
        {
        std::cout <<"GENbuffr: "<<cct.m_out_bufr.to_hex()<<"\n";		
        std::cout <<"GEN: gate: "<<cct.m_gate_ix<<" : "<< setBytes(current_zero_key).to_hex()<<" "<< setBytes(onev).to_hex() <<"  "<<current_gate->tag<<"\n";	
        }
      */
      
      
	}
  else if (current_gate->tag == TAG_INPUT_B) // this is an Eval input
    {
      __m128i a[2];
      
      Bytes tmp = cct.m_prng.rand_bits(Env::k());
      tmp.resize(16, 0);
      current_zero_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
      uint32_t evl_inp_ix = current_gate->wire1;
      
      tmp = (*cct.m_ot_keys)[2*evl_inp_ix+0];
      tmp.resize(16, 0);
      a[0] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
      
      tmp = (*cct.m_ot_keys)[2*evl_inp_ix+1];
      tmp.resize(16, 0);
      a[1] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
      
      // a[0] ^= zero_key; a[1] ^= zero_key ^ R;
      a[0] = _mm_xor_si128(a[0], current_zero_key);
      a[1] = _mm_xor_si128(a[1], _mm_xor_si128(current_zero_key, cct.m_R));
      
      // cct.m_out_bufr += a[0];
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), a[0]);
      cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
      
      // cct.m_out_bufr += a[1];
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), a[1]);
      cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
      
      cct.m_evl_inp_ix++; // after PCF compiler, this isn't really necessary
    }
  else
    {
#ifdef FREE_XOR
      if (current_gate->truth_table == 0x06) // if XOR gate
        {
          current_zero_key = _mm_xor_si128
            (
             *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire1)),
             *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire2))
             );
        }
      else
#endif
        {
          uint8_t bit;
          __m128i aes_key[2], aes_plaintext, aes_ciphertext;
          __m128i X[2], Y[2], Z[2];
          static Bytes tmp(16, 0);
          
          aes_plaintext = _mm_set1_epi64x(cct.m_gate_ix);
          
          X[0] = *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire1));
          Y[0] = *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire2));
          
          X[1] = _mm_xor_si128(X[0], cct.m_R); // X[1] = X[0] ^ R
          Y[1] = _mm_xor_si128(Y[0], cct.m_R); // Y[1] = Y[0] ^ R
          
          const uint8_t perm_x = _mm_extract_epi8(X[0], 0) & 0x01; // permutation bit for X
          const uint8_t perm_y = _mm_extract_epi8(Y[0], 0) & 0x01; // permutation bit for Y
          const uint8_t de_garbled_ix = (perm_y<<1)|perm_x; // wire1+2*wire2
          
          // encrypt the 0-th entry : (X[x], Y[y])
          aes_key[0] = _mm_load_si128(X+perm_x);
          aes_key[1] = _mm_load_si128(Y+perm_y);
          
          KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
          aes_ciphertext = _mm_and_si128(aes_ciphertext, cct.m_clear_mask); // clear extra bits so that only k bits left
          //bit = current_gate.m_table[de_garbled_ix];
          bit = (current_gate->truth_table>>(3-de_garbled_ix))&0x01;
          
#ifdef GRR
          // GRR technique: using zero entry's key as one of the output keys
          _mm_store_si128(Z+bit, aes_ciphertext);
          Z[1-bit] = _mm_xor_si128(Z[bit], cct.m_R);
          current_zero_key = _mm_load_si128(Z);
#else
          tmp = cct.m_prng.rand_bits(Env::k());
          tmp.resize(16, 0);
          Z[0] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
          Z[1] = _mm_xor_si128(Z[0], cct.m_R);
          
          aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
          _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);
          cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
#endif
          
          // encrypt the 1st entry : (X[1-x], Y[y])
          aes_key[0] = _mm_xor_si128(aes_key[0], cct.m_R);
          
          KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
          aes_ciphertext = _mm_and_si128(aes_ciphertext, cct.m_clear_mask);
          //bit = current_gate.m_table[0x01^de_garbled_ix];
          bit = (current_gate->truth_table>>(3-(0x01^de_garbled_ix)))&0x01;
          aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
          _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);
          cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
          
          // encrypt the 2nd entry : (X[x], Y[1-y])
          aes_key[0] = _mm_xor_si128(aes_key[0], cct.m_R);
          aes_key[1] = _mm_xor_si128(aes_key[1], cct.m_R);
          
          KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
          aes_ciphertext = _mm_and_si128(aes_ciphertext, cct.m_clear_mask);
          //bit = current_gate.m_table[0x02^de_garbled_ix];
          bit = (current_gate->truth_table>>(3-(0x02^de_garbled_ix)))&0x01;
          aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
          _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);
          cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
          
          // encrypt the 3rd entry : (X[1-x], Y[1-y])
          aes_key[0] = _mm_xor_si128(aes_key[0], cct.m_R);
          
          KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
          aes_ciphertext = _mm_and_si128(aes_ciphertext, cct.m_clear_mask);
          //bit = current_gate.m_table[0x03^de_garbled_ix];
          bit = (current_gate->truth_table>>(3-(0x03^de_garbled_ix)))&0x01;
          aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
          _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);
          cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
          
        }
      
      if (current_gate->tag == TAG_OUTPUT_A) // Gen output
        {
          cct.m_out_bufr.push_back(_mm_extract_epi8(current_zero_key, 0) & 0x01); // permutation bit
          cct.m_gen_out_ix++;
        }
      else if (current_gate->tag == TAG_OUTPUT_B) // Eval output
        {
          cct.m_out_bufr.push_back(_mm_extract_epi8(current_zero_key, 0) & 0x01); // permutation bit
          cct.m_evl_out_ix++;
        }
    }
  
  cct.m_gate_ix++;
  return &current_zero_key;
}

void * evl_next_malicious_gate(PCFState *st, PCFGate *current_gate){
  //garbled_circuit_m_t &cct = *reinterpret_cast<garbled_circuit_m_t*>(get_external_circuit(st));
  GarbledMal &cct = *reinterpret_cast<GarbledMal*>(get_external_circuit(st));
  
  static __m128i current_key;
  __m128i a;
  static Bytes tmp;
  
  // this function has no way of knowing or telling which of Gen's inputs we are pointing to, and is obviously a source of issues when we do not traverse them strictly in order
  if (current_gate->tag == TAG_INPUT_A) // Gen Input
    {
      //std::cout <<cct.m_gen_inp_mask.size()*8<<" " <<cct.m_gen_inp_ix <<" \n";
      
      //cct.tag_a_cnt++;
      //fprintf(stderr, "INPUT ALICE: %u\n",cct.tag_a_cnt);
      
      if(!(cct.m_gen_inp_ix < 8LL*cct.m_gen_inp_mask.size())){
        fprintf(stderr, "gen input index not < 8 * gen input mask size. %u !< %lu \n",cct.m_gen_inp_ix, 8*cct.m_gen_inp_mask.size());
        // fprintf(stderr,"gen input size is: %lu\n",cct.m_gen_inp.size());
      }
      
      uint8_t bit = cct.m_gen_inp_mask.get_ith_bit(cct.m_gen_inp_ix);
      Bytes::const_iterator it = cct.m_in_bufr_ix + bit*Env::key_size_in_bytes();
      
      uint32_t gen_inp_ix = current_gate->wire1;
      
      assert(cct.m_gen_inp_com.size() == cct.m_gen_inp_ix);
      
      //cct.m_gen_inp_com[gen_inp_ix] = Bytes(it, it+Env::key_size_in_bytes());
      cct.m_gen_inp_com.push_back(Bytes(it, it+Env::key_size_in_bytes()));
      
      it = cct.m_gen_inp_decom[cct.m_gen_inp_ix].begin();
      tmp.assign(it, it+Env::key_size_in_bytes());
      tmp.resize(16, 0);
      current_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
      
      /*    
      if(Env::is_root())
        {
        std::cout <<"EVLbuffr: "<<cct.m_in_bufr.to_hex()<<"\n";
        std::cout <<"EVL: gate: "<<cct.m_gate_ix<<" : "<< setBytes(current_key).to_hex() <<"\n";	
        }
      */
      cct.m_gen_inp_ix++;
    }
  else if (current_gate->tag == TAG_INPUT_B) // Eval input
    {
      uint32_t evl_inp_ix = current_gate->wire1;
      
      if(!(evl_inp_ix < 8LL *cct.m_evl_inp.size())){
        fprintf(stderr, "eval input index not < 8* circuit eval input size. %u !< %lu \n",evl_inp_ix, 8* cct.m_evl_inp.size());
      }
      assert(evl_inp_ix < 8LL * cct.m_evl_inp.size());
      uint8_t bit = cct.m_evl_inp.get_ith_bit(evl_inp_ix);
      Bytes::const_iterator it = cct.m_in_bufr_ix + bit*Env::key_size_in_bytes();
      
      tmp = (*cct.m_ot_keys)[evl_inp_ix];
      tmp.resize(16, 0);
      current_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
      
      tmp.assign(it, it+Env::key_size_in_bytes());
      tmp.resize(16, 0);
      a = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
      
      current_key = _mm_xor_si128(current_key, a);
      
      cct.m_in_bufr_ix += Env::key_size_in_bytes()*2;
      cct.m_evl_inp_ix++;
    }
  else
    {
#ifdef FREE_XOR
      if (current_gate->truth_table == 0x06)
        {
          current_key = _mm_xor_si128
            (
             *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire1)),
             *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire2))
             );
        }
      else
#endif
        {
          __m128i aes_key[2], aes_plaintext, aes_ciphertext;
          
          
          aes_plaintext = _mm_set1_epi64x(cct.m_gate_ix);
          
          aes_key[0] = *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire1));
          aes_key[1] = *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire2));
          
          const uint8_t perm_x = _mm_extract_epi8(aes_key[0], 0) & 0x01;
          const uint8_t perm_y = _mm_extract_epi8(aes_key[1], 0) & 0x01;
          
          KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
          aes_ciphertext = _mm_and_si128(aes_ciphertext, cct.m_clear_mask);
          uint8_t garbled_ix = (perm_y<<1)|perm_x;
          
#ifdef GRR
          if (garbled_ix == 0)
            {
              current_key = _mm_load_si128(&aes_ciphertext);
            }
          else
            {
              Bytes::const_iterator it = cct.m_in_bufr_ix+(garbled_ix-1)*Env::key_size_in_bytes();
              tmp.assign(it, it+Env::key_size_in_bytes());
              tmp.resize(16, 0);
              a = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
              current_key = _mm_xor_si128(aes_ciphertext, a);
            }
          cct.m_in_bufr_ix += 3*Env::key_size_in_bytes();
#else
          it = cct.m_in_bufr_ix + garbled_ix*Env::key_size_in_bytes();
          tmp.assign(it, it+Env::key_size_in_bytes());
          tmp.resize(16, 0);
          current_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
          current_key = _mm_xor_si128(current_key, aes_ciphertext);
          
          cct.m_in_bufr_ix += 4*Env::key_size_in_bytes();
#endif
        }
      /*
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire1)));
      std::cout << " " << current_gate->wire1 << " " << Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes()).to_hex();
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire2)));
      std::cout << " " << current_gate->wire2 << " " << Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes()).to_hex() << "\t";
      */

      if (current_gate->tag == TAG_OUTPUT_A) // Gen output
        {
          
          // fprintf(stderr,"Gen Output!\n");
          if (cct.m_gen_out.size()*8 <= cct.m_gen_out_ix)
            {
              // dynamically grown output array
              cct.m_gen_out.resize((cct.m_gen_out.size()+1)*2, 0);
            }
          
          uint8_t out_bit = _mm_extract_epi8(current_key, 0) & 0x01;
          out_bit ^= *cct.m_in_bufr_ix;
          cct.m_gen_out.set_ith_bit(cct.m_gen_out_ix, out_bit);
          cct.m_in_bufr_ix++;
          cct.m_gen_out_ix++;
        }
      else if (current_gate->tag == TAG_OUTPUT_B) // Eval output
        {
          // fprintf(stderr,"eval output!\n");
          if (cct.m_evl_out.size()*8 <= cct.m_evl_out_ix)
            {
              // dynamically grown output array
              cct.m_evl_out.resize((cct.m_evl_out.size()+1)*2, 0);
            }
          
          uint8_t out_bit = _mm_extract_epi8(current_key, 0) & 0x01;
          out_bit ^= *cct.m_in_bufr_ix;
          cct.m_evl_out.set_ith_bit(cct.m_evl_out_ix, out_bit);
          cct.m_in_bufr_ix++;
          
          cct.m_evl_out_ix++;
        }
    }
  //update_hash(cct, cct.m_in_bufr);
  cct.m_gate_ix++;
  
  return &current_key;
}

#endif
