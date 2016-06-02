

#include "garbling.h"



/**
   GARBLING ACCESSORY FUNCTIONS
 */



void print128_num(__m128i var)
{
  uint16_t *val = (uint16_t*) &var;
  fprintf(stderr,"Numerical: %x %x %x %x %x %x %x %x \n", 
          val[0], val[1], val[2], val[3], val[4], val[5], 
          val[6], val[7]);
}

/**
   This function doubles in GF2 by shifting "key" left one bit.
   the clear mask is included to ensure that all keys remain the same length
   (we don't want a key boundary overflowing!)
 */
void Double(__m128i & key, const __m128i & clear_mask){
  
  uint64_t * modifier = (uint64_t*) &key;
  uint64_t carry = modifier[0] & ((uint64_t) 0x8000000000000000) > 0 ? 1 : 0;
  modifier[0] = modifier[0]<<1;
  modifier[1] = (modifier[1]<<1) | carry;
  key = _mm_and_si128(key, clear_mask);

}


/**
   This function computes the permutation on an input key
   it is destructive of key so must make copies of the inputs first
   it returns H(K) = pi(L) xor L where L = 2key ^ tweak
 */
void H_Pi(__m128i & destination, __m128i &key, const __m128i & tweak,const __m128i & clear_mask, const AES_KEY_J & fixed_key){
  __m128i K; // ,K1;

  Double(key, clear_mask);
  K = _mm_xor_si128(key, tweak);

  KDF128((uint8_t*)&destination,(uint8_t*)&K, &fixed_key);
  destination = _mm_xor_si128(destination, K);
  destination = _mm_and_si128(destination,clear_mask);  
}


/**
   remember to copy the keys before they enter this function, because it's destructive to key1 and key2
 */
void H_Pi256(__m128i & destination, __m128i &key1, __m128i &key2,const __m128i & tweak, const __m128i & clear_mask, const AES_KEY_J & fixed_key){
  // takes two keys and computes a ciphertest, A la JustGarble
  __m128i K;
  
  Double(key1,clear_mask);
  Double(key2,clear_mask);
  Double(key2,clear_mask);
  
  K = _mm_xor_si128(key1,key2);
  K = _mm_xor_si128(K,tweak);

  KDF128((uint8_t*)&destination,(uint8_t*)&K, &fixed_key);
  destination = _mm_xor_si128(destination, K);
  destination = _mm_and_si128(destination,clear_mask);

}


void genHalfGatePair(__m128i& out_key, __m128i & key1, __m128i & key2, Bytes & out_bufr, const byte a1, const byte a2, const byte a3, const size_t keysize, const __m128i & clear_mask, const AES_KEY_J & fixed_key, const __m128i & R, const uint32_t j1, const uint32_t j2){
  // this function implements half-gate generation by
  // Zahur, Rosulek, and Evans

  assert(a1==0||a1==1);
  assert(a2==0||a2==1);
  assert(a3==0||a3==1);
  
  //double start = MPI_Wtime();

  const uint8_t perm_x = _mm_extract_epi8(key1,0) & 0x01;
  const uint8_t perm_y = _mm_extract_epi8(key2,0) & 0x01;

  __m128i Wa0, Wb0,Wa1,Wb1; // incoming wire keys
  __m128i Wg, We; // half gate wire keys

  // Tg and Te are the transmitted variables for the Generator and Evaluator, respectively 
  __m128i Tg, Te;
  // Ha and Hb will contain the results of applying the KDF
  __m128i Ha0, Ha1, Hb0, Hb1;

  // __m128i tmp; // temporarily stores the output of H that will be used again
  //Bytes tmp_bufr; // transfers keys to output buffer
  
  //uint32_t j1, j2;
  //j1 = increment_index();
  //j2 = increment_index();

  __m128i j1_128, j2_128;
  j1_128 = _mm_set1_epi64x(j1);
  j2_128 = _mm_set1_epi64x(j2);
  
  Wa0 = key1;
  Wb0 = key2;

  __m128i key1x = _mm_xor_si128(key1,R);

  Wa1 = _mm_xor_si128(Wa0,R);
  Wb1 = _mm_xor_si128(Wb0,R);

  __m128i H_Wa0j1, H_Wa1j1, H_Wb0j2, H_Wb1j2;
  H_Pi(H_Wa0j1, Wa0, j1_128, clear_mask, fixed_key);
  H_Pi(H_Wa1j1, Wa1, j1_128, clear_mask, fixed_key);
  H_Pi(H_Wb0j2, Wb0, j2_128, clear_mask, fixed_key);
  H_Pi(H_Wb1j2, Wb1, j2_128, clear_mask, fixed_key);
  

  // first half gate
  Tg = _mm_xor_si128(H_Wa0j1, H_Wa1j1);
  if(perm_y != a2){
    Tg = _mm_xor_si128(Tg,R);
  }

  Wg = perm_x? H_Wa1j1 : H_Wa0j1;
  if(((perm_x != a1) && (perm_y != a2)) != a3){
    Wg = _mm_xor_si128(Wg, R);
  }
  
  // second half gate
  We = perm_y ? H_Wb1j2 : H_Wb0j2;
  Te = _mm_xor_si128(H_Wb0j2,H_Wb1j2);
  Te = _mm_xor_si128(Te, a1? key1x :key1);

  // add Tg,Te to output buffer

  //  std::fill(out_bufr.begin(),out_bufr.end(),0);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&out_bufr[0]),Tg);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&out_bufr[keysize]),Te);

  // combine half gates:
  out_key = _mm_xor_si128(Wg,We);

  //benchmark_time += MPI_Wtime() - start;

}



void evlHalfGatePair(__m128i &current_key, __m128i & key1, __m128i & key2, Bytes & in_bufr,const size_t keysize,const __m128i & clear_mask,const AES_KEY_J & fixed_key,const uint32_t j1,const uint32_t j2){
  //assert(in_bufr.size() == 2*keysize);

  //double start = MPI_Wtime();

  // get the select bits
  byte sa,sb;
  sa = _mm_extract_epi8(key1,0) & 0x01;
  sb = _mm_extract_epi8(key2,0) & 0x01;

  // get the counter values
  //  uint32_t j1, j2;
  //j1 = increment_index();
  //j2 = increment_index();

  __m128i j1_128, j2_128;
  j1_128 = _mm_set1_epi64x(j1);
  j2_128 = _mm_set1_epi64x(j2);

  // fprintf(stdout,"Half Gate In Buffer: %s\n",in_bufr.to_hex().c_str());

  Bytes tmp_bufr;
  __m128i Tg, Te;
  

  // TG is always sent first, and then TE
  // where TG is the single row transmitted for the Generator's half-gate
  // and TE is the single row transmitted for the Evaluator's half-gate
  Tg = _mm_loadu_si128(reinterpret_cast<__m128i*>(&in_bufr[0]));
  Te = _mm_loadu_si128(reinterpret_cast<__m128i*>(&in_bufr[0]+keysize));
  Tg = _mm_and_si128(Tg,clear_mask);
  Te = _mm_and_si128(Te,clear_mask);


  __m128i Wa,Wb;

  Wa = key1;
  Wb = key2;

  __m128i H_Waj1,H_Wbj2; // output of hashes
  H_Pi(H_Waj1, Wa, j1_128, clear_mask, fixed_key);
  H_Pi(H_Wbj2, Wb, j2_128, clear_mask, fixed_key);


  __m128i Wg, We,tmp,tmpwe;
  __m128i xorTeKey1;

  Wg = H_Waj1;
  if(sa){
    Wg = _mm_xor_si128(H_Waj1,Tg);
  }
  
  We = H_Wbj2;
  xorTeKey1 = _mm_xor_si128(Te, key1);
  if(sb){
    We = _mm_xor_si128(H_Wbj2,xorTeKey1);
  }

  current_key = _mm_xor_si128(Wg,We);

}  


void genStandardGate(__m128i& current_key, __m128i & key1, __m128i & key2, Bytes & out_bufr, uint8_t truth_table,const __m128i & R,const size_t keysize,const __m128i & clear_mask,const AES_KEY_J & fixed_key ,const uint32_t j1){

  // X and Y are input, Z is output
  __m128i X[2], Y[2], Z[2];
  __m128i garble_ciphertext;
  __m128i garble_ciphertext1, garble_ciphertext2, garble_ciphertext3, garble_ciphertext4;
  
  uint8_t semantic_bit;
  uint8_t semantic_bit1, semantic_bit2, semantic_bit3, semantic_bit4;
    
  //double start = MPI_Wtime();

  // load the (zero-key) inputs from the PCF state container
  X[0] = key1;
  Y[0] = key2;
  // and XOR-complements
  X[1] = _mm_xor_si128(X[0], R); // X[1] = X[0] ^ R
  Y[1] = _mm_xor_si128(Y[0], R); // Y[1] = Y[0] ^ R
  

  // and get the permutation bits (tells if each zero-key has a 1 on the end)
  const uint8_t perm_x = _mm_extract_epi8(X[0],0) & 0x01;
  const uint8_t perm_y = _mm_extract_epi8(Y[0],0) & 0x01;
  const uint8_t de_garbled_ix = (perm_y << 1)|perm_x; // wire1 + 2*wire2
  
  
  // the last information for garbling
  // j1 passed to function
  //uint32_t j1;
  //j1 = increment_index();
  
  __m128i tweak;
  tweak = _mm_set1_epi64x(j1);


  
  // now run the key derivation function using the keys and the gate index  
  __m128i key1_1 = _mm_loadu_si128(X+perm_x);
  __m128i key2_1 = _mm_loadu_si128(Y+perm_y);

  __m128i key1_2 = _mm_loadu_si128(X+1-perm_x);;
  __m128i key2_2 = _mm_loadu_si128(Y+perm_y);

  __m128i key1_3 = _mm_loadu_si128(X+perm_x);
  __m128i key2_3 = _mm_loadu_si128(Y+1-perm_y);
  
  __m128i key1_4 = _mm_loadu_si128(X+1-perm_x);
  __m128i key2_4 = _mm_loadu_si128(Y+1-perm_y);
  

  H_Pi256(garble_ciphertext1, key1_1, key2_1, tweak, clear_mask, fixed_key);
  H_Pi256(garble_ciphertext2, key1_2, key2_2, tweak, clear_mask, fixed_key);
  H_Pi256(garble_ciphertext3, key1_3, key2_3, tweak, clear_mask, fixed_key);
  H_Pi256(garble_ciphertext4, key1_4, key2_4, tweak, clear_mask, fixed_key);


  semantic_bit1 = (truth_table >> (3-de_garbled_ix)) & 0x01;
  semantic_bit2 = (truth_table>>(3-(0x01^de_garbled_ix)))&0x01;
  semantic_bit3 = (truth_table>>(3-(0x02^de_garbled_ix)))&0x01;
  semantic_bit4 = (truth_table>>(3-(0x03^de_garbled_ix)))&0x01;


  _mm_store_si128(Z+semantic_bit1, garble_ciphertext1);
  Z[1 - semantic_bit1] = _mm_xor_si128(Z[semantic_bit1], R);
  current_key = _mm_loadu_si128(Z);

  
  garble_ciphertext2 = _mm_xor_si128(garble_ciphertext2, Z[semantic_bit2]);
  garble_ciphertext3 = _mm_xor_si128(garble_ciphertext3, Z[semantic_bit3]);
  garble_ciphertext4 = _mm_xor_si128(garble_ciphertext4, Z[semantic_bit4]);

 
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&out_bufr[0]),garble_ciphertext2);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&out_bufr[keysize]),garble_ciphertext3);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&out_bufr[2*keysize]),garble_ciphertext4);

 
  
  // current_key holds our output key, and it will be available to our calling function
  // the calling function will also be able to send the information in out_bufr

  //benchmark_time +=  MPI_Wtime() - start;

}


void evlStandardGate(__m128i& current_key, __m128i & key1, __m128i & key2, Bytes & in_bufr,const size_t keysize,const __m128i & clear_mask,const AES_KEY_J & fixed_key,const uint32_t j1){
  __m128i garble_key[2], aes_plaintext, garble_ciphertext;
  Bytes tmp;
  __m128i a;
  Bytes::const_iterator it;
  
  //double start = MPI_Wtime();

  //aes_plaintext = _mm_set1_epi64x(m_gate_index);
  //aes_plaintext = _mm_set1_epi64x(gate_index);


  garble_key[0] =  key1;
  garble_key[1] =  key2;
  
  const uint8_t perm_x = _mm_extract_epi8(garble_key[0], 0) & 0x01;
  const uint8_t perm_y = _mm_extract_epi8(garble_key[1], 0) & 0x01;
  
  
#ifndef AESNI
  // don't need this part
  //KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&garble_ciphertext, (uint8_t*)garble_key);     
  // garble_ciphertext = _mm_and_si128(garble_ciphertext, clear_mask);
#else
  // j1 provided as argument
  //uint32_t j1;
  //j1 = increment_index();
  
  __m128i tweak;
  tweak = _mm_set1_epi64x(j1);
  
  //  __m128i key1_in = garble_key[0];
  // __m128i key2_in = garble_key[1];
  //H_Pi256(garble_ciphertext, key1_in, key2_in, tweak, m_clear_mask, m_fixed_key);
  H_Pi256(garble_ciphertext, garble_key[0], garble_key[1], tweak, clear_mask, fixed_key);

#endif
  
  uint8_t garbled_ix = (perm_y<<1)|perm_x;
  
#ifdef GRR
  if (garbled_ix == 0) {
   current_key = _mm_load_si128(&garble_ciphertext);
  }
  else
    {
      it = in_bufr.begin() + (garbled_ix-1)*keysize;
      
      tmp.assign(it, it+keysize);
      tmp.resize(16, 0);
      a = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
      current_key = _mm_xor_si128(garble_ciphertext, a);
    }
#else
    it = in_bufr.begin() + (garbled_ix)*keysize;
    tmp.assign(it, it+keysize);
    tmp.resize(16, 0);
    current_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
    current_key = _mm_xor_si128(current_key, garble_ciphertext);
#endif
    
    // current key holds our output key, and it will be available to the calling function

    //benchmark_time +=  MPI_Wtime() - start;

}
