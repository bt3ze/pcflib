// garbled_circuit_m.cpp

#include <iostream>
#include "garbled_circuit_m.h"

const Bytes get_const_key(garbled_circuit_m_t &cct, byte c, byte b)
{
	assert(c == 0 || c == 1); // wire for constant 0 or 1
	assert(b == 0 || b == 1); // with bit value 0 or 1
	Bytes tmp(16);

	tmp.resize(16);
	if (b)
	{
		_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), _mm_xor_si128(cct.m_R, cct.m_const_wire[c]));
	}
	else
	{
		_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), cct.m_const_wire[c]);
	}

	tmp.resize(Env::key_size_in_bytes());
	return tmp;
}

void set_const_key(garbled_circuit_m_t &cct, byte c, const Bytes &key)
{
	assert(c == 0 || c == 1); // wire for constant 0 or 1
	Bytes tmp = key;
	tmp.resize(16);
	cct.m_const_wire[c] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
}

namespace
{ // already transfered to new garbledbase
const int CIRCUIT_HASH_BUFFER_SIZE = 1024*1024;
const int MAX_OUTPUT_SIZE = 1024;

  /*
void update_hash(garbled_circuit_m_t &cct, const Bytes &data)
{
	cct.m_bufr += data;

#ifdef RAND_SEED
	if (cct.m_bufr.size() > CIRCUIT_HASH_BUFFER_SIZE) // hash the circuit by chunks
	{
          //cct.m_hash.update(cct.m_bufr);
		cct.m_bufr.clear();
	}
#endif
}
  */
 
};

/**
   should be renamed initialize_generation_circuit
 */
void gen_init_circuit(garbled_circuit_m_t &cct, const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed)
{
	cct.m_ot_keys = &ot_keys;
	cct.m_gen_inp_mask = gen_inp_mask;
	cct.m_prng.seed_rand(seed);

	// R is a random k-bit string whose 0-th bit has to be 1
	Bytes tmp;
        
        // Eval does not do this part for evaluation circuits
	tmp = cct.m_prng.rand_bits(Env::k());
	tmp.set_ith_bit(0, 1);
	tmp.resize(16, 0);
	cct.m_R = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

	// pick zero-keys for constant wires
	tmp = cct.m_prng.rand_bits(Env::k());
	tmp.resize(16, 0);
	cct.m_const_wire[0] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

	tmp = cct.m_prng.rand_bits(Env::k());
	tmp.resize(16, 0);
	cct.m_const_wire[1] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));
        
	initialize_malicious_circuit(cct);
        
	cct.m_gen_inp_decom.clear();
}

/**
   initialize the circuit with ot keys, masked get input, and eval input
   clear the output buffers
   should be renamed initialize_evaluation_circuit
 */
void evl_init_circuit(garbled_circuit_m_t &cct, const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp)
{

  initialize_malicious_circuit(cct);

  cct.m_ot_keys = &ot_keys;
  cct.m_gen_inp_mask = masked_gen_inp; // this was overwritten in the last version
  cct.m_evl_inp = evl_inp;
  
  cct.m_evl_out.reserve(MAX_OUTPUT_SIZE);
  cct.m_evl_out.clear(); // will grow dynamically
  cct.m_gen_out.reserve(MAX_OUTPUT_SIZE);
  cct.m_gen_out.clear(); // will grow dynamically
  
  //cct.tag_a_cnt = 0;

  //	cct.m_bufr.reserve(CIRCUIT_HASH_BUFFER_SIZE);
  //cct.m_bufr.clear();
  //cct.m_hash.init();
  
  // clear input commitments
  cct.m_gen_inp_com.clear();
  //cct.m_gen_inp_decom.clear();
}


void *gen_next_gate_m(struct PCFState *st, struct PCFGate *current_gate)
{
	garbled_circuit_m_t &cct =
		*reinterpret_cast<garbled_circuit_m_t*>(get_external_circuit(st));

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
	
                //if(Env::is_root())
                // {
                //                   std::cout <<"GENbuffr: "<<cct.m_out_bufr.to_hex()<<"\n";		
                //                   std::cout <<"GEN: gate: "<<cct.m_gate_ix<<" : "<< setBytes(current_zero_key).to_hex()<<" "<< setBytes(onev).to_hex() <<"  "<<current_gate->tag<<"\n";	
                   //	}


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
/*
		Bytes tmp(16);
		__m128i XX;
		XX = *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire1));
	_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), XX);
std::cerr << " " << current_gate->wire1 << " (" << Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes()).to_hex();
	_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), _mm_xor_si128(XX, cct.m_R));
std::cerr << ", " << Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes()).to_hex() << ")";

		XX = *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire2));
	_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), XX);
std::cerr << " " << current_gate->wire2 << " (" << Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes()).to_hex();
	_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), _mm_xor_si128(XX, cct.m_R));
std::cerr << ", " << Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes()).to_hex() << ")\t";
*/
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
/*

Bytes tmp(16);
	_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), current_zero_key);
std::cerr << cct.m_gate_ix << ": (" << Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes()).to_hex();
	_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), _mm_xor_si128(current_zero_key, cct.m_R));
std::cerr << ", " << Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes()).to_hex() << ")" << std::endl;
}
*/
	cct.m_gate_ix++;
	return &current_zero_key;
}

void *evl_next_gate_m(struct PCFState *st, struct PCFGate *current_gate)
{
  garbled_circuit_m_t &cct = *reinterpret_cast<garbled_circuit_m_t*>(get_external_circuit(st));
  
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
      // assert(cct.m_gen_inp_ix < cct.m_gen_inp_mask.size());
      // HERE! THE BUG IS HERE!
      //fprintf(stderr,"gen input index: %u \t rank: %u\n",cct.m_gen_inp_ix, Env::group_rank());
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
//	_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire1)));
//std::cout << " " << current_gate->wire1 << " " << Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes()).to_hex();
//	_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), *reinterpret_cast<__m128i*>(get_wire_key(st, current_gate->wire2)));
//std::cout << " " << current_gate->wire2 << " " << Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes()).to_hex() << "\t";

		if (current_gate->tag == TAG_OUTPUT_A) // Gen output
		{
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
/*

	tmp.resize(16);
	_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), current_key);
std::cout << cct.m_gate_ix << ": " << Bytes(tmp.begin(), tmp.begin()+Env::key_size_in_bytes()).to_hex() << std::endl;
*/
	//update_hash(cct, cct.m_in_bufr);
	cct.m_gate_ix++;

	return &current_key;
}


/**
   @cct is the garbled circuit we're working on
   @row is a row of the 2-UHF
   @kx is the row number
 */
void gen_next_gen_inp_com(garbled_circuit_m_t &cct, const Bytes &row, size_t kx)
{
	Bytes tmp;

	__m128i out_key[2];
	tmp = cct.m_prng.rand_bits(Env::k());
	tmp.set_ith_bit(0, 0);
	tmp.resize(16, 0);
	out_key[0] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
	out_key[1] = _mm_xor_si128(out_key[0], cct.m_R);

        // this should always be true, since Gen commits to two wires for each input
	assert(cct.m_gen_inp_decom.size() % 2 == 0);

        // cct.m_gen_inp_decom.size() is twice Gen's input size for Gen,
        // exactly Gen's input size for Eval
        // this looks like it might be a bug
        // and m_gen_inp_decom[0].size() is the size of a single decommitment
	// it should be 2*key size
        // (and that should be made more explicit in the code)
        Bytes msg(cct.m_gen_inp_decom[0].size());
	for (size_t jx = 0; jx < cct.m_gen_inp_decom.size()/2; jx++)
          {
            if (row.get_ith_bit(jx))
              {
                byte bit = cct.m_gen_inp_mask.get_ith_bit(jx);
                msg ^= cct.m_gen_inp_decom[2*jx+bit];
                // xor all of the decommitments for which gen's input mask is a 1
                // i do not understand this yet
              }
          }
        
	__m128i in_key[2], aes_plaintext, aes_ciphertext;

        // this initializes both 64-bit halves of aes_plaintext to kx+10
        // why + 10??
	aes_plaintext = _mm_set1_epi64x((uint64_t)kx+10);

        // copy msg into tmp
	tmp.assign(msg.begin(), msg.begin()+Env::key_size_in_bytes());
        // and resize tmp so that it's 128-bits, using 0 to fill beyond its current size
	tmp.resize(16, 0);
        // in_key0 is tmp (which is message, which was set above)
	in_key[0] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
        // in_key1 is in_key[0] XORed with R (for free XOR trick)
	in_key[1] = _mm_xor_si128(in_key[0], cct.m_R);

        // send addresses to the key derivation function
	KDF128((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)&in_key[0]);
        // reduce to the number of bits we've selected w/ m_clear_mask
	aes_ciphertext = _mm_and_si128(aes_ciphertext, cct.m_clear_mask);
        // and mask the outkey with the ciphertext
	out_key[0] = _mm_xor_si128(out_key[0], aes_ciphertext);

        // send pointers to the key derivation function
	KDF128((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)&in_key[1]);
        // mask to the right number of bits
	aes_ciphertext = _mm_and_si128(aes_ciphertext, cct.m_clear_mask);
        // mask the out key with the ciphertext
	out_key[1] = _mm_xor_si128(out_key[1], aes_ciphertext);

        // this should be the permutation bit
	const byte bit = msg.get_ith_bit(0);

	tmp.resize(16);
	_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), out_key[  bit]);
	// insert out_key[bit] at the end
        cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());

	_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), out_key[1-bit]);
        // insert out_key[1-bit] at the very end
	cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());

	//cct.m_gen_inp_hash_ix++;
}

void evl_next_gen_inp_com(garbled_circuit_m_t &cct, const Bytes &row, size_t kx)
{
  // this should be 2*k() bits
  Bytes out(cct.m_gen_inp_decom[0].size());
  
  // this time, m_gen_inp_decom.size() ought to be m_gen_inp_cnt
  for (size_t jx = 0; jx < cct.m_gen_inp_decom.size(); jx++)
    { // evalute 2-UHF
      if (row.get_ith_bit(jx)) { out ^= cct.m_gen_inp_decom[jx]; }
    }
  
  byte bit = out.get_ith_bit(0); // permutation bit
  assert((bit == 0 || bit == 1));
  
  //static Bytes tmp;
  Bytes tmp;
  
  // strange formula for the iterator
  Bytes::iterator it = cct.m_in_bufr_ix + bit*Env::key_size_in_bytes();
  
  
  __m128i aes_key, aes_plaintext, aes_ciphertext, out_key;
  
  // copy out into tmp, then resize it to ensure 128 bits
  tmp.assign(out.begin(), out.begin()+Env::key_size_in_bytes());
  tmp.resize(16, 0);
  aes_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
  
  aes_plaintext = _mm_set1_epi64x((uint64_t)kx+10); // here again, +10. but why?
  
  KDF128((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)&aes_key);
  aes_ciphertext = _mm_and_si128(aes_ciphertext, cct.m_clear_mask);
  
  // for evl, we don't need to do two ciphertexts, just one
  tmp.assign(it, it+Env::key_size_in_bytes());
  tmp.resize(16, 0);
  out_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
  // mask the out key with aes_ciphertext
  out_key = _mm_xor_si128(out_key, aes_ciphertext);
  
  // and set the input hash bit?
  bit = _mm_extract_epi8(out_key, 0) & 0x01;
  cct.m_gen_inp_hash.set_ith_bit(kx, bit);
  
  // then increment the input buffer by 2?
  cct.m_in_bufr_ix += 2*Env::key_size_in_bytes();
        
}

