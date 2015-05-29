#ifndef GARBLED_BASE_H_
#define GARBLED_BASE_H_

#include <emmintrin.h>
#include <vector>

#include "Env.h"
#include "Prng.h"
#include "Aes.h"
#include "Bytes.h"

#ifdef __CPLUSPLUS
extern "C" {
#endif

#include "../pcflib.h"

void *copy_key(void *);
void delete_key(void *);

void *gen_next_gate(struct PCFState *st, struct PCFGate *gate);
void *evl_next_gate(struct PCFState *st, struct PCFGate *gate);

#ifdef __CPLUSPLUS
}
#endif

#define  _mm_extract_epi8(x, imm) \
	((((imm) & 0x1) == 0) ?   \
	_mm_extract_epi16((x), (imm) >> 1) & 0xff : \
	_mm_extract_epi16( _mm_srli_epi16((x), 8), (imm) >> 1))

// the maximum output for Gen and Evl is 1024 bits.
// this parameter can be changed, of course. for large circuits we must
const int MAX_OUTPUT_SIZE = 1024;

class GarbledBase {

public:
     GarbledBase();
     virtual ~GarbledBase() {};

     // virtual functions to be implemented by HBC and Malicious circuit objects
     virtual void gen_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed) = 0;
     virtual void evl_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp) = 0;

     virtual void initialize_circuit() = 0;

     __m128i             m_const_wire[2]; // keys for constant 0 and 1     
     
     struct PCFState    *m_st; // pointer to the PCF state
 

     // public methods
     inline void trim_output(){
         assert(m_gen_out.size() > 0);
         m_gen_out.resize((m_gen_out_ix+7)/8);
         m_evl_out.resize((m_evl_out_ix+7)/8);
     }

     inline void clear_and_replace_in_bufr(const Bytes &i_data){
       {
         m_in_bufr.clear();
         m_in_bufr += i_data;
         //assert(cct.m_i_bufr.size() > 0);
         m_in_bufr_ix = m_in_bufr.begin();
       }

     }

     inline const Bytes get_and_clear_out_bufr()
     {
       //static
       Bytes o_data;
       o_data.swap(m_out_bufr);
       m_out_bufr.clear();
       return o_data;
     }

     
     void set_const_key(byte c, const Bytes &key);
     const Bytes get_const_key(byte c, byte b);

     Bytes get_gen_out(){
       return m_gen_out;
     }
     Bytes get_evl_out(){
       return m_evl_out;
     }
     uint32_t get_gen_out_ix(){
       return m_gen_out_ix;
     }
     uint32_t get_evl_out_ix(){
       return m_evl_out_ix;
     }

protected:
     // these member functions are internal to the circuit object

     __m128i             m_R; // constant for free-XOR
     
     const std::vector<Bytes>  *m_ot_keys; // pointer to keys
  
     Prng                m_prng; // generator
  
     uint64_t            m_gate_ix;
     // the gate index is especially important because it is used as plaintext for AES encryption in masking wire keys
 
     __m128i             m_clear_mask; // this is used to ensure that the length of keys is correct (by setting the lowest S bits to 1 and ANDing)

     // input indices
     uint32_t            m_gen_inp_ix;
     uint32_t            m_evl_inp_ix;
     uint32_t            m_gen_out_ix;
     uint32_t            m_evl_out_ix;

     // inputs
     Bytes               m_gen_inp_mask;
     // Bytes               m_gen_inp;
     Bytes               m_evl_inp;
     Bytes               m_gen_out;
     Bytes               m_evl_out;
  
     Bytes               m_out_bufr; // out buffer
     Bytes               m_in_bufr; // in buffer
     Bytes::iterator     m_in_bufr_ix;
  
     

// gen next gate
// evl next gate
// gen next gen inp com
// evl next gen inp com

//private:
     
};


#endif

