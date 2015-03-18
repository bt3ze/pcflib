#ifndef GARBLED_BASE_H_
#define GARBLED_BASE_H_

#include <emmintrin.h>
#include <vector>

#include "Env.h"
#include "Prng.h"
//#include "Hash.h"
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

const int CIRCUIT_HASH_BUFFER_SIZE = 1024*1024;
const int MAX_OUTPUT_SIZE = 1024;

class GarbledBase {

public:
     GarbledBase();
     virtual ~GarbledBase();

     virtual void * gen_next_gate(struct PCFState *st, struct PCFGate *current_gate) = 0; 
     virtual void * evl_next_gate(struct PCFState *st, struct PCFGate *current_gate) = 0;

protected:

     Bytes               m_bufr;
     
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
  
     Bytes               m_out_bufr; // out buffer
     Bytes               m_in_bufr; // in buffer
     Bytes::iterator     m_in_bufr_ix;
  
     struct PCFState    *m_st; // pointer to the PCF state
     __m128i             m_const_wire[2]; // keys for constant 0 and 1     

     
     // methods

     inline void trim_output();
     inline void clear_and_replace_in_bufr(const Bytes &);
     inline const Bytes get_and_clear_out_bufr();
     
     void set_const_key(byte c, const Bytes &key);
     const Bytes get_const_key(byte c, byte b);

     virtual void gen_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed);
     virtual void evl_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp);

     virtual void initialize_circuit();

// gen next gate
// evl next gate
// gen next gen inp com
// evl next gen inp com

//private:
     
};


#endif

