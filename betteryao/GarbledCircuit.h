#ifndef GARBLED_CIRCUIT_H
#define GARBLED_CIRCUIT_H

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

const int MAX_OUTPUT_SIZE = 1024;

class GarbledCircuit {

public:

    GarbledCircuit();
    ~GarbledCircuit() {}
    void init_Generation_Circuit();
    void init_Evaluation_Circuit();

    void * gen_Next_Gate(PCFState *st, PCFGate *current_gate);
    void * evl_Next_Gate(PCFState *st, PCFGate *current_gate);

protected:

    // high-level call to garble a gate 
    void garble_Gate();
    // lower-level call to compute the KDF on a couple keys
    void garble_On_Keys();


    // circuit's free-XOR value
    __m128i  m_R;
    // enforces k-length keys
    __m128i  m_clear_mask;

    // hold values of constant wires 1 and 0
    __m128i m_const_wire[2];

    // these give access to Gen and Eval's input keys
    // (useful for both evaluation and generation circuits)
    const std::vector<Bytes>* gen_inputs;
    const std::vector<Bytes>* evl_inputs;

    // remember where are are, and also this is passed to the garbling function
    // garbled gate = H_{key1} ( H_{key2} ( gate index) )
    uint64_t m_gate_index;
    
    // pointer to the PCF State
    struct PCFState *m_st;
    
    

};



#endif

