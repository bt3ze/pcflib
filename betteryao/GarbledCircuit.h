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
    // should also include permutation bits in this one?
    void init_Generation_Circuit(const std::vector<Bytes> * gen_keys,
                                   const std::vector<Bytes> * evl_keys);
    void init_Evaluation_Circuit(const std::vector<Bytes> * gen_keys,
                                   const std::vector<Bytes> * evl_keys);

    void evaluate_Circuit();

// functions to:
// - transform Eval's input w/ k-probe matrix?
// - evaluate Gen's input consistency?

//    void * gen_Next_Gate(PCFState *st, PCFGate *current_gate);
//  void * evl_Next_Gate(PCFState *st, PCFGate *current_gate);

    void * gen_Next_Gate(PCFGate *current_gate);
    void * evl_Next_Gate(PCFGate *current_gate);

    void set_const_key(byte c, const Bytes &key);
    const Bytes get_const_key(byte c, byte b);

protected:



    // high-level call to garble a gate 
    void garble_Gate();


    // lower-level call to compute the KDF on a couple keys
    // should also include the wire index
    void * garble_On_Keys(void * x, void * y);

    // Input Key Accessor Functions
    // get-Input intended for Eval
    // get-Key intended for Gen
    Bytes get_Gen_Input(uint32_t idx);
    Bytes get_Evl_Input(uint32_t idx);
    Bytes get_Gen_Key(uint32_t idx, uint32_t parity);
    Bytes get_Evl_Key(uint32_t idx, uint32_t parity);
    // Access to current wire table
    Bytes get_Wire_Value(uint32_t idx);


    // circuit's free-XOR value
    __m128i  m_R;
    // enforces k-length keys
    __m128i  m_clear_mask;

    // hold values of constant wires 1 and 0
    __m128i m_const_wire[2];

    // these give access to Gen and Eval's input keys
    // (useful for both evaluation and generation circuits)
    const std::vector<Bytes> * m_gen_inputs;
    const std::vector<Bytes> * m_evl_inputs;

    // remember where are are, and also this is passed to the garbling function
    // garbled gate = H_{key1} ( H_{key2} ( gate index) )
    uint64_t m_gate_index;
    
    // pointer to the PCF State
    struct PCFState *m_st;
    
    

};



#endif

