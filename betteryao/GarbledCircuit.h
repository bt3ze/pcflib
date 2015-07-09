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

#include "macros.h"

#ifdef __CPLUSPLUS
}
#endif

// #include "macros.h"

#define  _mm_extract_epi8(x, imm) \
        ((((imm) & 0x1) == 0) ?   \
	_mm_extract_epi16((x), (imm) >> 1) & 0xff : \
	_mm_extract_epi16( _mm_srli_epi16((x), 8), (imm) >> 1))

const int MAX_OUTPUT_SIZE = 1024;

void print128_num(__m128i var);

void Double(__m128i & k, __m128i & mask);
void H_Pi(__m128i & destination, __m128i &key, __m128i & tweak, __m128i & clear_mask);
                        
class GarbledCircuit {

public:

    GarbledCircuit();
    ~GarbledCircuit() {}
    // should also include permutation bits in this one?
    void init_Generation_Circuit(const std::vector<Bytes> * gen_keys,
                                 const std::vector<Bytes> * evl_keys,
                                 Bytes & rand_seed,
                                 const Bytes & permutation_bits,
                                 const Bytes R,
                                 const Bytes & zero_key,
                                 const Bytes & one_key);
    void init_Evaluation_Circuit(const std::vector<Bytes> * gen_keys,
                                 const std::vector<Bytes> * evl_keys,
                                 const Bytes & evl_input,
                                 const Bytes & zero_key,
                                 const Bytes & one_key);
    void set_Gen_Circuit_Functions();
    void set_Evl_Circuit_Functions();

    void generate_Circuit();
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

    // pointer to the PCF State
    struct PCFState *m_st;
    // get the constant wires
    void * get_Const_Wire(uint32_t i);
 
    Bytes get_garbling_bufr();
    void set_garbling_bufr(Bytes buf);
    void clear_garbling_bufr();
    Bytes get_bob_out();
    Bytes get_alice_out();

    void trim_output_buffers();
protected:


    void set_Input_Keys(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys);

    // Gen and Eval must properly evaluate the k-probe matrix
    // as well as the hash on Gen's inputs
    // TODO: implement these 
    void evaluate_K_Probe_Matrix();
    void evalute_Gen_Inp_Hash();

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
    //Bytes get_Wire_Value(uint32_t idx);
    
    // Access to which input Gen or Eval chose
    uint32_t get_Input_Parity(uint32_t idx); 


    // things to help with garbling
    void generate_Random_Key(__m128i & destination);
 
    void generate_Alice_Input(PCFGate* current_Gate, __m128i &current_key);
    void evaluate_Alice_Input(PCFGate* current_Gate, __m128i &current_key);
    void generate_Bob_Input(PCFGate* current_Gate, __m128i &current_key);
    void evaluate_Bob_Input(PCFGate* current_Gate, __m128i &current_key);
    void generate_Gate(PCFGate* current_Gate, __m128i &current_key);
    void evaluate_Gate(PCFGate* current_Gate, __m128i &current_key);
    void evaluate_Alice_Output(PCFGate* current_Gate, __m128i &current_key);
    void evaluate_Bob_Output(PCFGate* current_Gate, __m128i &current_key);
    void generate_Alice_Output(PCFGate* current_Gate, __m128i &current_key);
    void generate_Bob_Output(PCFGate* current_Gate, __m128i &current_key);
    
    void genHalfGate(PCFGate* current_Gate, __m128i &current_key);
    void evlHalfGate(PCFGate* current_Gate, __m128i &current_key);

    void genHalfGatePair(__m128i& out_key, __m128i & key1, __m128i & key2, Bytes & out_bufr, byte a1, byte a2, byte a3); 
    void evlHalfGatePair(__m128i& current_key, __m128i & key1, __m128i & key2, Bytes & in_bufr); 
    
    void xor_Gate(PCFGate* current_gate, __m128i &current_key);


    uint32_t increment_index();

    // circuit's free-XOR value
    __m128i  m_R;
    // enforces k-length keys
    __m128i  m_clear_mask;
    
    // m_prng used for generating new wire keys
    Prng m_prng;
    // permutation bits important for generation circuits
    // (they will be secret to evaluation circuits)
    //    const Bytes * m_select_bits;
    Bytes m_select_bits;

    // hold values of constant wires 1 and 0
    __m128i m_const_wire[2];

    // these give access to Gen and Eval's input keys
    // (useful for both evaluation and generation circuits)
    const std::vector<Bytes> * m_gen_inputs;
    const std::vector<Bytes> * m_evl_inputs;

    // remember where are are, and also this is passed to the garbling function
    // garbled gate = H_{key1} ( H_{key2} ( gate index) )
    uint64_t m_gate_index;

    // this variable will be used by Gen to send messages to Eval
    // and by Eval to receive Gen's messages
    // it will contain garbling information for each gate sent
    Bytes m_garbling_bufr;
    
    // a couple of buffers that hold circuit output, used in the garbling functions
    Bytes m_alice_out;
    Bytes m_bob_out;
    uint32_t m_in_bufr_ix;
    uint32_t m_alice_out_ix;
    uint32_t m_bob_out_ix;
};



#endif

