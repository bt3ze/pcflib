#ifndef GARBLED_CIRCUIT_H
#define GARBLED_CIRCUIT_H

#include <emmintrin.h>
#include <vector>

#include "Env.h"
#include "Prng.h"
#include "Aes.h"
#include "Bytes.h"
#include "garbling.h"

//#ifdef __CPLUSPLUS
//extern "C" {
//#endif

#include "../pcflib.h"


/*
#ifdef __cplusplus
extern "C" {
#include "cthreadpool/thpool.h"
#endif
#ifdef __cplusplus
}
#endif
*/

//  void *copy_key(void *);
void copy_key(void *, void*);
void delete_key(void *);

void *gen_next_gate(struct PCFState *st, struct PCFGate *gate);
void *evl_next_gate(struct PCFState *st, struct PCFGate *gate);

#include "macros.h"

//#ifdef __CPLUSPLUS
//}
//#endif


#define MESSAGE_LIMIT 300

const int MAX_OUTPUT_SIZE = 1024;

void print128_num(__m128i var);

void Double(__m128i & k, __m128i & mask);
void H_Pi(__m128i & destination, __m128i &key, __m128i & tweak, __m128i & clear_mask);
             
void save_Key_to_128bit(const Bytes & key, __m128i & destination);
void append_m128i_to_Bytes(const __m128i & num, Bytes & buf);
void insert_to_garbling_bufr(const __m128i & num, Bytes & dest, Bytes & tmp_bufr, uint32_t pos);

          
#include <time.h>
#define BILN 1E9
static double copy_time = 0.0;
static double buffer_time = 0.0;
static int num_copies = 0;
static int num_buffers = 0;
static struct timespec copy_start, copy_end;

static struct timespec buffer_start, buffer_end;

static int num_b_cpy = 0;
static double b_cpy_time = 0.0;
static struct timespec better_cpy_start, better_cpy_end;

static int num_comms = 0;
static double comm_time = 0.0;
static struct timespec comm_start, comm_end;


class GarbledCircuit {

public:

    GarbledCircuit();
    ~GarbledCircuit() {}

    void init_Generation_Circuit(const std::vector<Bytes> * gen_keys,
                                 const std::vector<Bytes> * evl_keys,
                                 const uint32_t gen_inp_size,
                                 Bytes & rand_seed,
                                 const Bytes & permutation_bits,
                                 const Bytes R,
                                 const Bytes & zero_key,
                                 const Bytes & one_key);
    void init_Evaluation_Circuit(const std::vector<Bytes> * gen_keys,
                                 const std::vector<Bytes> * evl_keys,
                                 const uint32_t gen_inp_size,
                                 const Bytes & evl_input,
                                 const Bytes & zero_key,
                                 const Bytes & one_key);

    // need to set callbacks functions for the PCF interpreter
    void set_Gen_Circuit_Functions();
    void set_Evl_Circuit_Functions();
    void init_circuit_AES_key(__m128i & key);


    // Gen and Eval must properly evaluate the k-probe matrix
    // as well as the hash on Gen's inputs
    // TODO: implement these 
    void evaluate_K_Probe_Matrix(std::vector<Bytes> &matrix);
    void generate_K_Probe_Matrix(std::vector<Bytes> &matrix);
    void evl_next_hash_row(Bytes & row, Bytes & in_bufr);
    void gen_next_hash_row(Bytes & row, Bytes & out_bufr);


    void * gen_Next_Gate(PCFGate *current_gate);
    void * evl_Next_Gate(PCFGate *current_gate);

    void Garble_Circuit();
    void Evaluate_Circuit();


    // pointer to the PCF State
    struct PCFState *m_st;
    // get the constant wires
    void * get_Const_Wire(uint32_t i);
 
    void clear_garbling_bufr();
    Bytes get_bob_out();
    Bytes get_alice_out();
    Bytes get_hash_out(); // get output of 2UHF hash

    void trim_output_buffers();

    // for Gen output authenticity proof
    Bytes get_Gen_Output_Label(uint32_t idx);


    // timing communication
    double m_comm_time;

    void send_buffer();

protected:


    void set_Input_Keys(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys);



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
    
    // Access to which input Gen or Eval chose
    uint32_t get_Input_Parity(uint32_t idx); 




    // these are the individual functions for the types of gates
    // alice in, alice out, bob in, bob out, and general gate
 
    void generate_Alice_Input(PCFGate* current_Gate, __m128i &current_key, Bytes & garbling_bufr);
    void evaluate_Alice_Input(PCFGate* current_Gate, __m128i &current_key, Bytes & garbling_bufr);
    void generate_Bob_Input(PCFGate* current_Gate, __m128i &current_key);
    void evaluate_Bob_Input(PCFGate* current_Gate, __m128i &current_key);
    void evaluate_Alice_Output(PCFGate* current_Gate, __m128i &current_key, Bytes & garbling_bufr);
    void evaluate_Bob_Output(PCFGate* current_Gate, __m128i &current_key, Bytes & garbling_bufr);
    void generate_Alice_Output(PCFGate* current_Gate, __m128i &current_key, Bytes & garbling_bufr);
    void generate_Bob_Output(PCFGate* current_Gate, __m128i &current_key, Bytes & garbling_bufr);
    

    void generate_Gate(PCFGate* current_Gate, __m128i &current_key, Bytes & garbling_bufr);
    void evaluate_Gate(PCFGate* current_Gate, __m128i &current_key, Bytes & garbling_bufr);

 
    //void genHalfGatePair(__m128i& out_key, __m128i & key1, __m128i & key2, Bytes & out_bufr, byte a1, byte a2, byte a3); 
    
    //void genStandardGate(__m128i& current_key, __m128i & key1, __m128i & key2, Bytes & out_bufr,uint8_t truth_table);

    //void evlStandardGate(__m128i& current_key, __m128i & key1, __m128i & key2, Bytes & in_bufr);

    //void xor_Gate(__m128i & key1, __m128i & key2, __m128i &current_key);


    uint32_t increment_index();

    // circuit's free-XOR value
    __m128i  m_R;
    // convenient copy of m_R
    Bytes m_R_bytes;

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
    
    // these will be used to move things around, between buffers, and into _mm_128i fields
    //Bytes m_temp_bufr;
    //__m128i key_bufr[4];

    /**
       Parallel implementation additions
       we want to garble about MESSAGE_LIMIT gates at a time (or in a batch)
       and need a specific garbling address for each to ensure no data races hurt us
       also need a garbling buffer for each one
     */
    __m128i m_parallel_keys[MESSAGE_LIMIT];
    Bytes m_parallel_buffers[MESSAGE_LIMIT];
    uint32_t m_batch_index;

    void send_half_gate(const Bytes &buf);
    void send_full_gate(const Bytes &buf);
    void read_half_gate(Bytes &buf);
    void read_full_gate(Bytes &buf);

    // a couple of buffers that hold circuit output, used in the garbling functions
    Bytes m_alice_out;
    Bytes m_bob_out;
    uint32_t m_in_bufr_ix;
    uint32_t m_alice_out_ix;
    uint32_t m_bob_out_ix;

    // for Gen output authenticity proof
    std::vector<Bytes> m_gen_output_labels;

    // for Gen input consistency
    Bytes m_hash_out; // stores the hash function output
    uint32_t m_hash_row_idx;
    
    uint32_t m_gen_inp_size; //important for accessing output mask keys and hash keys

    // very important: used for encrypting
    AES_KEY_J m_fixed_key;


    uint32_t xor_gates;
    uint32_t half_gates;
    uint32_t other_gates;
    uint32_t total_gates;

    double xor_time;
    double hg_time; // half gate time
    double og_time; // other gate tim
    double garble_time; // other gate time

    struct timespec xor_start, xor_end;
    struct timespec half_start, half_end;
    struct timespec og_start, og_end;
    struct timespec garble_start, garble_end;


    Bytes m_message_queue;
    Bytes m_ciphertext_buff; // should always be two or three ciphertexts long

    // serial implementation!
    // follows the next few

    uint32_t m_message_limit;
    uint32_t m_messages_waiting; // provide similar functions
    uint32_t m_queue_index; // provide similar functions

    // these functions are for batch sending gates
    // they are all void because information will be returned by their calling functions
    void enqueue_messages(const Bytes & source, uint32_t num);
    void add_messages_to_queue(const Bytes & src, uint32_t num);
  

    void retrieve_buffer();
    void retrieve_ciphertexts(Bytes & buf, uint32_t num_ctexts);


    // parallel implementation:
    // the message buffer itself will have to be (3*keysize)*NUM_MESSAGES bytes
    // so to get the benefits of halfgates we'd need to actually eliminate all GRR gates (TODO)
    // the reason is that we will assign an index to every gate
    // and insert its garbled ciphertexts in a slotted position
    // in order to avoid parallel data races
    // once we've garbled the right number of gates,
    // we send off the whole buffer at once

    void insert_garbled_ciphertext(const Bytes & buf, const uint32_t idx, const uint32_t ctext_per_gate, const uint32_t keysize_bytes);
    void send_full_message_queue();
    
    void * gen_Next_Gate_Parallel(PCFGate * current_gate);

};


/*
We can send each gate separately to a thread to be evaluated

 */


typedef struct garble_gate_arg{
  PCFState * st;
  PCFGate * current_gate;
  __m128i * current_key;
  Bytes * garbling_bufr;
} garble_gate_arg;



void * garble_gate_parallel(void * arg);



#endif

