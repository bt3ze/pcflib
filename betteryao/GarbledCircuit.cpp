#ifndef GARBLED_CIRCUIT_CPP
#define GARBLED_CIRCUIT_CPP

#include "GarbledCircuit.h"
#include <stdio.h>
#include <iostream>
// notes: still need functions to load the circuit file and/or
// set the loaded circuit file to a field in the garbled circuit for access


void *copy_key(void *old_key)
{
  __m128i *new_key = 0;
  if (old_key != 0)
    {
      // first argument is size, second argument is allignment
      new_key = (__m128i*)_mm_malloc(sizeof(__m128i), sizeof(__m128i));
      *new_key = *reinterpret_cast<__m128i*>(old_key);
    }
  return new_key;
}

void delete_key(void *key)
{
  if (key != 0) _mm_free(key);
}

 
GarbledCircuit::GarbledCircuit(): m_gate_index(0) {

  // initialize the key Mask
  Bytes tmp(16);
  for(size_t ix = 0; ix< Env::k(); ix++){
    tmp.set_ith_bit(ix,1);
  }
  m_clear_mask = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));

}

// the big cahoona
// this is called by the BetterYao protocol object
// and runs the circuit till the end.
void GarbledCircuit::generate_Circuit(){
  // gen and eval versions
  set_key_copy_function(m_st, copy_key);
  set_key_delete_function(m_st, delete_key);
  set_callback(m_st,gen_next_gate);
  
  
  while(get_next_gate(m_st)){
    // do things
    
  }

}

// the big cahoona
// this is called by the BetterYao protocol object
// and runs the circuit till the end.
void GarbledCircuit::evaluate_Circuit(){
  // gen and eval versions
  set_key_copy_function(m_st, copy_key);
  set_key_delete_function(m_st, delete_key);

  set_callback(m_st,evl_next_gate);

  
  while(get_next_gate(m_st)){
    // do things
    
  }

}

/**
   These functions serve as intermediaries
   to get back to the Garbled Circuit object
   since we go through such crazy calling loops with 
   the pcf interpreter struct
 */
void * gen_next_gate(PCFState *st, PCFGate *current_gate){

  fprintf(stderr,"gen next gate\n");

  // returns void pointer, which is pointer to a key?
  // use this one to call the Garbled Circuit object again
  GarbledCircuit &cct = *reinterpret_cast<GarbledCircuit*>(get_external_circuit(st));

  fprintf(stderr,"gen next gate\n");
  // now, call the appropriate function from cct
  return cct.gen_Next_Gate(current_gate);
 }

void * evl_next_gate(PCFState *st, PCFGate *current_gate){
  // returns void pointer, which is pointer to a key?
  // use this one to call the Garbled Circuit object again
  GarbledCircuit &cct = *reinterpret_cast<GarbledCircuit*>(get_external_circuit(st));
  
  // now, call the appropriate function from cct
  return cct.evl_Next_Gate(current_gate);
 }


void GarbledCircuit::init_Generation_Circuit(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys, Bytes & circuit_seed, const Bytes * permutation_bits){
  
  Bytes tmp;
  
  m_prng.seed_rand(circuit_seed);
  m_permutation_bits = permutation_bits;

  // initialize the value of R
  tmp = m_prng.rand_bits(Env::k());
  tmp.resize(16,0);
  tmp.set_ith_bit(0,1);
  m_R = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  // initialize the constant wires
  tmp = m_prng.rand_bits(Env::k());
  tmp.resize(16,0);
  m_const_wire[0] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  tmp = m_prng.rand_bits(Env::k());
  tmp.resize(16,0);
  m_const_wire[1] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  // set the input keys
  set_Input_Keys(gen_keys, evl_keys);

}


void GarbledCircuit::set_Input_Keys(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys){
  m_gen_inputs = gen_keys;
  m_evl_inputs = evl_keys;
}

void GarbledCircuit::init_Evaluation_Circuit(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys){
  set_Input_Keys(gen_keys, evl_keys);

  //std::cout << "Gen first input key: " << (*gen_keys)[0].to_hex() << std::endl;
}


void * GarbledCircuit::get_Const_Wire(uint32_t i){
  assert(i == 0 || i == 1);
  // not sure i wanna reference it like this,
  // but seems OK for now.
  return m_const_wire+i;
}

void GarbledCircuit::garble_Gate(){


}


void * GarbledCircuit::garble_On_Keys(void * x_key, void * y_key){
  __m128i X[2], Y[2], Z[2];
  // __m128i aes_key[2], aes_plaintext, aes_ciphertext;
  Bytes tmp(16,0);
  
  

  
}

/**
   Input Key Accessor Functions 
   get-Key functions are intended for Gen (although they call the get-Input) functions
   get-Input functions are intended for Eval
 */

Bytes GarbledCircuit::get_Gen_Key(uint32_t input_idx, uint32_t parity){
  assert(parity == 0 || parity == 1);
  fprintf(stderr,"Get Gen Key: %u %u\n",input_idx,parity);
  return get_Gen_Input(2*input_idx + parity);
}

Bytes GarbledCircuit::get_Evl_Key(uint32_t input_idx, uint32_t parity){
  assert(parity == 0 || parity == 1);
  return get_Evl_Input(2*input_idx + parity);
}

Bytes GarbledCircuit::get_Gen_Input(uint32_t idx){
  fprintf(stderr,"Get Gen Input idx: %u\n",idx);
  Bytes ret = (*m_gen_inputs)[idx];
  fprintf(stderr, "Gen's Input Key: %s\n", ret.to_hex().c_str());
  return ret;
}

Bytes GarbledCircuit::get_Evl_Input(uint32_t idx){
  return (*m_evl_inputs)[idx];
}

/**
   Must be able to access the current wire values at indices from the wire table
 */

//Bytes GarbledCircuit::get_Wire_Value(uint32_t idx){}


void save_Key_to_128bit(Bytes & key, __m128i & destination){
  Bytes tmp = key;
  tmp.resize(16,0);
  destination = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
}

uint32_t GarbledCircuit::get_Input_Parity(uint32_t idx){
  fprintf(stderr,"get %u input parity: %u\n",idx,m_permutation_bits->get_ith_bit(idx));
  return m_permutation_bits->get_ith_bit(idx);
}

void GarbledCircuit::generate_Random_Key(__m128i & destination){
  Bytes tmp = m_prng.rand_bits(Env::k());
  tmp.resize(16,0);
  destination = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
}

void print128_num(__m128i var)
{
    uint16_t *val = (uint16_t*) &var;
    fprintf(stderr,"Numerical: %x %x %x %x %x %x %x %x \n", 
           val[0], val[1], val[2], val[3], val[4], val[5], 
           val[6], val[7]);
}

void * GarbledCircuit::gen_Next_Gate(PCFGate *current_gate){
  // somehow the PCFState should be available through the m_st pointer
  
  // first, we generate a new key that will be made available
  // to the pcf state at the end.
  // we manipulate the current_key as we generate the gate.

  static __m128i current_key; // must be static to return it
  
  if(current_gate->tag == TAG_INPUT_A){

    // Gen's input keys have already been generated and determined
    // by his input keys
    // here, we simply return the proper input key
    fprintf(stderr,"Alice/Gen Input Gate\n");

    uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index
    
    fprintf(stderr,"Next Wire Index: %u\n", gen_input_idx);
    
    //    Bytes gen_input = get_Gen_Key(gen_input_idx, get_Input_Parity(gen_input_idx));
    Bytes gen_input = get_Gen_Input(2*gen_input_idx + get_Input_Parity(gen_input_idx));

    fprintf(stderr,"Save to 128 Bit\n");
    save_Key_to_128bit(gen_input,current_key);
    
    return &current_key; // get_Gen_Key or create a random nonce


  }  else if (current_gate->tag == TAG_INPUT_B){

    // Here, Eval already has her input keys, but Gen doesn't know
    // which to use. So he will generate a new key, encrypt it
    // using both of Eval's inputs, and send the key to her.
    // She will decrypt the proper one and use it 
    Bytes tmp;
    generate_Random_Key(current_key);
    
    fprintf(stderr,"Eval Key: ");
    print128_num(current_key);
    fprintf(stderr,"\n");
    fprintf(stderr,"Bob/Evl Input Gate\n");

    return &current_key;

  } else if (current_gate->tag == TAG_OUTPUT_A) {

  } else if (current_gate->tag == TAG_OUTPUT_B){
    
  } else {
    // actual gate
  }
  
}

void * GarbledCircuit::evl_Next_Gate(PCFGate *current_gate){
  
 static __m128i current_key; // must be static to return it
 
  if(current_gate->tag == TAG_INPUT_A){
    // here, Gen's input keys have been sent to Eval
    // she needs only to assign the input key to the wire
    /*
    uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index
    Bytes gen_input = get_Gen_Input(gen_input_idx);
    save_Key_to_128bit(gen_input,current_key);
    
    return &current_key; // get_Gen_Key or create a random nonce
    */
    fprintf(stderr,"Alice/Gen Input Gate\n");

    uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index
    
    fprintf(stderr,"Next Wire Index: %u\n", gen_input_idx);

    Bytes gen_input = get_Gen_Input(gen_input_idx);
    fprintf(stderr,"Save to 128 Bit\n");
    save_Key_to_128bit(gen_input,current_key);
    
    return &current_key; // get_Gen_Key or create a random nonce

    
  } else if (current_gate->tag == TAG_INPUT_B){
    // here, Eval already knows which input key she wants to use
    // she selects it and assigns it to her wire value

    uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index
    Bytes evl_input = get_Evl_Input(gen_input_idx);
    save_Key_to_128bit(evl_input,current_key);
    
    return &current_key; // get_Gen_Key or create a random nonce
    
  } else if (current_gate->tag == TAG_OUTPUT_A) {

  } else if (current_gate->tag == TAG_OUTPUT_B){
    
  } else {
    // actual gate
  }
  
}


void GarbledCircuit::set_const_key(byte c, const Bytes &key)
{
  assert(c == 0 || c == 1); // wire for constant 0 or 1
  Bytes tmp = key;
  tmp.resize(16);
  m_const_wire[c] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
}


// TODO: figure out what this is intended to do
const Bytes GarbledCircuit::get_const_key(byte c, byte b)
{
  assert(c == 0 || c == 1); // wire for constant 0 or 1
  assert(b == 0 || b == 1); // with bit value 0 or 1
  Bytes tmp(16);
  
  tmp.resize(16);
  if (b)
    {
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), _mm_xor_si128(m_R, m_const_wire[c]));
    }
  else
    {
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), m_const_wire[c]);
    }
  
  tmp.resize(Env::key_size_in_bytes());
  return tmp;
}




#endif
