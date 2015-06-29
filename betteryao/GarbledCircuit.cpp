#ifndef GARBLED_CIRCUIT_CPP
#define GARBLED_CIRCUIT_CPP

#include "GarbledCircuit.h"

// notes: still need functions to load the circuit file and/or
// set the loaded circuit file to a field in the garbled circuit for access
 
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
void GarbledCircuit::evaluate_Circuit(){
  // gen and eval versions
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
  // returns void pointer, which is pointer to a key?
  // use this one to call the Garbled Circuit object again
  GarbledCircuit &cct = *reinterpret_cast<GarbledCircuit*>(get_external_circuit(st));
  
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


void GarbledCircuit::init_Generation_Circuit(const std::vector<Bytes>* gen_keys, const std::vector<Bytes> * evl_keys){
  Prng prng = Prng();
  Bytes tmp;
  
  // initialize the value of R
  tmp = prng.rand_bits(Env::k());
  tmp.resize(16,0);
  tmp.set_ith_bit(0,1);
  m_R = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  // initialize the constant wires
  tmp = prng.rand_bits(Env::k());
  tmp.resize(16,0);
  m_const_wire[0] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  tmp = prng.rand_bits(Env::k());
  tmp.resize(16,0);
  m_const_wire[1] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  // set the input keys
  m_gen_inputs = gen_keys;
  m_evl_inputs = evl_keys;

}


void GarbledCircuit::init_Evaluation_Circuit(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys){
  
  m_gen_inputs = gen_keys;
  m_evl_inputs = evl_keys;

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
  return get_Gen_Input(2*input_idx + parity);
}

Bytes GarbledCircuit::get_Evl_Key(uint32_t input_idx, uint32_t parity){
  assert(parity == 0 || parity == 1);
  return get_Evl_Input(2*input_idx + parity);
}

Bytes GarbledCircuit::get_Gen_Input(uint32_t idx){
  return m_gen_inputs[idx];
}

Bytes GarbledCircuit::get_Evl_Input(uint32_t idx){
  return m_evl_inputs[idx];
}

/**
   Must be able to access the current wire values at indices from the wire table
 */

Bytes GarbledCircuit::get_Wire_Value(uint32_t idx){
  

}


void * GarbledCircuit::gen_Next_Gate(PCFGate *current_gate){
  // somehow the PCFState should be available through the m_st pointer
  
  if(current_gate->tag == TAG_INPUT_A){
    
  } else if (current_gate->tag == TAG_INPUT_B){

  } else if (current_gate->tag == TAG_OUTPUT_A) {

  } else if (current_gate->tag == TAG_OUTPUT_B){
    
  } else {
    // actual gate
  }
  
}

void * GarbledCircuit::evl_Next_Gate(PCFGate *current_gate){
  if(current_gate->tag == TAG_INPUT_A){
    
  } else if (current_gate->tag == TAG_INPUT_B){

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
