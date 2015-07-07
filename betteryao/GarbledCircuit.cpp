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

 
GarbledCircuit::GarbledCircuit(): m_gate_index(0), m_bob_out_ix(0),m_alice_out_ix(0) {

  // initialize the key Mask
  Bytes tmp(16);
  for(size_t ix = 0; ix< Env::k(); ix++){
    tmp.set_ith_bit(ix,1);
  }
  m_clear_mask = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));

  m_alice_out.clear();
  m_bob_out.clear();

}

// the big cahoona
// this is called by the BetterYao protocol object
// and runs the circuit till the end.
void GarbledCircuit::generate_Circuit(){
  // gen and eval versions
  
  /*
  while(get_next_gate(m_st)){
    GEN_SEND(garbling_bufr);
  }
  */
}

// the big cahoona
// this is called by the BetterYao protocol object
// and runs the circuit till the end.
void GarbledCircuit::evaluate_Circuit(){
  // gen and eval versions

  // Bytes bufr;
  do {
    // do things?
    //    garbling_bufr = EVL_RECV();
  }  
  while(get_next_gate(m_st));

}

/**
   These functions serve as intermediaries
   to get back to the Garbled Circuit object
   since we go through such crazy calling loops with 
   the pcf interpreter struct
 */
void * gen_next_gate(PCFState *st, PCFGate *current_gate){

  //fprintf(stdout,"gen next gate\n");

  // returns void pointer, which is pointer to a key?
  // use this one to call the Garbled Circuit object again
  GarbledCircuit &cct = *reinterpret_cast<GarbledCircuit*>(get_external_circuit(st));

  //fprintf(stdout,"gen next gate\n");
  // now, call the appropriate function from cct

  //Bytes garbled_gate = Bytes(0);

  return cct.gen_Next_Gate(current_gate);
 }

void * evl_next_gate(PCFState *st, PCFGate *current_gate){
  // returns void pointer, which is pointer to a key?
  // use this one to call the Garbled Circuit object again
  GarbledCircuit &cct = *reinterpret_cast<GarbledCircuit*>(get_external_circuit(st));
  
  // now, call the appropriate function from cct
  return cct.evl_Next_Gate(current_gate);
 }


void GarbledCircuit::init_Generation_Circuit(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys, Bytes & circuit_seed, const Bytes & select_bits, const Bytes R, const Bytes & zero_key, const Bytes & one_key ){
  
  Bytes tmp;
  
  m_prng.seed_rand(circuit_seed);
  //  m_select_bits.resize()
  m_select_bits.insert(m_select_bits.begin(),select_bits.begin(),select_bits.begin() + select_bits.size());
  // m_select_bits = permutation_bits;

  fprintf(stdout, "Gen select bits: %s\n",select_bits.to_hex().c_str());

  // initialize the value of R
  //tmp = m_prng.rand_bits(Env::k());
  //tmp.resize(16,0);
  //tmp.set_ith_bit(0,1);
  
  // we initialize tmp with zeros for this so that 
  // the value of R doesn't get extra nonsense at the end
  // from whatever was sitting around in memory
  tmp.resize(16,0);
  tmp.insert(tmp.begin(),R.begin(),R.begin()+Env::k()/8);
  m_R = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));
  fprintf(stdout,"R value: %s\n",R.to_hex().c_str());
  //  print128_num(m_R);
  

  // initialize the constant wires
  tmp.clear();
  tmp.resize(Env::k());
  tmp.insert(tmp.end(),zero_key.begin(), zero_key.begin()+zero_key.size());
  tmp.resize(16,0);
  m_const_wire[0] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  tmp.clear();
  tmp.resize(Env::k());
  tmp.insert(tmp.end(),one_key.begin(), one_key.begin()+one_key.size());
  tmp.resize(16,0);
  m_const_wire[1] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));
  

  // set the input keys
  set_Input_Keys(gen_keys, evl_keys);

}

void GarbledCircuit::set_Gen_Circuit_Functions(){
  set_key_copy_function(m_st, copy_key);
  set_key_delete_function(m_st, delete_key);
  set_callback(m_st,gen_next_gate);
  
}


void GarbledCircuit::set_Evl_Circuit_Functions(){
  set_key_copy_function(m_st, copy_key);
  set_key_delete_function(m_st, delete_key);
  set_callback(m_st,evl_next_gate);
  
}


void GarbledCircuit::set_Input_Keys(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys){
  m_gen_inputs = gen_keys;
  m_evl_inputs = evl_keys;

  for(int i = 0; i < m_evl_inputs->size();i++){
    fprintf(stdout,"evl key %i: %s\n",i, (*m_evl_inputs)[i].to_hex().c_str());
  }
}

void GarbledCircuit::init_Evaluation_Circuit(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys, const Bytes & evl_input, const Bytes &zero_key, const Bytes & one_key){
  set_Input_Keys(gen_keys, evl_keys);

  //m_select_bits = evl_input;
  m_select_bits.insert(m_select_bits.begin(),evl_input.begin(),evl_input.begin() + evl_input.size());

  fprintf(stdout,"***Eval Input (select bits) *** %s\t%s\n",evl_input.to_hex().c_str(),m_select_bits.to_hex().c_str());
  //std::cout << "Gen first input key: " << (*gen_keys)[0].to_hex() << std::endl;

  Bytes tmp;
  // initialize the constant wires
  tmp.clear();
  tmp.resize(Env::k());
  tmp.insert(tmp.end(),zero_key.begin(), zero_key.begin()+zero_key.size());
  tmp.resize(16,0);
  m_const_wire[0] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));

  tmp.clear();
  tmp.resize(Env::k());
  tmp.insert(tmp.end(),one_key.begin(), one_key.begin()+one_key.size());
  tmp.resize(16,0);
  m_const_wire[1] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));


}


void * GarbledCircuit::get_Const_Wire(uint32_t i){
  assert(i == 0 || i == 1);
  return &(m_const_wire[i]);
}

/*
void GarbledCircuit::garble_Gate(){


}
*/

/*
void * GarbledCircuit::garble_On_Keys(void * x_key, void * y_key){
  __m128i X[2], Y[2], Z[2];
  // __m128i aes_key[2], aes_plaintext, aes_ciphertext;
  Bytes tmp(16,0);
  
}
*/

/**
   Input Key Accessor Functions 
   get-Key functions are intended for Gen (although they call the get-Input) functions
   get-Input functions are intended for Eval
 */

Bytes GarbledCircuit::get_Gen_Key(uint32_t input_idx, uint32_t parity){
  assert(parity == 0 || parity == 1);
  fprintf(stdout,"Get Gen Key: %u %u\n",input_idx,parity);
  return get_Gen_Input(2*input_idx + parity);
}

Bytes GarbledCircuit::get_Evl_Key(uint32_t input_idx, uint32_t parity){
  assert(parity == 0 || parity == 1);
  return get_Evl_Input(2*input_idx + parity);
}

Bytes GarbledCircuit::get_Gen_Input(uint32_t idx){
  fprintf(stdout,"Get Gen Input idx: %u\n",idx);
  Bytes ret = (*m_gen_inputs)[idx];
  //fprintf(stdout, "Gen's Input Key: %s\n", ret.to_hex().c_str());
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
  // Gen uses this to look up permutation bits
  // Eval uses this function to figure out what input key to decrpyt
  //   fprintf(stdout,"%s get %u input parity: %u\n",m_select_bits.to_hex().c_str(),idx,m_select_bits.get_ith_bit(idx));
  
 // fprintf(stdout,"HELP ME PLEASE\n");
  //for(int i = 0; i <= idx; i++){
  //  fprintf(stdout,"%x",m_select_bits.get_ith_bit(i));
  //}
  
 return m_select_bits.get_ith_bit(idx);

}

void GarbledCircuit::generate_Random_Key(__m128i & destination){
  Bytes tmp = m_prng.rand_bits(Env::k());
  tmp.resize(16,0);
  destination = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
}

void print128_num(__m128i var)
{
    uint16_t *val = (uint16_t*) &var;
    fprintf(stdout,"Numerical: %x %x %x %x %x %x %x %x \n", 
           val[0], val[1], val[2], val[3], val[4], val[5], 
           val[6], val[7]);
}


/**
  throughout garbling, Gen will maintain the zero-semantic keys for each wire
  and use them to generate all of the ciphertexts (find the 1-semantics by XOR with R)
  that are sent to Eval
 */
void * GarbledCircuit::gen_Next_Gate(PCFGate *current_gate){
  // the PCFState should be available through the m_st pointer
  
  static __m128i current_key; // must be static to return it
  
  if(current_gate->tag == TAG_INPUT_A){
    generate_Alice_Input(current_gate, current_key);
    return &current_key;
    
  }  else if (current_gate->tag == TAG_INPUT_B){
    generate_Bob_Input(current_gate, current_key);
    return &current_key;

  } else if (current_gate->tag == TAG_OUTPUT_A) {

    fprintf(stdout,"Alice Output!\n");
    generate_Gate(current_gate, current_key);
    generate_Alice_Output(current_gate,current_key);
    m_gate_index++;
    return &current_key;

  } else if (current_gate->tag == TAG_OUTPUT_B){

    fprintf(stdout,"Bob Output!\n");
    generate_Gate(current_gate,current_key);
    generate_Bob_Output(current_gate, current_key);
    m_gate_index++;
    return &current_key;

  } else {

    // actual gate
    generate_Gate(current_gate, current_key);
    m_gate_index++;
    return &current_key; // get_Gen_Key or create a random nonce
  }
}


void * GarbledCircuit::evl_Next_Gate(PCFGate *current_gate){
  
  static __m128i current_key; // must be static to return it
  
  if(current_gate->tag == TAG_INPUT_A){
    
    evaluate_Alice_Input(current_gate, current_key);
    return &current_key;
    
  } else if (current_gate->tag == TAG_INPUT_B){
    evaluate_Bob_Input(current_gate,current_key);
    return &current_key; // get_Gen_Key or create a random nonce
    
  } else if (current_gate->tag == TAG_OUTPUT_A) {

    fprintf(stdout,"Alice Output!\n");
    evaluate_Gate(current_gate,current_key);
    evaluate_Alice_Output(current_gate,current_key);
    m_gate_index++;
    return &current_key;

  } else if (current_gate->tag == TAG_OUTPUT_B){

    fprintf(stdout,"Bob Output!\n");
    evaluate_Gate(current_gate, current_key);
    evaluate_Bob_Output(current_gate, current_key);
    m_gate_index++;
    return &current_key;
    
} else {
    
    evaluate_Gate(current_gate,current_key);
    m_gate_index++;
    return &current_key;
  }
}


void GarbledCircuit::xor_Gate(PCFGate* current_gate, __m128i &current_key){
  
  __m128i key1 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1));
  __m128i key2 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2));

  current_key = _mm_xor_si128
    (
     key1,
     key2
     //*reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1)),
     //*reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2))
     );
  
  fprintf(stdout,"XOR Gate\n");
  fprintf(stdout,"input keys: %i %i (parity %x %x) \n",current_gate->wire1, current_gate->wire2, _mm_extract_epi8(key1,0)&0x01, _mm_extract_epi8(key2,0)&0x01);
  print128_num(key1);
  print128_num(key2);
  //  fprintf(stdout,"derived output key:\n");
  // print128_num(current_key);
  //  fprintf(stdout,"\n");
  
}

  uint32_t GarbledCircuit::increment_index(){
  return m_gate_index++;
}


void GarbledCircuit::generate_Bob_Input(PCFGate* current_gate, __m128i &current_key){
  // Gen's input keys have already been generated and determined
  // by his input keys
  // here, we return the proper input key encoding 0
  // to do this, we need to know both his permutation bits and his input bits
  fprintf(stdout,"\nBob/Gen Input Gate\n");
  
  uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index
  
  fprintf(stdout,"Next Wire Index: %u\n", gen_input_idx);
    
  //    Bytes gen_input = get_Gen_Key(gen_input_idx, get_Input_Parity(gen_input_idx));
  Bytes gen_input = get_Gen_Input(2*gen_input_idx + get_Input_Parity(gen_input_idx));
  
  Bytes gen_xor_input = get_Gen_Input(2*gen_input_idx + (1 ^ get_Input_Parity(gen_input_idx)));
    
  fprintf(stdout, "Gen Input Keys: %s %s::(%x) \n", gen_input.to_hex().c_str(), gen_xor_input.to_hex().c_str(), get_Input_Parity(gen_input_idx)); 
  
  fprintf(stdout,"Save to 128 Bit\n");
  save_Key_to_128bit(gen_input,current_key);
  __m128i xor_key = _mm_xor_si128(current_key,m_R);
  print128_num(current_key);
  print128_num(xor_key);
  fprintf(stdout,"\n\n");

  // we need to set the garbling buffer to something (or nothing)
  // because it will be sent by default between gates
  m_garbling_bufr = Bytes(0);
}

void GarbledCircuit::evaluate_Bob_Input(PCFGate* current_gate, __m128i &current_key){
    // here, Gen's input keys have been sent to Eval
  // she needs only to assign the input key to the wire

  fprintf(stdout,"\nBob/Gen Input Gate\n");
  
  uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index
  
  fprintf(stdout,"Next Wire Index: %u\n", gen_input_idx);
  
  Bytes gen_input = get_Gen_Input(gen_input_idx);
  
  fprintf(stdout,"Gen Input: %s\n", gen_input.to_hex().c_str());
  save_Key_to_128bit(gen_input,current_key);
  fprintf(stdout,"Parity: %x (%x)\n",_mm_extract_epi8(current_key,0)&0x01,gen_input.get_ith_bit(0));
  print128_num(current_key);
  // and current_key is available up the call stack
}

void GarbledCircuit::generate_Alice_Input(PCFGate* current_gate, __m128i &current_key){
    // Here, Eval already has her input keys, but Gen doesn't know
    // which to use. So he will generate a new key and its XOR-offser,
    // encrypt it using the respective Eval input keys, and send the keys to Eval.
    // She will decrypt the proper one and use it 
    
    __m128i output_keys[2];
    Bytes new_zero_key;
    Bytes tmp; // used for loading input keys
    new_zero_key = m_prng.rand_bits(Env::k());
    new_zero_key.resize(16,0);
    new_zero_key.set_ith_bit(0,0);

    // save this zero-key for future access
    current_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&new_zero_key[0]));

    fprintf(stdout, "\nAlice/Evl Input:: new current key:\n");
    print128_num(current_key);
    __m128i xor_key = _mm_xor_si128(current_key, m_R);
    fprintf(stdout, "Bob xor-input key:\n");
    print128_num(xor_key);

    // load the first output key
    tmp = get_Evl_Key(current_gate->wire1,0);
    tmp.resize(16,0);
    output_keys[0] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
    fprintf(stdout,"key1: %s\n",tmp.to_hex().c_str());
    //print128_num(output_keys[0]);
    //tmp.clear();

    // load the second output key
    tmp = get_Evl_Key(current_gate->wire1,1); // wire1 = wire2
    tmp.resize(16,0);
    output_keys[1] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
    fprintf(stdout,"key2: %s\n",tmp.to_hex().c_str());
    //print128_num(output_keys[1]);
    
    // encrypt the output keys so that Eval can get only one of them
    output_keys[0] = _mm_xor_si128(output_keys[0], current_key);
    output_keys[1] = _mm_xor_si128(output_keys[1], _mm_xor_si128(current_key,m_R));
    
    // now send the output keys to Eval for decryption
    // put the keys into the garbling buffer and send them
    m_garbling_bufr.clear();
    tmp.clear();
    tmp.resize(16,0);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), output_keys[0]);
    m_garbling_bufr.insert(m_garbling_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());

    tmp.clear();
    tmp.resize(16,0);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), output_keys[1]);
    m_garbling_bufr.insert(m_garbling_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
    
    // should have received garbling bufr from the high level protocol function
    assert(m_garbling_bufr.size() == 2*Env::key_size_in_bytes());
    
    // quick print check
    fprintf(stdout,"Alice/Eval Key: ");
    print128_num(current_key);
    fprintf(stdout,"Alice/Eval Xor: ");
    print128_num(xor_key);
    fprintf(stdout,"Garbling buffer transmission: %s\n",m_garbling_bufr.to_hex().c_str());
    fprintf(stdout,"\n");
    
}

void GarbledCircuit::evaluate_Alice_Input(PCFGate* current_gate, __m128i &current_key){
   // here, Eval already knows which input key she wants to use
  // she selects it and assigns it to her wire value
  
  uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index
  Bytes evl_input = get_Evl_Input(gen_input_idx);
  
  // must send 2 ciphertexts
  assert(m_garbling_bufr.size() == 2*Env::key_size_in_bytes());
  
  // find the input ciphertext that corresponds to Eval's chosen bit
  uint32_t bit = get_Input_Parity(current_gate->wire1);
  assert(bit == 0 || bit == 1);
  
  fprintf(stdout,"Evaluate Alice Input\n");

  fprintf(stdout,"\nreceived garbling buffer: %s\n",m_garbling_bufr.to_hex().c_str());
  fprintf(stdout,"selecting and decrypting Eval Input Key: %x (input wire %i)\n",bit,current_gate->wire1);
  
  // now select the one of two ciphertexts to decrypt
  Bytes encrypted_input;
  encrypted_input.insert(encrypted_input.end(),
                         m_garbling_bufr.begin() + bit * Env::key_size_in_bytes(),
                         m_garbling_bufr.begin() + Env::key_size_in_bytes() + bit * Env::key_size_in_bytes());
  
  fprintf(stdout,"decrypting encrypted buffer: %s\n",encrypted_input.to_hex().c_str()); 
  fprintf(stdout,"eval buffer decryption key: %s\n",evl_input.to_hex().c_str());

  assert(evl_input.size() == Env::key_size_in_bytes());
  assert(evl_input.size() == encrypted_input.size());
  
  evl_input = evl_input ^ encrypted_input;
  save_Key_to_128bit(evl_input,current_key);
  
  fprintf(stdout,"Alice/Eval Key: (parity %x)\n", _mm_extract_epi8(current_key,0)&0x01);
  // a quick print check
  print128_num(current_key);
  fprintf(stdout,"\n");

}

void GarbledCircuit::evaluate_Alice_Output(PCFGate* current_gate, __m128i &current_key){
    // make sure the output buffer is big enough
    if(m_alice_out.size()*8 <= m_alice_out_ix){
      m_alice_out.resize((m_alice_out.size()+1)*2,0); // grow by doubling for less work
    }

    uint8_t out_bit = _mm_extract_epi8(current_key, 0) & 0x01;
    // out_bit ^= *m_in_bufr_ix;
    m_alice_out.set_ith_bit(m_alice_out_ix, out_bit);
    // m_in_bufr_ix++;
    m_alice_out_ix++;

}
void GarbledCircuit::evaluate_Bob_Output(PCFGate* current_gate, __m128i &current_key){
    if (m_bob_out.size()*8 <= m_bob_out_ix)
      {
        // dynamically grow output array by doubling
        m_bob_out.resize((m_bob_out.size()+1)*2, 0);
      }
    
    uint8_t out_bit = _mm_extract_epi8(current_key, 0) & 0x01;
    
    
    //out_bit ^= m_in_bufr_ix;
    m_bob_out.set_ith_bit(m_bob_out_ix, out_bit);
    //m_in_bufr_ix++;
    m_bob_out_ix++;
    

}
void GarbledCircuit::generate_Alice_Output(PCFGate* current_gate, __m128i &current_key){
    
    //m_gen_out.push_back(_mm_extract_epi8(current_key,0)&0x01);
    //m_garbling_bufr.push_back(_mm_extract_epi8(current_key,0)&0x01);
    
    // gen moves through all of his output bits and saves the parity of
    // the current zero keys, to be transmitted to Evl all together after
    // protocol execution
    if(m_alice_out.size()*8 <= m_alice_out_ix){
      m_alice_out.resize((m_alice_out.size()+1)*2,0); // grow by doubling for less work
    }

    uint8_t out_bit = _mm_extract_epi8(current_key, 0) & 0x01;
    // gen out contains the permutation bits of the output keys
    m_alice_out.set_ith_bit(m_alice_out_ix, out_bit);
    m_alice_out_ix++;
}


void GarbledCircuit::generate_Bob_Output(PCFGate* current_gate, __m128i &current_key){
  //m_garbling_bufr.push_back(_mm_extract_epi8(current_key,0)&0x01);
  //    m_bob_out_ix++;
  
  if (m_bob_out.size()*8 <= m_bob_out_ix)
    {
     // dynamically grow output array by doubling
     m_bob_out.resize((m_bob_out.size()+1)*2, 0);
   }
  
  uint8_t out_bit = _mm_extract_epi8(current_key, 0) & 0x01;
  
  //out_bit ^= m_in_bufr_ix;
  m_bob_out.set_ith_bit(m_bob_out_ix, out_bit);
  //m_in_bufr_ix++;
  m_bob_out_ix++;
  
}


void GarbledCircuit::generate_Gate(PCFGate* current_gate, __m128i &current_key){
  Bytes tmp;
  
#ifdef FREE_XOR
  if(current_gate->truth_table == 0x06){ // if XOR gate
    xor_Gate(current_gate, current_key);
  } else {
#endif
    
    // we garble this gate
    // first, get the inputs
    
    // X and Y are input, Z is output
    __m128i X[2], Y[2], Z[2];
    // we need a couple of 128-bit variables
    __m128i aes_plaintext, aes_ciphertext;
    __m128i aes_key[2];
    
    Bytes gate_out_bufr;
    gate_out_bufr.resize(3*Env::key_size_in_bytes());
      
    uint8_t bit;
    
    // load the (zero-key) inputs from the PCF state container
    X[0] = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1));
    Y[0] = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2));
    // and derive their XOR-complements
    X[1] = _mm_xor_si128(X[0], m_R); // X[1] = X[0] ^ R
    Y[1] = _mm_xor_si128(Y[0], m_R); // Y[1] = Y[0] ^ R
    
    fprintf(stdout,"garbling input keys: %i %i (parity: %i %i) \n",
            current_gate->wire1,
            current_gate->wire2,
            _mm_extract_epi8(X[0],0) &0x01,
            _mm_extract_epi8(Y[0],0) &0x01);
    print128_num(X[0]);
    print128_num(X[1]);
    print128_num(Y[0]);
    print128_num(Y[1]);

    // and get the permutation bits (tells if each zero-key has a 1 on the end)
    const uint8_t perm_x = _mm_extract_epi8(X[0],0) & 0x01;
    const uint8_t perm_y = _mm_extract_epi8(Y[0],0) & 0x01;
    // construct the permutation index
    const uint8_t de_garbled_ix = (perm_y << 1)|perm_x; // wire1 + 2*wire2
    
    // now we are ready to garble

    // load the plaintext to be encrypted
    aes_plaintext = _mm_set1_epi64x(m_gate_index);
    // and load the garbling keys
    
    aes_key[0] = _mm_load_si128(X+perm_x);
    aes_key[1] = _mm_load_si128(Y+perm_y);
    //aes_key[0] = _mm_load_si128(X); // this is still the zero key
    //aes_key[1] = _mm_load_si128(Y); // this is still the zero key
    
    // now run the key derivation function using the keys and the gate index
    KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
    
    // clear extra bits at the front
    aes_ciphertext = _mm_and_si128(aes_ciphertext, m_clear_mask);
    // aes_ciphertext now contains the encryption of both zero-keys
    
    // this bit discovers what the semantics of the next wire will be
    bit = (current_gate->truth_table >> (3-de_garbled_ix)) & 0x01;
    
#ifdef GRR
    // GRR technique: using zero entry's key as one of the output keys
    // the output key is the encrypted gate index
    // this puts the zero key in Z[0] and the one key in Z[1]
    
    // we want the new zero-key to be in Z[0]
    _mm_store_si128(Z+bit, aes_ciphertext);
    // and other output is an offset
    Z[1 - bit] = _mm_xor_si128(Z[bit], m_R);
    //Z[1] = _mm_xor_si128(Z[0], m_R);
    // and load it into current_key to return it to the circuit's state container
    current_key = _mm_load_si128(Z);
    // this way we return the new zero-key
    
    /*
    fprintf(stdout,"output keys, destination %x: (parity %i %i) \n",current_gate->reswire,_mm_extract_epi8(Z[0],0) & 0x01,
      _mm_extract_epi8(Z[1],0) & 0x01);
    print128_num(Z[0]);
    //fprintf(stdout,"parity: %x\n",);
    print128_num(Z[1]);
    fprintf(stdout,"\n\n");
    //fprintf(stdout,"parity: %x\n",_mm_extract_epi8(Z[1],0) & 0x01);
    */

#else
    // practically, this code is obsolete. we should be using GRR

    // otherwise we generate a new random key for the 0th entry
    tmp = m_prng.rand_bits(Env::k());
    tmp.resize(16, 0);
    Z[0] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
    Z[1] = _mm_xor_si128(Z[0], m_R);
    
    // new ciphertext becomes the xor of this and the output key
    aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
    // we put it back in tmp so that we can communicate it to Evl
    
    _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);
    // now we have to save this key to send to Evl
    // cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
#endif
     
    // now we encrypt the next several truth table entries before sending them to Eval
    // the first entry is (X[1-x],Y[y])
    // so we take the XOR-offset of X
    // Bytes tmp;
    tmp.resize(16,0);
    
    // we just used X[0] as aex_key[0],
    // so now take the XOR with R to get X[1-x]
    aes_key[0] = _mm_xor_si128(aes_key[0], m_R);

    // and run through the KDF
    // key derivation function
    KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
    // and mask to remove unnecessary bits
    aes_ciphertext = _mm_and_si128(aes_ciphertext, m_clear_mask);

    // recover the permutation bit
    bit = (current_gate->truth_table>>(3-(0x01^de_garbled_ix)))&0x01;
    // encrypt the key using the ciphertext and load it into tmp to send
    aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);
    // and we add tmp to the buffer we're sending to Eval
    
    m_garbling_bufr.insert(m_garbling_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
    
    
    // encrypt the 2nd entry : (X[x], Y[1-y])
    aes_key[0] = _mm_xor_si128(aes_key[0], m_R);
    aes_key[1] = _mm_xor_si128(aes_key[1], m_R);
    
    KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
    aes_ciphertext = _mm_and_si128(aes_ciphertext, m_clear_mask);
    //bit = current_gate.m_table[0x02^de_garbled_ix];
    bit = (current_gate->truth_table>>(3-(0x02^de_garbled_ix)))&0x01;
    aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
    
    _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);
    
    m_garbling_bufr.insert(m_garbling_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
    
    //cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
    
    // encrypt the 3rd entry : (X[1-x], Y[1-y])
    aes_key[0] = _mm_xor_si128(aes_key[0], m_R);
    
    KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
    aes_ciphertext = _mm_and_si128(aes_ciphertext, m_clear_mask);
    //bit = current_gate.m_table[0x03^de_garbled_ix];
    bit = (current_gate->truth_table>>(3-(0x03^de_garbled_ix)))&0x01;
    aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);
    
    m_garbling_bufr.insert(m_garbling_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
    
    //m_garbling_bufr += tmp;
    //fprintf(stdout,"Input Gates to Garble Function: %x %x \n", current_gate->wire1, current_gate->wire2);
    //fprintf(stdout, "Gen Send Gate: %s\n", m_garbling_bufr.to_hex().c_str()); 
    
    //GEN_SEND(garbling_bufr);
    // m_garbling_bufr.clear();
    //cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
    
#ifdef FREE_XOR
  }
#endif
  
  
  fprintf(stdout,"derived keys, destination %i: (parity %x) \n", current_gate->reswire, _mm_extract_epi8(current_key,0) & 0x01);
  print128_num(current_key);
  __m128i xor_offset = _mm_xor_si128(current_key, m_R);
  //  fprintf(stdout,"derived xor offset: (parity %x)\n", _mm_extract_epi8(xor_offset,0) & 0x01 );
  print128_num(xor_offset);
    
  fprintf(stdout,"\n");
  
    
  /*
  fprintf(stdout,"derived keys, destination %x: (parity %i %i) \n",current_gate->reswire,_mm_extract_epi8(Z[0],0) & 0x01,
          _mm_extract_epi8(Z[1],0) & 0x01);
  print128_num(Z[0]);
  //fprintf(stdout,"parity: %x\n",);
  print128_num(Z[1]);
  fprintf(stdout,"\n\n");
  //fprintf(stdout,"parity: %x\n",_mm_extract_epi8(Z[1],0) & 0x01);
  */



}



void GarbledCircuit::evaluate_Gate(PCFGate* current_gate, __m128i &current_key){
   
#ifdef FREE_XOR
  if (current_gate->truth_table == 0x06)
    {
      xor_Gate(current_gate, current_key);
    } else {
#endif    
    
    // non-XOR gate
    __m128i aes_key[2], aes_plaintext, aes_ciphertext;
    Bytes tmp;
    __m128i a;
    Bytes::const_iterator it;
    
    //fprintf(stdout,"receive: \n");
    
    //fprintf(stdout,"EVL Receive Gate: %s\n",m_garbling_bufr.to_hex().c_str());
    
    //assert(m_garbling_bufr.size() == 3 * Env::key_size_in_bytes());
    
    aes_plaintext = _mm_set1_epi64x(m_gate_index);
    
    aes_key[0] = *reinterpret_cast<__m128i*>(get_wire_key(m_st,current_gate->wire1));
    aes_key[1] = *reinterpret_cast<__m128i*>(get_wire_key(m_st,current_gate->wire2));
    
    const uint8_t perm_x = _mm_extract_epi8(aes_key[0], 0) & 0x01;
    const uint8_t perm_y = _mm_extract_epi8(aes_key[1], 0) & 0x01;
    
    fprintf(stdout,"garbling input keys: %i %i (parity %i %i)\n",current_gate->wire1, current_gate->wire2,perm_x, perm_x);
    print128_num(aes_key[0]);
    print128_num(aes_key[1]);
    //    fprintf(stdout,"\n");
    
    KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
    aes_ciphertext = _mm_and_si128(aes_ciphertext, m_clear_mask);
    uint8_t garbled_ix = (perm_y<<1)|perm_x;
    
#ifdef GRR
    if (garbled_ix == 0)
      {
        current_key = _mm_load_si128(&aes_ciphertext);
      }
    else
      {
        //Bytes::const_iterator it = cct.m_in_bufr_ix+(garbled_ix-1)*Env::key_size_in_bytes();
        it = m_garbling_bufr.begin() + (garbled_ix-1)*Env::key_size_in_bytes();
        
        tmp.assign(it, it+Env::key_size_in_bytes());
        tmp.resize(16, 0);
        a = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
        current_key = _mm_xor_si128(aes_ciphertext, a);
      }
#else
    it = m_garbling_bufr.begin() + (garbled_ix)*Env::key_size_in_bytes();
    tmp.assign(it, it+Env::key_size_in_bytes());
    tmp.resize(16, 0);
    current_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
    current_key = _mm_xor_si128(current_key, aes_ciphertext);
#endif
    
#ifdef FREE_XOR
  }
#endif
  fprintf(stdout,"derived key, destination %i: (parity %x) \n", current_gate->reswire, _mm_extract_epi8(current_key,0) & 0x01);
  print128_num(current_key);
  fprintf(stdout,"\n");
  // current_key will be returned
}



  /*
  void GarbledCircuit::genHalfGate(PCFGate* current_Gate, __m128i &current_key){
  


}
  

  
void GarbledCircuit::evlHalfGate(PCFGate* current_gate, __m128i &current_key){
  Bytes TG, TE;

  assert(m_garbling_bufr.size() == 2*Env::key_size_in_bytes());
    
  // TG is always sent first, and then TE
  // where TG is the single row transmitted for the Generator's half-gate
  // and TE is the single row transmitted for the Evaluator's half-gate
  TG.insert(TG.end(),m_garbling_bufr.begin(),m_garbling_bufr.begin()+Env::key_size_in_bytes);
  TE.insert(TE.end(),m_garbling_bufr.begin()+Env::key_size_in_bytes(),m_garbling_bufr.begin()+2*Env::key_size_in_bytes());

  // load the inputs from the PCF state container
  
  __m128i Wa, Wb;
  
  Wa = *reinterpret_cast<__m128i*>(get_wire_key(m_st,current_gate->wire1));
  Wb = *reinterpret_cast<__m128i*>(get_wire_key(m_st,current_gate->wire2));

  Bytes WireA, WireB;
  
  WireA.resize(16,0);
  WireB.resize(16,0);
  
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&WireA[0]), Wa);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&WireB[0]), Wb);
  
  // now get the select bits

  byte sa,sb;
  sa = WireA.get_ith_bit(0); sb = WireB.get_ith_bit(0);

  uint32_t j, k;
  
  j = increment_index();
  k = increment_index();
  
  // ready for the KDF (or hash function)
  // here we use a single-key AES call
  
  
  
  
}
  */


void GarbledCircuit::set_const_key(byte c, const Bytes &key)
{
  assert(c == 0 || c == 1); // wire for constant 0 or 1
  Bytes tmp = key;
  tmp.resize(16,0);
  m_const_wire[c] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
}


Bytes GarbledCircuit::get_garbling_bufr(){
  //std::cout << "get garbling bufr: " << m_garbling_bufr.to_hex() << std::endl;
  //fprintf(stdout,"get garbling bufr: %s\n", m_garbling_bufr.to_hex().c_str()); 
  return m_garbling_bufr;
}

void GarbledCircuit::set_garbling_bufr(Bytes buf){
  m_garbling_bufr = buf;
  
  //fprintf(stdout,"setting garbling bufr: %s %s\n",buf.to_hex().c_str(),m_garbling_bufr.to_hex().c_str());
}

void GarbledCircuit::clear_garbling_bufr(){
  m_garbling_bufr.clear();
}

// TODO: figure out what this is intended to do
/*
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
*/

 Bytes GarbledCircuit::get_alice_out(){
  return m_alice_out;
}

 Bytes GarbledCircuit::get_bob_out(){
  return m_bob_out;
}

 void GarbledCircuit::trim_output_buffers(){
  if(m_alice_out_ix>0)
    m_alice_out.resize(m_alice_out_ix/8,0);
  if(m_bob_out_ix>0)
    m_bob_out.resize(m_bob_out_ix/8,0);
}

#endif
