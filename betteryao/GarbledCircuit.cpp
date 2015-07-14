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

  // returns void pointer, which is pointer to a key?
  // use this one to call the Garbled Circuit object again
  GarbledCircuit &cct = *reinterpret_cast<GarbledCircuit*>(get_external_circuit(st));

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
  m_select_bits.insert(m_select_bits.begin(),select_bits.begin(),select_bits.begin() + select_bits.size());

  fprintf(stdout, "Gen select bits: %s\n",select_bits.to_hex().c_str());

  // initialize the value of R
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
  //  fprintf(stdout,"Get Gen Key: %u %u\n",input_idx,parity);
  return get_Gen_Input(2*input_idx + parity);
}

Bytes GarbledCircuit::get_Evl_Key(uint32_t input_idx, uint32_t parity){
  assert(parity == 0 || parity == 1);
  return get_Evl_Input(2*input_idx + parity);
}

Bytes GarbledCircuit::get_Gen_Input(uint32_t idx){
  // fprintf(stdout,"Get Gen Input idx: %u\n",idx);
  Bytes ret = (*m_gen_inputs)[idx];
  //fprintf(stdout, "Gen's Input Key: %s\n", ret.to_hex().c_str());
  return ret;
}

Bytes GarbledCircuit::get_Evl_Input(uint32_t idx){
  if(idx < m_evl_inputs->size()){
    return (*m_evl_inputs)[idx];
  } else{
    Bytes tmp;
    tmp.resize(Env::k()/8,0);
    return tmp;
  }
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
  if(idx < m_select_bits.size()*8){
    return m_select_bits.get_ith_bit(idx);
  } else{
    return 0;
  }
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
    // fprintf(stdout, "Alice Input!");
    generate_Alice_Input(current_gate, current_key);
    return &current_key;
    
  }  else if (current_gate->tag == TAG_INPUT_B){
    //  fprintf(stdout, "Bob Input!");
    generate_Bob_Input(current_gate, current_key);
    return &current_key;

  } else if (current_gate->tag == TAG_OUTPUT_A) {
    
    // fprintf(stdout,"Alice Output!\n");
    generate_Gate(current_gate, current_key);
    generate_Alice_Output(current_gate,current_key);
    m_gate_index++;
    return &current_key;

  } else if (current_gate->tag == TAG_OUTPUT_B){

    // fprintf(stdout,"Bob Output!\n");
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

    // fprintf(stdout,"Alice Output!\n");
    evaluate_Gate(current_gate,current_key);
    evaluate_Alice_Output(current_gate,current_key);
    m_gate_index++;
    return &current_key;

  } else if (current_gate->tag == TAG_OUTPUT_B){

    // fprintf(stdout,"Bob Output!\n");
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

  current_key = _mm_xor_si128(key1,key2);
  
}

uint32_t GarbledCircuit::increment_index(){
  return m_gate_index++;
}


void GarbledCircuit::generate_Bob_Input(PCFGate* current_gate, __m128i &current_key){
  // Gen's input keys have already been generated and determined
  // by his input keys
  // here, we return the proper input key encoding 0
  // to do this, we need to know both his permutation bits and his input bits
  
  // fprintf(stdout,"\nBob/Gen Input Gate\n");
  
  uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index
  
  Bytes gen_input = get_Gen_Input(2*gen_input_idx + get_Input_Parity(gen_input_idx));

  save_Key_to_128bit(gen_input,current_key);
  
  // we need to set the garbling buffer to something (or nothing)
  // because it will be sent by default between gates
  m_garbling_bufr = Bytes(0);
}

void GarbledCircuit::evaluate_Bob_Input(PCFGate* current_gate, __m128i &current_key){
  // here, Gen's input keys have been sent to Eval
  // she needs only to assign the input key to the wire
  
  uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index
  Bytes gen_input = get_Gen_Input(gen_input_idx);

  save_Key_to_128bit(gen_input,current_key);
  
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

    
    //fprintf(stdout, "\nAlice/Evl Input:: new current key:\n");
    //print128_num(current_key);
    // __m128i xor_key = _mm_xor_si128(current_key, m_R);
    //fprintf(stdout, "Bob xor-input key:\n");
    //print128_num(xor_key);
    

    // load the first output key
    tmp = get_Evl_Key(current_gate->wire1,0);
    tmp.resize(16,0);
    output_keys[0] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
    // quick print check
    // fprintf(stdout,"key1: %s\n",tmp.to_hex().c_str());

    
    // load the second output key
    tmp = get_Evl_Key(current_gate->wire1,1); // wire1 = wire2
    tmp.resize(16,0);
    output_keys[1] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
    // quick print check
    // fprintf(stdout,"key2: %s\n",tmp.to_hex().c_str());
    
    
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
  
  
  //fprintf(stdout,"Evaluate Alice Input\n");
  //fprintf(stdout,"\nreceived garbling buffer: %s\n",m_garbling_bufr.to_hex().c_str());
  //fprintf(stdout,"selecting and decrypting Eval Input Key: %x (input wire %i)\n",bit,current_gate->wire1);
  
  
  // now select the one of two ciphertexts to decrypt
  Bytes encrypted_input;
  encrypted_input.insert(encrypted_input.end(),
                         m_garbling_bufr.begin() + bit * Env::key_size_in_bytes(),
                         m_garbling_bufr.begin() + Env::key_size_in_bytes() + bit * Env::key_size_in_bytes());
  

  assert(evl_input.size() == Env::key_size_in_bytes());
  assert(evl_input.size() == encrypted_input.size());
  
  evl_input = evl_input ^ encrypted_input;
  save_Key_to_128bit(evl_input,current_key);

}

void GarbledCircuit::evaluate_Alice_Output(PCFGate* current_gate, __m128i &current_key){
    // make sure the output buffer is big enough
    if(m_alice_out.size()*8 <= m_alice_out_ix){
      m_alice_out.resize((m_alice_out.size()+1)*2,0); // grow by doubling for less work
    }

    uint8_t out_bit = _mm_extract_epi8(current_key, 0) & 0x01;
    m_alice_out.set_ith_bit(m_alice_out_ix, out_bit);
    m_alice_out_ix++;

}
void GarbledCircuit::evaluate_Bob_Output(PCFGate* current_gate, __m128i &current_key){
    if (m_bob_out.size()*8 <= m_bob_out_ix)
      {
        // dynamically grow output array by doubling
        m_bob_out.resize((m_bob_out.size()+1)*2, 0);
      }
    
    uint8_t out_bit = _mm_extract_epi8(current_key, 0) & 0x01;
    m_bob_out.set_ith_bit(m_bob_out_ix, out_bit);
    m_bob_out_ix++;

}

void GarbledCircuit::generate_Alice_Output(PCFGate* current_gate, __m128i &current_key){
    
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
  if (m_bob_out.size()*8 <= m_bob_out_ix)
    {
     // dynamically grow output array by doubling
     m_bob_out.resize((m_bob_out.size()+1)*2, 0);
   }
  
  uint8_t out_bit = _mm_extract_epi8(current_key, 0) & 0x01;
  
  m_bob_out.set_ith_bit(m_bob_out_ix, out_bit);
  m_bob_out_ix++;
  
}


void GarbledCircuit::generate_Gate(PCFGate* current_gate, __m128i &current_key){
  Bytes tmp;
  
#ifdef FREE_XOR
  if(current_gate->truth_table == 0x06){ // if XOR gate
    xor_Gate(current_gate, current_key);
  } else {
#endif
    if(current_gate->truth_table == 0x01){ // AND Gate
      
      //fprintf(stdout,"AND GATE!\n");
      __m128i key1 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1));
      __m128i key2 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2));
      
      genHalfGatePair(current_key, key1, key2, m_garbling_bufr, 0, 0, 0);
    }
    
    else if(current_gate->truth_table == 0x07){ // OR Gate

      // fprintf(stdout,"OR Gate!");
      __m128i key1 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1));
      __m128i key2 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2));
      genHalfGatePair(current_key, key1, key2, m_garbling_bufr, 1, 1, 1);

    }
    
    else { 

      // here (most likely) we have a NOT gate or an XNOR gate 
      // the compiler's optimizer should do its best to remove them
      // but since they can't be garbled with half gates, we garble with GRR
      // TODO: implement fixed-key garbling for two keys as in JustGarble

      // fprintf(stdout,"Gate Type: %x\n",current_gate->truth_table);
      
      // first, get the inputs
      // X and Y are input, Z is output
      __m128i X[2], Y[2], Z[2];
      // we need a couple of 128-bit variables
      __m128i aes_plaintext, aes_ciphertext;
      __m128i aes_key[2];
    
      uint8_t bit;
      
      // load the (zero-key) inputs from the PCF state container
      X[0] = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1));
      Y[0] = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2));
      // and derive their XOR-complements
      X[1] = _mm_xor_si128(X[0], m_R); // X[1] = X[0] ^ R
      Y[1] = _mm_xor_si128(Y[0], m_R); // Y[1] = Y[0] ^ R
      
      // and get the permutation bits (tells if each zero-key has a 1 on the end)
      const uint8_t perm_x = _mm_extract_epi8(X[0],0) & 0x01;
      const uint8_t perm_y = _mm_extract_epi8(Y[0],0) & 0x01;
      // construct the permutation index
      const uint8_t de_garbled_ix = (perm_y << 1)|perm_x; // wire1 + 2*wire2
      
      // now we are ready to garble
      // first, load the plaintext to be encrypted
      aes_plaintext = _mm_set1_epi64x(m_gate_index);
      // and then load the garbling keys      
      aes_key[0] = _mm_load_si128(X+perm_x);
      aes_key[1] = _mm_load_si128(Y+perm_y);
      
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
      // and load it into current_key to return it to the circuit's state container
      current_key = _mm_load_si128(Z);
      // this way we return the new zero-key
      
      
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
      // and add tmp to the buffer we're sending to Eval
      m_garbling_bufr.insert(m_garbling_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
      
      
      // encrypt the 2nd entry : (X[x], Y[1-y])
      aes_key[0] = _mm_xor_si128(aes_key[0], m_R);
      aes_key[1] = _mm_xor_si128(aes_key[1], m_R);
    
      KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
      aes_ciphertext = _mm_and_si128(aes_ciphertext, m_clear_mask);
      bit = (current_gate->truth_table>>(3-(0x02^de_garbled_ix)))&0x01;
      aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
      
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);
      
      m_garbling_bufr.insert(m_garbling_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
      
      // encrypt the 3rd entry : (X[1-x], Y[1-y])
      aes_key[0] = _mm_xor_si128(aes_key[0], m_R);
      
      KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
      aes_ciphertext = _mm_and_si128(aes_ciphertext, m_clear_mask);
      bit = (current_gate->truth_table>>(3-(0x03^de_garbled_ix)))&0x01;
      aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);
      
      m_garbling_bufr.insert(m_garbling_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
      
    }
    
#ifdef FREE_XOR
  }
#endif
  
}


void GarbledCircuit::evaluate_Gate(PCFGate* current_gate, __m128i &current_key){
   
#ifdef FREE_XOR
  if (current_gate->truth_table == 0x06)
    {
      xor_Gate(current_gate, current_key);
    } else {
#endif    
    if(current_gate->truth_table == 0x01){ // AND Gate
      // fprintf(stdout,"AND GATE!\n");
      __m128i key1 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1));
      __m128i key2 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2));
      
      evlHalfGatePair(current_key, key1,key2,m_garbling_bufr,0);

    } 
    else if(current_gate->truth_table == 0x07){ // OR gate
      // fprintf(stdout,"OR GATE!\n");
      __m128i key1 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1));
      __m128i key2 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2));
      
      evlHalfGatePair(current_key, key1,key2,m_garbling_bufr,1);
      
    }
    else { 

      // here (most likely) we have a NOT gate or an XNOR gate 
      // the compiler's optimizer should do its best to remove them
      // but since they can't be garbled with half gates, we garble with GRR
      // TODO: implement fixed-key garbling for two keys as in JustGarble

    __m128i aes_key[2], aes_plaintext, aes_ciphertext;
    Bytes tmp;
    __m128i a;
    Bytes::const_iterator it;
    
    aes_plaintext = _mm_set1_epi64x(m_gate_index);
    
    aes_key[0] = *reinterpret_cast<__m128i*>(get_wire_key(m_st,current_gate->wire1));
    aes_key[1] = *reinterpret_cast<__m128i*>(get_wire_key(m_st,current_gate->wire2));
    
    const uint8_t perm_x = _mm_extract_epi8(aes_key[0], 0) & 0x01;
    const uint8_t perm_y = _mm_extract_epi8(aes_key[1], 0) & 0x01;
    
    //fprintf(stdout,"garbling input keys: %i %i (parity %i %i)\n",current_gate->wire1, current_gate->wire2,perm_x, perm_x);
    //print128_num(aes_key[0]);
    //print128_num(aes_key[1]);
    

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
    
    }
#ifdef FREE_XOR
  }
#endif

  // current_key will be returned
}



 
/**
   This function doubles in GF2 by shifting "key" left one bit.
   the clear mask is included to ensure that all keys remain the same length
   (we don't want a key boundary overflowing!)
 */
void Double(__m128i & key, __m128i & clear_mask){
  
  uint64_t * modifier = (uint64_t*) &key;
  uint64_t carry = modifier[0] & ((uint64_t) 0x8000000000000000) > 0 ? 1 : 0;
  modifier[0] = modifier[0]<<1;
  modifier[1] = (modifier[1]<<1) | carry;
  key = _mm_and_si128(key, clear_mask);

}


/**
   This function computes the permutation on an input key
   it is destructive of key so must make copies of the inputs first
   it returns H(K) = pi(L) xor L where L = 2key ^ tweak
 */
void H_Pi(__m128i & destination, __m128i &key, __m128i & tweak, __m128i & clear_mask){
  __m128i K,K1;
  
  Double(key, clear_mask);
  K = _mm_xor_si128(key,tweak);
  K1 = K;
  KDF128((uint8_t*)&K,(uint8_t*)&destination, (uint8_t*)&K1);
  destination = _mm_xor_si128(destination, K1);
  destination = _mm_and_si128(destination,clear_mask);
  
}


void GarbledCircuit::genHalfGatePair(__m128i& out_key, __m128i & key1, __m128i & key2, Bytes & out_bufr, byte a1, byte a2, byte a3){
  // this function implements half-gate generation by
  // Zahur, Rosulek, and Evans

  assert(a1==0||a1==1);
  assert(a2==0||a2==1);
  assert(a3==0||a3==1);
  
  const uint8_t perm_x = _mm_extract_epi8(key1,0) & 0x01;
  const uint8_t perm_y = _mm_extract_epi8(key2,0) & 0x01;

  __m128i Wa0, Wb0,Wa1,Wb1; // incoming wire keys
  __m128i Wg, We; // half gate wire keys

  // Tg and Te are the transmitted variables for the Generator and Evaluator, respectively 
  __m128i Tg, Te;
  // Ha and Hb will contain the results of applying the KDF
  __m128i Ha0, Ha1, Hb0, Hb1;

  // __m128i tmp; // temporarily stores the output of H that will be used again
  Bytes tmp_bufr; // transfers keys to output buffer
  
  uint32_t j1, j2;
  j1 = increment_index();
  j2 = increment_index();

  __m128i j1_128, j2_128;
  j1_128 = _mm_set1_epi64x(j1);
  j2_128 = _mm_set1_epi64x(j2);
  
  Wa0 = key1;
  Wb0 = key2;

  __m128i key1x = _mm_xor_si128(key1,m_R);

  Wa1 = _mm_xor_si128(Wa0,m_R);
  Wb1 = _mm_xor_si128(Wb0,m_R);

  __m128i H_Wa0j1, H_Wa1j1, H_Wb0j2, H_Wb1j2;
  H_Pi(H_Wa0j1, Wa0, j1_128, m_clear_mask);
  H_Pi(H_Wa1j1, Wa1, j1_128, m_clear_mask);
  H_Pi(H_Wb0j2, Wb0, j2_128, m_clear_mask);
  H_Pi(H_Wb1j2, Wb1, j2_128, m_clear_mask);
  

  // first half gate
  Tg = _mm_xor_si128(H_Wa0j1, H_Wa1j1);
  if(perm_y != a2){
    Tg = _mm_xor_si128(Tg,m_R);
  }

  Wg = perm_x? H_Wa1j1 : H_Wa0j1;
  if(((perm_x != a1) && (perm_y != a2)) != a3){
    Wg = _mm_xor_si128(Wg, m_R);
  }
  
  // second half gate
  We = perm_y ? H_Wb1j2 : H_Wb0j2;
  Te = _mm_xor_si128(H_Wb0j2,H_Wb1j2);
  Te = _mm_xor_si128(Te, a1? key1x :key1);

  // add Tg to output buffer
  tmp_bufr.resize(16,0);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp_bufr[0]),Tg);
  out_bufr.clear();
  out_bufr.insert(out_bufr.begin(),tmp_bufr.begin(),tmp_bufr.begin()+Env::key_size_in_bytes());

  // add Te to output buffer
  tmp_bufr.clear();
  tmp_bufr.resize(16,0);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp_bufr[0]),Te);
  out_bufr.insert(out_bufr.end(),tmp_bufr.begin(),tmp_bufr.begin()+Env::key_size_in_bytes());

  // combine half gates:
  out_key = _mm_xor_si128(Wg,We);

}
  


void GarbledCircuit::evlHalfGatePair(__m128i &current_key, __m128i & key1, __m128i & key2, Bytes & in_bufr, byte a1){
  assert(in_bufr.size() == 2*Env::key_size_in_bytes());


  // get the select bits

  byte sa,sb;
  sa = _mm_extract_epi8(key1,0) & 0x01;
  sb = _mm_extract_epi8(key2,0) & 0x01;

  // get the counter values
  uint32_t j1, j2;
  
  j1 = increment_index();
  j2 = increment_index();

  // fprintf(stdout,"Gate indices: %x %x\n", j1,j2);

  __m128i j1_128, j2_128;

  j1_128 = _mm_set1_epi64x(j1);
  j2_128 = _mm_set1_epi64x(j2);

  // fprintf(stdout,"Half Gate In Buffer: %s\n",in_bufr.to_hex().c_str());

  Bytes tmp_bufr;
  __m128i Tg, Te;
  
  // TG is always sent first, and then TE
  // where TG is the single row transmitted for the Generator's half-gate
  // and TE is the single row transmitted for the Evaluator's half-gate
  tmp_bufr.resize(16,0);
  tmp_bufr.insert(tmp_bufr.begin(),in_bufr.begin(), in_bufr.begin()+Env::key_size_in_bytes());
  Tg = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp_bufr[0]));
  
  tmp_bufr.clear();
  tmp_bufr.resize(16,0);
  tmp_bufr.insert(tmp_bufr.begin(),in_bufr.begin()+Env::key_size_in_bytes(),in_bufr.begin()+2*Env::key_size_in_bytes());
  Te = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp_bufr[0]));

  __m128i Wa,Wb;

  Wa = key1;
  Wb = key2;

  __m128i H_Waj1,H_Wbj2; // output of hashes
  H_Pi(H_Waj1, Wa, j1_128, m_clear_mask);
  H_Pi(H_Wbj2, Wb, j2_128, m_clear_mask);

  __m128i Wg, We,tmp,tmpwe;
  __m128i xorTeKey1;

  Wg = H_Waj1;
  if(sa){
    Wg = _mm_xor_si128(H_Waj1,Tg);
  }
  
  We = H_Wbj2;
  xorTeKey1 = _mm_xor_si128(Te, key1);
  if(sb){
    We = _mm_xor_si128(H_Wbj2,xorTeKey1);
  }

  current_key = _mm_xor_si128(Wg,We);

}  


void GarbledCircuit::set_const_key(byte c, const Bytes &key)
{
  assert(c == 0 || c == 1); // wire for constant 0 or 1
  Bytes tmp = key;
  tmp.resize(16,0);
  m_const_wire[c] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
}


Bytes GarbledCircuit::get_garbling_bufr(){
  return m_garbling_bufr;
}

void GarbledCircuit::set_garbling_bufr(Bytes buf){
  m_garbling_bufr = buf;
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
