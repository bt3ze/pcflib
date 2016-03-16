#ifndef GARBLED_CIRCUIT_CPP
#define GARBLED_CIRCUIT_CPP

#include "GarbledCircuit.h"
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <algorithm> 

/**
   ACCESSORY FUNCTIONS
 */

/*
void *copy_key(void *old_key)
{
  __m128i *new_key = 0;
  if (old_key != 0)
    {
      // first argument is size, second argument is allignment
      new_key = (__m128i*)_mm_malloc(sizeof(__m128i), sizeof(__m128i));
      *new_key = *reinterpret_cast<__m128i*>(old_key);
    }

  //print128_num(*new_key);
  return new_key;
}
*/

void copy_key(void* source_key, void * dest_key){
  //  __m128i *new_key = 0; 
 
  //fprintf(stdout,"copy key \n");

  if (source_key != 0)
    {
      // first argument is size, second argument is allignment
      //new_key = (__m128i*)_mm_malloc(sizeof(__m128i), sizeof(__m128i));

      num_copies++;
      clock_gettime(CLOCK_REALTIME, &copy_start);
      
      _mm_storeu_si128(reinterpret_cast<__m128i*>(dest_key),*reinterpret_cast<__m128i*>(source_key));
     
      clock_gettime(CLOCK_REALTIME, &copy_end);
      copy_time += ( copy_end.tv_sec - copy_start.tv_sec )
        + ( copy_end.tv_nsec - copy_start.tv_nsec )
        / BILN;

    } else{
    fprintf(stderr,"no copy\n");
  }
  
  //  return new_key;
 
}


void delete_key(void *key)
{
  // if (key != 0) _mm_free(key);
}


void save_Key_to_128bit(const Bytes & key, __m128i & destination){
  num_buffers++;
  clock_gettime(CLOCK_REALTIME, &buffer_start);

  Bytes tmp = key;

  tmp.resize(16,0);
  destination = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));

  //destination = _mm_loadu_si128(reinterpret_cast<__m128i*>(key[0]));

  clock_gettime(CLOCK_REALTIME, &buffer_end);
  buffer_time += ( buffer_end.tv_sec - buffer_start.tv_sec )
    + ( buffer_end.tv_nsec - buffer_start.tv_nsec )
    / BILN;
}

void append_m128i_to_Bytes(const __m128i & num, Bytes & dest){
  Bytes tmp;
  
  num_buffers++;
  clock_gettime(CLOCK_REALTIME, &buffer_start);
  

  tmp.resize(16,0);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]),num);
  dest.insert(dest.end(),tmp.begin(),tmp.begin()+Env::key_size_in_bytes());

  //m_temp_bufr.clear();
  //_mm_storeu_si128(reinterpret_cast<__m128i*>(&m_temp_bufr[0]),num);
  //dest.insert(dest.end(),m_temp_bufr.begin(),m_temp_bufr.begin()+Env::key_size_in_bytes());
  
  //  _mm_storeu_si128(reinterpret_cast<__m128i*>(&dest[0]),num);
  

  clock_gettime(CLOCK_REALTIME, &buffer_end);
  buffer_time += ( buffer_end.tv_sec - buffer_start.tv_sec )
    + ( buffer_end.tv_nsec - buffer_start.tv_nsec )
    / BILN;

}


void insert_to_garbling_bufr(const __m128i & num, Bytes & dest, Bytes & tmp_bufr, uint32_t pos){
  
  //std::fill(tmp_bufr.begin(),tmp_bufr.end(),0);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp_bufr[0]),num);
  dest.insert(dest.begin()+pos*Env::key_size_in_bytes(),tmp_bufr.begin(),tmp_bufr.begin()+Env::key_size_in_bytes());

  //_mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp_bufr[0]),num);
 

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

  //fprintf(stdout,"get next gate\n");
  GarbledCircuit &cct = *reinterpret_cast<GarbledCircuit*>(get_external_circuit(st));

  return cct.gen_Next_Gate(current_gate);
 }

void * evl_next_gate(PCFState *st, PCFGate *current_gate){
  // returns void pointer, which is pointer to a key
  // use this one to call the Garbled Circuit object again

  //fprintf(stdout,"evl next gate\n");
  GarbledCircuit &cct = *reinterpret_cast<GarbledCircuit*>(get_external_circuit(st));
  
  // now, call the appropriate function from cct
  return cct.evl_Next_Gate(current_gate);
 }

void print128_num(__m128i var)
{
  uint16_t *val = (uint16_t*) &var;
  fprintf(stderr,"Numerical: %x %x %x %x %x %x %x %x \n", 
          val[0], val[1], val[2], val[3], val[4], val[5], 
          val[6], val[7]);
}




/**
   INITIALIZING GARBLED CIRCUITS
 */

  
GarbledCircuit::GarbledCircuit(): m_gate_index(0), m_bob_out_ix(0),m_alice_out_ix(0),m_hash_row_idx(0), m_comm_time(0.0) {

  // initialize the key Mask
  Bytes tmp(16);
  for(size_t ix = 0; ix< Env::k(); ix++){
    tmp.set_ith_bit(ix,1);
  }
  m_clear_mask = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));

  m_alice_out.clear();
  m_bob_out.clear();

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


void GarbledCircuit::init_Generation_Circuit(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys, const uint32_t gen_inp_size, Bytes & circuit_seed, const Bytes & perm_bits, const Bytes R, const Bytes & zero_key, const Bytes & one_key){
  
  Bytes tmp;
  
  m_prng.seed_rand(circuit_seed);

  m_select_bits.insert(m_select_bits.begin(),perm_bits.begin(),perm_bits.begin() + perm_bits.size());

  // initialize the value of R
  // we initialize tmp with zeros for this so that 
  // the value of R doesn't get extra nonsense at the end
  // from whatever was sitting around in memory
  tmp.resize(16,0);
  tmp.insert(tmp.begin(),R.begin(),R.begin()+Env::k()/8);
  m_R = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&tmp[0]));
  
  // also initialize the Bytes version of R,
  // used for returning keys of proper parity for Gen
  m_R_bytes = R;

  // initialize the constant wires
  save_Key_to_128bit(zero_key,m_const_wire[0]);
  save_Key_to_128bit(one_key, m_const_wire[1]);
  

  m_gen_inp_size = gen_inp_size;
  set_Input_Keys(gen_keys, evl_keys);

  // for now, use zeroes and the aes key
  __m128i aes_key = m_const_wire[0];
  //aes_key = _mm_xor_si128(aes_key,aes_key);
  init_circuit_AES_key(aes_key);


  m_garbling_bufr.resize(256,0);
  m_temp_bufr.resize(256,0);
  
  m_message_queue.resize(MESSAGE_LIMIT*Env::key_size_in_bytes());
  m_message_limit = MESSAGE_LIMIT;
  m_messages_waiting = 0;
  m_queue_index = 0;


  xor_gates = half_gates = other_gates = total_gates = 0;
  xor_time = hg_time = og_time = garble_time = 0.0;
  
  

}

void GarbledCircuit::init_Evaluation_Circuit(const std::vector<Bytes> * gen_keys, const std::vector<Bytes> * evl_keys,const uint32_t gen_inp_size, const Bytes & evl_input, const Bytes &zero_key, const Bytes & one_key){
  set_Input_Keys(gen_keys, evl_keys);

  //m_select_bits = evl_input;
  m_select_bits.insert(m_select_bits.begin(),evl_input.begin(),evl_input.begin() + evl_input.size());

  save_Key_to_128bit(zero_key,m_const_wire[0]);
  save_Key_to_128bit(one_key, m_const_wire[1]);

  m_gen_inp_size = gen_inp_size;

  // for now, let the fixed aes key be zeroes.
  __m128i aes_key = m_const_wire[0];
  //  aes_key = _mm_xor_si128(aes_key,aes_key);
  init_circuit_AES_key(aes_key);

  m_garbling_bufr.resize(256,0);
  m_temp_bufr.resize(16,0);
   
  m_message_queue.resize(MESSAGE_LIMIT*Env::key_size_in_bytes());
  m_message_limit = MESSAGE_LIMIT;
  m_messages_waiting = 0;
  m_queue_index = 0;
  
  
  xor_gates = half_gates = other_gates = total_gates = 0;
  xor_time = hg_time = og_time = garble_time = 0.0;

  
}

/**
   this function sets our AES key for constant-key garbling
   the key must be fed to the function
 */
void GarbledCircuit::init_circuit_AES_key(__m128i &key){

  AES_set_encrypt_key((unsigned char*)&key, 128, &m_fixed_key);  

}


void * GarbledCircuit::get_Const_Wire(uint32_t i){
  assert(i == 0 || i == 1);
  return &(m_const_wire[i]);
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
  Bytes ret = (*m_gen_inputs)[idx];
  return ret;
}

Bytes GarbledCircuit::get_Evl_Input(uint32_t idx){
  if(idx < m_evl_inputs->size()){
    return (*m_evl_inputs)[idx];
  } else{
    fprintf(stdout,"Eval uninitialized wire: %i\n",idx);
    Bytes tmp;
    __m128i key = *((__m128i*)get_Const_Wire(0));
    append_m128i_to_Bytes(key, tmp );
    return tmp;
  }
}

uint32_t GarbledCircuit::get_Input_Parity(uint32_t idx){
  // Gen uses this to look up permutation bits
  // Eval uses this function to figure out what input key to decrpyt
  if(idx < m_select_bits.size()*8){
    return m_select_bits.get_ith_bit(idx);
  } else{
    fprintf(stdout,"input out of bounds: %x\n",idx);
    return 0;
  }
}

/**
   EVALUATION/GARBLING FUNCTION
 */
void GarbledCircuit::Garble_Circuit(){
  while(get_next_gate(m_st)){
  }
  //send_buffer(); // last send, in case there's something left over
}

void GarbledCircuit::Evaluate_Circuit(){
  do {
    
  } while(get_next_gate(m_st));
  //retrieve_buffer();
  fprintf(stderr,"finish evaluating\n");
}

/**
   GARBLING ACCESSORY FUNCTIONS
 */

 
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
void H_Pi(__m128i & destination, __m128i &key, __m128i & tweak, __m128i & clear_mask, AES_KEY_J & fixed_key){
  __m128i K; // ,K1;

  Double(key, clear_mask);
  K = _mm_xor_si128(key,tweak);

  KDF128((uint8_t*)&destination,(uint8_t*)&K, &fixed_key);
  destination = _mm_xor_si128(destination, K);
  destination = _mm_and_si128(destination,clear_mask);  
}


/**
   remember to copy the keys before they enter this function, because it's destructive to key1 and key2
 */
void H_Pi256(__m128i & destination, __m128i &key1, __m128i &key2, __m128i & tweak, __m128i & clear_mask, AES_KEY_J & fixed_key){
  // takes two keys and computes a ciphertest, A la JustGarble
  __m128i K;
  
  Double(key1,clear_mask);
  Double(key2,clear_mask);
  Double(key2,clear_mask);
  
  K = _mm_xor_si128(key1,key2);
  K = _mm_xor_si128(K,tweak);

  KDF128((uint8_t*)&destination,(uint8_t*)&K, &fixed_key);
  destination = _mm_xor_si128(destination, K);
  destination = _mm_and_si128(destination,clear_mask);

}



/**
  throughout garbling, Gen will maintain the zero-semantic keys for each wire
  and use them to generate all of the ciphertexts (find the 1-semantics by XOR with R)
  that are sent to Eval

  this function returns current_key to the pcf_state struct
 */
void * GarbledCircuit::gen_Next_Gate(PCFGate *current_gate){
  
  //clock_gettime(CLOCK_REALTIME, &bstart);
    

  // double start = MPI_Wtime();

  static __m128i current_key; // must be static to return it

  if(current_gate->tag == TAG_INTERNAL){
    // actual gate
    clear_garbling_bufr();
    generate_Gate(current_gate,current_key,m_garbling_bufr);
    m_gate_index++;
  
  } else if(current_gate->tag == TAG_INPUT_A){

    // fprintf(stdout, "Alice Input!");
    clear_garbling_bufr();
    generate_Alice_Input(current_gate, current_key, m_garbling_bufr);

    // send two ciphertexts

    // return &current_key;
    
  }  else if (current_gate->tag == TAG_INPUT_B){

    //  fprintf(stdout, "Bob Input!");
    generate_Bob_Input(current_gate, current_key);
    // return &current_key;

  } else if (current_gate->tag == TAG_OUTPUT_A) {
    
    // fprintf(stdout,"Alice Output!\n");
    clear_garbling_bufr();
    generate_Alice_Output(current_gate,current_key, m_garbling_bufr);
    m_gate_index++;
    // return &current_key;

  } else if (current_gate->tag == TAG_OUTPUT_B){

    // fprintf(stdout,"Bob Output!\n");
    
    clear_garbling_bufr();
    generate_Bob_Output(current_gate, current_key, m_garbling_bufr);
    
    Bytes tmp;
    append_m128i_to_Bytes(current_key,tmp);
    m_gen_output_labels.push_back(tmp); // for output authenticity proof
    
    m_gate_index++;
    //return &current_key;

  } else {

    fprintf(stdout,"read gate error\n");
    exit(0);
  
    //return &current_key; 
  }

  // benchmark_time +=  MPI_Wtime() - start;
  //  clock_gettime(CLOCK_REALTIME, &bend);
  //btime2 += ( bend.tv_sec - bstart.tv_sec )
  //  + ( bend.tv_nsec - bstart.tv_nsec )
  //  / BILN;

  return &current_key;
}


void * GarbledCircuit::evl_Next_Gate(PCFGate *current_gate){
  
  // double start = MPI_Wtime();

  //clock_gettime(CLOCK_REALTIME, &bstart);
   
  static __m128i current_key; // must be static to return it
  
  if(current_gate->tag == TAG_INTERNAL){
    evaluate_Gate(current_gate,current_key, m_garbling_bufr);
    m_gate_index++;

  } else if(current_gate->tag == TAG_INPUT_A){
    
    evaluate_Alice_Input(current_gate, current_key, m_garbling_bufr);
    //return &current_key;
    
  } else if (current_gate->tag == TAG_INPUT_B){
    
    evaluate_Bob_Input(current_gate,current_key);
    //return &current_key; 
    
  } else if (current_gate->tag == TAG_OUTPUT_A) {

    // fprintf(stdout,"Alice Output!\n");
    
    evaluate_Alice_Output(current_gate,current_key, m_garbling_bufr);

    m_gate_index++;
    //return &current_key;

  } else if (current_gate->tag == TAG_OUTPUT_B){

    // fprintf(stdout,"Bob Output!\n");

    //mask Bob's output with his output masking key
    
    evaluate_Bob_Output(current_gate, current_key, m_garbling_bufr);
 
    // save this for later, when we need to do gen's
    // output authenticity proof
    Bytes tmp;
    append_m128i_to_Bytes(current_key,tmp);
    m_gen_output_labels.push_back(tmp);

    m_gate_index++;
    //return &current_key;
    
} else {
    
    fprintf(stderr,"gate read error\n");
    exit(0);
    //return &current_key;
  }

  //  benchmark_time +=  MPI_Wtime() - start;
 

  //clock_gettime(CLOCK_REALTIME, &bend);
  //btime2 += ( bend.tv_sec - bstart.tv_sec )
  //  + ( bend.tv_nsec - bstart.tv_nsec )
  //  / BILN;

 return &current_key;

}


void GarbledCircuit::xor_Gate(__m128i & key1, __m128i & key2, __m128i &current_key){

  current_key = _mm_xor_si128(key1,key2);



}

uint32_t GarbledCircuit::increment_index(){
  return m_gate_index++;
}


void GarbledCircuit::generate_Bob_Input(PCFGate* current_gate, __m128i &current_key){
  // Gen's input keys have already been generated and determined
  // by his input keys; here, we return the proper input key encoding 0
  
  uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index 
 
  Bytes gen_input = get_Gen_Key(gen_input_idx,get_Input_Parity(gen_input_idx));
  
  save_Key_to_128bit(gen_input,current_key);
  
  // we need to set the garbling buffer to something (or nothing)
  // because it will be sent by default between gates
  // m_garbling_bufr = Bytes(0);

}

void GarbledCircuit::evaluate_Bob_Input(PCFGate* current_gate, __m128i &current_key){
  // here, Gen's input keys have been sent to Eval
  // she needs only to assign the input key to the wire
  
  uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index
  Bytes gen_input = get_Gen_Input(gen_input_idx);

  save_Key_to_128bit(gen_input,current_key);
  
  // and current_key is available up the call stack
}

void GarbledCircuit::generate_Alice_Input(PCFGate* current_gate, __m128i &current_key, Bytes & garbling_bufr){
    // Here, Eval already has her input keys, but Gen doesn't know
    // which to use. So he will generate a new key and its XOR-offset,
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


    // load the first output key
    tmp = get_Evl_Key(current_gate->wire1,0);
    tmp.resize(16,0);
    output_keys[0] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
    
    // load the second output key
    tmp = get_Evl_Key(current_gate->wire1,1); // wire1 = wire2
    tmp.resize(16,0);
    output_keys[1] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
    

    // encrypt the output keys so that Eval can get only one of them
    output_keys[0] = _mm_xor_si128(output_keys[0], current_key);
    output_keys[1] = _mm_xor_si128(output_keys[1], _mm_xor_si128(current_key,m_R));

    
    // send the output keys to Eval for decryption    
    append_m128i_to_Bytes(output_keys[0],garbling_bufr);
    append_m128i_to_Bytes(output_keys[1],garbling_bufr);

    assert(garbling_bufr.size() == 2*Env::key_size_in_bytes());
    send_half_gate(garbling_bufr);

}

void GarbledCircuit::evaluate_Alice_Input(PCFGate* current_gate, __m128i &current_key, Bytes & garbling_bufr){
  // here, Eval already knows which input key she wants to use
  // she selects it and assigns it to her wire value
  
  uint32_t gen_input_idx = current_gate->wire1; // wire1 holds the input index
  Bytes evl_input = get_Evl_Input(gen_input_idx);
  
  // receive 2 ciphertexts
  //garbling_bufr = read_half_gate();
  read_half_gate(garbling_bufr);

  assert(garbling_bufr.size() == 2*Env::key_size_in_bytes());
  
  // find the input ciphertext that corresponds to Eval's chosen bit
  uint32_t bit = get_Input_Parity(current_gate->wire1);
  assert(bit == 0 || bit == 1);  
  

  // select the one to decrypt
  Bytes encrypted_input;
  encrypted_input.insert(encrypted_input.end(),
                         garbling_bufr.begin() + bit * Env::key_size_in_bytes(),
                         garbling_bufr.begin() + Env::key_size_in_bytes() + bit * Env::key_size_in_bytes());
  

  assert(evl_input.size() == Env::key_size_in_bytes());
  
  evl_input = evl_input ^ encrypted_input;
  save_Key_to_128bit(evl_input,current_key);

}

// Alice's outputs
void GarbledCircuit::evaluate_Alice_Output(PCFGate* current_gate, __m128i &current_key, Bytes & garbling_bufr){

  evaluate_Gate(current_gate,current_key,garbling_bufr);

  //print128_num(current_key);
  
  if(m_alice_out.size()*8 <= m_alice_out_ix){
    m_alice_out.resize((m_alice_out.size()+1)*2,0); // grow by doubling
  }
  
  uint8_t out_bit = _mm_extract_epi8(current_key, 0) & 0x01;
  m_alice_out.set_ith_bit(m_alice_out_ix, out_bit);
  m_alice_out_ix++;
    
}

// Gen's outputs
void GarbledCircuit::evaluate_Bob_Output(PCFGate* current_gate, __m128i &current_key, Bytes & garbling_bufr){

if (m_bob_out.size()*8 <= m_bob_out_ix)
    {
      // dynamically grow output array by doubling
      m_bob_out.resize((m_bob_out.size()+1)*2, 0);
    }
  
  Bytes output_mask = get_Gen_Input(m_gen_inp_size + m_bob_out_ix);
  __m128i output_mask_key;
  save_Key_to_128bit(output_mask, output_mask_key);
  
  current_key = _mm_xor_si128(output_mask_key, current_key);

  evaluate_Gate(current_gate, current_key, garbling_bufr);
  
  uint8_t out_bit = _mm_extract_epi8(current_key,0)& 0x01;
  
  m_bob_out.set_ith_bit(m_bob_out_ix, out_bit);
  m_bob_out_ix++;

}

void GarbledCircuit::generate_Alice_Output(PCFGate* current_gate, __m128i &current_key, Bytes & garbling_bufr){
  
  
  generate_Gate(current_gate, current_key, garbling_bufr); 
  
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


void GarbledCircuit::generate_Bob_Output(PCFGate* current_gate, __m128i &current_key, Bytes & garbling_bufr){ 

  if (m_bob_out.size()*8 <= m_bob_out_ix)
    {
      // dynamically grow output array by doubling
      m_bob_out.resize((m_bob_out.size()+1)*2, 0);
    }
  

  Bytes output_mask = get_Gen_Key(m_gen_inp_size + m_bob_out_ix,get_Input_Parity(m_gen_inp_size+m_bob_out_ix));
  
  __m128i output_mask_key;
  save_Key_to_128bit(output_mask, output_mask_key);

  current_key = _mm_xor_si128(output_mask_key, current_key);
  
  generate_Gate(current_gate,current_key, garbling_bufr); 
  
  uint8_t out_bit = (_mm_extract_epi8(current_key, 0) & 0x01); 
    // output_mask.get_ith_bit(0) ^ m_select_bits.get_ith_bit(m_gen_inp_size+m_bob_out_ix);
//get_Input_Parity(m_gen_inp_size+m_bob_out_ix); 
  
  m_bob_out.set_ith_bit(m_bob_out_ix, out_bit );
  m_bob_out_ix++;
  
}



void GarbledCircuit::generate_Gate(PCFGate* current_gate, __m128i &current_key, Bytes & garbling_bufr){

  __m128i key1 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1));
  __m128i key2 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2));

  total_gates++;
  clock_gettime(CLOCK_REALTIME, &garble_start);    

#ifdef FREE_XOR
  if(current_gate->truth_table == 0x06){ // if XOR gate
  
    xor_gates++;
    clock_gettime(CLOCK_REALTIME, &xor_start);    
   
    xor_Gate(key1, key2, current_key);

    clock_gettime(CLOCK_REALTIME, &xor_end);
    xor_time += ( xor_end.tv_sec - xor_start.tv_sec )
      + ( xor_end.tv_nsec - xor_start.tv_nsec )
      / BILN;


  } else if (current_gate->truth_table == 0x01 || current_gate->truth_table == 0x07)   {
  
    //fprintf(stderr,"%i\n",current_gate->truth_table);

    half_gates++;
    clock_gettime(CLOCK_REALTIME, &half_start); 

    if(current_gate->truth_table == 0x01){ // AND Gate            
      genHalfGatePair(current_key, key1, key2, m_temp_bufr, 0, 0, 0);
      send_half_gate(m_temp_bufr);
    }
    
    else if(current_gate->truth_table == 0x07){ // OR Gate
      genHalfGatePair(current_key, key1, key2, m_temp_bufr, 1, 1, 1);
      send_half_gate(m_temp_bufr);
    }
    
    
    clock_gettime(CLOCK_REALTIME, &half_end);
    hg_time += ( half_end.tv_sec - half_start.tv_sec )
      + ( half_end.tv_nsec - half_start.tv_nsec )
      / BILN;
    
  } else { 

#endif
      
    //fprintf(stderr,"%i\n",current_gate->truth_table);


    // here (most likely) we have a NOT gate or an XNOR gate 
    // the compiler's optimizer should do its best to remove them
    // but since they can't be garbled with half gates, we garble with GRR
    // NOT or XNOR gates, however, might be a bit cryptographically dangerous
    // we also use this method for output gates
    
    other_gates++;
    clock_gettime(CLOCK_REALTIME, &og_start);
    
    genStandardGate(current_key, key1, key2, garbling_bufr, current_gate->truth_table);
    send_full_gate(garbling_bufr);

    clock_gettime(CLOCK_REALTIME, &og_end);
    og_time += ( og_end.tv_sec - og_start.tv_sec )
      + ( og_end.tv_nsec - og_start.tv_nsec )
      / BILN;

#ifdef FREE_XOR
    }
#endif

    
  clock_gettime(CLOCK_REALTIME, &garble_end);
  garble_time += ( garble_end.tv_sec - garble_start.tv_sec )
    + ( garble_end.tv_nsec - garble_start.tv_nsec )
    / BILN;

}


void GarbledCircuit::evaluate_Gate(PCFGate* current_gate, __m128i &current_key, Bytes & garbling_bufr){
   
  __m128i key1 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1));
  __m128i key2 = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2));

  total_gates++;
  clock_gettime(CLOCK_REALTIME, &garble_start);    

#ifdef FREE_XOR
  if (current_gate->truth_table == 0x06)
    {
      xor_gates++;
      clock_gettime(CLOCK_REALTIME, &xor_start);    

      xor_Gate(key1, key2, current_key);

      clock_gettime(CLOCK_REALTIME, &xor_end);
      xor_time += ( xor_end.tv_sec - xor_start.tv_sec )
        + ( xor_end.tv_nsec - xor_start.tv_nsec )
        / BILN;
      
    } else if(current_gate->truth_table == 0x01 || current_gate->truth_table == 0x07) {
  
    half_gates++;
    clock_gettime(CLOCK_REALTIME, &half_start); 
   
    if(current_gate->truth_table == 0x01){ // AND Gate
      // fprintf(stdout,"AND GATE!\n");
      read_half_gate(garbling_bufr);
      //garbling_bufr = read_half_gate();
      evlHalfGatePair(current_key, key1,key2, garbling_bufr);
      
    } 
    else if(current_gate->truth_table == 0x07){ // OR gate
      // fprintf(stdout,"OR GATE!\n");
      //garbling_bufr = read_half_gate();
      read_half_gate(garbling_bufr);
      evlHalfGatePair(current_key, key1,key2, garbling_bufr);
      
    }

        
    clock_gettime(CLOCK_REALTIME, &half_end);
    hg_time += ( half_end.tv_sec - half_start.tv_sec )
      + ( half_end.tv_nsec - half_start.tv_nsec )
      / BILN;

  }else { 

#endif


    other_gates++;
    clock_gettime(CLOCK_REALTIME, &og_start);
    
    // here (most likely) we have a NOT gate or an XNOR gate 
      // the compiler's optimizer should do its best to remove them
      // but since they can't be garbled with half gates, we garble with GRR
      // we also use this for output gates
    
    read_full_gate(garbling_bufr);
    //garbling_bufr = read_full_gate();
    evlStandardGate(current_key, key1, key2, garbling_bufr);      
      
    clock_gettime(CLOCK_REALTIME, &og_end);
    og_time += ( og_end.tv_sec - og_start.tv_sec )
      + ( og_end.tv_nsec - og_start.tv_nsec )
      / BILN;
    
#ifdef FREE_XOR

    }
#endif

  // current_key will be available to the calling function

  
  clock_gettime(CLOCK_REALTIME, &garble_end);
  garble_time += ( garble_end.tv_sec - garble_start.tv_sec )
    + ( garble_end.tv_nsec - garble_start.tv_nsec )
    / BILN;

}

void GarbledCircuit::genStandardGate(__m128i& current_key, __m128i & key1, __m128i & key2, Bytes & out_bufr, uint8_t truth_table){

  // X and Y are input, Z is output
  __m128i X[2], Y[2], Z[2];
  __m128i garble_ciphertext;
  
  uint8_t semantic_bit;
    
  //double start = MPI_Wtime();

  // load the (zero-key) inputs from the PCF state container
  X[0] = key1;
  Y[0] = key2;
  // and XOR-complements
  X[1] = _mm_xor_si128(X[0], m_R); // X[1] = X[0] ^ R
  Y[1] = _mm_xor_si128(Y[0], m_R); // Y[1] = Y[0] ^ R
  

  // and get the permutation bits (tells if each zero-key has a 1 on the end)
  const uint8_t perm_x = _mm_extract_epi8(X[0],0) & 0x01;
  const uint8_t perm_y = _mm_extract_epi8(Y[0],0) & 0x01;
  const uint8_t de_garbled_ix = (perm_y << 1)|perm_x; // wire1 + 2*wire2
  
  
  // the last information for garbling
  uint32_t j1;
  j1 = increment_index();
  
  __m128i tweak;
  tweak = _mm_set1_epi64x(j1);


  // now run the key derivation function using the keys and the gate index  
  __m128i key1_in = _mm_loadu_si128(X+perm_x);
  __m128i key2_in = _mm_loadu_si128(Y+perm_y);

  H_Pi256(garble_ciphertext, key1_in, key2_in, tweak, m_clear_mask, m_fixed_key);
  semantic_bit = (truth_table >> (3-de_garbled_ix)) & 0x01;
  

  // GRR technique: using zero entry's key as one of the output keys
  // the output key is the encrypted gate index
  // this puts the zero key in Z[0] and the one key in Z[1]
  
  // we want the new zero-key to be in Z[0]
  // and other output is an offset
  _mm_store_si128(Z+semantic_bit, garble_ciphertext);
  Z[1 - semantic_bit] = _mm_xor_si128(Z[semantic_bit], m_R);
  current_key = _mm_loadu_si128(Z);

  
  // encrypt the first entry: (X[1-x],Y[y]) or (key1 xor R, key2)
  key1_in = _mm_loadu_si128(X+1-perm_x); 
  key2_in = _mm_loadu_si128(Y+perm_y);
  H_Pi256(garble_ciphertext, key1_in, key2_in, tweak, m_clear_mask, m_fixed_key);

  semantic_bit = (truth_table>>(3-(0x01^de_garbled_ix)))&0x01;
  garble_ciphertext = _mm_xor_si128(garble_ciphertext, Z[semantic_bit]);
  append_m128i_to_Bytes(garble_ciphertext,out_bufr);  
  

  // encrypt the 2nd entry : (X[x], Y[1-y])
  key1_in = _mm_loadu_si128(X+perm_x);
  key2_in = _mm_loadu_si128(Y+1-perm_y);
  H_Pi256(garble_ciphertext, key1_in, key2_in, tweak, m_clear_mask, m_fixed_key);

  
  semantic_bit = (truth_table>>(3-(0x02^de_garbled_ix)))&0x01;
  garble_ciphertext = _mm_xor_si128(garble_ciphertext, Z[semantic_bit]);
  append_m128i_to_Bytes(garble_ciphertext,out_bufr);
  
  
  // encrypt the 3rd entry : (X[1-x], Y[1-y])
  key1_in = _mm_loadu_si128(X+1-perm_x);
  key2_in = _mm_loadu_si128(Y+1-perm_y);
  H_Pi256(garble_ciphertext, key1_in, key2_in, tweak, m_clear_mask, m_fixed_key);
   

  semantic_bit = (truth_table>>(3-(0x03^de_garbled_ix)))&0x01;
  garble_ciphertext = _mm_xor_si128(garble_ciphertext, Z[semantic_bit]);
  append_m128i_to_Bytes(garble_ciphertext,out_bufr);
  
  
  // current_key holds our output key, and it will be available to our calling function
  // the calling function will also be able to send the information in out_bufr

  //benchmark_time +=  MPI_Wtime() - start;

}


void GarbledCircuit::evlStandardGate(__m128i& current_key, __m128i & key1, __m128i & key2, Bytes & in_bufr){
  __m128i garble_key[2], aes_plaintext, garble_ciphertext;
  Bytes tmp;
  __m128i a;
  Bytes::const_iterator it;
  
  //double start = MPI_Wtime();

  aes_plaintext = _mm_set1_epi64x(m_gate_index);
  
  garble_key[0] =  key1;
  garble_key[1] =  key2;
  
  const uint8_t perm_x = _mm_extract_epi8(garble_key[0], 0) & 0x01;
  const uint8_t perm_y = _mm_extract_epi8(garble_key[1], 0) & 0x01;
  
  
#ifndef AESNI
  KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&garble_ciphertext, (uint8_t*)garble_key);     
  garble_ciphertext = _mm_and_si128(garble_ciphertext, m_clear_mask);
#else
  uint32_t j1;
  j1 = increment_index();
  
  __m128i tweak;
  tweak = _mm_set1_epi64x(j1);
  
  //  __m128i key1_in = garble_key[0];
  // __m128i key2_in = garble_key[1];
  //H_Pi256(garble_ciphertext, key1_in, key2_in, tweak, m_clear_mask, m_fixed_key);
  H_Pi256(garble_ciphertext, garble_key[0], garble_key[1], tweak, m_clear_mask, m_fixed_key);

#endif
  
  uint8_t garbled_ix = (perm_y<<1)|perm_x;
  
#ifdef GRR
  if (garbled_ix == 0) {
   current_key = _mm_load_si128(&garble_ciphertext);
   }
    else
      {
        it = in_bufr.begin() + (garbled_ix-1)*Env::key_size_in_bytes();
        
        tmp.assign(it, it+Env::key_size_in_bytes());
        tmp.resize(16, 0);
        a = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
        current_key = _mm_xor_si128(garble_ciphertext, a);
      }
#else
    it = in_bufr.begin() + (garbled_ix)*Env::key_size_in_bytes();
    tmp.assign(it, it+Env::key_size_in_bytes());
    tmp.resize(16, 0);
    current_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
    current_key = _mm_xor_si128(current_key, garble_ciphertext);
#endif
    
    // current key holds our output key, and it will be available to the calling function

    //benchmark_time +=  MPI_Wtime() - start;

}


void GarbledCircuit::genHalfGatePair(__m128i& out_key, __m128i & key1, __m128i & key2, Bytes & out_bufr, byte a1, byte a2, byte a3){
  // this function implements half-gate generation by
  // Zahur, Rosulek, and Evans

  assert(a1==0||a1==1);
  assert(a2==0||a2==1);
  assert(a3==0||a3==1);
  
  //double start = MPI_Wtime();

  const uint8_t perm_x = _mm_extract_epi8(key1,0) & 0x01;
  const uint8_t perm_y = _mm_extract_epi8(key2,0) & 0x01;

  __m128i Wa0, Wb0,Wa1,Wb1; // incoming wire keys
  __m128i Wg, We; // half gate wire keys

  // Tg and Te are the transmitted variables for the Generator and Evaluator, respectively 
  __m128i Tg, Te;
  // Ha and Hb will contain the results of applying the KDF
  __m128i Ha0, Ha1, Hb0, Hb1;

  // __m128i tmp; // temporarily stores the output of H that will be used again
  //Bytes tmp_bufr; // transfers keys to output buffer
  
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
  H_Pi(H_Wa0j1, Wa0, j1_128, m_clear_mask,m_fixed_key);
  H_Pi(H_Wa1j1, Wa1, j1_128, m_clear_mask,m_fixed_key);
  H_Pi(H_Wb0j2, Wb0, j2_128, m_clear_mask,m_fixed_key);
  H_Pi(H_Wb1j2, Wb1, j2_128, m_clear_mask,m_fixed_key);
  

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

  // add Tg,Te to output buffer

  //  std::fill(out_bufr.begin(),out_bufr.end(),0);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&out_bufr[0]),Tg);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&out_bufr[Env::key_size_in_bytes()]),Te);

  // combine half gates:
  out_key = _mm_xor_si128(Wg,We);

  //benchmark_time += MPI_Wtime() - start;

}
  


void GarbledCircuit::evlHalfGatePair(__m128i &current_key, __m128i & key1, __m128i & key2, Bytes & in_bufr){
  //assert(in_bufr.size() == 2*Env::key_size_in_bytes());

  //double start = MPI_Wtime();

  // get the select bits
  byte sa,sb;
  sa = _mm_extract_epi8(key1,0) & 0x01;
  sb = _mm_extract_epi8(key2,0) & 0x01;

  // get the counter values
  uint32_t j1, j2;
  j1 = increment_index();
  j2 = increment_index();

  __m128i j1_128, j2_128;
  j1_128 = _mm_set1_epi64x(j1);
  j2_128 = _mm_set1_epi64x(j2);

  // fprintf(stdout,"Half Gate In Buffer: %s\n",in_bufr.to_hex().c_str());

  Bytes tmp_bufr;
  __m128i Tg, Te;
  

  // TG is always sent first, and then TE
  // where TG is the single row transmitted for the Generator's half-gate
  // and TE is the single row transmitted for the Evaluator's half-gate
  Tg = _mm_loadu_si128(reinterpret_cast<__m128i*>(&in_bufr[0]));
  Te = _mm_loadu_si128(reinterpret_cast<__m128i*>(&in_bufr[0]+Env::key_size_in_bytes()));
  Tg = _mm_and_si128(Tg,m_clear_mask);
  Te = _mm_and_si128(Te,m_clear_mask);


  __m128i Wa,Wb;

  Wa = key1;
  Wb = key2;

  __m128i H_Waj1,H_Wbj2; // output of hashes
  H_Pi(H_Waj1, Wa, j1_128, m_clear_mask, m_fixed_key);
  H_Pi(H_Wbj2, Wb, j2_128, m_clear_mask, m_fixed_key);


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

  //benchmark_time += MPI_Wtime() - start;

}  


void GarbledCircuit::clear_garbling_bufr(){
  m_garbling_bufr.clear();
  //std::fill(m_garbling_bufr.begin(),m_garbling_bufr.end(),0);
}



// this function computes the 2-uhf matrix given by matrix
// including output gates
// collaboratively with Gen
// Evl does not learn the parity of the outputs though
// either i need to do output gates differently or 
// Eval needs that info from Gen
void GarbledCircuit::gen_next_hash_row(Bytes & row, Bytes & hash_bufr){
  

  if(m_hash_out.size()*8 <= m_hash_row_idx){
    m_hash_out.resize((m_hash_out.size()+1)*2,0);
  }

  // starter value for our row's key
  __m128i row_key;
  row_key = _mm_set_epi64x(0,0);
  
  Bytes next_key;
  __m128i next_key_128;
  
  __m128i output_key;
  
  for(int i = 0; i < m_gen_inputs->size()/2; i++){
    // should be approximately same size as row*8
    
    next_key = get_Gen_Key(i,get_Input_Parity(i));
    //next_key = get_Gen_Input(i);
    save_Key_to_128bit(next_key, next_key_128);
    
    if(row.get_ith_bit(i)==1){
      row_key = _mm_xor_si128(row_key, next_key_128);
    }    
  }
  
  hash_bufr.clear();
  
  // we use 5 for output keys
  genStandardGate(output_key, row_key, row_key, hash_bufr, 5);  
  
  // now get output bit
  // we probably need to garble our output gate for this
  if(_mm_extract_epi8(output_key,0)&0x01 == 1 ){
    m_hash_out.set_ith_bit(m_hash_row_idx,1);
  }

  m_hash_row_idx++;
 
}

void GarbledCircuit::evl_next_hash_row(Bytes & row, Bytes & in_bufr){

  if(m_hash_out.size()*8 <= m_hash_row_idx){
    m_hash_out.resize((m_hash_out.size()+1)*2,0);
  }
  
  __m128i row_key;
  row_key = _mm_set_epi64x(0,0);
  
  Bytes next_key;
  __m128i next_key_128;
  
  //assert(m_gen_inputs->size()>=row.size()*8);
  for(int i = 0; i < m_gen_inputs->size(); i++){
    // should be approximately same size as row*8
    next_key = get_Gen_Input(i);
    save_Key_to_128bit(next_key, next_key_128);
    if(row.get_ith_bit(i)==1){
      row_key = _mm_xor_si128(row_key, next_key_128);
    }
  }
    
  __m128i output_key;
  evlStandardGate(output_key, row_key, row_key, in_bufr);
  in_bufr.clear();

  // now get output bit
  // we probably need to garble our output gate for this
  if(_mm_extract_epi8(output_key,0)&0x01 == 1 ){
    m_hash_out.set_ith_bit(m_hash_row_idx,1);
  }

  m_hash_row_idx++;
  
}


void evaluate_K_Probe_Matrix(std::vector<Bytes> &matrix){
  // here, Eval derives keys with the proper semantics
  // to use as her inputs

}

void generate_K_Probe_Matrix(std::vector<Bytes> &matrix){
  // here, Gen derives keys for Eval with the proper semantics
  // that will be her "new" input keys

}


Bytes GarbledCircuit::get_alice_out(){
  //  fprintf(stdout,"benchmark time: %f\n",benchmark_time);
  //fprintf(stdout,"benchmark time2: %f\n",btime2);


  fprintf(stdout,"xor gates  : %i \t xor time  : %f\n",xor_gates,xor_time);
  fprintf(stdout,"half gates : %i \t hgate time: %f\n",half_gates,hg_time);
  fprintf(stdout,"other gates: %i \t other time: %f\n",other_gates,og_time);
  fprintf(stdout,"total gates: %i \t total time: %f\n",total_gates,garble_time);
  fprintf(stdout,"num copies : %i \t copy time : %f\n",num_copies, copy_time);
  fprintf(stdout,"num buffers: %i \tbuffer time: %f\n",num_buffers, buffer_time);
  fprintf(stdout,"num b cpy  : %i \t b_cpy time: %f\n",num_b_cpy, b_cpy_time);
  fprintf(stdout,"num comm   : %i \t comm time: %f\n",num_comms, comm_time);
  //fprintf(stdout,"num send   : %i \t send time: %f\n",num_sends, send_time);

  return m_alice_out;
}

Bytes GarbledCircuit::get_bob_out(){
  //  fprintf(stdout,"benchmark time: %f\n",benchmark_time);

  fprintf(stdout,"xor gates  : %i \t xor time  : %f\n",xor_gates,xor_time);
  fprintf(stdout,"half gates : %i \t hgate time: %f\n",half_gates,hg_time);
  fprintf(stdout,"other gates: %i \t other time: %f\n",other_gates,og_time);
  fprintf(stdout,"total gates: %i \t total time: %f\n",total_gates,garble_time);
  fprintf(stdout,"num copies : %i \t copy time : %f\n",num_copies, copy_time);
  fprintf(stdout,"num buffers: %i \tbuffer time: %f\n",num_buffers, buffer_time);
  fprintf(stdout,"num b cpy  : %i \t b_cpy time: %f\n",num_b_cpy, b_cpy_time);
  fprintf(stdout,"num comm   : %i \t comm time: %f\n",num_comms, comm_time);
  //fprintf(stdout,"num send   : %i \t send time: %f\n",num_sends, send_time);
  

  return m_bob_out;
}

Bytes GarbledCircuit::get_hash_out(){
  return m_hash_out;
}

 void GarbledCircuit::trim_output_buffers(){
  if(m_alice_out_ix>0)
    m_alice_out.resize(m_alice_out_ix/8,0);
  if(m_bob_out_ix>0)
    m_bob_out.resize(m_bob_out_ix/8,0);


}

Bytes GarbledCircuit::get_Gen_Output_Label(uint32_t idx){

  return m_gen_output_labels[idx];

}


void GarbledCircuit::send_half_gate(const Bytes &buf){
  //clock_t start_t;
  //start_t = clock();

  // std::cout << "enqueue half gate " << std::endl << buf.to_hex() << std::endl; 
  
  num_comms++;
  clock_gettime(CLOCK_REALTIME, &comm_start);

  // Env::remote()->write_2_ciphertexts(buf);
  enqueue_messages(buf,2);
  //m_comm_time += (double) (clock() - start_t)/CLOCKS_PER_SEC;

  
  clock_gettime(CLOCK_REALTIME, &comm_end);
  comm_time += ( comm_end.tv_sec - comm_start.tv_sec )
    + ( comm_end.tv_nsec - comm_start.tv_nsec )
    / BILN;

}

void GarbledCircuit::send_full_gate(const Bytes &buf){
  //clock_t start_t;
  //start_t = clock();
  
  //std::cout << "send full gate " << std::endl << buf.to_hex() << std::endl; 

#ifdef GRR
  enqueue_messages(buf,3);
  //Env::remote()->write_3_ciphertexts(buf);
#else
  enqueue_messages(buf,4);
  //Env::remote()->write_4_ciphertexts(buf);
#endif
  
  //  std::cout << "enqueue full gate " << std::endl << buf.to_hex() << std::endl; 

  //  m_comm_time += (double) (clock() - start_t)/CLOCKS_PER_SEC;

}

void GarbledCircuit::read_half_gate(Bytes & buf){
  //clock_t start_t;
  //start_t = clock();

  num_comms++;
  clock_gettime(CLOCK_REALTIME, &comm_start);
  
  //Bytes ret =  Env::remote()->read_2_ciphertexts();
  retrieve_ciphertexts(buf,2);

  
  
  clock_gettime(CLOCK_REALTIME, &comm_end);
  comm_time += ( comm_end.tv_sec - comm_start.tv_sec )
    + ( comm_end.tv_nsec - comm_start.tv_nsec )
    / BILN;
  

  //  std::cout << "read half gate " << std::endl << ret.to_hex() << std::endl; 

  // return ret;

    //m_comm_time += (double) (clock() - start_t)/CLOCKS_PER_SEC;

}

void GarbledCircuit::read_full_gate(Bytes & buf){
  //clock_t start_t;
  // start_t = clock();
  

#ifdef GRR
  retrieve_ciphertexts(buf,3);
  //Bytes ret = retrieve_ciphertexts(3);
  //return Env::remote()->read_3_ciphertexts();
#else
  retrieve_ciphertexts(buf,4);
//Bytes ret =  retrieve_ciphertexts(4);
  //return Env::remote()->read_4_ciphertexts();
#endif
  //std::cout << "read full gate " << std::endl << ret.to_hex() << std::endl; 
    
  //return ret;
  //  m_comm_time += (double) (clock() - start_t)/CLOCKS_PER_SEC;

}

void GarbledCircuit::enqueue_messages(const Bytes & source, uint32_t num){
  if(num + m_messages_waiting < m_message_limit){
    add_messages_to_queue(source,num);

  } else if(num + m_messages_waiting == m_message_limit) { 
    
    add_messages_to_queue(source,num);
    send_buffer();
    
  } else {
    int send_early = m_message_limit - m_messages_waiting;
    
    add_messages_to_queue(source, send_early);
    
    send_buffer();

    m_message_queue.insert(m_message_queue.begin(),
                           source.begin() + send_early*Env::key_size_in_bytes(), 
                           source.begin() + send_early*Env::key_size_in_bytes()
                           + (num - send_early)*Env::key_size_in_bytes());
    m_messages_waiting = (num - send_early);
  }
    
}
    
void GarbledCircuit::add_messages_to_queue(const Bytes & src, uint32_t num){
  m_message_queue.insert(m_message_queue.begin()+m_messages_waiting*Env::key_size_in_bytes(),src.begin(),src.begin()+num*Env::key_size_in_bytes());

  m_messages_waiting += num;
}

void GarbledCircuit::send_buffer(){
  // this is also like a flush
  //std::cout << "buffer send " << std::endl << m_message_queue.to_hex() << std::endl; 

  Env::remote()->write_n_ciphertexts(m_message_queue, m_message_limit);
  m_messages_waiting = 0;
  std::fill(m_message_queue.begin(),m_message_queue.end(),0);
  
  //std::cout << "sent" << std::endl;

}

void GarbledCircuit::retrieve_ciphertexts(Bytes & buf, uint32_t num_ctexts){
  //Bytes ret;
  //m_ciphertext_buff.clear();
  
  if( num_ctexts < m_messages_waiting ){
    // get num messages from the queue
    // advance the queue index / reduce the number waiting
    buf.insert(buf.begin(),
               m_message_queue.begin() + m_queue_index*Env::key_size_in_bytes(),
               m_message_queue.begin() + m_queue_index*Env::key_size_in_bytes() + num_ctexts*Env::key_size_in_bytes());
    m_queue_index += num_ctexts;
    m_messages_waiting -= num_ctexts;
    
    //assert(m_queue_index + m_messages_waiting == m_message_limit);
    
    // return m_ciphertext_buff;
    
  } else if (num_ctexts == m_messages_waiting ){
    // get the last messages
    // reset the queue index/ reset the number waiting
    // retrieve the next set of messages
    
    assert(num_ctexts + m_queue_index == m_message_limit);
    buf.insert(buf.begin(),             
               m_message_queue.begin() + m_queue_index*Env::key_size_in_bytes(),
               m_message_queue.begin() + m_queue_index*Env::key_size_in_bytes() + num_ctexts*Env::key_size_in_bytes());
    
    // fprintf(stdout,"case 2\n");
    retrieve_buffer();
    
    // assert(m_queue_index + m_messages_waiting == m_message_limit);

    // return m_ciphertext_buff;

  } else{  // num_ctexts > m_messages_waiting
    // get the last ones left
    // reset the queue index/ reset the number waiting
    // get the remaining ones we need and update appropriately
    uint32_t get_early = m_messages_waiting;
    uint32_t get_later = num_ctexts - get_early;

    if(get_early > 0){ // should always be true for this case
      buf.insert(buf.begin(),             
                 m_message_queue.begin() + m_queue_index*Env::key_size_in_bytes(),
                 m_message_queue.begin() + m_queue_index*Env::key_size_in_bytes() + get_early*Env::key_size_in_bytes());
      
      //   std::cout << "partial ciphertext" << m_ciphertext_buff.to_hex() << std::endl;
      // fprintf(stdout,"partial copy\n");
    }
    
    // fprintf(stdout,"case 3\n");
    retrieve_buffer();
    buf.insert(buf.begin() + get_early*Env::key_size_in_bytes(),
               m_message_queue.begin(),
               m_message_queue.begin() + get_later*Env::key_size_in_bytes());
    
    m_messages_waiting -= get_later;
    m_queue_index += get_later;
    
    //assert(m_queue_index + m_messages_waiting == m_message_limit);
    
    // return m_ciphertext_buff;
  }
  
}

void GarbledCircuit::retrieve_buffer(){
  // resets the queue index/ the number waiting
  // retrieves the next batch of messages from the wire
  
  std::fill(m_message_queue.begin(),m_message_queue.end(),0);
  m_queue_index = 0;
  m_messages_waiting = m_message_limit;
  Env::remote()->read_n_ciphertexts(m_message_queue,m_message_limit);
  
  // std::cout << "buffer retrieve " << std::endl << m_message_queue.to_hex() << std::endl; 

}

#endif
