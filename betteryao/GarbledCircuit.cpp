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

  // Bytes bufr;
  do {
    // do things?
    
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

  fprintf(stderr,"gen next gate\n");

  // returns void pointer, which is pointer to a key?
  // use this one to call the Garbled Circuit object again
  GarbledCircuit &cct = *reinterpret_cast<GarbledCircuit*>(get_external_circuit(st));

  fprintf(stderr,"gen next gate\n");
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
  
    // will return current_key;

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

#ifdef FREE_XOR
    if(current_gate->truth_table == 0x06){ // if XOR gate
      current_key = _mm_xor_si128
        (
         *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1)),
         *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2))
         );
    } else {
#endif

      // we garble this gate
      // first, get the inputs

      // X and Y are input, Z is output
      __m128i X[2], Y[2], Z[2];
      // we need a couple of 128-bit variables
      __m128i aes_plaintext, aes_ciphertext;
      __m128i aes_key[2];
      
      uint8_t bit;
      
      // load the inputs from the PCF state container
      X[0] = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire1));
      Y[0] = *reinterpret_cast<__m128i*>(get_wire_key(m_st, current_gate->wire2));
      // and derive their XOR-complements
      X[1] = _mm_xor_si128(X[0], m_R); // X[1] = X[0] ^ R
      Y[1] = _mm_xor_si128(Y[0], m_R); // Y[1] = Y[0] ^ R
      
      // and get the permutation bits
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
      
      // now run the key derivation function using the keys and the gate index
      KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);

      // clear extra bits at the front
      aes_ciphertext = _mm_and_si128(aes_ciphertext, m_clear_mask);
      
      // extract the permutation bit for this wire?
      bit = (current_gate->truth_table >> (3-de_garbled_ix)) & 0x01;

#ifdef GRR
      // GRR technique: using zero entry's key as one of the output keys
      // the output key is the encrypted gate index
      _mm_store_si128(Z+bit, aes_ciphertext);
      // and other output is an offset
      Z[1-bit] = _mm_xor_si128(Z[bit], m_R);
      // and load it into current_key to return it to the circuit's state container
      current_key = _mm_load_si128(Z);
      
#else
      // otherwise we generate a new random key for the 0th entry
      tmp = m_prng.rand_bits(Env::k());
      tmp.resize(16, 0);
      Z[0] = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
      Z[1] = _mm_xor_si128(Z[0], cct.m_R);

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
      Bytes tmp;
      tmp.resize(16,0);


      aes_key[0] = _mm_xor_si128(aes_key[0], m_R);
      // key derivation function
      KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
      // and then mask to remove unnecessary bits
      aes_ciphertext = _mm_and_si128(aes_ciphertext, m_clear_mask);
    
      // recover the permutation bit
      bit = (current_gate->truth_table>>(3-(0x01^de_garbled_ix)))&0x01;
      // and find the offset
      aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);
      // and we add tmp to the buffer we're sending to Eval

      garbling_bufr.insert(garbling_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());


      
      // encrypt the 2nd entry : (X[x], Y[1-y])
      aes_key[0] = _mm_xor_si128(aes_key[0], m_R);
      aes_key[1] = _mm_xor_si128(aes_key[1], m_R);
      
      KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
      aes_ciphertext = _mm_and_si128(aes_ciphertext, m_clear_mask);
      //bit = current_gate.m_table[0x02^de_garbled_ix];
      bit = (current_gate->truth_table>>(3-(0x02^de_garbled_ix)))&0x01;
      aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
      
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);

      garbling_bufr.insert(garbling_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
      
 //cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
      
      // encrypt the 3rd entry : (X[1-x], Y[1-y])
      aes_key[0] = _mm_xor_si128(aes_key[0], m_R);
      
      KDF256((uint8_t*)&aes_plaintext, (uint8_t*)&aes_ciphertext, (uint8_t*)aes_key);
      aes_ciphertext = _mm_and_si128(aes_ciphertext, m_clear_mask);
      //bit = current_gate.m_table[0x03^de_garbled_ix];
      bit = (current_gate->truth_table>>(3-(0x03^de_garbled_ix)))&0x01;
      aes_ciphertext = _mm_xor_si128(aes_ciphertext, Z[bit]);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(&tmp[0]), aes_ciphertext);

      garbling_bufr += tmp;

      GEN_SEND(garbling_bufr);
      garbling_bufr.clear();
      //cct.m_out_bufr.insert(cct.m_out_bufr.end(), tmp.begin(), tmp.begin()+Env::key_size_in_bytes());
      


      // as a last step, increment the gate index. it's important that we 
      // maintain this index.
      // Note: if we divorce the garbling from the processor running
      // each garbled circuit evaluation, we will remove this piece from
      // this location, and probably pass it as a parameter to the garbling
      // function
      
      m_gate_index++;
      
      return &current_key; // get_Gen_Key or create a random nonce


    }

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
    __m128i aes_key[2], aes_plaintext, aes_ciphertext;
    std::vector<Bytes> garbled_keys;
    Bytes tmp;
    __m128i a;
    Bytes::const_iterator it;

    garbling_bufr = EVL_RECV();
    
    aes_plaintext = _mm_set1_epi64x(m_gate_index);
    
    aes_key[0] = *reinterpret_cast<__m128i*>(get_wire_key(m_st,current_gate->wire1));
    aes_key[1] = *reinterpret_cast<__m128i*>(get_wire_key(m_st,current_gate->wire2));
    
    const uint8_t perm_x = _mm_extract_epi8(aes_key[0], 0) & 0x01;
    const uint8_t perm_y = _mm_extract_epi8(aes_key[1], 0) & 0x01;
    
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
        it = garbling_bufr.begin() + (garbled_ix-1)*Env::key_size_in_bytes();
        
        tmp.assign(it, it+Env::key_size_in_bytes());
        tmp.resize(16, 0);
        a = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
        current_key = _mm_xor_si128(aes_ciphertext, a);
      }
    //cct.m_in_bufr_ix += 3*Env::key_size_in_bytes();
#else
    //it = cct.m_in_bufr_ix + garbled_ix*Env::key_size_in_bytes();
    it = garbling_bufr.begin() + (garbled_ix)*Env::key_size_in_bytes();
    tmp.assign(it, it+Env::key_size_in_bytes());
    tmp.resize(16, 0);
    current_key = _mm_loadu_si128(reinterpret_cast<__m128i*>(&tmp[0]));
    current_key = _mm_xor_si128(current_key, aes_ciphertext);
    
    garbling_bufr.clear();

    //cct.m_in_bufr_ix += 4*Env::key_size_in_bytes();
#endif


  }
  m_gate_index++;
  return &current_key;
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
