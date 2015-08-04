#include "BetterYao5.h"
#include <algorithm>
#include "otextension.h"

#include <unistd.h>

#include <log4cxx/logger.h>
static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("BetterYao5.cpp"));

#define debug_evl_fprintf(a,b) if(!m_chks[ix]){ std::cout << a << "  " << b << "\t rank: " << Env::group_rank() << std::endl; }


BetterYao5::BetterYao5(EnvParams &params) : YaoBase(params)
{
  
  std::cout << "node load: " << Env::node_load() << std::endl;
  
  // Init variables
  m_gcs.resize(Env::node_load());
  for (size_t ix = 0; ix < m_gcs.size(); ix++)
    {
      m_gcs[ix] = GarbledCircuit();
    }

  get_and_size_inputs();

  
}

void BetterYao5::start()
{  
  // in HBC there is no commitment, 
  // we just use the function prototype from BY5
  fprintf(stdout,"Gen Generate Inputs\n");
  gen_generate_and_commit_to_inputs();
  fprintf(stdout,"Eval Input OT\n");
  eval_input_OT();
  fprintf(stdout,"Transfer Evaluation Circuit Info\n");
  transfer_evaluation_circuit_info();
  fprintf(stdout,"Garble\n");
  garble_and_check_circuits();
  fprintf(stdout,"Outputs\n");
  retrieve_outputs();
}



/**
 *
 GENERAL-PURPOSE RANDOMNESS GENERATIONS
 *
 */

/**
   this function generates the randomnesses that will be used to seed the circuit Prngs
   and will be transfered in the special circuit OTs
*/
void BetterYao5::generate_random_seeds(std::vector<Bytes> & seeds,uint32_t num_seeds){
 
  Bytes rand;
  G rand_elem;
  //Prng prng = Prng();
  
  seeds.clear();

  for(int i=0;i<num_seeds;i++){
    // these random elements will be used for OT,
    // so we generate them as random group elements

    
    // this version useful for ALSZ extension
    //rand = prng.rand_bits(Env::elm_size_in_bytes()*8);
    //std::cout << "prng rand element: " << rand.to_hex() << std::endl;
    //seeds.push_back(rand);
    
    rand_elem.random();
    seeds.push_back(rand_elem.to_bytes());
  }
}

/**
   given a vector of prngs and a vector of seeds,
   seed the prngs
 */
void BetterYao5::seed_prngs(std::vector<Prng> & prngs, std::vector<Bytes> & seeds){
  assert(prngs.size() == seeds.size());
  
  for(int i=0;i<prngs.size();i++){
    prngs[i].seed_rand(seeds[i]);
  }
}


/**
   inputs: reference to Prng, location where keys will be stores, number of keys
   the length of the key is determined by the security parameter K
   we use a prng for repeatable results
   note that even keys are given semantic value 0 and odd keys are given semantic value 1
*/
void BetterYao5::generate_input_keys(Prng & prng, std::vector<Bytes> & dest, uint32_t num_keys, uint32_t num_bits){
  Bytes key;
  int j;
  for(j=0;j<num_keys;j++){
    // generate the key (so that it is the right size)
    key = prng.rand_bits(num_bits);
    // set the permutation bit. evens are 0, odds are 1
    //    key.set_ith_bit(0, j%2);
    dest.push_back(key);
  }

}



Bytes BetterYao5::get_gen_full_input(){
  GEN_BEGIN
    return m_private_input + get_gen_output_mask();
    //  return m_private_input + get_gen_output_mask() + get_gen_input_randomness();
  GEN_END

    EVL_BEGIN
    //assert(0);
    return m_prng.rand_bits(get_gen_full_input_size());
    EVL_END
}

Bytes BetterYao5::get_gen_output_mask(){
  // return Bytes(0);
  return m_gen_output_mask;
}

Bytes BetterYao5::get_gen_input_randomness(){
  return Bytes(0);
  //  return m_gen_aux_random_input;
}

uint32_t BetterYao5::get_gen_full_input_size(){
  // assume Gen's output mask is the same size as his input mask
  // returns the number of bits of gen's fill input size
  

  // in hbc, still need input keys and output masks
  return 2*get_gen_inp_size();
  //  return 2*get_gen_inp_size() + 2*Env::k() + ceil_log_base_2(Env::k()); 
}

uint32_t BetterYao5::get_gen_inp_size(){
  return m_gen_inp_cnt;
  
}

// assume output size = input size for sake of number of keys
uint32_t BetterYao5::get_gen_output_size(){
  return get_gen_inp_size();
}





/**
   GENERATE GEN INPUT KEYS
 */

void BetterYao5::gen_generate_and_commit_to_inputs(){

  GEN_BEGIN
  
  // first, generate random seeds and seed the circuit prngs
  m_circuit_prngs.resize(Env::node_load());
  generate_random_seeds(m_circuit_seeds, Env::node_load());
  seed_prngs(m_circuit_prngs,m_circuit_seeds);

  // next, generate Gen's input keys
  m_gen_inp_keys.resize(Env::node_load());
  m_gen_inp_permutation_bits.resize(Env::node_load());
  m_R.resize(Env::node_load());
  m_gen_select_bits.resize(Env::node_load());
  m_evl_inp_keys.resize(Env::node_load());


  // also generate randomness for output
  m_gen_output_mask = m_circuit_prngs[0].rand_bits(m_private_input.size()*8);
  
  //  for(int j = 0; j < Env::node_load();j++){
  generate_gen_input_keys(0);
  generate_eval_input_keys(0);
    //}
  
  GEN_END
 
}


/**
   generate Gen input keys and permutation bits
   using the circuit prng
 */
void BetterYao5::generate_gen_input_keys(uint32_t circuit_num){
  // this function should be executed in parallel by all processors

  // first, generate all of the input keys for Gen
  // (K_{0},K_{1}) <-$ {0,1}^{2k}
  m_R[circuit_num] = m_circuit_prngs[circuit_num].rand_bits(Env::k());
  m_R[circuit_num].set_ith_bit(0,1);
 
  std::vector<Bytes> random_keys;

  // generate enough input keys for Gen's entire input
  // and then compute their XOR-offsets
  generate_input_keys(m_circuit_prngs[circuit_num],
                      random_keys,
                      get_gen_full_input_size(),
                      Env::k());
  
  // the proper XOR-offsets get pushed alongside the input keys
  for(int i = 0; i < get_gen_full_input_size();i++){
    random_keys[i].set_ith_bit(0,0);
    m_gen_inp_keys[circuit_num].push_back(random_keys[i]);
    m_gen_inp_keys[circuit_num].push_back(random_keys[i] ^ m_R[circuit_num]);
  }
  
  // next generate permutation bits for Gen's input keys
  // { pi_{i} } <-$ {0,1} for i in (0, gen_inputs ]  
  m_gen_inp_permutation_bits[circuit_num] = m_circuit_prngs[circuit_num].rand_bits(get_gen_full_input_size());
 
  // now we have a situation where we have a bunch of permutation bits
  // and a bunch of keys, ordered [0,1,0,1] etc.
  // the idea is that the permutation bit tells us if 
  // the 0 key encodes a 0 (0th bit 0) or a 1 (0th bit 1)
  // Gen also has his own input bits
  // so we have the following table to describe which key Gen sends
  // depending on his input key and the permutation bit
  //   Gen Input     Permutation    Selected
  //      0       |       0      |     0
  //      0       |       1      |     1
  //      1       |       0      |     1
  //      1       |       1      |     0
  // which is easily an XOR
 
  m_gen_select_bits[circuit_num] = m_gen_inp_permutation_bits[circuit_num] ^ get_gen_full_input();
  
}


/*
  Gen (or Eval) generates Gen's input keys for the ith circuit
  where i is given in circuit_num
  note that these keys do not need their permutation/selection bits set
  (so that their semantics can be derived), since the key will be
  used as a decryption key during an Eval input gate
*/
void BetterYao5::generate_eval_input_keys(uint32_t circuit_num){
    

  //Bytes rand_key;
  G rnd;
  for(int i = 0; i < get_evl_inp_count()*2;i++){

    // this implementation will replace the current one when
    // we link OT extensions
    // rand_key = m_circuit_prngs[circuit_num].rand_bits(Env::elm_size_in_bytes()*8);
    // rnd.from_bytes(rand_key);

    // m_evl_inp_keys[circuit_num].push_back(rand_key);
  
    rnd.random(m_circuit_prngs[circuit_num]);
    m_evl_inp_keys[circuit_num].push_back(rnd.to_bytes());
      
  }
  
}




void BetterYao5::eval_input_OT(){

  GEN_BEGIN
  assert(m_evl_inp_keys.size() == Env::node_load());
  m_evl_hashed_inp_keys.resize(Env::node_load());
  m_evl_hashed_inp_keys[0].resize(m_evl_inp_keys[0].size());

  ot_send_batch(m_evl_inp_keys);

  fprintf(stdout,"hash eval input keys");
  hash_eval_input_keys(m_evl_inp_keys[0],m_evl_hashed_inp_keys[0],Env::k());

  GEN_END

  EVL_BEGIN
    std::vector<std::vector<Bytes> > evl_receive_keys;
    assert(evl_receive_keys.size() == 0);
    assert(m_evl_received_keys.size() == 0);
    evl_receive_keys.resize(Env::node_load());
    m_evl_received_keys.resize(Env::node_load());
      
    // remember we require that the private input length be a multiple of 8
    ot_receive_batch(explode_vector(m_private_input, m_private_input.size()*8),evl_receive_keys);
    
    //for(int i = 0; i < Env::node_load(); i++){
    int i = 0;
    m_evl_received_keys[i].resize(evl_receive_keys[i].size());
    hash_eval_input_keys(evl_receive_keys[i],m_evl_received_keys[i],Env::k());
      //}

    EVL_END
    
}


void BetterYao5::hash_eval_input_keys(std::vector<Bytes> & source, std::vector<Bytes> & destination, uint32_t num_bits){
  assert(source.size()==destination.size());

  for(int j = 0; j < source.size(); j++){
    
    destination[j] = source[j].hash(num_bits);
  }
}




// evl needs:
// gen's input keys
// constant keys
void BetterYao5::transfer_evaluation_circuit_info(){

  GEN_BEGIN
    
    // constant keys
  m_const_0_keys.resize(1);
  m_const_1_keys.resize(1);
  m_const_0_keys[0] = m_circuit_prngs[0].rand_bits(Env::k());
  m_const_1_keys[0] = m_circuit_prngs[0].rand_bits(Env::k());

  GEN_SEND(m_const_0_keys[0]);
  GEN_SEND(m_const_1_keys[0]);

  std::vector<Bytes > gen_send_keys;

  // Gen has to figure out which keys to send
  assert(m_gen_inp_keys[0].size() == 2*get_gen_full_input_size());
  for(int i = 0; i < m_gen_inp_keys[0].size()/2;i++){

    gen_send_keys.push_back(m_gen_inp_keys[0][2*i+m_gen_select_bits[0].get_ith_bit(i)]);
  }

  assert(gen_send_keys.size() == get_gen_full_input_size());
  for(int i = 0; i < gen_send_keys.size();i++){
    GEN_SEND(gen_send_keys[i]);
  }

  // key generation seed
  G rnd;
  m_key_generation_seeds.resize(1);
  rnd.random();
  m_key_generation_seeds[0] = rnd.to_bytes();
  //GEN_SEND(m_key_generation_seeds[0]);


  GEN_END

  EVL_BEGIN
  
    // constant keys
  m_const_0_keys.resize(1);
  m_const_1_keys.resize(1);  
  m_const_0_keys[0] = EVL_RECV();
  m_const_1_keys[0] = EVL_RECV();
  
  // gen input keys
  m_gen_inp_keys.resize(1);
  for(int i = 0; i < get_gen_full_input_size(); i++){
    m_gen_inp_keys[0].push_back(EVL_RECV());
  }

  // key generation seed
  //  m_key_generation_seeds.resize(1);
  //m_key_generation_seeds[0] = EVL_RECV();
  
 
  EVL_END

}


void BetterYao5::garble_and_check_circuits(){
  initialize_circuits();
  evaluate_circuits();
}



void BetterYao5::initialize_circuits(){
  MPI_Barrier(m_mpi_comm);
  
  GEN_BEGIN
    int ix = 0;


  m_const_1_keys[ix] = m_const_1_keys[ix] ^ m_R[ix];
  
  m_gcs[ix].init_Generation_Circuit(&m_gen_inp_keys[ix],// gen input keys
                                    &m_evl_hashed_inp_keys[ix], // evl input keys
                                    get_gen_inp_size(),// gen input size
                                    m_key_generation_seeds[ix], // random seed
                                    m_gen_inp_permutation_bits[ix], // permutation bits
                                    m_R[ix], // circuit XOR offset
                                    m_const_0_keys[ix], m_const_1_keys[ix]);// constant keys
  


  m_gcs[ix].m_st = 
    load_pcf_file(Env::pcf_file(), m_gcs[ix].get_Const_Wire(0), m_gcs[ix].get_Const_Wire(1), copy_key);
  m_gcs[ix].m_st->alice_in_size = get_gen_full_input_size();
  m_gcs[ix].m_st->bob_in_size = get_evl_inp_count();
  
  set_external_circuit(m_gcs[ix].m_st, &m_gcs[ix]);
  
  m_gcs[ix].set_Gen_Circuit_Functions();
  
  
  fprintf(stderr,"Gen Input: \t%s\nGen Full Input: \t%s\nGen Permutation Bits: \t%s\n",
          m_private_input.to_hex().c_str(),
          get_gen_full_input().to_hex().c_str(),
          m_gen_inp_permutation_bits[ix].to_hex().c_str());
  

  GEN_END
  
  EVL_BEGIN

    //  if(Env::group_rank() == 0){
      // begin with one, just to make debugging easier
    // for(int ix = 0; ix < m_gcs.size();ix++){
    int ix = 0;    
       
          // evaluation circuit
  
  m_gcs[ix].init_Evaluation_Circuit(&m_gen_inp_keys[ix],// gen keys
                                    &m_evl_received_keys[ix],//evl keys
                                    get_gen_inp_size(),// gen inp size
                                    m_private_input, // private input
                                    m_const_0_keys[ix], //constant keys
                                    m_const_1_keys[ix]);
  
    
  m_gcs[ix].m_st = 
    load_pcf_file(Env::pcf_file(), m_gcs[ix].get_Const_Wire(0), m_gcs[ix].get_Const_Wire(1), copy_key);
    m_gcs[ix].m_st->alice_in_size = get_gen_full_input_size();
    m_gcs[ix].m_st->bob_in_size = get_evl_inp_count();
    
    set_external_circuit(m_gcs[ix].m_st, &m_gcs[ix]);
    
    m_gcs[ix].set_Evl_Circuit_Functions();
    
  EVL_END

    }


void BetterYao5::evaluate_circuits(){
  

  double start;
  double garble_time =0.0;
  double comm_time =0.0;

  GEN_BEGIN
  
  int ix = 0;
  Bytes bufr;

  start = MPI_Wtime();
  while(get_next_gate(m_gcs[ix].m_st)){
    garble_time += MPI_Wtime() - start;
    start = MPI_Wtime();

    bufr = m_gcs[ix].get_garbling_bufr();
    GEN_SEND(bufr);
    m_gcs[ix].clear_garbling_bufr();
    
    comm_time += MPI_Wtime() - start;
    start = MPI_Wtime();
  }
  garble_time += MPI_Wtime()-start;

  GEN_SEND(Bytes(0)); // redundant value to prevent Evl from hanging
  
  GEN_END
    
  EVL_BEGIN
  
  int ix =0;
  Bytes bufr;
  

  start = MPI_Wtime();
  do {

    garble_time += MPI_Wtime()-start;
    start = MPI_Wtime();

    bufr = EVL_RECV();
    m_gcs[ix].set_garbling_bufr(bufr);
    
    comm_time += MPI_Wtime()-start;
    start = MPI_Wtime();

  } while(get_next_gate(m_gcs[ix].m_st));

  garble_time += MPI_Wtime()-start;

  EVL_END

    fprintf(stdout,"garble: %f\n",garble_time);
    fprintf(stdout,"comm: %f\n",comm_time);
  

}


/**
   STEP 8: RETRIEVE OUTPUTS
 */

void BetterYao5::retrieve_outputs(){
  // this function currently outputs in the clear
  // TODO: gen output authenticity and masking gen's outputs

  GEN_BEGIN

  assert(m_gcs.size()>0);  

  int i = 0;
  
  m_gcs[i].trim_output_buffers();

  Bytes alice_out_parity = m_gcs[i].get_alice_out();
  GEN_SEND(alice_out_parity);
  

  Bytes bob_out_evaluated = GEN_RECV();
  
  //  if(bob_out_evaluated.size() > 0){
    // TODO: Gen ultimately won't send this information to Eval
  Bytes bob_out_parity = m_gcs[i].get_bob_out();
  GEN_SEND(bob_out_parity);
  Bytes bob_out = bob_out_evaluated ^ bob_out_parity;
  Bytes bob_mask = get_gen_output_mask();
  bob_mask.resize(bob_out_evaluated.size(),0);
  Bytes bob_decrypted = bob_mask ^ bob_out;
  
  
  std::cout << "bob out (parity): " << bob_out_parity.to_hex() << std::endl;
  std::cout << "bob out (parity ^ eval): " << bob_out.to_hex() << std::endl;
  std::cout << "bob out (mask): " << bob_mask.to_hex() << std::endl;
  std::cout << "bob out: " << bob_decrypted.to_hex() << std::endl; 
  //}
  
  
  GEN_END

  EVL_BEGIN
    
  int i = 0;

  m_gcs[i].trim_output_buffers();
  
  Bytes alice_out = m_gcs[i].get_alice_out();
  std::cout << "alice out (masked): " << alice_out.to_hex() << std::endl;
  
  Bytes alice_out_parity = EVL_RECV();
  alice_out = alice_out ^ alice_out_parity;
      
  // TODO: eval won't receive this info. Gen's outputs will 
  // be commincated properly
  Bytes bob_out = m_gcs[i].get_bob_out();
  EVL_SEND(bob_out);
  Bytes bob_out_parity = EVL_RECV();
  bob_out = bob_out ^ bob_out_parity;
  
  std::cout << "alice out (parity): " << alice_out_parity.to_hex() << std::endl;
  std::cout << "alice out  (final): " << alice_out.to_hex() << std::endl;      
  //  std::cout << "alice check hash: " << m_2UHF_hashes[i].to_hex() << std::endl;
  if(bob_out.size()>0){
    std::cout << "bob out  (final):" << bob_out.to_hex() << std::endl;
  }
 


  EVL_END

}



/**
 *
 OBLIVIOUS TRANSFER
 *
 */

void BetterYao5::ot_send(std::vector<Bytes> & sender_inputs){
  //OT_alsz_send(Env::ip_server(), Env::IP_SERVER_PORT + Env::node_amnt() + Env::group_rank(), sender_inputs.size(), Env::k(), Env::s(), sender_inputs, sender_inputs);

  ot_send_init();
  ot_send_random(sender_inputs);
}

void BetterYao5::ot_receive(Bytes selection_bits, std::vector<Bytes> & results_container){
  //OT_alsz_recv(Env::ip_server(), Env::IP_SERVER_PORT + Env::node_amnt() + Env::group_rank(), results_container.size(), Env::k(), Env::s(), selection_bits, results_container);

  ot_receive_init();
  ot_receive_random(selection_bits, results_container);
}


void BetterYao5::ot_send_batch(std::vector<std::vector<Bytes> > & sender_inputs){
  ot_send_init();
  for(int j = 0; j < sender_inputs.size(); j++){
    //if(Env::group_rank() == 0)
    //for(int i = 0; i < sender_inputs[j].size();i++){
      //  fprintf(stderr,"ot send input (%i): %s\n",i,sender_inputs[j][i].to_hex().c_str());
      //   }
    ot_send_random(sender_inputs[j]);
  }
}


void BetterYao5::ot_receive_batch(Bytes selection_bits, std::vector<std::vector<Bytes> > & results_container){
  ot_receive_init();
  for(int j = 0; j < results_container.size();j++){
    ot_receive_random(selection_bits, results_container[j]);
    //if(Env::group_rank() == 0)
      //for(int i = 0; i < results_container[j].size();i++){
      //   fprintf(stderr,"ot receive input (%i): %s\n",i,results_container[j][i].to_hex().c_str());
      //}
  }
}

void BetterYao5::ot_send_init(){
  double start;
  
  start = MPI_Wtime();
  std::vector<Bytes> bufr_chunks;
  Bytes bufr(Env::elm_size_in_bytes()*4);
  
  Z y, a;

  // step 1: ZKPoK of the CRS: g[0], h[0], g[1], h[1]
  if (Env::is_root())
    {
      
      //start = MPI_Wtime();
      GEN_BEGIN
      bufr = GEN_RECV();
      GEN_END
      EVL_BEGIN
      bufr = EVL_RECV();
      EVL_END
      // m_timer_com += MPI_Wtime() - start;
      
      // m_comm_sz += bufr.size();
    }
  
  // send g[0], g[1], h[0], h[1] to slave processes
  MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm);
  
  bufr_chunks = bufr.split(Env::elm_size_in_bytes());
  
  m_ot_g[0].from_bytes(bufr_chunks[0]);
  m_ot_g[1].from_bytes(bufr_chunks[1]);
  m_ot_h[0].from_bytes(bufr_chunks[2]);
  m_ot_h[1].from_bytes(bufr_chunks[3]);
        
  // group element pre-processing
  m_ot_g[0].fast_exp();
  m_ot_g[1].fast_exp();
  m_ot_h[0].fast_exp();
  m_ot_h[1].fast_exp();
  
}

void BetterYao5::ot_receive_init(){
  double start;

  start = MPI_Wtime();
  std::vector<Bytes> bufr_chunks;
  Bytes bufr(Env::elm_size_in_bytes()*4);

  Z y, a;
  // m_timer_gen += MPI_Wtime() - start;
  // m_timer_evl += MPI_Wtime() - start;

  // step 1: ZKPoK of the CRS: g[0], h[0], g[1], h[1]

  if(Env::is_root()){
    y.random();
    a.random();
    
    m_ot_g[0].random();
    m_ot_g[1] = m_ot_g[0]^y;          // g[1] = g[0]^y
    
    m_ot_h[0] = m_ot_g[0]^a;          // h[0] = g[0]^a
    m_ot_h[1] = m_ot_g[1]^(a + Z(1)); // h[1] = g[1]^(a+1)
    
    bufr.clear();
    bufr += m_ot_g[0].to_bytes();
    bufr += m_ot_g[1].to_bytes();
    bufr += m_ot_h[0].to_bytes();
    bufr += m_ot_h[1].to_bytes();
    // m_timer_evl += MPI_Wtime() - start;
  

    GEN_BEGIN
    GEN_SEND(bufr);
    GEN_END
    EVL_BEGIN
    EVL_SEND(bufr);
    EVL_END
  }

  // send g[0], g[1], h[0], h[1] to slave processes
  MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm);

  bufr_chunks = bufr.split(Env::elm_size_in_bytes());

  m_ot_g[0].from_bytes(bufr_chunks[0]);
  m_ot_g[1].from_bytes(bufr_chunks[1]);
  m_ot_h[0].from_bytes(bufr_chunks[2]);
  m_ot_h[1].from_bytes(bufr_chunks[3]);

  // group element pre-processing
  m_ot_g[0].fast_exp();
  m_ot_g[1].fast_exp();
  m_ot_h[0].fast_exp();
  m_ot_h[1].fast_exp();

}

void BetterYao5::ot_send_random(std::vector<Bytes> & send_inputs){
  double start;
  
  Bytes send, recv;
  std::vector<Bytes> recv_chunks;
  
  Z r, s[2], t[2];
  G gr, hr, X[2], Y[2];
    
  //  GEN_BEGIN // generator (OT sender)
  for (size_t bix = 0; bix < send_inputs.size()/2; bix++)
    {
      // Step 1: gr=g[b]^r, hr=h[b]^r, where b is the receiver's bit
      start = MPI_Wtime();
      recv.clear();
      
      // depending on who's calling,
      // receive
      GEN_BEGIN
        recv += GEN_RECV(); // receive gr, hr
      GEN_END
      EVL_BEGIN
        recv += EVL_RECV();
      EVL_END
        
        //      m_timer_com += MPI_Wtime() - start;
      
        //  m_comm_sz += recv.size();
        
      // Step 2: the evaluator computes X[0], Y[0], X[1], Y[1]
      start = MPI_Wtime();
      recv_chunks = recv.split(Env::elm_size_in_bytes());
      
      gr.from_bytes(recv_chunks[0]);
      hr.from_bytes(recv_chunks[1]);
      
      // retrieve keys from the supplied send_inputs vector
      // generate them from bytes rather than get them directly
      //Y[0] = send_inputs[bix*2]; // K[0]
      //Y[1] = send_inputs[bix*2+1]; // K[1]
      Y[0].from_bytes(send_inputs[bix*2]);
      Y[1].from_bytes(send_inputs[bix*2+1]);
      //if(Env::group_rank() == 0){
      //  fprintf(stderr,"new Y[0]: %s\n", Y[0].to_bytes().to_hex().c_str());
      //  fprintf(stderr,"new Y[1]: %s\n", Y[1].to_bytes().to_hex().c_str());
      //}

      s[0].random(); s[1].random();
      t[0].random(); t[1].random();
      
      // X[b] = ( g[b]^s[b] ) * ( h[b]^t[b] ) for b = 0, 1
      X[0] = m_ot_g[0]^s[0]; X[0] *= m_ot_h[0]^t[0];
      X[1] = m_ot_g[1]^s[1]; X[1] *= m_ot_h[1]^t[1];
      
      // Y[b] = ( gr^s[b] ) * ( hr^t[b] ) * K[b] for b = 0, 1
      Y[0] *= gr^s[0]; Y[0] *= hr^t[0];
      Y[1] *= gr^s[1]; Y[1] *= hr^t[1];
      
      send.clear();
      send += X[0].to_bytes();
      send += X[1].to_bytes();
      send += Y[0].to_bytes();
      send += Y[1].to_bytes();
      
      // depending on who's calling it,
      // send the buffer
      GEN_BEGIN
        GEN_SEND(send);
      GEN_END
      EVL_BEGIN
        EVL_SEND(send);
      EVL_END
    }
  
}

void BetterYao5::ot_receive_random(Bytes selection_bits, std::vector<Bytes> & results_container){
  double start;

  start = MPI_Wtime();
  Bytes send, recv;
  std::vector<Bytes> recv_chunks;

  Z r, s[2], t[2];
  G gr, hr, X[2], Y[2];

  // m_timer_gen += MPI_Wtime() - start;
  // m_timer_evl += MPI_Wtime() - start;

  for (size_t bix = 0; bix < selection_bits.size(); bix++)
    {
      // Step 1: gr=g[b]^r, hr=h[b]^r, where b is the receiver's bit
      start = MPI_Wtime();
      
      int bit_value = selection_bits[bix];

      r.random();
      
      gr = m_ot_g[bit_value]^r;
      hr = m_ot_h[bit_value]^r;
            
      send.clear();
      send += gr.to_bytes();
      send += hr.to_bytes();
      m_timer_evl += MPI_Wtime() - start;
            
      start = MPI_Wtime();
      
      GEN_BEGIN
      GEN_SEND(send);
      GEN_END
      EVL_BEGIN
      EVL_SEND(send);
      EVL_END
      

      // Step 2: the evaluator computes X[0], Y[0], X[1], Y[1]
      recv.clear();

      GEN_BEGIN
        recv += GEN_RECV();
      GEN_END
      EVL_BEGIN
        recv += EVL_RECV(); // receive X[0], Y[0], X[1], Y[1]
      EVL_END
      
      m_timer_com += MPI_Wtime() - start;
            
      m_comm_sz += send.size() + recv.size();
            
      // Step 3: the evaluator computes K = Y[b]/X[b]^r
      start = MPI_Wtime();
      recv_chunks = recv.split(Env::elm_size_in_bytes());
            
      X[bit_value].from_bytes(recv_chunks[    bit_value]); // X[b]
      Y[bit_value].from_bytes(recv_chunks[2 + bit_value]); // Y[b]
            
      // K = Y[b]/(X[b]^r)
      Y[bit_value] /= X[bit_value]^r;
      results_container.push_back(Y[bit_value].to_bytes());
      m_timer_evl += MPI_Wtime() - start;
    }
 
  // assert(results_container.size() = selection_bits.size()*8);
  // not quite the right assertion because number of bits need not be
  // an even multiple of 8.
}

