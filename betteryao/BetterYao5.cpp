#include "BetterYao5.h"
#include <algorithm>
#include "otextension.h"

#include <unistd.h>

#include <log4cxx/logger.h>
static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("BetterYao5.cpp"));

#define debug_evl_fprintf(a,b) if(!m_chks[ix]){ std::cout << a << "  " << b << "\t rank: " << Env::group_rank() << std::endl; }


BetterYao5::BetterYao5(EnvParams &params) : YaoBase(params), m_ot_bit_cnt(0)
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

// old start function
void BetterYao5::start()
{
  SS13();
}


// new start function
void BetterYao5::SS13(){

  modify_inputs();
  gen_generate_and_commit_to_inputs();
  agree_on_objective_circuit();
  gen_commit_to_io_labels();
  eval_input_OT();
  SS13_cut_and_choose();
  garble_and_check_circuits();
  retrieve_outputs();
  std::cout << "done protocol" << std::endl;
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
    key.set_ith_bit(0, j%2);
    dest.push_back(key);
  }


}


/**
 *
 COMMITMENTS
 *
 */


/**
   Generate commitments to informations stored in keys with a prng,
   then store them in the location given by commitments
 */
void BetterYao5::generate_commitments(Prng & prng, std::vector<Bytes> & source, std::vector<commitment_t> & dest){
  int i;
  for(i=0;i<source.size();i++){
    dest.push_back(make_commitment(prng,source[i]));
  }
}


/**
   Gen commits to inputs and transfers the commitments to Eval
 */
void BetterYao5::gen_send_evl_commitments(std::vector<commitment_t> & commits){
  
  Bytes committed;
  for(int j=0;j<commits.size();j++){
    committed = commit(commits[j]);
    GEN_SEND(committed); 
  }
}


/**
   Eval receives commitments from Gen and stores them in the container provided
 */
void BetterYao5::evl_receive_gen_commitments(std::vector<Bytes> & commits,uint32_t num_commitments){
  Bytes committed;
  for(int j=0;j< num_commitments;j++){
    committed = EVL_RECV();
    commits.push_back(committed);
  }
}


/**

   STEP 1: INPUT MODIFICATIONS

 */

void BetterYao5::modify_inputs(){
// Gen appends an output mask and extra randomness to his inputs
// Evl transforms her input in order to hide it from Gen
  GEN_BEGIN
  gen_generate_output_mask(m_prng);
  gen_generate_input_randomness(m_prng);
  GEN_END

  EVL_BEGIN
  choose_k_probe_resistant_matrix();
  evl_generate_new_input();
  EVL_END
}


void BetterYao5::gen_generate_output_mask(Prng & output_mask_prng){
  // assume for now that the output is at most the size of the input
  
  Bytes output_mask;

  if(Env::is_root()){
    // completely cover the input size (vector of bytes)
    output_mask = output_mask_prng.rand_bits(m_private_input.size()*8);
  }
  
   // send output mask to all processors
  output_mask.resize(m_private_input.size());
  MPI_Bcast(&output_mask[0],output_mask.size(),MPI_BYTE,0,m_mpi_comm);
  
  // save it to Gen's class instance
  m_gen_output_mask = Bytes(output_mask.begin(), output_mask.end());
  
}


void BetterYao5::gen_generate_input_randomness(Prng & input_prng){
  // compute ceil(lg_2(k)) and then generate some random bits 
  // these random bits are added as auxiliary inputs for Gen
  // in order to make the output of the 2-UHF appear random
  // to protect Gen's input privacy while enforcing
  // Gen's input consistency
  
  Bytes aux_input;
  
  uint32_t lg_k = ceil_log_base_2(Env::k());
  
  if(Env::is_root()){
    aux_input = input_prng.rand_bits(2*Env::k() + lg_k);
  }
  
  // send aux input to all processors so that the inputs are consistent
  aux_input.resize((2*Env::k()+lg_k+7)/8);
  MPI_Bcast(&aux_input[0], aux_input.size(), MPI_BYTE, 0, m_mpi_comm);
  
  // and save it to the class
  m_gen_aux_random_input = Bytes(aux_input.begin(), aux_input.end());
  
}

Bytes BetterYao5::get_gen_full_input(){
  GEN_BEGIN
  return m_private_input + get_gen_output_mask() + get_gen_input_randomness();
  GEN_END

    EVL_BEGIN
    //assert(0);
    return m_prng.rand_bits(get_gen_full_input_size());
    EVL_END
}

Bytes BetterYao5::get_gen_output_mask(){
  return m_gen_output_mask;
}

Bytes BetterYao5::get_gen_input_randomness(){
  return m_gen_aux_random_input;
}

uint32_t BetterYao5::get_gen_full_input_size(){
  // assume Gen's output mask is the same size as his input mask
  // returns the number of bits of gen's fill input size
  return 2*get_gen_inp_size() + 2*Env::k() + ceil_log_base_2(Env::k()); 
}

uint32_t BetterYao5::get_gen_inp_size(){
  std::cout << "gen input size: "<<m_gen_inp_cnt << std::endl;
  return m_gen_inp_cnt;
  
}

// assume output size = input size for sake of number of keys
uint32_t BetterYao5::get_gen_output_size(){
  return get_gen_inp_size();
}

void BetterYao5::choose_k_probe_resistant_matrix(){
  // this algorithm is given by SS13
  Prng choose_poly = Prng();
  uint32_t lg_4k,lg_4n,t,K,N;

  // initialize proper constants
  uint32_t k = Env::k();
  uint32_t n = m_private_input.size()*8; // multiply by 8, the byte width
  lg_4k = ceil_log_base_2(4*k);
  lg_4n = ceil_log_base_2(4*n);
  
  // find minimum t
  t = lg_4k > lg_4n ? lg_4k : lg_4n;
  while(1 << (t-1) > (k + (ceil_log_base_2(n)+n+k)/(t-1))){
    // while 2^(t-1) > k + (lg(n) + n + k)/(t-1)
    t--;
  }
  
  // set K and N
  K = ceil_log_base_2(n)+n+k;
  K = (K % t == 0) ? K : (K/t) + 1;
  N = K + k -1;
  
  

}


void BetterYao5::evl_generate_new_input(){
  // TODO: implement
  // transform input using k-probe-matrix
}


/**
   STEP 2: GEN COMMITS TO INPUTS
   
   first, Gen must generate the random seeds that he will use for the protocol
   and then seed his PRNGS

   then, he generates his input keys
   and sends commitments to the keys to Eval
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
  for(int j = 0; j < Env::node_load();j++){
    generate_gen_input_keys(j);
  }
  

  GEN_END
 
  // now an interactive step
  // in which Gen commits to his input keys
  commit_to_gen_input_keys();
}


/**
   generate Gen input keys and permutation bits
   using the circuit prng
 */
void BetterYao5::generate_gen_input_keys(uint32_t circuit_num){
  // this function should be executed in parallel by all processors

  //std::cout << "generate gen input keys\trank: " << Env::group_rank() << std::endl;
  
  // first, generate all of the input keys for Gen
  // (K_{0},K_{1}) <-$ {0,1}^{2k}
  m_R[circuit_num] = m_circuit_prngs[circuit_num].rand_bits(Env::k());
  m_R[circuit_num].set_ith_bit(0,1);
 
  std::vector<Bytes> random_keys;

  // generate enough input keys for Gen's entire input
  // and then compute their XOR-offsets
  generate_input_keys(m_circuit_prngs[circuit_num],
                      random_keys, //m_gen_inp_keys[circuit_num],
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

  // these lines are for debugging purposes, setting the permutation bits
  ////m_gen_inp_permutation_bits[circuit_num].resize(get_gen_full_input_size()/8+1);
  //  for(int i = 0; i < get_gen_full_input_size();i++){
  //  m_gen_inp_permutation_bits[circuit_num].set_ith_bit(i,0);
  //}

  std::cout << "gen permutation bits: " << m_gen_inp_permutation_bits[circuit_num].to_hex() << std::endl;
  
  m_gen_select_bits[circuit_num] = m_gen_inp_permutation_bits[circuit_num] ^ get_gen_full_input();
  
}


/**
   in this function, Gen commits to his input keys
   In Step 3 of the protocol, Gen must use a different prng
   to commit to his inputs, namely m_commitment_prngs
   rather than m_circuit_prngs
 */
void BetterYao5::commit_to_gen_input_keys(){
  // Gen uses m_commitment_prngs rather than one of m_circuit_prngs
  // since this randomness has to be independent of the randomness
  // used to generate the the keys and the wire labels
  
  int i,j;
  m_gen_inp_commitments.resize(Env::node_load());
  m_gen_received_input_commitments.resize(Env::node_load());

  GEN_BEGIN
    
  m_commitment_prngs.resize(Env::node_load());
  generate_random_seeds(m_commitment_seeds, Env::node_load());
  seed_prngs(m_commitment_prngs,m_commitment_seeds);

  for(j=0;j<Env::node_load();j++){
    generate_commitments(m_commitment_prngs[j],m_gen_inp_keys[j],m_gen_inp_commitments[j]);
  }

  GEN_END

    // now Gen sends Eval each set of commitments and Eval receives them  
  for(j=0;j<Env::node_load();j++){
    GEN_BEGIN
      assert(m_gen_inp_commitments[j].size() == 2 * get_gen_full_input_size());
      gen_send_evl_commitments(m_gen_inp_commitments[j]);
    GEN_END

    EVL_BEGIN
      // eval stores them in received_input_commitments
      evl_receive_gen_commitments(m_gen_received_input_commitments[j],get_gen_full_input_size()*2);
    assert(m_gen_received_input_commitments[j].size() == get_gen_full_input_size()*2);
    EVL_END
  }

}


/**
   STEP 3: AGREE ON OBJECTIVE CIRCUIT
 */

void BetterYao5::agree_on_objective_circuit(){
  std::cout << "Eval Announces K-Probe Matrix" << std::endl;
  eval_announce_k_probe_matrix();
 
  std::cout << "Collaboratively Choose 2-UHF" << std::endl;
  collaboratively_choose_2UHF();
}

void BetterYao5::eval_announce_k_probe_matrix(){
  // TODO: Implement
  // eval sends Gen the representation of the k-probe-matrix
  // as a binary array
  std::cout << "Eval Announce K Probe Unimplemented" << std::endl;
} 

void BetterYao5::collaboratively_choose_2UHF(){
  
  Bytes bufr;
  std::vector<Bytes> bufr_chunks;
  double start; // for timing

  reset_timers();

  // jointly pick a 2-UHF matrix
  // use gen's full input size (inputs + output xor pad + auxiliary randomness)
  // to come up with the proper size of a 2-UHF
  // number of bits is gen's input size * k
  bufr = flip_coins(Env::k()*((get_gen_full_input_size()+7)/8)); // only roots get the result

  // must resize buffer for all the subprocesses before sending it
  start = MPI_Wtime();
  bufr.resize(Env::k()*((get_gen_full_input_size()+7)/8));
  // comment out timers for now
  //  m_timer_evl += MPI_Wtime() - start;
  // m_timer_gen += MPI_Wtime() - start;
  
  //  start = MPI_Wtime();
  MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm);
  // m_timer_mpi += MPI_Wtime() - start;
 
  //  start = MPI_Wtime();
        
  // create m rows of length k
  m_2UHF_matrix = bufr.split(bufr.size()/Env::k());
  /*
  if(Env::group_rank() ==0){
    std::cout<<"k probe matrix: "<< std::endl;
    for(int i = 0; i < bufr.size();i++){
      std::cout << m_2UHF_matrix[i].to_hex() << std::endl;
    }
  }
  */
  // m_timer_evl += MPI_Wtime() - start;
  // m_timer_gen += MPI_Wtime() - start;

}


// interactive coin flipping protocol by ?
// Gen commits to some random coins
// Eval sends her random coins
// Gen decommits to his coins
// Eval checks the decommitment
// if checks pass, both accept the XOR
// of the two versions
// (timers commented)
Bytes BetterYao5::flip_coins(size_t len_in_bytes)
{
  double start;
  
  Bytes bufr;
  
  if (Env::is_root())
    {
      Bytes remote_coins, commitment, commit_value;
      
      //start = MPI_Wtime();
      Bytes coins = m_prng.rand_bits(len_in_bytes*8);	// Step 0: flip coins
      //m_timer_gen += MPI_Wtime() - start;
      //m_timer_evl += MPI_Wtime() - start;
      
      
      GEN_BEGIN
        //start = MPI_Wtime();
      commit_value = m_prng.rand_bits(Env::k()) + coins;	// Step 1: commit to coins
      commitment = commit_value.hash(Env::k());
      //m_timer_gen += MPI_Wtime() - start;
      
      //start = MPI_Wtime();
      GEN_SEND(commitment);
      remote_coins = GEN_RECV();     	// Step 2: receive alice's coins
      // Gen can decommit to coins only after receiving Alice's
      GEN_SEND(commit_value);		// Step 3: decommit to the coins
      //m_timer_com += MPI_Wtime() - start;
      GEN_END
        
      EVL_BEGIN
        //start = MPI_Wtime();
      commitment = EVL_RECV();	    	// Step 1: receive bob's commitment
      EVL_SEND(coins);	     	// Step 2: send coins to bob
      commit_value = EVL_RECV();
      //m_timer_com += MPI_Wtime() - start;
      
      start = MPI_Wtime();
      if (!(commit_value.hash(Env::k()) == commitment))		// Step 3: check bob's decommitment
        {
          LOG4CXX_FATAL(logger, "commitment to coins can't be properly opened");
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
      remote_coins = Bytes(commit_value.begin()+Env::k()/8, commit_value.end());
      //m_timer_evl += MPI_Wtime() - start;
      EVL_END
        
        // m_comm_sz = commitment.size() + remote_coins.size() + commit_value.size();
      
        //start = MPI_Wtime();
      
      coins ^= remote_coins;
      // combine randomnesses from both players
      
      bufr.swap(coins);
      //bufr = coins;
      
      //m_timer_evl += MPI_Wtime() - start;
      //m_timer_gen += MPI_Wtime() - start;
    }
  
  return bufr;
}


/**

   STEP 4: COMMIT TO INPUT AND OUTPUT LABELS

 */

void BetterYao5::gen_commit_to_io_labels(){
  // Gen has already generated his own input labels;
  // Now he constructs Eval's and then commits to both.
  // In the protocol, Gen must also commit to a method of constructing
  // his output labels, e.g. either generate output labels
  // immediately after output labels, or use a counter mode
  // stream cipher to generate output labels according to their 
  // gate indices (encrypt a nonce concatenated with the gate index).
  // This logic is embedded
  // TODO: fill out how we choose to do it.

  std::cout << "Gen Commits to Input/Output Labels" << std::endl;

  GEN_BEGIN
    
  m_evl_inp_keys.resize(Env::node_load());
  m_gen_inp_label_commitments.resize(Env::node_load());
  m_evl_inp_label_commitments.resize(Env::node_load());
  
  // order:
  // generate eval input keys
  // generate gen input label commitments
  // generate eval input label commitments

  m_evl_hashed_inp_keys.resize(Env::node_load());
  
  for(int i = 0; i < Env::node_load();i++){
    std::cout << "generate Eval Input Keys "<< i << std::endl;
    generate_eval_input_keys(i);
    
    std::cout << "hash eval input keys " << i << std::endl;
    m_evl_hashed_inp_keys[i].resize(m_evl_inp_keys[i].size());
    hash_eval_input_keys(m_evl_inp_keys[i],m_evl_hashed_inp_keys[i],Env::k());

    std::cout << "generate Gen's Input Label Commitments "<< i << std::endl;
    generate_gen_input_label_commitments(i);
    
    std::cout << "generate Eval's Input Label commitments "<< i << std::endl;

    //generate_eval_input_label_commitments(m_evl_hashed_inp_keys[i], m_evl_inp_label_commitments[i]);
    // generate eval input label commitments the manual way
    // using the hashed k-bit keys
    assert(m_evl_hashed_inp_keys[i].size() == 2*get_evl_inp_count());
    generate_commitments(m_circuit_prngs[i],m_evl_hashed_inp_keys[i],m_evl_inp_label_commitments[i]);
    
  }
  
  GEN_END
    
  // now, Gen commits to his and Eval's input labels
  // first his own labels, and then hers
  std:: cout << "commit to Gen input labels" << std::endl;
  commit_to_gen_input_labels();
  std:: cout << "commit to Eval input labels" << std::endl;
  commit_to_eval_input_labels();

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
    // we link ALSZ OT extensions
    // rand_key = m_circuit_prngs[circuit_num].rand_bits(Env::elm_size_in_bytes()*8);
    // rnd.from_bytes(rand_key);

    // m_evl_inp_keys[circuit_num].push_back(rand_key);

  
    rnd.random(m_circuit_prngs[circuit_num]);
    m_evl_inp_keys[circuit_num].push_back(rnd.to_bytes());
    
  
  }
  std::cout << " done generating Eval input\trank: " << Env::group_rank() << std::endl;

}


/*
  Gen (or Eval) generates Gen's input label commitments for the ith circuit
  where i is given in circuit_num
 */
void BetterYao5::generate_gen_input_label_commitments(uint32_t circuit_num){
  
  std::cout << "Generating Gen Input Label Commitments";

  assert(m_gen_inp_keys[circuit_num].size() == 2* get_gen_full_input_size());
  generate_commitments(m_circuit_prngs[circuit_num],m_gen_inp_keys[circuit_num],m_gen_inp_label_commitments[circuit_num]);
  
}


/*
  Gen (or Eval) generates Eval's input label commitments for the ith circuit
  where i is given in circuit_num
 */
// this function removed for debugging purposes while dealing with OT issues
//
//void BetterYao5::generate_eval_input_label_commitments(Prng & rng, std::vector<Bytes> & source, std::vector<Bytes> & dest){
//
//std::cout << "Generating Eval Input Label Commitments" << std::endl;
//
//assert(m_evl_inp_keys[circuit_num].size() == 2*get_evl_inp_count());
//generate_commitments(m_circuit_prngs[circuit_num],m_evl_inp_keys[circuit_num],m_evl_inp_label_commitments[circuit_num]);
//generate_commitments(rng, std::vector
//
//}


void BetterYao5::commit_to_gen_input_labels(){
  int j;

  EVL_BEGIN
  m_gen_received_label_commitments.resize(Env::node_load());
  EVL_END

  for(j=0;j<Env::node_load();j++){

  GEN_BEGIN
    assert(m_gen_inp_label_commitments[j].size() == 2*get_gen_full_input_size());
    gen_send_evl_commitments(m_gen_inp_label_commitments[j]);
  GEN_END
        
  EVL_BEGIN
    evl_receive_gen_commitments(m_gen_received_label_commitments[j],2*get_gen_full_input_size());
  EVL_END
    }
}

void BetterYao5::commit_to_eval_input_labels(){

  int j;

  EVL_BEGIN
  m_evl_received_label_commitments.resize(Env::node_load());
  EVL_END

  for(j=0;j<Env::node_load();j++){
    GEN_BEGIN
      assert(m_evl_inp_label_commitments[j].size() == 2*get_evl_inp_count());
    gen_send_evl_commitments(m_evl_inp_label_commitments[j]);
    GEN_END

    EVL_BEGIN
      
      evl_receive_gen_commitments(m_evl_received_label_commitments[j],2*get_evl_inp_count());
    assert(m_evl_received_label_commitments[j].size() == 2*get_evl_inp_count());
    EVL_END

    }
}

/**
   STEP 5: EVAL'S INPUT OTs
 */

void BetterYao5::eval_input_OT(){
  std::cout << "Eval Input OT" << std::endl;

  GEN_BEGIN
  assert(m_evl_inp_keys.size() == Env::node_load());
  
  ot_send_batch(m_evl_inp_keys);
  
  GEN_END

  EVL_BEGIN
    // TODO: make sure Eval asks for
    // her modified input keys
    // not just her original private input
    // (since at this point, Gen probably doesn't know
    // her real input keys)
    std::vector<std::vector<Bytes> > evl_receive_keys;
    assert(evl_receive_keys.size() == 0);
    assert(m_evl_received_keys.size() == 0);
    evl_receive_keys.resize(Env::node_load());
    m_evl_received_keys.resize(Env::node_load());
      
    // remember we require that the private input length be a multiple of 8
    ot_receive_batch(explode_vector(m_private_input, m_private_input.size()*8),evl_receive_keys);
    
    for(int i = 0; i < Env::node_load(); i++){
      m_evl_received_keys[i].resize(evl_receive_keys[i].size());
      hash_eval_input_keys(evl_receive_keys[i],m_evl_received_keys[i],Env::k());
    }

    fprintf(stderr,"done eval input ot: %i",Env::group_rank());
    EVL_END
    
}


void BetterYao5::hash_eval_input_keys(std::vector<Bytes> & source, std::vector<Bytes> & destination, uint32_t num_bits){
  for(int j = 0; j < source.size(); j++){
    destination[j] = source[j].hash(Env::k());
  }
}

/**
   STEP 6: CUT AND CHOOSE
 */
// Gen doesn't know which are check and which are evaluation circuits
void BetterYao5::SS13_cut_and_choose(){

  std::cout << "Cut and Choose \t rank: " << Env::group_rank() << std::endl;

  EVL_BEGIN
  evl_select_cut_and_choose_circuits();
  fprintf(stderr,"done selecting cut and choose\t rank: %i",Env::group_rank());
  EVL_END
  
  std::cout << "Special Circuit OT \t rank: " << Env::group_rank() <<  std::endl;
  special_circuit_ot();
  
  std::cout << "Transfer Masked Info W/ Prngs" << std::endl;
  transfer_evaluation_circuit_info();
  transfer_check_circuit_info();
  

}


void BetterYao5::evl_select_cut_and_choose_circuits(){

  Bytes coins;
  double start;
  
  // should be enough random bits to select choose or cut for every circuit
  Prng prng;
  
  fprintf(stderr,"evl select cut and choose circuits %i\n",Env::group_rank());

  if(Env::is_root()){
    coins = m_prng.rand_bits(Env::key_size_in_bytes());
  
    prng.seed_rand(coins);
    
    m_all_chks.assign(Env::s(),1);
    // FisherÐYates shuffle
    std::vector<uint16_t> indices(m_all_chks.size());
    for (size_t ix = 0; ix < indices.size(); ix++) { indices[ix] = ix; }
    
    // permute the indices
    // for debugging the garbling, set the 0th circuit to evaluation
    //    for (size_t ix = 0; ix < indices.size(); ix++)
    for(size_t ix=1; ix<indices.size();ix++)
      {
        int rand_ix = prng.rand_range(indices.size()-ix);
        std::swap(indices[ix], indices[ix+rand_ix]);
      }
    
    int num_of_evls;
    switch(m_all_chks.size())
      {
      case 0: case 1:
        LOG4CXX_FATAL(logger, "there aren't enough circuits for cut-and-choose");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        break;
      
      case 2: case 3:
        num_of_evls = 1;
        break;
      
      case 4:
        num_of_evls = 2;
        break;
      
      default:
        num_of_evls = m_all_chks.size()*2/5;
        break;
      }
  
    for (size_t ix = 0; ix < num_of_evls; ix++) {
      m_all_chks[indices[ix]] = 0;
    }
    m_timer_evl += MPI_Wtime() - start;

  }
  
  start = MPI_Wtime();
  // resize m_chks to prepare for scatter, which will tell whether each processor evaluates
  // check or evaluation circuits
  m_chks.resize(Env::node_load());
  m_timer_evl += MPI_Wtime() - start;
  
  start = MPI_Wtime();

  // MPI_Scatter sends a fragment of the array to each processor, which in turn allows each processor to only operate on its circuits
  MPI_Scatter(&m_all_chks[0], m_chks.size(), MPI_BYTE, &m_chks[0], m_chks.size(), MPI_BYTE, 0, m_mpi_comm);
  
  m_timer_mpi += MPI_Wtime() - start;
  
}

/**
   In this function, we perform "circuit OTs" where
   Eval selects whether a circuit is a check circuit or an evaluation circuit
   If evaluation circuit, Eval will receive the de-commitments to Gen's input keys
   and she can evaluate the circuit
   If check circuit, Eval will receive the random seeds that Gen used to generate
   the input keys (and subsequently the circuit itself) AND a nonce to seed the
   circuit's key-generating PRNG. Eval will then use this information to
   generate the circuit in lockstep with evaluating the other gates.
   To perform this circuit OT, we actually perform an OT over random seeds
   used to construct PRNGs that Gen uses to mask the information that he sends
   to Gen. If Eval has chosen the correct seed, she can decrypt the information
   Gen sends and retrieve the keys or seeds.
*/

void BetterYao5::special_circuit_ot(){
  // in this function, Gen and Evl perform a bunch of OTs in parallel,
  // Evl uses m_chks to select each seed, which will correspond to check or eval circuit
  // and each Gen processor generates m_otp_seeds for 2*his share of circuits
  // he will be oblivious to what Eval chooses
  
  GEN_BEGIN
      
  m_otp_seeds.resize(2*Env::node_load());
  m_otp_prngs.resize(2*Env::node_load());
  
  generate_random_seeds(m_otp_seeds,2*Env::node_load());
  
  seed_prngs(m_otp_prngs, m_otp_seeds);
  
  ot_send(m_otp_seeds);

  GEN_END

  EVL_BEGIN
  // m_otp_prngs will be sized appropriately by ot_receive.
  // must resize the prng vector here
  m_otp_prngs.resize(Env::node_load());
  
  ot_receive(m_chks, m_otp_seeds);
  
  // and now that Eval has the seeds, she seeds her prngs
  seed_prngs(m_otp_prngs,m_otp_seeds);
  EVL_END

}



/**
   check circuits are odd-numbered prngs
   (since we use m_chks in the OT, a 1 indicates both odd seed and check circuit)
 */

void BetterYao5::transfer_check_circuit_info(){
  
  // send the circuit generation seed used to create the input keys and permutation bits
 GEN_BEGIN
   assert(m_otp_prngs.size()/2 == Env::node_load());
   for(int i=0;i<m_otp_prngs.size()/2;i++){
    gen_send_masked_info(m_otp_prngs[2*i+1],m_circuit_seeds[i],
                         (int)(Env::elm_size_in_bytes()*8UL));
    //Env::k());
  }
  GEN_END

  EVL_BEGIN

  m_circuit_seeds.resize(m_chks.size());
  for(int i=0;i<m_chks.size();i++){
    if(m_chks[i]){ // check circuit
      evl_receive_masked_info(m_otp_prngs[i],m_circuit_seeds[i],
                              Env::elm_size_in_bytes()*8);
                              //Env::k());
    } else{
      evl_ignore_masked_info(1); // one Bytes segment sent
    }
  }
  EVL_END


    // also send a random seed for generating keys by the circuit
    // must generate it first
  GEN_BEGIN
  generate_random_seeds(m_key_generation_seeds, Env::node_load());
  for(int i = 0; i < m_otp_prngs.size()/2;i++){
    gen_send_masked_info(m_otp_prngs[2*i+1],m_key_generation_seeds[i],
                         (int)Env::elm_size_in_bytes()*8);
    //Env::k());
  }
    
  GEN_END

  EVL_BEGIN
    m_key_generation_seeds.resize(Env::node_load());
  for(int i = 0; i < m_chks.size();i++){
    if(m_chks[i]){
      evl_receive_masked_info(m_otp_prngs[i],m_key_generation_seeds[i],
                              Env::elm_size_in_bytes()*8);
                              //Env::k());
    } else{
      evl_ignore_masked_info(1); // one Bytes segment
    }
  }

  EVL_END

}


void BetterYao5::select_input_decommitments(std::vector<commitment_t> & source, std::vector<commitment_t> & dest, Bytes & permutation_bits, Bytes & input_bits){

  assert(permutation_bits.size() == input_bits.size());
  
  // there should be one key for each input
  // source should have twice the number of keys as select_bits has choices
  assert(source.size() == 2*get_gen_full_input_size());
  assert(source.size()/2 == (get_gen_full_input_size()%8 ==0 ? (permutation_bits.size()*8) : permutation_bits.size()*8 - 8 + ((source.size()/2)%8)));
  
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
  
  for(int i = 0; i < source.size()/2; i++){
    dest.push_back(source[2*i + (permutation_bits.get_ith_bit(i)
                                 ^ input_bits.get_ith_bit(i))]);
  }
}


/**
   For evaluation circuits, Gen must decommit to the inputs he chooses to send Eval
   Eval receives Gen's commitments. For check circuits, Eval simply ignores this information
 */
void BetterYao5::transfer_evaluation_circuit_info(){
  std::cout << "transfer evaluation circuit info" << std::endl;

  GEN_BEGIN 

  assert(m_otp_prngs.size() == 2*Env::node_load());

  std::vector<commitment_t> gen_decommit_inputs;
  std::vector<commitment_t> gen_decommit_labels;
  Bytes gen_input = get_gen_full_input();
  m_const_0_keys.resize(Env::node_load());
  m_const_1_keys.resize(Env::node_load());

  for(int i =0; i < Env::node_load(); i++){
    // first, select Gen's actual inputs from all those he generates
    // this uses the permutation bit to figure out if he should select the
    // first or second key allocated to each input
    select_input_decommitments(m_gen_inp_commitments[i],
                               gen_decommit_inputs,
                               m_gen_inp_permutation_bits[i],
                               gen_input);
    select_input_decommitments(m_gen_inp_label_commitments[i],
                               gen_decommit_labels,
                               m_gen_inp_permutation_bits[i],
                               gen_input);

    // gen sends Eval the decommitments to his inputs
    // this is (X1 U X2), where X1 is his set of commitments from Step 3 (Gen's Input Commitments)
    // and X2 is his set of commitments from Step 5 (Commit to I/O labels)
    gen_decommit_and_send_masked_vector(m_otp_prngs[2*i],gen_decommit_inputs);
    gen_decommit_and_send_masked_vector(m_otp_prngs[2*i],gen_decommit_labels);


    // GEN also sends Eval two keys for the constant wires
    Bytes const0,const1;
    m_const_0_keys[i] = m_circuit_prngs[i].rand_bits(Env::k());
    m_const_1_keys[i] = m_circuit_prngs[i].rand_bits(Env::k());
    m_const_0_keys[i].set_ith_bit(0,0);
    m_const_1_keys[i].set_ith_bit(0,1);
    gen_send_masked_info(m_otp_prngs[2*i],m_const_0_keys[i],Env::k());
    gen_send_masked_info(m_otp_prngs[2*i],m_const_1_keys[i],Env::k());

  }

  GEN_END

  EVL_BEGIN
    
  assert(m_chks.size() == Env::node_load());
  assert(m_chks.size() == m_gen_received_input_commitments.size());
  
  // prepare the vectors Eval will use to store Gen's decommitments
  m_cc_recv_gen_inp_commitments.resize(Env::node_load());
  m_cc_recv_gen_inp_label_commitments.resize(Env::node_load());
  
  // prepare to receive constant keys for the circuit
  m_const_0_keys.resize(Env::node_load());
  m_const_1_keys.resize(Env::node_load());


  for(int i = 0; i < Env::node_load();i++){

    assert(m_gen_received_label_commitments[i].size() == 2*get_gen_full_input_size());
    assert(m_gen_received_input_commitments[i].size() == 2*get_gen_full_input_size());

    if(!m_chks[i]){ // evaluation circuit
      // Eval receives Gen's decommitments
      // these sizes should be get_gen_full_input_size()
      
      std::cout << "received input commitments size: " 
                << m_gen_received_input_commitments[i].size()
                << "\treceived label commitments size: "
                << m_gen_received_label_commitments[i].size()
                << "\trank: " << Env::group_rank() << "\tcircuit: " << i << std::endl;
      
      evl_receive_masked_vector(m_otp_prngs[i],m_cc_recv_gen_inp_commitments[i],2*Env::k(),m_gen_received_input_commitments[i].size()/2);
      evl_receive_masked_vector(m_otp_prngs[i],m_cc_recv_gen_inp_label_commitments[i],2*Env::k(),m_gen_received_label_commitments[i].size()/2);

      assert(m_cc_recv_gen_inp_commitments[i].size() == get_gen_full_input_size());
      assert(m_cc_recv_gen_inp_label_commitments[i].size() == get_gen_full_input_size());

      evl_receive_masked_info(m_otp_prngs[i],m_const_0_keys[i],Env::k());
      evl_receive_masked_info(m_otp_prngs[i],m_const_1_keys[i],Env::k());


    } else {
      // these sizes should be get_gen_full_input_size()
      evl_ignore_masked_info(m_gen_received_input_commitments[i].size()/2); // decommitted inputs
      evl_ignore_masked_info(m_gen_received_label_commitments[i].size()/2); // decommitted labels

      // ignore the constant keys
      evl_ignore_masked_info(1);
      evl_ignore_masked_info(1);
        
    }
  }

  EVL_END   

}

/**
   Gen decommits to a set of commitments and sends the information to Eval
   he provides the Prng that generates the masks and the vector of commitments.
 */
void BetterYao5::gen_decommit_and_send_masked_vector(Prng & mask_generator, std::vector<commitment_t> & vec){ //, uint32_t chunk_size){
  for(int i = 0; i < vec.size();i++){
    assert(decommit(vec[i]).size() > 0);
    gen_send_masked_info(mask_generator,decommit(vec[i]),decommit(vec[i]).size()*8);
  }
}

/**
   gen sends a piece of information to Eval
   he must know how long it is (in bits) to compute
   an XOR mask with his Prng
 */
void BetterYao5::gen_send_masked_info(Prng & mask_generator, Bytes info, uint32_t chunk_size){
  assert(chunk_size == info.size()*8);
  Bytes send = info ^ mask_generator.rand_bits(chunk_size);
  GEN_SEND(send);

}

/**
   Eval receives an entire vector of information from Gen
   and decrypts with bits from the Prng she provides
   she must also know the size of the "chunks" or vector pieces she will receive
   as well as the vector length
 */
void BetterYao5::evl_receive_masked_vector(Prng & mask_generator, std::vector<Bytes> & destination, uint32_t chunk_size, uint32_t len){
  destination.resize(len);
  for(int i = 0; i < len; i++){
    evl_receive_masked_info(mask_generator, destination[i], chunk_size);
  }
}

/**
   Eval receives information from Gen and decrypts it with bits from the prgn she provides
 */
void BetterYao5::evl_receive_masked_info(Prng & mask_generator, Bytes & destination, uint32_t chunk_size){
  destination = EVL_RECV();
  destination ^= mask_generator.rand_bits(chunk_size);
  assert(destination.size()*8 == chunk_size);
}

void BetterYao5::evl_ignore_masked_info(uint32_t len){
  Bytes recv;
  for(int i = 0; i < len; i++){
    recv = EVL_RECV();
  }
  // and then do nothing because Eval can't decrypt this information anyway
}

/**
   STEP 7: GARBLE AND CHECK
*/


void BetterYao5::garble_and_check_circuits(){
  std::cout << "garble and check circuits" << std::endl;
  
  
  //MPI_Barrier(m_mpi_comm);

  EVL_BEGIN

  m_circuit_prngs.resize(Env::node_load());
  m_gen_inp_keys.resize(Env::node_load());
  m_evl_inp_keys.resize(Env::node_load());
  m_evl_hashed_inp_keys.resize(Env::node_load());
  m_gen_inp_label_commitments.resize(Env::node_load());
  m_evl_inp_label_commitments.resize(Env::node_load());
  m_gen_inp_permutation_bits.resize(Env::node_load());
  m_R.resize(Env::node_load());
  m_gen_select_bits.resize(Env::node_load());

  for(int i = 0; i < Env::node_load();i++){
    if(m_chks[i]){ // this is a check circuit
      evl_regenerate_circuits(i);
      evl_check_commitment_regeneration(i);
    } else{
      evl_check_garbled_circuit_commitments(i);
      evl_set_inp_keys(i);
    }
  }
  
  // prepare for garbling checks
  m_2UHF_hashes.resize(Env::node_load());

  EVL_END
    // now that checks are done, we can garble!
    
    initialize_circuits();
    evaluate_circuits();
  

}


void BetterYao5::evl_inputs_transform(std::vector<Bytes> &source, std::vector<Bytes> &dest){
  G a;
  for(int i = 0; i < source.size();i++){
    a.from_bytes(source[i]);
    dest.push_back(a.to_bytes());
  }
}


// this function only called on check circuits
// (circuit_num is the index of a check circuit
void BetterYao5::evl_regenerate_circuits(uint32_t circuit_num){
  // first, Evl seeds the circuit prngs for which she has received seeds
  //TODO: FILL IN GENERATING CONSTANT KEYS

  m_circuit_prngs[circuit_num].seed_rand(m_circuit_seeds[circuit_num]);

  // next, Evl must regenerate the commitments that Gen sent her
  
  // order: of generating things with circuit seed rho:
  // 1. generate Gen input keys and permutation bits
  // 2. generate eval input keys
  // 3. generate gen input label commitments
  // 4. generate eval input label commitments

  std::cout << "generate Gen Input Keys" << std::endl;
  generate_gen_input_keys(circuit_num);
  
  std::cout << "generate Eval Input Keys" << std::endl;
  generate_eval_input_keys(circuit_num);
  
  // these have to be transformed by changing them to and from 
  // G elements a couple times, to simulate what would happen in OT.
  std::vector<Bytes> evl_transform_inputs;
  evl_inputs_transform(m_evl_inp_keys[circuit_num],evl_transform_inputs);

  std::cout << "Hash Eval Input Keys" << std::endl;
  // these have to be hashed to become the correct length for an input key
  assert(m_evl_hashed_inp_keys[circuit_num].size() == 0);
  
  m_evl_hashed_inp_keys[circuit_num].resize(m_evl_inp_keys[circuit_num].size());
  hash_eval_input_keys(evl_transform_inputs,m_evl_hashed_inp_keys[circuit_num], Env::k());
  //hash_eval_input_keys(m_evl_inp_keys[circuit_num],m_evl_hashed_inp_keys[circuit_num], Env::k());
  assert(m_evl_inp_keys[circuit_num].size() == m_evl_hashed_inp_keys[circuit_num].size());

  std::cout << "generate Gen's Input Label Commitments" << std::endl;
  generate_gen_input_label_commitments(circuit_num);

  std::cout << "generate Eval's Input Label commitments" << std::endl;
  //  generate_eval_input_label_commitments(circuit_num);
  
  generate_commitments(m_circuit_prngs[circuit_num],m_evl_hashed_inp_keys[circuit_num],m_evl_inp_label_commitments[circuit_num]);
  assert(m_evl_hashed_inp_keys[circuit_num].size() == m_evl_inp_label_commitments[circuit_num].size());
  
}



void BetterYao5::evl_check_commitment_regeneration(uint32_t circuit_num){
  bool verify = true;

  // check consistency of the gen input label commitments
  assert(m_gen_received_label_commitments[circuit_num].size() == m_gen_inp_label_commitments[circuit_num].size());
  verify &= check_received_commitments_vs_generated(m_gen_received_label_commitments[circuit_num],m_gen_inp_label_commitments[circuit_num],0);

  if(!verify)
    fprintf(stderr,"check of gen received input labels/commitments failed\n");
  else 
    fprintf(stderr,"check of gen received input labels/commitments passed\n");

  // TODO: will need to come back to eval's input label commitments
  // currently failing because of the changes that happen to the representation
  // of random strings converted to G elements
  
  // check consistency of the eval input label commitments
  // assert(m_evl_received_label_commitments[circuit_num].size() == m_evl_inp_label_commitments[circuit_num].size());
  // verify &= check_received_commitments_vs_generated(m_evl_received_label_commitments[circuit_num], m_evl_inp_label_commitments[circuit_num],1);

  //if(!verify)
  //  fprintf(stderr,"check of eval received input labels/commitments failed\n");
  //else
  //  fprintf(stderr,"check of eval received input labels/commitments passed\n");
  
  if(!verify){
    std::cout << "no verify" << std::endl;
    LOG4CXX_FATAL(logger, "Eval's Check of Gen's Check Circuits Failed");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  } else{
    std::cout << "verified" << std::endl;
  }
  
}

void BetterYao5::evl_check_garbled_circuit_commitments(uint32_t circuit_num){
  // the cc_received vectors should be half the size of the previously received vectors
  // since the cc_received vectors should only have the input keys Gen uses
  assert(m_cc_recv_gen_inp_commitments[circuit_num].size() == get_gen_full_input_size());
  assert(m_cc_recv_gen_inp_label_commitments[circuit_num].size() == get_gen_full_input_size());
  assert(m_cc_recv_gen_inp_commitments[circuit_num].size() == m_gen_received_input_commitments[circuit_num].size()/2);
  assert(m_cc_recv_gen_inp_label_commitments[circuit_num].size() == m_gen_received_label_commitments[circuit_num].size()/2);

  bool verify = true;
  for(int i = 0; i < m_cc_recv_gen_inp_commitments[circuit_num].size(); i++){
    
    // verify that the decommitment we just received is consistent with one that we got in Gen's original input commitments
    verify &= (m_cc_recv_gen_inp_commitments[circuit_num][i].hash(Env::k()) == m_gen_received_input_commitments[circuit_num][2*i] || m_cc_recv_gen_inp_commitments[circuit_num][i].hash(Env::k()) == m_gen_received_input_commitments[circuit_num][2*i+1]);

    // verify that the decommitment we just received is consistent with one that we got in Gen's io commitments
    verify &= (m_cc_recv_gen_inp_label_commitments[circuit_num][i].hash(Env::k()) == m_gen_received_label_commitments[circuit_num][2*i] || m_cc_recv_gen_inp_label_commitments[circuit_num][i].hash(Env::k()) == m_gen_received_label_commitments[circuit_num][2*i+1]);
    
    // verify that the two decommitments we just received are consistent with each other
    verify &= (reconstruct_commitment(m_cc_recv_gen_inp_commitments[circuit_num][i]).msg == reconstruct_commitment(m_cc_recv_gen_inp_label_commitments[circuit_num][i]).msg);

  }

  if(!verify){
    std::cout << "no verify" << std::endl;
    LOG4CXX_FATAL(logger, "Eval's Check of Gen's Evaluation Circuits Failed");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  std::cout << "done verify" << std::endl;
}


// set eval's input keys for use in circuit evaluation
void BetterYao5::evl_set_inp_keys(uint32_t circuit_num){
  for(int i = 0; i < get_gen_full_input_size();i++){
    m_gen_inp_keys[circuit_num].push_back(reconstruct_commitment(m_cc_recv_gen_inp_commitments[circuit_num][i]).msg);
  }
  
  // the following loop is not strictly necessary but makes the code cleaner
  // we coud just use m_evl_received_keys
  for(int i = 0; i < get_evl_inp_count();i++){
    m_evl_inp_keys[circuit_num].push_back(m_evl_received_keys[circuit_num][i]);
  }

}

bool BetterYao5::check_received_commitments_vs_generated(std::vector<Bytes> & received, std::vector<commitment_t> & generated, uint32_t print_t){
  assert(received.size() == generated.size());

  std::cout << "checking received vs generated circuits " << std::endl;

  bool verify = true;
  for(int i = 0; i < received.size();i++){
    verify &= (received[i] == decommit(generated[i]).hash(Env::k()));

    // debugging output (flagged)
    if(print_t){
      fprintf(stderr,"received commitment: %s\tregenerated commitment: %s\thashed decommitment: %s\trank: %i\n",
              received[i].to_hex().c_str(),
              decommit(generated[i]).to_hex().c_str(),
              decommit(generated[i]).hash(Env::k()).to_hex().c_str(),
              Env::group_rank());
      }
  }
  return verify;
}


void BetterYao5::initialize_circuits(){
  MPI_Barrier(m_mpi_comm);
  GEN_BEGIN
  if(Env::group_rank() == 0){
      // begin with one, just to make debugging easier
    for(int ix = 0; ix < m_gcs.size();ix++){
      
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
      
      fprintf(stderr,"Gen Input: \t%s\nGen Full Input: \t%s\nGen Permutation Bits: \t%s\nGen select bits:\t%s\n",
              m_private_input.to_hex().c_str(),
              get_gen_full_input().to_hex().c_str(),
              m_gen_inp_permutation_bits[ix].to_hex().c_str(),
              m_gen_select_bits[ix].to_hex().c_str());
      
    }
  }
  GEN_END
  
  EVL_BEGIN

  if(Env::group_rank() == 0){
    // begin with one, just to make debugging easier
    for(int ix = 0; ix < m_gcs.size();ix++){
      
      if(m_chks[ix]){
        // if check circuit, then we will generate it
        // TODO: implement this
      }
      else if(!m_chks[ix]){
        
        // if evaluation circuit, we will evaluate it
        m_gcs[ix].init_Evaluation_Circuit(&m_gen_inp_keys[ix],// gen keys
                                          &m_evl_received_keys[ix],//evl keys
                                          get_gen_inp_size(),// gen inp size
                                          m_private_input, // private input
                                          m_const_0_keys[ix], //output keys
                                          m_const_1_keys[ix]);
        
        m_gcs[ix].m_st = 
          load_pcf_file(Env::pcf_file(), m_gcs[ix].get_Const_Wire(0), m_gcs[ix].get_Const_Wire(1), copy_key);
        m_gcs[ix].m_st->alice_in_size = get_gen_full_input_size();
        m_gcs[ix].m_st->bob_in_size = get_evl_inp_count();
        
        set_external_circuit(m_gcs[ix].m_st, &m_gcs[ix]);
        
        m_gcs[ix].set_Evl_Circuit_Functions();
        
      }
    }
  }

  EVL_END


}

void BetterYao5::generate_2UHF(){
  if(Env::group_rank()==0){
    for(int ix = 0; ix < Env::node_load();ix++){
      m_gcs[ix].generate_Gen_Inp_Hash(m_2UHF_matrix);
      GEN_SEND(m_gcs[ix].get_hash_out());
    }
  }
}

void BetterYao5::evaluate_2UHF(){
  if(Env::group_rank()==0){
    for(int ix = 0; ix < Env::node_load();ix++){
      m_gcs[ix].evaluate_Gen_Inp_Hash(m_2UHF_matrix);
      Bytes hash_parity = EVL_RECV();
      m_2UHF_hashes[ix] = hash_parity ^ m_gcs[ix].get_hash_out();
    }
  }
}

void BetterYao5::evaluate_circuits(){
  
  GEN_BEGIN
    if(Env::group_rank()==0){
      for(int ix = 0; ix < Env::node_load();ix++){
        Bytes bufr;
        while(get_next_gate(m_gcs[ix].m_st)){
          bufr = m_gcs[ix].get_garbling_bufr();
          GEN_SEND(bufr);
          m_gcs[ix].clear_garbling_bufr();
        }
        GEN_SEND(Bytes(0)); // redundant value to prevent Evl from hanging
        
      }
    }

  GEN_END
    
    EVL_BEGIN

    if(Env::group_rank()==0){
      for(int ix = 0; ix < Env::node_load(); ix++){ 
        Bytes bufr;
        do {
          bufr = EVL_RECV();
          m_gcs[ix].set_garbling_bufr(bufr);
          
        } while(get_next_gate(m_gcs[ix].m_st));
      }
    }
    EVL_END
}


/**
   STEP 8: RETRIEVE OUTPUTS
 */

void BetterYao5::retrieve_outputs(){
  // this function currently outputs in the clear
  // TODO: gen output authenticity and masking gen's outputs

  GEN_BEGIN

  assert(m_gcs.size()>0);  

  for(int i = 0; i < m_gcs.size();i++){
    m_gcs[i].trim_output_buffers();

    Bytes alice_out_parity = m_gcs[i].get_alice_out();
    GEN_SEND(alice_out_parity);
    
    Bytes bob_out_parity = m_gcs[i].get_bob_out();
    Bytes bob_evaluated_out = GEN_RECV();
    Bytes bob_out = bob_out_parity ^ bob_evaluated_out;
    //GEN_SEND(bob_out_parity);
    std::cout << "bob out (masked): " << bob_evaluated_out.to_hex() << std::endl;
    std::cout << "bob out (parity): " << bob_out_parity.to_hex() << std::endl;
    //    bob_out = bob_out ^ bob_out_parity;
    std::cout << "bob out  (final):" << bob_out.to_hex() << std::endl;

    }
    
  GEN_END

  EVL_BEGIN
    std::cout << "Retrieve Outputs" << std::endl;
  
  
  assert(m_gcs.size()>0);
  for(int i = 0; i < m_gcs.size();i++){
    m_gcs[i].trim_output_buffers();

    Bytes alice_out = m_gcs[i].get_alice_out();
    Bytes alice_out_parity = EVL_RECV();
    
    Bytes bob_out = m_gcs[i].get_bob_out();
    EVL_SEND(bob_out);
    //    Bytes bob_out_parity = EVL_RECV();
   
    std::cout << "alice out (masked): " << alice_out.to_hex() << std::endl;
    std::cout << "alice out (parity): " << alice_out_parity.to_hex() << std::endl;
    alice_out = alice_out ^ alice_out_parity;
    std::cout << "alice out  (final): " << alice_out.to_hex() << std::endl;
   
  }


  EVL_END

}

void BetterYao5::gen_output_auth_proof(){
  
  // step 1: Gen chooses random nonce r in {0,1}^k and encrypts it with each gate key u_a^j
  // then sends Evl the encryptions of r under each key

  // step 2: Evl commits to the decryption of the nonce with the output key for this bit from circuit m

  // step 3: Gen receives Evl's commitment and decommits to his output keys (which technically were shared inputs)

  // step 4: Eval checks Gen's decommitments:
  //         making sure the output keys decommit properly, and
  //         all of the decryptions using the keys produce the same nonce
  // if so Evl accepts, and then she decommits to r'

  // step 5: Gen checks that com(r') properly decommits to r. If so, he accepts.
  
}



/**
   UTILITY FUNCTIONS
 */


uint32_t ceil_log_base_2(uint32_t k){
  uint32_t tmp = 2;
  uint32_t lg_k = 1;

  while(tmp<k){
    lg_k++;
    tmp = tmp << 1;
  }
  return lg_k;
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
    if(Env::group_rank() == 0)
      for(int i = 0; i < sender_inputs[j].size();i++){
         fprintf(stderr,"ot send input (%i): %s\n",i,sender_inputs[j][i].to_hex().c_str());
      }
    ot_send_random(sender_inputs[j]);
  }
}

void BetterYao5::ot_receive_batch(Bytes selection_bits, std::vector<std::vector<Bytes> > & results_container){
  ot_receive_init();
  for(int j = 0; j < results_container.size();j++){
    ot_receive_random(selection_bits, results_container[j]);
    if(Env::group_rank() == 0)
      for(int i = 0; i < results_container[j].size();i++){
         fprintf(stderr,"ot receive input (%i): %s\n",i,results_container[j][i].to_hex().c_str());
     }
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
  
  //  m_ot_out.clear();
  //  m_ot_out.reserve(2*m_ot_bit_cnt); // the receiver only uses half of it

  
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

