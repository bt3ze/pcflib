#include "BetterYao5.h"
// #include "garbled_circuit.h"
#include <algorithm>

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
      //m_gcs[ix] = GarbledMal();
      m_gcs[ix] = GarbledCircuit();
    }

  // m_gen_inp_hash.resize(Env::node_load());
  // m_gen_inp_masks.resize(Env::node_load());
  // m_gen_inp_decom.resize(Env::node_load());
  // m_gen_inp_decom.resize(m_gen_inp_cnt);
  
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
  Prng prng = Prng();

  seeds.clear();

  for(int i=0;i<num_seeds;i++){
    // these random elements will be used for OT,
    // so we generate them as random group elements

    /*
    // this version does not work
    rand = prng.rand_bits(Env::elm_size_in_bytes()*8);
    std::cout << "prng rand element: " << rand.to_hex() << std::endl;
    seeds.push_back(rand);
    */

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
*/
void BetterYao5::generate_input_keys(Prng & prng, std::vector<Bytes> & dest, uint32_t num_keys){
  Bytes key;
  int j;
  for(j=0;j<num_keys;j++){
    // generate the key (so that it is the right size)
    key = prng.rand_bits(Env::k());
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
  return m_private_input + get_gen_output_mask() + get_gen_input_randomness();
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
  return m_gen_inp_cnt + m_gen_inp_cnt + 2*Env::k() + ceil_log_base_2(Env::k()); 
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
  std::cout << "generate random circuit seeds" << std::endl;
  m_circuit_prngs.resize(Env::node_load());
  generate_random_seeds(m_circuit_seeds, Env::node_load());
  seed_prngs(m_circuit_prngs,m_circuit_seeds);
  
  std::cout << "m_circuit_prngs seeded\t rank: " << Env::group_rank() << std::endl; 

  // next, generate Gen's input keys
  m_gen_inp_keys.resize(Env::node_load());
  m_gen_inp_permutation_bits.resize(Env::node_load());
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

  std::cout << "generate gen input keys\trank: " << Env::group_rank() << std::endl;

  // first, generate all of the input keys for Gen
  // (K_{0},K_{1}) <-$ {0,1}^{2k}
  generate_input_keys(m_circuit_prngs[circuit_num],
                      m_gen_inp_keys[circuit_num],
                      get_gen_full_input_size()*2);
  
  
  // next generate permutation bits for Gen's input keys
  // { pi_{i} } <-$ {0,1} for i in (0, gen_inputs ]  
  m_gen_inp_permutation_bits[circuit_num] = m_circuit_prngs[circuit_num].rand_bits(get_gen_full_input_size());
  
  // now set Gen's permutation bits
  // permutation bit 1 --> will be sent as the second in key pair to Eval
  for(int i = 0; i < 2*get_gen_full_input_size();i+=2){
    // this permutes the keys based on the permutation bit
    if(m_gen_inp_permutation_bits[circuit_num].get_ith_bit(i/2)==1){
      // swap the two keys, then set their bits so that they look
      // like nothing's up (which lets us commit to them in this order)
      // but the secret permutation bits generated above hold the key to their true values.
      std::swap(m_gen_inp_keys[circuit_num][i],m_gen_inp_keys[circuit_num][i+1]);
      m_gen_inp_keys[circuit_num][i].set_ith_bit(0,0);
      m_gen_inp_keys[circuit_num][i+1].set_ith_bit(0,1);
    }
  }
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
  //  m_timer_evl += MPI_Wtime() - start;
  // m_timer_gen += MPI_Wtime() - start;
  
  //  start = MPI_Wtime();
  MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm);
  // m_timer_mpi += MPI_Wtime() - start;
 
  //  start = MPI_Wtime();
        
  // create m rows of length k
  m_2UHF_matrix = bufr.split(bufr.size()/Env::k());
  
  // m_timer_evl += MPI_Wtime() - start;
  // m_timer_gen += MPI_Wtime() - start;

  /*
  if(Env::is_root()){
    std::cout << "agree on 2-UHF: " << std::endl;
    int i;
    for(i=0;i<m_2UHF_matrix.size();i++){
      std::cout << m_2UHF_matrix[i].to_hex() << std::endl;
    }
  }
  */
}


// interactive coin flipping protocol by ?
// Gen commits to some random coins
// Eval sends her random coins
// Gen decommits to his coins
// Eval checks the decommitment
// if checks pass, both accept the XOR
// of the two versions
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

  for(int i = 0; i < Env::node_load();i++){
    std::cout << "generate Eval Input Keys "<< i << std::endl;
    generate_eval_input_keys(i);
    
    std::cout << "generate Gen's Input Label Commitments "<< i << std::endl;
    generate_gen_input_label_commitments(i);
    
    std::cout << "generate Eval's Input Label commitments "<< i << std::endl;
    generate_eval_input_label_commitments(i);

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
 */
void BetterYao5::generate_eval_input_keys(uint32_t circuit_num){
    
  generate_input_keys(m_circuit_prngs[circuit_num],m_evl_inp_keys[circuit_num],get_evl_inp_count()*2);
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
void BetterYao5::generate_eval_input_label_commitments(uint32_t circuit_num){
  
  std::cout << "Generating Eval Input Label Commitments" << std::endl;
  
assert(m_evl_inp_keys[circuit_num].size() == 2*get_evl_inp_count());
  generate_commitments(m_circuit_prngs[circuit_num],m_evl_inp_keys[circuit_num],m_evl_inp_label_commitments[circuit_num]);
  

}

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
  
  for(int i = 0; i < m_evl_inp_keys.size(); i++){
    std::cout << "eval input OT: " << m_evl_inp_keys[i].size() << "\trank: " << Env::group_rank() << std::endl;
  }
  ot_send_batch(m_evl_inp_keys);

  GEN_END

  EVL_BEGIN
    // TODO: make sure Eval asks for
    // her modified input keys
    // not just her original private input
    // (since at this point, Gen probably doesn't know
    // her real input keys)
    assert(m_evl_received_keys.size() == 0);
    m_evl_received_keys.resize(Env::node_load());

    std::cout << "eval input OT: " << explode_vector(m_private_input,m_private_input.size()*8).size() << "\trank: " << Env::group_rank() << std::endl;
      
    // remember we require that the private input length be a multiple of 8
    ot_receive_batch(explode_vector(m_private_input, m_private_input.size()*8),m_evl_received_keys);
    
  EVL_END

}


/**
   STEP 6: CUT AND CHOOSE
 */
// Gen doesn't know which are check and which are evaluation circuits
void BetterYao5::SS13_cut_and_choose(){

  std::cout << "Cut and Choose" << std::endl;

  EVL_BEGIN
  evl_select_cut_and_choose_circuits();
  EVL_END
  
  std::cout << "Special Circuit OT" << std::endl;
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
  
  if(Env::is_root()){
    coins = m_prng.rand_bits(Env::key_size_in_bytes());
  
    std::cout << "cut-and-choose coins: " << coins.to_hex() << std::endl; 

    prng.seed_rand(coins);
    
    m_all_chks.assign(Env::s(),1);
    // FisherÐYates shuffle
    std::vector<uint16_t> indices(m_all_chks.size());
    for (size_t ix = 0; ix < indices.size(); ix++) { indices[ix] = ix; }
    
    // permute the indices
    for (size_t ix = 0; ix < indices.size(); ix++)
      {
        int rand_ix = prng.rand_range(indices.size()-ix);
        std::swap(indices[ix], indices[ix+rand_ix]);
      }
    
    int num_of_evls;
    std::cout << "m_all_chks.size: " << m_all_chks.size() << std::endl;
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
   the input keys (and subsequently the circuit itself). Eval will then use this
   information to generate the circuit in lockstep with evaluating the other gates.
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
  
  std::cout << "circuit OT: " << m_otp_seeds.size() << "\trank: " << Env::group_rank() << std::endl;
  ot_send(m_otp_seeds);

  GEN_END

  EVL_BEGIN
  // m_otp_prngs will be sized appropriately by ot_receive.
  // must resize the prng vector here
  m_otp_prngs.resize(Env::node_load());
  
  std::cout << "circuit OT: " << m_chks.size() << "\trank: " << Env::group_rank() << std::endl;
  
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
  GEN_BEGIN
  for(int i=0;i<m_otp_prngs.size()/2;i++){
    gen_send_masked_info(m_otp_prngs[2*i+1],m_circuit_seeds[i],(int)(Env::elm_size_in_bytes()*8UL));
    std::cout << "m_circuit_seeds " << m_circuit_seeds[i].to_hex()
              << "\trank: " << Env::group_rank() << std::endl;
  }
  GEN_END

  EVL_BEGIN

    // std::cout << "m_circuit_seeds size: " << m_circuit_seeds.size() << std::endl;
  m_circuit_seeds.resize(m_chks.size());
  for(int i=0;i<m_chks.size();i++){
    if(m_chks[i]){ // check circuit
      evl_receive_masked_info(m_otp_prngs[i],m_circuit_seeds[i],Env::elm_size_in_bytes()*8);
      
      std::cout << "m_circuit_seeds " << m_circuit_seeds[i].to_hex()
                << "\trank: " << Env::group_rank() << std::endl;
    } else{
      evl_ignore_masked_info(1); // one Bytes segment sent
    }
  }
  EVL_END

}


void BetterYao5::select_input_decommitments(std::vector<commitment_t> & source, std::vector<commitment_t> & dest, Bytes & permutation_bits, Bytes input_bits){
  // source should have twice the number of keys as select_bits has choices

  std::cout << "select input decommitments\trank: " << Env::group_rank() << std::endl;
  assert(permutation_bits.size() == input_bits.size());

  std::cout << " source/2: " << source.size()/2 << "\tperm: " << permutation_bits.size()*8 
            << "\trank: " << Env::group_rank() <<"\t: "
            << permutation_bits.size()*8 - 8 + (source.size()/2)%8 << std::endl;
  
  // there should be one key for each input
  assert(source.size() == 2*get_gen_full_input_size());
  assert(source.size()/2 == (get_gen_full_input_size()%8 ==0 ? (permutation_bits.size()*8) : permutation_bits.size()*8 - 8 + ((source.size()/2)%8)));
  
  // now, select the proper input key
  // which is given by m_gen_inp_permutation_bits[i] XOR m_gen_inp[i]
  
  std::cout << "index: ";
  for(int i = 0; i < source.size()/2; i++){
    std::cout << i << "\t";
    if(permutation_bits.get_ith_bit(i)==0){
      
      dest.push_back(source[2*i+  (permutation_bits.get_ith_bit(i)^input_bits.get_ith_bit(i))]);
    } else{
      dest.push_back(source[2*i+(1^(permutation_bits.get_ith_bit(i)^input_bits.get_ith_bit(i)))]);
    }
  }
  std::cout << "done selecting decommitments" << std::endl;
}

/**
   For evaluation circuits, Gen must decommit to the inputs he chooses to send Eval
   TODO: make sure he only decommits to half of the keys
   Eval receives Gen's commitments. For check circuits, Eval simply ignores this information
   TODO: fix the sizes that Eval expects (only half the number of keys)
 */
void BetterYao5::transfer_evaluation_circuit_info(){
  std::cout << "transfer evaluation circuit info" << std::endl;

  GEN_BEGIN 

  assert(m_otp_prngs.size() == 2*Env::node_load());

  std::vector<commitment_t> gen_decommit_inputs;
  std::vector<commitment_t> gen_decommit_labels;

  for(int i =0; i < Env::node_load(); i++){
    // TODO: fix chunk size
    // first, select Gen's actual inputs from all those he generates
    // this uses the permutation bit to figure out if he should select the
    // first or second key allocated to each input
    select_input_decommitments(m_gen_inp_commitments[i],
                               gen_decommit_inputs,
                               m_gen_inp_permutation_bits[i],
                               get_gen_full_input());
    select_input_decommitments(m_gen_inp_label_commitments[i],
                               gen_decommit_labels,
                               m_gen_inp_permutation_bits[i],
                               get_gen_full_input());
    

    // gen sends Eval the decommitments to his inputs
    // this is (X1 U X2), where X1 is his set of commitments from Step 3 (Gen's Input Commitments)
    // and X2 is his set of commitments from Step 5 (Commit to I/O labels)
    gen_decommit_and_send_masked_vector(m_otp_prngs[2*i],gen_decommit_labels);
    gen_decommit_and_send_masked_vector(m_otp_prngs[2*i],gen_decommit_labels);
  }

  GEN_END

  EVL_BEGIN
    
  assert(m_chks.size() == Env::node_load());
  assert(m_chks.size() == m_gen_received_input_commitments.size());
  
  // prepare the vectors Eval will use to store Gen's decommitments
  m_cc_recv_gen_inp_commitments.resize(Env::node_load());
  m_cc_recv_gen_inp_label_commitments.resize(Env::node_load());

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

    } else {
      // these sizes should be get_gen_full_input_size()
      evl_ignore_masked_info(m_gen_received_input_commitments[i].size()); // decommitted inputs
      evl_ignore_masked_info(m_gen_received_label_commitments[i].size()); // decommitted labels
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
    //std::cout << " eval receive vector index: " << i << std::endl;
    evl_receive_masked_info(mask_generator, destination[i], chunk_size);
  }
}

/**
   Eval receives information from Gen and decrypts it with bits from the prgn she provides
 */
void BetterYao5::evl_receive_masked_info(Prng & mask_generator, Bytes & destination, uint32_t chunk_size){
  //  Bytes recv;
  destination = EVL_RECV();
  //  std::cout << "dest size: " << destination.size()*8 << "\tchunk size: " << chunk_size << std::endl;
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
  m_gen_inp_label_commitments.resize(Env::node_load());
  m_evl_inp_label_commitments.resize(Env::node_load());
  m_gen_inp_permutation_bits.resize(Env::node_load());

  for(int i = 0; i < Env::node_load();i++){
    if(m_chks[i]){ // this is a check circuit
      evl_regenerate_circuits(i);
      evl_check_commitment_regeneration(i);
    } else{
      evl_check_garbled_circuit_commitments(i);
    }
  }
  
  EVL_END
    // now that checks are done, we can garble!
    

}

// this function only called on check circuits
// (circuit_num is the index of a check circuit
void BetterYao5::evl_regenerate_circuits(uint32_t circuit_num){
  // first, Evl seeds the circuit prngs for which she has received seeds
  
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

  std::cout << "generate Gen's Input Label Commitments" << std::endl;
  generate_gen_input_label_commitments(circuit_num);

  std::cout << "generate Eval's Input Label commitments" << std::endl;
  generate_eval_input_label_commitments(circuit_num);
  
}



void BetterYao5::evl_check_commitment_regeneration(uint32_t circuit_num){
  bool verify = true;

  // check consistency of the gen input label commitments
  assert(m_gen_received_label_commitments[circuit_num].size() == m_gen_inp_label_commitments[circuit_num].size());
  verify &= check_received_commitments_vs_generated(m_gen_received_label_commitments[circuit_num],m_gen_inp_label_commitments[circuit_num]);
      
  // check consistency of the eval input label commitments
  std::cout << "number of evl label commitments: "
            << m_evl_received_label_commitments[circuit_num].size()
            << "\t number of generated evl input labels: " 
            << m_evl_inp_label_commitments[circuit_num].size()
            << "\trank: " << Env::group_rank()
            << std::endl;
  
  assert(m_evl_received_label_commitments[circuit_num].size() == m_evl_inp_label_commitments[circuit_num].size());
  verify &= check_received_commitments_vs_generated(m_evl_received_label_commitments[circuit_num], m_evl_inp_label_commitments[circuit_num]);
  
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
  // TODO: fix these assertions (divide the right hand sides by 2)
  // will also need to fix the verify statements below to ensure
  // that each cc_received entry is consistent with one of two earlier received commitments
  assert(m_cc_recv_gen_inp_commitments[circuit_num].size() == get_gen_full_input_size());
  assert(m_cc_recv_gen_inp_label_commitments[circuit_num].size() == get_gen_full_input_size());
  assert(m_cc_recv_gen_inp_commitments[circuit_num].size() == m_gen_received_input_commitments[circuit_num].size()/2);
  assert(m_cc_recv_gen_inp_label_commitments[circuit_num].size() == m_gen_received_label_commitments[circuit_num].size()/2);

  std::cout << "checking garbled circuit commitments" << std::endl;
  
  bool verify = true;
  std::cout << "gen inp commitment size: " << m_cc_recv_gen_inp_commitments[circuit_num].size() << std::endl;
  for(int i = 0; i < m_cc_recv_gen_inp_commitments[circuit_num].size(); i++){
    //verify &= (m_cc_recv_gen_inp_commitments[circuit_num][i].hash(Env::k()) == m_gen_received_input_commitments[circuit_num][2*i] || m_cc_recv_gen_inp_commitments[circuit_num][i].hash(Env::k()) == m_gen_received_input_commitments[circuit_num][2*i+1]);
    //verify &= (m_cc_recv_gen_inp_label_commitments[circuit_num][i].hash(Env::k()) == m_gen_received_label_commitments[circuit_num][i] || m_cc_recv_gen_inp_label_commitments[circuit_num][i].hash(Env::k()) == m_gen_received_label_commitments[circuit_num][2*i+1]);
    //    verify &= (m_cc_recv_gen_inp_label_commitments[circuit_num][i].hash(Env::k()) == m_cc_recv_gen_inp_commitments[circuit_num][i].hash(Env::k()));
  }


  
  if(!verify){
    std::cout << "no verify" << std::endl;
    LOG4CXX_FATAL(logger, "Eval's Check of Gen's Evaluation Circuits Failed");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  std::cout << "done verify" << std::endl;
}

bool BetterYao5::check_received_commitments_vs_generated(std::vector<Bytes> & received, std::vector<commitment_t> & generated){
  assert(received.size() == generated.size());

  std::cout << "checking received vs generated circuits " << std::endl;

  bool verify = true;
  for(int i = 0; i < received.size();i++){
    //verify &= verify_commitment(generated[i],received[i]);
    //verify &= verify_commitment(received[i],generated[i]);
    verify &= (received[i] == decommit(generated[i]).hash(Env::k()));

    /*
    std::cout << "received commitment: "
              << received[i].to_hex()
              << "\tregenerated commitment: "
              << decommit(generated[i]).to_hex()
              << "\thashed decommitment: "
              << decommit(generated[i]).hash(Env::k()).to_hex()
              << "\trank: " << Env::group_rank()
              << std::endl;
    */
  }
  //  verify = false;
  return verify;
}

/**
   STEP 8: RETRIEVE OUTPUTS
 */

void BetterYao5::retrieve_outputs(){
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
  ot_send_init();
  ot_send_random(sender_inputs);
}

void BetterYao5::ot_receive(Bytes selection_bits, std::vector<Bytes> & results_container){
  ot_receive_init();
  ot_receive_random(selection_bits, results_container);
}


void BetterYao5::ot_send_batch(std::vector<std::vector<Bytes> > & sender_inputs){
  ot_send_init();
  for(int j = 0; j < sender_inputs.size(); j++){
    ot_send_random(sender_inputs[j]);
  }
}

void BetterYao5::ot_receive_batch(Bytes selection_bits, std::vector<std::vector<Bytes> > & results_container){
  ot_receive_init();
  for(int j = 0; j < results_container.size();j++){
    ot_receive_random(selection_bits, results_container[j]);
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
      Y[0] = send_inputs[bix*2]; // K[0]
      Y[1] = send_inputs[bix*2+1]; // K[1]

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

  //  m_ot_out.clear();
  // m_ot_out.reserve(2*m_ot_bit_cnt); // the receiver only uses half of it
  // m_timer_gen += MPI_Wtime() - start;
  // m_timer_evl += MPI_Wtime() - start;

  //assert(m_ot_recv_bits.size() >= fit_to_byte_containers(m_ot_bit_cnt));
  //assert(m_chks.size() == Env::node_load());
  
  for (size_t bix = 0; bix < selection_bits.size(); bix++)
    {
      // Step 1: gr=g[b]^r, hr=h[b]^r, where b is the receiver's bit
      start = MPI_Wtime();
      
      //int bit_value = m_ot_recv_bits.get_ith_bit(bix);
      //int bit_value = m_chks[bix];
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
      //      m_ot_out.push_back(Y[bit_value].to_bytes().hash(Env::k()));
      results_container.push_back(Y[bit_value].to_bytes());
      m_timer_evl += MPI_Wtime() - start;
    }
 
  // assert(results_container.size() = selection_bits.size()*8);
  // not quite the right assertion because number of bits need not be
  // an even multiple of 8.
}


//TODO: IMPLEMENT PROPERLY!!!


/**
   THE FUNCTIONS BELOW ARE LEGACY FUNCTIONS
   USED BY BETTERYAO4 (WITH SOME MODIFICATION)
   THE PROTOCOL WORKS, BUT THERE ARE MANY BUGS
   INCLUDING ONE THAT ALLOWS US TO COMMENT OUT
   ALL OF EVL_NEXT_GEN_INP_COMMITMENT
   AND RUN THE PROTOCOL AS IF NOTHING CHANGED
 */



// this function uses an interactive coin flipping protocol between
// Gen and Eval to agree on which circuits will be check circuits
// and which will be evaluation circuits
// it is from the SS11 protocol, but not included in SS13, as Eval
// performs the random selection on her own
/*
void BetterYao5::cut_and_choose()
{
  reset_timers();
  
  double start;
  
  Bytes coins = flip_coins(Env::key_size_in_bytes()); // only roots get the result
  // this collaborative part of cut and choose is not a part of SS13

  if (Env::is_root())
    {
      Prng prng;
      start = MPI_Wtime();
      prng.seed_rand(coins); // use the coins to generate more random bits, deterministically for both sides
      
      // make 60-40 check-vs-evaluation circuit ratio
      m_all_chks.assign(Env::s(), 1);
      
      // FisherÐYates shuffle
      std::vector<uint16_t> indices(m_all_chks.size());
      for (size_t ix = 0; ix < indices.size(); ix++) { indices[ix] = ix; }
      
      // starting from 1 since the 0-th circuit is always evaluation-circuit
      for (size_t ix = 1; ix < indices.size(); ix++)
        {
          int rand_ix = prng.rand_range(indices.size()-ix);
          std::swap(indices[ix], indices[ix+rand_ix]);
        }
      
      int num_of_evls;
      std::cout << "m_all_chks.size: " << m_all_chks.size() << std::endl;
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
      m_timer_gen += MPI_Wtime() - start;
    }
  
  start = MPI_Wtime();
  m_chks.resize(Env::node_load());
  m_timer_evl += MPI_Wtime() - start;
  m_timer_gen += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  MPI_Scatter(&m_all_chks[0], m_chks.size(), MPI_BYTE, &m_chks[0], m_chks.size(), MPI_BYTE, 0, m_mpi_comm);
// distributes m_all_chks to all subprocesses
  // because the prng is seeded with the same (collaboratively tossed) coins, these arrays will be equivalent
  m_timer_mpi += MPI_Wtime() - start;
  
  step_report("cut-&-check");
}
*/

/*
void BetterYao5::cut_and_choose2()
{
  reset_timers();
  
  cut_and_choose2_ot();
  cut_and_choose2_precomputation();
  
  for (size_t ix = 0; ix < m_gcs.size(); ix++)
    {
      cut_and_choose2_evl_circuit(ix);
      cut_and_choose2_chk_circuit(ix);
    
      fprintf(stderr," cut-and-choose2-evl-check: %lu\trank: %i\n",ix,Env::group_rank());
    }
  
  fprintf(stderr,"done cut and choose 2\t rank: %i \n",Env::group_rank());
  step_report("cut-'n-chk2");
}
*/

//
// this function computes OTs between Gen and Eval on random inputs
// the randomness they exchange will be used to seed the prngs which
// Gen and Eval use to generate the circuits and key generators
// Outputs: m_circuit_prngs[]
//
/*
void BetterYao5::cut_and_choose2_ot()
{

  double start;
  m_ot_bit_cnt = Env::node_load();
     
  // initialize OT by carrying out first couple of rounds
  // Eval selects random elements from a cyclic group
  // and sends them to Gen, who then preprocesses
  ot_init();
  
  // Eval selects random inputs for each OT that will happen
  // and sends them to Gen, who masks his own inputs and
  // returns the masks to Eval
  ot_random();
  
  // Gen's m_ot_out has 2*Env::node_load() seeds and
  // Evl's m_ot_out has   Env::node_load() seeds according to m_chks.
  // however, m_chks does not follow SS13, and I want to
  // update the coin flipping/random selection of check circuits
  
  // now that Gen and Eval share some randomness, they will
  // seed the prngs with the random outputs from the OTs
  // these prngs will generate the following in the protocol:
  // at the very least, the prngs are used to mask (and decrypt) 
  // all of the information sent by Gen to Eval about the check/evaluation
  // circuits. Eval will only be able to decrypt information that Gen
  // sends if she has the correct seed for the prng, and this enforces
  // that although Gen sends all of the information check and evaluation
  // information for every circuit, Eval can only decrypt either the check
  // or evaluation information, meaning each circuit to her is one or the other
  // 

  GEN_BEGIN
    start = MPI_Wtime();
  // seeds m_circuit_prngs with seeds from m_ot_out
  seed_m_prngs(2*Env::node_load(), m_ot_out);
  m_timer_gen += MPI_Wtime() - start;
  GEN_END
    
  EVL_BEGIN
    start = MPI_Wtime();
  // seeds m_circuit_prngs with seeds from m_ot_out
  seed_m_prngs(Env::node_load(), m_ot_out);
  
  m_timer_evl += MPI_Wtime() - start;
  
  EVL_END
    
    }
*/

/**
   this function seeds m_circuit_prngs with the seeds provided.
   the number of prngs is provided as a bounds check
   (rather than just using m_circuit_prngs.size())
 */
/*
void BetterYao5::seed_m_prngs(size_t num_circuit_prngs, std::vector<Bytes> seeds){
  assert(seeds.size() == num_circuit_prngs); // this actually not really necessary if assertion is programmatically true. We could just use the size of the seeds, but it is a useful assertion
  m_circuit_prngs.resize(num_circuit_prngs);
  for(size_t ix = 0; ix < num_circuit_prngs; ix++){
    m_circuit_prngs[ix].seed_rand(seeds[ix]);
  }
  std::cout << "prngs seeded " << std::endl;
}
*/


// cut-and-choose2-precomputation
// in this function, Gen (Eval does not execute anything)
// initializes all of the circuits that he will need for the protocol
// and sets a bunch of pointers and constants
// the last part runs through Gen's inputs (and input decommitments)
// and is a source of crashing
// it seems that this (badly) needs to be rewritten
// Outputs: m_circuit_seeds[], m_gen_inp_masks[], m_gcs[].m_gen_inp_decom
//
/*
extern "C" {
void finalize(PCFState *st);
}

void BetterYao5::cut_and_choose2_precomputation()
{
  double start;
  
  GEN_BEGIN

    std::cout << "cut-and-choose2-precomp" << std::endl;
    start = MPI_Wtime();
  for (size_t ix = 0; ix < m_gcs.size(); ix++)
    {
      
      // straight from static circuit implementation
      m_circuit_seeds[ix] = m_prng.rand_bits(Env::k());

      // input masks are the length of the input
      m_gen_inp_masks[ix] = m_prng.rand_bits(m_gen_inp_cnt);
      
      m_gcs[ix].initialize_gen_circuit(m_ot_keys[ix], m_gen_inp_masks[ix], m_circuit_seeds[ix]);

      m_gcs[ix].m_st = 
        load_pcf_file(Env::pcf_file(), m_gcs[ix].m_const_wire, m_gcs[ix].m_const_wire+1, copy_key);
      m_gcs[ix].m_st->alice_in_size = m_gen_inp_cnt;
      m_gcs[ix].m_st->bob_in_size = m_evl_inp_cnt;
      
      set_external_circuit(m_gcs[ix].m_st, &m_gcs[ix]);
      set_key_copy_function(m_gcs[ix].m_st, copy_key);
      set_key_delete_function(m_gcs[ix].m_st, delete_key);
      set_callback(m_gcs[ix].m_st, gen_next_malicious_gate);

      fprintf(stderr, "generate decommitments: index %lu \trank %x\n",ix, Env::group_rank());
      
      //      std::cout << "generate decommitments? " << ix << std::endl;
      std::cout << "gen input count: " << m_gen_inp_cnt << std::endl; // this is the one we care about right now.
      std::cout << "evl input count: " << m_evl_inp_cnt << std::endl;
   
      // this loop either runs through all of Gen's inputs 
      // or terminates when the circuit is out of gates
      // this second check is obviously an error
      // since that should never happen before Gen's inputs are exhausted
      // 
      //while ((m_gcs[ix].m_gen_inp_decom.size()/2 < m_gen_inp_cnt))
      while ((m_gcs[ix].get_gen_decommitments().size()/2 < m_gen_inp_cnt))
        {
          if(!get_next_gate(m_gcs[ix].m_st)){
            fprintf(stderr,"error: circuit completed before finding all of Gen's inputs");
            break;
          }
          m_gcs[ix].get_and_clear_out_bufr();
          // discard the garbled gates for now
          
        }

      std::cout << "done decommitting" << std::endl;
      
      // why do we need this finalize? 
      finalize(m_gcs[ix].m_st);
    }
  m_timer_gen += MPI_Wtime() - start;
  
  std::cout << "end cut-and-choose2-precomp" << std::endl;
  GEN_END
}
*/

// gen does this for every one of the circuits, with index taken as a parameter
// Gen sends Eval his masked input, encrypted with the output of a prng
// he also sends Eval the keys that will be used for constants 1 and 0
// (they should not need to be transmitted every time with the gates)
// finally, he sends the decommitments to his own input keys (why now?)
// note: Gen's even-indexed m_circuit_prngs are for evaluation circuits,
//       and   odd-indexed  m_circuit_prngs are for check circuits
// question: why does Gen need to send Evl his masked input?
//           shouldn't he just send his input keys?
/*
void BetterYao5::cut_and_choose2_evl_circuit(size_t ix)
{
  double start;
  
  Bytes bufr;
  
  fprintf(stderr,"Begin sending things for evaluation circuit\n");

  // send masked gen inp
  GEN_BEGIN
    start = MPI_Wtime();

  assert(m_gen_inp_masks[ix].size() == m_private_input.size());

  // mask gen's input
  bufr = m_gen_inp_masks[ix] ^ m_private_input;

  debug_evl_fprintf("Gen's decrypted inp mask xor input ", bufr.to_hex());

  // and then encrypt it again
  // seriously, why do we need gen's input mask? Evl shouldn't see
  // any form of Gen's input except his input keys
  // which makes this kind of ridiculous:
  // Gen send's Evl masked input which is also masked by a prng
  // (seed for which was exchanged in OT)
  bufr ^= m_circuit_prngs[2*ix+0].rand_bits(bufr.size()*8); // encrypt message
  
  debug_evl_fprintf("Gen send inp mask xor input ", bufr.to_hex());

  m_timer_gen += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  GEN_SEND(bufr);
  m_timer_com += MPI_Wtime() - start;
  
  GEN_END
    
    EVL_BEGIN
    // fprintf(stderr,"cut-choose-2-evl: %lu \t rank: %i\n",ix,Env::group_rank());
    // fprintf(stderr,"cut-and-choose-2-evl: m_chks[ix] %i\trank: %i\n",m_chks[ix],Env::group_rank());
  
  start = MPI_Wtime();
  bufr = EVL_RECV();
  // debug_evl_fprintf("Evl Receive inp mask xor input ", bufr.to_hex());
  m_timer_com += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  // eval can only decrypt if this is an evaluation circuit,
  // and she gets Gen's masked input
  if (!m_chks[ix]) // evaluation circuit
    {
      bufr ^= m_circuit_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
      m_gen_inp_masks[ix] = bufr;
      // fprintf(stderr,"step 1\t rank %i\n",Env::group_rank());
    }
  // debug_evl_fprintf("Evl decrypted inp mask xor input ", bufr.to_hex());
  m_timer_evl += MPI_Wtime() - start;
  EVL_END
    
    m_comm_sz += bufr.size();
  
  // if this is an evaluation circuit, Gen and Eval
  // must agree on some constant keys
  
  // send constant keys m_gcs[ix].m_const_wire
  GEN_BEGIN
    start = MPI_Wtime();
  
  bufr = m_gcs[ix].get_const_key(0, 0) + m_gcs[ix].get_const_key(1, 1);
  // debug_evl_fprintf("Gen const keys", bufr.to_hex());

  bufr ^= m_circuit_prngs[2*ix+0].rand_bits(bufr.size()*8); // encrypt message
  // debug_evl_fprintf("Gen encrypted const keys", bufr.to_hex());
  
  m_timer_gen += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  GEN_SEND(bufr);
  m_timer_com += MPI_Wtime() - start;
  GEN_END
    
    EVL_BEGIN
    start = MPI_Wtime();
  bufr = EVL_RECV();
  //debug_evl_fprintf("Evl encrypted const keys", bufr.to_hex());
  
  m_timer_com += MPI_Wtime() - start;
    
  start = MPI_Wtime();
  if (!m_chks[ix]) // evaluation circuit
    {
      bufr ^= m_circuit_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
      //debug_evl_fprintf("Evl const keys", bufr.to_hex());
      
      std::vector<Bytes> bufr_chunks = bufr.split(Env::key_size_in_bytes());
      m_gcs[ix].set_const_key(0, bufr_chunks[0]);
      m_gcs[ix].set_const_key(1, bufr_chunks[1]);
      //set_const_key(m_gcs[ix], 0, bufr_chunks[0]);
      //set_const_key(m_gcs[ix], 1, bufr_chunks[1]);
      //fprintf(stderr,"step 2\trank %i\n",Env::group_rank());
    }
  m_timer_evl += MPI_Wtime() - start;
  EVL_END
    
  m_comm_sz += bufr.size();
    

  // next, Gen decommits to his input keys

  // send m_gcs[ix].m_gen_inp_decom
  GEN_BEGIN
    //assert(m_gcs[ix].m_gen_inp_decom.size() == 2*m_gen_inp_cnt);
    assert(m_gcs[ix].get_gen_decommitments().size() == 2*m_gen_inp_cnt);
  GEN_END
    
    EVL_BEGIN
    if (!m_chks[ix]) {
      m_gcs[ix].resize_gen_decommitments(m_gen_inp_cnt);
      //m_gcs[ix].get_gen_decommitments().resize(m_gen_inp_cnt); 
      // fprintf(stderr,"step 3\trank: %i\n",Env::group_rank());
    }
  EVL_END
    
    for (size_t jx = 0; jx < m_gen_inp_cnt; jx++)
      {
        GEN_BEGIN
          start = MPI_Wtime();
        // this bit selects which decommitment to select and send
        
        byte bit = m_private_input.get_ith_bit(jx) ^ m_gen_inp_masks[ix].get_ith_bit(jx);
        //bufr = m_gcs[ix].m_gen_inp_decom[2*jx+bit];
        
        bufr = m_gcs[ix].get_gen_decommitments()[2*jx+bit];
        //debug_evl_fprintf("Gen ith input decommitment", bufr.to_hex());
        bufr ^= m_circuit_prngs[2*ix+0].rand_bits(bufr.size()*8); // encrypt message
        // debug_evl_fprintf("Gen ith input decommitment encrypted", bufr.to_hex());
        
        // remember, even m_circuit_prngs are for evaluation circuits
        // and odd m_circuit_prngs are for check circuits
        m_timer_gen += MPI_Wtime() - start;
        
        start = MPI_Wtime();
        GEN_SEND(bufr);
        m_timer_com += MPI_Wtime() - start;
        GEN_END
          
          EVL_BEGIN
          start = MPI_Wtime();
        bufr = EVL_RECV();
        m_timer_com += MPI_Wtime() - start;
 
        //debug_evl_fprintf("Evl ith input decommitment encrypted", bufr.to_hex());
       
        start = MPI_Wtime();
        if (!m_chks[ix]) // evaluation circuit
          {
            //fprintf(stderr,"step 4 start\t rank: %i\n",Env::group_rank());
            //!! BUG HERE
            // eval receives decommitment and decrypts it
            // she now has the decommitments to Gen's input keys
            bufr ^= m_circuit_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
            //debug_evl_fprintf("Evl ith input decommitment", bufr.to_hex());
           //fprintf(stderr,"step 4 mid\t rank: %i\n",Env::group_rank());
            //m_gcs[ix].m_gen_inp_decom[jx] = bufr;
            //m_gcs[ix].get_gen_decommitments()[jx] = bufr;
            m_gcs[ix].set_gen_decommitment(jx,bufr);
            //fprintf(stderr,"step 4\trank: %i\n",Env::group_rank());
          }
        m_timer_evl += MPI_Wtime() - start;
        EVL_END
          
          m_comm_sz += bufr.size();
      }
  fprintf(stderr,"end cut-and-choose2 evl %lu\trank: %i\n",ix,Env::group_rank());
     
}
*/

// Gen must send Evl info for check circuits
// he sends the masks for his own input,
// random seeds used to generate:
//     <something, probably the rest of the circuit)>,
// all of his OT keys for Eval's inputs
// and all of his own input decommitments
// Gen always sends these, but Eval can only decrypt if 
// her prng was properly seeded by the OT result
// note: Gen's even-indexed m_circuit_prngs are for evaluation circuits,
//       and   odd-indexed  m_circuit_prngs are for check circuits
/*
void BetterYao5::cut_and_choose2_chk_circuit(size_t ix)
{
  double start;
  
  Bytes bufr;
  std::vector<Bytes> bufr_chunks;
  
  // send m_gen_inp_masks[ix]
  GEN_BEGIN
    start = MPI_Wtime();
 
  // input mask was generated randomly in precomp
  // and is the length of Gen's inputs (m_gen_inp_cnt)
  bufr = m_gen_inp_masks[ix];
  assert(bufr.size() == m_gen_inp_cnt/8);
  
  // encrypt the input mask with some randomness
  // use the 2*ix+1 prng for check circuits (checks get odd values)
  bufr ^= m_circuit_prngs[2*ix+1].rand_bits(bufr.size()*8);
  m_timer_gen += MPI_Wtime() - start;
  
  // and send encrypted input mask to Evl
  start = MPI_Wtime();
  GEN_SEND(bufr);
  m_timer_com += MPI_Wtime() - start;
  GEN_END
    

    // Eval gets the input mask,
    // and if this is a check circuit,
    // then she decrypts with a properly seeded generator
    EVL_BEGIN
    fprintf(stderr,"cut-and-chooose2-check: %lu\trank: %i\n",ix,Env::group_rank());
    start = MPI_Wtime();
  bufr = EVL_RECV();
  m_timer_com += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  if (m_chks[ix]) // check circuit
    {
      bufr ^= m_circuit_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
      m_gen_inp_masks[ix] = bufr;
    }
  m_timer_evl += MPI_Wtime() - start;
  EVL_END
    
    m_comm_sz += bufr.size();
  

  // next, Gen sends random seeds to Eval
  
  // send m_circuit_seeds[ix]
  GEN_BEGIN
    start = MPI_Wtime();
  bufr = m_circuit_seeds[ix];
  bufr ^= m_circuit_prngs[2*ix+1].rand_bits(bufr.size()*8); // encrypt message
  m_timer_gen += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  GEN_SEND(bufr);
  m_timer_com += MPI_Wtime() - start;
  GEN_END
    
    EVL_BEGIN
    start = MPI_Wtime();
  bufr = EVL_RECV();
  m_timer_com += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  if (m_chks[ix]) // check circuit
    {
      bufr ^= m_circuit_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
      m_circuit_seeds[ix] = bufr;
    }
  m_timer_evl += MPI_Wtime() - start;
  EVL_END

    m_comm_sz += bufr.size();


  // send m_ot_keys[ix]
  // these are the ot keys used for eval's inputs
  GEN_BEGIN
    assert(m_ot_keys[ix].size() == 2*m_evl_inp_cnt);
  GEN_END

    EVL_BEGIN
    if (m_chks[ix]) { m_ot_keys[ix].resize(2*m_evl_inp_cnt); }
  EVL_END

    for (size_t jx = 0; jx < 2*m_evl_inp_cnt; jx++)
      {
        GEN_BEGIN
          start = MPI_Wtime();
        bufr = m_ot_keys[ix][jx];
        bufr ^= m_circuit_prngs[2*ix+1].rand_bits(bufr.size()*8); // encrypt message
        m_timer_gen += MPI_Wtime() - start;

        start = MPI_Wtime();
        GEN_SEND(bufr);
        m_timer_com += MPI_Wtime() - start;
        GEN_END

          EVL_BEGIN
          start = MPI_Wtime();
        bufr = EVL_RECV();
        m_timer_com += MPI_Wtime() - start;

        start = MPI_Wtime();
        if (m_chks[ix]) // check circuit
          {
            bufr ^= m_circuit_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
            m_ot_keys[ix][jx] = bufr;
          }
        m_timer_evl += MPI_Wtime() - start;
        EVL_END

          m_comm_sz += bufr.size();
      }

  // send m_gcs[ix].m_gen_inp_decom
  // this is all of the keys to gen's inputs
  GEN_BEGIN
    // assert(m_gcs[ix].m_gen_inp_decom.size() == 2*m_gen_inp_cnt);
    assert(m_gcs[ix].get_gen_decommitments().size() == 2*m_gen_inp_cnt);
  GEN_END

    EVL_BEGIN
    if (m_chks[ix]) { m_gen_inp_decom[ix].resize(2*m_gen_inp_cnt); }
  EVL_END

    for (size_t jx = 0; jx < 2*m_gen_inp_cnt; jx++)
      {
        GEN_BEGIN
          start = MPI_Wtime();
        bufr = m_gcs[ix].get_gen_decommitments()[jx]; //m_gen_inp_decom[jx];
        // why use 2*ix+1?
        bufr ^= m_circuit_prngs[2*ix+1].rand_bits(bufr.size()*8); // encrypt message
        m_timer_gen += MPI_Wtime() - start;

        start = MPI_Wtime();
        GEN_SEND(bufr);
        m_timer_com += MPI_Wtime() - start;
        GEN_END

          EVL_BEGIN
          start = MPI_Wtime();
        bufr = EVL_RECV();
        m_timer_com += MPI_Wtime() - start;

        start = MPI_Wtime();
        if (m_chks[ix]) // check circuit
          {
            bufr ^= m_circuit_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
            m_gen_inp_decom[ix][jx] = bufr;
          }
        m_timer_evl += MPI_Wtime() - start;
        EVL_END
          
          m_comm_sz += bufr.size();
      }
  fprintf(stderr,"end cut-and-choose2-check %lu\trank: %i\n",ix,Env::group_rank());
}
*/

// in this consistency check, 
// first, Gen and Evl agree on a 2-UHF
// note: in order to ensure Gen's input consistency
//       agreement on the 2-UHF function must happen
//       after Gen commits to his inputs
// then, gen creates his input commitments
// and sends them to Evl, who evaluates them
// finally, Evl saves the input hashes from all of the circuits
/*
void BetterYao5::consistency_check()
{
  // std::cout << "const check start" << std::endl;

  reset_timers();

  Bytes bufr;
  std::vector<Bytes> bufr_chunks;

  double start;

  // jointly pick a 2-UHF matrix
  // m_gen_inp_cnt is almost the right number to use here, since
  // the number of random bits should equal k * gen's input size,
  // but in the current protocol Gen does not supplement his inputs
  // by selecting random ciphertext "c" (used for masking his output)
  // and he doesn't add the extra 2k+lg(k) bits
  bufr = flip_coins(Env::k()*((m_gen_inp_cnt+7)/8)); // only roots get the result

  // must resize buffer for all the subprocesses before sending it
  start = MPI_Wtime();
  bufr.resize(Env::k()*((m_gen_inp_cnt+7)/8));
  m_timer_evl += MPI_Wtime() - start;
  m_timer_gen += MPI_Wtime() - start;

  start = MPI_Wtime();
  MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm);
  m_timer_mpi += MPI_Wtime() - start;

  start = MPI_Wtime();
        
  // create m rows of length k
  m_2UHF_matrix = bufr.split(bufr.size()/Env::k());
  
  m_timer_evl += MPI_Wtime() - start;
  m_timer_gen += MPI_Wtime() - start;

  // std::cout << "agree on UHF" << std::endl;

  // now everyone agrees on the UHF given by m_2UHF_matrix
  for (size_t ix = 0; ix < m_gcs.size(); ix++)
    for (size_t kx = 0; kx < m_2UHF_matrix.size(); kx++) // kx ranges from 0 to m1
      {
        // in this loop, 
        GEN_BEGIN
        start = MPI_Wtime();
        m_gcs[ix].gen_next_gen_inp_com(m_2UHF_matrix[kx], kx);
        // this doesn't really generate an input commitment. we shouldn't call it that.
        bufr = m_gcs[ix].get_and_clear_out_bufr();
        m_timer_gen += MPI_Wtime() - start;
              
        start = MPI_Wtime();
        GEN_SEND(bufr);
        m_timer_com += MPI_Wtime() - start;
              
        GEN_END
                
        EVL_BEGIN
        start = MPI_Wtime();
        bufr = EVL_RECV();
        m_timer_com += MPI_Wtime() - start;

        if (!m_chks[ix]) // evaluation circuit
          {
            start = MPI_Wtime();
            m_gcs[ix].clear_and_replace_in_bufr(bufr);
            m_gcs[ix].evl_next_gen_inp_com(m_2UHF_matrix[kx], kx);
            m_timer_evl += MPI_Wtime() - start;
          }
        EVL_END
                
          m_comm_sz += bufr.size();
      }

  // std::cout << "EVL check hashes" << std::endl;

  EVL_BEGIN
    for (size_t ix = 0; ix < m_gcs.size(); ix++)
      if (!m_chks[ix]) // if evaluation circuit, save this hash
        {
          std::cout << "save Gen input hash: " << m_gcs[ix].get_gen_inp_hash().to_hex() << std::endl;
          m_gen_inp_hash[ix] = m_gcs[ix].get_gen_inp_hash();//m_gen_inp_hash;
        }
  EVL_END
          
    step_report("const-check");
}
*/

/*
void BetterYao5::circuit_evaluate()
{
	reset_timers();

	double start;

	int verify = 1;
	Bytes bufr;

        std::cout << "begin circuit evaluate" << std::endl;

	for (size_t ix = 0; ix < m_gcs.size(); ix++)
	{
          GEN_BEGIN
            start = MPI_Wtime();

          m_gcs[ix].gen_init_circuit(m_ot_keys[ix], m_gen_inp_masks[ix], m_circuit_seeds[ix]);   
          std::cout << ix << " gen const 0: " << m_gcs[ix].get_const_key( 0, 0).to_hex() << std::endl;
          std::cout << ix << " gen const 1: " << m_gcs[ix].get_const_key( 1, 1).to_hex() << std::endl;
          m_timer_gen += MPI_Wtime() - start;
          GEN_END
              
              EVL_BEGIN
              std::cout << "eval start evaluating" << std::endl;
            
            start = MPI_Wtime();
            if (m_chks[ix]) // check-circuits
              {
                m_gcs[ix].gen_init_circuit(m_ot_keys[ix], m_gen_inp_masks[ix], m_circuit_seeds[ix]);
              }
            else // evaluation-circuits
              {
                m_gcs[ix].evl_init_circuit(m_ot_keys[ix], m_gen_inp_masks[ix], m_private_input);
              }
            m_timer_evl += MPI_Wtime() - start;

            //std::cout << ix << " evl const 0: " << get_const_key(m_gcs[ix], 0, 0).to_hex() << std::endl;
            //std::cout << ix << " evl const 1: " << get_const_key(m_gcs[ix], 1, 1).to_hex() << std::endl;
            EVL_END
              
              std::cout << "load pcf file" << std:: endl;
            start = MPI_Wtime();
            m_gcs[ix].m_st = 
              load_pcf_file(Env::pcf_file(), m_gcs[ix].m_const_wire, m_gcs[ix].m_const_wire+1, copy_key);
            m_gcs[ix].m_st->alice_in_size = m_gen_inp_cnt;
            m_gcs[ix].m_st->bob_in_size = m_evl_inp_cnt;

            // the next three should have been done by Gen during precomputation
            // important for Eval (and all the other mpi cores)
            set_external_circuit(m_gcs[ix].m_st, &m_gcs[ix]);
            set_key_copy_function(m_gcs[ix].m_st, copy_key);
            set_key_delete_function(m_gcs[ix].m_st, delete_key);

            m_timer_gen += MPI_Wtime() - start;
            m_timer_evl += MPI_Wtime() - start;
            
            GEN_BEGIN // generate and send the circuit gate-by-gate
              
              std::cout << "gen generate and send" << std::endl;
            start = MPI_Wtime();
            //set_callback(m_gcs[ix].m_st, gen_next_gate_m);
            set_callback(m_gcs[ix].m_st, gen_next_malicious_gate);
            
            while (get_next_gate(m_gcs[ix].m_st))
              {
                bufr = m_gcs[ix].get_and_clear_out_bufr();
                m_timer_gen += MPI_Wtime() - start;
                
                start = MPI_Wtime();
                std::cout << "Gen Send Gate: " << bufr.to_hex() << std::endl; 
                GEN_SEND(bufr);
                m_timer_com += MPI_Wtime() - start;
                
                m_comm_sz += bufr.size();
                
                start = MPI_Wtime(); // start m_timer_gen
              }
            m_timer_gen += MPI_Wtime() - start;
            
            GEN_SEND(Bytes(0)); // a redundant value to prevent the evaluator from hanging
            GEN_END
              
              EVL_BEGIN // receive and evaluate the circuit gate-by-gate
              
              std::cout << "eval receive and evaluate" << std::endl;
            if (m_chks[ix]) // check circuit
              {
                fprintf(stderr,"eval check circuit: index %lu, rank %i\n",ix, Env::group_rank());
                start = MPI_Wtime();
                set_callback(m_gcs[ix].m_st, gen_next_malicious_gate);
                //set_callback(m_gcs[ix].m_st, gen_next_gate_m);
                while (get_next_gate(m_gcs[ix].m_st))
                  {
                    
                    bufr = m_gcs[ix].get_and_clear_out_bufr();
                    m_timer_evl += MPI_Wtime() - start;
                    
                     start = MPI_Wtime();
                    Bytes recv = EVL_RECV();
                    m_timer_com += MPI_Wtime() - start;
                    
                    m_comm_sz += bufr.size();
                    
                    start = MPI_Wtime(); // start m_timer_evl
                    verify &= (bufr == recv);
                  }
                m_timer_gen += MPI_Wtime() - start;
                
                EVL_RECV(); // a redundant value to prevent the evlauator from hanging
              }
            else // evaluation circuit
              {
                fprintf(stderr,"eval evaluation circuit: index %lu, rank %i\n",ix, Env::group_rank());
                start = MPI_Wtime();
                set_callback(m_gcs[ix].m_st, evl_next_malicious_gate);
                std::cout<< "enter do loop" << std::endl;
                do {
                  // std::cout << "timing stuff" << std::endl;
                  m_timer_evl += MPI_Wtime() - start;
                  
                  start = MPI_Wtime();
                  bufr = EVL_RECV();
                  std::cout << "Eval receive gate: " << bufr.to_hex() << std::endl;
                  m_timer_com += MPI_Wtime() - start;
                  
                  // std::cout << "buffer size add" << std::endl;
                  
                  m_comm_sz += bufr.size();
                  
                  start = MPI_Wtime();
                  m_gcs[ix].clear_and_replace_in_bufr(bufr);
                  
                  //std::cout << "received" << std::endl;
                  //fprintf(stderr, "received");
                  
                  // std::cout << "got a gate" << std::endl;                  
                } while (get_next_gate(m_gcs[ix].m_st));
                
                std::cout << "complete" << std::endl;
                m_timer_evl += MPI_Wtime() - start;
              }
            
            EVL_END
              }

	EVL_BEGIN // check the hash of all the garbled circuits
          int all_verify = 0;
        
        std::cout << "evl check hash of garbled circuits" << std::endl;
        start = MPI_Wtime();
        MPI_Reduce(&verify, &all_verify, 1, MPI_INT, MPI_LAND, 0, m_mpi_comm);
        m_timer_mpi += MPI_Wtime() - start;
        
        start = MPI_Wtime();
        if (Env::is_root() && !all_verify)
          {
            LOG4CXX_FATAL(logger, "Verification failed");
            //MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
          }
        
        Bytes gen_inp_hash;
        
        for (size_t ix = 0; ix < m_gcs.size(); ix++)
          {
            // check the commitments associated with the generator's input wires
            if (m_chks[ix]) // check circuit
              {
                for (size_t jx = 0; jx < m_gen_inp_cnt*2; jx++)
                  {
                    if (!(m_gen_inp_decom[ix][jx] == m_gcs[ix].get_gen_decommitments()[jx]))//.m_gen_inp_decom[jx]))
                      {
                        LOG4CXX_FATAL(logger, "Commitment Verification Failure (check circuit)");
                        //MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                      }
                  }
              }
            else // evaluation circuit
              {
                
                if(!m_gcs[ix].pass_check())
                  {
                    LOG4CXX_FATAL(logger, "Commitment Verification Failure (evaluation circuit)");
                    //MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                  }
                
                if (gen_inp_hash.size() == 0)
                  {
                    std::cout << "gen inp hash is empty?" << std::endl;
                    gen_inp_hash = m_gen_inp_hash[ix];
                  }
                else if (!(gen_inp_hash == m_gen_inp_hash[ix]))
                  {
                    LOG4CXX_FATAL(logger, "Generator Input Hash Inconsistent Failure (evaluation circuit)");
                    //MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                  }
              }
            
            m_gcs[ix].trim_output();
            //trim_output(m_gcs[ix]);
          }
        m_timer_evl += MPI_Wtime() - start;
	EVL_END
          
          step_report("circuit-evl");
        
	if (m_gcs[0].get_evl_out_ix() != 0) //m_evl_out_ix != 0)
          proc_evl_out();
        
        if (m_gcs[0].get_gen_out_ix() != 0) //m_gen_out_ix != 0)
          proc_gen_out();
}
*/

/*
void BetterYao5::proc_evl_out()
{
  fprintf(stderr,"begin evl out\n");
  EVL_BEGIN
    reset_timers();
  
  double start;
  Bytes send, recv;
  
  start = MPI_Wtime();
  for (size_t ix = 0; ix < m_gcs.size(); ix++) // fill zeros for uniformity (convenient for MPIs)
    {
      send += m_gcs[ix].get_evl_out(); // m_evl_out;
    }
  
  if (Env::is_root() == 0)
    {
      fprintf(stderr,"node amount %i\n", Env::node_amnt());
      recv.resize(send.size()*Env::node_amnt());
    }
  m_timer_evl += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  MPI_Gather(&send[0], send.size(), MPI_BYTE, &recv[0], send.size(), MPI_BYTE, 0, m_mpi_comm);
  m_timer_mpi += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  if (Env::is_root())
    {
      size_t chks_total = 0;
      for (size_t ix = 0; ix < m_all_chks.size(); ix++)
        chks_total += m_all_chks[ix];
      
      // find majority by locating the median of output from evaluation-circuits
      
      // this line is OBVIOUSLY wrong. evl_out_cnt was never anything but 0
      // std::vector<Bytes> vec = recv.split((Env::circuit().evl_out_cnt()+7)/8); 
      std::vector<Bytes> vec = recv.split(0);
      
      size_t median_ix = (chks_total+vec.size())/2;
      std::nth_element(vec.begin(), vec.begin()+median_ix, vec.end());
      
      m_evl_out = *(vec.begin()+median_ix);
    }
  m_timer_evl += MPI_Wtime() - start;
  
  step_report_no_sync("chk-evl-out");
  EVL_END
}


void BetterYao5::proc_gen_out()
{
  // THIS SIMPLY TAKES THE OUTPUT OF THE FIRST CIRCUIT!
  // NO MAJORITY OPERATION IS DONE
  // if Gen cheats, that's his problem?
  // not if Evl tries to mess with the first circuit

	reset_timers();

	// TODO: implement Ki08
	//m_gen_out = m_gcs[0].get_gen_out();//.m_gen_out;

	EVL_BEGIN
          m_gen_out = m_gcs[0].get_gen_out();//m_gen_out;
        EVL_SEND(m_gen_out);
	EVL_END

	GEN_BEGIN
        m_gen_out = GEN_RECV();
        GEN_END

	step_report("chk-gen-out");
}
*/
/** 
    witness indestinguishable proof of Gen's output authenticity
    Gen has output keys {u_0^j,u_1^j}for j in all circuits
    Eval knows the index m of the majority circuit (representing other circuits)
    and the random key v corresponding to Gen's output wire of semantic value a
    they share the security parameter, statistical parameter, commitments to Gen's output keys, and Gen's output a, which presumably Evl tells Gen in the clear (protected by their secure channel)
*/



//
// Implementation of "Two-Output Secure Computation with Malicious Adversaries"
// by abhi shelat and Chih-hao Shen from EUROCRYPT'11 (Protocol 2)
//
// The evaluator (sender) generates m_ot_bit_cnt pairs of k-bit random strings, and
// the generator (receiver) has input m_ot_bits and will receive output m_ot_out.
//
/*
void BetterYao5::ot_init()
{
  double start;
  
  start = MPI_Wtime();
  std::vector<Bytes> bufr_chunks;
  Bytes bufr(Env::elm_size_in_bytes()*4);
  
  Z y, a;
  m_timer_gen += MPI_Wtime() - start;
  m_timer_evl += MPI_Wtime() - start;
  
  // step 1: ZKPoK of the CRS: g[0], h[0], g[1], h[1]
  if (Env::is_root())
    {
      EVL_BEGIN // evaluator (OT receiver)
        start = MPI_Wtime();
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
      m_timer_evl += MPI_Wtime() - start;
      
      start = MPI_Wtime(); // send to Gen
      EVL_SEND(bufr);
      m_timer_com += MPI_Wtime() - start;
      EVL_END
        
        GEN_BEGIN // generator (OT sender)
        start = MPI_Wtime();
      bufr = GEN_RECV();
      m_timer_com += MPI_Wtime() - start;
      GEN_END
        
        m_comm_sz += bufr.size();
    }
  
  // send g[0], g[1], h[0], h[1] to slave processes
  start = MPI_Wtime();
  MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm);
  m_timer_mpi += MPI_Wtime() - start;
  
  start = MPI_Wtime();
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
  m_timer_gen += MPI_Wtime() - start;
  m_timer_evl += MPI_Wtime() - start;
}


void BetterYao5::ot_random()
{
	double start;

	start = MPI_Wtime();
        Bytes send, recv;
        std::vector<Bytes> recv_chunks;
        
        Z r, s[2], t[2];
        G gr, hr, X[2], Y[2];
        
        m_ot_out.clear();
        m_ot_out.reserve(2*m_ot_bit_cnt); // the receiver only uses half of it
        m_timer_gen += MPI_Wtime() - start;
        m_timer_evl += MPI_Wtime() - start;

	EVL_BEGIN // evaluator (OT receiver)
          assert(m_chks.size() == Env::node_load());
        assert(m_ot_bit_cnt == Env::node_load());
        
	//for (size_t bix = 0; bix < m_ot_bit_cnt; bix++)
        for(size_t bix = 0; bix < m_chks.size(); bix++)
        {
            // Step 1: gr=g[b]^r, hr=h[b]^r, where b is the receiver's bit
            start = MPI_Wtime();
            int bit_value = m_chks[bix];
            
            r.random();
            
            gr = m_ot_g[bit_value]^r;
            hr = m_ot_h[bit_value]^r;
            
            send.clear();
            send += gr.to_bytes();
            send += hr.to_bytes();
            m_timer_evl += MPI_Wtime() - start;
            
            start = MPI_Wtime();
            EVL_SEND(send);
            
            // Step 2: the evaluator computes X[0], Y[0], X[1], Y[1]
            recv.clear();
            recv += EVL_RECV(); // receive X[0], Y[0], X[1], Y[1]
            m_timer_com += MPI_Wtime() - start;
            
            m_comm_sz += send.size() + recv.size();
            
            // Step 3: the evaluator computes K = Y[b]/X[b]^r
            start = MPI_Wtime();
            recv_chunks = recv.split(Env::elm_size_in_bytes());
            
            X[bit_value].from_bytes(recv_chunks[    bit_value]); // X[b]
            Y[bit_value].from_bytes(recv_chunks[2 + bit_value]); // Y[b]
            
            // K = Y[b]/(X[b]^r)
            Y[bit_value] /= X[bit_value]^r;
            m_ot_out.push_back(Y[bit_value].to_bytes().hash(Env::k()));
            m_timer_evl += MPI_Wtime() - start;
          }
        
        assert(m_ot_out.size() == m_ot_bit_cnt);
	EVL_END

	GEN_BEGIN // generator (OT sender)
		for (size_t bix = 0; bix < m_ot_bit_cnt; bix++)
		{
			// Step 1: gr=g[b]^r, hr=h[b]^r, where b is the receiver's bit
			start = MPI_Wtime();
				recv.clear();
				recv += GEN_RECV(); // receive gr, hr
			m_timer_com += MPI_Wtime() - start;

			m_comm_sz += recv.size();

			// Step 2: the evaluator computes X[0], Y[0], X[1], Y[1]
			start = MPI_Wtime();
				recv_chunks = recv.split(Env::elm_size_in_bytes());

				gr.from_bytes(recv_chunks[0]);
				hr.from_bytes(recv_chunks[1]);

				Y[0].random(); Y[1].random(); // K[0], K[1] sampled at random
                                
				m_ot_out.push_back(Y[0].to_bytes().hash(Env::k()));
				m_ot_out.push_back(Y[1].to_bytes().hash(Env::k()));

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
			m_timer_gen += MPI_Wtime() - start;

			start = MPI_Wtime();
				GEN_SEND(send);
			m_timer_com += MPI_Wtime() - start;

			m_comm_sz += send.size();
		}

		assert(m_ot_out.size() == 2*m_ot_bit_cnt);
	GEN_END
}


*/
