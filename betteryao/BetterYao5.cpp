#include "BetterYao5.h"
// #include "garbled_circuit.h"
#include <algorithm>

#include <unistd.h>

#include <log4cxx/logger.h>
static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("BetterYao5.cpp"));


BetterYao5::BetterYao5(EnvParams &params) : YaoBase(params), m_ot_bit_cnt(0)
{

  std::cout << "node load: " << Env::node_load() << std::endl;
  // Init variables
  m_rnds.resize(Env::node_load());
  m_gcs.resize(Env::node_load());
  
  for (size_t ix = 0; ix < m_gcs.size(); ix++)
    {
      //initialize_malicious_circuit(m_gcs[ix]);
      m_gcs[ix] = new GarbledMal();
      m_gcs[ix].initialize_circuit();
  }
  
  m_gen_inp_hash.resize(Env::node_load());
  m_gen_inp_masks.resize(Env::node_load());
  m_gen_inp_decom.resize(Env::node_load());
  
  get_and_size_inputs();
}


void BetterYao5::start()
{
  oblivious_transfer();
  cut_and_choose();
  cut_and_choose2();
  consistency_check();
  circuit_evaluate();
  final_report();
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
      
      start = MPI_Wtime();
      Bytes coins = m_prng.rand_bits(len_in_bytes*8);	// Step 0: flip coins
      m_timer_gen += MPI_Wtime() - start;
      m_timer_evl += MPI_Wtime() - start;
      
      
      GEN_BEGIN
      start = MPI_Wtime();
      commit_value = m_prng.rand_bits(Env::k()) + coins;	// Step 1: commit to coins
      commitment = commit_value.hash(Env::k());
      m_timer_gen += MPI_Wtime() - start;
      
      start = MPI_Wtime();
      GEN_SEND(commitment);
      remote_coins = GEN_RECV();     	// Step 2: receive alice's coins
      // Gen can decommit to coins only after receiving Alice's
      GEN_SEND(commit_value);		// Step 3: decommit to the coins
      m_timer_com += MPI_Wtime() - start;
      GEN_END
        
      EVL_BEGIN
      start = MPI_Wtime();
      commitment = EVL_RECV();	    	// Step 1: receive bob's commitment
      EVL_SEND(coins);	     	// Step 2: send coins to bob
      commit_value = EVL_RECV();
      m_timer_com += MPI_Wtime() - start;
      
      start = MPI_Wtime();
      if (!(commit_value.hash(Env::k()) == commitment))		// Step 3: check bob's decommitment
        {
          LOG4CXX_FATAL(logger, "commitment to coins can't be properly opened");
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
      remote_coins = Bytes(commit_value.begin()+Env::k()/8, commit_value.end());
      m_timer_evl += MPI_Wtime() - start;
      EVL_END
        
        m_comm_sz = commitment.size() + remote_coins.size() + commit_value.size();
      
      start = MPI_Wtime();
      
      coins ^= remote_coins;
      // combine randomnesses from both players
      
      bufr.swap(coins);
      //bufr = coins;
      
      m_timer_evl += MPI_Wtime() - start;
      m_timer_gen += MPI_Wtime() - start;
    }
  
  return bufr;
}


// this function uses an interactive coin flipping protocol between
// Gen and Eval to agree on which circuits will be check circuits
// and which will be evaluation circuits
// it is from the SS11 protocol, but not included in SS13, as Eval
// performs the random selection on her own
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


void BetterYao5::cut_and_choose2()
{
  reset_timers();
  
  cut_and_choose2_ot();
  cut_and_choose2_precomputation();
  
  for (size_t ix = 0; ix < m_gcs.size(); ix++)
    {
      cut_and_choose2_evl_circuit(ix);
      cut_and_choose2_chk_circuit(ix);
    }
  
  step_report("cut-'n-chk2");
}

//
// this function computes OTs between Gen and Eval on random inputs
// the randomness they exchange will be used to seed the prngs which
// Gen and Eval use to generate the circuits and key generators
// Outputs: m_prngs[]
//
void BetterYao5::cut_and_choose2_ot()
{

  double start;
  m_ot_bit_cnt = Env::node_load();
  
  EVL_BEGIN
  start = MPI_Wtime();
  m_ot_recv_bits.resize((m_ot_bit_cnt+7)/8);
  for (size_t ix = 0; ix < m_chks.size(); ix++)
    {
      m_ot_recv_bits.set_ith_bit(ix, m_chks[ix]);
    }
  m_timer_evl += MPI_Wtime() - start;
  EVL_END
   
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
  // seeds m_prngs with seeds from m_ot_out
  seed_m_prngs(2*Env::node_load(), m_ot_out);
  m_timer_gen += MPI_Wtime() - start;
  GEN_END
    
  EVL_BEGIN
    start = MPI_Wtime();
  // seeds m_prngs with seeds from m_ot_out
  seed_m_prngs(Env::node_load(), m_ot_out);
  
  m_timer_evl += MPI_Wtime() - start;
  
  EVL_END
    
    }

/**
   this function seeds m_prngs with the seeds provided.
   the number of prngs is provided as a bounds check
   (rather than just using m_prngs.size())
 */
void BetterYao5::seed_m_prngs(size_t num_prngs, std::vector<Bytes> seeds){
  assert(seeds.size() == num_prngs); // this actually not really necessary if assertion is programmatically true. We could just use the size of the seeds, but it is a useful assertion
  m_prngs.resize(num_prngs);
  for(size_t ix = 0; ix < num_prngs; ix++){
    m_prngs[ix].seed_rand(seeds[ix]);
  }
  std::cout << "prngs seeded " << std::endl;
}



// cut-and-choose2-precomputation
// in this function, Gen (Eval does not execute anything)
// initializes all of the circuits that he will need for the protocol
// and sets a bunch of pointers and constants
// the last part runs through Gen's inputs (and input decommitments)
// and is a source of crashing
// it seems that this (badly) needs to be rewritten
// Outputs: m_rnds[], m_gen_inp_masks[], m_gcs[].m_gen_inp_decom
//

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
      m_rnds[ix] = m_prng.rand_bits(Env::k());

      // input masks are the length of the input
      m_gen_inp_masks[ix] = m_prng.rand_bits(m_gen_inp_cnt);
      
      // gen_init_circuit(m_gcs[ix], m_ot_keys[ix], m_gen_inp_masks[ix], m_rnds[ix]);
      m_gcs[ix].initialize_gen_circuit(m_ot_keys[ix], m_gen_inp_masks[ix], m_rnds[ix]);


      m_gcs[ix].m_st = 
        load_pcf_file(Env::pcf_file(), m_gcs[ix].m_const_wire, m_gcs[ix].m_const_wire+1, copy_key);
      m_gcs[ix].m_st->alice_in_size = m_gen_inp_cnt;
      m_gcs[ix].m_st->bob_in_size = m_evl_inp_cnt;
      
      set_external_circuit(m_gcs[ix].m_st, &m_gcs[ix]);
      set_key_copy_function(m_gcs[ix].m_st, copy_key);
      set_key_delete_function(m_gcs[ix].m_st, delete_key);
      set_callback(m_gcs[ix].m_st, gen_next_gate);

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
          //get_and_clear_out_bufr(m_gcs[ix]); // discard the garbled gates for now
          
          //if(m_gcs[ix].m_gen_inp_decom.size() % 16 == 0){
            //  last_decom_size = m_gcs[ix].m_gen_inp_decom.size();
          //  std::cout << "decommitments: " << m_gcs[ix].m_gen_inp_decom.size() << std::endl;
          // }
        }

      std::cout << "done decommitting" << std::endl;
      
      // this runs through Gen's inputs. but why?
      // std::cout << "finalize:..." << std::endl;

      // why do we need this finalize? 
      finalize(m_gcs[ix].m_st);
    }
  m_timer_gen += MPI_Wtime() - start;
  
  std::cout << "end cut-and-choose2-precomp" << std::endl;
  GEN_END
}


// gen does this for every one of the circuits, with index taken as a parameter
// Gen sends Eval his masked input, encrypted with the output of a prng
// he also sends Eval the keys that will be used for constants 1 and 0
// (they should not need to be transmitted every time with the gates)
// finally, he sends the decommitments to his own input keys (why now?)
// note: Gen's even-indexed m_prngs are for evaluation circuits,
//       and   odd-indexed  m_prngs are for check circuits
// question: why does Gen need to send Evl his masked input?
//           shouldn't he just send his input keys?
void BetterYao5::cut_and_choose2_evl_circuit(size_t ix)
{
  double start;
  
  Bytes bufr;
  
  // send masked gen inp
  GEN_BEGIN
    start = MPI_Wtime();

  //assert(m_gen_inp_masks[ix].size() == m_gen_inp.size());
  assert(m_gen_inp_masks[ix].size() == m_private_input.size());
  // mask gen's input
  //  bufr = m_gen_inp_masks[ix] ^ m_gen_inp;
  bufr = m_gen_inp_masks[ix] ^ m_private_input;
  // and then encrypt it again
  bufr ^= m_prngs[2*ix+0].rand_bits(bufr.size()*8); // encrypt message
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
  // eval can only decrypt if this is an evaluation circuit,
  // and she gets Gen's masked input
  if (!m_chks[ix]) // evaluation circuit
    {
      bufr ^= m_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
      m_gen_inp_masks[ix] = bufr;
    }
  m_timer_evl += MPI_Wtime() - start;
  EVL_END
    
    m_comm_sz += bufr.size();
  
  // if this is an evaluation circuit, Gen and Eval
  // must agree on some constant keys
  
  // send constant keys m_gcs[ix].m_const_wire
  GEN_BEGIN
    start = MPI_Wtime();
  
//bufr = get_const_key(m_gcs[ix], 0, 0) + get_const_key(m_gcs[ix], 1, 1);
  bufr = m_gcs[ix].get_const_key(0, 0) + m_gcs[ix].get_const_key(1, 1);
  
  bufr ^= m_prngs[2*ix+0].rand_bits(bufr.size()*8); // encrypt message
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
  if (!m_chks[ix]) // evaluation circuit
    {
      bufr ^= m_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
      std::vector<Bytes> bufr_chunks = bufr.split(Env::key_size_in_bytes());
      m_gcs[ix].set_const_key(0, bufr_chunks[0]);
      m_gcs[ix].set_const_key(1, bufr_chunks[0]);
      //set_const_key(m_gcs[ix], 0, bufr_chunks[0]);
      //set_const_key(m_gcs[ix], 1, bufr_chunks[1]);
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
    if (!m_chks[ix]) { m_gcs[ix].get_gen_decommitments().resize(m_gen_inp_cnt); }
  EVL_END
    
    for (size_t jx = 0; jx < m_gen_inp_cnt; jx++)
      {
        GEN_BEGIN
          start = MPI_Wtime();
        // this bit selects which decommitment to select and send
        //byte bit = m_gen_inp.get_ith_bit(jx) ^ m_gen_inp_masks[ix].get_ith_bit(jx);
        byte bit = m_private_input.get_ith_bit(jx) ^ m_gen_inp_masks[ix].get_ith_bit(jx);
        //bufr = m_gcs[ix].m_gen_inp_decom[2*jx+bit];
        bufr = m_gcs[ix].get_gen_decommitments()[2*jx+bit];
        bufr ^= m_prngs[2*ix+0].rand_bits(bufr.size()*8); // encrypt message
        // remember, even m_prngs are for evaluation circuits
        // and odd m_prngs are for check circuits
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
        if (!m_chks[ix]) // evaluation circuit
          {
            // eval receives decommitment and decrypts it
            // she now has the decommitments to Gen's input keys
            bufr ^= m_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
            //m_gcs[ix].m_gen_inp_decom[jx] = bufr;
            m_gcs[ix].get_gen_decommitments()[jx] = bufr;
          }
        m_timer_evl += MPI_Wtime() - start;
        EVL_END
          
          m_comm_sz += bufr.size();
      }
}


// Gen must send Evl info for check circuits
// he sends the masks for his own input,
// random seeds used to generate:
//     <something, probably the rest of the circuit)>,
// all of his OT keys for Eval's inputs
// and all of his own input decommitments
// Gen always sends these, but Eval can only decrypt if 
// her prng was properly seeded by the OT result
// note: Gen's even-indexed m_prngs are for evaluation circuits,
//       and   odd-indexed  m_prngs are for check circuits
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
  bufr ^= m_prngs[2*ix+1].rand_bits(bufr.size()*8);
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
    start = MPI_Wtime();
  bufr = EVL_RECV();
  m_timer_com += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  if (m_chks[ix]) // check circuit
    {
      bufr ^= m_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
      m_gen_inp_masks[ix] = bufr;
    }
  m_timer_evl += MPI_Wtime() - start;
  EVL_END
    
    m_comm_sz += bufr.size();
  

  // next, Gen sends random seeds to Eval
  
  // send m_rnds[ix]
  GEN_BEGIN
    start = MPI_Wtime();
  bufr = m_rnds[ix];
  bufr ^= m_prngs[2*ix+1].rand_bits(bufr.size()*8); // encrypt message
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
      bufr ^= m_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
      m_rnds[ix] = bufr;
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
        bufr ^= m_prngs[2*ix+1].rand_bits(bufr.size()*8); // encrypt message
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
            bufr ^= m_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
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
        bufr ^= m_prngs[2*ix+1].rand_bits(bufr.size()*8); // encrypt message
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
            bufr ^= m_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
            m_gen_inp_decom[ix][jx] = bufr;
          }
        m_timer_evl += MPI_Wtime() - start;
        EVL_END
          
          m_comm_sz += bufr.size();
      }
}


// in this consistency check, 
// first, Gen and Evl agree on a 2-UHF
// note: in order to ensure Gen's input consistency
//       agreement on the 2-UHF function must happen
//       after Gen commits to his inputs
// then, gen creates his input commitments
// and sends them to Evl, who evaluates them
// finally, Evl saves the input hashes from all of the circuits
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

	start = MPI_Wtime();
		bufr.resize(Env::k()*((m_gen_inp_cnt+7)/8));
        m_timer_evl += MPI_Wtime() - start;
	m_timer_gen += MPI_Wtime() - start;

	start = MPI_Wtime();
                MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm);
	m_timer_mpi += MPI_Wtime() - start;

	start = MPI_Wtime();
        
        // create m rows of length k
        m_matrix = bufr.split(bufr.size()/Env::k());
        
        m_timer_evl += MPI_Wtime() - start;
	m_timer_gen += MPI_Wtime() - start;

        // std::cout << "agree on UHF" << std::endl;

	// now everyone agrees on the UHF given by m_matrix
	for (size_t ix = 0; ix < m_gcs.size(); ix++)
		for (size_t kx = 0; kx < m_matrix.size(); kx++)
	{
		GEN_BEGIN
                  start = MPI_Wtime();
                //gen_next_gen_inp_com(m_gcs[ix], m_matrix[kx], kx);
                m_gcs[ix].gen_next_gen_inp_com(m_matrix[kx], kx);
                bufr = m_gcs[ix].get_and_clear_out_bufr();
                //bufr = get_and_clear_out_bufr(m_gcs[ix]);
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
                                //clear_and_replace_in_bufr(m_gcs[ix], bufr);
                                m_gcs[ix].evl_next_gen_inp_com(m_matrix[kx], kx);
                                //evl_next_gen_inp_com(m_gcs[ix], m_matrix[kx], kx);
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
                  m_gen_inp_hash[ix] = m_gcs[ix].get_gen_inp_hash();//m_gen_inp_hash;
		}
	EVL_END

	step_report("const-check");
}

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
          //gen_init_circuit(m_gcs[ix], m_ot_keys[ix], m_gen_inp_masks[ix], m_rnds[ix]);   
          m_gcs[ix].gen_init_circuit(m_ot_keys[ix], m_gen_inp_masks[ix], m_rnds[ix]);   
            /*
            std::cout << ix << " gen const 0: " << get_const_key(m_gcs[ix], 0, 0).to_hex() << std::endl;
            std::cout << ix << " gen const 1: " << get_const_key(m_gcs[ix], 1, 1).to_hex() << std::endl;
            */
            std::cout << ix << " gen const 0: " << m_gcs[ix].get_const_key( 0, 0).to_hex() << std::endl;
            std::cout << ix << " gen const 1: " << m_gcs[ix].get_const_key( 1, 1).to_hex() << std::endl;
            m_timer_gen += MPI_Wtime() - start;
            GEN_END
              
              EVL_BEGIN
              std::cout << "eval start evaluating" << std::endl;
            
            start = MPI_Wtime();
            if (m_chks[ix]) // check-circuits
              {
                m_gcs[ix].gen_init_circuit(m_ot_keys[ix], m_gen_inp_masks[ix], m_rnds[ix]);
                //gen_init_circuit(m_gcs[ix], m_ot_keys[ix], m_gen_inp_masks[ix], m_rnds[ix]);
                //std::cout << "check ";
              }
            else // evaluation-circuits
              {
                // evl_init_circuit(m_gcs[ix], m_ot_keys[ix], m_gen_inp_masks[ix], m_evl_inp);
                m_gcs[ix].evl_init_circuit(m_ot_keys[ix], m_gen_inp_masks[ix], m_private_input);
                //std::cout << "evl ";
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
            set_callback(m_gcs[ix].m_st, gen_next_gate);
            
            while (get_next_gate(m_gcs[ix].m_st))
              {
                //bufr = get_and_clear_out_bufr(m_gcs[ix]);
                bufr = m_gcs[ix].get_and_clear_out_bufr();
                m_timer_gen += MPI_Wtime() - start;
                
                start = MPI_Wtime();
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
                set_callback(m_gcs[ix].m_st, gen_next_gate);
                //set_callback(m_gcs[ix].m_st, gen_next_gate_m);
                while (get_next_gate(m_gcs[ix].m_st))
                  {
                    
                    bufr = m_gcs[ix].get_and_clear_out_bufr();
                    //bufr = get_and_clear_out_bufr(m_gcs[ix]);
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
                set_callback(m_gcs[ix].m_st, evl_next_gate);
                std::cout<< "enter do loop" << std::endl;
                do {
                  // std::cout << "timing stuff" << std::endl;
                  m_timer_evl += MPI_Wtime() - start;
                  
                  start = MPI_Wtime();
                  bufr = EVL_RECV();
                  m_timer_com += MPI_Wtime() - start;
                  
                  // std::cout << "buffer size add" << std::endl;

                  m_comm_sz += bufr.size();
                  
                  start = MPI_Wtime();
                  m_gcs[ix].clear_and_replace_in_bufr(bufr);
                  //clear_and_replace_in_bufr(m_gcs[ix], bufr);
                  
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
                                  //if (!pass_check(m_gcs[ix]))
					{
						LOG4CXX_FATAL(logger, "Commitment Verification Failure (evaluation circuit)");
						//MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
					}

					if (gen_inp_hash.size() == 0)
					{
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



void BetterYao5::proc_evl_out()
{
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
	reset_timers();

	// TODO: implement Ki08
	m_gen_out = m_gcs[0].get_gen_out();//.m_gen_out;

	EVL_BEGIN
          m_gen_out = m_gcs[0].get_gen_out();//m_gen_out;
        EVL_SEND(m_gen_out);
	EVL_END

	GEN_BEGIN
        m_gen_out = GEN_RECV();
        GEN_END

	step_report("chk-gen-out");
}


//
// Implementation of "Two-Output Secure Computation with Malicious Adversaries"
// by abhi shelat and Chih-hao Shen from EUROCRYPT'11 (Protocol 2)
//
// The evaluator (sender) generates m_ot_bit_cnt pairs of k-bit random strings, and
// the generator (receiver) has input m_ot_bits and will receive output m_ot_out.
//
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

			start = MPI_Wtime();
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
		assert(m_ot_recv_bits.size() >= ((m_ot_bit_cnt+7)/8));

		for (size_t bix = 0; bix < m_ot_bit_cnt; bix++)
		{
			// Step 1: gr=g[b]^r, hr=h[b]^r, where b is the receiver's bit
			start = MPI_Wtime();
				int bit_value = m_ot_recv_bits.get_ith_bit(bix);

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
