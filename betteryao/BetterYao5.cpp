#ifndef BETTERYAO5_CPP_
#define BETTERYAO5_CPP_

#include "BetterYao5.h"
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
      initialize_circuit_mal(m_gcs[ix]);
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


void BetterYao5::cut_and_choose()
{
  reset_timers();
  
  double start;
  
  Bytes coins = flip_coins(Env::key_size_in_bytes()); // only roots get the result
  // this is the collaborative part of the cut-and-choose
  
  if (Env::is_root())
    {
      Prng prng;
      start = MPI_Wtime();
      prng.seed_rand(coins); // use the coins to generate more random bits
      
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
  
  MPI_Barrier(m_mpi_comm);

  start = MPI_Wtime();
  m_chks.resize(Env::node_load());
  m_timer_evl += MPI_Wtime() - start;
  m_timer_gen += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  MPI_Scatter(&m_all_chks[0], m_chks.size(), MPI_BYTE, &m_chks[0], m_chks.size(), MPI_BYTE, 0, m_mpi_comm); // distributes m_all_chks to all subprocesses
  // because the prng is seeded with the same (collaboratively tossed) coins, these arrays will be equivalent
  m_timer_mpi += MPI_Wtime() - start;
  
  step_report("cut-&-check");
}

void BetterYao5::cut_and_choose2()
{
  reset_timers();
  
  cut_and_choose2_ot();
  cut_and_choose2_precomputation(); // this is empty for Evl, includes things for Gen
  
  for (size_t ix = 0; ix < m_gcs.size(); ix++)
    {
      cut_and_choose2_evl_circuit(ix);
      cut_and_choose2_chk_circuit(ix);
    }
  
  step_report("cut-'n-chk2");
}

void BetterYao5::seed_m_prngs(size_t num_prngs, std::vector<Bytes> seeds){
  assert(seeds.size() == num_prngs); // this actually not really necessary if assertion is programmatically true. could just use the size of the seeds.
  
  m_prngs.resize(num_prngs);
  for(size_t ix = 0; ix < num_prngs; ix++){
    m_prngs[ix].seed_rand(seeds[ix]);
  }

  std::cout << "prngs seeded " << std::endl;

}



#endif /* BETTERYAO5_CPP_ */
