#ifndef _BETTERYAOEVL_CPP_
#define _BETTERYAOEVL_CPP_

#include "YaoBase.h"
#include "BetterYao5.h"
#include "BetterYaoEvl.h"
// #include "garbled_circuit.h"
#include <algorithm>

#include <unistd.h>

#include <log4cxx/logger.h>
static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("BetterYaoEvl.cpp"));


BetterYaoEvl::BetterYaoEvl(EnvParams &params) : BetterYao5(params)
{
}


void BetterYaoEvl::start()
{
  oblivious_transfer();
  //  MPI_Barrier(m_mpi_comm);
  cut_and_choose();
  //  MPI_Barrier(m_mpi_comm);
  cut_and_choose2();
  //  MPI_Barrier(m_mpi_comm);
  consistency_check();
  //  MPI_Barrier(m_mpi_comm);
  circuit_evaluate();
  //  MPI_Barrier(m_mpi_comm);
  final_report();
}


Bytes BetterYaoEvl::flip_coins(size_t len_in_bytes)
{
  double start;
  
  Bytes bufr;
  
  //  if (Env::is_root())
  //  {
      Bytes remote_coins, commitment, commit_value;
      
      start = MPI_Wtime();
      Bytes coins = m_prng.rand_bits(len_in_bytes*8);	// Step 0: flip coins
      m_timer_gen += MPI_Wtime() - start;
      m_timer_evl += MPI_Wtime() - start;
      
      
      //      EVL_BEGIN
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
      // EVL_END
      
      m_comm_sz = commitment.size() + remote_coins.size() + commit_value.size();
      
      start = MPI_Wtime();
      
      coins ^= remote_coins;
      // combine randomnesses from both players
      
      bufr.swap(coins);
      //bufr = coins;
      
      m_timer_evl += MPI_Wtime() - start;
      //m_timer_gen += MPI_Wtime() - start;
      //}
  
  return bufr;
}

//
// Outputs: m_prngs[]
//
void BetterYaoEvl::cut_and_choose2_ot()
{
  double start;
  m_ot_bit_cnt = Env::node_load();
  
  EVL_BEGIN
    start = MPI_Wtime();
  m_ot_recv_bits.resize(fit_to_byte_containers(m_ot_bit_cnt));
  
  for (size_t ix = 0; ix < m_chks.size(); ix++)
    {
      m_ot_recv_bits.set_ith_bit(ix, m_chks[ix]);
    }
  m_timer_evl += MPI_Wtime() - start;
  EVL_END
    
  ot_init();
  ot_random();
  
  // Gen's m_ot_out has 2*Env::node_load() seeds and
  // Evl's m_ot_out has   Env::node_load() seeds according to m_chks.
  
  EVL_BEGIN
    start = MPI_Wtime();
  seed_m_prngs(Env::node_load(), m_ot_out);
  /*
  m_prngs.resize(Env::node_load());
  for (size_t ix = 0; ix < m_prngs.size(); ix++)
    {
      m_prngs[ix].seed_rand(m_ot_out[ix]);
    }
  */
  m_timer_evl += MPI_Wtime() - start;
  EVL_END
    
    }

//
// Outputs: m_rnds[], m_gen_inp_masks[], m_gcs[].m_gen_inp_decom
//

extern "C" {
void finalize(PCFState *st);
}

void BetterYaoEvl::cut_and_choose2_precomputation()
{
  // this function empty for eval
  
}

void BetterYaoEvl::cut_and_choose2_evl_circuit(size_t ix)
{
  double start;
  
  Bytes bufr;
  
  // send masked gen inp
    
  EVL_BEGIN
    start = MPI_Wtime();
  bufr = EVL_RECV();
  m_timer_com += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  if (!m_chks[ix]) // evaluation circuit
    {
      bufr ^= m_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
      m_gen_inp_masks[ix] = bufr;
    }
  m_timer_evl += MPI_Wtime() - start;
  EVL_END
    
    m_comm_sz += bufr.size();
  
  EVL_BEGIN
    start = MPI_Wtime();
  bufr = EVL_RECV();
  m_timer_com += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  if (!m_chks[ix]) // evaluation circuit
    {
      bufr ^= m_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
      std::vector<Bytes> bufr_chunks = bufr.split(Env::key_size_in_bytes());
      set_const_key(m_gcs[ix], 0, bufr_chunks[0]);
      set_const_key(m_gcs[ix], 1, bufr_chunks[1]);
    }
  m_timer_evl += MPI_Wtime() - start;
  EVL_END
    
    m_comm_sz += bufr.size();
  
    
  EVL_BEGIN
    if (!m_chks[ix]) { m_gcs[ix].m_gen_inp_decom.resize(m_gen_inp_cnt); }
  EVL_END
    
    for (size_t jx = 0; jx < m_gen_inp_cnt; jx++)
      {
        
        EVL_BEGIN
          start = MPI_Wtime();
        bufr = EVL_RECV();
        m_timer_com += MPI_Wtime() - start;
        
        start = MPI_Wtime();
        if (!m_chks[ix]) // evaluation circuit
          {
            bufr ^= m_prngs[ix].rand_bits(bufr.size()*8); // decrypt message
            m_gcs[ix].m_gen_inp_decom[jx] = bufr;
          }
        m_timer_evl += MPI_Wtime() - start;
        EVL_END
          
          m_comm_sz += bufr.size();
      }
}


void BetterYaoEvl::cut_and_choose2_chk_circuit(size_t ix)
{
	double start;

	Bytes bufr;
        std::vector<Bytes> bufr_chunks;

	// send m_gen_inp_masks[ix]

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



	EVL_BEGIN
		if (m_chks[ix]) { m_ot_keys[ix].resize(2*m_evl_inp_cnt); }
	EVL_END

	for (size_t jx = 0; jx < 2*m_evl_inp_cnt; jx++)
	{

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
        
	
	EVL_BEGIN
          if (m_chks[ix]) { m_gen_inp_decom[ix].resize(2*m_gen_inp_cnt); }
	EVL_END
          
	for (size_t jx = 0; jx < 2*m_gen_inp_cnt; jx++)
	{
          
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


void BetterYaoEvl::consistency_check()
{
  // std::cout << "const check start" << std::endl;
  
  reset_timers();
  
  Bytes bufr;
  
  double start;
  
  // jointly pick a 2-UHF matrix
  if(Env::is_root()){
    bufr = flip_coins(Env::k()*fit_to_byte_containers(m_gen_inp_cnt));
    // only roots get the result
  }
  
  start = MPI_Wtime();
  bufr.resize(Env::k()*fit_to_byte_containers(m_gen_inp_cnt));
  
  m_timer_evl += MPI_Wtime() - start;
  //	m_timer_gen += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm);
  m_timer_mpi += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  m_matrix = bufr.split(bufr.size()/Env::k());
  m_timer_evl += MPI_Wtime() - start;
  //m_timer_gen += MPI_Wtime() - start;
  
  // std::cout << "agree on UHF" << std::endl;
  
  // now everyone agrees on the UHF given by m_matrix
  for (size_t ix = 0; ix < m_gcs.size(); ix++)
    {
      for (size_t kx = 0; kx < m_matrix.size(); kx++)
        {
          
          EVL_BEGIN
            start = MPI_Wtime();
          bufr = EVL_RECV();
          m_timer_com += MPI_Wtime() - start;
          
          if (!m_chks[ix]) // evaluation circuit
            {
              start = MPI_Wtime();
              clear_and_replace_in_bufr(m_gcs[ix], bufr);
              evl_next_gen_inp_com(m_gcs[ix], m_matrix[kx], kx);
              m_timer_evl += MPI_Wtime() - start;
            }
          EVL_END
            
            m_comm_sz += bufr.size();
        }
    }
  
  // std::cout << "EVL check hashes" << std::endl;
  
  EVL_BEGIN
    for (size_t ix = 0; ix < m_gcs.size(); ix++)
      if (!m_chks[ix])
        {
          m_gen_inp_hash[ix] = m_gcs[ix].m_gen_inp_hash;
        }
  EVL_END
    
    step_report("const-check");
}

void BetterYaoEvl::circuit_evaluate()
{
	reset_timers();

	double start;

	int verify = 1;
	Bytes bufr;

        std::cout << "begin circuit evaluate" << std::endl;

        MPI_Barrier(m_mpi_comm);

	for (size_t ix = 0; ix < m_gcs.size(); ix++)
	{
          
          EVL_BEGIN
            std::cout << "eval start evaluating" << std::endl;
          
          start = MPI_Wtime();
          if (m_chks[ix]) // check-circuits
            {
              gen_init_circuit(m_gcs[ix], m_ot_keys[ix], m_gen_inp_masks[ix], m_rnds[ix]);
              //std::cout << "check ";
            }
          else // evaluation-circuits
            {
              evl_init_circuit(m_gcs[ix], m_ot_keys[ix], m_gen_inp_masks[ix], m_evl_inp);
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
            
            
            set_external_circuit(m_gcs[ix].m_st, &m_gcs[ix]);
            set_key_copy_function(m_gcs[ix].m_st, copy_key);
            set_key_delete_function(m_gcs[ix].m_st, delete_key);
            //            m_timer_gen += MPI_Wtime() - start;
            m_timer_evl += MPI_Wtime() - start;
            
            EVL_BEGIN // receive and evaluate the circuit gate-by-gate
              
              std::cout << "eval receive and evaluate" << std::endl;
            
            if (m_chks[ix]) // check circuit
              {
                std::cout << "eval check circuit" << m_chks[ix] << std::endl;
                start = MPI_Wtime();
                set_callback(m_gcs[ix].m_st, gen_next_gate_m);
                while (get_next_gate(m_gcs[ix].m_st))
                  {
                    
                    bufr = get_and_clear_out_bufr(m_gcs[ix]);
                    m_timer_evl += MPI_Wtime() - start;
                    
                    start = MPI_Wtime();
                    Bytes recv = EVL_RECV();
                    m_timer_com += MPI_Wtime() - start;
                    
                    m_comm_sz += bufr.size();
                    
                    start = MPI_Wtime(); // start m_timer_evl
                    verify &= (bufr == recv);
                    //FLAG: what's up with this check?
                    // how will bufr = recv. what are we doing to bufr other than getting the out bufr from the circuit?
                    // must check what happens in that method, and see how m_gcs has been altered
                    // out buffer will be set by the gen_next_gate_m callback function
                  }
                m_timer_gen += MPI_Wtime() - start;
                
                EVL_RECV(); // a redundant value to prevent the evlauator from hanging
           
                }
                
                  else // evaluation circuit
              {
                std::cout << "eval evaluate circuit" << std::endl;
                start = MPI_Wtime();
                set_callback(m_gcs[ix].m_st, evl_next_gate_m);
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
                  clear_and_replace_in_bufr(m_gcs[ix], bufr);
                  
                  //std::cout << "received" << std::endl;
                  //fprintf(stderr, "received");
                  
                  //std::cout << "got a gate" << std::endl;                  
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
                    if (!(m_gen_inp_decom[ix][jx] == m_gcs[ix].m_gen_inp_decom[jx]))
                      {
                        LOG4CXX_FATAL(logger, "Commitment Verification Failure (check circuit)");
                        //MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                      }
                  }
              }
            else // evaluation circuit
              {
                if (!pass_check(m_gcs[ix]))
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
            
            trim_output(m_gcs[ix]);
          }
        m_timer_evl += MPI_Wtime() - start;
	EVL_END
          
          MPI_Barrier(m_mpi_comm);
          
          step_report("circuit-evl");

	if (m_gcs[0].m_evl_out_ix != 0)
          proc_evl_out();
        
        if (m_gcs[0].m_gen_out_ix != 0)
          proc_gen_out();
}



void BetterYaoEvl::proc_evl_out()
{
  EVL_BEGIN
    reset_timers();
  
  double start;
  Bytes send, recv;
  
  start = MPI_Wtime();
  for (size_t ix = 0; ix < m_gcs.size(); ix++) // fill zeros for uniformity (convenient for MPIs)
    {
      send += m_gcs[ix].m_evl_out;
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

void BetterYaoEvl::proc_gen_out()
{
  reset_timers();

  // TODO: implement Ki08
  m_gen_out = m_gcs[0].m_gen_out;

  EVL_BEGIN
    m_gen_out = m_gcs[0].m_gen_out;
  EVL_SEND(m_gen_out);
  EVL_END

    step_report("chk-gen-out");
}


//
// Implementation of "Two-Output Secure Computation with Malicious Adversaries"
// by abhi shelat and Chih-hao Shen from EUROCRYPT'11 (Protocol 2)
//
// The evaluator (sender) generates m_ot_bit_cnt pairs of k-bit random strings, and
// the generator (receiver) has input m_ot_bits and will receive output m_ot_out.
//
void BetterYaoEvl::ot_init()
{
  double start;

  start = MPI_Wtime();
  std::vector<Bytes> bufr_chunks;
  Bytes bufr(Env::elm_size_in_bytes()*4);

  Z y, a;
  // m_timer_gen += MPI_Wtime() - start;
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


void BetterYaoEvl::ot_random()
{
  double start;

  start = MPI_Wtime();
  Bytes send, recv;
  std::vector<Bytes> recv_chunks;

  Z r, s[2], t[2];
  G gr, hr, X[2], Y[2];

  m_ot_out.clear();
  m_ot_out.reserve(2*m_ot_bit_cnt); // the receiver only uses half of it
  // m_timer_gen += MPI_Wtime() - start;
  m_timer_evl += MPI_Wtime() - start;

  EVL_BEGIN // evaluator (OT receiver)
    assert(m_ot_recv_bits.size() >= fit_to_byte_containers(m_ot_bit_cnt));
  
        
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
          
    }

#endif