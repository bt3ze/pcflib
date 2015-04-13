#ifndef _BETTERYAOGEN_CPP_
#define _BETTERYAOGEN_CPP_

#include "BetterYaoGen.h"
// #include "garbled_circuit.h"
#include <algorithm>

#include <unistd.h>

#include <log4cxx/logger.h>
static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("BetterYaoGen.cpp"));


BetterYaoGen::BetterYaoGen(EnvParams &params) : BetterYao5(params)
{
}


void BetterYaoGen::start()
{
  oblivious_transfer();
  MPI_Barrier(m_mpi_comm);
  cut_and_choose();
  MPI_Barrier(m_mpi_comm);
  cut_and_choose2();
  MPI_Barrier(m_mpi_comm);
  consistency_check();
  MPI_Barrier(m_mpi_comm);
  circuit_evaluate();
  MPI_Barrier(m_mpi_comm);
  final_report();
}


Bytes BetterYaoGen::flip_coins(size_t len_in_bytes)
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

//
// Outputs: m_prngs[]
//
void BetterYaoGen::cut_and_choose2_ot()
{
  double start;
  m_ot_bit_cnt = Env::node_load();
      
  ot_init();
  ot_random();
  
  // Gen's m_ot_out has 2*Env::node_load() seeds and
  // Evl's m_ot_out has   Env::node_load() seeds according to m_chks.

  GEN_BEGIN
    start = MPI_Wtime();
  seed_m_prngs(2*Env::node_load(), m_ot_out);
  /*
  m_prngs.resize(2*Env::node_load());
  for (size_t ix = 0; ix < m_prngs.size(); ix++)
    {
      m_prngs[ix].seed_rand(m_ot_out[ix]);
    }
  */
  m_timer_gen += MPI_Wtime() - start;
  GEN_END
    
    std::cout << "end-cut-n-choose2-ot " << std::endl;
    
    }


//
// Outputs: m_rnds[], m_gen_inp_masks[], m_gcs[].m_gen_inp_decom
//

extern "C" {
void finalize(PCFState *st);
}

void BetterYaoGen::cut_and_choose2_precomputation()
{
  double start;
  
  GEN_BEGIN
    start = MPI_Wtime();
  for (size_t ix = 0; ix < m_gcs.size(); ix++)
    {
      m_rnds[ix] = m_prng.rand_bits(Env::k());
      
      m_gen_inp_masks[ix] = m_prng.rand_bits(m_gen_inp_cnt);
      //m_gen_inp_masks[ix]=m_prgn.rand(Env::k());
      
      gen_init_circuit(m_gcs[ix], m_ot_keys[ix], m_gen_inp_masks[ix], m_rnds[ix]);
      
      m_gcs[ix].m_st = 
        load_pcf_file(Env::pcf_file(), m_gcs[ix].m_const_wire, m_gcs[ix].m_const_wire+1, copy_key);
      m_gcs[ix].m_st->alice_in_size = m_gen_inp_cnt;
      m_gcs[ix].m_st->bob_in_size = m_evl_inp_cnt;
      
      set_external_circuit(m_gcs[ix].m_st, &m_gcs[ix]);
      set_key_copy_function(m_gcs[ix].m_st, copy_key);
      set_key_delete_function(m_gcs[ix].m_st, delete_key);
      set_callback(m_gcs[ix].m_st, gen_next_gate_m);
      
      assert(m_gcs[ix].m_gen_inp_decom.size() <= 2*m_gen_inp_cnt);
      assert(m_gcs[ix].m_gen_inp_decom.size() <= m_gen_inp_cnt);
      //      assert(m_gcs[ix].m_gen_inp_decom.size() == m_gen_inp_cnt);
      std::cout << "decom size: " << m_gcs[ix].m_gen_inp_decom.size() << " \t gen inp count: " << m_gen_inp_cnt << std::endl;
      // the following one fails.
      // so m_gen_inp_decom.size() is less than twice the input count
      //assert(m_gcs[ix].m_gen_inp_decom.size() == 2*m_gen_inp_cnt);
      
        // i am not sure what the purpose of this is.
        // it runs through the commitments
      while ((m_gcs[ix].m_gen_inp_decom.size()/2 < m_gen_inp_cnt) && get_next_gate(m_gcs[ix].m_st))
        {
          get_and_clear_out_bufr(m_gcs[ix]); // discard the garbled gates for now
        }
      assert(m_gcs[ix].m_gen_inp_decom.size() <= 2*m_gen_inp_cnt);
      assert(m_gcs[ix].m_gen_inp_decom.size() <= m_gen_inp_cnt);
      //      assert(m_gcs[ix].m_gen_inp_decom.size() == m_gen_inp_cnt);
      std::cout << "decom size: " << m_gcs[ix].m_gen_inp_decom.size() << " \t gen inp count: " << m_gen_inp_cnt << std::endl;
      
      
//      std::cout << "prngs seeded " << std::endl;

      std::cout << "gen circuits initialized" << std::endl;
      MPI_Barrier(m_mpi_comm);

      finalize(m_gcs[ix].m_st);

      std::cout << "end cut-n-choose2-precomp" << std::endl;

    }
  m_timer_gen += MPI_Wtime() - start;
  GEN_END
}

void BetterYaoGen::cut_and_choose2_evl_circuit(size_t ix)
{
  double start;
  
  Bytes bufr;
  
  // send masked gen inp
  GEN_BEGIN
    start = MPI_Wtime();
  bufr = m_gen_inp_masks[ix] ^ m_gen_inp;
  bufr ^= m_prngs[2*ix+0].rand_bits(bufr.size()*8); // encrypt message
  m_timer_gen += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  GEN_SEND(bufr);
  m_timer_com += MPI_Wtime() - start;
  GEN_END
    
    m_comm_sz += bufr.size();
  
  // send constant keys m_gcs[ix].m_const_wire
  GEN_BEGIN
    start = MPI_Wtime();
  bufr = get_const_key(m_gcs[ix], 0, 0) + get_const_key(m_gcs[ix], 1, 1);
  bufr ^= m_prngs[2*ix+0].rand_bits(bufr.size()*8); // encrypt message
  m_timer_gen += MPI_Wtime() - start;
  
  start = MPI_Wtime();
  GEN_SEND(bufr);
  m_timer_com += MPI_Wtime() - start;
  GEN_END
    
    m_comm_sz += bufr.size();
  
  // send m_gcs[ix].m_gen_inp_decom
  GEN_BEGIN
    assert(m_gcs[ix].m_gen_inp_decom.size() == 2*m_gen_inp_cnt);
  GEN_END
    
    
    for (size_t jx = 0; jx < m_gen_inp_cnt; jx++)
      {
        GEN_BEGIN
          start = MPI_Wtime();
        byte bit = m_gen_inp.get_ith_bit(jx) ^ m_gen_inp_masks[ix].get_ith_bit(jx);
        bufr = m_gcs[ix].m_gen_inp_decom[2*jx+bit];
        bufr ^= m_prngs[2*ix+0].rand_bits(bufr.size()*8); // encrypt message
        m_timer_gen += MPI_Wtime() - start;
        
        start = MPI_Wtime();
        GEN_SEND(bufr);
        m_timer_com += MPI_Wtime() - start;
        GEN_END
         
          m_comm_sz += bufr.size();
      }
}


void BetterYaoGen::cut_and_choose2_chk_circuit(size_t ix)
{
	double start;

	Bytes bufr;
        std::vector<Bytes> bufr_chunks;

	// send m_gen_inp_masks[ix]
	GEN_BEGIN
		start = MPI_Wtime();
			bufr = m_gen_inp_masks[ix];
			bufr ^= m_prngs[2*ix+1].rand_bits(bufr.size()*8); // encrypt message
		m_timer_gen += MPI_Wtime() - start;

		start = MPI_Wtime();
			GEN_SEND(bufr);
		m_timer_com += MPI_Wtime() - start;
	GEN_END

	m_comm_sz += bufr.size();

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

	m_comm_sz += bufr.size();


	// send m_ot_kesy[ix]
	GEN_BEGIN
		assert(m_ot_keys[ix].size() == 2*m_evl_inp_cnt);
	GEN_END

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

		m_comm_sz += bufr.size();
	}

	// send m_gcs[ix].m_gen_inp_decom
	GEN_BEGIN
		assert(m_gcs[ix].m_gen_inp_decom.size() == 2*m_gen_inp_cnt);
	GEN_END

	for (size_t jx = 0; jx < 2*m_gen_inp_cnt; jx++)
	{
		GEN_BEGIN
			start = MPI_Wtime();
				bufr = m_gcs[ix].m_gen_inp_decom[jx];
				bufr ^= m_prngs[2*ix+1].rand_bits(bufr.size()*8); // encrypt message
			m_timer_gen += MPI_Wtime() - start;

			start = MPI_Wtime();
				GEN_SEND(bufr);
			m_timer_com += MPI_Wtime() - start;
		GEN_END

		m_comm_sz += bufr.size();
	}
}


void BetterYaoGen::consistency_check()
{
  // std::cout << "const check start" << std::endl;

	reset_timers();

	Bytes bufr;

	double start;

	// jointly pick a 2-UHF matrix
	bufr = flip_coins(Env::k()*fit_to_byte_containers(m_gen_inp_cnt));
         // only roots get the result

	start = MPI_Wtime();
        bufr.resize(Env::k()*fit_to_byte_containers(m_gen_inp_cnt));
        
                //   m_timer_evl += MPI_Wtime() - start;
	m_timer_gen += MPI_Wtime() - start;

	start = MPI_Wtime();
        MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm);
	m_timer_mpi += MPI_Wtime() - start;
        
	start = MPI_Wtime();
        m_matrix = bufr.split(bufr.size()/Env::k());
                //m_timer_evl += MPI_Wtime() - start;
	m_timer_gen += MPI_Wtime() - start;

        // std::cout << "agree on UHF" << std::endl;

	// now everyone agrees on the UHF given by m_matrix
	for (size_t ix = 0; ix < m_gcs.size(); ix++)
          {
            for (size_t kx = 0; kx < m_matrix.size(); kx++)
              {
		GEN_BEGIN
                  start = MPI_Wtime();
                gen_next_gen_inp_com(m_gcs[ix], m_matrix[kx], kx);
                bufr = get_and_clear_out_bufr(m_gcs[ix]);
                m_timer_gen += MPI_Wtime() - start;
                
                start = MPI_Wtime();
                GEN_SEND(bufr);
                m_timer_com += MPI_Wtime() - start;
                
                GEN_END
                  
                  m_comm_sz += bufr.size();
              }
          }
        // std::cout << "EVL check hashes" << std::endl;

	step_report("const-check");
}

void BetterYaoGen::circuit_evaluate()
{
	reset_timers();

	double start;

	int verify = 1;
	Bytes bufr;

        std::cout << "begin circuit evaluate" << std::endl;


        MPI_Barrier(m_mpi_comm);
	for (size_t ix = 0; ix < m_gcs.size(); ix++)
	{
          GEN_BEGIN
            start = MPI_Wtime();
            gen_init_circuit(m_gcs[ix], m_ot_keys[ix], m_gen_inp_masks[ix], m_rnds[ix]);
          
            std::cout << ix << " gen const 0: " << get_const_key(m_gcs[ix], 0, 0).to_hex() << std::endl;
            std::cout << ix << " gen const 1: " << get_const_key(m_gcs[ix], 1, 1).to_hex() << std::endl;
            m_timer_gen += MPI_Wtime() - start;
          GEN_END
                            
            std::cout << "load pcf file" << std:: endl;
          start = MPI_Wtime();
          m_gcs[ix].m_st = 
            load_pcf_file(Env::pcf_file(), m_gcs[ix].m_const_wire, m_gcs[ix].m_const_wire+1, copy_key);
          m_gcs[ix].m_st->alice_in_size = m_gen_inp_cnt;
          m_gcs[ix].m_st->bob_in_size = m_evl_inp_cnt;
          
          
          set_external_circuit(m_gcs[ix].m_st, &m_gcs[ix]);
          set_key_copy_function(m_gcs[ix].m_st, copy_key);
          set_key_delete_function(m_gcs[ix].m_st, delete_key);
          m_timer_gen += MPI_Wtime() - start;
          //          m_timer_evl += MPI_Wtime() - start;
          
          GEN_BEGIN // generate and send the circuit gate-by-gate
            
            std::cout << "gen generate and send" << std::endl;
          start = MPI_Wtime();
          set_callback(m_gcs[ix].m_st, gen_next_gate_m);
          while (get_next_gate(m_gcs[ix].m_st))
            {
              bufr = get_and_clear_out_bufr(m_gcs[ix]);
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
            
            }

        //fprintf(stderr,"end gen circuit eval\n");
        std::cout << "end gen circuit eval" << std::endl;
        MPI_Barrier(m_mpi_comm);
                
	step_report("circuit-evl");
        
	if (m_gcs[0].m_evl_out_ix != 0)
          proc_evl_out();
        
        if (m_gcs[0].m_gen_out_ix != 0)
          proc_gen_out();
}



void BetterYaoGen::proc_evl_out()
{
  // this is entirely an evl function. gen doesn't output eval's inputs?
  // gen should never learn Evl's inputs
}

void BetterYaoGen::proc_gen_out()
{
	reset_timers();

	// TODO: implement Ki08
	m_gen_out = m_gcs[0].m_gen_out;

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
void BetterYaoGen::ot_init()
{
	double start;

	start = MPI_Wtime();
        std::vector<Bytes> bufr_chunks;
        Bytes bufr(Env::elm_size_in_bytes()*4);
        
        Z y, a;
	m_timer_gen += MPI_Wtime() - start;
        //	m_timer_evl += MPI_Wtime() - start;

	// step 1: ZKPoK of the CRS: g[0], h[0], g[1], h[1]
	if (Env::is_root())
	{
          
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


void BetterYaoGen::ot_random()
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
	// m_timer_evl += MPI_Wtime() - start;

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


#endif
