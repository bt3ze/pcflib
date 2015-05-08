#ifndef YAO_GEN_CPP_
#define YAO_GEN_CPP_

#include "YaoGen.h"

#include <log4cxx/logger.h>
static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("YaoGen.cpp"));


Yao::Yao(EnvParams &params) : YaoBase(params), m_gcs(0)
{
	if (Env::s() != 1)
	{
		LOG4CXX_FATAL(logger, "s has to be 1 in the honest-but-curious setting");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

	m_rnds.resize(1);
	m_gen_inp_masks.resize(1);
	m_gcs.resize(1);

        get_and_size_inputs();

	std::cout << "gencount: "<< m_gen_inp_cnt <<"\n";
	std::cout << "evlcount: "<< m_evl_inp_cnt <<"\n";
}


void Yao::start()
{
	oblivious_transfer();
	circuit_evaluate();
	final_report();
}


void Yao::circuit_evaluate()
{
	reset_timers();

	double start;

	Bytes bufr;
        std::vector<Bytes> bufr_chunks;

	G M;

        // step 1:
        // gen generates some random bit and random input mask
        // gen sends masked input to evl
        // gen initializes circuit with 
	GEN_BEGIN
		start = MPI_Wtime();
			m_rnds[0] = m_prng.rand_bits(Env::k());
			m_gen_inp_masks[0] = m_prng.rand_bits(m_gen_inp_cnt);
                     
		m_timer_gen += MPI_Wtime() - start;

		start = MPI_Wtime();
			GEN_SEND(m_gen_inp_masks[0] ^ m_gen_inp); // send the masked gen_inp
		m_timer_com += MPI_Wtime() - start;

		start = MPI_Wtime();
			gen_init_circuit(m_gcs[0], m_ot_keys[0], m_gen_inp_masks[0], m_rnds[0]);
			m_gcs[0].m_gen_inp = m_gen_inp;
		m_timer_gen += MPI_Wtime() - start;
	GEN_END


	m_comm_sz += m_gen_inp_masks[0].size();

	GEN_BEGIN
		for (size_t ix = 0; ix < 2; ix++)
		{
			start = MPI_Wtime();
				bufr = get_const_key(m_gcs[0], ix, ix);
			m_timer_gen += MPI_Wtime() - start;

			start = MPI_Wtime();
				GEN_SEND(bufr);	
			m_timer_com += MPI_Wtime() - start;
		
			m_comm_sz += bufr.size();
		}
	GEN_END

	m_gcs[0].m_st = load_pcf_file(Env::pcf_file(), m_gcs[0].m_const_wire, m_gcs[0].m_const_wire+1, copy_key);
        m_gcs[0].m_st->alice_in_size = m_gen_inp_cnt;
        m_gcs[0].m_st->bob_in_size = m_evl_inp_cnt;

	set_external_circuit(m_gcs[0].m_st, &m_gcs[0]);
	set_key_copy_function(m_gcs[0].m_st, copy_key);
	set_key_delete_function(m_gcs[0].m_st, delete_key);


	step_report("pre-cir-evl");
	reset_timers();

	GEN_BEGIN // generate and send the circuit gate-by-gate
		set_callback(m_gcs[0].m_st, gen_next_gate);
		start = MPI_Wtime();
			while (get_next_gate(m_gcs[0].m_st))
			{
                          bufr = get_and_clear_out_bufr(m_gcs[0]);
                          m_timer_gen += MPI_Wtime() - start;

                          start = MPI_Wtime();
                          //assert(bufr.size() > 0);
                          GEN_SEND(bufr);
                          m_timer_com += MPI_Wtime() - start;

                          m_comm_sz += bufr.size();

                          start = MPI_Wtime(); // start m_timer_gen
			}
		m_timer_gen += MPI_Wtime() - start;

		GEN_SEND(Bytes(0)); // a redundant value to prevent the evlauator from hanging
	GEN_END

	std::cout << "gencount: "<<m_gcs[0].m_gen_inp_ix<<"\n";
	std::cout << "evlcount: "<<m_gcs[0].m_evl_inp_ix<<"\n";
	
	step_report("circuit-evl");

	trim_output(m_gcs[0]);

	if (m_gcs[0].m_evl_out_ix != 0)
		proc_evl_out();

	if (m_gcs[0].m_gen_out_ix != 0)
		proc_gen_out();
}


void Yao::proc_evl_out()
{
	reset_timers();

        // this is really an Evl method
        
	step_report("chk-evl-out");
}

void Yao::proc_gen_out()
{
	reset_timers();

	double start;

	GEN_BEGIN
		start = MPI_Wtime();
			m_gen_out = GEN_RECV();
		m_timer_com += MPI_Wtime() - start;
	GEN_END

	m_comm_sz += m_gen_out.size();

	step_report("chk-gen-out");
}


#endif
