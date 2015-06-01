#ifndef YAOBASE_H_
#define YAOBASE_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "mpi.h"

#include "Env.h"
#include "NetIO.h"

#include "../pcflib.h"

#ifdef GEN_CODE

	// User mode (code for the generator)
	#define GEN_BEGIN
	#define GEN_END
	#define EVL_BEGIN     if (0) {
	#define EVL_END       }
	#define GEN_SEND(d)   Env::remote()->write_bytes(d)
	#define EVL_RECV()    Env::remote()->read_bytes()
	#define EVL_SEND(d)   Env::remote()->write_bytes(d)
	#define GEN_RECV()    Env::remote()->read_bytes()

#elif defined EVL_CODE

	// User mode (code for the evaluator)
	#define GEN_BEGIN     if (0) {
	#define GEN_END       }
	#define EVL_BEGIN
	#define EVL_END
	#define GEN_SEND(d)   Env::remote()->write_bytes(d)
	#define EVL_RECV()    Env::remote()->read_bytes()
	#define EVL_SEND(d)   Env::remote()->write_bytes(d)
	#define GEN_RECV()    Env::remote()->read_bytes()

#else

	// Simulation mode
	#define GEN_BEGIN     if (!Env::is_evl()) {
	#define GEN_END       }
	#define GEN_SEND(d)   send_data(Env::world_rank()+1, (d))
	#define GEN_RECV()    recv_data(Env::world_rank()+1)
	#define EVL_BEGIN     if ( Env::is_evl()) {
	#define EVL_END       }
	#define EVL_SEND(d)   send_data(Env::world_rank()-1, (d))
	#define EVL_RECV()    recv_data(Env::world_rank()-1)

#endif


class YaoBase {
public:
        YaoBase(EnvParams &params);
        virtual ~YaoBase();

        virtual void start() = 0;
        

private:
	void init_cluster(EnvParams &params);
	void init_network(EnvParams &params);
	void init_environ(EnvParams &params);
	void init_private(EnvParams &params);

protected:
	// subroutines for the communication in the Simulation mode
	Bytes recv_data(int src_node);
	void send_data(int dst_node, const Bytes &data);

	// subroutines for profiling
	void reset_timers();
	void step_report(std::string step_name);
	void step_report_no_sync(std::string step_name);
	void final_report();

        // subroutines for the protocol
        void oblivious_transfer();
        void get_and_size_inputs();
        void size_inputs();
        void get_inputs(); 

protected:
	// variables for MPI
	MPI_Comm            m_mpi_comm;

	// variables for profiling
	double              m_timer_gen;
	double              m_timer_evl;
	double              m_timer_com; // inter-cluster communication
	double              m_timer_mpi; // intra-cluster communication

	uint64_t            m_comm_sz;

        std::vector<double>      m_timer_cmp_vec;
        std::vector<double>      m_timer_mpi_vec;
        std::vector<double>      m_timer_cmm_vec;

        std::vector<std::string> m_step_name_vec;
        std::vector<uint64_t>    m_comm_sz_vec;

	// variables for Yao protocol
	Bytes               m_gen_out;
        Bytes               m_evl_out;

        Bytes               m_private_input;
        

	Prng                m_prng;

        // variables for IKNP03 OT-extension implementation
        // or SS11 committing OT implementation
	G                               m_ot_g[2];
	G                               m_ot_h[2];
        std::vector<std::vector<Bytes> > m_ot_keys; // ot output

        // variables for input counts
	uint32_t               m_gen_inp_cnt;
	uint32_t               m_evl_inp_cnt;

};

#endif
