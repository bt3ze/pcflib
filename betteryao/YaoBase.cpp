#include <unistd.h>
#include <iomanip>

#include <netdb.h>
#include <arpa/inet.h>
#include <mpi.h>
#include <cstring>
#include <vector>

#include "YaoBase.h"

#include <log4cxx/logger.h>
static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("YaoBase.cpp"));

using std::vector;


inline std::string get_IP()
{
	// get local IP and display
	char hostname[1024];
	gethostname(hostname, 1024);
	struct hostent *host = gethostbyname(hostname);
	const std::string local_IP = inet_ntoa(*((struct in_addr *)host->h_addr_list[0]));
	return local_IP;
}


YaoBase::YaoBase(EnvParams &params)
{

  LOG4CXX_INFO(logger, "Begin YaoBase");

	init_cluster(params);
#if defined EVL_CODE || defined GEN_CODE
	init_network(params); // no need in simulation mode
#endif
	init_environ(params);
	init_private(params);

	// [WARNING] no access to params after this point. use Env instead.

	// display local IP
	EVL_BEGIN
		LOG4CXX_INFO(logger, "EVL (" << Env::group_rank() << ") is at " << get_IP());
	EVL_END
xo
	GEN_BEGIN
		LOG4CXX_INFO(logger, "GEN (" << Env::group_rank() << ") is at " << get_IP());
	GEN_END

	if (!Env::is_root())
		return;

	LOG4CXX_INFO(logger, "========================================================");
	LOG4CXX_INFO(logger, "Starting Yao protocol");
	LOG4CXX_INFO(logger, "========================================================");
}

void YaoBase::init_cluster(EnvParams &params)
{
	// inquire world info
	MPI_Comm_rank(MPI_COMM_WORLD, &params.wrld_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &params.node_amnt);

#if defined GEN_CODE || defined EVL_CODE
	params.node_rank = params.wrld_rank;
	m_mpi_comm = MPI_COMM_WORLD;
#else
	if (params.node_amnt % 2 != 0)
	{
		LOG4CXX_FATAL(logger, "statistical parameter s needs to be an even number in the Simulation mode");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

	// divide the world into two groups of the same size (Simulation mode)
	MPI_Comm_split(MPI_COMM_WORLD, params.wrld_rank % 2, params.wrld_rank, &m_mpi_comm);

	// inquire info of the new group
	MPI_Comm_rank(m_mpi_comm, &params.node_rank);
	MPI_Comm_size(m_mpi_comm, &params.node_amnt);
#endif
}


void YaoBase::init_network(EnvParams &params)
{
#if !(defined GEN_CODE || defined EVL_CODE) // just in case (shouldn't reach this point to begin with)
	return;
#endif

	const int IP_SERVER_PORT = params.port_base;
	const int PORT = params.port_base + params.node_rank+1;
	Bytes send, recv;

	// get local IP
	char hostname[1024];
	gethostname(hostname, 1024);
	struct hostent *host = gethostbyname(hostname);
	const std::string local_ip = inet_ntoa(*((struct in_addr *)host->h_addr_list[0]));

	EVL_BEGIN
		// collect IPs from slaves and send them the evaluator via IP server
		send.resize(sizeof(struct in_addr));
		memcpy(&send[0], host->h_addr_list[0], send.size());

		recv.resize(sizeof(struct in_addr)*params.node_amnt); // only used by node 0
		MPI_Gather(&send[0], send.size(), MPI_BYTE, &recv[0], send.size(), MPI_BYTE, 0, m_mpi_comm);

		if (params.node_rank == 0)
		{
			//ServerSocket ip_exchanger(Env::IP_SERVER_PORT);
			ServerSocket ip_exchanger(IP_SERVER_PORT);
			Socket *sock = ip_exchanger.accept();
			sock->write_bytes(recv); // send slaves' IPs to remote
		}

		LOG4CXX_INFO(logger, "EVL (" << params.node_rank << ":" << local_ip << ") is listening at port " << PORT);
		params.server = new ServerSocket(PORT);
		params.remote = params.server->accept();
		LOG4CXX_INFO(logger, "EVL (" << params.node_rank << ":" << local_ip << ") is connected at port " << PORT);
	EVL_END

	GEN_BEGIN
		// receive IPs from the generator via IP server and forward them to slaves
		if (params.node_rank == 0)
		{
			//ClientSocket ip_exchanger(params.ipserve_addr, Env::IP_SERVER_PORT);
			ClientSocket ip_exchanger(params.ipserve_addr, IP_SERVER_PORT);
			send = ip_exchanger.read_bytes(); // receive the generator's slaves' IPs
		}

		recv.resize(sizeof(struct in_addr));
		MPI_Scatter(&send[0], recv.size(), MPI_BYTE, &recv[0], recv.size(), MPI_BYTE, 0, m_mpi_comm);

		std::string remote_ip = inet_ntoa(*((struct in_addr *)&recv[0]));
		LOG4CXX_INFO(logger, "GEN (" << params.node_rank << ":" << local_ip << ") is connecting (" <<  remote_ip << ") at port " << PORT);
		params.remote = new ClientSocket(remote_ip.c_str(), PORT);
		LOG4CXX_INFO(logger, "GEN (" << params.node_rank << ":" << local_ip << ") succeeded connecting");
	GEN_END
}


void YaoBase::init_environ(EnvParams &params)
{
	if (params.secu_param % 8 != 0 || params.secu_param > 128)
	{
		LOG4CXX_FATAL(logger, "security parameter k needs to be a multiple of 8 less than or equal to 128");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

	if (params.stat_param % params.node_amnt != 0)
	{
		LOG4CXX_FATAL(logger, "statistical  parameter s needs to be a multiple of cluster size");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

	// # of copies of a circuit each node is responsible for
	params.node_load = params.stat_param/params.node_amnt;

	//if (!params.circuit.load_binary(params.circuit_file))
	//{
	//	LOG4CXX_FATAL(logger, "circuit parsing failed");
	//	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	//}

	Env::init(params);

        /*
          // SS13 does not use claw-free collections... remove
          // for this version

	// synchronize claw-free collections
	ClawFree claw_free;
	claw_free.init();
	Bytes bufr(claw_free.size_in_bytes());

	if (Env::is_root())
	{
		EVL_BEGIN
			bufr = claw_free.to_bytes();
			EVL_SEND(bufr);
		EVL_END

		GEN_BEGIN
			bufr = GEN_RECV();
		GEN_END
	}

	// synchronize claw-free collections to the root evaluator's
	MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm);
	Env::claw_free_from_bytes(bufr);
        */
}


void YaoBase::init_private(EnvParams &params)
{
	static byte MASK[8] = { 0xFF, 0x01, 0x03, 0x07, 0x0F, 0x1F, 0x3F, 0x7F};

	std::ifstream private_file(params.private_file);
	std::string input;

	if (!private_file.is_open())
	{
		LOG4CXX_FATAL(logger, "file open failed: " << params.private_file);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

	EVL_BEGIN // evaluator
          //private_file >> input;          // 1st line is the evaluator's input
        //m_evl_inp.from_hex(input);
        //m_private_input.from_hex(input);
	//m_evl_inp.resize((Env::circuit().evl_inp_cnt()+7)/8);
		//m_evl_inp.back() &= MASK[Env::circuit().evl_inp_cnt()%8];
	EVL_END

	GEN_BEGIN // generator
          //private_file >> input >> input; // 2nd line is the generator's input
        //m_gen_inp.from_hex(input);
        // m_private_input.from_hex(input);
	//m_gen_inp.resize((Env::circuit().gen_inp_cnt()+7)/8);
		//m_gen_inp.back() &= MASK[Env::circuit().gen_inp_cnt()%8];
	GEN_END

	private_file.close();
}


YaoBase::~YaoBase()
{
	Env::destroy();

	int res;
	MPI_Comm_compare(m_mpi_comm, MPI_COMM_WORLD, &res);
	if (res == MPI_UNEQUAL) MPI_Comm_free(&m_mpi_comm); // MPI_COMM_WORLD can't be freed
}




inline std::string print_longlong(uint64_t l)
{
	char buf[8];
	std::string str;

	while (l >= 1000)
	{
		sprintf(buf, ",%03d", (int)(l%1000));
		str = buf + str;

		l = l / 1000LL;
	}

	sprintf(buf, "%d", (int)l);
	str = buf + str;

	return str;
}


void YaoBase::reset_timers()
{
    m_timer_gen = m_timer_evl = m_timer_mpi = m_timer_com = 0;
    m_comm_sz = 0;
}


void YaoBase::step_report(std::string step_name)
{
  double start = MPI_Wtime();
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  m_timer_mpi += MPI_Wtime() - start;
  
  step_report_no_sync(step_name);
}


void YaoBase::step_report_no_sync(std::string step_name)
{
	uint64_t all_comm_sz = 0LL;

	MPI_Reduce(&m_comm_sz, &all_comm_sz, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, m_mpi_comm);

	if (!Env::is_root())
		return;

	m_timer_mpi_vec.push_back(m_timer_mpi);
	m_timer_cmm_vec.push_back(m_timer_com);
	m_step_name_vec.push_back(step_name);
	m_comm_sz_vec.push_back(all_comm_sz);

	EVL_BEGIN
		m_timer_cmp_vec.push_back(m_timer_evl);
		LOG4CXX_INFO(logger, "EVL finish " << step_name << "");
	EVL_END

	GEN_BEGIN
		m_timer_cmp_vec.push_back(m_timer_gen);
		LOG4CXX_INFO(logger, "GEN finish " << step_name << "");
	GEN_END
}


void YaoBase::final_report()
{
	if (!Env::is_root())
          {
            std::cout << "end child process" << std::endl;
            return;
          }
	LOG4CXX_INFO(logger, "========================================================");
	LOG4CXX_INFO(logger, "Yao protocol completed");
	LOG4CXX_INFO(logger, "========================================================");

	std::string name;

	EVL_BEGIN
          name = "EVL";
        // LOG4CXX_INFO(logger, "EVL  input: " << m_evl_inp.to_hex() << "");
        
        LOG4CXX_INFO(logger, "EVL  input: " << m_private_input.to_hex() << "");
        LOG4CXX_INFO(logger, "EVL output: " << m_evl_out.to_hex() << "");
	EVL_END

	GEN_BEGIN
          name = "GEN";
        // LOG4CXX_INFO(logger, "GEN  input: " << m_gen_inp.to_hex() << "");
        LOG4CXX_INFO(logger, "GEN  input: " << m_private_input.to_hex() << "");
        LOG4CXX_INFO(logger, "GEN output: " << m_gen_out.to_hex() << "");
	GEN_END

	for (size_t i = 0; i < m_comm_sz_vec.size(); i++)
	{
		LOG4CXX_INFO
		(
			logger,
			name << " in " << m_step_name_vec[i] << "> " << std::fixed <<
			"  cmp:"  << std::setw(12) << std::setprecision(4) << m_timer_cmp_vec[i] <<
			", cmm:"  << std::setw(12) << std::setprecision(4) << m_timer_cmm_vec[i] <<
			", mpi:"  << std::setw(12) << std::setprecision(4) << m_timer_mpi_vec[i] <<
			", size:" << std::setw(16) << std::setprecision(4) << print_longlong(m_comm_sz_vec[i])
		);
	}
}


void YaoBase::get_and_size_inputs(){
  size_inputs();
  get_inputs();
}

void YaoBase::get_inputs(){
  char * raw_input_bufr;
  std::string input_string;
  
  GEN_BEGIN
    
    raw_input_bufr = (char*)malloc(sizeof(char)*m_gen_inp_cnt);
    raw_input_bufr = get_bob_input(m_gen_inp_cnt,Env::private_file());
    fprintf(stderr, "Bob's (1) input is: %s",raw_input_bufr);
    
    m_private_input.from_char_hex(raw_input_bufr, m_gen_inp_cnt);
    std::cout << "Bob's (2) input is: " << m_private_input.to_hex() << std::endl;
    // std::cout << "Bob's (3) input is: " << m_gen_inp.to_hex() << std::endl;
    
    GEN_END
      
    EVL_BEGIN
    
    raw_input_bufr = (char*)malloc(sizeof(char)*m_evl_inp_cnt);
    raw_input_bufr = get_alice_input(m_evl_inp_cnt,Env::private_file());
    fprintf(stderr,"Alice's (1) input is: %s",raw_input_bufr);
    
    m_private_input.from_char_hex(raw_input_bufr,m_evl_inp_cnt);
    std::cout << "Alice's (2) input is: " << m_private_input.to_hex() << std::endl;
    // std::cout << "Alice's (3) input is: " << m_evl_inp.to_hex() << std::endl;
    
    EVL_END

}

void YaoBase::size_inputs(){
  	static byte MASK[8] = { 0xFF, 0x01, 0x03, 0x07, 0x0F, 0x1F, 0x3F, 0x7F};

	m_gen_inp_cnt = read_alice_length(Env::private_file());
	m_evl_inp_cnt = read_bob_length(Env::private_file());

        EVL_BEGIN
          //m_evl_inp.resize((m_evl_inp_cnt+7)/8);
          //m_evl_inp.back() &= MASK[m_evl_inp_cnt%8];
          m_private_input.resize((m_evl_inp_cnt+7)/8);
          m_private_input.back() &= MASK[m_evl_inp_cnt%8];
        EVL_END

          GEN_BEGIN
	m_private_input.resize((m_gen_inp_cnt+7)/8);
	m_private_input.back() &= MASK[m_gen_inp_cnt%8];
        GEN_END
        
}

/*
void YaoBase::oblivious_transfer()
{
	reset_timers();

	double start; // time marker

	Bytes send, recv, bufr(Env::elm_size_in_bytes()*4);
	std::vector<Bytes> bufr_chunks, recv_chunks;

	G X[2], Y[2], gr, hr;
	Z s[2], t[2],  y,  a,  r;

        //std::cout << "OT Step 1" << std::endl;

	// step 1: generating the CRS: g[0], h[0], g[1], h[1]
	if (Env::is_root())
	{
          EVL_BEGIN
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
          
          start = MPI_Wtime(); // send to Gen's root process
          EVL_SEND(bufr);
          m_timer_com += MPI_Wtime() - start;
          EVL_END
            
            GEN_BEGIN
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
        
        // pre-processing
        m_ot_g[0].fast_exp();
        m_ot_g[1].fast_exp();
        m_ot_h[0].fast_exp();
        m_ot_h[1].fast_exp();
        
        // allocate memory for m_keys
        m_ot_keys.resize(Env::node_load());
        // in HBC, Env::node_load() should be 1
        // in malicious, this will vary
        for (size_t ix = 0; ix < m_ot_keys.size(); ix++)
          {
            // hbc only executes once
            m_ot_keys[ix].reserve(m_evl_inp_cnt*2);
          }
	m_timer_evl += MPI_Wtime() - start;
	m_timer_gen += MPI_Wtime() - start;
        
        // std::cout << "OT Step 2" << std::endl;
        
	// Step 2: ZKPoK of (g[0], g[1], h[0], h[1])
	// TODO

        // std::cout << "OT Step 3" << std::endl;

	// Step 3: gr=g[b]^r, hr=h[b]^r, where b is the evaluator's bit
	if (Env::is_root())
	{
          EVL_BEGIN
            start = MPI_Wtime();
          bufr.clear(); bufr.reserve(Env::exp_size_in_bytes()*m_evl_inp_cnt);
          send.clear(); send.reserve(Env::elm_size_in_bytes()*m_evl_inp_cnt*2);
          for (size_t bix = 0; bix < m_evl_inp_cnt; bix++)
            {
              r.random();
              bufr += r.to_bytes();  // to be shared with slave evaluators
              
              //byte bit_value = m_evl_inp.get_ith_bit(bix);
              byte bit_value = m_private_input.get_ith_bit(bix);
              
              send += (m_ot_g[bit_value]^r).to_bytes(); // gr
              send += (m_ot_h[bit_value]^r).to_bytes(); // hr
            }
          m_timer_evl += MPI_Wtime() - start;
          
          start = MPI_Wtime();
          EVL_SEND(send); // send (gr, hr)'s
          m_timer_com += MPI_Wtime() - start;
          
          m_comm_sz += send.size();
          EVL_END
            
            GEN_BEGIN
            start = MPI_Wtime();
          bufr = GEN_RECV(); // receive (gr, hr)'s
          m_timer_com += MPI_Wtime() - start;
          
          m_comm_sz += bufr.size();
          GEN_END
        }
        
	EVL_BEGIN // forward rs to slave evaluators
          start = MPI_Wtime();
        bufr.resize(Env::exp_size_in_bytes()*m_evl_inp_cnt);
        m_timer_evl += MPI_Wtime() - start;
        
        start = MPI_Wtime();
        MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm); // now every evaluator has r's
        m_timer_mpi += MPI_Wtime() - start;
        
        start = MPI_Wtime();
        bufr_chunks = bufr.split(Env::exp_size_in_bytes());
        // bufr_chunks now holds each value of r
        m_timer_evl += MPI_Wtime() - start;
	EVL_END

        GEN_BEGIN // forward (gr, hr)s to slave generators
        start = MPI_Wtime();
        bufr.resize(Env::elm_size_in_bytes()*m_evl_inp_cnt*2);
        m_timer_gen += MPI_Wtime() - start;
        
        start = MPI_Wtime();
        MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm); // now every generator has (gr, hr)s
        m_timer_mpi += MPI_Wtime() - start;
        
        start = MPI_Wtime();
        bufr_chunks = bufr.split(Env::elm_size_in_bytes());
        // bufr_chunks now hold all the rs for Gen
        m_timer_gen += MPI_Wtime() - start;
	GEN_END
          
          //  std::cout << "OT Step 4" << std::endl;
        
	// Step 4: the generator computes X[0], Y[0], X[1], Y[1]
	GEN_BEGIN
          for (size_t bix = 0; bix < m_evl_inp_cnt; bix++)
            {
              start = MPI_Wtime();
              gr.from_bytes(bufr_chunks[2*bix+0]);
              hr.from_bytes(bufr_chunks[2*bix+1]);
              
              if (m_ot_keys.size() > 2)
                {
                  gr.fast_exp();
                  hr.fast_exp();
                }
              m_timer_gen += MPI_Wtime() - start;
              
              for (size_t cix = 0; cix < m_ot_keys.size(); cix++)
                {
                  start = MPI_Wtime();
                  Y[0].random(); // K[0]
                  Y[1].random(); // K[1]
                  
                  m_ot_keys[cix].push_back(Y[0].to_bytes().hash(Env::k()));
                  m_ot_keys[cix].push_back(Y[1].to_bytes().hash(Env::k()));
                  
                  s[0].random(); s[1].random();
                  t[0].random(); t[1].random();
                  
                  // X[b] = ( g[b]^s[b] ) * ( h[b]^t[b] ), where b = 0, 1
                  X[0] = m_ot_g[0]^s[0]; X[0] *= m_ot_h[0]^t[0];
                  X[1] = m_ot_g[1]^s[1]; X[1] *= m_ot_h[1]^t[1];
                  
                  // Y[b] = ( gr^s[b] ) * ( hr^t[b] ) * K[b], where b = 0, 1
                  Y[0] *= gr^s[0]; Y[0] *= hr^t[0];
                  Y[1] *= gr^s[1]; Y[1] *= hr^t[1];
                  
                  send.clear();
                  send += X[0].to_bytes(); send += X[1].to_bytes();
                  send += Y[0].to_bytes(); send += Y[1].to_bytes();
                  m_timer_gen += MPI_Wtime() - start;
                  
                  start = MPI_Wtime();
                  GEN_SEND(send);
                  m_timer_com += MPI_Wtime() - start;
                  
                  m_comm_sz += send.size();
                }
            }
        
        for (size_t ix = 0; ix < m_ot_keys.size(); ix++)
          {
            assert(m_ot_keys[ix].size() == m_evl_inp_cnt*2);
          }
	GEN_END

          // std::cout << "OT Step 5" << std::endl;
          
	// Step 5: the evaluator computes K = Y[b]/X[b]^r
	EVL_BEGIN
          
          // std::cout << "OT 5 part 1" << std::endl;
          for (size_t bix = 0; bix < m_evl_inp_cnt; bix++)
            {
              start = MPI_Wtime();
              //int bit_value = m_evl_inp.get_ith_bit(bix);
              int bit_value = m_private_input.get_ith_bit(bix);
              r.from_bytes(bufr_chunks[bix]);
              m_timer_evl += MPI_Wtime() - start;
              
              for (size_t cix = 0; cix < m_ot_keys.size(); cix++)
                {
                  start = MPI_Wtime();
                  recv = EVL_RECV(); // receive X[0], X[1], Y[0], Y[1]
                  m_timer_com += MPI_Wtime() - start;
                  
                  m_comm_sz += recv.size();
                  
                  start = MPI_Wtime();
                  recv_chunks = recv.split(Env::elm_size_in_bytes());
                  
                  X[bit_value].from_bytes(recv_chunks[    bit_value]); // X[b]
                  Y[bit_value].from_bytes(recv_chunks[2 + bit_value]); // Y[b]
                  
                  // K = Y[b]/(X[b]^r)
                  Y[bit_value] /= X[bit_value]^r;
                  m_ot_keys[cix].push_back(Y[bit_value].to_bytes().hash(Env::k()));
                  m_timer_evl += MPI_Wtime() - start;
                }
            }
        
        // std::cout << "OT 5 part 2" << std::endl;
        
        for (size_t ix = 0; ix < m_ot_keys.size(); ix++)
          {
            assert(m_ot_keys[ix].size() == m_evl_inp_cnt);
          }
	EVL_END
          
        step_report("ob-transfer");
}
*/
