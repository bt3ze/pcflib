void Yao::oblivious_transfer()
{
	reset_timers();

	double start; // time marker

	Bytes send, recv, bufr(Env::elm_size_in_bytes()*4);
	std::vector<Bytes> bufr_chunks, recv_chunks;

	G X[2], Y[2], gr, hr;
	Z s[2], t[2],  y,  a,  r;

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
		//m_ot_keys.resize(Env::node_load());
		m_ot_keys.resize(1);
		//for (size_t ix = 0; ix < m_ot_keys.size(); ix++)
		//{
		//	m_ot_keys[ix].reserve(Env::circuit().evl_inp_cnt()*2);
		//}
		m_ot_keys[0].reserve(m_evl_inp_cnt*2);
	m_timer_evl += MPI_Wtime() - start;
	m_timer_gen += MPI_Wtime() - start;

	// Step 2: ZKPoK of (g[0], g[1], h[0], h[1])
	// TODO

	// Step 3: gr=g[b]^r, hr=h[b]^r, where b is the evaluator's bit
	if (Env::is_root())
	{
		EVL_BEGIN
			start = MPI_Wtime();

				bufr.clear();
				bufr.reserve(Env::exp_size_in_bytes()*m_evl_inp_cnt);

				send.clear();
				send.reserve(Env::elm_size_in_bytes()*m_evl_inp_cnt*2);
				for (size_t bix = 0; bix < m_evl_inp_cnt; bix++)
				{
					r.random();
					bufr += r.to_bytes();  // to be shared with slave evaluators

					byte bit_value = m_evl_inp.get_ith_bit(bix);
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
		m_timer_evl += MPI_Wtime() - start;
	EVL_END

	GEN_BEGIN // forward (gr, hr)s to slave generators
		start = MPI_Wtime();
			//bufr.resize(Env::elm_size_in_bytes()*Env::circuit().evl_inp_cnt()*2);
			bufr.resize(Env::elm_size_in_bytes()*m_evl_inp_cnt*2);
		m_timer_gen += MPI_Wtime() - start;

		start = MPI_Wtime();
			MPI_Bcast(&bufr[0], bufr.size(), MPI_BYTE, 0, m_mpi_comm); // now every Bob has bufr
		m_timer_mpi += MPI_Wtime() - start;

		start = MPI_Wtime();
			bufr_chunks = bufr.split(Env::elm_size_in_bytes());
		m_timer_gen += MPI_Wtime() - start;
	GEN_END

	// Step 4: the generator computes X[0], Y[0], X[1], Y[1]
	GEN_BEGIN
		//for (size_t bix = 0; bix < Env::circuit().evl_inp_cnt(); bix++)
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

	// Step 5: the evaluator computes K = Y[b]/X[b]^r
	EVL_BEGIN
		for (size_t bix = 0; bix < m_evl_inp_cnt; bix++)
		{
			start = MPI_Wtime();
				int bit_value = m_evl_inp.get_ith_bit(bix);
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

		for (size_t ix = 0; ix < m_ot_keys.size(); ix++)
		{
			assert(m_ot_keys[ix].size() == m_evl_inp_cnt);
		}
	EVL_END
	step_report("ob-transfer");
}
