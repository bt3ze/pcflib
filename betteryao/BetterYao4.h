#ifndef BETTERYAO4_H_
#define BETTERYAO4_H_

#include "YaoBase.h"
//#include "GarbledCct3.h"
#include "garbled_circuit_m.h"
#include <vector>

class BetterYao4 : public YaoBase
{
public:
	BetterYao4(EnvParams &params);
	virtual ~BetterYao4() {}

	virtual void start();

	void oblivious_transfer();
	void cut_and_choose();
	void cut_and_choose2();
	void consistency_check();
	void circuit_evaluate();

private:
	void ot_init();
	void ot_random(); // sender has m pairs of l-bit strings, and receiver has m bits
	void cut_and_choose2_ot();
	void cut_and_choose2_precomputation();
	void cut_and_choose2_evl_circuit(size_t ix);
	void cut_and_choose2_chk_circuit(size_t ix);

	Bytes flip_coins(size_t len);

	void proc_gen_out();
	void proc_evl_out();

	// variables for IKNP03 OT-extension implementation
	G                               m_ot_g[2];
	G                               m_ot_h[2];

	size_t                          m_ot_bit_cnt;
	Bytes                           m_ot_recv_bits;
        std::vector<Bytes>                   m_ot_send_pairs;
        std::vector<Bytes>                   m_ot_out;

        std::vector<std::vector<Bytes> >          m_ot_keys; // ot output

	// variables for cut-and-choose
	Bytes                           m_chks;
	Bytes                           m_all_chks;

	// variables for Yao protocol
        std::vector<Bytes>                   m_gen_inp_masks;
        std::vector<Bytes>                   m_coms;
        std::vector<Bytes>                   m_rnds;
        std::vector<garbled_circuit_m_t>     m_gcs; 

        // variables for Gen's input check
        std::vector<Bytes>                   m_gen_inp_hash;
        std::vector<std::vector<Bytes> >      m_gen_inp_decom;
        std::vector<Bytes>                   m_matrix;

        std::vector<Prng>		     m_prngs;
	uint32_t                        m_gen_inp_cnt;
	uint32_t                        m_evl_inp_cnt;
};


#endif /* BETTERYAO4_H_ */
