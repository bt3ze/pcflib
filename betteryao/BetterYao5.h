#ifndef BETTERYAO5_H_
#define BETTERYAO5_H_

#include "YaoBase.h"
#include "garbled_circuit_m.h"

#include <vector>


class BetterYao5 : public YaoBase
{
public:
	BetterYao5(EnvParams &params);
	virtual ~BetterYao5() {}

	virtual void start();

	//void oblivious_transfer();
	virtual void consistency_check() = 0;
	virtual void circuit_evaluate() = 0;

protected:
	void cut_and_choose();
	void cut_and_choose2();

        
	virtual void ot_init() = 0;
	virtual void ot_random() = 0; // sender has m pairs of l-bit strings, and receiver has m bits
	virtual void cut_and_choose2_ot() = 0;
	virtual void cut_and_choose2_precomputation() = 0;
	virtual void cut_and_choose2_evl_circuit(size_t ix) = 0;
	virtual void cut_and_choose2_chk_circuit(size_t ix) = 0;

	virtual Bytes flip_coins(size_t len) = 0;

	virtual void proc_gen_out() = 0;
	virtual void proc_evl_out() = 0;

        // useful utility functions
        // seed the m_prngs with a vector of seeds (usuallu provided by m_ot_out, but that should be independent of the function
        void seed_m_prngs(size_t num_prngs, std::vector<Bytes> seeds);

	size_t                          m_ot_bit_cnt;
	Bytes                           m_ot_recv_bits;
        std::vector<Bytes>                   m_ot_out;
        
	// variables for cut-and-choose
        Bytes                           m_chks;
        Bytes                           m_all_chks;
        // useful because these are vectors of 8-bit integers
        // they fulfill the function of booleans

	// variables for Yao protocol
        std::vector<Bytes>                   m_gen_inp_masks;
        std::vector<Bytes>                   m_rnds;
        std::vector<garbled_circuit_m_t>     m_gcs; 

        // variables for Gen's input check
        std::vector<Bytes>                   m_gen_inp_hash;
        std::vector<std::vector<Bytes> >      m_gen_inp_decom;
        std::vector<Bytes>                   m_matrix;

        std::vector<Prng>		     m_prngs;

};


#endif /* BETTERYAO5_H_ */
