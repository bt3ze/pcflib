#ifndef BETTERYAO5_H_
#define BETTERYAO5_H_

#include "YaoBase.h"
//#include "garbled_circuit_m.h"
#include "GarbledBase.h"
#include "GarbledMal.h"

#include <vector>

class BetterYao5 : public YaoBase
{
public:
	BetterYao5(EnvParams &params);
	virtual ~BetterYao5() {}

	void start();

	//void oblivious_transfer();
	void cut_and_choose();
	void cut_and_choose2();
	void consistency_check();
	void circuit_evaluate();

protected:
	void ot_init();
	void ot_random(); // sender has m pairs of l-bit strings, and receiver has m bits
	void cut_and_choose2_ot();
	void cut_and_choose2_precomputation();
	void cut_and_choose2_evl_circuit(size_t ix);
	void cut_and_choose2_chk_circuit(size_t ix);

	Bytes flip_coins(size_t len);
        void seed_m_prngs(size_t num_prngs, std::vector<Bytes> seeds);

	void proc_gen_out();
	void proc_evl_out();
        
        /**
           a couple of functions used as gadgets in the protocol
           for enforcing privacy or security properties
        */
        // Eval proves Gen's output authenticity
        void gen_output_auth_proof();
        
        // Eval protects her inputs by generating
        // a k-probe-resistant matrix
        void choose_k_probe_resistant_matrix();
        
        // Eval determines her protocol input based on
        // her private input and the matrix
        void evl_generate_new_input();
        
        /**
          gen must generate two auxiliary inputs:
          a mask for his output and extra randomness
          that is necessary for the 2-UHF
        */
        // calls the two other functions
        void gen_generate_aux_inputs();
        // generates Gen's output mask 
        // (referred to as e)
        void gen_generate_output_mask();
        // generates 2k+lg(k) extra random bits
        void gen_generate_input_randomness();
        
        /**
           Gen has a couple of inputs segments, and therefore gets
           extra input accessors
         */

        // get Gen's output mask
        Bytes get_gen_output_mask();
        // get Gen's extra inputs
        Bytes get_gen_input_randomness();
        // get all of Gen's inputs concatenated
        // private_inp || e || rho_{2k+lg(k)}
        Bytes get_gen_full_input();

	size_t                          m_ot_bit_cnt;
	Bytes                           m_ot_recv_bits;
        std::vector<Bytes>              m_ot_out;

	// variables for cut-and-choose
	Bytes                           m_chks;
	Bytes                           m_all_chks;

	// variables for Yao protocol
        std::vector<Bytes>                   m_gen_inp_masks;
        std::vector<Bytes>                   m_coms;
        std::vector<Bytes>                   m_rnds;
        //std::vector<garbled_circuit_m_t>     m_gcs; 
        std::vector<GarbledMal>             m_gcs;

        // variables for Gen's input check
        std::vector<Bytes>                   m_gen_inp_hash;
        std::vector<std::vector<Bytes> >      m_gen_inp_decom;
        std::vector<Bytes>                   m_matrix;
        
        std::vector<Prng>		     m_prngs;

        // Gen needs a couple of extra variables for his own inputs
        Bytes m_gen_output_mask;
        Bytes m_gen_aux_random_input;

};


#endif /* BETTERYAO5_H_ */
