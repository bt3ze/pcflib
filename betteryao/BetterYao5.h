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

        // old protocol:
	void start();
        //void oblivious_transfer();
        void cut_and_choose();
	void cut_and_choose2();
	void consistency_check();
	void circuit_evaluate();

        // new protocol:
        void SS13();
	
protected:
        /**
           new protocol functions
         */

        // highest level functions
        // roughly correspond to steps in the SS13 protocol
        void modify_inputs();
        void gen_generate_and_commit_to_inputs();
        void agree_on_objective_circuit();
        void gen_commit_to_io_labels();
        void eval_input_OT();
        void SS13_cut_and_choose();
        void garble_and_check_circuits();
        void retrieve_outputs();

        // intermediate level functions
        void gen_generate_input_keys();
        void gen_commit_to_inputs();

        void evl_select_cut_and_choose_circuits();

        /**
           old and new protocol functions
         */

	void ot_init();
	void ot_random(); // sender has m pairs of l-bit strings, and receiver has m bits
        Bytes flip_coins(size_t len);
        void seed_m_prngs(size_t num_prngs, std::vector<Bytes> seeds);

        /**
           old protocol functions
         */
	void cut_and_choose2_ot();
	void cut_and_choose2_precomputation();
	void cut_and_choose2_evl_circuit(size_t ix);
	void cut_and_choose2_chk_circuit(size_t ix);

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

        /**
           new variables
         */
        // Gen needs a couple of extra variables for his own inputs
        Bytes m_gen_output_mask;

        // Gen must generate extra random input
        // to make the output of the 2-UHF appear random
        Bytes m_gen_aux_random_input;

        // list of Gen's input commitments
        std::vector<Bytes>   m_gen_commitments;

        // placeholder for the k-probe-resistant matrix we'll need
        std::vector<Bytes>                   m_k_probe_resistant_matrix;
        
        /**
           new and old variables
        */
        // access to the garbled circuits (important!)
        std::vector<GarbledMal>             m_gcs;
        

        // variables for cut-and-choose
	Bytes                           m_chks;
	Bytes                           m_all_chks;

        // Gen's input decommitments
        std::vector<std::vector<Bytes> >      m_gen_inp_decom;

        // m_rnd_seeds holds the seeds that Gen uses for his circuit Prngs
        // these are the "randomness" that Gen uses to construct his circuits
        // and are very important to the protocol
        std::vector<Bytes>                   m_rnd_seeds;
        
        // m_prngs are used to "extend the OTs",
        // using the random seeds sent in the OTs to compute
        // masks (or one-time pads) over the information
        // that Gen must obliviously send to Eval;
        // Eval must only have access to 1/2 of the information
        // (depending on check or evaluation circuit)
        std::vector<Prng>		     m_prngs;

        // variables for Gen's input check
        // this is the vector of Gen's input hashes, which must all be consistent
        std::vector<Bytes>                   m_gen_inp_hash;
        
        // this matrix defines the 2-UHF that is used to enforce
        // Gen's input consistency
        // each entry in the array is a row (or column?)
        std::vector<Bytes>                   m_2UHF_matrix;
        

        /**
           old variables
         */

        // m_ot_bit_cnt is the number of circuit OTs that each processor has to do
        // (set to node_load)
	size_t                          m_ot_bit_cnt;

        // m_ot_out contains the output of the OTs for Gen and Evl
        // and serves as the container for the seeds
        // for the prngs (m_prngs) that will mask the information
        // that Gen sends to Eval
        // so that Eval can only decrypt if she chose the seed in the OT
        // contemplating changing the name to m_seeds
        std::vector<Bytes>              m_ot_out;


        // i do not understand what this is used for
        // or why it is used as it is
        std::vector<Bytes>              m_gen_inp_masks;
        

};

uint32_t ceil_log_base_2(uint32_t k);


#endif /* BETTERYAO5_H_ */
