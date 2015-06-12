#ifndef BETTERYAO5_H_
#define BETTERYAO5_H_

#include "YaoBase.h"
//#include "garbled_circuit_m.h"
#include "GarbledBase.h"
#include "GarbledMal.h"

#include <vector>


// useful function for modifying Gen and Evl inputs
uint32_t ceil_log_base_2(uint32_t k); 

//Bytes make_commitment(Bytes commit_value);

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
        
        /**
           STEP 1: MODIFY INPUT KEYS
           gen must generate two auxiliary inputs:
          a mask for his output and extra randomness
          that is necessary for the 2-UHF
        */

        // generates Gen's output mask (referred to as e)
        void gen_generate_output_mask(Prng &);

        // generates 2k+lg(k) extra random bits
        void gen_generate_input_randomness(Prng &);
        
        // Eval protects her inputs by generating
        // a k-probe-resistant matrix
        void choose_k_probe_resistant_matrix();
        
        // Eval determines her protocol input based on
        // her private input and the matrix
        void evl_generate_new_input();
        
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

        uint32_t get_gen_full_input_size();

        /**
           STEP 2: GEN INPUT COMMITMENTS
         */
        
        void generate_random_seeds();
        void seed_prngs(std::vector<Prng> & prngs, std::vector<Bytes> & seeds);

        // Gen 
        void generate_gen_input_keys();

        void generate_commitments(Prng & rng, std::vector<Bytes> & keys, std::vector<commitment_t> & commitments);
        
        void commit_to_gen_input_keys();

        void generate_input_keys(Prng & rng, std::vector<Bytes> & keys, uint32_t num_keys);

        void gen_send_evl_commitments(std::vector<commitment_t> &commits);
        void evl_receive_gen_commitments(std::vector<Bytes> &commits, uint32_t num_commitments);
        /**
           STEP 3: AGREE ON OBJECTIVE CIRCUIT
         */

        // Objective Circuit Agreement
        void eval_announce_k_probe_matrix();
        void collaboratively_choose_2UHF();


        /**
           STEP 4: COMMITMENT TO INPUT (AND OUTPUT) LABELS
           note: Gen commits to a method of generating output labels,
                 not the labels themselves
         */

        // gen_commit_to_io_labels declared above
        void generate_eval_input_keys();
        

        void generate_gen_input_commitments();
        void generate_eval_input_commitments();

        void commit_to_gen_input_labels();
        void commit_to_eval_input_labels();

        /**
           STEP 5: EVAL'S INPUT OTS
         */

        // unimplemented
        //eval_input_ot delared above
        
        /**
           STEP 6: CUT AND CHOOSE
         */

        // Cut and Choose
        void evl_select_cut_and_choose_circuits();

        // the following are legacy functions
        // will need new ones that perform the special OT
        // where Eval gets seeds for check circuits
        // and inputs for evaluation circuits
        void ot_init();
        void ot_random(); // sender has m pairs of l-bit strings, and receiver has m bits
        Bytes flip_coins(size_t len);
        void seed_m_prngs(size_t num_prngs, std::vector<Bytes> seeds);

        void cut_and_choose2_ot();
	void cut_and_choose2_precomputation();
	void cut_and_choose2_evl_circuit(size_t ix);
	void cut_and_choose2_chk_circuit(size_t ix);


        /**
           STEP 7: GARBLE CIRCUITS
           AND CHECK CIRCUIT CONSISTENCY

        */
        // unimplemented

        /**
           STEP 8: RETRIEVE OUTPUTS
         */


        // Eval proves Gen's output authenticity
        void gen_output_auth_proof();
        
        // these are legacy functions for this 
	void proc_gen_out();
	void proc_evl_out();
        
        

        /**
 
          new variables

 
         */

        /**
           Gen's inputs
         */

        // Gen needs a couple of extra variables for his own inputs
        Bytes m_gen_output_mask;

        // Gen must generate extra random input
        // to make the output of the 2-UHF appear random
        Bytes m_gen_aux_random_input;

        // Gen also has m_private_input
        // both parties have m_gen_inp_cnt
        // and m_evl_inp_cnt
        
        
        /**
           new variables
         */

        std::vector<std::vector<Bytes> > m_gen_inp_keys;
        std::vector<std::vector<Bytes> > m_evl_inp_keys;

        // list of Gen's input commitments
        std::vector<std::vector<commitment_t> >   m_gen_inp_commitments;
        std::vector<std::vector<Bytes> >          m_gen_committed_inputs;
        
        std::vector<std::vector<commitment_t> >   m_evl_inp_commitments;
        std::vector<std::vector<Bytes> >          m_evl_committed_inputs;
        
        

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
        
        // this vector tracks Gen's permutation bits
        // 
        std::vector<Bytes>     m_gen_inp_permutation_bits;

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


        // i think this was used for 
        // permutation bits
        // it is replaced with m_gen_inp_permutation_bits
        std::vector<Bytes>              m_gen_inp_masks;
        

};


#endif /* BETTERYAO5_H_ */
