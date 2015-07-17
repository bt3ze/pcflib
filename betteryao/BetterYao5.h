#ifndef BETTERYAO5_H_
#define BETTERYAO5_H_

#include "YaoBase.h"
#include "Bytes.h"
#include "GarbledCircuit.h"
#include "otextension.h"

#include <vector>


// useful function for modifying Gen and Evl inputs
uint32_t ceil_log_base_2(uint32_t k); 

//Bytes make_commitment(Bytes commit_value);

class BetterYao5 : public YaoBase
{
public:
	BetterYao5(EnvParams &params);
	virtual ~BetterYao5() {}

	void start();
	
protected:
        /**
           HIGH LEVEL PROTOCOL FUNCTIONS
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
           GENERAL-PURPOSE RANDOMNESS GENERATION FUNCTIONS
           generate random seeds, keys, seed prngs, etc.
         */
        void generate_random_seeds(std::vector<Bytes> & seeds, uint32_t num_seeds);
        void seed_prngs(std::vector<Prng> & prngs, std::vector<Bytes> & seeds);
        void generate_input_keys(Prng & rng, std::vector<Bytes> & keys, uint32_t num_keys,uint32_t num_bits);



        /**
         *
         COMMITMENTS: GENERATION AND TRANSFER
         *
         */
        void gen_send_evl_commitments(std::vector<commitment_t> &commits);
        void evl_receive_gen_commitments(std::vector<Bytes> &commits, uint32_t num_commitments);
        void generate_commitments(Prng & rng, std::vector<Bytes> & keys, std::vector<commitment_t> & commitments);
        

        /**
         *
         OBLIVIOUS TRANSFER
         *
         */

        // these two high-level functions 
        // submit and select a vector for OT
        void ot_send(std::vector<Bytes> & sender_inputs);
        void ot_receive(Bytes selection_bits, std::vector<Bytes> & results_container);

        // these two functions are for sending groups of keys
        // they perform multiple parallel OTs
        void ot_send_batch(std::vector<std::vector<Bytes> > & sender_inputs);
        void ot_receive_batch(Bytes selection_bits, std::vector<std::vector<Bytes> > & results_container);


        // core OT functions
        void ot_send_init();
        void ot_receive_init();
        void ot_send_random(std::vector<Bytes> & send_inputs);
        void ot_receive_random(Bytes selection_bits, std::vector<Bytes> & results_container);


        
        /**
           STEP 1: MODIFY INPUT KEYS
           gen must generate two auxiliary inputs:
           a mask for his output and extra randomness
           that is necessary for the 2-UHF
        */

        // generates Gen's output mask
        void gen_generate_output_mask(Prng &);

        // generates 2k+lg(k) extra random bits
        void gen_generate_input_randomness(Prng &);
        
        // Eval protects her inputs by generating
        // a k-probe-resistant matrix
        void choose_k_probe_resistant_matrix();
        
        // Eval determines her protocol input based on
        // her private input and the matrix
        void evl_generate_new_input();
        
        // TODO: change this to account for eval's 
        // modified inputs
        uint32_t get_evl_inp_count(){
          return m_evl_inp_cnt; 
        }

        
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
        uint32_t get_gen_inp_size();
        // for now, the following just calles get_gen_inp_size()
        uint32_t get_gen_output_size();

        
        /**
           STEP 2: GEN INPUT COMMITMENTS
         */
        
        void generate_gen_input_keys(uint32_t circuit_num);
        void commit_to_gen_input_keys();


        /**
           STEP 3: AGREE ON OBJECTIVE CIRCUIT
         */

        // Objective Circuit Agreement

        // K Probe
        void eval_announce_k_probe_matrix();
        // choose 2UFH
        void collaboratively_choose_2UHF();
        // interactive coin flipping protocol
        Bytes flip_coins(size_t len);


        /**
           STEP 4: COMMITMENT TO INPUT (AND OUTPUT) LABELS
           note: Gen commits to a method of generating output labels,
                 not the labels themselves
         */

        // these functions take a circuit_num
        // which is an index to a processor's
        // list of circuit evaluations
        void generate_eval_input_keys(uint32_t circuit_num);
        void generate_gen_input_label_commitments(uint32_t circuit_num);

        // gen_commit_to_io_labels declared above


        void commit_to_gen_input_labels();
        void commit_to_eval_input_labels();



        /**
           STEP 5: EVAL'S INPUT OTS
         */

        // Eval uses ot methods here, declared above in auxiliary functions
        
        // after Eval receives her input keys
        // the seeds need to be hashed to the right size
        // (we use much longer seeds than necessary
        // because of our old OT)
        void hash_eval_input_keys(std::vector<Bytes> & source, std::vector<Bytes> & destination, uint32_t num_bits);
        
        
        /**
           STEP 6: CUT AND CHOOSE
         */

        // Cut and Choose
        void evl_select_cut_and_choose_circuits();
        void special_circuit_ot();
        
        // Gen must select his input decommitments
        void select_input_decommitments(std::vector<commitment_t>& source, std::vector<commitment_t>& dest, Bytes & select_bits);

        // Gen and Eval transfer information in the "circuit OT"
        void transfer_evaluation_circuit_info();
        void transfer_check_circuit_info();

        // transferred circuit information is masked by Prngs
        // so that Eval can only decrypt what she asked for
        // these functions help encrypt and decrypt information
        
        // send and receive a whole vector
        void gen_decommit_and_send_masked_vector(Prng & mask_generator, std::vector<commitment_t> & vec);
        void evl_receive_masked_vector(Prng & mask_generator, std::vector<Bytes> & destination, uint32_t chunk_size, uint32_t len);
                
        // send and receive one message
        void gen_send_masked_info(Prng & mask_generator, Bytes info, uint32_t chunk_size);
        void evl_receive_masked_info(Prng & mask_generator, Bytes & destinatuion, uint32_t chunk_size);
        
        // len is passed to tell Evl how many messages to ignore
        // if Evl cannot decrypt, she ignores
        void evl_ignore_masked_info(uint32_t len);
                

      
        /**
           STEP 7: GARBLE CIRCUITS
           AND CHECK CIRCUIT CONSISTENCY

        */
        
        // these functions help Evl regenerate her check circuits
        void evl_regenerate_circuits(uint32_t circuit_num);
        void evl_check_garbled_circuit_commitments(uint32_t circuit_num);
        void evl_check_commitment_regeneration(uint32_t circuit_num);
        bool check_received_commitments_vs_generated(std::vector<Bytes> & received, std::vector<commitment_t> & generated,uint32_t p);
        void evl_set_inp_keys(uint32_t circuit_num);

        void evl_inputs_transform(std::vector<Bytes> &source, std::vector<Bytes> &dest);


        // Gen and Evl compute their 2UHF circuits
        void generate_2UHF();
        void evaluate_2UHF();
        
        // and initialize/evaluate circuits
        void initialize_circuits();
        void evaluate_circuits();


        /**
           STEP 8: RETRIEVE OUTPUTS
         */


        void post_evaluation_checks(); // TODO

        // Eval proves Gen's output authenticity
        void gen_output_auth_proof();
                  


        /**
           PROTOCOL VARIABLES
         */

        /**
           Gen's inputs
         */

        // Gen needs a couple of extra variables for his own inputs
        Bytes m_gen_output_mask;

        // Gen must generate extra random input
        // to make the output of the 2-UHF appear random
        Bytes m_gen_aux_random_input;

        // both parties also have m_private_input
        // both parties have m_gen_inp_cnt
        // and m_evl_inp_cnt
        
        
        /**
           vectors of keys
         */

        std::vector<std::vector<Bytes> > m_gen_inp_keys;
        std::vector<std::vector<Bytes> > m_evl_inp_keys;
        std::vector<std::vector<Bytes> > m_evl_inp_ot_keys;

        std::vector<std::vector<Bytes> > m_evl_received_keys;


        // Gen uses this vector to store Eval's hashed input keys
        // since he inputs them to the OT in longer form
        // (size of group element)
        // and he uses the k-bit hash as the key itself
        std::vector<std::vector<Bytes> > m_evl_hashed_inp_keys;


        /**
           vectors of commitments
         */

        // these input commitments are used in step 3
        std::vector<std::vector<commitment_t> >   m_gen_inp_commitments;
        std::vector<std::vector<Bytes> >          m_gen_received_input_commitments;
        

        // these input label commitments are used in
        // Step 5: Gen's I/O Commitments
        // right before Eval's input OTs
        std::vector<std::vector<commitment_t> >          m_gen_inp_label_commitments;
        std::vector<std::vector<commitment_t> >          m_evl_inp_label_commitments;
        



        /**
           eval's received commitments
         */

        // when Gen sends his input label commitments in Step 4/5,
        // Evl puts them in here
        std::vector<std::vector<Bytes> >          m_gen_received_label_commitments;
        std::vector<std::vector<Bytes> >          m_evl_received_label_commitments;
        // Gen also commits to a way of generating output labels
        // but that's embedded in the logic of the code



        
        // Eval puts the decommitments that she receives during the Cut and Choose
        // in these containers
        // the first, gen_inp_commitments, corresponds to what she received in Step 2 (Step 3 in SS13)
        // the second, gen_inp_label_commitments, corresponds to step 4 (step 5 in SS13)
        std::vector<std::vector<Bytes> >                   m_cc_recv_gen_inp_commitments; 
        std::vector<std::vector<Bytes> >                   m_cc_recv_gen_inp_label_commitments; 



        // *********************************************************
        // this protocol implementation requires THREE sets of PRNGS
        // 1) random generators that create Gen's input keys
        //    and are used for Gen's I/O label commitments
        //    (Gen's input labels and Eval's input labels)
        //    these information will be checked in the cut-and-choose
        // 2) random generators to construct Gen's input commitments
        //    (in step 3)
        //    which will be decommitted as the choose in cut-and-choose
        // 3) random generators used to mask the information
        //    sent in the special OT that we use for cut-and-choose
        //    Eval receives 1-of-2 of these seeds,
        //    and uses the seeded generator to then decrypt the information
        //    she chose for each circuit, whether check information
        //   or evaluation information
        // *********************************************************

        // we maintain the random seeds for each set of PRNGs
        // because we will need to access them
        // (technically, we can do without them,
        // but I think it makes the code cleaner to store and differentiate) 
        std::vector<Bytes>                   m_circuit_seeds;
        std::vector<Prng>		     m_circuit_prngs;
        
        std::vector<Bytes>                   m_commitment_seeds;
        std::vector<Prng>                    m_commitment_prngs;

        // the otp-prng (or one-time pad prng)
        // is #3 above, and is used in our "special OT"
        // that we use to transmit information after the cut-and-choose
        // is selected
        // on Gen's side, each circuit will have 2 prngs (and 2 seeds)
        // Eval will get one of the seeds to decrypt as chooses
        std::vector<Bytes>                   m_otp_seeds;
        std::vector<Prng>                    m_otp_prngs;
        
        
        //   The next set of seeds is for generation circuits
        //   the (randomly generated) seed is used to seed the circuit
        //   object's PRNG to generate new wire keys.
        //   For lack of a better place to put it, these are generated
        //  during circuit information transfer
        
        std::vector<Bytes>                   m_key_generation_seeds;


        /**
           CUT AND CHOOSE
           these variables track which circuits are 
           evaluation circuits and which are check circuits
         */
	
        Bytes                           m_chks;
        
        // the following variable must be a class member
        // because if it is scoped in a function
        // it might go away before everyone is done
        Bytes                           m_all_chks;


        /**
           GARBLED CIRCUIT ACCESSORS
           and variables that the circuits need
         */       

        // access to the garbled circuits (important!)
        std::vector<GarbledCircuit>           m_gcs;

        // R keys
        std::vector<Bytes> m_R;
        
        // placeholder for the k-probe-resistant matrix we'll need
        std::vector<Bytes>                   m_k_probe_resistant_matrix;

        // this matrix defines the 2-UHF that is used to enforce
        // Gen's input consistency
        // each entry in the array is a row the size of gen's inputs
        // there are K rows (column length = k)
        std::vector<Bytes>                   m_2UHF_matrix;
        
        // this vector holds all of Evl's input hashes
        std::vector<Bytes>                   m_2UHF_hashes;
        
        // this vector tracks Gen's permutation bits
        std::vector<Bytes>     m_gen_inp_permutation_bits;

        // and these will track Gen's select bits
        // which are defined to be the XOR of the permutation and input bits
        // and tell the circuit which alice key represents 0
        std::vector<Bytes> m_gen_select_bits;
        
        std::vector<Bytes> m_const_0_keys;
        std::vector<Bytes> m_const_1_keys;

};


#endif /* BETTERYAO5_H_ */
