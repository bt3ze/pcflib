#ifndef GARBLED_MAL_H
#define GARBLED_MAL_H

#include "GarbledBase.h"
#include <vector>

class GarbledMal: public GarbledBase
{
 public:
  
  GarbledMal();
  virtual ~GarbledMal() {}
  

  /**
     accessor functions
   */
  std::vector<Bytes> get_gen_decommitments (){
    return m_gen_inp_decom;
  }

  std::vector<Bytes> get_gen_commitments(){
    return m_gen_inp_com;
  }

  Bytes get_gen_inp_hash(){
    return m_gen_inp_hash;
  }

  
  /**
     pass_check is inlined for efficiency.
     it ensures that Gen's input decommitments are consistent with his commitments.
   */
  inline bool pass_check(){
    assert(m_gen_inp_decom.size() == m_gen_inp_com.size());
    
    bool pass_chk = true;
    for (size_t ix = 0; ix < m_gen_inp_decom.size(); ix++)
      {
        pass_chk &= (m_gen_inp_decom[ix].hash(Env::k()) == m_gen_inp_com[ix]);
      }
    return pass_chk;
  }

  /**
     these initialization functions set up the circuit to be either 
     a generator circuit or evaluator circuit
     setting OY keys, seeding the PRNG, and saving Gen's masked input
     the generator circuits will also, using the PRNG, select constant keys
     and choose the value R for free-XOR
  */
  // I would like to eliminate the first two of the following (since they just call the latter two anyway)
  virtual void gen_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed);

  virtual void evl_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp);

  void initialize_gen_circuit(const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed);

  void initialize_eval_circuit(const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp);

  
  /**
     the following functions generate and evaluate Gen's input commitments
     I'm not sold on the current format, since they accept a row and an index to the row
     which are supposed to correspond to some hashing or input obscuring matrix
     but that should not really be necessary
   */
  void evl_next_gen_inp_com(const Bytes &row, size_t kx);
  void gen_next_gen_inp_com(const Bytes &row, size_t kx);

  // protected:

  // the following are specific to malicious circuit
  /**
     we need vectors for Gen's input commitments and decommitments
     for both Gen and Evl -- Gen to store them and Eval to check them
  */
  std::vector<Bytes>  m_gen_inp_com; // commitments?
  std::vector<Bytes>  m_gen_inp_decom; // decommitments?
  Bytes               m_gen_inp_hash;
 
 
};


void * gen_next_malicious_gate(struct PCFState *st, struct PCFGate *current_gate); 
void * evl_next_malicious_gate(struct PCFState *st, struct PCFGate *current_gate);


#endif
