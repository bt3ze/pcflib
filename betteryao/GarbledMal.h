#ifndef GARBLED_MAL_H
#define GARBLED_MAL_H

#include "GarbledBase.h"
#include <vector>

class GarbledMal: public GarbledBase
{
 public:
  GarbledMal();
  virtual ~GarbledMal() {}

  virtual void * gen_next_gate(struct PCFState *st, struct PCFGate *current_gate); 
  virtual void * evl_next_gate(struct PCFState *st, struct PCFGate *current_gate);

  std::vector<Bytes> get_gen_decommitments (){
    return this->m_gen_inp_decom;
  }

  Bytes get_gen_inp_hash(){
    return m_gen_inp_hash;
  }

   inline bool pass_check();

  virtual void initialize_circuit();
  virtual void gen_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed);
  virtual void evl_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp);


  void initialize_gen_circuit(const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed);
  void initialize_eval_circuit(const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp);

  void evl_next_gen_inp_com(const Bytes &row, size_t kx);
  void gen_next_gen_inp_com(const Bytes &row, size_t kx);


 protected:

  // specific to malicious circuit
  std::vector<Bytes>  m_gen_inp_com; // commitments?
  std::vector<Bytes>  m_gen_inp_decom; // decommitments?
  Bytes               m_gen_inp_hash;

 

};

#endif
