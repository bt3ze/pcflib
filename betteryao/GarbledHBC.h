#ifndef GARBLED_HBC_H_
#define GARBLED_HBC_H_

#include "GarbledBase.h"

class GarbledHBC: public GarbledBase
{
  
 public:
  
  GarbledHBC();
  virtual ~GarbledHBC() {}

    virtual void initialize_circuit();

  virtual void gen_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &gen_inp_mask, const Bytes &seed);
  virtual void evl_init_circuit(const std::vector<Bytes> &ot_keys, const Bytes &masked_gen_inp, const Bytes &evl_inp);


};

void * gen_next_gate(struct PCFState *st, struct PCFGate *current_gate); 
void * evl_next_gate(struct PCFState *st, struct PCFGate *current_gate);


#endif
