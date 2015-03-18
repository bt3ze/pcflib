#ifndef GARBLED_MAL_H
#define GARBLED_MAL_H

#include "GarbledBase.h"

class GarbledMal: public GarbledBase
{
 public:
  GarbledMal();
  ~GarbledMal();

  virtual void * gen_next_gate(struct PCFState *st, struct PCFGate *current_gate); 
  virtual void * evl_next_gate(struct PCFState *st, struct PCFGate *current_gate);

 protected:
  Bytes m_gen_inp_hash;

};

#endif
