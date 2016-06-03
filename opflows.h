#ifndef __OPFLOWS_H
#define __OPFLOWS_H


#include "pcflib.h"


void bits_flow(struct PCFState * st, struct PCFOP * op);

void gate_flow(struct PCFState * st, struct PCFOP * op);

void const_flow(struct PCFState * st, struct PCFOP * op);

// add, sub, mul
void arith_flow(struct PCFState * st, struct PCFOP * op);

void initbase_flow(struct PCFState * st, struct PCFOP * op);

void clear_flow(struct PCFState * st, struct PCFOP * op);

void copy_flow(struct PCFState * st, struct PCFOP * op);

void mkprt_flow(struct PCFState * st, struct PCFOP * op);

void copy_indir_flow(struct PCFState * st, struct PCFOP * op);

void indir_copy_flow(struct PCFState * st, struct PCFOP * op);

void call_flow(struct PCFState * st, struct PCFOP * op);

void ret_flow(struct PCFState * st, struct PCFOP * op);

void branch_flow(struct PCFState * st, struct PCFOP * op);

void label_flow(struct PCFState * st, struct PCFOP * op);

void join_flow(struct PCFState * st, struct PCFOP * op);


#endif


