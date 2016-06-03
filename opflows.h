#ifndef __OPFLOWS_H
#define __OPFLOWS_H


#include "pcflib.h"


void bits_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void gate_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void const_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

// add, sub, mul
void arith_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void initbase_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void clear_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void copy_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void mkprt_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void copy_indir_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void indir_copy_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void call_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void ret_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void branch_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void label_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);

void join_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by);


#endif


