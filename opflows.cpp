#include "opflows.h"
#include "opdefs.h"

void bits_flow(struct PCFState * st, struct PCFOP * op)
/*
(BITS :dest (dest-list) :op1 wire1)  
       preds: next to use wire1
       succs: last to use (d for d in dest-list)
 */
{
  // update ownership
  struct bits_op_data * data = (bits_op_data *)op->data;
  uint32_t s_idx = data->source + st->base;  
  s_idx++;
}

void gate_flow(struct PCFState * st, struct PCFOP * op)
/*
(GATE :dest dst :op1 op1 :op2 op2 :truth-table tt)
       preds: last to use op1,op2
       succs: next to use dst
 */
{
  struct PCFGate * data = (struct PCFGate*)op->data;
  uint32_t op1idx = data->wire1 + st->base;
  uint32_t op2idx = data->wire2 + st->base;
  uint32_t destidx = data->reswire + st->base;
  uint32_t sum = op1idx + op2idx + destidx;
  sum++;
}

void const_flow(struct PCFState * st, struct PCFOP * op)
/*
(CONST :DEST dst :OP1 op1 )
       preds: last to use dst
       succs: next to use dst
 */
{
  struct const_op_data * data = (const_op_data *) op->data;
  uint32_t idx = data->dest + st->base;
  idx++;
}

// add, sub, mul
void arith_flow(struct PCFState * st, struct PCFOP * op)
/*
(ADD :dest dst :op1 op1 :op2 op2)
      preds: last to use op1,op2
      succs: next to use dst
 */
{
  struct arith_op_data * data = (arith_op_data *) op->data;
  uint32_t sum = data->dest + data->op1 + data->op2;
  sum++;

}

void initbase_flow(struct PCFState * st, struct PCFOP * op)
/*
(INITBASE :BASE base)
       preds: none/previous instruction
       succs: label immediately following
 */
{

}

void clear_flow(struct PCFState * st, struct PCFOP * op)
/*
(CLEAR :LOCALSIZE lsize ) // sets the first lsize 
       preds: label immediately preceding
       succs: next to use anything in lsize
 */
{
  struct clear_op_data * data = (struct clear_op_data*)op->data;
  uint32_t i;
  for(i = st->base; i < st->base + data->localsize; i++)
    {
      //      st->wires[i].value = 0;
      //st->copy_key(st->constant_keys[0],st->wires[i].keydata);
      //      copy_key1(st->constant_keys[0],st->wires[i].keydata);
      
      st->wires[i].flags = KNOWN_WIRE;
    }

}

void copy_flow(struct PCFState * st, struct PCFOP * op)
/*
(COPY :DEST dst :OP1 op1 :OP2 len )
      preds: last to use op1 through op1+len
      succs: next to use dest through dest+len
 */
{

}

void mkprt_flow(struct PCFState * st, struct PCFOP * op)
/*
(MKPTR :DEST dst )
      preds: last to use dst
      succs: next to use dst
 */
{

}

void copy_indir_flow(struct PCFState * st, struct PCFOP * op)
/*
(COPY-INDIR :DEST dst :OP1 src :OP2 len )
      preds: last to use (src, src+len)
      succs: all following instructions (non dynamic)
this is a tough one, because it's hard to know what the destination will be
      "block"
 */
{

}

void indir_copy_flow(struct PCFState * st, struct PCFOP * op)
/*
(INDIR-COPY :DEST dst :OP1 src :OP2 len )
      preds: all previous instructions (non dynamic)
      succs: all following instructions (non dynamic)
      "block"
 */
{

}

void call_flow(struct PCFState * st, struct PCFOP * op)
/*
(CALL :NEWBASE nbase :FNAME name )
     preds: all previous instructions (since the stack will be set up immediately before)
     succs: function called with name (and instruction immediately after)
     "block"
 */
{

}

void ret_flow(struct PCFState * st, struct PCFOP * op)
/*
(RET :VALUE wireptr )
     preds: all previous instructions
     succs: instruction following the call
     "block"
 */
{

}

void branch_flow(struct PCFState * st, struct PCFOP * op)
/*
(BRANCH :CND cnd-wire :TARG lbl )
     preds: last to use cnd wire
     succs: lbl pointing to, instruction following (pick)
 */
{

}

void label_flow(struct PCFState * st, struct PCFOP * op)
/*
(LABEL :STR lblname )
     preds: prev instruction, branches pointing to it 
            in case branches point to it, only need to check that lower-
            numbered predecessors have passed. otherwise must add to ready queue after every branch instruction (if chosen)
     succs: next instruction
 */
{

}

void join_flow(struct PCFState * st, struct PCFOP * op)
/*
(JOIN :DEST dst-wire :OP1 (list of input wires) )
    preds: last to use input wires
    succs: next to use dst-wire
 */
{

}
