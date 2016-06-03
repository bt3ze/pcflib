#include "opflows.h"
#include "opdefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>


/**
   in these definitions, we don't need to care much for the state of the stack
   since we are only interested in the owner of a wire within each function call
   to figure out which instructions are dependent on others
   however, this means we must use a separate table for each function call
 */
void bits_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by)
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

void gate_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
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

void const_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
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
void arith_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
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

void initbase_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
/*
(INITBASE :BASE base)
       preds: none/previous instruction
       succs: label immediately following
 */
{

}

void clear_flow(struct PCFState * st, struct PCFOP * op, uint32_t * owned_by)
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

void copy_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
/*
(COPY :DEST dst :OP1 op1 :OP2 len )
      preds: last to use op1 through op1+len
      succs: next to use dest through dest+len
 */
{
    struct copy_op_data * data = (struct copy_op_data*)op->data;
    uint32_t dest = data->dest + st->base;
    uint32_t source = data->source + st->base;
    uint32_t i;
    for(i = 0; i < data->width; i++)
      {
        source += dest;
        // do something
      }
}

void mkprt_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
/*
(MKPTR :DEST dst )
      preds: last to use dst
      succs: next to use dst
 */
{
  uint32_t idx = *((uint32_t*)op->data);
  idx++;

}

void copy_indir_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
/*
(COPY-INDIR :DEST dst :OP1 src :OP2 len )
      preds: last to use (src, src+len)
      succs: all following instructions (non dynamic)
this is a tough one, because it's hard to know what the destination will be
      "block"
 */
{
  struct copy_op_data * data = (struct copy_op_data*)op->data;
  uint32_t dest = data->dest + st->base;
  uint32_t source = st->wires[data->source + st->base].value;
  uint32_t i;
  //assert(data->width > 0);
  for(i = 0; i < data->width; i++)
    {
      // do something (if possible!)
      source += dest;
    }
}

void indir_copy_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
/*
(INDIR-COPY :DEST dst :OP1 src :OP2 len )
      preds: all previous instructions (non dynamic)
      succs: all following instructions (non dynamic)
      "block"
 */
{
    struct copy_op_data * data = (struct copy_op_data*)op->data;
  uint32_t dest = st->wires[data->dest + st->base].value;
  uint32_t source = data->source + st->base;
  uint32_t i;
  //assert(data->width > 0);
  for(i = 0; i < data->width; i++)
    {
      // do something!
      source += dest;
    }
}

void call_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
/*
(CALL :NEWBASE nbase :FNAME name )
     preds: all previous instructions (since the stack will be set up immediately before)
     succs: function called with name (and instruction immediately after)
     "block"
 */
{
  struct call_op_data * data = (struct call_op_data*)op->data;
  //ENTRY * ent = data->target;
  //ENTRY * r = 0;

  if(strcmp(data->target->key, "alice") == 0)
    {
    }
  else if(strcmp(data->target->key, "bob") == 0)
    {
    }
  else if(strcmp(data->target->key, "output_alice") == 0)
    {
    }
  else if(strcmp(data->target->key, "output_bob") == 0)
    {
    }
  else
    {
      //struct activation_record * newtop = (activation_record *)  malloc(sizeof(struct activation_record));
      
    }
}

void ret_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
/*
(RET :VALUE wireptr )
     preds: all previous instructions
     succs: instruction following the call
     "block"
 */
{
  //struct activation_record * rec = st->call_stack;
  
  //if(st->call_stack == 0)
  //  st->done = -1;
  //else
  //  {
  //    st->call_stack = rec->rest;
  //    st->PC = rec->ret_pc;
  //    st->base = rec->base;
  //    free(rec);
  //  }
}

void branch_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
/*
(BRANCH :CND cnd-wire :TARG lbl )
     preds: last to use cnd wire
     succs: lbl pointing to, instruction following (pick)
 */
{
  struct branch_op_data * data = (struct branch_op_data *)op->data;
  ENTRY * ent, * r;
  ent = data->target;
  if(hsearch_r(*ent, FIND, &r, st->labels) == 0)
    {
      fprintf(stderr, "Problem searching hash table for %s: %s\n", ent->key, strerror(errno));
      abort();
    }
 
  uint32_t * target = ( uint32_t *)r->data;
  
  (*target)++;
  //if(st->wires[data->cnd_wire + st->base].value == 1)
  //  st->PC = *target;
}

void label_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
/*
(LABEL :STR lblname )
     preds: prev instruction, branches pointing to it 
            in case branches point to it, only need to check that lower-
            numbered predecessors have passed. otherwise must add to ready queue after every branch instruction (if chosen)
     succs: next instruction
 */
{
  // point to previous instruction
  // branches to be handled elsewhere

}

void join_flow(struct PCFState * st, struct PCFOP * op,uint32_t * owned_by)
/*
(JOIN :DEST dst-wire :OP1 (list of input wires) )
    preds: last to use input wires
    succs: next to use dst-wire
 */
{
  struct join_op_data * data = (join_op_data *)op->data;
    int32_t i = 0;
  uint32_t cval = 0;
  
  for(i = data->nsources-1; i >= 0; i--)
    {
         cval = cval << 1;
    }
}
