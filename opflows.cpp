#include "opflows.h"
#include "opdefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <vector>
#include <assert.h>

/**
   in these definitions, we don't need to care much for the state of the stack
   since we are only interested in the owner of a wire within each function call
   to figure out which instructions are dependent on others
   however, this means we must use a separate table for each function call


   whenever an instruction touches any data, whether source or destination,
   it must mark that data as owned
   so that data are not used before ready
   or pulled out from beneath us before we use them
   ^ the above may not be strictly necessarly. i think we may be able to get away
   with only claiming successors, since future operations really should not ever 
   have to pull out data before it's used
   can come back to this later though
   
   
 */

//#define declare_swap_vars(idx,pred) uint32_t a,b;

// the following two macros deal with integers as indices
#define set_owner(new_owner,idx) owned_by[idx]=new_owner;
#define get_owner(old_owner,idx) old_owner = owned_by[idx];


// the following two macros deal with ops 
/*
#define add_succ(pred_op, succ_op)            \
  pred_op.succs.push_back(succ_op->idx);
*/
void add_succ(PCFOP * pred_op, PCFOP * succ_op){
  pred_op->succs.push_back(succ_op->idx);
}

void add_pred (PCFOP * succ_op, PCFOP* pred_op){
  succ_op->preds.push_back(pred_op->idx);
}

/*
#define add_pred(succ_op,pred_op)          \   
  succ_op.preds.push_back(pred_op.idx);
*/

#define DEBUG_OUTPUT 1

/*
#define exchange_ownership_and_pointers(pred,dest,idx)        \
  declare_swap_vars(pred,idx)                                 \
  get_prev_owner(pred,idx)                                    \
  set_owner(dest,idx)
*/

void bits_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
/*
(BITS :dest (dest-list) :op1 wire1)  
       succs: next to use wire1
       preds: last to use (d for d in dest-list)
n successors
1 predecessor
 */
{
  // update ownership
  struct bits_op_data * data = (bits_op_data *)op->data;
  uint32_t ndests = data->ndests;

  if(DEBUG_OUTPUT){
    fprintf(stdout,"BITS: %i\n",op->idx);    
  }

  for(uint32_t i = 0; i < ndests;i++){

    uint32_t old_owner;
    get_owner(old_owner, data->dests[i]);
    
    if(old_owner != (uint32_t)-1){
     
      //assert(old_owner != 0);

      if(DEBUG_OUTPUT){
        fprintf(stdout,"BITS succs \t dest: %i\t owner: %i\n", data->dests[i], old_owner);    
        fprintf(stdout,"BITS preds \t pred idx: %i\t pred type %i \t succ id: %i\t succ type %i\n",st->ops[old_owner].idx, st->ops[old_owner].type, op->idx, op->type);        
      }
      
      //add_succ(&st->ops[old_owner],op);
      //add_pred(op,&st->ops[old_owner]);

      op->preds.push_back(old_owner);
      st->ops[old_owner].succs.push_back(op->idx);
      

    }
     // get the ith dest, set its ownership
  }

  for(uint32_t i =0; i < ndests; i++){
    set_owner(op->idx,data->dests[i]);
  }
  
  // claim the source
  //owned_by[data->source] = op->idx;
  //set_owner(op->idx,data->source);


  if(DEBUG_OUTPUT){
  for(unsigned int j = 0; j < op->preds.size();j++){
    fprintf(stdout,"%i\t",op->preds[j]);
  }
  fprintf(stdout,"\n");
  }
}


void join_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
/*
(JOIN :DEST dst-wire :OP1 (list of input wires) )
    preds: last to use input wires
    succs: next to use dst-wire
 */
{
  struct join_op_data * data = (join_op_data *)op->data;
  uint32_t nsources = data->nsources;
  for(uint32_t i = 0; i < nsources;i++){

    if(DEBUG_OUTPUT){
      fprintf(stdout,"JOIN: %i\n",op->idx);    
    }

    uint32_t old_owner;
    //    fprintf(stdout,"old owner: %i \n",old_owner);    
    get_owner(old_owner, data->sources[i]);

    if(old_owner != (uint32_t)-1){

      if(DEBUG_OUTPUT){
        fprintf(stdout,"JOIN\t old_owner: %i \t old_type %i\t op_idx: %i\t op_type %i \n",st->ops[old_owner].idx,st->ops[old_owner].type, op->idx, op->type);
      }

      add_succ(&st->ops[old_owner],op);      
      add_pred(op,&st->ops[old_owner]);
    }


    // get the ith source
    // claim it
    set_owner(op->idx,data->sources[i]);
  }
  // claim the destination
  //  owned_by[data->dest] = op->idx;
  set_owner(op->idx,data->dest);
  
}


void gate_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
/*
(GATE :dest dst :op1 op1 :op2 op2 :truth-table tt)
       preds: last to use op1,op2
       succs: next to use dst
 */
{
  struct PCFGate * data = (struct PCFGate*)op->data;
  //owned_by[data->reswire] = op->idx;

  uint32_t pred1, pred2;
  get_owner(pred1, data->wire1);
  get_owner(pred2, data->wire2);

  if(DEBUG_OUTPUT){
    fprintf(stdout,"GATE\t op_idx: %i \t pred1: %i\t wire1: %i\t pred2: %i\t wire2 %i \t outwire: %i\n", op->idx, pred1, data->wire1, pred2, data->wire2, data->reswire); 
  }

  //assert(pred1 > -1);
  //assert(pred2 > -1);
  //add_pred(op,&st->ops[pred1]);
  //add_pred(op,&st->ops[pred2]);
  //add_succ(&st->ops[pred1],op);
  //add_succ(&st->ops[pred2],op);
  op->preds.push_back(pred1);
  op->preds.push_back(pred2);
  st->ops[pred1].succs.push_back(op->idx);
  st->ops[pred2].succs.push_back(op->idx);


  set_owner(op->idx,data->reswire);
  // the gate only really owns the result wire, not the inputs. they might be reused.
  //  set_owner(op->idx,data->wire1);
  //set_owner(op->idx,data->wire2);

  //uint32_t op1idx = data->wire1;
  //uint32_t op2idx = data->wire2;

}

void const_flow(struct PCFState * st, struct PCFOP * op,int32_t * owned_by)
/*
(CONST :DEST dst :OP1 op1 )
       preds: last to use dst
       succs: next to use dst
 */
{
  struct const_op_data * data = (const_op_data *) op->data;
  
  uint32_t old_owner;
  get_owner(old_owner,data->dest);
  
  if(DEBUG_OUTPUT){
    fprintf(stdout,"CONST: idx %i \t dst: %i \t prev: %i \n",op->idx, data->dest, old_owner);
  }
  
  if(old_owner != (uint32_t) -1){
    add_succ(&st->ops[old_owner],op);      
    add_pred(op,&st->ops[old_owner]);
  }

  set_owner(op->idx,data->dest);

}

// add, sub, mul
void arith_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
/*
(ADD :dest dst :op1 op1 :op2 op2)
      preds: last to use op1,op2
      succs: next to use dst
 */
{
  struct arith_op_data * data = (arith_op_data *) op->data;


  uint32_t pred1, pred2;
  get_owner(pred1, data->op1);
  get_owner(pred2, data->op2);
  op->preds.push_back(pred1);
  op->preds.push_back(pred2);
  st->ops[pred1].succs.push_back(op->idx);
  st->ops[pred2].succs.push_back(op->idx);

  set_owner(op->idx, data->dest);
  //  set_owner(op->idx, data->op1);
  //set_owner(op->idx, data->op2);
}

void initbase_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
/*
(INITBASE :BASE base)
       preds: none/previous instruction
       succs: label immediately following
 */
{
  // we will not concern ourselves all that much with ownership of the base memory
  // it will have to be initialized in order and through copy_indir and indir_copy
  // so we block anyway on them.
  // if our data dependency analysis becomes more sophisticated, this must change
  
}

void clear_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
/*
(CLEAR :LOCALSIZE lsize ) // sets the first lsize 
       preds: label immediately preceding
       succs: next to use anything in lsize
 */
{
  struct clear_op_data * data = (struct clear_op_data*)op->data;
  uint32_t i;
  //  uint32_t pred;
  for(i = 0; i < data->localsize; i++)
    {
      // this only shows up at the beginning of a function, so
      // there should be no predecessors at the beginning of a function
      // successors will set themselves later
 
      // set new ownership
      //uint32_t old_owner;
      //get_owner(old_owner, i);
      set_owner(op->idx,i);
    }

}

void copy_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
/*
(COPY :DEST dst :OP1 op1 :OP2 len )
      preds: last to use op1 through op1+len
      succs: next to use dest through dest+len
 */
{
  struct copy_op_data * data = (struct copy_op_data*)op->data;
  uint32_t dest = data->dest;
  uint32_t source = data->source;
  uint32_t i;
  //  uint32_t pred;

  // we need to mark both the source and the destination
  // as owned by this instruction to make sure
  // that nothing changes from beneath us

  if(DEBUG_OUTPUT){
    fprintf(stdout,"COPY: idx: %i \t src: %i \t dst: %i \t width: %i\n",op->idx, source,dest,data->width); 
  }

  for(i = 0; i < data->width; i++)
    {

      uint32_t old_owner;
      get_owner(old_owner, source+i);
      
      if(old_owner != (uint32_t) -1){

        if(DEBUG_OUTPUT){
          fprintf(stdout,"source: %i\t owner: %i\t \n",source+i, old_owner);
        }

        add_succ(&st->ops[old_owner],op);      
        add_pred(op,&st->ops[old_owner]);
  
      }
      //set_owner(op->idx,source+i);
    }
  for(i = 0; i < data->width; i++)
    {
      uint32_t old_owner;
      get_owner(old_owner, dest+i);
      if(DEBUG_OUTPUT){
        fprintf(stdout,"dest: %i\t old owner: %i\t new owner: %i\n",dest+i, old_owner, op->idx);
      }      
      if(old_owner != (uint32_t) -1){  
        
        add_succ(&st->ops[old_owner],op);      
        add_pred(op,&st->ops[old_owner]);
      }
      set_owner(op->idx,dest+i);
    }
}

void mkprt_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
/*
(MKPTR :DEST dst )
      preds: last to use dst
      succs: next to use dst
 */
{
  uint32_t dest = *((uint32_t*)op->data);

  uint32_t old_owner;
  get_owner(old_owner,dest);
  
  if(DEBUG_OUTPUT){
    fprintf(stdout,"MKPTR: idx %i \t dst: %i \t prev: %i \n",op->idx, dest, old_owner);
  }
  
  if(old_owner != (uint32_t) -1){
    add_succ(&st->ops[old_owner],op);      
    add_pred(op,&st->ops[old_owner]);
  }

  set_owner(op->idx,dest);


}

void copy_indir_flow(struct PCFState * st, struct PCFOP * op,int32_t * owned_by)
/*
(COPY-INDIR :DEST dst :OP1 src :OP2 len )
      preds: last to use (src, src+len)
      succs: all following instructions (non dynamic)
this is a tough one, because it's hard to know what the destination will be
      "block"
 */
{
  struct copy_op_data * data = (struct copy_op_data*)op->data;
  uint32_t i;
  uint32_t dest = data->dest;

  if(DEBUG_OUTPUT){
    fprintf(stdout,"COPY-INDIR \t dest: %i\t source: %i\t width %i \n", data->dest, data->source, data->width);
  }

  for(i = 0; i < data->width; i++)
    {
      
      uint32_t old_owner;
      get_owner(old_owner, dest+i);
      if(old_owner != (uint32_t) -1){      

        if(DEBUG_OUTPUT){
          fprintf(stdout,"dest: %i\t old owner: %i\t new owner: %i\n",dest+i, old_owner, op->idx);
        }
        add_succ(&st->ops[old_owner],op);
        add_pred(op,&st->ops[old_owner]);
      }
      set_owner(op->idx,dest+i);
      // I now own the destination addresses
      
    }

  // very hard to figure out how to claim the source addresses
}

void indir_copy_flow(struct PCFState * st, struct PCFOP * op,int32_t * owned_by)
/*
(INDIR-COPY :DEST dst :OP1 src :OP2 len )
      preds: all previous instructions (non dynamic)
      succs: all following instructions (non dynamic)
      "block"
 */
{
  struct copy_op_data * data = (struct copy_op_data*)op->data;
  uint32_t source = data->source;
  uint32_t i;
  for(i = 0; i < data->width; i++)
    {
      uint32_t old_owner;
      get_owner(old_owner, source+i);
      
      if(old_owner != (uint32_t) -1){
        add_succ(&st->ops[old_owner],op);
        add_pred(op,&st->ops[old_owner]);
      }
      //set_owner(op->idx,source+i);
    }

  // but very hard to figure out how to claim the destinations
}

void call_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
/*
(CALL :NEWBASE nbase :FNAME name )
     preds: all previous instructions (since the stack will be set up immediately before)
     succs: function called with name (and instruction immediately after, which will always depend on data owned by the call instruction if there is a return)
     note that we may have issues with global data if it's manipulated within function calls. copy-indir and indir-copy really control those accesses
     "block"
 */
{
  struct call_op_data * data = (struct call_op_data*)op->data;
  //ENTRY * ent = data->target;
  //ENTRY * r = 0;

  // update ownership
  if(strcmp(data->target->key, "alice") == 0)
    {
      for(uint32_t i = 0; i < 32; i++){
        //owned_by[data->newbase + i] = op->idx;
        set_owner(op->idx,data->newbase+i);
      }
    }
  else if(strcmp(data->target->key, "bob") == 0)
    {
      for(uint32_t i = 0; i < 32; i++){
        //owned_by[data->newbase + i] = op->idx;
        set_owner(op->idx,data->newbase+i);
      }
    }
  else if(strcmp(data->target->key, "output_alice") == 0)
    {
      for(uint32_t i = 0; i < 32; i++){
        //owned_by[data->newbase -i-1] = op->idx;
        set_owner(op->idx,data->newbase-i-1);
      }
    }
  else if(strcmp(data->target->key, "output_bob") == 0)
    {
      for(uint32_t i = 0; i < 32; i++){
        //owned_by[data->newbase -i-1] = op->idx;
        set_owner(op->idx,data->newbase-i-1);
      }
    }
  else
    {
      for(uint32_t i = 0; i < 32; i++){
        set_owner(op->idx,data->newbase+i);
        //owned_by[data->newbase + i] = op->idx;
      }
    }
}

void ret_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
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

void branch_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
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
  
  if(DEBUG_OUTPUT){
    fprintf(stdout,"BRANCH:\t target: %i\t idx: %i\n",*target,op->idx);
  }

  // update branching and label targets
  add_succ(&st->ops[*target],op);
  add_pred(op,&st->ops[*target]);
  
  add_succ(&st->ops[op->idx+1],op);
  add_pred(op, &st->ops[op->idx+1]);

  // and update the condition wire
  uint32_t old_owner;
  get_owner(old_owner,data->cnd_wire);

  add_succ(&st->ops[old_owner],op);
  add_pred(op,&st->ops[old_owner]);

  //owned_by[data->cnd_wire] = op->idx;
  set_owner(op->idx, data->cnd_wire);

  if(DEBUG_OUTPUT){
    for(unsigned int j = 0; j < op->preds.size();j++){
      fprintf(stdout,"%i \t",op->preds[j]);
    }
    fprintf(stdout,"\n");
  }
}

void label_flow(struct PCFState * st, struct PCFOP * op, int32_t * owned_by)
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

  // do we need to set it as a dependent on the previous instruction?
  // unsure
}

