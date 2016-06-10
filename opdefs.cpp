#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <search.h>
#include <errno.h>
#include <string.h>

#include "opdefs.h"

#include <emmintrin.h>
#include <vector>


void copy_key1(void* source_key, void * dest_key){
  //  __m128i *new_key = 0; 
 
  //fprintf(stdout,"copy key \n");

  if (source_key != 0)
    {
      // first argument is size, second argument is allignment
      //new_key = (__m128i*)_mm_malloc(sizeof(__m128i), sizeof(__m128i));

      //num_copies++;
      //clock_gettime(CLOCK_REALTIME, &copy_start);
      
      _mm_storeu_si128(reinterpret_cast<__m128i*>(dest_key),*reinterpret_cast<__m128i*>(source_key));
     
      //clock_gettime(CLOCK_REALTIME, &copy_end);
      //copy_time += ( copy_end.tv_sec - copy_start.tv_sec )
      //  + ( copy_end.tv_nsec - copy_start.tv_nsec )
      //  / BILN;

    } else{
    fprintf(stderr,"no copy\n");
  }
  
  //  return new_key;
 
}


struct PCFGate _gate;
#define DEBUG_OUTPUT 0

void clear_op(struct PCFState * st, struct PCFOP  * op)
{
  struct clear_op_data * data = (struct clear_op_data*)op->data;
  uint32_t i;
  for(i = st->base; i < st->base + data->localsize; i++)
    {
      st->wires[i].value = 0;
      
      //st->copy_key(st->constant_keys[0],st->wires[i].keydata);
      copy_key1(st->constant_keys[0],st->wires[i].keydata);

      
      st->wires[i].flags = KNOWN_WIRE;
    }
}

void nop(struct PCFState * st, struct PCFOP * op)
{}

void initbase_op(struct PCFState * st, struct PCFOP * op)
{
  ENTRY ent, * r;
  ent.key = (char*) "main";
  if(hsearch_r(ent, FIND, &r, st->labels) == 0)
        {
          fprintf(stderr, "Problem searching hash table for main: %s\n", strerror(errno));
          abort();
        }
 
  uint32_t * target = (uint32_t *) r->data;
  st->PC = *target;

  st->base += *((uint32_t*)op->data);
}

void mkptr_op(struct PCFState * st, struct PCFOP * op)
{
  uint32_t idx = *((uint32_t*)op->data);
  //assert(st->wires[idx + st->base].flags == KNOWN_WIRE);

  
  //if(DEBUG_OUTPUT)
  // fprintf(stderr, "Making pointer %d -> %d\n", st->wires[idx + st->base].value, st->wires[idx+st->base].value + st->base);
  
  // .. now we clear keydata by copying 0 key into it rather than freeing
  //st->copy_key(st->constant_keys[0], st->wires[idx+st->base].keydata);
  copy_key1(st->constant_keys[0], st->wires[idx+st->base].keydata);
 
  st->wires[idx + st->base].value = st->wires[idx+st->base].value + st->base;
  
}

void const_op(struct PCFState * st, struct PCFOP * op)
{
  struct const_op_data * data = (const_op_data *) op->data;
  uint32_t idx = data->dest + st->base;
  st->wires[idx].value = data->value;

  if(data->value >= 2) {
    //st->copy_key(st->constant_keys[0], st->wires[idx].keydata);
    copy_key1(st->constant_keys[0], st->wires[idx].keydata);
  }
  else
    {
      //assert((data->value == 0) || (data->value == 1));
      
      //st->copy_key(st->constant_keys[data->value],st->wires[idx].keydata);
      copy_key1(st->constant_keys[data->value],st->wires[idx].keydata);

    }
  st->wires[idx].flags = KNOWN_WIRE;
}

void add_op(struct PCFState * st, struct PCFOP* op)
{
  struct arith_op_data * data = (arith_op_data *) op->data;
  //assert(st->wires[st->base + data->op1].flags == KNOWN_WIRE);
  //assert(st->wires[st->base + data->op2].flags == KNOWN_WIRE);
  st->wires[st->base + data->dest].flags = KNOWN_WIRE;

  st->wires[st->base + data->dest].value =
    st->wires[st->base + data->op1].value +
    st->wires[st->base + data->op2].value;

  //st->copy_key(st->constant_keys[0], st->wires[data->dest + st->base].keydata);
  copy_key1(st->constant_keys[0], st->wires[data->dest + st->base].keydata);
  

}

void mul_op(struct PCFState * st, struct PCFOP* op)
{
  struct arith_op_data * data =(arith_op_data *) op->data;
  //assert(st->wires[st->base + data->op1].flags == KNOWN_WIRE);
  //assert(st->wires[st->base + data->op2].flags == KNOWN_WIRE);
  st->wires[st->base + data->dest].flags = KNOWN_WIRE;
  st->wires[st->base + data->dest].value =
    st->wires[st->base + data->op1].value *
    st->wires[st->base + data->op2].value;

  //st->copy_key(st->constant_keys[0], st->wires[data->dest + st->base].keydata);
  copy_key1(st->constant_keys[0], st->wires[data->dest + st->base].keydata);

}

void join_op(struct PCFState * st, struct PCFOP * op)
{
  struct join_op_data * data = (join_op_data *)op->data;

  //  uint32_t d_idx = data->dest + st->base;
  int32_t i = 0;
  uint32_t cval = 0;
  
  for(i = data->nsources-1; i >= 0; i--)
    {
      cval = cval << 1;
      //assert(st->wires[data->sources[i] + st->base].flags == KNOWN_WIRE);
      //assert(st->wires[data->sources[i] + st->base].value < 2);
      cval += st->wires[data->sources[i] + st->base].value;
    }

  st->wires[data->dest + st->base].flags = KNOWN_WIRE;

  //st->copy_key(st->constant_keys[0], st->wires[data->dest + st->base].keydata);
  copy_key1(st->constant_keys[0], st->wires[data->dest + st->base].keydata);
    
  st->wires[data->dest+st->base].value = cval;
}

void bits_op(struct PCFState * st, struct PCFOP * op)
{
  //if(DEBUG_OUTPUT){
  //  fprintf(stderr,"Bits \n");
  //} 
  struct bits_op_data * data = (bits_op_data *)op->data;
  uint32_t s_idx = data->source + st->base;
  uint32_t i = 0, cval;

  //assert(st->wires[s_idx].flags == KNOWN_WIRE);

  cval = st->wires[s_idx].value;
  
  //if(DEBUG_OUTPUT){
  // fprintf(stderr,"cval: %i\treldest: %i\tabsdest: %i\n",cval, data->source, s_idx);
  //}

  for(i = 0; i < data->ndests; i++)
    {
      st->wires[data->dests[i] + st->base].value = (cval & 0x01);
      st->wires[data->dests[i] + st->base].flags = KNOWN_WIRE;

      //st->copy_key(st->constant_keys[cval & 0x01], st->wires[data->dests[i] + st->base].keydata);
      copy_key1(st->constant_keys[cval & 0x01], st->wires[data->dests[i] + st->base].keydata);

      cval = cval >> 1;

    }

  //  if(DEBUG_OUTPUT){
  // fprintf(stderr,"end Bits \n");
  //} 
}


void call_op (struct PCFState * st, struct PCFOP * op)
{
  struct call_op_data * data = (struct call_op_data*)op->data;
  
  ENTRY * ent = data->target;
  ENTRY * r = 0;

  if(strcmp(data->target->key, "alice") == 0)
    {

      uint32_t i = 0, idx = 0;

      // Get the argument to this function.
      // the inputs will be included in the 32 wires preceding the "call" op's newbase,
      // since they are passed as the only argument to alice().
      // so step backwards through them, shifting the accumulator.
      // remember that the LSB is the lowest wire
      
      if(st->inp_i == 0)
        {
          for(i = 1; i <= 32; i++)
            {

              // shift over to make room for next bit
              idx = idx << 1;

              // each wire holding the address must have a wire value (0,1)
              //assert(st->wires[st->base + data->newbase - i].value < 2);
               // an input MUST be a known wire.
              //assert(st->wires[st->base + data->newbase - i].flags == KNOWN_WIRE);
              // gather bit into index accumulator
              idx += st->wires[st->base + data->newbase - i].value;
            }
          
          // store this index to identify which input to grab
          // st->inp_idx = idx*32;
          st->inp_idx = idx;
        }
      // idx (similarly, st->idx) now has the argument that was passed to alice()

      // fprintf(stderr,"alice call: %u, (newbase: %u)\n ",idx, data->newbase);

      // this loop with inp_i ensures that we call the function 
      // to retrieve outputs 32 times (for a whole word)
      // it iterates this once for each bit in the word
      // and then resets the input index (inp_i) for the next input
      if(st->inp_i < 32)
        { 
          // fprintf(stderr,"inp_i < 32 \t inp_i=%x\n",st->inp_i);
          i = st->inp_i;
          st->inp_i++;
          st->input_g.wire1 = st->inp_idx + i;
          st->input_g.wire2 = st->inp_idx + i;
          st->input_g.reswire = st->base + data->newbase + i;
          st->input_g.truth_table = 5; //0101?
          st->input_g.tag = TAG_INPUT_A;

          // if input wire already has a value, delete it
          //st->copy_key(st->constant_keys[0], st->wires[st->input_g.reswire].keydata);
          copy_key1(st->constant_keys[0], st->wires[st->input_g.reswire].keydata);
          


          if(st->inp_idx + i < st->alice_in_size)
            { // if alice provided a value for this wire, then get it
              // fprintf(stderr,"A\t");
              
              //st->copy_key(st->callback(st, &st->input_g), st->wires[st->input_g.reswire].keydata);
              st->curgate = &st->input_g;
              
              //copy_key1(st->callback(st, &st->input_g), st->wires[st->input_g.reswire].keydata);
              static garble_cb_arg * cb_arg;
              cb_arg->st = st;
              cb_arg->gate = &st->input_g;
              st->callback(cb_arg);
            }
          else 
            { // otherwise, fill with zeros
            
              //fprintf(stdout,"filling zeros for uninitialized wire");
              //st->copy_key(st->constant_keys[0],st->wires[st->input_g.reswire].keydata);
              copy_key1(st->constant_keys[0],st->wires[st->input_g.reswire].keydata);

            }
          
          st->wires[st->input_g.reswire].flags = UNKNOWN_WIRE;
          // Not yet done with function call
          // so decrement the program counter to go through again
          st->PC--;
          // fprintf(stderr,"End Loop\n");
        }
      else
        { // now we have iterated 32 times for an input
          // so we reset inp_i for the next time
          st->inp_i = 0;
        } 
      
      // fprintf(stderr,"End Alice Call\n");
    }
  else if(strcmp(data->target->key, "bob") == 0)
    {

      uint32_t i = 0, idx = 0;
      // Get the argument to this function
      // we only need to do the check and input retrieval once
      if(st->inp_i == 0)
        {
          for(i = 1; i <= 32; i++)
            {
              idx = idx << 1;
              //assert(st->wires[st->base + data->newbase - i].value < 2);
              //assert(st->wires[st->base + data->newbase - i].flags == KNOWN_WIRE);
              idx += st->wires[st->base + data->newbase - i].value;
            }
          //st->inp_idx = idx * 32;
          st->inp_idx = idx;
        }
 
      // fprintf(stderr,"bob call: %x\n",idx);

      if(st->inp_i < 32)
        {
          // fprintf(stderr,"inp_i < 32 \t inp_i=%x\n",st->inp_i);
          i = st->inp_i;
          st->inp_i++;
          st->input_g.wire1 = st->inp_idx+i;
          st->input_g.wire2 = st->inp_idx+i;
          st->input_g.reswire = st->base + data->newbase + i;
          st->input_g.truth_table = 5;
          st->input_g.tag = TAG_INPUT_B;
          

          if(st->inp_idx + i < st->bob_in_size)
            {
              //st->copy_key(st->callback(st, &st->input_g), st->wires[st->input_g.reswire].keydata);
              //copy_key1(st->callback(st, &st->input_g), st->wires[st->input_g.reswire].keydata);
              
              // questionable if we even need this line
              // st->curgate = &st->input_g;
              
              static garble_cb_arg * cb_arg;
              cb_arg->st = st;
              cb_arg->gate = &st->input_g;
              st->callback(cb_arg);
            }
          else
            {
         
              //st->copy_key(st->constant_keys[0],st->wires[st->input_g.reswire].keydata);
              copy_key1(st->constant_keys[0],st->wires[st->input_g.reswire].keydata);
              
            }
          st->wires[st->input_g.reswire].flags = UNKNOWN_WIRE;
          // Not yet done with function call
          st->PC--;
        }
      else
        {
          st->inp_i = 0;
        }


    }
  else if(strcmp(data->target->key, "output_alice") == 0)
    {
      // fprintf(stdout,"output alice opdef: %x\n",st->inp_i);
      uint32_t i;
      if(st->inp_i < 32)
        {
          i = st->inp_i;
          st->inp_i++;
          st->input_g.wire1 = st->base + data->newbase - (32 - i);
          st->input_g.wire2 = st->base + data->newbase - (32 - i);
          st->input_g.reswire = st->base + data->newbase - (32 - i);
          st->input_g.truth_table = 5;
          st->input_g.tag = TAG_OUTPUT_A;
          st->curgate = &st->input_g;
          st->PC--;
          static garble_cb_arg * cb_arg;
          cb_arg->st = st;
          cb_arg->gate = &st->input_g;
          st->callback(cb_arg);
          //st->callback(st, &st->input_g);

        }
      else
        st->inp_i = 0;
    }
  else if(strcmp(data->target->key, "output_bob") == 0)
    {
      uint32_t i;
      if(st->inp_i < 32)
        {
          i = st->inp_i;
          st->inp_i++;
          st->input_g.wire1 = st->base + data->newbase - (32 - i);
          st->input_g.wire2 = st->base + data->newbase - (32 - i);
          st->input_g.reswire = st->base + data->newbase - (32 - i);
          st->input_g.truth_table = 5;
          st->input_g.tag = TAG_OUTPUT_B;
          st->curgate = &st->input_g;
          st->PC--;
          static garble_cb_arg * cb_arg;
          cb_arg->st = st;
          cb_arg->gate = &st->input_g;
          st->callback(cb_arg);
          //st->callback(st, &st->input_g);
          
        }
      else
        st->inp_i = 0;
    }
  else
    {
      struct activation_record * newtop = (activation_record *)  malloc(sizeof(struct activation_record));
      //check_alloc(newtop);
      if(newtop == 0){
        fprintf(stderr, "Failed to allocate memory: %s\n", strerror(errno));
        exit(-1);
      }

      newtop->rest = st->call_stack;
      newtop->ret_pc = st->PC;
      newtop->base = st->base;
      st->call_stack = newtop;

      if(hsearch_r(*ent, FIND, &r, st->labels) == 0)
        {
          fprintf(stderr, "Problem searching hash table for %s: %s\n", ent->key, strerror(errno));
          abort();
        }
 
      long * target = (long *) r->data;
      st->PC = *target;
      st->base += data->newbase;
    }
}

void branch_op(struct PCFState * st, struct PCFOP * op)
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

  //assert(st->wires[data->cnd_wire + st->base].flags == KNOWN_WIRE);
  //assert(st->wires[data->cnd_wire + st->base].value < 2);

  //if(DEBUG_OUTPUT)
  //  fprintf(stderr, "Branch: to %s on w %x\n", data->target->key, st->wires[data->cnd_wire + st->base].value);

  uint32_t old_PC = st->PC;
  if(st->wires[data->cnd_wire + st->base].value == 1)
    st->PC = *target;
  
  if(old_PC > *target){
    for(uint32_t j = *target; j < old_PC; j++){
      st->ops[j].num_exec--;
    }
  }
}

void gate_op(struct PCFState * st, struct PCFOP * op)
{
  //fprintf(stderr,"gate op\n");
  struct PCFGate * data = (struct PCFGate*)op->data;
  uint32_t op1idx = data->wire1 + st->base;
  uint32_t op2idx = data->wire2 + st->base;
  uint32_t destidx = data->reswire + st->base;
  uint8_t bits[4];
  int8_t i = 0;
  uint8_t tab = data->truth_table;
  //void * tmp = st->wires[destidx].keydata;

  //fprintf(stdout,"wire1: %u, wire2: %u\tabsolute1: %u, absolute2: %u\n",
  //        data->wire1, data->wire2, op1idx, op2idx);
  //fprintf(stdout,"wire1: %u\twire2: %u\n",
  //        st->wires[op1idx].keydata, st->wires[op2idx].keydata);


  for(i = 0; i < 4; i++)
    {
      bits[i] = tab & 0x01;
      tab = tab >> 1;
    }

  ////assert(st->curgate == 0);
  ////assert(data->truth_table < 16);
  
  if((st->wires[op1idx].flags != KNOWN_WIRE) || (st->wires[op2idx].flags != KNOWN_WIRE))
    {
      // Time for the callback
      ////assert((st->wires[op1idx].keydata != 0) && (st->wires[op2idx].keydata != 0)); // //asserts no null pointers, but asserting not known constant wires is unclear. should be taken care of by the KNOWN_WIRE checks
 
      st->curgate = &_gate;

      st->curgate->wire1 = op1idx;
      st->curgate->wire2 = op2idx;
      st->curgate->reswire = destidx;


#ifdef BETTERYAO
      // construct the truth table according to some method that should be dependent on the
      // input gates' permutation bits
      tab = 
        ((data->truth_table & 1) << 3) |
        ((data->truth_table & 2) << 1) |
        ((data->truth_table & 4) >> 1) |
        ((data->truth_table & 8) >> 3);
#else
      tab = data->truth_table;
#endif

      st->curgate->truth_table = tab;
      st->curgate->tag = TAG_INTERNAL;
      st->wires[destidx].flags = UNKNOWN_WIRE;

      
      clock_gettime(CLOCK_REALTIME, &(st->requestStart2));
      
      // note: we will be replacing this line with just the callback,
      // and make sure that the copy key is evaluated as part of the callback
      //copy_key1(st->callback(st, st->curgate),st->wires[destidx].keydata);

      //st->callback(st,st->curgate);
      static garble_cb_arg * cb_arg;
      cb_arg->st = st;
      cb_arg->gate = &st->input_g;
      st->callback(cb_arg);

      clock_gettime(CLOCK_REALTIME, &(st->requestEnd2));
      st->accum2 += ( st->requestEnd2.tv_sec - st->requestStart2.tv_sec )
        + ( st->requestEnd2.tv_nsec - st->requestStart2.tv_nsec )
        / BILLION;
      
    }
  else
    {
      // Check that we are dealing only with bits
      ////assert((st->wires[op1idx].value < 2) && (st->wires[op2idx].value < 2));
      ////assert(((st->wires[op1idx].value) + (2*(st->wires[op2idx].value))) < 4);
      
      //if((bits[(st->wires[op1idx].value) + (2*(st->wires[op2idx].value))]) >= 2){
      //  fprintf(stderr, "Problem!\n");
      // }
      
      //      //assert((bits[(st->wires[op1idx].value) + (2*(st->wires[op2idx].value))]) < 2);

      //st->copy_key(st->constant_keys[bits[(st->wires[op1idx].value) + (2*(st->wires[op2idx].value))]], st->wires[destidx].keydata );
      copy_key1(st->constant_keys[bits[(st->wires[op1idx].value) + (2*(st->wires[op2idx].value))]], st->wires[destidx].keydata );

      
      st->wires[destidx].value = bits[(st->wires[op1idx].value) + (2*(st->wires[op2idx].value))];
      st->wires[destidx].flags = KNOWN_WIRE;
    }

 

}

void copy_op(struct PCFState * st, struct PCFOP * op)
{
  struct copy_op_data * data = (struct copy_op_data*)op->data;
  uint32_t i;
  uint32_t dest = data->dest + st->base;
  uint32_t source = data->source + st->base;

  //assert(data->width > 0);
  for(i = 0; i < data->width; i++)
    {
      

      if(st->wires[source+i].keydata != 0) {
        //st->copy_key(st->wires[source+i].keydata, st->wires[dest+i].keydata );
        copy_key1(st->wires[source+i].keydata, st->wires[dest+i].keydata );
      }
      else {
        
        //st->copy_key(st->constant_keys[0], st->wires[dest+i].keydata);
        copy_key1(st->constant_keys[0], st->wires[dest+i].keydata);
      }

      st->wires[dest+i].value = st->wires[source+i].value;
      st->wires[dest+i].flags = st->wires[source+i].flags;
    }
}

void indir_copy_op(struct PCFState * st, struct PCFOP * op)
{
  struct copy_op_data * data = (struct copy_op_data*)op->data;
  uint32_t dest = st->wires[data->dest + st->base].value;
  uint32_t source = data->source + st->base;
  uint32_t i;
  //assert(data->width > 0);
  for(i = 0; i < data->width; i++)
    {

      //if(DEBUG_OUTPUT) {
      //  if(st->wires[source+i].keydata != 0) {
      //    fprintf(stderr, "Copying value: %d %x\n", dest, *((uint32_t*)st->wires[source+i].keydata));
      //  }      
      //}

      //if(st->wires[source+i].keydata != 0) {
      //st->copy_key(st->wires[source+i].keydata,st->wires[dest+i].keydata);
      copy_key1(st->wires[source+i].keydata,st->wires[dest+i].keydata);

      //}
      
      //else {
      //st->copy_key(st->constant_keys[0], st->wires[dest+i].keydata);
      //  copy_key1(st->constant_keys[0], st->wires[dest+i].keydata);

      //st->wires[dest+i].keydata = 0;
      //}
      
      st->wires[dest+i].value = st->wires[source+i].value;
      st->wires[dest+i].flags = st->wires[source+i].flags;
    }
}

void copy_indir_op(struct PCFState * st, struct PCFOP * op)
{
  struct copy_op_data * data = (struct copy_op_data*)op->data;
  uint32_t dest = data->dest + st->base;
  uint32_t source = st->wires[data->source + st->base].value;
  uint32_t i;
  //assert(data->width > 0);
  for(i = 0; i < data->width; i++)
    {
      //if(DEBUG_OUTPUT) {
      //  if(st->wires[source+i].keydata != 0) {
      //    fprintf(stderr, "Copying value: %d %x\n", dest, *((uint32_t*)st->wires[source+i].keydata));
      //  }
      // }


      //if(st->wires[source+i].keydata != 0) {
      //st->copy_key(st->wires[source+i].keydata, st->wires[dest+i].keydata);
      
      copy_key1(st->wires[source+i].keydata, st->wires[dest+i].keydata);
     
      //  } else { 
      //        
      //st->copy_key(st->constant_keys[0], st->wires[dest+i].keydata);
      //  copy_key1(st->constant_keys[0], st->wires[dest+i].keydata);
        
      // }

      st->wires[dest+i].value = st->wires[source+i].value;
      st->wires[dest+i].flags = st->wires[source+i].flags;
    }

}

void ret_op(struct PCFState * st, struct PCFOP * op)
{
  struct activation_record * rec = st->call_stack;

  if(st->call_stack == 0)
    st->done = -1;
  else
    {
      st->call_stack = rec->rest;
      st->PC = rec->ret_pc;
      st->base = rec->base;
      free(rec);
    }
}
