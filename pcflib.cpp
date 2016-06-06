#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include <search.h>
#include <assert.h>

#include "pcflib.h"
#include "opdefs.h"
#include "opflows.h"
#include "opgen.h"

void check_alloc(void * ptr)
{
 if(ptr == 0)
    {
      fprintf(stderr, "Failed to allocate memory: %s\n", strerror(errno));
      exit(-1);
    }
}



/**
   this function returns a PCFState object
   it accepts a filename, two keys, and a function used to copy keys
   when this funciton is called, the two keys 
   (for some reason yet to be discerned) are offset by 1
*/ 

PCFState * load_pcf_file(const char * fname, void * key0, void * key1, void (*copy_key)(void*,void*))
{
  FILE * input;
  PCFState * ret; 
  char line[LINE_MAX];
  uint32_t icount = 0;
  uint32_t i = 0;

  uint32_t num_wires = NUM_WIRES;

  ret = (PCFState*)malloc(sizeof(struct PCFState));
  check_alloc(ret);

  //ret->alice_outputs = 0;
  //ret->bob_outputs = 0;
  ret->inp_i = 0;
  ret->inp_idx = 0; // should not be strictly necessary because will be set before first use

  ret->accum=0.0;
  ret->accum2=0.0;

  fprintf(stdout,"allocate const keys:\n");
  
  ret->constant_keys[0] =(void*) malloc(4*sizeof(uint32_t));
  ret->constant_keys[1] =(void*) malloc(4*sizeof(uint32_t));
  fprintf(stdout,"set up const keys:\n");

  copy_key(key0, ret->constant_keys[0]);
  fprintf(stdout,"set const key 0\n");
  copy_key(key1, ret->constant_keys[1]);
  fprintf(stdout,"set const key 1\n");
  ret->copy_key = copy_key;
  fprintf(stdout,"set copy key function\n");
  ret->call_stack = 0;
  ret->done = 0;
#ifdef __APPLE__
// TODO:fix this ????
#else
  ret->labels = (struct hsearch_data *)malloc(sizeof(struct hsearch_data));
  check_alloc(ret->labels);
  memset(ret->labels, 0, sizeof(struct hsearch_data));
#endif

  ret->wires = (struct wire *)malloc(num_wires * sizeof(struct wire)); // !note here the limit on size of the wire table
  check_alloc(ret->wires);
  
  fprintf(stdout,"allocate keys\n");
  for(i = 0; i < num_wires; i++)
    {
      //fprintf(stdout,"%i",i);
      ret->wires[i].flags = KNOWN_WIRE;
      ret->wires[i].value = 0;
      ret->wires[i].keydata = (void *)malloc(4*sizeof(uint32_t));
      check_alloc(ret->wires[i].keydata);
      
      //      copy_key(key0,ret->wires[i].keydata);
      //ret->wires[i].keydata = copy_key(key0);
    }
  fprintf(stdout,"keys allocated\n");
  

  
  ret->done = 0;
  ret->base = 1;
  ret->PC = 0;

  fprintf(stderr, "%s\n", fname);
  input = fopen(fname, "r");
  if(input == 0)
    {
      fprintf(stderr, "%s: %s\n", fname, strerror(errno));
      assert(0);
    }

  while(!feof(input))
    {
      fgets(line, LINE_MAX-1, input);
      icount++;
    }

#ifdef __APPLE__
  if(hcreate(icount) == 0)
#else
  if(hcreate_r(icount, ret->labels) == 0)
#endif
    {
      fprintf(stderr, "Unable to allocate hash table: %s\n", strerror(errno));
      abort();
      //      exit(-1);
    }

  ret->icount = icount;
  ret->ops = (PCFOP*)malloc(icount * sizeof(PCFOP));
  check_alloc(ret->ops);

  assert(fseek(input, 0, SEEK_SET) == 0);

  icount = 0;


  // read all of the ops into memory
  while(!feof(input))
    {
      PCFOP * op;
      fgets(line, LINE_MAX-1, input);
      op = read_instr(ret, line, icount);
      op->idx = icount;

      ret->ops[icount] = *op;
      
      //free(op);
      // is this a good idea after setting the pointer to it?
      // read_instr uses malloc, and then ret sets the shallow copy
      // but I really think free is a mistake
      
      icount++;
    }

  fclose(input);

  ret->wires[0].value = 1;

  ret->copy_key(ret->constant_keys[1],ret->wires[0].keydata);
  //  ret->wires[0].keydata = ret->copy_key(ret->constant_keys[1]);
  ret->wires[0].flags = KNOWN_WIRE;

  
  // at this point, we can construct a graph of the circuit using dependency analysis
  build_tree(ret);

  // and now done initialization

  return ret;
}

void finalize(PCFState * st)
{

  //  fprintf(stderr, "finalize\n");
  
  // uint32_t i = 0;
  
  // still need some way to delete the keys and free memory
  //for(i = 0; i < 200000; i++)
  //  {
  //    if(st->wires[i].keydata != 0)
  //      st->delete_key(st->wires[i].keydata);
  //  } 
  uint32_t i = 0;
  for(i=0;i< NUM_WIRES;i++){
    free(st->wires[i].keydata);
  }
  //  free(st->wires);
  free(st);
  
  //  fprintf(stderr,"done finalize\n");
  
  fprintf(stderr,"accumulator: %f\n",st->accum);
  fprintf(stderr,"accumulator2: %f\n",st->accum2);

}

void apply_flow(struct PCFState *st, struct PCFOP * op, int32_t * table){
  switch (op->type){
  case GATE_OP:
    gate_flow(st,op,table);
    break;
  case BITS_OP:
    bits_flow(st,op,table);
    break;
  case CONST_OP:
    const_flow(st,op,table);
    break;
  case ADD_OP:
  case SUB_OP:
  case MUL_OP:
    arith_flow(st,op,table);
    break;
  case INITBASE_OP:
    initbase_flow(st,op,table);
    break;
  case CLEAR_OP:
    clear_flow(st,op,table);
    break; 
  case COPY_OP:
    copy_flow(st,op,table);
    break; 
  case MKPTR_OP:
    mkprt_flow(st,op,table);
    break;
  case COPY_INDIR_OP:
    copy_indir_flow(st,op,table);
    break;
  case INDIR_COPY_OP:
    indir_copy_flow(st,op,table);
    break;
  case CALL_OP:
    call_flow(st,op,table);
    break;
  case RET_OP:
    ret_flow(st,op,table);
    break;
  case BRANCH_OP:
    branch_flow(st,op,table);
    break;
  case LABEL_OP:
    label_flow(st,op,table);
    break;
  case JOIN_OP:
    join_flow(st,op,table);
    break;
  default:
    fprintf(stdout,"error determining op type! %x\n",op->type);
    break;
  }
  
}


PCFState * build_tree(struct PCFState *st){
  uint32_t i;

  // table good for the whole circuit, in the worst case
  // must be cleared for each function though
  int32_t * owned_by = (int32_t *)malloc(sizeof(int32_t*)*st->icount);
  for(int32_t j = 0; j < (int32_t)st->icount; j++){
    owned_by[j]=-1; // note that -1 means "unowned"
  }

  clock_gettime(CLOCK_REALTIME, &(st->requestStart2));

  // run through the list of ops
  // building dependencies
  for(i=0; i< st->icount; i++){
    apply_flow(st,&st->ops[i],owned_by);
  }

  
  clock_gettime(CLOCK_REALTIME, &(st->requestEnd2));

  
  fprintf(stdout,"\nbuild time: %f\n", ( st->requestEnd2.tv_sec - st->requestStart2.tv_sec )
          + ( st->requestEnd2.tv_nsec - st->requestStart2.tv_nsec )
          / BILLION);

  free(owned_by); // don't need it after the graph has been built

  return st;
}


void evaluate_circuit(struct PCFState *st){

  // will need here a couple of threads
  // because no internet, use a threading abstraction for now
  
  // execute circuit on main thread,
  // dispatching gates in parallel
  

#ifndef __APPLE__
  clock_gettime(CLOCK_REALTIME, &(st->requestStart)); 
#endif

  uint32_t i = 0;
  i++;
  while(st->done != 0)
    {
      if(st->curgate ==0){
        //dispatch
      } else{
        // take the top off of the ready queue
        // execute it
        // add successors to ready queue, if ready
        // continue

        // in fact, I think this can be evaluated in parallel as well
        // note that whenever we get an indir_copy the data ownership analysis fails,
        // so we must think about some sequentiality
        // or figure out how to bottleneck until that operation is done 
        // it can only execute after all those before it
        // and others can only execute after it is done

        /*
         PCFOP op = readyQueue.pop();
         // this pop should be atomic and blocking
         // so only one thread at a time can get it
         
         if(op != 0){
           //
           

         } else{
           // nothing on the ready queue
           // either an error or we just need to wait
           // figure that out later
         }
        */

      }
      // note that first I am implementing the parallel garbler,
      // then I can make adaptations to communicate properly with the evaluator
      // but want to do first construct the parallel infrastructure to run in sequence
    }



#ifndef __APPLE__
  clock_gettime(CLOCK_REALTIME, &(st->requestEnd));
  st->accum += ( st->requestEnd.tv_sec - st->requestStart.tv_sec )
    + ( st->requestEnd.tv_nsec - st->requestStart.tv_nsec )
    / BILLION;
#endif

}

// this function will change for parallel garbler
// and really the whole execution paradigm
struct PCFGate * get_next_gate(struct PCFState * st)
{
  
  //fprintf(stderr, "get next gate");
  //  std::cout << "get next gate" << std::endl;
#ifndef __APPLE__
  clock_gettime(CLOCK_REALTIME, &(st->requestStart)); 
#endif

  st->curgate = 0;
  //fprintf(stderr,"program counter: %u\n",st->PC);
  while((st->curgate == 0) && (st->done == 0))
  {
    // if curgate is 0, why are we executing things?
    // that check is an indicator (a null pointer)
    // that we are on a non-gate instruction


      // call the instruction's op function
      st->ops[st->PC].op(st, &st->ops[st->PC]);
      // then increment the program counter
      st->PC++;
      // and verify that the PC has not eclipsed the instruction count
      //assert((st->PC < st->icount));
    }

  //  if((st->curgate == 0) || (st->done != 0))
  // this seems redundant, since st->curgate should only be 0 after the above loop if st->done is nonzero
  // might want to look this up in an old version of PCF - June 3 2016
 
#ifndef __APPLE__
  clock_gettime(CLOCK_REALTIME, &(st->requestEnd));
  st->accum += ( st->requestEnd.tv_sec - st->requestStart.tv_sec )
    + ( st->requestEnd.tv_nsec - st->requestStart.tv_nsec )
    / BILLION;
#endif

  if(st->done != 0)
  {
      // no more gates and the PCFState thinks it's done
    //      fprintf(stderr,"enter finalize\n");
    finalize(st);
    return 0;
  }
  else {
    // return the next gate
    return st->curgate;
  }
}


void * get_wire_key(struct PCFState * st, uint32_t idx)
{
  return st->wires[idx].keydata;
}
/*
void set_wire_key(struct PCFState * st, uint32_t idx, void * kd)
{
  setWireKey(&st->wires[idx], kd);
}
*/

/*
wire * getWire(struct PCFState * st, uint32_t idx)
{
  return &st->wires[idx];
}
*/

void set_external_circuit(struct PCFState * st, void * ec)
{
  st->external_circuit = ec;
}

void * get_external_circuit(struct PCFState * st)
{
  return st->external_circuit;
}

void set_callback(struct PCFState * st, void* (*callback)(struct PCFState *,struct PCFGate*))
{
  st->callback = callback;
}

//void set_key_copy_function(struct PCFState * st, void *(*f)(void*))
void set_key_copy_function(struct PCFState * st, void (*f)(void*,void*))
{
  st->copy_key = f;
}

void set_key_delete_function(struct PCFState * st, void (*f)(void*))
{
  st->delete_key = f;
}

uint32_t read_alice_length(const char * fname)
{
  uint32_t alice_in = -1;
  char line[LINE_MAX];
  FILE * cnf;
  cnf = fopen(fname, "r");
  // Assume that the inputs are all on one line -- probably a bad assumption, but whatever
  fgets(line, LINE_MAX-1, cnf);
  fclose(cnf);
  alice_in = (strlen(line) - 1) * 4;
  fprintf(stderr,"alice length: %u\n",alice_in);
  return alice_in;
}

uint32_t read_bob_length(const char * fname)
{
  uint32_t bob_in = 32;
  FILE * cnf;
  char line[LINE_MAX];
  cnf = fopen(fname, "r");
  // Bob's input is on the second line
  //
  // The input file format should be like this:
  //
  // 0xALICEINPUTSINHEX
  // 0x0000000000000000
  //
  // or
  //
  // 0x0000000000000000
  // 0xBOBINPUTSINHEX00
  fgets(line, LINE_MAX-1, cnf);
  fgets(line, LINE_MAX-1, cnf);
  fclose(cnf);
  bob_in = (strlen(line) - 1) * 4;
  fprintf(stderr,"bob length: %u\n",bob_in);
  return bob_in;
}

char* get_alice_input(uint32_t length, const char* fname){
  FILE * cnf;
  char * alice_in = (char*) malloc(sizeof(char)*length);
  cnf = fopen(fname,"r");
  // alice's input is on the first line.
  // read length characters into alice_in
  fgets(alice_in, length, cnf);
  fclose(cnf);
  return alice_in;
}

char* get_bob_input(uint32_t length, const char* fname){
  FILE * cnf;
  char * bob_in = (char*)malloc(sizeof(char)*length);
  char discard[LINE_MAX];
  cnf = fopen(fname,"r");
  // bob's input is on the second line
  // read length characters into bob_in
  fgets(discard, LINE_MAX-1,cnf);
  fgets(bob_in, length, cnf);
  fclose(cnf);
  return bob_in;
}
