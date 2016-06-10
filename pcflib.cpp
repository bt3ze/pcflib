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

//#ifdef __cplusplus
//extern "C" {
//#endif

//#include "cthreadpool/thpool.h"

//#ifdef __cplusplus
//}
//#endif

//#define _parallel_exec 1
#define _serial_exec 1


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
      op->num_exec = 0;

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

  fprintf(stdout,"build begin\n");

  // table good for the whole circuit, in the worst case
  // must be cleared for each function though
  int32_t * owned_by = (int32_t *)malloc(sizeof(int32_t*)*NUM_WIRES);
  
  for(int32_t j = 0; j < (int32_t)NUM_WIRES; j++){
    owned_by[j]=0; // note that -1 means "unowned"
    // or should we really set the initial owner of every wire to 0 -- the first instruction?
  }

  clock_gettime(CLOCK_REALTIME, &(st->requestStart2));

  fprintf(stdout,"apply flow\n");

  // run through the list of ops
  // building dependencies
  for(i=0; i< st->icount; i++){
    apply_flow(st,&st->ops[i],owned_by);
  }

  
  clock_gettime(CLOCK_REALTIME, &(st->requestEnd2));

  
  for(i=0; i< st->icount; i++){
    //st->ops[i].num_exec= 0;
    st->ops[i].preds.clear();
    st->ops[i].succs.clear();
    apply_flow(st,&st->ops[i],owned_by);
  }

  fprintf(stdout,"\nbuild time: %f\n",
          ( st->requestEnd2.tv_sec - st->requestStart2.tv_sec )
          + ( st->requestEnd2.tv_nsec - st->requestStart2.tv_nsec )
          / BILLION);
  
  /*
  for(i=0; i < st->icount;i++){
    PCFOP * op = &st->ops[i];
    fprintf(stdout,"idx: %i\ttype: %i\n preds:",i,op->type);
    for(unsigned int j = 0; j < op->preds.size(); j++){
      fprintf(stdout,"%i\t",op->preds[j]);
    }
    fprintf(stdout,"succs:\t");
    for(unsigned int j = 0; j < op->succs.size(); j++){
      fprintf(stdout,"%i\t",op->succs[j]);
    }
    fprintf(stdout,"\n");
  }
  */
  free(owned_by); // don't need it after the graph has been built

  return st;
}

void * execute_op(void * arg){
  PCFState * st = ((exec_op_arg *)(arg))->st;
  uint32_t PC = ((exec_op_arg *)(arg))->PC;
  st->ops[PC].op(st, &st->ops[PC]);
  st->ops[PC].num_exec++;
  return st;
}


void * execute_queue(void * arg){
  // this is what the ready op manager queue should do
  queue_op_arg * ar = (queue_op_arg *)arg;
  PCFState * st = ar->st;
  uint32_t PC = ar->PC;
  PCFOP * op = ar->op;
  threadpool * qpool = ar->qpool;
  threadpool * wpool = ar->wpool;


  uint32_t num_exec = op->num_exec;
  //bool valid = true; qq
  
  for(uint32_t j = 0; j < op->succs.size();j++){
    if(st->ops[op->succs[j]].num_exec <= num_exec){
      thpool_add_work(*qpool, execute_queue,arg);
      return st;
    }
  }
  static exec_op_arg warg;
  warg.PC = PC;
  warg.st = st;
  thpool_add_work(*wpool, execute_op,(void *)&warg);
  return st;

}

void evaluate_circuit(struct PCFState *st){

  // will need here a couple of threads
  // because no internet, use a threading abstraction for now
  
  // execute circuit on main thread,
  // dispatching gates in parallel
  

  // note: we can implement the gate queue as a heap
  // and execute gates if all of the dependencies have been met
  // notice that if we arrive at an instruction in the heap
  // then it should be necessary that all of its predecessors have been executed 
  // since we start with a topological sort
  // some checking may still be necessary
  // but the idea is:
  // when finished executing an instruction, add all of its successors to the heap
  // make sure that the heap has unique entries (don't add a gate address twice)
  // ... unfortunately, a heap might not do such a great job of exploiting parallelism
  // ... still under design


  // must allocate threads
  // _and_ a garbling buffer for each thread

  // if a branch target is earlier than a branch, then must decrement the number of times each op has been executed in between them by 1


#ifndef __APPLE__
  clock_gettime(CLOCK_REALTIME, &(st->requestStart)); 
#endif


#ifdef _parallel_exec

  threadpool queuepool = thpool_init(1);
  threadpool workpool = thpool_init(1);

    // some kind of new instruction queue here
    // some kind of ready instruction queue here


  while(st->done == 0)
    {
      static exec_op_arg arg;
      arg.PC = st->PC;
      arg.st = st;
      
      //   at this point, add checks for blocking instructions
      // if (blocking) wait
      uint32_t optype = st->ops[st->PC].type;
       // if(optype == BRANCH_OP || optype == COPY_INDIR_OP || optype == INDIR_COPY_OP || optype == CALL_OP){
       //thpool_wait(queuepool);
       //thpool_wait(workpool);
       // }
      
       //thpool_add_work(queuepool, execute_op,(void *)&arg);

      if(optype == COPY_OP || optype == BITS_OP || optype == JOIN_OP){
        thpool_add_work(queuepool, execute_op,(void *)&arg);
      } 
      else {
        if(optype == GATE_OP ||  optype == COPY_INDIR_OP || optype == INDIR_COPY_OP){
          thpool_wait(queuepool);
          thpool_add_work(queuepool, execute_op,(void *)&arg);
        } else{

          if(optype == BRANCH_OP || optype == CALL_OP){
            thpool_wait(queuepool);
            thpool_wait(workpool);
          }
          
          
        //     Execute  (&arg);
          execute_op(&arg);
        }
      }
      st->PC++;
      

      //  }
      // note that first I am implementing the parallel garbler,
      // then I can make adaptations to communicate properly with the evaluator
      // but want to do first construct the parallel infrastructure to run in sequence
    }
  // thpool_wait(queuepool);
#endif

#ifdef _serial_exec
   while(st->done == 0)
    {
      static exec_op_arg arg;
      arg.PC = st->PC;
      arg.st = st;
      execute_op(&arg);
      st->PC++;
    }
#endif


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

//void set_callback(struct PCFState * st, void* (*callback)(struct PCFState *,struct PCFGate*))

void set_callback(struct PCFState * st, void* (*callback)(struct garble_cb_arg *))
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
