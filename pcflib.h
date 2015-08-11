#ifndef __PCFLIB_H
#define __PCFLIB_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <pthread.h>
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <search.h>

  struct PCFState;

  /* This should be immutable -- we should be able to safely make a shallow copy. */
typedef struct PCFOP {
  void * data;
  void (* op)(struct PCFState*, struct PCFOP*);
}PCFOP;

typedef struct PCFGate {
  uint32_t wire1;
  uint32_t wire2;
  uint32_t reswire;

  /* The truth table is laid out like this:
     reswire = truth_table{wire1 + 2 * wire2}

     In other words, the truth table bits b0-b3 will be:

     wire1 | wire2 | reswire
      0    |  0    |  b0
      1    |  0    |  b1
      0    |  1    |  b2
      1    |  1    |  b3
   */
  uint8_t truth_table;

  // The tag
  uint32_t tag;
}PCFGate;


  // gates are assigned the following values as their "tags"
  // input and output retrievals are also considered gates
  enum {TAG_INTERNAL = 0, TAG_INPUT_A = 2, TAG_INPUT_B = 4, TAG_OUTPUT_A = 6, TAG_OUTPUT_B = 8};

  // wires are assigned these values as their "flags"
  // wires can either be known to all parties, or unknown
  // if unknown, the information they hold must not be leaked
enum {KNOWN_WIRE = 0, UNKNOWN_WIRE = 1};

typedef struct wire {
  uint32_t value;
  uint8_t flags;

  /* The key(s) associated with this wire, which will be set by some
     external function. */
  void * keydata;
} wire;

  void * getWireKey(struct wire *);
  void setWireKey(struct wire *, void *);


  // the activation record is used as one would expect
  // it's a linked list that contains return addresses and base pointers
  // of all the function calls performed by the circuit
struct activation_record {
  uint32_t ret_pc;
  uint32_t base;
  struct activation_record * rest;
};


  // check_alloc is used frequently to ensure that memory has been allocated properly from the heap. if not, the program quits
  void check_alloc(void * ptr);

typedef struct PCFState {
  // base pointer on the circuit's emulated stack
  uint32_t base;
  // program counter on the circuit's emulated machine
  uint32_t PC;

  // pointer to the wire table
  wire * wires;
  uint32_t wire_table_size;

  // pointer to the ops
  // come back to this one
  PCFOP * ops;

  // instruction count?
  // appears to be a counter for number of instructions
  // as a safety mechanism to ensure PC hasn't gone awry
  uint32_t icount; 

  // the labels data structure holds the addresses
  // of all of the functions and is used to look them up
  // during function calls
  struct hsearch_data * labels;

  // input size counters
  uint32_t alice_in_size;
  uint32_t bob_in_size;

  // these two are indices that deal with how inputs are collected
  int32_t inp_i; // not sure specific purpose
  uint32_t inp_idx; // not sure specific purpose


  // pointer to the currently evaluated gates
  PCFGate * curgate;

  // pointer to two constants: 0 and 1
  void * constant_keys[2];

  // PCFState contains an input gate that it uses when processing inputs
  PCFGate input_g;

  struct activation_record * call_stack;

  uint8_t done; // flag for whether the computation is done

  /* This is the state of the program that is utilizing this library,
     which may be needed by the callback function to generate the keys
     for the wires in the circuit. */
  void * external_circuit; // used to hold a reference to one of the garbled circuits (m_gcs)

  /* This function is called when a gate should be emitted.  It should
     create the appropriate keys for the output wire of the gate, and
     return them. */
  void * (*callback)(struct PCFState *, struct PCFGate*);

  /* The function that will be used to make copies of the keys
     associated with a wire. */
  void * (*copy_key)(void *);

  /* The function that will be used to delete keys when a wire is
     destroyed. */
  void (*delete_key)(void*);
}PCFState;

  enum {ALICE = 0, BOB = 1};

  void set_external_circuit(struct PCFState *, void *);
  void * get_external_circuit(struct PCFState *);
  void set_key_delete_function(struct PCFState *, void (*)(void*));
  void set_key_copy_function(struct PCFState *, void *(*)(void*));
  void set_callback(struct PCFState *, void* (*)(struct PCFState *, struct PCFGate *));
  PCFGate * get_next_gate(PCFState *);
  PCFState * load_pcf_file(const char *, void *, void *, void *(*)(void*));

  void set_constant_keys(PCFState *, void *, void*);

  uint32_t get_input_size(PCFState *, uint32_t);

  wire * getWire(struct PCFState *, uint32_t);
  void * get_wire_key(struct PCFState *, uint32_t);
  void set_wire_key(struct PCFState *, uint32_t, void *);

  uint32_t read_alice_length(const char *);
  uint32_t read_bob_length(const char *); 

  char * get_alice_input(uint32_t,const char *);
  char * get_bob_input(uint32_t,const char *);

  /*
    obsolete:
    void reinitialize(PCFState *);
    PCFState * copy_pcf_state(struct PCFState *);
    void make_internal_thread(PCFState * st);
  */

#ifdef __cplusplus
}
#endif
#endif //__PCFLIB_H
