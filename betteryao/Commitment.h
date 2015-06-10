#ifndef COMMITMENT_H_
#define COMMITMENT_H_

#include "Bytes.h"
#include "Hash.h"
#include "Prng.h"


struct commitment_t {
  Bytes r;
  Bytes msg;
};

Bytes commit(commitment_t com);

Bytes decommit(commitment_t com);

commitment_t make_commitment(Prng rnd, Bytes msg);

commitment_t make_commitment(Bytes r,Bytes msg);

commitment_t reconstruct_commitment(Bytes decom);

uint32_t verify_commitment(commitment_t com, commitment_t decom);


#endif
