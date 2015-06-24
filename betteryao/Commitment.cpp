#ifndef COMMITMENT_CPP_
#define COMMITMENT_CPP_

#include "Commitment.h"
#include "Env.h"
#include "Hash.h"
#include "Prng.h"


commitment_t make_commitment(Prng &rnd, Bytes msg){
  Bytes r = rnd.rand_bits(Env::k());
  commitment_t comm;
  comm.r = r;
  comm.msg = msg;
  return comm;
}

commitment_t make_commitment(Bytes rnd, Bytes msg){
  commitment_t comm;
  comm.r = rnd;
  comm.msg = msg;
  return comm;
}

commitment_t reconstruct_commitment(Bytes decom){
  // should be some room left for the message
  assert(decom.size() > Env::k()/8);
  commitment_t comm;
  comm.r = Bytes(&decom[0],&decom[Env::k()/8]);
  comm.msg = Bytes(&decom[Env::k()/8],&decom[decom.size()]);
  return comm;

}

Bytes commit(commitment_t comm){
  return (comm.r + comm.msg).hash(Env::k());
}

Bytes decommit(commitment_t comm){
  return comm.r + comm.msg;
}

std::vector<Bytes> decommit_to_vector(std::vector<commitment_t> & vec){
  std::vector<Bytes> ret;
  for(int i = 0; i < vec.size();i++){
    ret.push_back(decommit(vec[i]));
  }
}

bool verify_commitment(commitment_t & com, Bytes & decom){
  return decom == decommit(com).hash(Env::k());
}



#endif
