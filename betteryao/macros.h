#ifndef MACROS_H
#define MACROS_H

#include "mpi.h"
#include "Bytes.h"
#include <vector>

#include <stdint.h>
#include <limits.h>


Bytes recv_data(int src_node);
void send_data(int dst_node, const Bytes &data);


#ifdef GEN_CODE

	// User mode (code for the generator)
	#define GEN_BEGIN
	#define GEN_END
	#define EVL_BEGIN     if (0) {
	#define EVL_END       }
	#define GEN_SEND(d)   Env::remote()->write_bytes(d)
	#define EVL_RECV()    Env::remote()->read_bytes()
	#define EVL_SEND(d)   Env::remote()->write_bytes(d)
	#define GEN_RECV()    Env::remote()->read_bytes()

#elif defined EVL_CODE

	// User mode (code for the evaluator)
	#define GEN_BEGIN     if (0) {
	#define GEN_END       }
	#define EVL_BEGIN
	#define EVL_END
	#define GEN_SEND(d)   Env::remote()->write_bytes(d)
	#define EVL_RECV()    Env::remote()->read_bytes()
	#define EVL_SEND(d)   Env::remote()->write_bytes(d)
	#define GEN_RECV()    Env::remote()->read_bytes()

#else

/*
	// Simulation mode
	#define GEN_BEGIN     if (!Env::is_evl()) {
	#define GEN_END       }
	#define GEN_SEND(d)   send_data(Env::world_rank()+1, (d))
	#define GEN_RECV()    recv_data(Env::world_rank()+1)
	#define EVL_BEGIN     if ( Env::is_evl()) {
	#define EVL_END       }
	#define EVL_SEND(d)   send_data(Env::world_rank()-1, (d))
	#define EVL_RECV()    recv_data(Env::world_rank()-1)
*/
/*
	#define GEN_SEND(d)   Env::remote()->write_bytes(d)
	#define EVL_RECV()    Env::remote()->read_bytes()
	#define EVL_SEND(d)   Env::remote()->write_bytes(d)
	#define GEN_RECV()    Env::remote()->read_bytes()
*/

#endif

#endif
