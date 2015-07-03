#ifndef MACROS_CPP
#define MACROS_CPP

#include "macros.h"

Bytes recv_data(int src_node)
{
	MPI_Status status;

	uint32_t comm_sz;
	MPI_Recv(&comm_sz, 1, MPI_INT, src_node, 0, MPI_COMM_WORLD, &status);

	Bytes recv(comm_sz);
	MPI_Recv(&recv[0], recv.size(), MPI_BYTE, src_node, 0, MPI_COMM_WORLD, &status);

	return recv;
}


void send_data(int dst_node, const Bytes &data)
{
  assert(data.size() < INT_MAX);
  //assert(data.size() < MY_INT32_MAX);
  
  uint32_t comm_sz = data.size();
  MPI_Send(&comm_sz, 1, MPI_INT, dst_node, 0, MPI_COMM_WORLD);
  
  MPI_Send(const_cast<byte*>(&data[0]), data.size(), MPI_BYTE, dst_node, 0, MPI_COMM_WORLD);
}

#endif
