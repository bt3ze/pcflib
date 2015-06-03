#!/bin/sh

mpirun -n 4 ./gen5 24 8 $1 $2 127.0.0.1 5000 1
