#!/bin/sh

mpirun -n 1 ./gen5hbc 80 1 $1 $2 127.0.0.1 5000 0
