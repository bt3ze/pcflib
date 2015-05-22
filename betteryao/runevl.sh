#!/bin/sh

mpirun -n 4 ./evl 24 4 $1 $2 127.0.0.1 5000 1
