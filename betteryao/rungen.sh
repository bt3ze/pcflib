#!/bin/sh

mpirun -n 2 ./splitgen 24 2 $1 $2 127.0.0.1 5000 1