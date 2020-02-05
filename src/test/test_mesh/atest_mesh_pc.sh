#!/bin/bash

TESTNAME="mesh test"
EXE=./$1/TEST_MESH_MAIN
echo $TESTNAME $EXE
mpirun -n 6 $EXE
