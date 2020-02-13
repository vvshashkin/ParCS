#!/bin/bash

TESTNAME="global diagnostics test"
EXE=./$1/TEST_GLOBAL_DIAG_MAIN
echo $TESTNAME $EXE
mpirun -n 12 $EXE
