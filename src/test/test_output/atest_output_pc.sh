#!/bin/bash

TESTNAME="output test"
EXE=./$1/TEST_OUTPUT_MAIN
echo $TESTNAME "(" $EXE")"
mpirun -n 48 $EXE
