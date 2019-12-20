#!/bin/bash

FoBiS.py build -s ./src -m Makefile -compiler intel -fc "mpiifort" -lflags " -traceback -init=snan -init=arrays -check all -ftrapuv" -cflags " -c -traceback -init=snan -init=arrays -check all -ftrapuv"

#build all fortran programs from src
#please note that this is a dirty trick, we need a better solution
#maybe to gig in FoBiS.py options/code
echo "all: \$(addprefix \$(DEXE),\$(EXES))" >> Makefile
