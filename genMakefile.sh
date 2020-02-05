#!/bin/bash

FoBiS.py build -s ./src -m Makefile -compiler intel -fc "mpiifort" -lflags " -traceback -init=snan -init=arrays -check all -ftrapuv" -cflags " -c -traceback -init=snan -init=arrays -check all -ftrapuv"

FoBiS.py build -s ./src -m Makefile.opt -compiler intel -fc "mpiifort" -lflags " -traceback -O3 -ipo" -cflags " -c -traceback -O3 -ipo" --obj_dir ./obj_opt --mod_dir ./mod_opt

#build all fortran programs from src
#please note that this is a dirty trick, we need a better solution
#maybe to gig in FoBiS.py options/code
echo "all: \$(addprefix \$(DEXE),\$(EXES))" >> Makefile
