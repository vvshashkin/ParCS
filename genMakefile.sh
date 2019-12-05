#!/bin/bash

FoBiS.py build -s ./src -m Makefile -compiler intel -fc "mpiifort" -lflags " -traceback -init=snan -init=arrays -check all -ftrapuv" -cflags " -c -traceback -init=snan -init=arrays -check all -ftrapuv"
