#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = obj/
DMOD    = mod/
DEXE    = ./
LIBS    =
FC      = mpiifort
OPTSC   =  -c -traceback -init=snan -init=arrays -check all -ftrapuv -module mod
OPTSL   =  -traceback -init=snan -init=arrays -check all -ftrapuv -module mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)MAIN: $(MKDIRS) $(DOBJ)main.o
	@rm -f $(filter-out $(DOBJ)main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) MAIN
$(DEXE)FIVEPOINTFILTER_MAIN: $(MKDIRS) $(DOBJ)fivepointfilter_main.o
	@rm -f $(filter-out $(DOBJ)fivepointfilter_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) FIVEPOINTFILTER_MAIN

#compiling rules
$(DOBJ)buffer_mod.o: src/buffer_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grid_function_mod.o: src/grid_function_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)topology_mod.o: src/topology_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)tile_mod.o: src/tile_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_profile_mod.o: src/exchange_profile_mod.f90 \
	$(DOBJ)partition_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)partition_mod.o: src/partition_mod.f90 \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_mod.o: src/exchange_mod.f90 \
	$(DOBJ)exchange_profile_mod.o \
	$(DOBJ)buffer_mod.o \
	$(DOBJ)grid_function_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parcs_mpi_mod.o: src/ParCS_mpi_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_factory_mod.o: src/exchange_factory_mod.f90 \
	$(DOBJ)exchange_mod.o \
	$(DOBJ)exchange_profile_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)topology_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_mod.o: src/test/test_mod.f90 \
	$(DOBJ)tile_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)grid_function_mod.o \
	$(DOBJ)exchange_mod.o \
	$(DOBJ)exchange_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)main.o: src/test/main.f90 \
	$(DOBJ)test_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)fivepointfilter_mod.o: src/test/FivePointFilter/FivePointFilter_mod.f90 \
	$(DOBJ)grid_function_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)fivepointfilter_main.o: src/test/FivePointFilter/FivePointFilter_main.f90 \
	$(DOBJ)exchange_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)grid_function_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)fivepointfilter_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
