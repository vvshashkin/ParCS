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
$(DEXE)FIVEPOINTFILTER_MAIN: $(MKDIRS) $(DOBJ)fivepointfilter_main.o
	@rm -f $(filter-out $(DOBJ)fivepointfilter_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) FIVEPOINTFILTER_MAIN
$(DEXE)TEST_HALO_MAIN: $(MKDIRS) $(DOBJ)test_halo_main.o
	@rm -f $(filter-out $(DOBJ)test_halo_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_HALO_MAIN
$(DEXE)TEST_MESH_MAIN: $(MKDIRS) $(DOBJ)test_mesh_main.o
	@rm -f $(filter-out $(DOBJ)test_mesh_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_MESH_MAIN
$(DEXE)TEST_EXCH_MAIN: $(MKDIRS) $(DOBJ)test_exch_main.o
	@rm -f $(filter-out $(DOBJ)test_exch_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_EXCH_MAIN
$(DEXE)TEST_METRIC_MAIN: $(MKDIRS) $(DOBJ)test_metric_main.o
	@rm -f $(filter-out $(DOBJ)test_metric_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_METRIC_MAIN

#compiling rules
$(DOBJ)tile_mod.o: src/tile_mod.f90
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

$(DOBJ)buffer_mod.o: src/buffer_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)halo_mod.o: src/halo_mod.f90 \
	$(DOBJ)grid_function_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_profile_mod.o: src/exchange_profile_mod.f90 \
	$(DOBJ)partition_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mesh_mod.o: src/mesh_mod.f90 \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grid_function_mod.o: src/grid_function_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parcs_mpi_mod.o: src/ParCS_mpi_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)const_mod.o: src/const_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)topology_mod.o: src/topology_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mesh_factory_mod.o: src/mesh_factory_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)topology_mod.o
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

$(DOBJ)ecs_geometry_mod.o: src/equiang_cs/ecs_geometry_mod.f90 \
	$(DOBJ)const_mod.o \
	$(DOBJ)topology_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_factory_mod.o: src/equiang_cs/ecs_halo_factory_mod.f90 \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_mod.o: src/equiang_cs/ecs_halo_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)grid_function_mod.o
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

$(DOBJ)test_halo_main.o: src/test/test_halo/test_halo_main.f90 \
	$(DOBJ)test_halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_halo_mod.o: src/test/test_halo/test_halo_mod.f90 \
	$(DOBJ)grid_function_mod.o \
	$(DOBJ)exchange_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)mesh_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)ecs_halo_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_mesh_main.o: src/test/test_mesh/test_mesh_main.f90 \
	$(DOBJ)test_mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_mesh_mod.o: src/test/test_mesh/test_mesh_mod.f90 \
	$(DOBJ)partition_mod.o \
	$(DOBJ)mesh_factory_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_exch_main.o: src/test/test_exch/test_exch_main.f90 \
	$(DOBJ)test_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_mod.o: src/test/test_exch/test_mod.f90 \
	$(DOBJ)grid_function_mod.o \
	$(DOBJ)exchange_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)exchange_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_metric_main.o: src/test/test_metric/test_metric_main.f90 \
	$(DOBJ)test_metric_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_metric_mod.o: src/test/test_metric/test_metric_mod.f90 \
	$(DOBJ)grid_function_mod.o \
	$(DOBJ)exchange_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)ecs_geometry_mod.o
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
