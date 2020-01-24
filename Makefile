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
$(DEXE)TEST_NAMELIST: $(MKDIRS) $(DOBJ)test_namelist.o
	@rm -f $(filter-out $(DOBJ)test_namelist.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_NAMELIST
$(DEXE)TEST_OUTPUT_MAIN: $(MKDIRS) $(DOBJ)test_output_main.o
	@rm -f $(filter-out $(DOBJ)test_output_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_OUTPUT_MAIN
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
$(DEXE)TEST_TS: $(MKDIRS) $(DOBJ)test_ts.o
	@rm -f $(filter-out $(DOBJ)test_ts.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_TS
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

$(DOBJ)namelist_read_mod.o: src/namelist_read_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mesh_mod.o: src/mesh_mod.f90 \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grid_function_mod.o: src/grid_function_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_abstract_mod.o: src/stvec_abstract_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parcs_mpi_mod.o: src/ParCS_mpi_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ouputer_factory_mod.o: src/ouputer_factory_mod.f90 \
	$(DOBJ)master_process_outputer_mod.o \
	$(DOBJ)exchange_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)master_process_outputer_mod.o: src/master_process_outputer_mod.f90 \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)grid_function_mod.o \
	$(DOBJ)exchange_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_abstract_mod.o: src/operator_abstract_mod.f90 \
	$(DOBJ)stvec_abstract_mod.o
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
	$(DOBJ)topology_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)ecs_geometry_mod.o \
	$(DOBJ)ecs_halo_factory_mod.o
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

$(DOBJ)outputer_abstract_mod.o: src/outputer_abstract_mod.f90 \
	$(DOBJ)grid_function_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_iomega_mod.o: src/iomega_model/stvec_iomega_mod.f90 \
	$(DOBJ)stvec_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_iomega_mod.o: src/iomega_model/operator_iomega_mod.f90 \
	$(DOBJ)operator_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)stvec_iomega_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_vec_a_mod.o: src/equiang_cs/ecs_halo_vec_a_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)grid_function_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_geometry_mod.o: src/equiang_cs/ecs_geometry_mod.f90 \
	$(DOBJ)topology_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_factory_mod.o: src/equiang_cs/ecs_halo_factory_mod.f90 \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_vec_a_factory_mod.o: src/equiang_cs/ecs_halo_vec_a_factory_mod.f90 \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)ecs_halo_vec_a_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)ecs_geometry_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_mod.o: src/equiang_cs/ecs_halo_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)grid_function_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)explicit_eul1_mod.o: src/time_schemes/explicit_Eul1_mod.f90 \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)timescheme_abstract_mod.o: src/time_schemes/timescheme_abstract_mod.f90 \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)rk4_mod.o: src/time_schemes/rk4_mod.f90 \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_namelist.o: src/test/test_namelist/test_namelist.f90 \
	$(DOBJ)test_namelist_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_namelist_mod.o: src/test/test_namelist/test_namelist_mod.f90 \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_output_main.o: src/test/test_output/test_output_main.f90 \
	$(DOBJ)test_output_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_output_mod.o: src/test/test_output/test_output_mod.f90 \
	$(DOBJ)grid_function_mod.o \
	$(DOBJ)exchange_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)master_process_outputer_mod.o \
	$(DOBJ)ouputer_factory_mod.o
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
	$(DOBJ)ecs_halo_factory_mod.o \
	$(DOBJ)ecs_halo_vec_a_factory_mod.o
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

$(DOBJ)test_rk4.o: src/test/test_time_steping/test_rk4.f90 \
	$(DOBJ)stvec_iomega_mod.o \
	$(DOBJ)operator_iomega_mod.o \
	$(DOBJ)explicit_eul1_mod.o \
	$(DOBJ)rk4_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_ts.o: src/test/test_time_steping/test_ts.f90 \
	$(DOBJ)test_rk4.o
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
	$(DOBJ)mesh_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)topology_mod.o
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
all: $(addprefix $(DEXE),$(EXES))
