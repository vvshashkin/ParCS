#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = obj/
DMOD    = mod/
DEXE    = ./
LIBS    =
FC      = mpiifort
OPTSC   =  -c -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -module mod
OPTSL   =  -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -lmkl_intel_lp64 -lmkl_core -lmkl_gf_lp64 -lmkl_sequential -lmkl_lapack95_lp64 -module mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)SWLIN: $(MKDIRS) $(DOBJ)swlin.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)swlin.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) SWLIN
$(DEXE)NHLIN: $(MKDIRS) $(DOBJ)nhlin.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)nhlin.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) NHLIN
$(DEXE)TEST_CMD_LINE: $(MKDIRS) $(DOBJ)test_cmd_line.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_cmd_line.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_CMD_LINE
$(DEXE)TEST_NAMELIST: $(MKDIRS) $(DOBJ)test_namelist.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_namelist.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_NAMELIST
$(DEXE)TEST_GRID_FIELD: $(MKDIRS) $(DOBJ)test_grid_field.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_grid_field.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_GRID_FIELD
$(DEXE)TEST_EXCH_MAIN: $(MKDIRS) $(DOBJ)test_exch_main.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_exch_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_EXCH_MAIN
$(DEXE)TEST_PANELED_OUTPUT: $(MKDIRS) $(DOBJ)test_paneled_output.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_paneled_output.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_PANELED_OUTPUT
$(DEXE)TEST_DOMAIN_MAIN: $(MKDIRS) $(DOBJ)test_domain_main.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_domain_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_DOMAIN_MAIN
$(DEXE)TEST_HALO_MAIN: $(MKDIRS) $(DOBJ)test_halo_main.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_halo_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_HALO_MAIN
$(DEXE)TEST_TS: $(MKDIRS) $(DOBJ)test_ts.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_ts.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_TS
$(DEXE)TEST_GLOBAL_DIAG_MAIN: $(MKDIRS) $(DOBJ)test_global_diag_main.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_global_diag_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_GLOBAL_DIAG_MAIN
$(DEXE)TEST_METRIC_MAIN: $(MKDIRS) $(DOBJ)test_metric_main.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_metric_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_METRIC_MAIN
$(DEXE)TEST_MESH_MAIN: $(MKDIRS) $(DOBJ)test_mesh_main.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_mesh_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_MESH_MAIN

#compiling rules
$(DOBJ)grid_function_mod.o: src/grid_function_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)namelist_read_mod.o: src/namelist_read_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)tile_mod.o: src/tile_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grid_field_factory_mod.o: src/grid_field_factory_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mesh_factory_mod.o: src/mesh_factory_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)metric_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mesh_mod.o: src/mesh_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grid_field_mod.o: src/grid_field_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)container_abstract_mod.o: src/container_abstract_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)avost.o: src/avost.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)global_diag_mod.o: src/global_diag_mod.f90 \
	$(DOBJ)container_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)const_mod.o: src/const_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)cmd_args_mod.o: src/cmd_args_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_abstract_mod.o: src/operator_abstract_mod.f90 \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)container_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_abstract_mod.o: src/stvec_abstract_mod.f90 \
	$(DOBJ)container_abstract_mod.o
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

$(DOBJ)ecs_geometry_mod.o: src/equiang_cs/ecs_geometry_mod.f90 \
	$(DOBJ)topology_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_metric_factory_mod.o: src/equiang_cs/ecs_metric_factory_mod.f90 \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)ecs_metric_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_metric_mod.o: src/equiang_cs/ecs_metric_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_factory_mod.o: src/equiang_cs/ecs_halo_factory_mod.f90 \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_mod.o: src/equiang_cs/ecs_halo_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_vec_a_mod.o: src/equiang_cs/ecs_halo_vec_a_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)auxhs.o: src/aux/auxhs.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_halo_c_mod.o: src/parallel/exchange_halo_C_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)buffer_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)buffer_mod.o: src/parallel/buffer_mod.f90 \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)partition_mod.o: src/parallel/partition_mod.f90 \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parcomm_factory_mod.o: src/parallel/parcomm_factory_mod.f90 \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_gather_mod.o: src/parallel/exchange_gather_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)buffer_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_factory_mod.o: src/parallel/exchange_factory_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)exchange_halo_c_mod.o \
	$(DOBJ)exchange_gather_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_abstract_mod.o: src/parallel/exchange_abstract_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parcomm_mod.o: src/parallel/parcomm_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_halo_mod.o: src/parallel/exchange_halo_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)buffer_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)metric_factory_mod.o: src/metric/metric_factory_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)ecs_metric_mod.o \
	$(DOBJ)ecs_metric_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)metric_mod.o: src/metric/metric_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)cubed_sphere_topology_mod.o: src/topology/cubed_sphere_topology_mod.f90 \
	$(DOBJ)topology_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)topology_mod.o: src/topology/topology_mod.f90 \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)topology_factory_mod.o: src/topology/topology_factory_mod.f90 \
	$(DOBJ)topology_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)halo_mod.o: src/halo/halo_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)halo_factory_mod.o: src/halo/halo_factory_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)ecs_halo_factory_mod.o \
	$(DOBJ)halo_a_default_mod.o \
	$(DOBJ)exchange_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)halo_a_default_mod.o: src/halo/halo_A_default_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)swlin.o: src/models/linear_shallow_water/swlin.f90 \
	$(DOBJ)swlin_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)swlin_initial_cond_mod.o: src/models/linear_shallow_water/swlin_initial_cond_mod.f90 \
	$(DOBJ)stvec_swlin_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)swlin_operator_factory_mod.o: src/models/linear_shallow_water/swlin_operator_factory_mod.f90 \
	$(DOBJ)partition_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)ecs_halo_factory_mod.o \
	$(DOBJ)hor_difops_basic_mod.o \
	$(DOBJ)parameters_swlin_mod.o \
	$(DOBJ)operator_swlin_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_swlin_factory_mod.o: src/models/linear_shallow_water/stvec_swlin_factory_mod.f90 \
	$(DOBJ)stvec_swlin_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parameters_swlin_mod.o: src/models/linear_shallow_water/parameters_swlin_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)mesh_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)swlin_output_mod.o: src/models/linear_shallow_water/swlin_output_mod.f90 \
	$(DOBJ)grid_function_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)stvec_swlin_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)swlin_mod.o: src/models/linear_shallow_water/swlin_mod.f90 \
	$(DOBJ)parameters_swlin_mod.o \
	$(DOBJ)stvec_swlin_mod.o \
	$(DOBJ)operator_swlin_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)diag_swlin_mod.o \
	$(DOBJ)cmd_args_mod.o \
	$(DOBJ)namelist_read_mod.o \
	$(DOBJ)swlin_output_mod.o \
	$(DOBJ)swlin_initial_cond_mod.o \
	$(DOBJ)swlin_operator_factory_mod.o \
	$(DOBJ)stvec_swlin_factory_mod.o \
	$(DOBJ)rk4_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_swlin_mod.o: src/models/linear_shallow_water/stvec_swlin_mod.f90 \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)diag_swlin_mod.o: src/models/linear_shallow_water/diag_swlin_mod.f90 \
	$(DOBJ)global_diag_mod.o \
	$(DOBJ)stvec_swlin_mod.o \
	$(DOBJ)parameters_swlin_mod.o \
	$(DOBJ)container_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_swlin_mod.o: src/models/linear_shallow_water/operator_swlin_mod.f90 \
	$(DOBJ)operator_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)stvec_swlin_mod.o \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)parameters_swlin_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)hor_difops_abstract_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nhlin_mod.o: src/models/linear_NH/NHlin_mod.f90 \
	$(DOBJ)parameters_nhlin_mod.o \
	$(DOBJ)stvec_nhlin_mod.o \
	$(DOBJ)operator_nhlin_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)diag_nhlin_mod.o \
	$(DOBJ)cmd_args_mod.o \
	$(DOBJ)namelist_read_mod.o \
	$(DOBJ)nhlin_output_mod.o \
	$(DOBJ)nhlin_initial_cond_mod.o \
	$(DOBJ)tscheme_nhlin_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parameters_nhlin_mod.o: src/models/linear_NH/parameters_NHlin_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)mesh_factory_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)tscheme_nhlin_mod.o: src/models/linear_NH/tscheme_NHlin_mod.f90 \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)operator_nhlin_mod.o \
	$(DOBJ)operator_nhlin_split_mod.o \
	$(DOBJ)parameters_nhlin_mod.o \
	$(DOBJ)stvec_nhlin_mod.o \
	$(DOBJ)rk4_mod.o \
	$(DOBJ)exp_krylov_mod.o \
	$(DOBJ)ars343.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nhlin_initial_cond_mod.o: src/models/linear_NH/NHlin_initial_cond_mod.f90 \
	$(DOBJ)stvec_nhlin_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parameters_nhlin_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_nhlin_mod.o: src/models/linear_NH/stvec_NHlin_mod.f90 \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)grid_function_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)diag_nhlin_mod.o: src/models/linear_NH/diag_NHlin_mod.f90 \
	$(DOBJ)global_diag_mod.o \
	$(DOBJ)stvec_nhlin_mod.o \
	$(DOBJ)parameters_nhlin_mod.o \
	$(DOBJ)container_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_nhlin_split_mod.o: src/models/linear_NH/operator_NHlin_split_mod.f90 \
	$(DOBJ)operator_nhlin_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)stvec_nhlin_mod.o \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)parameters_nhlin_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)hor_difops_abstract_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)grid_function_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_nhlin_mod.o: src/models/linear_NH/operator_NHlin_mod.f90 \
	$(DOBJ)operator_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)stvec_nhlin_mod.o \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)parameters_nhlin_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)hor_difops_abstract_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)ecs_halo_factory_mod.o \
	$(DOBJ)hor_difops_basic_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nhlin_output_mod.o: src/models/linear_NH/NHlin_output_mod.f90 \
	$(DOBJ)partition_mod.o \
	$(DOBJ)grid_function_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)stvec_nhlin_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nhlin.o: src/models/linear_NH/NHlin.f90 \
	$(DOBJ)nhlin_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_iomega_mod.o: src/models/iomega_model/operator_iomega_mod.f90 \
	$(DOBJ)operator_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)stvec_iomega_mod.o \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)parameters_iomega_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_iomega_mod.o: src/models/iomega_model/stvec_iomega_mod.f90 \
	$(DOBJ)stvec_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parameters_iomega_mod.o: src/models/iomega_model/parameters_iomega_mod.f90 \
	$(DOBJ)container_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mpi_paneled_outputer_mod.o: src/outputer/mpi_paneled_outputer_mod.f90 \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)partition_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)outputer_factory_mod.o: src/outputer/outputer_factory_mod.f90 \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)master_paneled_outputer_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)mpi_paneled_outputer_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)master_paneled_outputer_mod.o: src/outputer/master_paneled_outputer_mod.f90 \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)outputer_abstract_mod.o: src/outputer/outputer_abstract_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)partition_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)domain_mod.o: src/domain/domain_mod.f90 \
	$(DOBJ)topology_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)domain_factory_mod.o: src/domain/domain_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)topology_factory_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)metric_factory_mod.o \
	$(DOBJ)mesh_factory_mod.o \
	$(DOBJ)parcomm_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hor_difops_abstract_mod.o: src/differential_operators/hor_difops_abstract_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hor_difops_basic_mod.o: src/differential_operators/hor_difops_basic_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ars343.o: src/time_schemes/ars343.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)explicit_eul1_mod.o: src/time_schemes/explicit_Eul1_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)rk4_mod.o: src/time_schemes/rk4_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)timescheme_abstract_mod.o: src/time_schemes/timescheme_abstract_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exp_krylov_mod.o: src/time_schemes/exp_krylov_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exp_taylor_mod.o: src/time_schemes/exp_taylor_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_cmd_line.o: src/test/test_cmd_line/test_cmd_line.f90 \
	$(DOBJ)cmd_args_test_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)cmd_args_test_mod.o: src/test/test_cmd_line/cmd_args_test_mod.f90 \
	$(DOBJ)cmd_args_mod.o
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

$(DOBJ)test_grid_field_mod.o: src/test/test_grid_field/test_grid_field_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_grid_field.o: src/test/test_grid_field/test_grid_field.f90 \
	$(DOBJ)test_grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_exch_main.o: src/test/test_exch/test_exch_main.f90 \
	$(DOBJ)test_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_mod.o: src/test/test_exch/test_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)exchange_halo_c_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_paneled_output.o: src/test/test_paneled_output/test_paneled_output.f90 \
	$(DOBJ)test_paneled_output_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_paneled_output_mod.o: src/test/test_paneled_output/test_paneled_output_mod.f90 \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)mesh_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_domain_mod.o: src/test/test_domain/test_domain_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_domain_main.o: src/test/test_domain/test_domain_main.f90 \
	$(DOBJ)test_domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_halo_main.o: src/test/test_halo/test_halo_main.f90 \
	$(DOBJ)test_ecs_halo_mod.o \
	$(DOBJ)test_halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_halo_mod.o: src/test/test_halo/test_halo_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_ecs_halo_mod.o: src/test/test_halo/test_ecs_halo_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_ts_mod.o: src/test/test_time_steping/test_ts_mod.f90 \
	$(DOBJ)stvec_iomega_mod.o \
	$(DOBJ)operator_iomega_mod.o \
	$(DOBJ)parameters_iomega_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)explicit_eul1_mod.o \
	$(DOBJ)rk4_mod.o \
	$(DOBJ)exp_taylor_mod.o \
	$(DOBJ)exp_krylov_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_ts.o: src/test/test_time_steping/test_ts.f90 \
	$(DOBJ)test_ts_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_gl_diag_mod.o: src/test/test_global_diagnostics/test_gl_diag_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parameters_swlin_mod.o \
	$(DOBJ)stvec_swlin_mod.o \
	$(DOBJ)global_diag_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_global_diag_main.o: src/test/test_global_diagnostics/test_global_diag_main.f90 \
	$(DOBJ)test_gl_diag_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_metric_mod.o: src/test/test_metric/test_metric_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)mesh_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_metric_main.o: src/test/test_metric/test_metric_main.f90 \
	$(DOBJ)test_metric_mod.o \
	$(DOBJ)test_metric_class_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_metric_class_mod.o: src/test/test_metric/test_metric_class_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)metric_factory_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)topology_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_mesh_mod.o: src/test/test_mesh/test_mesh_mod.f90 \
	$(DOBJ)partition_mod.o \
	$(DOBJ)mesh_factory_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_mesh_main.o: src/test/test_mesh/test_mesh_main.f90 \
	$(DOBJ)test_mesh_mod.o
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
