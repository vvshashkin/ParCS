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
$(DEXE)TEST_DOMAIN_MAIN: $(MKDIRS) $(DOBJ)test_domain_main.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_domain_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_DOMAIN_MAIN
$(DEXE)TEST_PANELED_OUTPUT: $(MKDIRS) $(DOBJ)test_paneled_output.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_paneled_output.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_PANELED_OUTPUT
$(DEXE)TEST_NAMELIST: $(MKDIRS) $(DOBJ)test_namelist.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_namelist.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_NAMELIST
$(DEXE)TEST_TS: $(MKDIRS) $(DOBJ)test_ts.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_ts.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_TS
$(DEXE)TEST_EXCH_MAIN: $(MKDIRS) $(DOBJ)test_exch_main.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_exch_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_EXCH_MAIN
$(DEXE)TEST_CMD_LINE: $(MKDIRS) $(DOBJ)test_cmd_line.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_cmd_line.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_CMD_LINE
$(DEXE)TEST_DIFFOPS_ALL: $(MKDIRS) $(DOBJ)test_diffops_all.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_diffops_all.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_DIFFOPS_ALL
$(DEXE)TEST_DIFFOPS_CONV: $(MKDIRS) $(DOBJ)test_diffops_conv.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_diffops_conv.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_DIFFOPS_CONV
$(DEXE)TEST_DIFFOPS: $(MKDIRS) $(DOBJ)test_diffops.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_diffops.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_DIFFOPS
$(DEXE)TEST_LAPLACE_SPECTRE: $(MKDIRS) $(DOBJ)test_laplace_spectre.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_laplace_spectre.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_LAPLACE_SPECTRE
$(DEXE)TEST_HALO_MAIN: $(MKDIRS) $(DOBJ)test_halo_main.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_halo_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_HALO_MAIN
$(DEXE)TEST_GRID_FIELD: $(MKDIRS) $(DOBJ)test_grid_field.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_grid_field.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_GRID_FIELD
$(DEXE)TEST_METRIC_MAIN: $(MKDIRS) $(DOBJ)test_metric_main.o \
	$(DOBJ)avost.o \
	$(DOBJ)auxhs.o
	@rm -f $(filter-out $(DOBJ)test_metric_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_METRIC_MAIN

#compiling rules
$(DOBJ)tile_mod.o: src/tile_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)global_diag_mod.o: src/global_diag_mod.f90 \
	$(DOBJ)container_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grid_field_mod.o: src/grid_field_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_mod.o: src/stvec_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)namelist_read_mod.o: src/namelist_read_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_abstract_mod.o: src/operator_abstract_mod.f90 \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)container_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)avost.o: src/avost.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grid_field_factory_mod.o: src/grid_field_factory_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)const_mod.o: src/const_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mesh_factory_mod.o: src/mesh_factory_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)metric_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_mod.o: src/operator_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mesh_mod.o: src/mesh_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)tiles_mod.o: src/tiles_mod.f90 \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sph_coords_mod.o: src/sph_coords_mod.f90 \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)cmd_args_mod.o: src/cmd_args_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)container_abstract_mod.o: src/container_abstract_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_abstract_mod.o: src/stvec_abstract_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)buffer_mod.o: src/parallel/buffer_mod.f90 \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_abstract_mod.o: src/parallel/exchange_abstract_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)partition_mod.o: src/parallel/partition_mod.f90 \
	$(DOBJ)tile_mod.o \
	$(DOBJ)tiles_mod.o \
	$(DOBJ)parcomm_mod.o
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

$(DOBJ)parcomm_factory_mod.o: src/parallel/parcomm_factory_mod.f90 \
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
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_gather_mod.o
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

$(DOBJ)parcomm_mod.o: src/parallel/parcomm_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)partition_factory_mod.o: src/parallel/partition_factory_mod.f90 \
	$(DOBJ)partition_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_gather_mod.o: src/parallel/exchange_gather_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)buffer_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)metric_mod.o: src/metric/metric_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)metric_factory_mod.o: src/metric/metric_factory_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)ecs_metric_mod.o \
	$(DOBJ)ecs_metric_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)outputer_factory_mod.o: src/outputer/outputer_factory_mod.f90 \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)master_paneled_outputer_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mpi_paneled_outputer_mod.o: src/outputer/mpi_paneled_outputer_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)outputer_abstract_mod.o: src/outputer/outputer_abstract_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)master_paneled_outputer_mod.o: src/outputer/master_paneled_outputer_mod.f90 \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)tiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)halo_mod.o: src/halo/halo_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)halo_factory_mod.o: src/halo/halo_factory_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)ecs_halo_factory_mod.o \
	$(DOBJ)ecs_halo_vec_a_factory_mod.o \
	$(DOBJ)ecs_halo_vec_c_factory_mod.o \
	$(DOBJ)ecs_ah_vec_sync_factory_mod.o \
	$(DOBJ)halo_a_default_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)halo_c_default_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)halo_c_default_mod.o: src/halo/halo_C_default_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)halo_a_default_mod.o: src/halo/halo_A_default_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_fields_mod.o: src/test_fields/test_fields_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)latlon_functions_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)latlon_functions_mod.o: src/test_fields/latlon_functions_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)timescheme_factory_mod.o: src/time_schemes/timescheme_factory_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)explicit_eul1_mod.o \
	$(DOBJ)rk4_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)rk4_mod.o: src/time_schemes/rk4_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)timescheme_mod.o: src/time_schemes/timescheme_mod.f90 \
	$(DOBJ)operator_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)explicit_eul1_mod.o: src/time_schemes/explicit_Eul1_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)domain_factory_mod.o: src/domain/domain_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)topology_factory_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)metric_factory_mod.o \
	$(DOBJ)mesh_factory_mod.o \
	$(DOBJ)parcomm_factory_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)domain_mod.o: src/domain/domain_mod.f90 \
	$(DOBJ)topology_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)curl_factory_mod.o: src/differential_operators/curl_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)curl_div_based_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_a2_mod.o: src/differential_operators/div_a2_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_grad_mod.o: src/differential_operators/abstract_grad_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_contra_ah_sbp_mod.o: src/differential_operators/grad_contra_ah_sbp_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_factory_mod.o: src/differential_operators/div_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)div_c2_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)div_c_sbp42_mod.o \
	$(DOBJ)div_a2_mod.o \
	$(DOBJ)div_ah2_mod.o \
	$(DOBJ)div_ah_sbp_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_contra_ah2_mod.o: src/differential_operators/grad_contra_ah2_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_contra_a2_mod.o: src/differential_operators/grad_contra_a2_mod.f90 \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_ah2_mod.o: src/differential_operators/div_ah2_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_curl_mod.o: src/differential_operators/abstract_curl_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)curl_div_based_mod.o: src/differential_operators/curl_div_based_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_ah_sbp_mod.o: src/differential_operators/div_ah_sbp_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_c_sbp42_mod.o: src/differential_operators/div_c_sbp42_mod.f90 \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_c2_mod.o: src/differential_operators/div_c2_mod.f90 \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_contra_c_sbp42_mod.o: src/differential_operators/grad_contra_c_sbp42_mod.f90 \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)sbp_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_contra_c2_ecs_mod.o: src/differential_operators/grad_contra_c2_ecs_mod.f90 \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sbp_mod.o: src/differential_operators/sbp_mod.f90 \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_div_mod.o: src/differential_operators/abstract_div_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_factory_mod.o: src/differential_operators/grad_factory_mod.f90 \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grad_contra_c2_ecs_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)grad_contra_c_sbp42_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)grad_contra_a2_mod.o \
	$(DOBJ)grad_contra_ah2_mod.o \
	$(DOBJ)ecs_metric_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)grad_contra_ah_sbp_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)auxhs.o: src/aux/auxhs.f
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
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_metric_mod.o: src/equiang_cs/ecs_metric_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_vec_c_mod.o: src/equiang_cs/ecs_halo_vec_c_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_ah_vec_sync_mod.o: src/equiang_cs/ecs_Ah_vec_sync_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_ah_vec_sync_factory_mod.o: src/equiang_cs/ecs_Ah_vec_sync_factory_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)ecs_ah_vec_sync_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)ecs_metric_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_geometry_mod.o: src/equiang_cs/ecs_geometry_mod.f90 \
	$(DOBJ)topology_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_vec_c_factory_mod.o: src/equiang_cs/ecs_halo_vec_c_factory_mod.f90 \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)ecs_halo_vec_c_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)ecs_metric_mod.o
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

$(DOBJ)ecs_halo_vec_a_factory_mod.o: src/equiang_cs/ecs_halo_vec_a_factory_mod.f90 \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)ecs_halo_vec_a_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)ecs_halo_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)ecs_metric_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_vec_a_mod.o: src/equiang_cs/ecs_halo_vec_a_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_mod.o: src/equiang_cs/ecs_halo_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_metric_factory_mod.o: src/equiang_cs/ecs_metric_factory_mod.f90 \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)ecs_metric_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_iomega_mod.o: src/models/iomega_model/stvec_iomega_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parameters_iomega_mod.o: src/models/iomega_model/parameters_iomega_mod.f90 \
	$(DOBJ)container_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_iomega_mod.o: src/models/iomega_model/operator_iomega_mod.f90 \
	$(DOBJ)operator_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_iomega_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_domain_main.o: src/test/test_domain/test_domain_main.f90 \
	$(DOBJ)test_domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_domain_mod.o: src/test/test_domain/test_domain_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_paneled_output.o: src/test/test_paneled_output/test_paneled_output.f90 \
	$(DOBJ)test_paneled_output_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_paneled_output_mod.o: src/test/test_paneled_output/test_paneled_output_mod.f90 \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_namelist_mod.o: src/test/test_namelist/test_namelist_mod.f90 \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_namelist.o: src/test/test_namelist/test_namelist.f90 \
	$(DOBJ)test_namelist_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_ts.o: src/test/test_time_steping/test_ts.f90 \
	$(DOBJ)test_ts_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_ts_mod.o: src/test/test_time_steping/test_ts_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_iomega_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_iomega_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_mod.o: src/test/test_exch/test_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)exchange_halo_c_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_exch_main.o: src/test/test_exch/test_exch_main.f90 \
	$(DOBJ)test_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)cmd_args_test_mod.o: src/test/test_cmd_line/cmd_args_test_mod.f90 \
	$(DOBJ)cmd_args_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_cmd_line.o: src/test/test_cmd_line/test_cmd_line.f90 \
	$(DOBJ)cmd_args_test_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_diffops_all.o: src/test/test_diff_ops/test_diffops_all.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_diffops_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_diffops_conv.o: src/test/test_diff_ops/test_diffops_conv.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_diffops_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_diffops.o: src/test/test_diff_ops/test_diffops.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_diffops_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_diffops_mod.o: src/test/test_diff_ops/test_diffops_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)curl_factory_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_laplace_spectre.o: src/test/test_diff_ops/test_laplace_spectre.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_diffops_mod.o \
	$(DOBJ)cmd_args_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_halo_main.o: src/test/test_halo/test_halo_main.f90 \
	$(DOBJ)test_ecs_halo_mod.o \
	$(DOBJ)test_ecs_halo_c_mod.o \
	$(DOBJ)test_halo_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_ecs_halo_c_mod.o: src/test/test_halo/test_ecs_halo_c_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)mesh_mod.o
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
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_grid_field.o: src/test/test_grid_field/test_grid_field.f90 \
	$(DOBJ)test_grid_field_mod.o \
	$(DOBJ)parcomm_mod.o
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

$(DOBJ)test_metric_class_mod.o: src/test/test_metric/test_metric_class_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)metric_factory_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)topology_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_metric_main.o: src/test/test_metric/test_metric_main.f90 \
	$(DOBJ)test_metric_mod.o \
	$(DOBJ)test_metric_class_mod.o \
	$(DOBJ)parcomm_mod.o
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

$(DOBJ)timescheme_abstract_mod.o: src/time_schemes_old/timescheme_abstract_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ars343.o: src/time_schemes_old/ars343.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exp_taylor_mod.o: src/time_schemes_old/exp_taylor_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exp_krylov_mod.o: src/time_schemes_old/exp_krylov_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)stvec_abstract_mod.o \
	$(DOBJ)timescheme_abstract_mod.o \
	$(DOBJ)operator_abstract_mod.o
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
