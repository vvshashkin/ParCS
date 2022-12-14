#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = obj/
DMOD    = mod/
DEXE    = ./
LIBS    =
FC      = mpiifort
OPTSC   =  -c -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict -module mod
OPTSL   =  -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict  -module mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)RH4_WAVE_MAIN: $(MKDIRS) $(DOBJ)rh4_wave_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)rh4_wave_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) RH4_WAVE_MAIN
$(DEXE)TS2_MAIN: $(MKDIRS) $(DOBJ)ts2_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)ts2_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TS2_MAIN
$(DEXE)FORCING_TEST_MAIN: $(MKDIRS) $(DOBJ)forcing_test_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)forcing_test_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) FORCING_TEST_MAIN
$(DEXE)TS5_MAIN: $(MKDIRS) $(DOBJ)ts5_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)ts5_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TS5_MAIN
$(DEXE)BAROTROPIC_INST_MAIN: $(MKDIRS) $(DOBJ)barotropic_inst_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)barotropic_inst_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) BAROTROPIC_INST_MAIN
$(DEXE)ELDRED_TEST_MAIN: $(MKDIRS) $(DOBJ)eldred_test_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)eldred_test_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) ELDRED_TEST_MAIN
$(DEXE)NH_MAIN: $(MKDIRS) $(DOBJ)nh_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)nh_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) NH_MAIN
$(DEXE)TEST_TS: $(MKDIRS) $(DOBJ)test_ts.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_ts.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_TS
$(DEXE)TEST_VERTICAL_OPERATORS_MAIN: $(MKDIRS) $(DOBJ)test_vertical_operators_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_vertical_operators_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_VERTICAL_OPERATORS_MAIN
$(DEXE)TEST_CMD_LINE: $(MKDIRS) $(DOBJ)test_cmd_line.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_cmd_line.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_CMD_LINE
$(DEXE)TEST_NAMELIST: $(MKDIRS) $(DOBJ)test_namelist.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_namelist.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_NAMELIST
$(DEXE)TEST_EXCH_MAIN: $(MKDIRS) $(DOBJ)test_exch_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_exch_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_EXCH_MAIN
$(DEXE)TEST_GRID_FIELD: $(MKDIRS) $(DOBJ)test_grid_field.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_grid_field.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_GRID_FIELD
$(DEXE)TEST_DIFFOPS_CONV: $(MKDIRS) $(DOBJ)test_diffops_conv.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_diffops_conv.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_DIFFOPS_CONV
$(DEXE)TEST_DIFFOPS: $(MKDIRS) $(DOBJ)test_diffops.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_diffops.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_DIFFOPS
$(DEXE)TEST_DIFFOPS_ALL: $(MKDIRS) $(DOBJ)test_diffops_all.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_diffops_all.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_DIFFOPS_ALL
$(DEXE)TEST_DIFFOPS_3D: $(MKDIRS) $(DOBJ)test_diffops_3d.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_diffops_3d.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_DIFFOPS_3D
$(DEXE)TEST_LAPLACE_SPECTRE: $(MKDIRS) $(DOBJ)test_laplace_spectre.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_laplace_spectre.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_LAPLACE_SPECTRE
$(DEXE)TEST_HALO_MAIN: $(MKDIRS) $(DOBJ)test_halo_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_halo_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_HALO_MAIN
$(DEXE)TEST_REGRID: $(MKDIRS) $(DOBJ)test_regrid.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_regrid.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_REGRID
$(DEXE)TEST_LATLON_OUTPUT: $(MKDIRS) $(DOBJ)test_latlon_output.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_latlon_output.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_LATLON_OUTPUT
$(DEXE)TEST_VERTICAL_PROFILES_MAIN: $(MKDIRS) $(DOBJ)test_vertical_profiles_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_vertical_profiles_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_VERTICAL_PROFILES_MAIN
$(DEXE)TEST_METRIC_MAIN: $(MKDIRS) $(DOBJ)test_metric_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_metric_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_METRIC_MAIN
$(DEXE)TEST_DOMAIN_MAIN: $(MKDIRS) $(DOBJ)test_domain_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_domain_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_DOMAIN_MAIN
$(DEXE)TEST_ADVECTION_MAIN: $(MKDIRS) $(DOBJ)test_advection_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_advection_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_ADVECTION_MAIN
$(DEXE)TEST_HORDIFF_MAIN: $(MKDIRS) $(DOBJ)test_hordiff_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_hordiff_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_HORDIFF_MAIN
$(DEXE)TEST_VERTICAL_TRANSFORM_MAIN: $(MKDIRS) $(DOBJ)test_vertical_transform_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_vertical_transform_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_VERTICAL_TRANSFORM_MAIN
$(DEXE)TEST_QUADRATURE_MAIN: $(MKDIRS) $(DOBJ)test_quadrature_main.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_quadrature_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_QUADRATURE_MAIN
$(DEXE)TEST_PANELED_OUTPUT: $(MKDIRS) $(DOBJ)test_paneled_output.o \
	$(DOBJ)avost.o
	@rm -f $(filter-out $(DOBJ)test_paneled_output.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_PANELED_OUTPUT

#compiling rules
$(DOBJ)stvec_mod.o: src/stvec_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grid_field_mod.o: src/grid_field_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)namelist_read_mod.o: src/namelist_read_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)tiles_mod.o: src/tiles_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)container_abstract_mod.o: src/container_abstract_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mesh_mod.o: src/mesh_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)const_mod.o: src/const_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sph_coords_mod.o: src/sph_coords_mod.f90 \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)global_diag_mod.o: src/global_diag_mod.f90 \
	$(DOBJ)container_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_mod.o: src/operator_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)key_value_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grid_field_factory_mod.o: src/grid_field_factory_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)avost.o: src/avost.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mesh_factory_mod.o: src/mesh_factory_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tiles_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)orography_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vec_math_mod.o: src/vec_math_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)cmd_args_mod.o: src/cmd_args_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_abstract_mod.o: src/stvec_abstract_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)tile_mod.o: src/tile_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_mod.o: src/config_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_vertical_transform_mod.o: src/metric/abstract_vertical_transform_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_metric_mod.o: src/metric/config_metric_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vertical_transform_factory_mod.o: src/metric/vertical_transform_factory_mod.f90 \
	$(DOBJ)abstract_vertical_transform_mod.o \
	$(DOBJ)vertical_transform_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)metric_mod.o: src/metric/metric_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)orography_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)metric_factory_mod.o: src/metric/metric_factory_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_metric_mod.o \
	$(DOBJ)ecs_metric_mod.o \
	$(DOBJ)ecs_metric_factory_mod.o \
	$(DOBJ)shallow_atm_metric_mod.o \
	$(DOBJ)vertical_transform_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vertical_transform_mod.o: src/metric/vertical_transform_mod.f90 \
	$(DOBJ)abstract_vertical_transform_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)shallow_atm_metric_mod.o: src/metric/shallow_atm_metric_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)abstract_vertical_transform_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)orography_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)key_value_mod.o: src/stuff/key_value_mod.f90
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

$(DOBJ)stvec_iomega_mod.o: src/models/iomega_model/stvec_iomega_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_swm_mod.o: src/models/shallow_water/stvec/stvec_swm_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_swm_factory_mod.o: src/models/shallow_water/stvec/stvec_swm_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)rh4_wave_mod.o: src/models/shallow_water/test/RH4_wave/RH4_wave_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_swm_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)rh4_wave_main.o: src/models/shallow_water/test/RH4_wave/RH4_wave_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)rh4_wave_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ts2_mod.o: src/models/shallow_water/test/ts2/ts2_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_swm_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)config_ts2_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_ts2_mod.o: src/models/shallow_water/test/ts2/config_ts2_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ts2_main.o: src/models/shallow_water/test/ts2/ts2_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)ts2_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)random_friction_mod.o: src/models/shallow_water/test/forcing_test/random_friction_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)forcing_test_mod.o: src/models/shallow_water/test/forcing_test/forcing_test_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_swm_mod.o \
	$(DOBJ)operator_swm_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)random_friction_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)forcing_test_main.o: src/models/shallow_water/test/forcing_test/forcing_test_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)forcing_test_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ts5_mod.o: src/models/shallow_water/test/ts5/ts5_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_swm_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)config_ts5_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ts5_main.o: src/models/shallow_water/test/ts5/ts5_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)ts5_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_ts5_mod.o: src/models/shallow_water/test/ts5/config_ts5_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)barotropic_instability_u_mod.o: src/models/shallow_water/test/barotropic_instability/barotropic_instability_u_mod.f90 \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)barotropic_inst_main.o: src/models/shallow_water/test/barotropic_instability/barotropic_inst_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)barotropic_inst_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)barotropic_inst_mod.o: src/models/shallow_water/test/barotropic_instability/barotropic_inst_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_swm_mod.o \
	$(DOBJ)config_barotropic_inst_mod.o \
	$(DOBJ)operator_swm_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)curl_factory_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o \
	$(DOBJ)barotropic_instability_u_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_barotropic_inst_mod.o: src/models/shallow_water/test/barotropic_instability/config_barotropic_inst_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)eldred_test_main.o: src/models/shallow_water/test/Eldred_test/Eldred_test_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)elsred_test_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)elsred_test_mod.o: src/models/shallow_water/test/Eldred_test/Elsred_test_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_swm_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)operator_swm_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)random_friction_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_swm_mod.o: src/models/shallow_water/config/config_swm_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)config_domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_swm_factory_mod.o: src/models/shallow_water/operator/operator_swm_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)config_swm_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)operator_swm_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)curl_factory_mod.o \
	$(DOBJ)coriolis_factory_mod.o \
	$(DOBJ)ke_factory_mod.o \
	$(DOBJ)massflux_factory_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)quadrature_factory_mod.o \
	$(DOBJ)hordiff_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)operator_adv_swm_mod.o \
	$(DOBJ)vector_advection_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_adv_swm_mod.o: src/models/shallow_water/operator/operator_adv_swm_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)key_value_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_swm_diff_mod.o: src/models/shallow_water/operator/operator_swm_diff_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_swm_mod.o: src/models/shallow_water/operator/operator_swm_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)abstract_ke_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)coriolis_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_swm_diff_factory_mod.o: src/models/shallow_water/operator/operator_swm_diff_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)config_swm_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)hordiff_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nh_model_mod.o: src/models/NH/nh_model_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)stvec_nh_mod.o \
	$(DOBJ)abstract_postprocessing_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nh_main.o: src/models/NH/NH_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)cmd_args_mod.o \
	$(DOBJ)namelist_read_mod.o \
	$(DOBJ)nh_model_mod.o \
	$(DOBJ)nh_model_config_mod.o \
	$(DOBJ)nh_model_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nh_model_config_mod.o: src/models/NH/nh_model_config_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)config_domain_mod.o \
	$(DOBJ)config_postnh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)config_nh_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nh_model_factory_mod.o: src/models/NH/nh_model_factory_mod.f90 \
	$(DOBJ)nh_model_mod.o \
	$(DOBJ)nh_model_config_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_nh_factory_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)nh_testcases_mod.o \
	$(DOBJ)postprocessing_nh_factory_mod.o \
	$(DOBJ)nh_operator_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nh_operator_factory_mod.o: src/models/NH/operators/nh_operator_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)config_mod.o \
	$(DOBJ)config_nh_operator_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)ptheta_linear_nh_oper_mod.o \
	$(DOBJ)grad_3d_factory_mod.o \
	$(DOBJ)div_3d_factory_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)interpolator_w2uv_factory_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)vertical_test_field_mod.o \
	$(DOBJ)const_n_profile_mod.o \
	$(DOBJ)advection3d_oper_mod.o \
	$(DOBJ)scalar_advection_factory_mod.o \
	$(DOBJ)solid_rotation_wind_field_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)nonlin_nh_oper_mod.o \
	$(DOBJ)mixvec_transform_factory_mod.o \
	$(DOBJ)vector_advection3d_factory_mod.o \
	$(DOBJ)coriolis_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nonlin_nh_oper_mod.o: src/models/NH/operators/nonlin_nh_oper_mod.f90 \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_grad_3d_mod.o \
	$(DOBJ)abstract_div_3d_mod.o \
	$(DOBJ)abstract_mixvec_transform_mod.o \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)abstract_scalar_advection3d_mod.o \
	$(DOBJ)abstract_vector_advection3d_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_nh_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ptheta_linear_nh_oper_mod.o: src/models/NH/operators/Ptheta_linear_nh_oper_mod.f90 \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_grad_3d_mod.o \
	$(DOBJ)abstract_div_3d_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_nh_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)advection3d_oper_mod.o: src/models/NH/operators/advection3d_oper_mod.f90 \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_nh_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_scalar_advection3d_mod.o \
	$(DOBJ)test_fieds_3d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_nh_operator_mod.o: src/models/NH/operators/config_nh_operator_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)config_advection_3d_mod.o \
	$(DOBJ)config_mixvec_transform_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_nh_factory_mod.o: src/models/NH/stvec/stvec_nh_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_nh_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_nh_mod.o: src/models/NH/stvec/stvec_nh_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)postprocessing_nh_factory_mod.o: src/models/NH/postprocessing/postprocessing_nh_factory_mod.f90 \
	$(DOBJ)abstract_postprocessing_mod.o \
	$(DOBJ)simple_postnh_mod.o \
	$(DOBJ)config_postnh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)simple_postnh_mod.o: src/models/NH/postprocessing/simple_postnh_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)abstract_postprocessing_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)stvec_nh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_postprocessing_mod.o: src/models/NH/postprocessing/abstract_postprocessing_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_postnh_mod.o: src/models/NH/postprocessing/config_postnh_mod.f90 \
	$(DOBJ)config_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nh_testcases_mod.o: src/models/NH/testcases/nh_testcases_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)gw_testcase_mod.o \
	$(DOBJ)nh_solid_rotation_testcase_mod.o \
	$(DOBJ)advection3d_testcases_mod.o \
	$(DOBJ)straka_testcase_mod.o \
	$(DOBJ)hot_bubble_testcase_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)advection3d_testcases_mod.o: src/models/NH/testcases/advection3d_testcases_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_nh_mod.o \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)straka_testcase_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)gw_testcase_mod.o: src/models/NH/testcases/GW_testcase_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_nh_mod.o \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)solid_rotation_fields_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hot_bubble_testcase_mod.o: src/models/NH/testcases/hot_bubble_testcase_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_nh_mod.o \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)straka_testcase_mod.o: src/models/NH/testcases/Straka_testcase_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_nh_mod.o \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)nh_solid_rotation_testcase_mod.o: src/models/NH/testcases/nh_solid_rotation_testcase_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_nh_mod.o \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)solid_rotation_fields_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)topology_factory_mod.o: src/topology/topology_factory_mod.f90 \
	$(DOBJ)topology_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)topology_mod.o: src/topology/topology_mod.f90 \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)cubed_sphere_topology_mod.o: src/topology/cubed_sphere_topology_mod.f90 \
	$(DOBJ)topology_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)quadrature_factory_mod.o: src/quadrature/quadrature_factory_mod.f90 \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)default_quadrature_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_quadrature_mod.o \
	$(DOBJ)sbp_operators_collection_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_quadrature_mod.o: src/quadrature/abstract_quadrature_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sbp_quadrature_mod.o: src/quadrature/sbp_quadrature_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_quadrature_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)default_quadrature_mod.o: src/quadrature/default_quadrature_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_quadrature_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)halo_ah_scalar_sync_mod.o: src/halo/halo_Ah_scalar_sync_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)halo_c_default_mod.o: src/halo/halo_C_default_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
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
	$(DOBJ)halo_ah_scalar_sync_mod.o \
	$(DOBJ)halo_c_default_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)halo_mod.o: src/halo/halo_mod.f90 \
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

$(DOBJ)adv_z_mod.o: src/differential_operators/vertical/adv_z_mod.f90 \
	$(DOBJ)abstract_adv_z_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)identity_vertical_operator_mod.o: src/differential_operators/vertical/identity_vertical_operator_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sbp_vertical_operator_mod.o: src/differential_operators/vertical/sbp_vertical_operator_mod.f90 \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_vertical_operator_mod.o: src/differential_operators/vertical/abstract_vertical_operator_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)adv_z_factory_mod.o: src/differential_operators/vertical/adv_z_factory_mod.f90 \
	$(DOBJ)abstract_adv_z_mod.o \
	$(DOBJ)adv_z_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_adv_z_mod.o: src/differential_operators/vertical/abstract_adv_z_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vertical_operator_factory_mod.o: src/differential_operators/vertical/vertical_operator_factory_mod.f90 \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)identity_vertical_operator_mod.o \
	$(DOBJ)sbp_vertical_operator_mod.o \
	$(DOBJ)sbp_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_div_mod.o: src/differential_operators/horizontal/divergence/abstract_div_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_c2_mod.o: src/differential_operators/horizontal/divergence/div_c2_mod.f90 \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_factory_mod.o: src/differential_operators/horizontal/divergence/div_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)div_c2_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)div_c_sbp42_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)div_a2_mod.o \
	$(DOBJ)div_ah2_mod.o \
	$(DOBJ)div_ah_sbp_mod.o \
	$(DOBJ)div_ch_sbp_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_ch_sbp_mod.o: src/differential_operators/horizontal/divergence/div_ch_sbp_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_ah_sbp_mod.o: src/differential_operators/horizontal/divergence/div_ah_sbp_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_c_sbp42_mod.o: src/differential_operators/horizontal/divergence/div_c_sbp42_mod.f90 \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_a2_mod.o: src/differential_operators/horizontal/divergence/div_a2_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_ah2_mod.o: src/differential_operators/horizontal/divergence/div_ah2_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_ke_mod.o: src/differential_operators/horizontal/kinetic_energy/abstract_KE_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ke_factory_mod.o: src/differential_operators/horizontal/kinetic_energy/KE_factory_mod.f90 \
	$(DOBJ)abstract_ke_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)ke_colocated_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)ke_cgrid_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ke_cgrid_mod.o: src/differential_operators/horizontal/kinetic_energy/KE_Cgrid_mod.f90 \
	$(DOBJ)abstract_ke_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ke_colocated_mod.o: src/differential_operators/horizontal/kinetic_energy/KE_colocated_mod.f90 \
	$(DOBJ)abstract_ke_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_c2_ecs_mod.o: src/differential_operators/horizontal/gradient/grad_c2_ecs_mod.f90 \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_factory_mod.o: src/differential_operators/horizontal/gradient/grad_factory_mod.f90 \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grad_c2_ecs_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)grad_c_sbp42_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)grad_a2_mod.o \
	$(DOBJ)grad_ah_sbp_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)grad_ch_sbp_mod.o \
	$(DOBJ)grad_ch_halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_grad_mod.o: src/differential_operators/horizontal/gradient/abstract_grad_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_ah_sbp_mod.o: src/differential_operators/horizontal/gradient/grad_ah_sbp_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_ch_halo_mod.o: src/differential_operators/horizontal/gradient/grad_ch_halo_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)sbp_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_c_sbp42_mod.o: src/differential_operators/horizontal/gradient/grad_c_sbp42_mod.f90 \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_a2_mod.o: src/differential_operators/horizontal/gradient/grad_a2_mod.f90 \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_ch_sbp_mod.o: src/differential_operators/horizontal/gradient/grad_ch_sbp_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vector_advection_factory_mod.o: src/differential_operators/horizontal/vec_advection/vector_advection_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)v_nabla_mod.o \
	$(DOBJ)hor_christofel_factory_mod.o \
	$(DOBJ)vector_advection_ah_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)v_nabla_sbp_factory_mod.o \
	$(DOBJ)vector_advection_c_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)abstract_v_nabla_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vector_advection_c_mod.o: src/differential_operators/horizontal/vec_advection/vector_advection_C_mod.f90 \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)abstract_hor_christofel_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_vector_advection_mod.o: src/differential_operators/horizontal/vec_advection/abstract_vector_advection_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vector_advection_ah_mod.o: src/differential_operators/horizontal/vec_advection/vector_advection_Ah_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)abstract_hor_christofel_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_perp_factory_mod.o: src/differential_operators/horizontal/grad_perp/grad_perp_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_perp_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grad_perp_c_sbp_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)halo_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_perp_c_sbp_mod.o: src/differential_operators/horizontal/grad_perp/grad_perp_c_sbp_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_perp_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_grad_perp_mod.o: src/differential_operators/horizontal/grad_perp/abstract_grad_perp_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_coriolis_mod.o: src/differential_operators/horizontal/coriolis/abstract_coriolis_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)coriolis_factory_mod.o: src/differential_operators/horizontal/coriolis/coriolis_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)coriolis_cgrid_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)interpolator_w2h_factory_mod.o \
	$(DOBJ)coriolis_cgrid_noncons_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)coriolis_colocated_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)coriolis_cgrid_noncons_mod.o: src/differential_operators/horizontal/coriolis/coriolis_Cgrid_noncons_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)abstract_co2contra_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)coriolis_cgrid_mod.o: src/differential_operators/horizontal/coriolis/coriolis_Cgrid_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)interpolator_w2h_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)coriolis_colocated_mod.o: src/differential_operators/horizontal/coriolis/coriolis_colocated_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_massflux_mod.o: src/differential_operators/horizontal/massflux/abstract_massflux_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)massflux_cgrid_mod.o: src/differential_operators/horizontal/massflux/massflux_Cgrid_mod.f90 \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)massflux_factory_mod.o: src/differential_operators/horizontal/massflux/massflux_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)massflux_colocated_mod.o \
	$(DOBJ)massflux_cgrid_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)sbp_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)massflux_colocated_mod.o: src/differential_operators/horizontal/massflux/massflux_colocated_mod.f90 \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)co2contra_ch_mod.o: src/differential_operators/horizontal/co2contra/co2contra_Ch_mod.f90 \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)co2contra_colocated_mod.o: src/differential_operators/horizontal/co2contra/co2contra_colocated_mod.f90 \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)co2contra_factory_mod.o: src/differential_operators/horizontal/co2contra/co2contra_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)co2contra_colocated_mod.o \
	$(DOBJ)co2contra_cgrid_mod.o \
	$(DOBJ)co2contra_ch_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)co2contra_cgrid_mod.o: src/differential_operators/horizontal/co2contra/co2contra_Cgrid_mod.f90 \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_co2contra_mod.o: src/differential_operators/horizontal/co2contra/abstract_co2contra_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hor_christofel_shallow_mod.o: src/differential_operators/horizontal/Christofel/hor_Christofel_shallow_mod.f90 \
	$(DOBJ)abstract_hor_christofel_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hor_christofel_factory_mod.o: src/differential_operators/horizontal/Christofel/hor_Christofel_factory_mod.f90 \
	$(DOBJ)abstract_hor_christofel_mod.o \
	$(DOBJ)hor_christofel_shallow_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_hor_christofel_mod.o: src/differential_operators/horizontal/Christofel/abstract_hor_Christofel_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)laplace_ah_sbp63_narrow_mod.o: src/differential_operators/horizontal/laplace/laplace_Ah_sbp63_narrow_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_laplace_mod.o: src/differential_operators/horizontal/laplace/abstract_laplace_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)divgrad_laplace_mod.o: src/differential_operators/horizontal/laplace/divgrad_laplace_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_co2contra_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)laplace_ah_sbp42_narrow_mod.o: src/differential_operators/horizontal/laplace/laplace_Ah_sbp42_narrow_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)laplace_factory_mod.o: src/differential_operators/horizontal/laplace/laplace_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)divgrad_laplace_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)laplace_ch_halo_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)laplace_ah_sbp21_narrow_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)laplace_ah_sbp42_narrow_mod.o \
	$(DOBJ)laplace_ah_sbp63_narrow_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)laplace_ah_sbp21_narrow_mod.o: src/differential_operators/horizontal/laplace/laplace_Ah_sbp21_narrow_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)laplace_ch_halo_mod.o: src/differential_operators/horizontal/laplace/laplace_ch_halo_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)curl_c_sbp_mod.o: src/differential_operators/horizontal/curl/curl_c_sbp_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_curl_mod.o: src/differential_operators/horizontal/curl/abstract_curl_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)curl_factory_mod.o: src/differential_operators/horizontal/curl/curl_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)curl_div_based_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)curl_c_sbp_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)halo_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)curl_div_based_mod.o: src/differential_operators/horizontal/curl/curl_div_based_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)v_nabla_factory_mod.o: src/differential_operators/horizontal/v_dot_nabla/v_nabla_factory_mod.f90 \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)v_nabla_mod.o \
	$(DOBJ)v_nabla_sbp_factory_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_v_nabla_mod.o: src/differential_operators/horizontal/v_dot_nabla/abstract_v_nabla_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)v_nabla_ah_sbp_mod.o: src/differential_operators/horizontal/v_dot_nabla/v_nabla_Ah_sbp_mod.f90 \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)v_nabla_sbp_factory_mod.o: src/differential_operators/horizontal/v_dot_nabla/v_nabla_sbp_factory_mod.f90 \
	$(DOBJ)v_nabla_ah_sbp_mod.o \
	$(DOBJ)sbp_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)v_nabla_mod.o: src/differential_operators/horizontal/v_dot_nabla/v_nabla_mod.f90 \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolator_w2h_factory_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_w2h_factory_mod.f90 \
	$(DOBJ)interpolator_w2h_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolator_uv2p_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_uv2p_mod.f90 \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolator_q2uv_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_q2uv_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolator_w2h_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_w2h_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_interpolators2d_mod.o: src/differential_operators/horizontal/interpolator_2d/abstract_interpolators2d_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolator_p2uv_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_p2uv_mod.f90 \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolator_uv2q_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_uv2q_mod.f90 \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolator2d_factory_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator2d_factory_mod.f90 \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)interpolator_p2uv_mod.o \
	$(DOBJ)interpolator_uv2p_mod.o \
	$(DOBJ)interpolator_q2uv_mod.o \
	$(DOBJ)interpolator_uv2q_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hordiff_factory_mod.o: src/differential_operators/horizontal/hordiff/hordiff_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)hordiff_cgrid_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)curl_factory_mod.o \
	$(DOBJ)grad_perp_factory_mod.o \
	$(DOBJ)hordiff_scalar_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)laplace_factory_mod.o \
	$(DOBJ)hordiff_colocated_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hordiff_no_metric_mod.o: src/differential_operators/horizontal/hordiff/hordiff_no_metric_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)vec_math_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hordiff_cgrid_mod.o: src/differential_operators/horizontal/hordiff/hordiff_Cgrid_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_grad_perp_mod.o \
	$(DOBJ)abstract_co2contra_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hordiff_colocated_mod.o: src/differential_operators/horizontal/hordiff/hordiff_colocated_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)vec_math_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hordiff_scalar_mod.o: src/differential_operators/horizontal/hordiff/hordiff_scalar_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_laplace_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_hordiff_mod.o: src/differential_operators/horizontal/hordiff/abstract_hordiff_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_div_3d_mod.o: src/differential_operators/3d/divergence/abstract_div_3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_3d_hor_vert_mod.o: src/differential_operators/3d/divergence/div_3d_hor_vert_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_3d_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)vec_math_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_3d_factory_mod.o: src/differential_operators/3d/divergence/div_3d_factory_mod.f90 \
	$(DOBJ)abstract_div_3d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)div_3d_hor_vert_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_3d_factory_mod.o: src/differential_operators/3d/gradient/grad_3d_factory_mod.f90 \
	$(DOBJ)abstract_grad_3d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grad_3d_hor_vert_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_3d_hor_vert_mod.o: src/differential_operators/3d/gradient/grad_3d_hor_vert_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_3d_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_grad_3d_mod.o: src/differential_operators/3d/gradient/abstract_grad_3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vector_advection3d_factory_mod.o: src/differential_operators/3d/advection/vector_advection3d_factory_mod.f90 \
	$(DOBJ)abstract_vector_advection3d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_mod.o \
	$(DOBJ)shallow_atm_vecadv_mod.o \
	$(DOBJ)vector_advection_factory_mod.o \
	$(DOBJ)adv_z_factory_mod.o \
	$(DOBJ)interpolator_w2uv_factory_mod.o \
	$(DOBJ)scalar_advection_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)config_advection_3d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_advection_3d_mod.o: src/differential_operators/3d/advection/config_advection_3d_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)advection_p3d_mod.o: src/differential_operators/3d/advection/advection_p3d_mod.f90 \
	$(DOBJ)abstract_scalar_advection3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)abstract_adv_z_mod.o \
	$(DOBJ)exchange_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_scalar_advection3d_mod.o: src/differential_operators/3d/advection/abstract_scalar_advection3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_vector_advection3d_mod.o: src/differential_operators/3d/advection/abstract_vector_advection3d_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)advection_w3d_mod.o: src/differential_operators/3d/advection/advection_w3d_mod.f90 \
	$(DOBJ)abstract_scalar_advection3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)abstract_adv_z_mod.o \
	$(DOBJ)exchange_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)shallow_atm_vecadv_mod.o: src/differential_operators/3d/advection/shallow_atm_vecadv_mod.f90 \
	$(DOBJ)abstract_vector_advection3d_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)abstract_scalar_advection3d_mod.o \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_adv_z_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)scalar_advection_factory_mod.o: src/differential_operators/3d/advection/scalar_advection_factory_mod.f90 \
	$(DOBJ)abstract_scalar_advection3d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_mod.o \
	$(DOBJ)v_nabla_factory_mod.o \
	$(DOBJ)adv_z_factory_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)interpolator_uv2w_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)advection_p3d_mod.o \
	$(DOBJ)config_advection_3d_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)advection_w3d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolators_w2uv_mod.o: src/differential_operators/3d/interpolators/interpolators_w2uv_mod.f90 \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolators_uv2w_mod.o: src/differential_operators/3d/interpolators/interpolators_uv2w_mod.f90 \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolator_w2uv_factory_mod.o: src/differential_operators/3d/interpolators/interpolator_w2uv_factory_mod.f90 \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)interpolators_w2uv_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_interpolators3d_mod.o: src/differential_operators/3d/interpolators/abstract_interpolators3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolator_uv2w_factory_mod.o: src/differential_operators/3d/interpolators/interpolator_uv2w_factory_mod.f90 \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)interpolators_uv2w_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mixvec_transform_staggered_mod.o: src/differential_operators/3d/mixvec_transform/mixvec_transform_staggered_mod.f90 \
	$(DOBJ)abstract_mixvec_transform_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mixvec_transform_hor_colocated_mod.o: src/differential_operators/3d/mixvec_transform/mixvec_transform_hor_colocated_mod.f90 \
	$(DOBJ)abstract_mixvec_transform_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mixvec_transform_colocated_mod.o: src/differential_operators/3d/mixvec_transform/mixvec_transform_colocated_mod.f90 \
	$(DOBJ)abstract_mixvec_transform_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_mixvec_transform_mod.o: src/differential_operators/3d/mixvec_transform/config_mixvec_transform_mod.f90 \
	$(DOBJ)config_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_mixvec_transform_mod.o: src/differential_operators/3d/mixvec_transform/abstract_mixvec_transform_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mixvec_transform_factory_mod.o: src/differential_operators/3d/mixvec_transform/mixvec_transform_factory_mod.f90 \
	$(DOBJ)abstract_mixvec_transform_mod.o \
	$(DOBJ)mixvec_transform_colocated_mod.o \
	$(DOBJ)mixvec_transform_hor_colocated_mod.o \
	$(DOBJ)mixvec_transform_staggered_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)config_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)config_mixvec_transform_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_co2contra_3d_mod.o: src/differential_operators/3d/co2contra/abstract_co2contra_3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)co2contra_3d_cgrid_mod.o: src/differential_operators/3d/co2contra/co2contra_3d_Cgrid_mod.f90 \
	$(DOBJ)abstract_co2contra_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)co2contra_3d_colocated_mod.o: src/differential_operators/3d/co2contra/co2contra_3d_colocated_mod.f90 \
	$(DOBJ)abstract_co2contra_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)co2contra_3d_h_colocated_mod.o: src/differential_operators/3d/co2contra/co2contra_3d_h_colocated_mod.f90 \
	$(DOBJ)abstract_co2contra_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)co2contra_3d_factory_mod.o: src/differential_operators/3d/co2contra/co2contra_3d_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_co2contra_3d_mod.o \
	$(DOBJ)co2contra_3d_colocated_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)co2contra_3d_cgrid_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)co2contra_3d_h_colocated_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sbp_factory_mod.o: src/differential_operators/sbp_operators/sbp_factory_mod.f90 \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_operators_collection_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sbp_norm_mod.o: src/differential_operators/sbp_operators/sbp_norm_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sbp_operator_mod.o: src/differential_operators/sbp_operators/sbp_operator_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sbp_operators_collection_mod.o: src/differential_operators/sbp_operators/sbp_operators_collection_mod.f90
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

$(DOBJ)test_vertical_operators_main.o: src/test/test_vertical_operators/test_vertical_operators_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_vertical_operators_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_vertical_operators_mod.o: src/test/test_vertical_operators/test_vertical_operators_mod.f90 \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)config_domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)vertical_test_field_mod.o \
	$(DOBJ)const_n_profile_mod.o \
	$(DOBJ)abstract_vertical_profile_mod.o \
	$(DOBJ)vertical_div_test_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_cmd_line.o: src/test/test_cmd_line/test_cmd_line.f90 \
	$(DOBJ)cmd_args_test_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)cmd_args_test_mod.o: src/test/test_cmd_line/cmd_args_test_mod.f90 \
	$(DOBJ)cmd_args_mod.o
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

$(DOBJ)test_mod.o: src/test/test_exch/test_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)exchange_halo_c_mod.o \
	$(DOBJ)exchange_halo_ch_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_exch_main.o: src/test/test_exch/test_exch_main.f90 \
	$(DOBJ)test_mod.o \
	$(DOBJ)parcomm_mod.o
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

$(DOBJ)test_diffops_conv.o: src/test/test_diff_ops/test_diffops_conv.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_diffops_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_diffops_3d_mod.o: src/test/test_diff_ops/test_diffops_3d_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)config_mod.o \
	$(DOBJ)config_domain_mod.o \
	$(DOBJ)config_orography_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)grad_3d_factory_mod.o \
	$(DOBJ)abstract_grad_3d_mod.o \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grad3d_test_field_mod.o \
	$(DOBJ)div_3d_factory_mod.o \
	$(DOBJ)abstract_div_3d_mod.o \
	$(DOBJ)div3d_test_field_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)abstract_co2contra_3d_mod.o \
	$(DOBJ)co2contra_3d_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)interpolator_w2uv_factory_mod.o \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)interpolator_uv2w_factory_mod.o \
	$(DOBJ)mixvec_transform_factory_mod.o \
	$(DOBJ)abstract_mixvec_transform_mod.o \
	$(DOBJ)config_mixvec_transform_mod.o \
	$(DOBJ)abstract_scalar_advection3d_mod.o \
	$(DOBJ)scalar_advection_factory_mod.o \
	$(DOBJ)config_advection_3d_mod.o \
	$(DOBJ)solid_rotation_wind_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_diffops.o: src/test/test_diff_ops/test_diffops.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_diffops_mod.o \
	$(DOBJ)key_value_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_diffops_all.o: src/test/test_diff_ops/test_diffops_all.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_diffops_mod.o \
	$(DOBJ)key_value_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_diffops_mod.o: src/test/test_diff_ops/test_diffops_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)laplace_factory_mod.o \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)grad_perp_factory_mod.o \
	$(DOBJ)abstract_grad_perp_mod.o \
	$(DOBJ)curl_factory_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)ke_factory_mod.o \
	$(DOBJ)abstract_ke_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)coriolis_factory_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)vector_advection_factory_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)quadrature_factory_mod.o \
	$(DOBJ)abstract_quadrature_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_diffops_3d.o: src/test/test_diff_ops/test_diffops_3d.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_diffops_3d_mod.o \
	$(DOBJ)key_value_mod.o
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
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_generic_halo_mod.o
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

$(DOBJ)test_generic_halo_mod.o: src/test/test_halo/test_generic_halo_mod.f90 \
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

$(DOBJ)test_regrid_mod.o: src/test/test_regrid/test_regrid_mod.f90 \
	$(DOBJ)abstract_regridder_mod.o \
	$(DOBJ)regrid_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)latlon_functions_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_regrid.o: src/test/test_regrid/test_regrid.f90 \
	$(DOBJ)test_regrid_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_latlon_output_mod.o: src/test/test_latlon_output/test_latlon_output_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_fields_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_latlon_output.o: src/test/test_latlon_output/test_latlon_output.f90 \
	$(DOBJ)test_latlon_output_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_vertical_profiles_mod.o: src/test/test_vertical_profiles/test_vertical_profiles_mod.f90 \
	$(DOBJ)abstract_vertical_profile_mod.o \
	$(DOBJ)const_n_profile_mod.o \
	$(DOBJ)isotermal_profile_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_vertical_profiles_main.o: src/test/test_vertical_profiles/test_vertical_profiles_main.f90 \
	$(DOBJ)test_vertical_profiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_metric_class_mod.o: src/test/test_metric/test_metric_class_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)metric_factory_mod.o \
	$(DOBJ)config_metric_mod.o \
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
	$(DOBJ)config_domain_mod.o \
	$(DOBJ)config_orography_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)tile_mod.o
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
	$(DOBJ)config_domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)config_orography_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_advection_main.o: src/test/test_advection/test_advection_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_solid_rotation_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_advection_factory_mod.o: src/test/test_advection/operator_advection_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_advection_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)massflux_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_solid_rotation_mod.o: src/test/test_advection/test_solid_rotation_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_advection_mod.o \
	$(DOBJ)stvec_advection_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_advection_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)vec_math_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_advection_mod.o: src/test/test_advection/stvec_advection_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_advection_factory_mod.o: src/test/test_advection/stvec_advection_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_advection_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_advection_mod.o: src/test/test_advection/operator_advection_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)stvec_advection_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_hordiff_main.o: src/test/test_hordiff/test_hordiff_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_hordiff_scalar_mod.o \
	$(DOBJ)test_hordiff_vector_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_hordiff_scalar_mod.o: src/test/test_hordiff/test_hordiff_scalar_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)hordiff_factory_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)quadrature_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_hordiff_vector_mod.o: src/test/test_hordiff/test_hordiff_vector_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)hordiff_factory_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)halo_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_vertical_transform_mod.o: src/test/test_vertical_transform/test_vertical_transform_mod.f90 \
	$(DOBJ)abstract_vertical_transform_mod.o \
	$(DOBJ)vertical_transform_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_vertical_transform_main.o: src/test/test_vertical_transform/test_vertical_transform_main.f90 \
	$(DOBJ)test_vertical_transform_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_quadrature_mod.o: src/test/test_quadrature/test_quadrature_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)quadrature_factory_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_quadrature_main.o: src/test/test_quadrature/test_quadrature_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_quadrature_mod.o
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

$(DOBJ)test_paneled_output.o: src/test/test_paneled_output/test_paneled_output.f90 \
	$(DOBJ)test_paneled_output_mod.o \
	$(DOBJ)parcomm_mod.o
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

$(DOBJ)ecs_metric_mod.o: src/equiang_cs/ecs_metric_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)orography_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_factory_mod.o: src/equiang_cs/ecs_halo_factory_mod.f90 \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)ecs_halo_xy_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_halo_xy_mod.o: src/equiang_cs/ecs_halo_xy_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)parcomm_mod.o \
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

$(DOBJ)ecs_ah_vec_sync_mod.o: src/equiang_cs/ecs_Ah_vec_sync_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)ecs_metric_factory_mod.o: src/equiang_cs/ecs_metric_factory_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)ecs_metric_mod.o \
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

$(DOBJ)ecs_halo_vec_c_mod.o: src/equiang_cs/ecs_halo_vec_c_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
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

$(DOBJ)parcomm_factory_mod.o: src/parallel/parcomm_factory_mod.f90 \
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

$(DOBJ)exchange_halo_ch_mod.o: src/parallel/exchange_halo_Ch_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)buffer_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)partition_factory_mod.o: src/parallel/partition_factory_mod.f90 \
	$(DOBJ)partition_mod.o
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

$(DOBJ)exchange_abstract_mod.o: src/parallel/exchange_abstract_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)buffer_mod.o: src/parallel/buffer_mod.f90 \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)partition_mod.o: src/parallel/partition_mod.f90 \
	$(DOBJ)tiles_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parcomm_mod.o: src/parallel/parcomm_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_factory_mod.o: src/parallel/exchange_factory_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)tiles_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)exchange_halo_ch_mod.o \
	$(DOBJ)exchange_halo_c_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_gather_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)exchange_gather_mod.o: src/parallel/exchange_gather_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)buffer_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)domain_mod.o: src/domain/domain_mod.f90 \
	$(DOBJ)topology_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)orography_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)domain_factory_mod.o: src/domain/domain_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)topology_factory_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)metric_factory_mod.o \
	$(DOBJ)orography_mod.o \
	$(DOBJ)orography_factory_mod.o \
	$(DOBJ)config_domain_mod.o \
	$(DOBJ)mesh_factory_mod.o \
	$(DOBJ)parcomm_factory_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_domain_mod.o: src/domain/config_domain_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)config_metric_mod.o \
	$(DOBJ)config_orography_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)config_orography_mod.o: src/domain/orography/config_orography_mod.f90 \
	$(DOBJ)config_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)orography_factory_mod.o: src/domain/orography/orography_factory_mod.f90 \
	$(DOBJ)orography_mod.o \
	$(DOBJ)config_mod.o \
	$(DOBJ)config_orography_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)orography_test_field_mod.o \
	$(DOBJ)schar_orography_field_mod.o \
	$(DOBJ)test_fieds_3d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)orography_mod.o: src/domain/orography/orography_mod.f90 \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)master_paneled_outputer_mod.o: src/outputer/master_paneled_outputer_mod.f90 \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mpi_paneled_outputer_mod.o: src/outputer/mpi_paneled_outputer_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)latlon_outputer_mod.o: src/outputer/latlon_outputer_mod.f90 \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tiles_mod.o \
	$(DOBJ)abstract_regridder_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)outputer_factory_mod.o: src/outputer/outputer_factory_mod.f90 \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)master_paneled_outputer_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)latlon_outputer_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)parcomm_factory_mod.o \
	$(DOBJ)regrid_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)outputer_abstract_mod.o: src/outputer/outputer_abstract_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)regrid_factory_mod.o: src/regridders/regrid_factory_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)abstract_regridder_mod.o \
	$(DOBJ)regrid_to_latlon_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_regridder_mod.o: src/regridders/abstract_regridder_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)regrid_to_latlon_mod.o: src/regridders/regrid_to_latlon_mod.f90 \
	$(DOBJ)abstract_regridder_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)latlon_functions_mod.o: src/test_fields/latlon_functions_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_fields_mod.o: src/test_fields/test_fields_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o \
	$(DOBJ)latlon_functions_mod.o \
	$(DOBJ)barotropic_instability_u_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vertical_test_field_mod.o: src/test_fields/test_fields_3d/vertical_test_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)abstract_vertical_profile_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)schar_orography_field_mod.o: src/test_fields/test_fields_3d/Schar_orography_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)orography_test_field_mod.o: src/test_fields/test_fields_3d/orography_test_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_fieds_3d_mod.o: src/test_fields/test_fields_3d/test_fieds_3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div3d_test_field_mod.o: src/test_fields/test_fields_3d/div3d_test_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vertical_div_test_field_mod.o: src/test_fields/test_fields_3d/vertical_div_test_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad3d_test_field_mod.o: src/test_fields/test_fields_3d/grad3d_test_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solid_rotation_fields_factory_mod.o: src/test_fields/test_fields_3d/solid_rotation/solid_rotation_fields_factory_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)solid_rotation_therm_mod.o \
	$(DOBJ)solid_rotation_wind_field_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solid_rotation_therm_mod.o: src/test_fields/test_fields_3d/solid_rotation/solid_rotation_therm_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solid_rotation_wind_field_mod.o: src/test_fields/test_fields_3d/solid_rotation/solid_rotation_wind_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)isotermal_profile_mod.o: src/test_fields/vertical_thermodynamic_profiles/isotermal_profile_mod.f90 \
	$(DOBJ)abstract_vertical_profile_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)abstract_vertical_profile_mod.o: src/test_fields/vertical_thermodynamic_profiles/abstract_vertical_profile_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)const_n_profile_mod.o: src/test_fields/vertical_thermodynamic_profiles/const_N_profile_mod.f90 \
	$(DOBJ)abstract_vertical_profile_mod.o \
	$(DOBJ)const_mod.o
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

$(DOBJ)rk4_mod.o: src/time_schemes/rk4_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o
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
