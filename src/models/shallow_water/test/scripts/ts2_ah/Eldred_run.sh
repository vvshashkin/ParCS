#!/bin/bash

EXE=$1/ELDRED_TEST_MAIN
Nprocs=$2

NAMELIST_TEMPLATE="
&domain\n
    N = %%%N,\n
    Nz = 1,\n
    topology_type    = 'cube',\n
    staggering_type  = 'Ah',\n
    metric_type      = 'ecs'\n
/\n
&metric\n
/\n
&shallow_water_model\n
    swm_op_type       = 'vector_invariant'\n
    v_components_type = 'covariant'\n
    div_op_name       = '%%%divergence',\n
    grad_op_name      = '%%%gradient',\n
    coriolis_op_name  = 'coriolis_colocated',\n
    curl_op_name      = '%%%curl',\n
    KE_op_name        = 'KE_colocated',\n
    co2contra_op_name = 'co2contra_colocated',\n
    massflux_op_name  = 'massflux_colocated',\n
    quadrature_name   = '%%%quadrature',\n
    diff_time_scheme  = 'explicit_Eul1'\n
    uv_diff_coeff     =  0.04,\n
    hordiff_uv_name   = 'hordiff_vec_xyz_Ah',\n
    h_diff_coeff      =  0.01,\n
    hordiff_h_name    = 'hordiff_scalar_Ah',\n
    dt=%%%dt,\n
    tau_write = 86400.0,\n
    tau_diagnostics = 3600.0\n
    simulation_time_days  = 1400.0,\n
    simulation_time_hours = 0.0,\n
    simulation_time_min   = 0.0,\n
    simulation_time_sec   = 0.0,\n
/"

gen_namelist(){
    N=$1
    DT=$2
    if [[ $3 -eq "21" ]]; then
        div="divergence_ah2"
        grad="gradient_ah21_sbp_ecs"
        curl="curl_"$div
        quad='SBP_Ah21_quadrature'
    elif  [[ $3 -eq "42" ]]; then
        div="divergence_ah42_sbp"
        grad="gradient_ah42_sbp_ecs"
        curl="curl_"$div
        quad='SBP_Ah42_quadrature'
    elif  [[ $3 -eq "43" ]]; then
        div="divergence_ah43_sbp"
        grad="gradient_ah43_sbp_ecs"
        curl="curl_"$div
        quad='SBP_Ah63_quadrature'
    elif  [[ $3 -eq "63" ]]; then
        div="divergence_ah63_sbp"
        grad="gradient_ah63_sbp_ecs"
        curl="curl_"$div
        quad='SBP_Ah63_quadrature'
	fi
	echo -e $NAMELIST_TEMPLATE |
             sed "s/%%%divergence/$div/" |
             sed "s/%%%gradient/$grad/"  |
             sed "s/%%%curl/$curl/"      |
             sed "s/%%%quadrature/$quad/" |
             sed "s/%%%N/$N/" |
             sed "s/%%%dt/$DT/"
}

run_Eldred(){
	gen_namelist $1 $2 $3  > namelist_swm
	mpirun -n $Nprocs $EXE &> swm_N$1_dt$2_Ah$3.out
    mv h.dat h_N$1_Ah$3.dat
    mv u.dat u_N$1_Ah$3.dat
    mv v.dat v_N$1_Ah$3.dat
    mv div.dat div_N$1_Ah$3.dat
    mv curl.dat curl_N$1_Ah$3.dat
}
run_Eldred 096 300 21
run_Eldred 096 300 42
run_Eldred 096 300 43
run_Eldred 096 300 63
