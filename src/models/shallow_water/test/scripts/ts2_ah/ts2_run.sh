#!/bin/bash

EXE=$1/TS2_MAIN
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
    uv_diff_coeff     =  0.03,\n
    hordiff_uv_name   = 'hordiff_vec_xyz_Ah',\n
    h_diff_coeff      =  0.01,\n
    hordiff_h_name    = 'hordiff_scalar_Ah',\n
    dt=%%%dt,\n
    tau_write = 86400.0,\n
    tau_diagnostics = 3600.0\n
    simulation_time_days  = 10.0,\n
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

run_ts2(){
	gen_namelist $1 $2 $3  > namelist_swm
	mpirun -n $Nprocs $EXE &> swm_N$1_dt$2_Ah$3.out
    grep "l2_h" swm_N$1_dt$2_Ah$3.out > errors_N$1_dt$2_Ah$3.txt
    mv h.dat h_N$1_Ah$3.dat
}
run_ts2 020 800 21
run_ts2 040 400 21
run_ts2 080 200 21
run_ts2 160 100 21
run_ts2 020 800 42
run_ts2 040 400 42
run_ts2 080 200 42
run_ts2 160 100 42
run_ts2 020 800 43
run_ts2 040 400 43
run_ts2 080 200 43
run_ts2 160 100 43
run_ts2 020 800 63
run_ts2 040 400 63
run_ts2 080 200 63
run_ts2 160 100 63

