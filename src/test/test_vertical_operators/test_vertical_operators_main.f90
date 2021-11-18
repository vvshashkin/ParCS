program test_vertical_operators_main

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment, parcomm_global
use test_vertical_operators_mod, only : test_vertical_gradient_operator, test_vertical_div_operator

implicit none

call init_global_parallel_enviroment()

call test_vertical_gradient_operator(20,"eta_diff_p2w_sbp21","CharneyPhilips")
call test_vertical_gradient_operator(20,"eta_diff_p2w_sbp42","CharneyPhilips")
call test_vertical_gradient_operator(20,"eta_diff_sbp21","None")
call test_vertical_gradient_operator(20,"eta_diff_sbp42","None")
call test_vertical_gradient_operator(20,"eta_diff_sbp63","None")

call test_vertical_div_operator(20,"eta_diff_w2p_sbp21","CharneyPhilips")
call test_vertical_div_operator(20,"eta_diff_w2p_sbp42","CharneyPhilips")

call deinit_global_parallel_enviroment()

end program test_vertical_operators_main
