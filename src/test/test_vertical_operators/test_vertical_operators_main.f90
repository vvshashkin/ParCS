program test_vertical_operators_main

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment, parcomm_global
use test_vertical_operators_mod, only : test_vertical_gradient_operator

implicit none

call init_global_parallel_enviroment()

call test_vertical_gradient_operator()

call deinit_global_parallel_enviroment()

end program test_vertical_operators_main
