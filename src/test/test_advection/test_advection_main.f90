program advection_test

use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment, parcomm_global

use test_advection_mod, only : test_advection

call init_global_parallel_enviroment()


call test_advection()


call deinit_global_parallel_enviroment()

end program advection_test
