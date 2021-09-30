program swm_test

use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment, parcomm_global

use test_swm_mod, only : test_swm

call init_global_parallel_enviroment()


call test_swm()


call deinit_global_parallel_enviroment()

end program swm_test
