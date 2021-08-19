program main

use test_mod,    only : test_A_halo_exchange, test_D_halo_exchange, test_halo_vec_C_exchange
use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_A_halo_exchange()
    call test_D_halo_exchange()
    call test_halo_vec_C_exchange()

    call deinit_global_parallel_enviroment()

end
