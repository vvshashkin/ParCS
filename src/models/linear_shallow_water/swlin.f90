program swlin_main

    use swlin_mod,   only : init_swlin_model, run_swlin_model
    use parcomm_mod, only : init_global_parallel_enviroment, &
                            deinit_global_parallel_enviroment
    implicit none

    call init_global_parallel_enviroment()

    call init_swlin_model()
    call run_swlin_model()

    call deinit_global_parallel_enviroment()

end program swlin_main
