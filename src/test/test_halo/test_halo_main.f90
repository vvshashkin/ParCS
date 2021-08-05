program test_halo_main

    use test_ecs_halo_mod,   only : test_ecs_halo
    use test_ecs_halo_c_mod, only : test_ecs_cvec_halo
    use test_halo_mod,       only : test_halo
    use parcomm_mod,         only : init_global_parallel_enviroment, &
                                    deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    !call test_halo()
    !call test_ecs_halo
    call test_ecs_cvec_halo

    call deinit_global_parallel_enviroment()

end
