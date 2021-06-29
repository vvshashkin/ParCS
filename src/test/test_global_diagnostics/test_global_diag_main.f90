program test_gl_diag_ecs

    use test_gl_diag_mod, only : test_gl_diag
    use parcomm_mod,      only : init_global_parallel_enviroment, &
                                 deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_gl_diag()

    call deinit_global_parallel_enviroment()

end
