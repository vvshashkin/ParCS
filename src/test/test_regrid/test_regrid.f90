program main

    use test_regrid_mod, only : test_regrid
    use parcomm_mod,     only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_regrid('A')
    call test_regrid('Ah')

    call deinit_global_parallel_enviroment()

end
