program mesh_test

    use test_mesh_mod, only : test_mesh
    use parcomm_mod,   only : init_global_parallel_enviroment, &
                              deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_mesh()

    call deinit_global_parallel_enviroment()

end program mesh_test
