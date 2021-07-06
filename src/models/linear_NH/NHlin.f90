program NHlin_main

    use parcomm_mod, only : init_global_parallel_enviroment, &
                            deinit_global_parallel_enviroment
    use NHlin_mod,  only : init_NHlin_model, run_NHlin_model
    implicit none

    real(kind=8)    t1, t2, t3

    call init_global_parallel_enviroment()

    t1 = mpi_wtime()
    call init_NHlin_model()
    t2 = mpi_wtime()
    call run_NHlin_model()
    t3 = mpi_wtime()

    call deinit_global_parallel_enviroment()

    print *, "Time (with/without init)(s):", t3-t1, t3-t2

end program NHlin_main
