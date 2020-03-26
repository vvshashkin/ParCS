program NHlin_main

    use mpi
    use NHlin_mod, only : init_NHlin_model, run_NHlin_model
    implicit none

    integer(kind=4) ierr, myid
    real(kind=8)    t1, t2, t3

    call mpi_init(ierr)
    call MPI_comm_rank(mpi_comm_world , myid, ierr)

    t1 = mpi_wtime()
    call init_NHlin_model()
    t2 = mpi_wtime()
    call run_NHlin_model()
    t3 = mpi_wtime()

    call mpi_finalize(ierr)
    if(myid == 0) &
        print *, "Time (with/without init)(s):", t3-t1, t3-t2

end program NHlin_main
