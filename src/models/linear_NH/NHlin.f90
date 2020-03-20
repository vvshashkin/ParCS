program NHlin_main

    use mpi
    use NHlin_mod, only : init_NHlin_model, run_NHlin_model
    implicit none

    integer(kind=4) ierr

    call mpi_init(ierr)

    call init_NHlin_model()
    call run_NHlin_model()

    call mpi_finalize(ierr)

end program NHlin_main
