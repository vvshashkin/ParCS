program swlin_main

    use mpi
    use swlin_mod, only : init_swlin_model, run_swlin_model
    implicit none

    integer(kind=4) ierr

    call mpi_init(ierr)

    call init_swlin_model()
    call run_swlin_model()

    call mpi_finalize(ierr)

end program swlin_main
