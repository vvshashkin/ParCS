program swlin_main

    use mpi
    use swlin_mod, only : swlin_model_main
    implicit none

    integer(kind=4) ierr

    call mpi_init(ierr)

    print *, "linear shallow water model"

    call swlin_model_main()

    call mpi_finalize(ierr)

end program swlin_main
