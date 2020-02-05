program test_cmd_line

    use mpi
    use cmd_args_test_mod, only : test_cmd_args

    implicit none

    integer(kind=4) ierr

    call MPI_init(ierr)

    call test_cmd_args()

    call MPI_finalize(ierr)

end program test_cmd_line
