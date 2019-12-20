program output_test

    use test_output_mod, only : test_output
    use mpi

    call MPI_init(ierr)

    call test_output()

    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize(ierr)

end program output_test
