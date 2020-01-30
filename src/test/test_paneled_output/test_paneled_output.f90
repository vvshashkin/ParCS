program paneled_output_test

    use test_paneled_output_mod, only : test_paneled_output
    use mpi

    call MPI_init(ierr)

    call test_paneled_output()

    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize(ierr)

end program paneled_output_test