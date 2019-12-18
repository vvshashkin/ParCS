program mesh_test

    use test_mesh_mod, only : test_mesh
    use mpi

    call MPI_init(ierr)

    call test_mesh()


    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize(ierr)



end program mesh_test
