program main

use test_mod, only : test_A_halo_exchange, test_halo_vec_C_exchange
use mpi

call MPI_init(ierr)

call test_A_halo_exchange()
call test_halo_vec_C_exchange()

! call test_gather_exchange()

call mpi_barrier(mpi_comm_world, ierr)
call mpi_finalize(ierr)

end
