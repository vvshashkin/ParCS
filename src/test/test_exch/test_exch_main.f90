program main

use test_mod, only : test_cross_halo_exchange, test_full_halo_exchange
use mpi

call MPI_init(ierr)

call test_cross_halo_exchange()
call test_full_halo_exchange()

call mpi_barrier(mpi_comm_world, ierr)
call mpi_finalize(ierr)

end
