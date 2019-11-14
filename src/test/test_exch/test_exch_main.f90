program main

use test_mod, only : test_exchange
use mpi

call MPI_init(ierr)

call test_exchange()

call mpi_barrier(mpi_comm_world, ierr)
call mpi_finalize(ierr)

end
