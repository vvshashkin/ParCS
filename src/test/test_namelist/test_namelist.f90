program main

use test_namelist_mod, only : test_namelist
use mpi

call MPI_init(ierr)

call test_namelist

call mpi_barrier(mpi_comm_world, ierr)
call mpi_finalize(ierr)

end
