program test_halo_ecs

use test_halo_mod, only : test_ecs_halo
use mpi

call MPI_init(ierr)

call test_ecs_halo

call mpi_barrier(mpi_comm_world, ierr)
call mpi_finalize(ierr)

end
