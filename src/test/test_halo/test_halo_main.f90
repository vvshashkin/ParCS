program test_halo_main

use test_ecs_halo_mod, only : test_ecs_halo
use test_halo_mod, only : test_halo

use mpi

call MPI_init(ierr)

!call test_halo()
call test_ecs_halo

call mpi_finalize(ierr)

end
