program test_gl_diag_ecs

use test_gl_diag_mod, only : test_gl_diag
use mpi

call MPI_init(ierr)

call test_gl_diag()

call mpi_barrier(mpi_comm_world, ierr)
call mpi_finalize(ierr)

end
