program main

use test_metric_mod, only : test_metric
use mpi

call MPI_init(ierr)

call test_metric

call mpi_barrier(mpi_comm_world, ierr)
call mpi_finalize(ierr)

end
