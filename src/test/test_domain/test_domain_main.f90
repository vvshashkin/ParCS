program domain_test

use test_domain_mod, only : test_domain
use mpi

implicit none

integer(kind=4) :: ierr

call MPI_init(ierr)

call test_domain()

call mpi_barrier(mpi_comm_world, ierr)
call mpi_finalize(ierr)

end program domain_test
