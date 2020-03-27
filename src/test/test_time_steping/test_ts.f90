program test_ts
use test_rk4_mod, only: test_rk4

use mpi

implicit none

integer(kind=4) ierr

call mpi_init(ierr)

call test_rk4()

call mpi_finalize(ierr)

end program test_ts
