program test_grid_field_main

use test_grid_field_mod, only : test_grid_field

use mpi

call MPI_init(ierr)

call test_grid_field()

call mpi_finalize(ierr)

end
