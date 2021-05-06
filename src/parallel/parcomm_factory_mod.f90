module parcomm_factory_mod

use parcomm_mod, only : parcomm_t
use mpi

implicit none

contains

subroutine create_parcomm(parcomm)

    type(parcomm_t), intent(out) :: parcomm

    integer(kind=4) :: ierr

    parcomm%comm_w = mpi_comm_world
    parcomm%myid = parcomm%get_mpi_rank()
    parcomm%np   = parcomm%get_mpi_proc_number()

end subroutine create_parcomm

end module parcomm_factory_mod
