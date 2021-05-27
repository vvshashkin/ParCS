module parcomm_mod

use mpi

implicit none

type, public :: parcomm_t
    integer(kind=4)                :: myid, np
    integer(kind=mpi_integer_kind) :: comm_w
contains
    procedure, public :: get_mpi_rank
    procedure, public :: get_mpi_proc_number
    procedure, public :: barrier
    procedure, public :: print
end type parcomm_t

contains

function get_mpi_rank(this) result(myid)

    class(parcomm_t), intent(in) :: this
    integer(kind=4) :: myid, ierr

    call MPI_comm_rank(this%comm_w, myid, ierr)

end function get_mpi_rank

function get_mpi_proc_number(this) result(np)

    class(parcomm_t), intent(in) :: this
    integer(kind=4) :: np, ierr

    call MPI_comm_size(this%comm_w, np, ierr)

end function get_mpi_proc_number

subroutine barrier(this)

    class(parcomm_t), intent(in) :: this
    integer(kind=4) :: np, ierr

    call mpi_barrier(this%comm_w, ierr)

end subroutine barrier

subroutine print(this, message)

    class(parcomm_t), intent(in) :: this
    character(len=*), intent(in) :: message

    if (this%myid == 0) print*, message

end subroutine print

end module parcomm_mod
