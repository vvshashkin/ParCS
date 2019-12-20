module master_process_outputer_mod

use outputer_abstract_mod, only : outputer_t
use grid_function_mod,     only : grid_function_t
use exchange_mod,          only : exchange_t
use mpi

    implicit none

type, public, extends(outputer_t) :: master_process_outputer_t

    type(exchange_t) :: gather_exch
    integer(kind=4)  :: master_id

contains

    procedure, public :: write => master_process_write

end type master_process_outputer_t

contains

subroutine master_process_write(this, f, ts, te, file_name)

    class(master_process_outputer_t), intent(inout) :: this
    integer(kind=4),                  intent(in)    :: ts, te
    type(grid_function_t),            intent(inout) :: f(ts:te)
    character(*),                     intent(in)    :: file_name

    integer(kind=4) :: t, i, j, k, myid, ierr

    call this%gather_exch%do(f, lbound(f, dim=1), ubound(f, dim=1))

    call mpi_comm_rank(mpi_comm_world, myid, ierr)

    if (myid == this%master_id) then

        open(this%out_stream , file = ''//trim(file_name)//'.txt')

        do t = ts, te
            do k = f(t)%ks, f(t)%ke
                do j = f(t)%js, f(t)%je
                    do i = f(t)%is, f(t)%ie
                        write(this%out_stream, *) f(t)%p(i,j,k)
                    end do
                end do
            end do
        end do

        close(this%out_stream)

    end if

end subroutine master_process_write

end module master_process_outputer_mod
