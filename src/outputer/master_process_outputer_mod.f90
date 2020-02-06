module master_process_outputer_mod

use outputer_abstract_mod, only : outputer_t
use grid_function_mod,     only : grid_function_t
use mesh_mod,              only : mesh_t
use exchange_abstract_mod, only : exchange_t
use mpi

    implicit none

type, public, extends(outputer_t) :: master_process_outputer_t

    class(exchange_t),allocatable :: gather_exch
    integer(kind=4)               :: master_id
    character(len=:), allocatable :: write_type !'bin' or 'txt'
    integer(kind=4)               :: rec_num = 0

contains

    procedure, public :: write => master_process_write

end type master_process_outputer_t

contains

subroutine master_process_write(this, f, mesh, ts, te, file_name, rec_num)

    class(master_process_outputer_t), intent(inout) :: this
    integer(kind=4),                  intent(in)    :: ts, te
    type(grid_function_t),            intent(inout) :: f(ts:te)
    type(mesh_t),                     intent(in)    :: mesh(ts:te)
    character(*),                     intent(in)    :: file_name
    integer(kind=4),                  intent(in), &
                                      optional      :: rec_num

    integer(kind=4) :: t, i, j, k, myid, ierr, code

    call this%gather_exch%do(f,ts,te)

    call mpi_comm_rank(mpi_comm_world, myid, ierr)

    if (myid == this%master_id) then

        if (this%write_type == 'txt') then

            if (this%rec_num == 0) then
                open(this%out_stream , file = ''//trim(file_name)//'.txt', status = 'replace')
            else
                open(this%out_stream , file = ''//trim(file_name)//'.txt', action = 'write', position = 'append', status = 'old')
            end if

            do t = ts, te
                do k = f(t)%ks, f(t)%ke
                    do j = f(t)%js, f(t)%je
                        do i = f(t)%is, f(t)%ie
                            this%rec_num = this%rec_num + 1
                            write(this%out_stream, *) f(t)%p(i,j,k)
                        end do
                    end do
                end do
            end do

            close(this%out_stream)

        else if (this%write_type == 'bin') then

            open(this%out_stream , file = ''//trim(file_name)//'.bin', access="direct", recl = 1)

            do t = ts, te
                do k = f(t)%ks, f(t)%ke
                    do j = f(t)%js, f(t)%je
                        do i = f(t)%is, f(t)%ie
                            this%rec_num = this%rec_num + 1
                            write(this%out_stream, rec = this%rec_num) real(f(t)%p(i,j,k), 4)
                        end do
                    end do
                end do
            end do

            close(this%out_stream)

        else

            print*, "Error in master_process_write!!! Wrong write_type argument value!"
            call mpi_abort(mpi_comm_world, code, ierr)
            stop

        end if

    end if

end subroutine master_process_write

end module master_process_outputer_mod
