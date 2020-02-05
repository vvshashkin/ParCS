module master_paneled_outputer_mod

use outputer_abstract_mod, only : outputer_t
use grid_function_mod,     only : grid_function_t
use mesh_mod,              only : mesh_t
use exchange_mod,          only : exchange_t
use mpi

    implicit none

type, public, extends(outputer_t) :: master_paneled_outputer_t

    type(exchange_t)              :: gather_exch
    integer(kind=4)               :: master_id
    character(len=:), allocatable :: write_type !'bin' or 'txt'
    integer(kind=4)               :: rec_num = 1

contains

    procedure, public :: write => master_paneled_write

end type master_paneled_outputer_t

contains

subroutine master_paneled_write(this, f, mesh, ts, te, file_name, rec_num)

    class(master_paneled_outputer_t), intent(inout) :: this
    integer(kind=4),                  intent(in)    :: ts, te
    type(grid_function_t),            intent(inout) :: f(ts:te)
    type(mesh_t),                     intent(in)    :: mesh(ts:te)
    character(*),                     intent(in)    :: file_name
    integer(kind=4),                  intent(in), &
                                      optional      :: rec_num

    real(kind=4), allocatable :: buffer(:,:,:)
    integer(kind=4) irec
    integer(kind=4) panel_s,panel_e, nx, reclen, klev, ks, ke

    integer(kind=4) :: t, i, j, k, myid, ierr, code

    call this%gather_exch%do(f, lbound(f, dim=1), ubound(f, dim=1))

    call mpi_comm_rank(mpi_comm_world, myid, ierr)

    if (myid == this%master_id) then

            if(present(rec_num)) then
                irec = rec_num
            else
                irec = this%rec_num
                this%rec_num = this%rec_num+1
            end if

            nx = mesh(ts)%nx
            ks = mesh(ts)%ks
            ke = mesh(ts)%ke
            klev = ke-ks+1
            panel_s = minval(mesh(ts:te)%panel_ind)
            panel_e = maxval(mesh(ts:te)%panel_ind)
            reclen  = (panel_e-panel_s+1)*nx*nx
            allocate(buffer(nx,nx,panel_s:panel_e))

            open(newunit=this%out_stream, file = trim(file_name), &
                 access="direct", recl = reclen)
            do k = ks, ke
                do t = ts, te
                    do j = f(t)%js, f(t)%je
                        do i = f(t)%is, f(t)%ie
                            buffer(i,j,f(t)%panel_ind) = real(f(t)%p(i,j,k),4)
                        end do
                    end do
                end do

                write(this%out_stream, rec = (irec-1)*klev+k-ks+1) buffer
            end do
            deallocate(buffer)
            close(this%out_stream)

    end if

end subroutine master_paneled_write

end module master_paneled_outputer_mod
