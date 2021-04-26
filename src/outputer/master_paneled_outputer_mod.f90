module master_paneled_outputer_mod

use outputer_abstract_mod, only : outputer_t
use partition_mod,         only : partition_t
use grid_field_mod,        only : grid_field_t
use exchange_abstract_mod, only : exchange_t
use mpi

    implicit none

type, public, extends(outputer_t) :: master_paneled_outputer_t

    class(exchange_t) , allocatable :: gather_exch
    type(grid_field_t)              :: buf
    integer(kind=4)                 :: master_id
    integer(kind=4)                 :: rec_num = 1

contains

    procedure, public :: write => master_paneled_write

end type master_paneled_outputer_t

contains

subroutine master_paneled_write(this, f, partition, file_name, rec_num)

    class(master_paneled_outputer_t), intent(inout) :: this
    type(grid_field_t),               intent(inout) :: f
    type(partition_t),                intent(in)    :: partition
    character(*),                     intent(in)    :: file_name
    integer(kind=4),                  intent(in), &
                                      optional      :: rec_num

    real(kind=4), allocatable :: buffer(:,:,:)
    integer(kind=4) irec, reclen
    integer(kind=4) :: is, ie, js, je, ks, ke
    integer(kind=4) :: t, i, j, k, k0
    integer(kind=4) :: myid, ierr


    call mpi_comm_rank(mpi_comm_world, myid, ierr)

    if (myid == this%master_id) then
        do t = partition%ts, partition%te
            ks = partition%tile(t)%ks; ke = partition%tile(t)%ke
            js = partition%tile(t)%js; je = partition%tile(t)%je
            is = partition%tile(t)%is; ie = partition%tile(t)%ie
            this%buf%block(t)%p(is:ie,js:je,ks:ke) = f%block(t)%p(is:ie,js:je,ks:ke)
        end do
        call this%gather_exch%do(this%buf)
    else
        call this%gather_exch%do(f)
    end if

    if (myid == this%master_id) then
            irec = 1
            if (present(rec_num)) irec = rec_num

            reclen  = partition%num_panels*partition%Nh*partition%Nh
            allocate(buffer(partition%Nh,partition%Nh,1:partition%num_panels))

            open(newunit=this%out_stream, file = trim(file_name), &
                 access="direct", recl = reclen)

            k0 = partition%tile(1)%ks
            print *, "K0, 4to eto???????????????????????????????????", k0 !!!?????????????????????????
            do k = k0, partition%Nz
                do t = 1, partition%num_panels*partition%num_tiles
                    do j = partition%tile(t)%js, partition%tile(t)%je
                        do i = partition%tile(t)%is, partition%tile(t)%ie
                            buffer(i,j,partition%tile(t)%panel_number) = real(this%buf%block(t)%p(i,j,k),4)
                        end do
                    end do
                end do

                write(this%out_stream, rec = (irec-1)*(partition%Nz-k0+1)+(k-k0+1)) buffer
            end do
            deallocate(buffer)
            close(this%out_stream)
    end if

end subroutine master_paneled_write

end module master_paneled_outputer_mod
