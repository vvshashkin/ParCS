module exchange_mod

use exchange_profile_mod, only : exchange_profile_t
use buffer_mod,           only : buffer_t
use grid_function_mod,    only : grid_function_t
implicit none


type, public :: exchange_t

    type(exchange_profile_t), allocatable :: profile
    type(buffer_t),           allocatable :: send_buff(:), recv_buff(:)
    integer(kind=4),          allocatable :: mpi_send_req(:), mpi_recv_req(:)

contains

    procedure, public :: do   => do_exchange

end type exchange_t

contains

subroutine do_exchange(this, f, tile_s, tile_e)

    use mpi

    class(exchange_t),     intent(inout) :: this
    integer(kind=4),       intent(in)    :: tile_s, tile_e
    type(grid_function_t), intent(inout) :: f(tile_s:tile_e)
    integer(kind=4) :: ierr, myid
    integer(kind=4) :: i, ind, ind_recv

    this%mpi_send_req = MPI_REQUEST_NULL
    this%mpi_recv_req = MPI_REQUEST_NULL

call MPI_comm_rank(mpi_comm_world , myid, ierr)


    do i = 1, this%profile%exch_num
        call MPI_irecv(this%recv_buff(i)%p, this%profile%recv_pts_num(i), mpi_real8, this%profile%exchg_proc_id(i), MPI_ANY_TAG, mpi_comm_world, this%mpi_recv_req(i), ierr)
    end do

    do i = 1, this%profile%exch_num

        call pack_to_buf(f(this%profile%recv_tile_ind(i)), this%send_buff(i)%p,    &
             this%profile%send_is(i), this%profile%send_ie(i), &
             this%profile%send_js(i), this%profile%send_je(i), &
             this%profile%send_ks(i), this%profile%send_ke(i), &
             this%profile%first_dim_index(i),                  &
             this%profile%send_i_step(i),                      &
             this%profile%send_j_step(i),                      &
             this%profile%send_pts_num(i) )

        call MPI_isend(this%send_buff(i)%p, this%profile%send_pts_num(i), mpi_real8, this%profile%exchg_proc_id(i), 0, mpi_comm_world, this%mpi_send_req(i), ierr)

    end do

    do ind = 1, this%profile%exch_num
        call mpi_waitany(this%profile%exch_num, this%mpi_recv_req, i, mpi_status_ignore, ierr)

        call unpack_from_buf(f(this%profile%recv_tile_ind(i)), this%recv_buff(i)%p,      &
             this%profile%recv_is(i), this%profile%recv_ie(i), &
             this%profile%recv_js(i), this%profile%recv_je(i), &
             this%profile%recv_ks(i), this%profile%recv_ke(i), &
             this%profile%recv_pts_num(i) )
    end do

end subroutine do_exchange

subroutine unpack_from_buf(f, buf, is, ie, js, je, ks, ke, pts_num)

    type(grid_function_t), intent(inout) :: f
    integer(kind=4),       intent(in)    :: pts_num
    real(kind=8),          intent(inout) :: buf(pts_num)
    integer(kind=4),       intent(in)    :: is, ie, js, je, ks, ke

    integer(kind=4) :: ind, i, j ,k, idx

    idx = 0

    do k = ks, ke
        do j = js, je
            do i = is, ie
                idx = idx + 1
                f.p(i,j,k) = buf(idx)
            end do
        end do
    end do

    if (idx /= pts_num) print*, 'Error in unpacking!'

end subroutine unpack_from_buf


subroutine pack_to_buf(f, buf, is, ie, js, je, ks, ke, first_dim_index, send_i_step, send_j_step, pts_num)

    type(grid_function_t), intent(in)    :: f
    integer(kind=4),       intent(in)    :: pts_num
    real(kind=8),          intent(inout) :: buf(pts_num)
    integer(kind=4),       intent(in)    :: is, ie, js, je, ks, ke
    integer(kind=4),       intent(in)    :: send_i_step, send_j_step
    character(len=1),      intent(in)    :: first_dim_index

    integer(kind=4) :: ind, i, j ,k, idx

    idx = 0

    if (first_dim_index == 'i') then
        if (send_j_step == 1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do j = js, je
                    do i = is, ie
                        idx = idx + 1
                        buf(idx) = f.p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == 1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do j = js, je
                    do i =ie,  is, -1
                        idx = idx + 1
                        buf(idx) = f.p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do j = je, js, -1
                    do i = is, ie
                        idx = idx + 1
                        buf(idx) = f.p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do j = je, js, -1
                    do i =ie,  is, -1
                        idx = idx + 1
                        buf(idx) = f.p(i,j,k)
                    end do
                end do
            end do
        end if

    else if (first_dim_index == 'j') then

        if (send_j_step == 1 .and. send_i_step ==1 ) then
            do k = ks, ke
                do i = is, ie
                    do j = js, je
                        idx = idx + 1
                        buf(idx) = f.p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == 1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do i =ie,  is, -1
                    do j = js, je
                        idx = idx + 1
                        buf(idx) = f.p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == 1 ) then
            do k = ks, ke
                do i = is, ie
                    do j = je, js, -1
                        idx = idx + 1
                        buf(idx) = f.p(i,j,k)
                    end do
                end do
            end do
        else if (send_j_step == -1 .and. send_i_step == -1 ) then
            do k = ks, ke
                do i =ie,  is, -1
                    do j = je, js, -1
                        idx = idx + 1
                        buf(idx) = f.p(i,j,k)
                    end do
                end do
            end do
        end if

        if ( pts_num /= idx ) print*, 'Error in packing!'

    else
        print*, 'Error in packing!!!'
    end if

end subroutine pack_to_buf
end module exchange_mod
