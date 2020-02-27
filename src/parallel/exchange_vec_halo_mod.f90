module exchange_vec_halo_mod

use grid_function_mod,         only : grid_function_t
use exchange_vec_abstract_mod, only : exchange_vec_t
use buffer_mod,                only : buffer_t, pack_to_buf_vec, unpack_from_buf_vec
use mpi

implicit none

type, extends(exchange_vec_t), public :: exchange_vec_2D_halo_t

    type(buffer_t), allocatable :: send_buff(:), recv_buff(:)

    integer(kind=4) :: mpi_message_type = mpi_real8

    integer(kind=4) :: recv_number, send_number

    integer(kind=4), allocatable, dimension(:) :: mpi_send_req, mpi_recv_req

    integer(kind=4), allocatable, dimension(:) :: recv_points_num, send_points_num
    integer(kind=4), allocatable, dimension(:) :: recv_to_tile_ind, send_from_tile_ind
    integer(kind=4), allocatable, dimension(:) :: send_to_proc_id, recv_from_proc_id
    integer(kind=4), allocatable, dimension(:) :: send_tag, recv_tag

    integer(kind=4), allocatable, dimension(:) :: send_is, send_ie, send_js, send_je, send_ks, send_ke
    integer(kind=4), allocatable, dimension(:) :: recv_is, recv_ie, recv_js, recv_je, recv_ks, recv_ke

    integer(kind=4) :: ts, te

    integer(kind=4),  allocatable :: send_i_step(:), send_j_step(:)
    character(len=1), allocatable :: first_dim_index(:)

contains

    procedure, public:: do => do_halo_exchange

end type exchange_vec_2D_halo_t

contains

subroutine do_halo_exchange(this, u, v, ts, te)

    class(exchange_vec_2D_halo_t), intent(inout) :: this
    integer(kind=4),        intent(in)           :: ts, te
    type(grid_function_t),  intent(inout)        :: u(ts:te), v(ts:te)

    integer(kind=4) :: ierr, myid
    integer(kind=4) :: i, ind, ind_recv

    this%mpi_send_req = MPI_REQUEST_NULL
    this%mpi_recv_req = MPI_REQUEST_NULL

    do i = 1, this%recv_number
        call MPI_irecv(this%recv_buff(i)%p,        &
                       2*this%recv_points_num(i),    &
                       this%mpi_message_type,      &
                       this%recv_from_proc_id(i),  &
                       this%recv_tag(i),   &
                       mpi_comm_world,             &
                       this%mpi_recv_req(i),       &
                       ierr)
    end do

    do i = 1, this%send_number
        call pack_to_buf_vec(                  &
             u(this%send_from_tile_ind(i)),    &
             v(this%send_from_tile_ind(i)),    &
             this%send_buff(i)%p,              &
             this%send_is(i), this%send_ie(i), &
             this%send_js(i), this%send_je(i), &
             this%send_ks(i), this%send_ke(i), &
             this%first_dim_index(i),          &
             this%send_i_step(i),              &
             this%send_j_step(i),              &
             this%send_points_num(i) )
    end do

    do i = 1, this%send_number
        call MPI_isend(this%send_buff(i)%p,      &
                       2*this%send_points_num(i),  &
                       this%mpi_message_type,    &
                       this%send_to_proc_id(i),  &
                       this%send_tag(i),         &
                       mpi_comm_world,           &
                       this%mpi_send_req(i),     &
                       ierr )
    end do

    do ind = 1, this%recv_number
        call mpi_waitany(this%recv_number, this%mpi_recv_req, i, mpi_status_ignore, ierr)

        call unpack_from_buf_vec(              &
             u(this%recv_to_tile_ind(i)),      &
             v(this%recv_to_tile_ind(i)),      &
             this%recv_buff(i)%p,              &
             this%recv_is(i), this%recv_ie(i), &
             this%recv_js(i), this%recv_je(i), &
             this%recv_ks(i), this%recv_ke(i), &
             this%recv_points_num(i) )

    end do

    call mpi_waitall(this%send_number, this%mpi_send_req, mpi_statuses_ignore, ierr)

end subroutine do_halo_exchange

end module exchange_vec_halo_mod
