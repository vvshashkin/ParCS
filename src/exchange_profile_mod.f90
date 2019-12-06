module exchange_profile_mod

use partition_mod, only : partition_t
use tile_mod,      only : tile_t

implicit none

type, public :: exchange_profile_t

    integer(kind=4),  allocatable :: send_is(:), send_ie(:), send_js(:), send_je(:), send_ks(:), send_ke(:)
    integer(kind=4),  allocatable :: recv_is(:), recv_ie(:), recv_js(:), recv_je(:), recv_ks(:), recv_ke(:)

    integer(kind=4),  allocatable :: send_pts_num(:), recv_pts_num(:)
    integer(kind=4),  allocatable :: send_i_step(:), send_j_step(:)
    integer(kind=4),  allocatable :: send_tile_ind(:), recv_tile_ind(:)
    integer(kind=4),  allocatable :: send_proc_id(:), recv_proc_id(:)
    integer(kind=4)               :: send_exch_num, recv_exch_num
    character(len=1), allocatable :: first_dim_index(:)

contains

    procedure, public :: check

end type exchange_profile_t


contains

subroutine check(this)

!This routine perfoms validation check of the exchange_profile instance

    class(exchange_profile_t), intent(in) :: this

    integer :: ind
    logical :: test(21)

    test(1)  = allocated(this%send_is)
    test(2)  = allocated(this%send_ie)
    test(3)  = allocated(this%send_js)
    test(4)  = allocated(this%send_je)
    test(5)  = allocated(this%send_ks)
    test(6)  = allocated(this%send_ke)
    test(7)  = allocated(this%recv_is)
    test(8)  = allocated(this%recv_ie)
    test(9)  = allocated(this%recv_js)
    test(10) = allocated(this%recv_je)
    test(11) = allocated(this%recv_ks)
    test(12) = allocated(this%recv_ke)
    test(13) = allocated(this%send_pts_num)
    test(14) = allocated(this%recv_pts_num)
    test(15) = allocated(this%send_i_step)
    test(16) = allocated(this%send_j_step)
    test(17) = allocated(this%send_tile_ind)
    test(18) = allocated(this%recv_tile_ind)
    test(19) = allocated(this%send_proc_id)
    test(20) = allocated(this%recv_proc_id)
    test(21) = allocated(this%first_dim_index)

    if ( not(all(test)) ) then
        print*, 'Error in exchange_profile check!!! Some arrays are not allocated!'
        stop
    end if

    test(1) = all(this%send_is <= this%send_ie)
    test(2) = all(this%send_js <= this%send_je)
    test(3) = all(this%send_ks <= this%send_ke)
    test(4) = all(this%recv_is <= this%recv_ie)
    test(5) = all(this%recv_js <= this%recv_je)
    test(6) = all(this%recv_ks <= this%recv_ke)

    if ( not(all(test(1:6))) ) then
        print*, 'Error in exchange_profile check!!! Problem with send/recv index validity!'
        stop
    end if

    test(1) = all(this%send_pts_num>0)
    test(2) = all(this%recv_pts_num>0)

    if ( not(all(test(1:2))) ) then
        print*, 'Error in exchange_profile check!!! Problem with send/recv points number!'
        stop
    end if

end subroutine check


end module exchange_profile_mod
