module exchange_profile_mod

use partition_mod, only : partition_t
use tile_mod,      only : tile_t

implicit none

type, public :: exchange_profile_t

    integer(kind=4),  allocatable :: send_is(:), send_ie(:), send_js(:), send_je(:), send_ks(:), send_ke(:)
    integer(kind=4),  allocatable :: recv_is(:), recv_ie(:), recv_js(:), recv_je(:), recv_ks(:), recv_ke(:)
    integer(kind=4),  allocatable :: send_pts_num(:), recv_pts_num(:)
    integer(kind=4),  allocatable :: send_i_step(:), send_j_step(:)
    integer(kind=4),  allocatable :: send_tile_ind(:), recv_tile_ind(:),  exchg_proc_id(:)
    integer(kind=4)               :: exch_num
    character(len=1), allocatable :: first_dim_index(:)

contains

end type exchange_profile_t


contains


! subroutine init_2d_cross_halo_exchange(this, partition, halo_width, myid)
!
!     class(exchange_profile_t), intent(inout) :: this
!     type(partition_t), target, intent(in)    :: partition
!     integer(kind=4),           intent(in)    :: halo_width !halo width
!     integer(kind=4),           intent(in)    :: myid
!
!     type(tile_t), pointer :: local_tile, remote_tile
!
!     integer(kind=4),  dimension(6*partition%num_tiles) :: recv_is, recv_ie, recv_js, recv_je, recv_ks, recv_ke, &
!                                                           send_is, send_ie, send_js, send_je, send_ks, send_ke, &
!                                                           send_i_step, send_j_step, exchg_tile_ind
!     integer(kind=4),  dimension(6*partition%num_tiles) :: send_pts_num, recv_pts_num
!     character(len=1), dimension(6*partition%num_tiles) :: first_dim_index
!
!
!     integer(kind=4) :: local_ind, remote_ind, exch_num
!     logical :: is_intersection
!
!     if (halo_width > partition%npoints) then
!         write(*,*) 'Error! Halo is too wide! Abort!'
!         stop
!     end if
!
!     exch_num = 0
!
!     do local_ind = 1, 6*partition%num_tiles
!
!         if (partition%proc_map(local_ind) /= myid) cycle
!
!         local_tile => partition%tile(local_ind)
!
!         do remote_ind = 1, 6*partition%num_tiles
!
!             remote_tile => partition%tile(remote_ind)
!
!             call find_cross_halo_intersection(local_tile, remote_tile, halo_width, partition%npoints, is_intersection , &
!                                               recv_is(exch_num + 1), recv_ie(exch_num + 1)                            , &
!                                               recv_js(exch_num + 1), recv_je(exch_num + 1)                            , &
!                                               send_is(exch_num + 1), send_ie(exch_num + 1)                            , &
!                                               send_js(exch_num + 1), send_je(exch_num + 1) )
!
!             if (is_intersection) then
!                 exch_num = exch_num + 1
!                 exchg_tile_ind(exch_num) = remote_ind
!
!                 send_ks(exch_num) = local_tile%ks
!                 send_ke(exch_num) = local_tile%ke
!
!                 send_pts_num(exch_num) = (send_ie(exch_num) - send_is(exch_num) + 1)* &
!                                          (send_je(exch_num) - send_js(exch_num) + 1)* &
!                                          (send_ke(exch_num) - send_ks(exch_num) + 1)
!
!                 recv_pts_num(exch_num) = (recv_ie(exch_num) - recv_is(exch_num) + 1)* &
!                                          (recv_je(exch_num) - recv_js(exch_num) + 1)* &
!                                          (recv_ke(exch_num) - recv_ks(exch_num) + 1)
!             end if
!
!         end do
!
!     end do
!
!     print*, myid, exchg_tile_ind(1:exch_num)
!
!     allocate(this%send_is, source = send_is(1:exch_num))
!     allocate(this%send_ie, source = send_ie(1:exch_num))
!
!     allocate(this%send_js, source = send_js(1:exch_num))
!     allocate(this%send_je, source = send_je(1:exch_num))
!
!     allocate(this%send_ks, source = send_ks(1:exch_num))
!     allocate(this%send_ke, source = send_ke(1:exch_num))
!
!     allocate(this%recv_is, source = recv_is(1:exch_num))
!     allocate(this%recv_ie, source = recv_ie(1:exch_num))
!
!     allocate(this%recv_js, source = recv_js(1:exch_num))
!     allocate(this%recv_je, source = recv_je(1:exch_num))
!
!     allocate(this%recv_ks, source = recv_ks(1:exch_num))
!     allocate(this%recv_ke, source = recv_ke(1:exch_num))
!
!     allocate(this%send_i_step, source = send_i_step(1:exch_num))
!     allocate(this%send_j_step, source = send_j_step(1:exch_num))
!
!     allocate(this%first_dim_index, source = first_dim_index(1:exch_num))
!
!     allocate(this%exchg_tile_ind, source = exchg_tile_ind(1:exch_num))
!
!     allocate(this%send_pts_num, source = send_pts_num(1:exch_num))
!     allocate(this%recv_pts_num, source = recv_pts_num(1:exch_num))
!
! end subroutine init_2d_cross_halo_exchange
!
! subroutine find_cross_halo_intersection(l_tile, r_tile, halo_width, npoints, is_intersection, &
!                                         recv_is, recv_ie, recv_js, recv_je,                   &
!                                         send_is, send_ie, send_js, send_je )
!
!     use topology_mod, only : transform_index
!
!     type(tile_t),    intent(in)  :: l_tile, r_tile
!     integer(kind=4), intent(in)  :: halo_width, npoints
!     logical,         intent(out) :: is_intersection
!     integer(kind=4), intent(out) :: recv_is, recv_ie, recv_js, recv_je, &
!                                     send_is, send_ie, send_js, send_je
!
!     integer(kind=4)  :: r_tile_is, r_tile_ie, r_tile_js, r_tile_je
!     integer(kind=4)  :: is_r, ie_r, js_r, je_r
!     integer(kind=4)  :: is_l, ie_l, js_l, je_l
!     integer(kind=4)  :: vertex_i(4), vertex_j(4)
!     integer(kind=4)  :: i_step, j_step
!     character(len=1) :: fisrt_dim
!     is_intersection = .true.
!
!     if (l_tile%panel_number == r_tile%panel_number) then
!         r_tile_is = r_tile%is; r_tile_ie = r_tile%ie;
!         r_tile_js = r_tile%js; r_tile_je = r_tile%je;
!     else
!         call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%is, r_tile%js, vertex_i(1), vertex_j(1))
!         call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%is, r_tile%je, vertex_i(2), vertex_j(2))
!         call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%ie, r_tile%js, vertex_i(3), vertex_j(3))
!         call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%ie, r_tile%je, vertex_i(4), vertex_j(4), i_step, j_step, fisrt_dim)
!
!         r_tile_is = minval(vertex_i); r_tile_ie = maxval(vertex_i)
!         r_tile_js = minval(vertex_j); r_tile_je = maxval(vertex_j)
!
!     end if
!
!     !left recv halo
!     is_l = l_tile%is-halo_width; ie_l = l_tile%is-1
!     js_l = l_tile%js           ; je_l = l_tile%je
!
!     ! if (r_tile%panel_number==4) then
!     !     print*, is_l, ie_l, js_l, je_l
!     ! end if
!
!     recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
!     recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)
!
!     !left send halo
!     if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then
!
!         is_l = l_tile%is; ie_l = l_tile%ie
!         js_l = l_tile%js; je_l = l_tile%je
!
!         is_r = r_tile_ie+1; ie_r = r_tile_ie+halo_width
!         js_r = r_tile_js  ; je_r = r_tile_je
!
!         send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
!         send_js = max(js_l, js_r); send_je = min(je_l, je_r)
!
!         return
!     end if
!
!     !right recv halo
!     is_l = l_tile%ie+1; ie_l = l_tile%ie+halo_width
!     js_l = l_tile%js  ; je_l = l_tile%je
!
!     recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
!     recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)
!
!     !right send halo
!     if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then
!
!         is_l = l_tile%is; ie_l = l_tile%ie
!         js_l = l_tile%js; je_l = l_tile%je
!
!         is_r = r_tile_is - halo_width; ie_r = r_tile_is - 1
!         js_r = r_tile_js             ; je_r = r_tile_je
!
!         send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
!         send_js = max(js_l, js_r); send_je = min(je_l, je_r)
!
!         return
!     end if
!
!     !bottom recv halo
!     is_l = l_tile%is           ; ie_l = l_tile%ie
!     js_l = l_tile%js-halo_width; je_l = l_tile%js-1
!
!     recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
!     recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)
!
!     !bottom send halo
!     if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then
!
!         is_l = l_tile%is; ie_l = l_tile%ie
!         js_l = l_tile%js; je_l = l_tile%je
!
!         is_r = r_tile_is    ; ie_r = r_tile_ie
!         js_r = r_tile_je + 1; je_r = r_tile_je + halo_width
!
!         send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
!         send_js = max(js_l, js_r); send_je = min(je_l, je_r)
!
!         return
!     end if
!
!     !top recv halo
!     is_l = l_tile%is  ; ie_l = l_tile%ie
!     js_l = l_tile%je+1; je_l = l_tile%je+halo_width
!
!     recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
!     recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)
!
!     !top send halo
!     if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then
!
!         is_l = l_tile%is; ie_l = l_tile%ie
!         js_l = l_tile%js; je_l = l_tile%je
!
!         is_r = r_tile_is             ; ie_r = r_tile_ie
!         js_r = r_tile_js - halo_width; je_r = r_tile_js - 1
!
!         send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
!         send_js = max(js_l, js_r); send_je = min(je_l, je_r)
!
!         return
!     end if
!
!     is_intersection = .false.
!
! end subroutine find_cross_halo_intersection

end module exchange_profile_mod
