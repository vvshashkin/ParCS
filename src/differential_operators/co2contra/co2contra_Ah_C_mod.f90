module co2contra_ah_c_mod

use abstract_co2contra_mod, only : co2contra_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use mesh_mod,               only : tile_mesh_t
use domain_mod,             only : domain_t
use exchange_abstract_mod,  only : exchange_t
use parcomm_mod,            only : parcomm_global
use halo_mod,               only : halo_vec_t
use sbp_operator_mod,       only : sbp_operator_t
use interpolator_w2v_mod,   only : interpolator_w2v_t
use interpolator_v2w_mod,   only : interpolator_v2w_t

implicit none

type, extends(co2contra_operator_t), public :: co2contra_ah_c_sbp_t
    character(len=:), allocatable :: operator_name
    type(interpolator_w2v_t) :: interp_w2v_op
    type(interpolator_v2w_t) :: interp_v2w_op
    type(grid_field_t) :: uw, vw
    contains
        procedure :: transform    => transform_co2contra_ah_c_sbp
!        procedure :: transform2co => transform_contra2co_ah_c_sbp
end type co2contra_ah_c_sbp_t

contains

subroutine transform_co2contra_ah_c_sbp(this, u_contra, v_contra, u_cov, v_cov, domain)
    class(co2contra_ah_c_sbp_t), intent(inout)  :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u_cov, v_cov
    !output:
    type(grid_field_t),           intent(inout) :: u_contra, v_contra

    integer(kind=4) :: t

    !WORKAROUND
    call this%interp_v2w_op%interp_v2w_Ah_C(this%uw, this%vw, u_cov, v_cov, domain)

    do t = domain%partition%ts, domain%partition%te
        call start_co2contra_transform_at_w_tile(this%uw%tile(t), this%vw%tile(t), &
                                       domain%mesh_o%tile(t))
    end do

    !u_contra contains v part at u, v_contra contains u part at v
    !WORKAROUND
    call this%interp_w2v_op%interp_w2v_Ah_C(u_contra, v_contra, this%uw, this%vw, domain)

    do t = domain%partition%ts, domain%partition%te
        call finalize_co2contra_transform_at_uv_tile(u_contra%tile(t), v_contra%tile(t), &
                   u_cov%tile(t), v_cov%tile(t), &
                domain%mesh_y%tile(t), domain%mesh_x%tile(t))
    end do

end subroutine transform_co2contra_ah_c_sbp
subroutine start_co2contra_transform_at_w_tile(uw, vw, mesh_w)

    type(tile_field_t),     intent(inout) :: uw, vw
    type(tile_mesh_t),      intent(in)    :: mesh_w

    real(kind=8)    :: u, v
    integer(kind=4) :: i, j, k

    do k = mesh_w%ks,mesh_w%ke
        do j = mesh_w%js, mesh_w%je
            do i = mesh_w%is, mesh_w%ie
                u = uw%p(i,j,k)*mesh_w%G(i,j)*mesh_w%Qi(2,i,j)
                v = vw%p(i,j,k)*mesh_w%G(i,j)*mesh_w%Qi(2,i,j)
                !Change components for (u,v) to (v,u) interpolation
                uw%p(i,j,k) = v
                vw%p(i,j,k) = u
            end do
        end do
    end do
end subroutine start_co2contra_transform_at_w_tile
subroutine finalize_co2contra_transform_at_uv_tile(u_contra, v_contra, u_cov, v_cov, mesh_u, mesh_v)

    type(tile_field_t),     intent(inout) :: u_contra, v_contra, u_cov, v_cov
    type(tile_mesh_t),      intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k

    do k = mesh_u%ks,mesh_u%ke
        do j = mesh_u%js, mesh_u%je
            do i = mesh_u%is, mesh_u%ie
                u_contra%p(i,j,k) = mesh_u%Qi(1,i,j)*u_cov%p(i,j,k)+u_contra%p(i,j,k)/mesh_u%G(i,j)
            end do
        end do
    end do
    do k = mesh_v%ks,mesh_v%ke
        do j = mesh_v%js, mesh_v%je
            do i = mesh_v%is, mesh_v%ie
                v_contra%p(i,j,k) = mesh_v%Qi(3,i,j)*v_cov%p(i,j,k)+v_contra%p(i,j,k)/mesh_v%G(i,j)
            end do
        end do
    end do

end subroutine finalize_co2contra_transform_at_uv_tile
! subroutine transform_contra2co_ah_c_sbp(this, u_cov, v_cov, u_contra, v_contra, domain)
!     class(co2contra_ah_c_sbp_t), intent(inout) :: this
!     type(domain_t),               intent(in)    :: domain
!     type(grid_field_t),           intent(inout) :: u_contra, v_contra
!     !output:
!     type(grid_field_t),           intent(inout) :: u_cov, v_cov
!
!     integer(kind=4) :: t
!
!     call this%interp_v2h_op%interp_v2h(this%uh, this%vh, u_contra, v_contra, domain)
!
!     do t = domain%partition%ts, domain%partition%te
!         call start_contra2co_transform_at_h_tile(this%uh%tile(t), this%vh%tile(t),&
!                                                  domain%mesh_p%tile(t))
!     end do
!
!     !u_cov contains v part at u, v_cov contains u part at v
!     call this%interp_h2v_op%interp_h2v(u_cov, v_cov, this%uh, this%vh, domain)
!
!     do t = domain%partition%ts, domain%partition%te
!         call finalize_contra2co_transform_at_uv_tile(u_cov%tile(t), v_cov%tile(t), &
!                    u_contra%tile(t), v_contra%tile(t), &
!                    domain%mesh_u%tile(t), domain%mesh_v%tile(t))
!     end do
!
! end subroutine transform_contra2co_ah_c_sbp
! subroutine start_contra2co_transform_at_h_tile(uh, vh, mesh_o)
!
!     type(tile_field_t),     intent(inout) :: uh, vh
!     type(tile_mesh_t),      intent(in)    :: mesh_o
!
!     real(kind=8)    :: u, v
!     integer(kind=4) :: i, j, k
!
!     do k = mesh_o%ks,mesh_o%ke
!         do j = mesh_o%js, mesh_o%je
!             do i = mesh_o%is, mesh_o%ie
!                 u = uh%p(i,j,k)*mesh_o%G(i,j)*mesh_o%Q(2,i,j)
!                 v = vh%p(i,j,k)*mesh_o%G(i,j)*mesh_o%Q(2,i,j)
!                 !Change components for (u,v) to (v,u) interpolation
!                 uh%p(i,j,k) = v
!                 vh%p(i,j,k) = u
!             end do
!         end do
!     end do
! end subroutine start_contra2co_transform_at_h_tile
! subroutine finalize_contra2co_transform_at_uv_tile(u_cov, v_cov, u_contra, v_contra, mesh_u, mesh_v)
!
!     type(tile_field_t),     intent(inout) :: u_contra, v_contra, u_cov, v_cov
!     type(tile_mesh_t),      intent(in)    :: mesh_u, mesh_v
!
!     integer(kind=4) :: i, j, k
!
!     do k = mesh_u%ks,mesh_u%ke
!         do j = mesh_u%js, mesh_u%je
!             do i = mesh_u%is, mesh_u%ie
!                 u_cov%p(i,j,k) = mesh_u%Q(1,i,j)*u_contra%p(i,j,k)+u_cov%p(i,j,k)/mesh_u%G(i,j)
!             end do
!         end do
!     end do
!     do k = mesh_v%ks,mesh_v%ke
!         do j = mesh_v%js, mesh_v%je
!             do i = mesh_v%is, mesh_v%ie
!                 v_cov%p(i,j,k) = mesh_v%Q(3,i,j)*v_contra%p(i,j,k)+v_cov%p(i,j,k)/mesh_v%G(i,j)
!             end do
!         end do
!     end do
!
! end subroutine finalize_contra2co_transform_at_uv_tile

end module co2contra_ah_c_mod
