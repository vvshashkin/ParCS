!Object to interpolate vector values to virtual points beyond face edge
!non-staggered (A) grid. The algorithm is as follows:
!1)transform adjacent face's vector components to this face components
!2)make halo-interpolation for u,v as for scalars (non-staggered grid!)

module ecs_halo_vec_a_mod

use halo_mod,     only : halo_vec_t
use ecs_halo_mod, only : ecs_halo_t

implicit none

type, extends(halo_vec_t) :: ecs_halo_vec_t
  integer(kind=4) n          !number of real grid points along cubed-sphere face edge
  integer(kind=4) panel_ind
  integer(kind=4) ish, ieh, jsh, jeh
  integer(kind=4) halo_width !number of rows in halo-zone
  logical lhalo(4)   !flags to perform haloprocedures at specific edge and corner

  !Transformation matrices source panel vect ->target panel vect
  real(kind=8)   , allocatable :: TM1(:,:,:) !lower edge
  real(kind=8)   , allocatable :: TM2(:,:,:) !upper edge
  real(kind=8)   , allocatable :: TM3(:,:,:) !left edge
  real(kind=8)   , allocatable :: TM4(:,:,:) !right edge

  !will pointer below be better (memory efficiency)?
  class(ecs_halo_t), allocatable :: halo !scalar halo to interpolate individual
                                         !vector components after transformation
  contains
  procedure, public :: interpv  => ext_halo_vec_a_ecs

end type ecs_halo_vec_t

contains

subroutine ext_halo_vec_a_ecs(this, u, v, halo_width)
!interpolate source face vectors at halo zones to target face virtual points
use grid_function_mod, only: grid_function_t

class(ecs_halo_vec_t), intent(in)    :: this
type(grid_function_t), intent(inout) :: u, v
integer(kind=4),       intent(in)    :: halo_width
!locals
integer(kind=4) i, j, k, ks, ke
real(kind=8)    zu, zv

ks = u%ks - u%nvk
ke = u%ke + u%nvk

!Apply transform matrix
if(this%lhalo(1)) then
    do k = ks, ke
        do j=1,halo_width
            do i=this%ish, this%ieh
                zu = u%p(i,1-j,k)
                zv = v%p(i,1-j,k)
                !if(this%panel_ind == 1 .and. j==1 .and. k==1) then
                !    print '(i4,6f15.7)', i, u%p(i,1,k), v%p(i,1,k), zu, zv, &
                !                         this%TM1(1,i,j)*zu + this%TM1(2,i,j)*zv, &
                !                         this%TM1(3,i,j)*zu + this%TM1(4,i,j)*zv
                !end if
                u%p(i,1-j,k) = this%TM1(1,i,j)*zu + this%TM1(2,i,j)*zv
                v%p(i,1-j,k) = this%TM1(3,i,j)*zu + this%TM1(4,i,j)*zv
            end do
        end do
    end do
end if
if(this%lhalo(2)) then
    do k = ks, ke
        do j=1,halo_width
            do i=this%ish, this%ieh
                zu = u%p(i,this%n+j,k)
                zv = v%p(i,this%n+j,k)
                u%p(i,this%n+j,k) = this%TM2(1,i,j)*zu + this%TM2(2,i,j)*zv
                v%p(i,this%n+j,k) = this%TM2(3,i,j)*zu + this%TM2(4,i,j)*zv
            end do
        end do
    end do
end if
if(this%lhalo(3)) then
    do k = ks, ke
        do j=1,halo_width
            do i=this%ish, this%ieh
                zu = u%p(1-j,i,k)
                zv = v%p(1-j,i,k)
                u%p(1-j,i,k) = this%TM3(1,i,j)*zu + this%TM3(2,i,j)*zv
                v%p(1-j,i,k) = this%TM3(3,i,j)*zu + this%TM3(4,i,j)*zv
            end do
        end do
    end do
end if
if(this%lhalo(4)) then
    do k = ks, ke
        do j=1,halo_width
            do i=this%ish, this%ieh
                zu = u%p(this%n+j,i,k)
                zv = v%p(this%n+j,i,k)
                u%p(this%n+j,i,k) = this%TM4(1,i,j)*zu + this%TM4(2,i,j)*zv
                v%p(this%n+j,i,k) = this%TM4(3,i,j)*zu + this%TM4(4,i,j)*zv
            end do
        end do
    end do
end if



call this%halo%interp(u,halo_width)
call this%halo%interp(v,halo_width)

end subroutine ext_halo_vec_a_ecs

end module ecs_halo_vec_a_mod
