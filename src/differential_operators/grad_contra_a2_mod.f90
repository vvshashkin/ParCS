module grad_contra_a2_mod

use abstract_grad_mod,  only : grad_operator_t
use grid_field_mod,     only : grid_field_t, tile_field_t
use mesh_mod,           only : mesh_t, tile_mesh_t
use halo_mod,           only : halo_t
use domain_mod,         only : domain_t

implicit none

type, public, extends(grad_operator_t) :: grad_contra_a2_t
    class(halo_t), allocatable :: halo_procedure
contains
    procedure, public :: calc_grad => calc_grad_contra_a2
end type grad_contra_a2_t

contains

subroutine calc_grad_contra_a2(this, gx, gy, f, domain, multiplier)
    class(grad_contra_a2_t), intent(inout) :: this
    type(grid_field_t),      intent(inout) :: gx
    type(grid_field_t),      intent(inout) :: gy
    type(grid_field_t),      intent(inout) :: f
    type(domain_t),          intent(in)    :: domain
    real(kind=8), optional, intent(in)     :: multiplier

    integer(kind=4) :: t
    integer(kind=4), parameter :: halo_width=1
    real(kind=8)    :: mult_loc

    mult_loc = 1.0_8
    if (present(multiplier)) mult_loc = multiplier

    call this%halo_procedure%get_halo_scalar(f,domain,halo_width)

    do t = domain%partition%ts, domain%partition%te
        call calc_grad_on_tile(gx%tile(t), gy%tile(t), f%tile(t), domain%mesh_p%tile(t), mult_loc)
    end do

end subroutine calc_grad_contra_a2

subroutine calc_grad_on_tile(gx, gy, f, mesh, multiplier)

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(in)    :: f
    type(tile_mesh_t),      intent(in)    :: mesh
    real(kind=8),           intent(in)    :: multiplier

    real(kind=8) :: hx, mult_loc
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k
    real(kind=8) fdx, fdy

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx

    do k = ks, ke
        do j = js, je
            do i = is, ie
                fdx = 0.5_8*(f%p(i+1,j,k)-f%p(i-1,j,k))/hx
                fdy = 0.5_8*(f%p(i,j+1,k)-f%p(i,j-1,k))/hx
                gx%p(i,j,k) = (mesh%Qi(1,i,j)*fdx + mesh%Qi(2,i,j)*fdy)*multiplier
                gy%p(i,j,k) = (mesh%Qi(3,i,j)*fdy + mesh%Qi(2,i,j)*fdx)*multiplier
            end do
        end do
    end do

end subroutine calc_grad_on_tile

end module grad_contra_a2_mod