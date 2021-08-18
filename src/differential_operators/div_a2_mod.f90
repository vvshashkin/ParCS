module div_a2_mod

use domain_mod,         only : domain_t
use abstract_div_mod,   only : div_operator_t
use grid_field_mod,     only : grid_field_t, tile_field_t
use halo_mod,           only : halo_vec_t
use parcomm_mod,        only : parcomm_global

implicit none

type, public, extends(div_operator_t) :: div_a2_t
    class(halo_vec_t), allocatable :: halo_procedure
    character(:),      allocatable :: subtype
contains
    procedure, public :: calc_div => calc_div_a2
end type div_a2_t

contains

subroutine calc_div_a2(this, div, u, v, domain, multiplier)
    class(div_a2_t),        intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: u, v
    real(kind=8), optional, intent(in)    :: multiplier
    !output
    type(grid_field_t),     intent(inout) :: div

    integer(kind=4), parameter :: halo_width = 1
    integer(kind=4) :: t
    real(kind=8)    :: mult_loc

    mult_loc = 1.0_8
    if (present(multiplier)) mult_loc = multiplier

    call this%halo_procedure%get_halo_vector(u,v,domain,halo_width)
    select case(this%subtype)
    case("cons")
        do t = domain%partition%ts, domain%partition%te
            call calc_div_on_tile_cons(div%tile(t), u%tile(t), v%tile(t), &
                                       domain%mesh_p%tile(t), domain%partition%Nh, mult_loc)
        end do
    case("default")
        do t = domain%partition%ts, domain%partition%te
            call calc_div_on_tile(div%tile(t), u%tile(t), v%tile(t),  &
                                  domain%mesh_p%tile(t), mult_loc)
        end do
    case default
        call parcomm_global%abort("div_a2_mod, calc_div_a2, unknown subtype: "// this%subtype)
    end select


end subroutine calc_div_a2

subroutine calc_div_on_tile(div, u, v, mesh, multiplier)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),     intent(inout) :: div
    type(tile_field_t),     intent(in)    :: u, v
    type(tile_mesh_t),      intent(in)    :: mesh
    real(kind=8),           intent(in)    :: multiplier

    real(kind=8)    :: hx, mult_loc
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx
    do k = ks, ke
        do j = js, je
            do i = is, ie
                div%p(i,j,k) = (mesh%G(i+1,j)*u%p(i+1,j,k)-mesh%G(i-1,j)*u%p(i-1,j,k) +  &
                                mesh%G(i,j+1)*v%p(i,j+1,k)-mesh%G(i,j-1)*v%p(i,j-1,k))/  &
                                (2._8*mesh%G(i,j)*hx)*multiplier
            end do
        end do
    end do

end subroutine calc_div_on_tile

subroutine calc_div_on_tile_cons(div, u, v, mesh, nx, multiplier)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),     intent(inout) :: div
    type(tile_field_t),     intent(in)    :: u, v
    type(tile_mesh_t),      intent(in)    :: mesh
    integer(kind=4),        intent(in)    :: nx
    real(kind=8),           intent(in)    :: multiplier

    real(kind=8)    :: hx, mult_loc
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k
    integer(kind=4) :: ip1, im1, jp1, jm1

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx
    do k = ks, ke
        do j = js, je
            jm1 = max(1,j-1)
            jp1 = min(nx,j+1)
            do i = is, ie
                im1 = max(1,i-1)
                ip1 = min(nx,i+1)    
                div%p(i,j,k) = (mesh%G(ip1,j)*u%p(i+1,j,k)-mesh%G(im1,j)*u%p(i-1,j,k) +  &
                                mesh%G(i,jp1)*v%p(i,j+1,k)-mesh%G(i,jm1)*v%p(i,j-1,k))/  &
                                (2._8*mesh%G(i,j)*hx)*multiplier
            end do
        end do
    end do

end subroutine calc_div_on_tile_cons

end module div_a2_mod
