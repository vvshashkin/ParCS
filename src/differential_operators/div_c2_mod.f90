module div_c2_mod

use abstract_div_mod,   only : div_operator_t
use domain_mod,         only : domain_t
use mesh_mod,           only : mesh_t, tile_mesh_t
use grid_field_mod,     only : grid_field_t, tile_field_t
use exchange_halo_mod,  only : exchange_t

implicit none

type, public, extends(div_operator_t) :: div_c2_t
    class(exchange_t), allocatable :: exch_halo
contains
    procedure, public :: calc_div => calc_div_c2
end type div_c2_t

contains

subroutine calc_div_c2(this, div, u, v, domain, multiplier)

    class(div_c2_t),        intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: u, v
    real(kind=8), optional, intent(in)    :: multiplier
    !out put
    type(grid_field_t),     intent(inout) :: div

    integer(kind=4) :: i, j, k, t
    real(kind=8) mult_loc, hx

    mult_loc = 1.0_8
    if(present(multiplier)) mult_loc=multiplier

    call this%exch_halo%do_vec(u,v,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_div_on_tile(div%tile(t), u%tile(t), v%tile(t),            &
                              domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                              domain%mesh_o%tile(t), mult_loc)
    end do

end subroutine calc_div_c2

subroutine calc_div_on_tile(div, u, v, mesh_u, mesh_v, mesh_p, multiplier)

    type(tile_field_t),  intent(inout) :: div
    type(tile_field_t),  intent(in)    :: u, v
    type(tile_mesh_t),   intent(in)    :: mesh_u, mesh_v, mesh_p
    real(kind=8), optional, intent(in) :: multiplier

    real(kind=8) :: hx
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k

    is = mesh_p%is; ie = mesh_p%ie
    js = mesh_p%js; je = mesh_p%je
    ks = mesh_p%ks; ke = mesh_p%ke

    hx = mesh_p%hx

    do k = ks, ke
        do j = js, je
            do i = is, ie
                div%p(i,j,k) = (mesh_u%G(i+1,j)*u%p(i+1,j,k)-mesh_u%G(i,j)*u%p(i,j,k)  +  &
                                mesh_v%G(i,j+1)*v%p(i,j+1,k)-mesh_v%G(i,j)*v%p(i,j,k))/  &
                                (mesh_p%G(i,j)*hx)*multiplier
            end do
        end do
    end do

end subroutine calc_div_on_tile

end module div_c2_mod
