module div_ah2_mod

use domain_mod,             only : domain_t
use abstract_div_mod,       only : div_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use halo_mod,               only : halo_vec_t
use exchange_abstract_mod,  only : exchange_t
use parcomm_mod,            only : parcomm_global

implicit none

type, public, extends(div_operator_t) :: div_ah2_t
    class(exchange_t), allocatable     :: exch_halo
contains
    procedure, public :: calc_div => calc_div_ah2
end type div_ah2_t

contains

subroutine calc_div_ah2(this, div, u, v, domain, multiplier)
    class(div_ah2_t),        intent(inout) :: this
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

    call this%exch_halo%do_vec(u,v,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        div%tile(t)%p = 0.0_8
        call calc_div_on_tile(div%tile(t), u%tile(t), v%tile(t),  &
                              domain%mesh_xy%tile(t), mult_loc)
    end do


end subroutine calc_div_ah2

subroutine calc_div_on_tile(div, u, v, mesh, multiplier)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),     intent(inout) :: div
    type(tile_field_t),     intent(in)    :: u, v
    type(tile_mesh_t),      intent(in)    :: mesh
    real(kind=8),           intent(in)    :: multiplier

    real(kind=8)    :: hx, mult_loc
    integer(kind=4) :: ks, ke, js, je, is, ie, n, i, j, k

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx

    do k = ks, ke
        if(js==1 .and. is == 1) then
            div%p(1,1,k) = mesh%G(2,1)*(u%p(2,1,k)-u%p(-1,1,k)+ &
                                        v%p(1,2,k)-v%p(1,-1,k)+ &
                                        v%p(0,2,k)+u%p(2,0,k))/ &
                                             (3.0_8*mesh%G(1,1)*hx)*multiplier
        end if
        n = mesh%nx
        if(js==1 .and. ie == n+1) then
            div%p(n+1,1,k) = mesh%G(n,1)*(u%p(n+3,1,k)-u%p(n,1,k)+ &
                                        v%p(n+1,2,k)-v%p(n+1,-1,k)+ &
                                        v%p(n+2,2,k)-u%p(n,0,k))/ &
                                             (3.0_8*mesh%G(n+1,n+1)*hx)*multiplier
        end if
        if(js == 1) then
            j=1
            do i = max(is,2), min(ie,mesh%nx)
                div%p(i,j,k) = (mesh%G(i+1,j)*u%p(i+1,j,k)-mesh%G(i-1,j)*u%p(i-1,j,k) +  &
                                mesh%G(i+1,j)*u%p(i+1,j-1,k)-mesh%G(i-1,j)*u%p(i-1,j-1,k) +  &
                                2._8*mesh%G(i,j+1)*(v%p(i,j+1,k)-v%p(i,j-2,k)))/  &
                                (mesh%G(i,j)*hx)*multiplier*0.25_8
            end do
        end if
        do j = max(js,2), min(je,mesh%ny)
            do i = max(is,2), min(ie,mesh%nx)
                div%p(i,j,k) = (mesh%G(i+1,j)*u%p(i+1,j,k)-mesh%G(i-1,j)*u%p(i-1,j,k) +  &
                                mesh%G(i,j+1)*v%p(i,j+1,k)-mesh%G(i,j-1)*v%p(i,j-1,k))/  &
                                (2._8*mesh%G(i,j)*hx)*multiplier
            end do
        end do
    end do

end subroutine calc_div_on_tile

end module div_ah2_mod
