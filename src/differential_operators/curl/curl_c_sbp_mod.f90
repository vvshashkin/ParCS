module curl_c_sbp_mod

use grid_field_mod,        only : grid_field_t, tile_field_t
use domain_mod,            only : domain_t
use abstract_curl_mod,     only : curl_operator_t
use sbp_operator_mod,      only : sbp_operator_t
use exchange_abstract_mod, only : exchange_t

implicit none

type, public, extends(curl_operator_t) :: curl_c_sbp_t
    class(exchange_t),     allocatable :: exchange
    class(sbp_operator_t), allocatable :: sbp_diff
contains
    procedure, public :: calc_curl
end type curl_c_sbp_t

contains

subroutine calc_curl(this, curl, u, v, domain)

    class(curl_c_sbp_t),   intent(inout) :: this
    type(domain_t),        intent(in)    :: domain
    type(grid_field_t),    intent(inout) :: u, v !covariant components
    type(grid_field_t),    intent(inout) :: curl

    integer(kind=4) :: t

    call this%exchange%do_vec(u,v,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_curl_c_sbp_tile(curl%tile(t), u%tile(t), v%tile(t),   &
                                  this%sbp_diff, domain%mesh_x%tile(t), &
                                  domain%mesh_y%tile(t), domain%mesh_xy%tile(t), &
                                  domain%mesh_xy%scale)
    end do

end subroutine calc_curl

subroutine calc_curl_c_sbp_tile(curl, u, v, sbp_diff, mesh_x, mesh_y, mesh_xy, scale)
    use tile_mod,       only : tile_t
    use mesh_mod,       only : tile_mesh_t
    use grid_field_mod, only : grid_field_t

    type(tile_field_t),    intent(inout)    :: u, v
    class(sbp_operator_t), intent(in)    :: sbp_diff
    type(tile_mesh_t),     intent(in)    :: mesh_x, mesh_y, mesh_xy
    real(kind=8),          intent(in)    :: scale

    type(tile_field_t),    intent(inout) :: curl

    integer(kind=4) :: i, j, k
    integer(kind=4) :: is, ie, js, je, ks, ke
    real(kind=8) :: dvx(mesh_xy%is:mesh_xy%ie,mesh_xy%js:mesh_xy%je,1)
    real(kind=8) :: duy(mesh_xy%is:mesh_xy%ie,mesh_xy%js:mesh_xy%je,1)
    type(tile_t) :: dxdy_tile

    is = mesh_xy%is; ie = mesh_xy%ie
    js = mesh_xy%js; je = mesh_xy%je
    ks = mesh_xy%ks; ke = mesh_xy%ke

    dxdy_tile = tile_t(is=is, ie=ie, js=js, je=je, ks=1, ke=1)

    do k=ks, ke
        dxdy_tile%ks = k; dxdy_tile%ke = k

        call sbp_diff%apply(dvx, dxdy_tile, dxdy_tile, mesh_xy%nx+1, 'x', v)

        call sbp_diff%apply(duy, dxdy_tile, dxdy_tile, mesh_xy%ny+1, 'y', u)

        do j=js, je
            do i=is, ie
                curl%p(i,j,k) = (dvx(i,j,1)-duy(i,j,1)) / (mesh_xy%G(i,j)*mesh_xy%hx*scale)
            end do
        end do
    end do

end subroutine calc_curl_c_sbp_tile


end module curl_c_sbp_mod
