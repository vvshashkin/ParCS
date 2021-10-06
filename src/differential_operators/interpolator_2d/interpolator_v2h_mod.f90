module interpolator_v2h_mod

use grid_field_mod,         only : grid_field_t, tile_field_t
use domain_mod,             only : domain_t
use exchange_abstract_mod,  only : exchange_t
use mesh_mod,               only : tile_mesh_t
use sbp_operator_mod,       only : sbp_operator_t

implicit none

type, public :: interpolator_v2h_t
    class(exchange_t),     allocatable :: exchange
    class(sbp_operator_t), allocatable :: sbp_interp_v2h
contains
    procedure, public :: interp_v2h
end type interpolator_v2h_t

contains

subroutine interp_v2h(this, u_h, v_h, u, v, domain)

    class(interpolator_v2h_t), intent(inout) :: this
    type(grid_field_t),        intent(inout) :: u_h, v_h, u, v
    type(domain_t),            intent(in)    :: domain

    integer(kind=4) :: t

    !WORKAROUND
    !THIS WORKS PROPERLY ONLY WHEN THERE ARE NO CROSS PANEL EXCANGES
    !SINCE U, V ARE NOT NECESARILY VECTOR QUANTITIES
    call this%exchange%do_vec(u, v, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_v2h_tile(u_h%tile(t), v_h%tile(t), u%tile(t), v%tile(t), &
                    this%sbp_interp_v2h, domain%mesh_u%tile(t), &
                    domain%mesh_v%tile(t), domain%mesh_p%tile(t))
    end do

end subroutine interp_v2h

subroutine interp_v2h_tile(u_h, v_h, u, v, sbp_interp_v2h, &
                           mesh_x, mesh_y, mesh_o)

    use tile_mod, only : tile_t

    type(tile_field_t),     intent(inout) :: u_h, v_h
    type(tile_field_t),     intent(inout) :: u, v
    class(sbp_operator_t),  intent(in)    :: sbp_interp_v2h
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_o

    type(tile_t) :: work

    work = tile_t(is = mesh_o%is, ie = mesh_o%ie, &
                  js = mesh_o%js, je = mesh_o%je, &
                  ks = mesh_o%ks, ke = mesh_o%ke)

        call sbp_interp_v2h%apply(u_h, work, mesh_o%globnx, 'x', u)
        call sbp_interp_v2h%apply(v_h, work, mesh_o%globny, 'y', v)

end subroutine interp_v2h_tile

end module interpolator_v2h_mod
