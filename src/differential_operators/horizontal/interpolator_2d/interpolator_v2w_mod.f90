module interpolator_v2w_mod

use grid_field_mod,         only : grid_field_t, tile_field_t
use domain_mod,             only : domain_t
use exchange_abstract_mod,  only : exchange_t
use mesh_mod,               only : tile_mesh_t
use sbp_operator_mod,       only : sbp_operator_t

implicit none

type, public :: interpolator_v2w_t
    class(exchange_t),     allocatable :: exchange
    class(sbp_operator_t), allocatable :: sbp_interp_v2w
contains
    procedure, public :: interp_v2w
    procedure, public :: interp_v2w_Ah_C
end type interpolator_v2w_t

contains

subroutine interp_v2w(this, u_w, v_w, u, v, domain)

    class(interpolator_v2w_t), intent(inout) :: this
    type(grid_field_t),        intent(inout) :: u_w, v_w, u, v
    type(domain_t),            intent(in)    :: domain

    integer(kind=4) :: t


end subroutine interp_v2w

subroutine interp_v2w_Ah_C(this, u_w, v_w, u, v, domain)

    class(interpolator_v2w_t), intent(inout) :: this
    type(grid_field_t),        intent(inout) :: u_w, v_w, u, v
    type(domain_t),            intent(in)    :: domain

    integer(kind=4) :: t

    !WORKAROUND
    !THIS WORKS PROPERLY ONLY WHEN THERE ARE NO CROSS PANEL EXCANGES
    !SINCE U, V ARE NOT NECESARILY VECTOR QUANTITIES
    call this%exchange%do_vec(v, u, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_v2w_tile(u_w%tile(t), v_w%tile(t), u%tile(t), v%tile(t), &
                    this%sbp_interp_v2w, domain%mesh_o%tile(t))
    end do

end subroutine interp_v2w_Ah_C

subroutine interp_v2w_tile(u_w, v_w, u, v, sbp_interp_v2w, mesh_o)

    use tile_mod, only : tile_t

    type(tile_field_t),     intent(inout) :: u_w, v_w
    type(tile_field_t),     intent(inout) :: u, v
    class(sbp_operator_t),  intent(in)    :: sbp_interp_v2w
    type(tile_mesh_t),      intent(in)    :: mesh_o

    type(tile_t) :: work

    work = tile_t(is = mesh_o%is, ie = mesh_o%ie, &
                  js = mesh_o%js, je = mesh_o%je, &
                  ks = mesh_o%ks, ke = mesh_o%ke)

    call sbp_interp_v2w%apply(u_w, work, mesh_o%ny, 'y', u)
    call sbp_interp_v2w%apply(v_w, work, mesh_o%nx, 'x', v)

end subroutine interp_v2w_tile

end module interpolator_v2w_mod
