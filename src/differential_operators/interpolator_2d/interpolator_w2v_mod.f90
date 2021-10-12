module interpolator_w2v_mod

use grid_field_mod,         only : grid_field_t, tile_field_t
use domain_mod,             only : domain_t
use exchange_abstract_mod,  only : exchange_t
use mesh_mod,               only : tile_mesh_t
use sbp_operator_mod,       only : sbp_operator_t

implicit none

type, public :: interpolator_w2v_t
    class(exchange_t),     allocatable :: exchange
    class(sbp_operator_t), allocatable :: sbp_interp_w2v
contains
    procedure, public :: interp_w2v
    procedure, public :: interp_w2v_Ah_C
end type interpolator_w2v_t

contains

subroutine interp_w2v(this, u, v, u_w, v_w, domain)

    class(interpolator_w2v_t), intent(inout) :: this
    type(grid_field_t),        intent(inout) :: u_w, v_w, u, v
    type(domain_t),            intent(in)    :: domain

    integer(kind=4) :: t

end subroutine interp_w2v
subroutine interp_w2v_Ah_C(this, u, v, u_w, v_w, domain)

    class(interpolator_w2v_t), intent(inout) :: this
    type(grid_field_t),        intent(inout) :: u_w, v_w, u, v
    type(domain_t),            intent(in)    :: domain

    integer(kind=4) :: t

    call this%exchange%do_vec(u_w, v_w, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_w2v_tile(u%tile(t), v%tile(t), u_w%tile(t), v_w%tile(t), &
                    this%sbp_interp_w2v, domain%mesh_y%tile(t), &
                    domain%mesh_x%tile(t))
    end do

end subroutine interp_w2v_Ah_C

subroutine interp_w2v_tile(u, v, u_w, v_w, sbp_interp_w2v, &
                           mesh_u, mesh_v)

    use tile_mod, only : tile_t

    type(tile_field_t),     intent(inout) :: u_w, v_w
    type(tile_field_t),     intent(inout) :: u, v
    class(sbp_operator_t),  intent(in)    :: sbp_interp_w2v
    type(tile_mesh_t),      intent(in)    :: mesh_u, mesh_v

    type(tile_t) :: work

    work = tile_t(is = mesh_u%is, ie = mesh_u%ie, &
                  js = mesh_u%js, je = mesh_u%je, &
                  ks = mesh_u%ks, ke = mesh_u%ke)
    print*, mesh_u%ny
    call sbp_interp_w2v%apply(u, work, mesh_u%ny, 'y', u_w)

    work = tile_t(is = mesh_v%is, ie = mesh_v%ie, &
                  js = mesh_v%js, je = mesh_v%je, &
                  ks = mesh_v%ks, ke = mesh_v%ke)
    print*, mesh_v%nx
    call sbp_interp_w2v%apply(v, work, mesh_v%nx, 'x', v_w)

end subroutine interp_w2v_tile

end module interpolator_w2v_mod
