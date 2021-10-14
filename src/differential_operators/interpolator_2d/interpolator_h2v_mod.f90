module interpolator_h2v_mod

use grid_field_mod,         only : grid_field_t, tile_field_t
use domain_mod,             only : domain_t
use exchange_abstract_mod,  only : exchange_t
use mesh_mod,               only : tile_mesh_t
use sbp_operator_mod,       only : sbp_operator_t

implicit none

type, public :: interpolator_h2v_t
    class(exchange_t),     allocatable :: exchange
    class(exchange_t),     allocatable :: exchange_Ah
    class(sbp_operator_t), allocatable :: sbp_interp_h2v
contains
    procedure, public :: interp_h2v
    procedure, public :: interp_h2v_Ah_C
end type interpolator_h2v_t

contains

subroutine interp_h2v(this, u, v, u_h, v_h, domain)

    class(interpolator_h2v_t), intent(inout) :: this
    type(grid_field_t),        intent(inout) :: u_h, v_h, u, v
    type(domain_t),            intent(in)    :: domain

    integer(kind=4) :: t

    call this%exchange%do_vec(u_h, v_h, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_h2v_tile(u%tile(t), v%tile(t), u_h%tile(t), v_h%tile(t), &
                    this%sbp_interp_h2v, domain%mesh_u%tile(t), &
                    domain%mesh_v%tile(t), domain%mesh_p%tile(t))
    end do

end subroutine interp_h2v
subroutine interp_h2v_Ah_C(this, u, v, u_h, v_h, domain)

    class(interpolator_h2v_t), intent(inout) :: this
    type(grid_field_t),        intent(inout) :: u_h, v_h, u, v
    type(domain_t),            intent(in)    :: domain

    integer(kind=4) :: t

    call this%exchange_Ah%do_vec(u_h, v_h, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call interp_h2v_tile(u%tile(t), v%tile(t), u_h%tile(t), v_h%tile(t), &
                    this%sbp_interp_h2v, domain%mesh_y%tile(t), &
                    domain%mesh_x%tile(t), domain%mesh_p%tile(t))
    end do

end subroutine interp_h2v_Ah_C

subroutine interp_h2v_tile(u, v, u_h, v_h, sbp_interp_h2v, &
                           mesh_x, mesh_y, mesh_o)

    use tile_mod, only : tile_t

    type(tile_field_t),     intent(inout) :: u_h, v_h
    type(tile_field_t),     intent(inout) :: u, v
    class(sbp_operator_t),  intent(in)    :: sbp_interp_h2v
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_o

    type(tile_t) :: work

    work = tile_t(is = mesh_x%is, ie = mesh_x%ie, &
                  js = mesh_x%js, je = mesh_x%je, &
                  ks = mesh_x%ks, ke = mesh_x%ke)
    call sbp_interp_h2v%apply(u, work, mesh_x%nx, 'x', u_h)

    work = tile_t(is = mesh_y%is, ie = mesh_y%ie, &
                  js = mesh_y%js, je = mesh_y%je, &
                  ks = mesh_y%ks, ke = mesh_y%ke)
    call sbp_interp_h2v%apply(v, work, mesh_y%ny, 'y', v_h)

end subroutine interp_h2v_tile

end module interpolator_h2v_mod
