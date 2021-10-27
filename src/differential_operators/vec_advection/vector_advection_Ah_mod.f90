module vector_advection_Ah_mod

use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use sbp_operator_mod,              only : sbp_operator_t
use halo_mod,                      only : halo_vec_t
use exchange_abstract_mod,         only : exchange_t
use abstract_vector_advection_mod, only : vector_advection_operator_t

implicit none

type, public, extends(vector_advection_operator_t) :: vector_advection_Ah_t
    class(sbp_operator_t), allocatable :: sbp_op
    class(halo_vec_t),     allocatable :: sync_edges
    class(exchange_t),     allocatable :: exch_uv_interior
contains
    procedure :: calc_vec_advection
end type vector_advection_Ah_t

contains

subroutine calc_vec_advection(this, u_tend, v_tend, u, v, ut, vt, domain)
    class(vector_advection_Ah_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u,  v!covariant components
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    integer(kind=4) :: t

    !call u_tend%assign(0.0_8,domain%mesh_u)
    !call v_tend%assign(0.0_8,domain%mesh_v)

    do t = domain%mesh_xy%ts, domain%mesh_xy%te
        call calc_vec_advection_tile(u_tend%tile(t), v_tend%tile(t),               &
                                     u%tile(t), v%tile(t), ut%tile(t), vt%tile(t), &
                                     domain%mesh_xy%tile(t), domain%mesh_xy%scale, &
                                     this%sbp_op)
    end do
end subroutine calc_vec_advection

subroutine calc_vec_advection_tile(u_tend, v_tend, u, v, ut, vt, mesh, scale, sbp_op)

    use mesh_mod, only : tile_mesh_t
    use tile_mod, only : tile_t

    type(tile_field_t),     intent(in)    :: u, v, ut, vt
    type(tile_mesh_t),      intent(in)    :: mesh
    real(kind=8),           intent(in)    :: scale
    class(sbp_operator_t),  intent(in)    :: sbp_op

    type(tile_field_t),     intent(inout) :: u_tend, v_tend

    u_tend%p = 0.0_8
    v_tend%p = 0.0_8

end subroutine calc_vec_advection_tile

end module vector_advection_Ah_mod
