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
    class(halo_vec_t),     allocatable :: sync_edges_cov
    class(halo_vec_t),     allocatable :: sync_edges_contra
    class(exchange_t),     allocatable :: exch_uv_interior
contains
    procedure :: calc_vec_advection
    procedure :: calc_vec_advection_contra
end type vector_advection_Ah_t

contains

subroutine calc_vec_advection(this, u_tend, v_tend, u, v, ut, vt, domain)
    class(vector_advection_Ah_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u,  v!covariant components
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    integer(kind=4) :: t

    do t = domain%mesh_xy%ts, domain%mesh_xy%te
        call calc_vec_advection_tile(u_tend%tile(t), v_tend%tile(t),               &
                                     u%tile(t), v%tile(t), ut%tile(t), vt%tile(t), &
                                     domain%mesh_xy%tile(t), domain%mesh_xy%scale, &
                                     this%sbp_op)
    end do

    call this%sync_edges_cov%get_halo_vector(u_tend,v_tend,domain,0)
end subroutine calc_vec_advection

subroutine calc_vec_advection_tile(u_tend, v_tend, u, v, ut, vt, mesh, scale, sbp_op)

    use mesh_mod, only : tile_mesh_t
    use tile_mod, only : tile_t

    type(tile_field_t),     intent(in)    :: u, v, ut, vt
    type(tile_mesh_t),      intent(in)    :: mesh
    real(kind=8),           intent(in)    :: scale
    class(sbp_operator_t),  intent(in)    :: sbp_op

    type(tile_field_t),     intent(inout) :: u_tend, v_tend

    real(kind=8)    :: Dx(mesh%is:mesh%ie,mesh%js:mesh%je,1)
    real(kind=8)    :: Dy(mesh%is:mesh%ie,mesh%js:mesh%je,1)
    type(tile_t)    :: dxdy_tile
    integer(kind=4) :: i, j, k
    integer(kind=4) :: is, ie, js, je, ks, ke
    real(kind=8)    :: hx

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    dxdy_tile = tile_t(is = is, ie=ie, js=js, je=je,ks = 1, ke=1)

    hx = mesh%hx

    do k = ks, ke

        dxdy_tile%ks = k; dxdy_tile%ke = k

        call sbp_op%apply(Dx, dxdy_tile, dxdy_tile, mesh%nx, 'x', u)
        call sbp_op%apply(Dy, dxdy_tile, dxdy_tile, mesh%ny, 'y', u)

        do j = js, je
            do i = is, ie
                u_tend%p(i,j,k) =-(ut%p(i,j,k)*Dx(i,j,1)+vt%p(i,j,k)*Dy(i,j,1)) / (hx*scale)+&
                                  (u%p(i,j,k)*ut%p(i,j,k)*mesh%T(1,1,1,i,j)+ &
                                   v%p(i,j,k)*ut%p(i,j,k)*mesh%T(1,1,2,i,j)+ &
                                   u%p(i,j,k)*vt%p(i,j,k)*mesh%T(1,2,1,i,j)+ &
                                   v%p(i,j,k)*vt%p(i,j,k)*mesh%T(1,2,2,i,j)) / scale
            end do
        end do

        call sbp_op%apply(Dx, dxdy_tile, dxdy_tile, mesh%nx, 'x', v)
        call sbp_op%apply(Dy, dxdy_tile, dxdy_tile, mesh%ny, 'y', v)

        do j = js, je
            do i = is, ie
                v_tend%p(i,j,k) =-(ut%p(i,j,k)*Dx(i,j,1)+vt%p(i,j,k)*Dy(i,j,1)) / (hx*scale)+&
                                  (u%p(i,j,k)*ut%p(i,j,k)*mesh%T(2,1,1,i,j)+ &
                                   v%p(i,j,k)*ut%p(i,j,k)*mesh%T(2,1,2,i,j)+ &
                                   u%p(i,j,k)*vt%p(i,j,k)*mesh%T(2,2,1,i,j)+ &
                                   v%p(i,j,k)*vt%p(i,j,k)*mesh%T(2,2,2,i,j)) / scale
            end do
        end do

    end do

end subroutine calc_vec_advection_tile

subroutine calc_vec_advection_contra(this, u_tend, v_tend, ut, vt, domain)
    class(vector_advection_Ah_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    integer(kind=4) :: t

    do t = domain%mesh_xy%ts, domain%mesh_xy%te
        call calc_vec_advection_contra_tile(u_tend%tile(t), v_tend%tile(t),               &
                                            ut%tile(t), vt%tile(t),                       &
                                            domain%mesh_xy%tile(t), domain%mesh_xy%scale, &
                                            this%sbp_op)
    end do

    call this%sync_edges_contra%get_halo_vector(u_tend,v_tend,domain,0)
end subroutine calc_vec_advection_contra

subroutine calc_vec_advection_contra_tile(u_tend, v_tend, ut, vt, mesh, scale, sbp_op)

    use mesh_mod, only : tile_mesh_t
    use tile_mod, only : tile_t

    type(tile_field_t),     intent(in)    :: ut, vt
    type(tile_mesh_t),      intent(in)    :: mesh
    real(kind=8),           intent(in)    :: scale
    class(sbp_operator_t),  intent(in)    :: sbp_op

    type(tile_field_t),     intent(inout) :: u_tend, v_tend

    real(kind=8)    :: Dx(mesh%is:mesh%ie,mesh%js:mesh%je,1)
    real(kind=8)    :: Dy(mesh%is:mesh%ie,mesh%js:mesh%je,1)
    type(tile_t)    :: dxdy_tile
    integer(kind=4) :: i, j, k
    integer(kind=4) :: is, ie, js, je, ks, ke
    real(kind=8)    :: hx

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    dxdy_tile = tile_t(is = is, ie=ie, js=js, je=je,ks = 1, ke=1)

    hx = mesh%hx

    do k = ks, ke

        dxdy_tile%ks = k; dxdy_tile%ke = k

        call sbp_op%apply(Dx, dxdy_tile, dxdy_tile, mesh%nx, 'x', ut)
        call sbp_op%apply(Dy, dxdy_tile, dxdy_tile, mesh%ny, 'y', ut)

        do j = js, je
            do i = is, ie
                u_tend%p(i,j,k) =-(ut%p(i,j,k)*Dx(i,j,1)+vt%p(i,j,k)*Dy(i,j,1)) / (hx*scale)-&
                                  (ut%p(i,j,k)*ut%p(i,j,k)*mesh%T(1,1,1,i,j)+ &
                                   2.0_8*ut%p(i,j,k)*vt%p(i,j,k)*mesh%T(1,2,1,i,j)+ &
                                   vt%p(i,j,k)*vt%p(i,j,k)*mesh%T(2,2,1,i,j)) / scale
            end do
        end do

        call sbp_op%apply(Dx, dxdy_tile, dxdy_tile, mesh%nx, 'x', vt)
        call sbp_op%apply(Dy, dxdy_tile, dxdy_tile, mesh%ny, 'y', vt)

        do j = js, je
            do i = is, ie
                v_tend%p(i,j,k) =-(ut%p(i,j,k)*Dx(i,j,1)+vt%p(i,j,k)*Dy(i,j,1)) / (hx*scale)-&
                                  (ut%p(i,j,k)*ut%p(i,j,k)*mesh%T(1,1,2,i,j)+ &
                                   2.0_8*vt%p(i,j,k)*ut%p(i,j,k)*mesh%T(1,2,2,i,j)+ &
                                   vt%p(i,j,k)*vt%p(i,j,k)*mesh%T(2,2,2,i,j)) / scale
            end do
        end do

    end do

end subroutine calc_vec_advection_contra_tile

end module vector_advection_Ah_mod
