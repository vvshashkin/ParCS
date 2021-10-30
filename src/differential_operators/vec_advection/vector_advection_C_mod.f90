module vector_advection_C_mod

use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use halo_mod,                      only : halo_vec_t
use interpolator_v2h_mod,          only : interpolator_v2h_t
use interpolator_h2v_mod,          only : interpolator_h2v_t
use abstract_vector_advection_mod, only : vector_advection_operator_t
use parcomm_mod,                   only : parcomm_global

implicit none

type, public, extends(vector_advection_operator_t) :: vector_advection_C_t
    type(grid_field_t) :: u_at_v, v_at_u, uh, vh
    type(interpolator_v2h_t) :: interp_v2h_op
    type(interpolator_h2v_t) :: interp_h2v_op
    class(halo_vec_t), allocatable :: halo_uv, tendency_edge_sync
contains
    procedure :: calc_vec_advection
    procedure :: calc_vec_advection_contra
end type vector_advection_C_t

contains

subroutine calc_vec_advection(this, u_tend, v_tend, u, v, ut, vt, domain)
    class(vector_advection_C_t),  intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u,  v!covariant components
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    call parcomm_global%abort("vector_advection_c_t, calc_vec_advection not implemented")

end subroutine calc_vec_advection

subroutine calc_vec_advection_contra(this, u_tend, v_tend, ut, vt, domain)
    class(vector_advection_C_t),  intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    integer(kind=4) :: t

    call this%interp_v2h_op%interp_v2h(this%uh, this%vh, ut, vt, domain)
    call this%interp_h2v_op%interp_h2v(this%v_at_u, this%u_at_v, this%vh, this%uh, domain)

    call this%halo_uv%get_halo_vector(ut, vt, domain, 3)

    do t = domain%mesh_o%ts, domain%mesh_o%te
        call calc_advection_1comp_contra_tile(u_tend%tile(t), ut%tile(t), ut%tile(t), this%v_at_u%tile(t), &
                                              domain%mesh_u%scale, domain%mesh_u%tile(t),"u")
        call calc_advection_1comp_contra_tile(v_tend%tile(t), vt%tile(t), this%u_at_v%tile(t), vt%tile(t), &
                                              domain%mesh_v%scale, domain%mesh_v%tile(t),"v")
    end do

    call this%tendency_edge_sync%get_halo_vector(u_tend, v_tend, domain, 1)
end subroutine calc_vec_advection_contra

subroutine calc_advection_1comp_contra_tile(uvt_tend, uvt, ut, vt, scale, mesh, component)

    use mesh_mod, only : tile_mesh_t
    use tile_mod, only : tile_t

    type(tile_field_t),     intent(in)    :: uvt, ut, vt
    real(kind=8),           intent(in)    :: scale
    type(tile_mesh_t),      intent(in)    :: mesh
    character(len=*),       intent(in)    :: component

    type(tile_field_t),     intent(inout) :: uvt_tend

    integer(kind=4) :: i, j, k
    integer(kind=4) :: component_num
    integer(kind=4) :: is, ie, js, je, ks, ke
    real(kind=8)    :: hx
    real(kind=8)    :: dx, dy, zl, zr

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx

    if(component=="u") then
        component_num = 1
    else if(component=="v") then
        component_num = 2
    else
        call parcomm_global%abort("vector_advection_C_mod, unknown vector component type:"//component)
    end if

    do k = ks, ke
        do j = js, je
            do i = is, ie
                zl = .5_8+sign(.5_8,ut%p(i,j,k))
                zr = 1._8-zl
                dx = (zl*(uvt%p(i,j,k)-uvt%p(i-1,j,k))+zr*(uvt%p(i+1,j,k)-uvt%p(i,j,k)))
                zl = .5_8+sign(.5_8,vt%p(i,j,k))
                zr = 1._8-zl
                dy = (zl*(uvt%p(i,j,k)-uvt%p(i,j-1,k))+zr*(uvt%p(i,j+1,k)-uvt%p(i,j,k)))
                uvt_tend%p(i,j,k) =-(ut%p(i,j,k)*dx+vt%p(i,j,k)*dy) / (hx*scale)-&
                                    (ut%p(i,j,k)*ut%p(i,j,k)*mesh%T(1,1,component_num,i,j)+ &
                                     2.0_8*ut%p(i,j,k)*vt%p(i,j,k)*mesh%T(1,2,component_num,i,j)+ &
                                     vt%p(i,j,k)*vt%p(i,j,k)*mesh%T(2,2,component_num,i,j)) / scale
            end do
        end do
    end do

end subroutine calc_advection_1comp_contra_tile

end module vector_advection_C_mod
