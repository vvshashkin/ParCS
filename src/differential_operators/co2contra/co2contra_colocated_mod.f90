module co2contra_colocated_mod

use abstract_co2contra_mod, only : co2contra_operator_t
use grid_field_mod,         only : grid_field_t
use domain_mod,             only : domain_t

implicit none

type, extends(co2contra_operator_t), public :: co2contra_colocated_t

contains

procedure :: transform => transform_co2contra_colocated

end type co2contra_colocated_t

contains

subroutine transform_co2contra_colocated(this, u_contra, v_contra, u_cov, v_cov, domain)
    class(co2contra_colocated_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u_cov, v_cov
    !output:
    type(grid_field_t),           intent(inout) :: u_contra, v_contra

    integer(kind=4) :: t

    do t=domain%mesh_u%ts, domain%mesh_u%te
        call transform_co2contra_colocated_tile(u_contra%tile(t), v_contra%tile(t),&
                                                u_cov%tile(t), v_cov%tile(t), domain%mesh_u%tile(t))
    end do

end subroutine transform_co2contra_colocated

subroutine transform_co2contra_colocated_tile(u_contra, v_contra, u_cov, v_cov, mesh)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(in)    :: u_cov, v_cov
    type(tile_mesh_t),  intent(in)    :: mesh
    !output
    type(tile_field_t), intent(inout) :: u_contra, v_contra

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks, ke
        do j=js, je
            do i=is,ie
                u_contra%p(i,j,k) = mesh%Qi(1,i,j)*u_cov%p(i,j,k)+mesh%Qi(2,i,j)*v_cov%p(i,j,k)
                v_contra%p(i,j,k) = mesh%Qi(2,i,j)*u_cov%p(i,j,k)+mesh%Qi(3,i,j)*v_cov%p(i,j,k)
            end do
        end do
    end do

end subroutine transform_co2contra_colocated_tile

end module co2contra_colocated_mod
