module co2contra_Cgrid_mod

use abstract_co2contra_mod, only : co2contra_operator_t
use grid_field_mod,         only : grid_field_t
use domain_mod,             only : domain_t
use exchange_abstract_mod,  only : exchange_t

implicit none

type, extends(co2contra_operator_t), public :: co2contra_c_sbp21_t
    class(exchange_t), allocatable :: exchange_inner
    contains
        procedure :: transform => transform_co2contra_c_sbp21
end type co2contra_c_sbp21_t

contains

subroutine transform_co2contra_c_sbp21(this, u_contra, v_contra, u_cov, v_cov, domain)
    class(co2contra_c_sbp21_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u_cov, v_cov
    !output:
    type(grid_field_t),           intent(inout) :: u_contra, v_contra

    integer(kind=4) :: t

    call this%exchange_inner%do_vec(u_cov, v_cov, domain%parcomm)

    do t=domain%mesh_o%ts, domain%mesh_o%te
        call transform_co2contra_c_sbp21_tile(u_contra%tile(t), v_contra%tile(t),&
                                              u_cov%tile(t), v_cov%tile(t),      &
                                              domain%mesh_u%tile(t), domain%mesh_v%tile(t), &
                                              domain%mesh_o%tile(t))
    end do

end subroutine transform_co2contra_c_sbp21

subroutine transform_co2contra_c_sbp21_tile(u_contra, v_contra, u_cov, v_cov, mesh_u, mesh_v, mesh_o)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(in)    :: u_cov, v_cov
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v, mesh_o
    !output
    type(tile_field_t), intent(inout) :: u_contra, v_contra

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke
    real(kind=8)    :: u_at_o(mesh_u%is:mesh_u%ie,mesh_u%js-1:mesh_u%je)
    real(kind=8)    :: v_at_o(mesh_u%is-1:mesh_u%ie,mesh_u%js:mesh_u%je)
    real(kind=8)    :: v_at_u, u_at_v

    ks = mesh_o%ks; ke = mesh_o%ke

    do k=ks, ke

        is = mesh_u%is; ie = mesh_u%ie
        js = mesh_u%js; je = mesh_u%je

        do j=js,je
            do i=max(is-1,1),min(ie,mesh_o%nx)
                v_at_o(i,j) = 0.5_8*(v_cov%p(i,j,k)+v_cov%p(i,j+1,k))*mesh_o%G(i,j)*mesh_o%Qi(2,i,j)
            end do

            if(is == 1) then
                u_contra%p(1,j,k) = mesh_u%Qi(1,1,j)*u_cov%p(1,j,k)+v_at_o(1,j) / mesh_u%G(1,j)
            end if
            do i = max(is,2), min(ie,mesh_o%nx)
                v_at_u = 0.5_8*(v_at_o(i,j)+v_at_o(i-1,j))
                u_contra%p(i,j,k) = mesh_u%Qi(1,i,j)*u_cov%p(i,j,k)+v_at_u / mesh_u%G(i,j)
            end do
            if(ie == mesh_o%nx+1) then
                u_contra%p(ie,j,k) = mesh_u%Qi(1,ie,j)*u_cov%p(ie,j,k)+v_at_o(ie-1,j) / mesh_u%G(ie,j)
            end if
        end do

        is = mesh_v%is; ie = mesh_v%ie
        js = mesh_v%js; je = mesh_v%je

        do j=max(js-1,1),min(je,mesh_o%ny)
            do i = is,ie
                u_at_o(i,j) = 0.5_8*(u_cov%p(i+1,j,k)+u_cov%p(i,j,k))*mesh_o%G(i,j)*mesh_o%Qi(2,i,j)
            end do
        end do

        if(js == 1) then
            do i=is, ie
                v_contra%p(i,1,k) = mesh_v%Qi(3,i,1)*v_cov%p(i,1,k)+u_at_o(i,1) / mesh_v%G(i,1)
            end do
        end if
        do j = max(js,2), min(je,mesh_o%ny)
            do i = is, ie
                u_at_v = 0.5_8*(u_at_o(i,j)+u_at_o(i,j-1))
                v_contra%p(i,j,k) = mesh_v%Qi(3,i,j)*v_cov%p(i,j,k)+u_at_v / mesh_v%G(i,j)
            end do
        end do
        if(je == mesh_o%ny+1) then
            do i = is, ie
                v_contra%p(i,je,k) = mesh_v%Qi(3,i,je)*v_cov%p(i,je,k)+u_at_o(i,je-1) / mesh_v%G(i,je)
            end do
        end if

    end do

end subroutine transform_co2contra_c_sbp21_tile

end module co2contra_Cgrid_mod
