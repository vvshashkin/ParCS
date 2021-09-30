module ke_unstaggered_mod

use abstract_KE_mod, only : KE_operator_t
use grid_field_mod,  only : grid_field_t, tile_field_t
use domain_mod,      only : domain_t
use mesh_mod,        only : tile_mesh_t

implicit none

type, public, extends(KE_operator_t) :: ke_unstaggered_t
contains
    procedure, public :: calc_KE
end type ke_unstaggered_t

contains

subroutine calc_KE(this, KE, u, v, domain)
    class(ke_unstaggered_t), intent(inout) :: this
    type(domain_t),          intent(in)    :: domain
    type(grid_field_t),      intent(inout) :: u, v
    type(grid_field_t),      intent(inout) :: KE

    integer(kind=4) :: t

    do t = domain%mesh_p%ts, domain%mesh_p%ts
        call calc_KE_on_tile(KE%tile(t), u%tile(t), v%tile(t), domain%mesh_p%tile(t))
    end do

end subroutine calc_KE

subroutine calc_KE_on_tile(KE, u, v, mesh_p)

    type(tile_field_t), intent(in)    :: u, v
    type(tile_field_t), intent(inout) :: KE
    type(tile_mesh_t),  intent(in)    :: mesh_p

    integer(kind=4) :: i, j, k
    real(kind=8)    :: u_contra, v_contra

    do k = mesh_p%ks, mesh_p%ke
        do j = mesh_p%js, mesh_p%je
            do i = mesh_p%is, mesh_p%ie
                !transform to contravariant components
                u_contra = mesh_p%Qi(1,i,j)*u%p(i,j,k)+mesh_p%Qi(2,i,j)*v%p(i,j,k)
                v_contra = mesh_p%Qi(2,i,j)*u%p(i,j,k)+mesh_p%Qi(3,i,j)*v%p(i,j,k)
                KE%p = u%p(i,j,k)*u_contra + v%p(i,j,k)*v_contra
            end do
        end do
    end do
end subroutine calc_KE_on_tile

end module ke_unstaggered_mod
