module div_c2_mod

use abstract_div_mod,   only : div_operator_t
use grid_function_mod,  only : grid_function_t
use mesh_mod,           only : mesh_t
use partition_mod,      only : partition_t
!use timer_parallel_mod,        only : t_bstab_matvec
implicit none

type, public, extends(div_operator_t) :: div_c2_t
contains
    procedure, public :: calc_div => calc_div_c2
end type div_c2_t

contains

subroutine calc_div_c2(this, u, v, div, mesh, partition, multiplier)

    class(div_c2_t),        intent(inout) :: this
    type(partition_t),      intent(in)    :: partition
    type(grid_function_t),  intent(inout) :: u(partition%ts:partition%te), &
                                             v(partition%ts:partition%te)
    type(grid_function_t),  intent(inout) :: div(partition%ts:partition%te)
    type(mesh_t),           intent(in)    :: mesh(partition%ts:partition%te)
    real(kind=8), optional, intent(in)    :: multiplier

    integer(kind=4) :: i, j, k, t
    real(kind=8) mult_loc, hx

    do t = partition%ts, partition%te
        call calc_div_on_tile(u(t), v(t), div(t), mesh(t), multiplier)
    end do


end subroutine calc_div_c2

subroutine calc_div_on_tile(u, v, div, mesh, multiplier)

    type(grid_function_t),  intent(in)    :: u, v
    type(grid_function_t),  intent(inout) :: div
    type(mesh_t),           intent(in)    :: mesh
    real(kind=8), optional, intent(in)    :: multiplier

    real(kind=8) :: hx, mult_loc
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k

    mult_loc = 1.0_8
    if(present(multiplier)) mult_loc = multiplier

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx
    do k = ks, ke
        do j = js, je
            do i = is, ie
                div%p(i,j,k) = (mesh%Gu(i,j)*u%p(i,j,k)-mesh%Gu(i-1,j  )*u%p(i-1,j  ,k) +  &
                                mesh%Gv(i,j)*v%p(i,j,k)-mesh%Gv(i  ,j-1)*v%p(i  ,j-1,k))/  &
                                (mesh%G(i,j)*hx)*mult_loc
            end do
        end do
    end do

end subroutine calc_div_on_tile

end module div_c2_mod
