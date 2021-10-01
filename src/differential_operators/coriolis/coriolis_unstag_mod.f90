module coriolis_unstag_mod

use grid_field_mod,        only : grid_field_t, tile_field_t
use domain_mod,            only : domain_t
use abstract_coriolis_mod, only : coriolis_operator_t

implicit none

type, public, extends(coriolis_operator_t) :: coriolis_unstag_t
    type(grid_field_t) :: f !coriolis parameter
contains
    procedure, public :: calc_coriolis
end type coriolis_unstag_t

contains

subroutine calc_coriolis(this, cor_u, cor_v, ut, vt, domain)
    class(coriolis_unstag_t), intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),     intent(inout) :: cor_u, cor_v

    integer(kind=4) :: t

    do t = domain%partition%ts, domain%partition%te
        call calc_coriolis_on_tile(ut%tile(t), vt%tile(t), this%f%tile(t), &
                                   cor_u%tile(t), cor_v%tile(t), &
                                   domain%mesh_u%tile(t), domain%mesh_v%tile(t))
    end do


end subroutine calc_coriolis

subroutine calc_coriolis_on_tile(ut, vt, f, cor_u, cor_v, mesh_u, mesh_v)

    !coriolis(\vec{u}) = f*\vec{u}^\perp

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t), intent(in)    :: ut, vt, f
    type(tile_field_t), intent(inout) :: cor_u, cor_v
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k
    real(kind=8)    :: u_covariant, v_covariant

!This implementation works only for unstaggered case, so mesh_u==mesh_v==mesh_p
    do k = mesh_u%ks, mesh_u%ke
        do j = mesh_u%js, mesh_u%je
            do i = mesh_u%is, mesh_u%ie
                cor_u%p(i,j,k) =  f%p(i,j,1)*vt%p(i,j,k)*mesh_u%G(i,j)
                cor_v%p(i,j,k) = -f%p(i,j,1)*ut%p(i,j,k)*mesh_u%G(i,j)
            end do
        end do
    end do


end subroutine calc_coriolis_on_tile

end module coriolis_unstag_mod
