module coriolis_factory_mod

use domain_mod,             only : domain_t
use abstract_coriolis_mod,  only : coriolis_operator_t
use parcomm_mod,            only : parcomm_global
use grid_field_mod,         only : grid_field_t

implicit none

private

public :: create_coriolis

contains

subroutine create_coriolis(coriolis_op, coriolis_op_name, domain)

    class(coriolis_operator_t), allocatable, intent(out) :: coriolis_op
    character(len=*),                        intent(in)  :: coriolis_op_name
    type(domain_t),                          intent(in)  :: domain

    select case(coriolis_op_name)

    case("coriolis_A_Ah")
        call create_coriolis_unstaggered(coriolis_op, domain)
    case default
        call parcomm_global%abort("Unknown coriolis operator: "//coriolis_op_name)
    end select

end subroutine create_coriolis

subroutine create_coriolis_unstaggered(coriolis_op, domain)

    use coriolis_unstag_mod,    only : coriolis_unstag_t
    use grid_field_factory_mod, only : create_grid_field

    class(coriolis_operator_t), allocatable, intent(out) :: coriolis_op
    type(domain_t),                          intent(in)  :: domain

    type(coriolis_unstag_t), allocatable :: coriolis_unstag

    allocate(coriolis_unstag)

    call create_grid_field(coriolis_unstag%f, 0, 0, domain%mesh_u)

    call calc_coriolis_parameter(coriolis_unstag%f, domain%mesh_p)

    call move_alloc(coriolis_unstag, coriolis_op)

end subroutine create_coriolis_unstaggered

subroutine calc_coriolis_parameter(f, mesh)

    use mesh_mod, only : mesh_t
    use sph_coords_mod, only : cart2sph

    type(grid_field_t), intent(inout) :: f
    type(mesh_t),       intent(in)  :: mesh

    integer(kind=4) :: t, i, j, k
    real(kind=8) :: pr, n(3), lam, phi

    do t = mesh%ts, mesh%te
        do j = mesh%tile(t)%js, mesh%tile(t)%je
            do i = mesh%tile(t)%is, mesh%tile(t)%ie
                n = [mesh%tile(t)%rx(i,j),mesh%tile(t)%ry(i,j),mesh%tile(t)%rz(i,j)]
                pr = dot_product(mesh%rotation_axis, n) / sqrt(dot_product(n,n))
                f%tile(t)%p(i,j,1) = 2*mesh%omega*pr
            end do
        end do
    end do

end subroutine calc_coriolis_parameter

end module coriolis_factory_mod
