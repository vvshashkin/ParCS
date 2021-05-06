module test_domain_mod

use domain_mod,     only : domain_t
use grid_field_mod, only : grid_field_t

implicit none

contains

subroutine test_domain()

    use domain_factory_mod,     only : create_ecs_global_domain
    use grid_field_factory_mod, only : create_grid_field

    type(domain_t) :: domain
    type(grid_field_t) :: f, f2

    integer(kind=4)  :: nh=100, nz=10, halo_width=50
    character(len=1) :: hor_grid_type = 'C'

    call create_ecs_global_domain(domain, hor_grid_type, nh, nz)

    call create_grid_field(f , 2, 0, domain%mesh_p)
    call create_grid_field(f2, 2, 0, domain%mesh_p)

    call  f%assign(1.0_8, domain%mesh_p)
    call f2%assign(1.0_8, domain%mesh_p)

    call f%update(f2, 1.0_8, domain%mesh_p)

    print*, domain%mesh_u%tile(domain%mesh_u%ts)%ie

end subroutine test_domain

end module test_domain_mod
