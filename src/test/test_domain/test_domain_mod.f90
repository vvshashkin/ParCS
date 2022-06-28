module test_domain_mod

use domain_mod,             only : domain_t
use grid_field_mod,         only : grid_field_t
use config_domain_mod,      only : config_domain_t
use domain_factory_mod,     only : create_domain
use grid_field_factory_mod, only : create_grid_field

implicit none

contains

subroutine test_domain_simple()

    type(domain_t) :: domain
    type(grid_field_t) :: f, f2

    integer(kind=4)  :: nh=32, nz=10
    character(len=1) :: hor_grid_type = 'C'

    call create_domain(domain, "cube", hor_grid_type, nh, nz)

    call create_grid_field(f , 2, 0, domain%mesh_p)
    call create_grid_field(f2, 2, 0, domain%mesh_p)

    call  f%assign(1.0_8, domain%mesh_p)
    call f2%assign(1.0_8, domain%mesh_p)

    call f%update(1.0_8, f2, domain%mesh_p)

    print*, domain%mesh_u%tile(domain%mesh_u%ts)%ie

    print *, "domain test passed"

end subroutine test_domain_simple

subroutine test_domain_config(horizontal_staggering, vertical_staggering)

    character(len=*), intent(in) :: horizontal_staggering, vertical_staggering

    type(domain_t)        :: domain
    type(config_domain_t) :: config_domain
    type(grid_field_t)    :: f, f2

    config_domain%N  = 32
    config_domain%Nz = 10
    config_domain%staggering_type     = horizontal_staggering
    config_domain%vertical_staggering = vertical_staggering
    config_domain%metric_type         = "shallow_atmosphere_metric"
    config_domain%topology_type       = "cube"
    config_domain%h_top               = 1.0_8
    config_domain%is_orographic_curvilinear = .true.
    config_domain%orography_name = "test_orography"
    call config_domain%config_metric%set_defaults()

    call create_domain(domain, config_domain)

    call create_grid_field(f , 2, 0, domain%mesh_p)
    call create_grid_field(f2, 2, 0, domain%mesh_p)

    call  f%assign(1.0_8, domain%mesh_p)
    call f2%assign(1.0_8, domain%mesh_p)

    call f%update(1.0_8, f2, domain%mesh_p)

    print*, domain%mesh_u%tile(domain%mesh_u%ts)%ie

    print *, "domain test passed"

end subroutine test_domain_config

end module test_domain_mod
