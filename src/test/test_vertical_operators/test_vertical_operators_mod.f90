module test_vertical_operators_mod

use test_fields_3d_mod,       only: scalar_field3d_t
use vertical_test_field_mod,  only: vertical_ExnerP_t
use const_N_profile_mod,      only: const_N_profile_t
use domain_mod,               only: domain_t
use domain_factory_mod,       only: create_domain
use config_domain_mod,        only: config_domain_t
use grid_field_mod,           only: grid_field_t
use grid_field_factory_mod,   only: create_grid_field

implicit none

contains

subroutine test_vertical_gradient_operator()
    integer, parameter :: nh = 8, nz = 10
    real(kind=8), parameter :: h_top = 30e3_8

    class(scalar_field3d_t), allocatable :: scalar_gen
    type(config_domain_t) :: config_domain
    type(domain_t)     :: domain
    type(grid_field_t) :: p, pz_true, pz

    scalar_gen = vertical_ExnerP_t(t0=300.0_8, p0=1e5_8, &
                                   vert_profile=const_N_profile_t(N=0.01))

    config_domain%N  = nh
    config_domain%Nz = nz
    config_domain%staggering_type     = "C"
    config_domain%vertical_staggering = "CharneyPhilips"
    config_domain%metric_type         = "shallow_atmosphere_metric"
    config_domain%topology_type       = "cube"
    config_domain%h_top = h_top
    call config_domain%config_metric%set_defaults()
    config_domain%config_metric%vertical_scale = h_top

    call create_domain(domain, config_domain)

    call create_grid_field(p,0,0,domain%mesh_p)
    call create_grid_field(pz_true,0,0,domain%mesh_n)

    call scalar_gen%get_scalar_field(p,domain%mesh_p,0)
    call scalar_gen%grad%get_vertical_component(pz_true,domain%mesh_n,0,"covariant")

    print *, p%tile(1)%p(1,1,2)-p%tile(1)%p(1,1,1)
    print *, pz_true%tile(1)%p(1,1,2)*3e3_8

end subroutine test_vertical_gradient_operator

end module test_vertical_operators_mod
