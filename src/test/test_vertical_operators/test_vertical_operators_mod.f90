module test_vertical_operators_mod

use abstract_vertical_operator_mod, only: vertical_operator_t
use vertical_operator_factory_mod,  only: create_vertical_operator

use test_fields_3d_mod,       only: scalar_field3d_t
use vertical_test_field_mod,  only: vertical_ExnerP_t
use const_N_profile_mod,      only: const_N_profile_t
use domain_mod,               only: domain_t
use domain_factory_mod,       only: create_domain
use config_domain_mod,        only: config_domain_t
use grid_field_mod,           only: grid_field_t
use grid_field_factory_mod,   only: create_grid_field
use vec_math_mod,             only : l2norm

implicit none

contains

subroutine test_vertical_gradient_operator(nz,vertical_grad_name,vertical_staggering)

    integer(kind=4),  intent(in) :: nz
    character(len=*), intent(in) :: vertical_grad_name, vertical_staggering

    integer, parameter :: nh = 8
    real(kind=8), parameter :: h_top = 30e3_8

    class(scalar_field3d_t), allocatable :: scalar_gen
    type(config_domain_t) :: config_domain
    type(domain_t)     :: domain
    type(grid_field_t) :: p, pz_true, pz

    ! integer(kind=4) :: k

    class(vertical_operator_t), allocatable :: vert_grad

    scalar_gen = vertical_ExnerP_t(t0=300.0_8, p0=1e5_8, &
                                   vert_profile=const_N_profile_t(N=0.01))

    config_domain%N  = nh
    config_domain%Nz = nz
    config_domain%staggering_type     = "C"
    config_domain%vertical_staggering = vertical_staggering
    config_domain%metric_type         = "shallow_atmosphere_metric"
    config_domain%topology_type       = "cube"
    config_domain%h_top = h_top
    call config_domain%config_metric%set_defaults()
    config_domain%config_metric%vertical_scale = h_top

    call create_domain(domain, config_domain)

    call create_grid_field(p,0,0,domain%mesh_p)
    call create_grid_field(pz,0,0,domain%mesh_n)
    call create_grid_field(pz_true,0,0,domain%mesh_n)

    call scalar_gen%get_scalar_field(p,domain%mesh_p,0)
    call scalar_gen%grad%get_vertical_component(pz_true,domain%mesh_n,0,"covariant")

    call create_vertical_operator(vert_grad,vertical_grad_name)
    call vert_grad%apply(pz,p,domain)

    call pz%update(-1.0_8,pz_true,domain%mesh_n)
    print *, vertical_grad_name
    print *, "rel l2 error:", l2norm(pz, domain%mesh_n,domain%parcomm) / &
                              l2norm(pz_true, domain%mesh_n,domain%parcomm)
    print *, "rel linf error:", pz%maxabs(domain%mesh_n,domain%parcomm) / &
                                pz_true%maxabs(domain%mesh_n,domain%parcomm)

end subroutine test_vertical_gradient_operator

end module test_vertical_operators_mod
