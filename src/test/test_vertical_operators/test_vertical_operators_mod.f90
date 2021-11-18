module test_vertical_operators_mod

use abstract_vertical_operator_mod, only: vertical_operator_t
use vertical_operator_factory_mod,  only: create_vertical_operator

use test_fields_3d_mod,       only: scalar_field3d_t, vector_field3d_t
use domain_mod,               only: domain_t
use domain_factory_mod,       only: create_domain
use config_domain_mod,        only: config_domain_t
use grid_field_mod,           only: grid_field_t
use grid_field_factory_mod,   only: create_grid_field
use vec_math_mod,             only : l2norm

implicit none

contains

subroutine test_vertical_gradient_operator(nz,vertical_grad_name,vertical_staggering)

    use vertical_test_field_mod,  only: vertical_ExnerP_t
    use const_N_profile_mod,      only: const_N_profile_t

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
    call create_grid_field(pz,0,0,domain%mesh_w)
    call create_grid_field(pz_true,0,0,domain%mesh_w)

    call scalar_gen%get_scalar_field(p,domain%mesh_p,0)
    call scalar_gen%grad%get_vertical_component(pz_true,domain%mesh_w,0,"covariant")

    call create_vertical_operator(vert_grad,vertical_grad_name)
    call vert_grad%apply(pz,p,domain)

    call pz%update(-1.0_8,pz_true,domain%mesh_w)
    print *, vertical_grad_name
    print *, "rel l2 error:", l2norm(pz, domain%mesh_w,domain%parcomm) / &
                              l2norm(pz_true, domain%mesh_w,domain%parcomm)
    print *, "rel linf error:", pz%maxabs(domain%mesh_w,domain%parcomm) / &
                                pz_true%maxabs(domain%mesh_w,domain%parcomm)

end subroutine test_vertical_gradient_operator

subroutine test_vertical_div_operator(nz,vertical_div_name,vertical_staggering)

    use vertical_div_test_field_mod,  only: simple_divergent_w_t

    integer(kind=4),  intent(in) :: nz
    character(len=*), intent(in) :: vertical_div_name, vertical_staggering

    integer, parameter :: nh = 8
    real(kind=8), parameter :: h_top = 30e3_8

    class(vector_field3d_t), allocatable :: vector_gen
    type(config_domain_t) :: config_domain
    type(domain_t)     :: domain
    type(grid_field_t) :: w, w_div,w_div_true

    class(vertical_operator_t), allocatable :: vert_div

    vector_gen = simple_divergent_w_t(h_top)

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

    call create_grid_field(w,0,0,domain%mesh_w)
    call create_grid_field(w_div,0,0,domain%mesh_p)
    call create_grid_field(w_div_true,0,0,domain%mesh_p)

    call vector_gen%get_vertical_component(w,domain%mesh_w,0,"contravariant")
    call vector_gen%div%get_scalar_field(w_div_true,domain%mesh_p,0)

    call create_vertical_operator(vert_div,vertical_div_name)
    call vert_div%apply(w_div,w,domain)

    call w_div%update(-1.0_8,w_div_true,domain%mesh_p)
    print *, vertical_div_name
    print *, "rel l2 error:", l2norm(w_div, domain%mesh_p,domain%parcomm) / &
                              l2norm(w_div_true, domain%mesh_p,domain%parcomm)
    print *, "rel linf error:", w_div%maxabs(domain%mesh_p,domain%parcomm) / &
                                w_div_true%maxabs(domain%mesh_p,domain%parcomm)

end subroutine test_vertical_div_operator

end module test_vertical_operators_mod
