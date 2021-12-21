module test_diffops_3d_mod

use domain_mod,             only : domain_t
use domain_factory_mod,     only : create_domain
use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field
use parcomm_mod,            only : parcomm_global
use vec_math_mod,           only : l2norm

use key_value_mod,          only : key_value_r8_t

implicit none

private
public :: test_scalar_advection_3d

contains

function test_scalar_advection_3d(Nh, Nz, advection_oper_name, hor_advection_oper_name, &
                                          vert_advection_oper_name, points_type, &
                                          horizontal_staggering, vertical_staggering) result(errs)
    use const_mod, only : Earth_radii

    use abstract_scalar_advection3d_mod, only : scalar_advection3d_t
    use scalar_advection_factory_mod,    only : create_scalar_advection3d_operator
    use config_domain_mod,               only : config_domain_t

    use mesh_mod,                        only : mesh_t

    use test_fields_3d_mod,              only : scalar_field3d_t, vector_field3d_t
    use grad3d_test_field_mod,           only : grad3d_test_input_t, grad3d_test_out_t

    use solid_rotation_wind_field_mod,   only : solid_rotation_wind_field_t

    integer(kind=4),  intent(in) :: Nh, Nz
    character(len=*), intent(in) :: advection_oper_name, hor_advection_oper_name, &
                                    vert_advection_oper_name, points_type, &
                                    horizontal_staggering, vertical_staggering

    type(key_value_r8_t) :: errs

    class(scalar_advection3d_t), allocatable :: advection_op
    type(config_domain_t)                    :: config_domain
    type(domain_t), target                   :: domain

    type(grad3d_test_input_t) :: scalar_generator
    type(grad3d_test_out_t)   :: grad_generator
    class(vector_field3d_t), allocatable    :: wind_generator

    type(grid_field_t) :: f, fx, fy, fz, v_nabla_f
    type(grid_field_t) :: u, v, eta_dot
    type(grid_field_t) :: u_p, v_p, eta_dot_p
    type(mesh_t), pointer :: mesh

    integer(kind=4) :: halo_width = 8
    real(kind=8), parameter :: h_top = 30e3, T0 = 300._8

    config_domain%N  = nh
    config_domain%Nz = nz
    config_domain%staggering_type     = horizontal_staggering
    config_domain%vertical_staggering = vertical_staggering
    config_domain%metric_type         = "shallow_atmosphere_metric"
    config_domain%topology_type       = "cube"
    config_domain%h_top = h_top
    call config_domain%config_metric%set_defaults()
    config_domain%config_metric%vertical_scale = h_top
    config_domain%config_metric%scale = Earth_radii

    call create_domain(domain, config_domain)

    select case(points_type)
    case("p")
        mesh => domain%mesh_p
    case("w")
        mesh => domain%mesh_w
    case default
        call parcomm_global%abort("unsupported points type in test_scalar_advection_3d: "//&
                                  points_type)
    end select


    call create_grid_field(f,  halo_width, 0, mesh)
    call create_grid_field(v_nabla_f,  0, 0, mesh)

    call create_grid_field(fx, 0, 0, mesh)
    call create_grid_field(fy, 0, 0, mesh)
    call create_grid_field(fz, 0, 0, mesh)
    call create_grid_field(u, halo_width, 0, domain%mesh_u)
    call create_grid_field(v, halo_width, 0, domain%mesh_v)
    call create_grid_field(eta_dot, halo_width, 0, domain%mesh_w)
    call create_grid_field(u_p, 0, 0, mesh)
    call create_grid_field(v_p, 0, 0, mesh)
    call create_grid_field(eta_dot_p, 0, 0, mesh)

    scalar_generator = grad3d_test_input_t(h_top=h_top,T0=T0)
    grad_generator   = grad3d_test_out_t(h_top=h_top,T0=T0)
    wind_generator   = solid_rotation_wind_field_t(w_max = 0.1_8, u0 = 1.0_8, t= 0.0_8)

    call scalar_generator%get_scalar_field(f, mesh, 0)
    call grad_generator%get_vector_field(fx, fy, fz, mesh, mesh, mesh, &
                                         0,"covariant")
    call wind_generator%get_vector_field(u,v,eta_dot,domain%mesh_u,domain%mesh_v, &
                                         domain%mesh_w,0, "contravariant")
    call wind_generator%get_vector_field(u_p,v_p,eta_dot_p, mesh, mesh, mesh, &
                                                       0, "contravariant")

    !call grad_3d%calc_grad(fx, fy, fz, f, domain)
    call create_scalar_advection3d_operator(advection_op,advection_oper_name, &
                                            hor_advection_oper_name,          &
                                            vert_advection_oper_name,domain)
    call advection_op%calc_adv3d(v_nabla_f,f,u,v,eta_dot,domain)

    allocate(errs%keys(2), errs%values(2))
    errs%keys(1)%str = "scalar advection 3d L_inf"
    errs%keys(2)%str = "scalar advection 3d L2"

    call fx%assign_prod(-1.0_8,fx,u_p,mesh)
    call fy%assign_prod(-1.0_8,fy,v_p,mesh)
    call fz%assign_prod(-1.0_8,fz,eta_dot_p,mesh)

    call f%assign(1.0_8,fx,mesh)
    call f%update(1.0_8,fy,1.0_8,fz,mesh)
    call v_nabla_f%update(-1.0_8,f,mesh)

    errs%values(1) = v_nabla_f%maxabs(mesh, domain%parcomm) / f%maxabs(mesh, domain%parcomm)
    errs%values(2) = l2norm(v_nabla_f, mesh, domain%parcomm) / l2norm(f, mesh, domain%parcomm)
end function test_scalar_advection_3d

end module test_diffops_3d_mod
