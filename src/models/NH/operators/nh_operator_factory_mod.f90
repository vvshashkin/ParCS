module nh_operator_factory_mod

use domain_mod,             only : domain_t
use operator_mod,           only : operator_t
use config_mod,             only : config_t
use config_nh_operator_mod, only : config_Ptheta_linear_t
use parcomm_mod,            only : parcomm_global

implicit none

contains

subroutine create_nh_operator(nh_operator, config, domain)
    class(operator_t), allocatable, intent(out) :: nh_operator
    class(config_t),   intent(in)  :: config
    type(domain_t),    intent(in)  :: domain

    select type(config)
    class is (config_Ptheta_linear_t)
        call create_Ptheta_linear_nh_operator(nh_operator, config, domain)
    class default
        call parcomm_global%abort("unknown operator config type in create_nh_operator")
    end select
end subroutine create_nh_operator

subroutine create_Ptheta_linear_nh_operator(nh_operator,config,domain)

    use Ptheta_linear_nh_oper_mod,     only: Ptheta_linear_nh_operator_t
    use grad_3d_factory_mod,           only: create_grad_3d_operator
    use div_3d_factory_mod,            only: create_div_3d_operator
    use co2contra_factory_mod,         only: create_co2contra_operator
    use interpolator_w2uv_factory_mod, only: create_w2uv_interpolator
    use vertical_operator_factory_mod, only: create_vertical_operator
    use grid_field_factory_mod,        only: create_grid_field
    use vertical_test_field_mod,       only: vertical_ExnerP_t, vertical_ExnerP_grad_t, &
                                             vertical_theta_t,  vertical_theta_grad_t
    use const_N_profile_mod,           only: const_N_profile_t

    class(operator_t), allocatable,  intent(out) :: nh_operator
    class(config_Ptheta_linear_t),   intent(in)  :: config
    type(domain_t),                  intent(in)  :: domain

    type(Ptheta_linear_nh_operator_t), allocatable :: operator
    type(vertical_ExnerP_t)      :: P0_generator
    type(vertical_ExnerP_grad_t) :: dP0dz_generator
    type(vertical_theta_t)       :: theta0_generator
    type(vertical_theta_grad_t)  :: dtheta0dz_generator

    integer(kind=4), parameter :: halo_width_xy = 4

    allocate(operator)

    call create_grad_3d_operator(operator%grad_op, domain, &
                                 config%grad_hor_part_name, config%grad_vert_part_name)

    call create_div_3d_operator(operator%div_op, domain, &
                                config%div_hor_part_name, config%div_vert_part_name)

    operator%co2contra_op = create_co2contra_operator(domain, config%co2contra_operator_name)

    call create_w2uv_interpolator(operator%w2uv_op,config%w2uv_operator_name,domain)

    call create_vertical_operator(operator%w2p_op, config%w2p_operator_name)

    call create_grid_field(operator%P0, 0, 0, domain%mesh_p)
    call create_grid_field(operator%dP0_dz, 0, 0, domain%mesh_p)
    call create_grid_field(operator%dP0_dz_w, 0, 0, domain%mesh_w)
    call create_grid_field(operator%theta0, 0, 0, domain%mesh_w)
    call create_grid_field(operator%theta0_u, 0, 0, domain%mesh_u)
    call create_grid_field(operator%theta0_v, 0, 0, domain%mesh_v)
    call create_grid_field(operator%dtheta0_dz, 0, 0, domain%mesh_w)
    call create_grid_field(operator%w, 0, 0, domain%mesh_w)
    call create_grid_field(operator%wp, 0, 0, domain%mesh_p)
    call create_grid_field(operator%div3, 0, 0, domain%mesh_p)
    call create_grid_field(operator%grad_x, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(operator%grad_y, halo_width_xy, 0, domain%mesh_v)
    call create_grid_field(operator%grad_z, 0, 0, domain%mesh_w)
    !call create_grid_field(operator%grad_x_contra, 0, 0, domain%mesh_u)
    !call create_grid_field(operator%grad_y_contra, 0, 0, domain%mesh_v)

    select case(config%background_type)
    case("Nb_const")
        P0_generator%t0 = 300.0
        P0_generator%p0 = 1e5
        P0_generator%vert_profile = const_N_profile_t(N=config%Nb)
        dP0dz_generator%t0 = 300.0
        dP0dz_generator%p0 = 1e5
        dP0dz_generator%vert_profile = const_N_profile_t(N=config%Nb)
        theta0_generator%t0 = 300.0
        theta0_generator%p0 = 1e5
        theta0_generator%vert_profile = const_N_profile_t(N=config%Nb)
        dtheta0dz_generator%t0 = 300.0
        dtheta0dz_generator%p0 = 1e5
        dtheta0dz_generator%vert_profile = const_N_profile_t(N=config%Nb)
    case default
        call parcomm_global%abort("create_Ptheta_linear_nh_operator, "//&
                                  "unknown type of background profile: "//&
                                  config%background_type)
    end select

    call P0_generator%get_scalar_field(operator%P0,domain%mesh_p,0)
    call dP0dz_generator%get_vertical_component(operator%dP0_dz,domain%mesh_p,0,"covariant")
    call dP0dz_generator%get_vertical_component(operator%dP0_dz_w,domain%mesh_w,0,"covariant")
    call theta0_generator%get_scalar_field(operator%theta0,domain%mesh_w,0)
    call theta0_generator%get_scalar_field(operator%theta0_u,domain%mesh_u,0)
    call theta0_generator%get_scalar_field(operator%theta0_v,domain%mesh_v,0)
    call dtheta0dz_generator%get_vertical_component(operator%dtheta0_dz,domain%mesh_w,0,"covariant")
    call move_alloc(operator, nh_operator)

end subroutine create_Ptheta_linear_nh_operator

end module nh_operator_factory_mod
