module vector_advection3d_factory_mod

use abstract_vector_advection3d_mod, only : vector_advection3d_t
use domain_mod,                      only : domain_t
use parcomm_mod,                     only : parcomm_global
use config_mod,                      only : config_t

implicit none

contains

subroutine create_vector_advection3d_operator(vec_adv_op, vec_adv_op_name, config, domain)

    class(vector_advection3d_t), allocatable, intent(out) :: vec_adv_op
    !Input:
    character(len=*), intent(in) :: vec_adv_op_name
    class(config_t),  intent(in) :: config
    type(domain_t),   intent(in) :: domain

    select case(vec_adv_op_name)
    case("shallow_atm_staggered_vector_advection")
        call create_staggered_shallow_atm_vector_advection3d(vec_adv_op, config, domain)
    case default
        call parcomm_global%abort("create_vector_advection3d_operator error - "// &
                                  "unknown vec_adv_op_name: "// vec_adv_op_name)
    end select
end subroutine create_vector_advection3d_operator

subroutine create_staggered_shallow_atm_vector_advection3d(vec_adv_op,config,domain)

    use shallow_atm_vecadv_mod,        only : shallow_atm_staggered_vecadv_t
    use vector_advection_factory_mod,  only : create_vector_advection_operator
    use adv_z_factory_mod,             only : create_adv_z_operator
    use interpolator_w2uv_factory_mod, only : create_w2uv_interpolator
    use scalar_advection_factory_mod,  only : create_scalar_advection3d_operator
    use grid_field_factory_mod,        only : create_grid_field
    use config_advection_3d_mod,       only : config_vector_advection_3d_t

    class(vector_advection3d_t), allocatable, intent(out) :: vec_adv_op
    !Input:
    class(config_t),  intent(in) :: config
    type(domain_t),   intent(in) :: domain

    type(shallow_atm_staggered_vecadv_t), allocatable :: operator

    allocate(operator)

    select type(config)
    class is (config_vector_advection_3d_t)

    call create_vector_advection_operator(operator%uv_hor_advection_op, &
                                              config%uv_hor_advection_oper_name, domain)
    call create_adv_z_operator(operator%uv_z_advec_op,config%uv_ver_advection_oper_name)
    call create_w2uv_interpolator(operator%interp_w2uv, config%w2uv_operator_name,       &
                                  config%w2uv_hor_part_name, config%w2uv_vert_part_name, domain)

    call create_scalar_advection3d_operator(operator%w_advection3d_op, config%w_advection_oper_name, &
                                            config%w_advection_config, domain)

    call create_grid_field(operator%eta_dot_u,0,0,domain%mesh_u)
    call create_grid_field(operator%eta_dot_v,0,0,domain%mesh_v)
    call create_grid_field(operator%u_tend_z,0,0,domain%mesh_u)
    call create_grid_field(operator%v_tend_z,0,0,domain%mesh_v)

    class default
        call parcomm_global%abort("create_staggered_shallow_atm_vector_advection3d "// &
                                  "config type error")
    end select

    call move_alloc(operator,vec_adv_op)
end subroutine create_staggered_shallow_atm_vector_advection3d

end module vector_advection3d_factory_mod
