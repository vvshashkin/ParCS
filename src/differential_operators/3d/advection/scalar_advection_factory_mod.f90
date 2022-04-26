module scalar_advection_factory_mod

use abstract_scalar_advection3d_mod, only : scalar_advection3d_t
use domain_mod,                      only : domain_t
use parcomm_mod,                     only : parcomm_global
use config_mod,                      only : config_t
use v_nabla_factory_mod,             only : create_v_nabla_hor_operator
use adv_z_factory_mod,               only : create_adv_z_operator
use halo_factory_mod,                only : create_halo_procedure
use grid_field_factory_mod,          only : create_grid_field
use interpolator_uv2w_factory_mod,   only : create_uv2w_interpolator
use interpolator2d_factory_mod,      only : create_vec2vec_interpolator2d
use vertical_operator_factory_mod,   only : create_vertical_operator

implicit none

contains

subroutine create_scalar_advection3d_operator(adv_op, scalar_advection_op_name, &
                                                                  config, domain)
    class(scalar_advection3d_t), allocatable, intent(out) :: adv_op
    character(len=*), intent(in) :: scalar_advection_op_name
    class(config_t),  intent(in) :: config
    type(domain_t),   intent(in) :: domain

    select case(scalar_advection_op_name)
    case("advection_p_staggered")
        call create_p_3d_advection(adv_op, config, domain)
    case("advection_w_staggered")
        call create_w_3d_advection(adv_op, config, domain)
    case default
        call parcomm_global%abort("create_scalar_advection3d_operator, unknown "//&
                                  "scalar_advection_op_name: "// scalar_advection_op_name)
    end select
end subroutine create_scalar_advection3d_operator

subroutine create_p_3d_advection(adv_op, config, domain)

    use advection_p_3d_mod,      only : advection_p_C3d_t
    use config_advection_3d_mod, only : config_p_advection_t

    class(scalar_advection3d_t), allocatable, intent(out) :: adv_op
    class(config_t),  intent(in) :: config
    type(domain_t),   intent(in) :: domain

    type(advection_p_C3d_t), allocatable :: adv_p3d
    integer(kind=4) :: halo_width

    allocate(adv_p3d)

    select type(config)
    class is (config_p_advection_t)

    call create_vec2vec_interpolator2d(adv_p3d%interp_uv2p_op, config%uv2p_operator_name, domain)
    call create_vertical_operator(adv_p3d%interp_w2p_op, config%w2p_operator_name)

    call create_v_nabla_hor_operator(adv_p3d%v_nabla_op,halo_width, &
                                     config%hor_advection_oper_name)

    adv_p3d%halo_width = halo_width
    call create_halo_procedure(adv_p3d%halo_f,domain,max(halo_width,2),config%p_halo)

    call create_adv_z_operator(adv_p3d%adv_z, config%z_advection_oper_name)
    call create_grid_field(adv_p3d%up, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%vp, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%eta_dot_p, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%f_tend_z, 0, 0, domain%mesh_p)

    class default
        call parcomm_global%abort("create_p_3d_advection type error")
    end select

    call move_alloc(adv_p3d, adv_op)

end subroutine create_p_3d_advection

subroutine create_w_3d_advection(adv_op, config, domain)

    use advection_w_3d_mod,      only : advection_w_C3d_t
    use config_advection_3d_mod, only : config_w_advection_t


    class(scalar_advection3d_t), allocatable, intent(out) :: adv_op
    class(config_t),  intent(in) :: config
    type(domain_t),   intent(in) :: domain

    type(advection_w_C3d_t), allocatable :: adv_w3d
    integer(kind=4) :: halo_width

    allocate(adv_w3d)

    select type(config)
    class is (config_w_advection_t)

    call create_uv2w_interpolator(adv_w3d%interp_uv2w_op, config%uv2w_operator_name, &
                                  config%uv2w_hor_part_name, config%uv2w_vert_part_name, &
                                  domain)

    call create_v_nabla_hor_operator(adv_w3d%v_nabla_op,halo_width,config%hor_advection_oper_name)

    adv_w3d%halo_width = halo_width
    call create_halo_procedure(adv_w3d%halo_f,domain,max(halo_width,2),config%w_halo)

    call create_adv_z_operator(adv_w3d%adv_z, config%z_advection_oper_name)
    call create_grid_field(adv_w3d%uw, 0, 0, domain%mesh_w)
    call create_grid_field(adv_w3d%vw, 0, 0, domain%mesh_w)
    call create_grid_field(adv_w3d%f_tend_z, 0, 0, domain%mesh_w)

    class default
        call parcomm_global%abort("create_w_3d_advection type error")
    end select

    call move_alloc(adv_w3d, adv_op)

end subroutine create_w_3d_advection

end module scalar_advection_factory_mod
