module scalar_advection_factory_mod

use abstract_scalar_advection3d_mod, only : scalar_advection3d_t
use domain_mod,                      only : domain_t
use parcomm_mod,                     only : parcomm_global

implicit none

contains

subroutine create_scalar_advection3d_operator(adv_op, scalar_advection_op_name, &
                                              hor_adv_op_name, z_adv_op_name, domain)
    class(scalar_advection3d_t), allocatable, intent(out) :: adv_op
    character(len=*), intent(in) :: scalar_advection_op_name, hor_adv_op_name, &
                                    z_adv_op_name
    type(domain_t),   intent(in) :: domain

    select case(scalar_advection_op_name)
    case("advection_p_staggered")
        call create_p_3d_advection(adv_op, hor_adv_op_name, z_adv_op_name, domain)
    case default
        call parcomm_global%abort("create_scalar_advection3d_operator, unknown "//&
                                  "scalar_advection_op_name: "// scalar_advection_op_name)
    end select
end subroutine create_scalar_advection3d_operator

subroutine create_p_3d_advection(adv_op, hor_adv_op_name, z_adv_op_name, domain)

    use v_nabla_mod,        only : v_nabla_c2_operator_t, v_nabla_c4_operator_t,   &
                                   v_nabla_up1_operator_t, v_nabla_up3_operator_t, &
                                   v_nabla_up4_operator_t
    use adv_z_factory_mod,  only : create_adv_z_operator
    use halo_factory_mod,   only : create_halo_procedure
    use advection_p_3d_mod, only : advection_p_C3d_t

    use interpolator_v2h_factory_mod,  only : create_v2h_interpolator
    use vertical_operator_factory_mod, only : create_vertical_operator
    use grid_field_factory_mod,        only : create_grid_field

    class(scalar_advection3d_t), allocatable, intent(out) :: adv_op
    character(len=*), intent(in) :: hor_adv_op_name, z_adv_op_name
    type(domain_t),   intent(in) :: domain

    type(advection_p_C3d_t), allocatable :: adv_p3d
    integer(kind=4) :: halo_width

    allocate(adv_p3d)

    call create_v2h_interpolator(adv_p3d%interp_uv2p_op, "W42_stagered_interp_c2i", domain)
    call create_vertical_operator(adv_p3d%interp_w2p_op, "vertical_interp_w2p_sbp42")

    select case(hor_adv_op_name)
    case("c2")
        adv_p3d%v_nabla_op = v_nabla_c2_operator_t()
        halo_width = 1
    case("c4")
        adv_p3d%v_nabla_op = v_nabla_c4_operator_t()
        halo_width = 2
    case("up1")
        adv_p3d%v_nabla_op = v_nabla_up1_operator_t()
        halo_width = 1
    case("up3")
        adv_p3d%v_nabla_op = v_nabla_up3_operator_t()
        halo_width = 2
    case("up4")
        adv_p3d%v_nabla_op = v_nabla_up4_operator_t()
        halo_width = 3
    case default
        call parcomm_global%abort("create_p_3d_advection, unknown horizontal advection operator:"//&
                                  hor_adv_op_name)
    end select

    adv_p3d%halo_width = halo_width
    call create_halo_procedure(adv_p3d%halo_f,domain,max(halo_width,2),"ECS_O")

    call create_adv_z_operator(adv_p3d%adv_z, z_adv_op_name)
    call create_grid_field(adv_p3d%up, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%vp, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%eta_dot_p, 0, 0, domain%mesh_p)
    call create_grid_field(adv_p3d%f_tend_z, 0, 0, domain%mesh_p)

    call move_alloc(adv_p3d, adv_op)

end subroutine create_p_3d_advection
end module scalar_advection_factory_mod
