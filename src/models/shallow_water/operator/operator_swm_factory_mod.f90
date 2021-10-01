module operator_swm_factory_mod

use domain_mod,     only : domain_t
use operator_mod,   only : operator_t
use config_swm_mod, only : config_swm_t

implicit none

contains

subroutine create_swm_operator(operator, swm_config, domain)

    use operator_swm_mod,       only : operator_swm_t
    use div_factory_mod,        only : create_div_operator
    use grad_factory_mod,       only : create_grad_operator
    use curl_factory_mod,       only : create_curl_operator_div_based
    use coriolis_factory_mod,   only : create_coriolis
    use KE_factory_mod,         only : create_KE_operator
    use massflux_factory_mod,   only : create_massflux_operator
    use co2contra_factory_mod,  only : create_co2contra_operator

    use grid_field_factory_mod, only : create_grid_field

    type(domain_t),                 intent(in)  :: domain
    type(config_swm_t),             intent(in)  :: swm_config
    class(operator_t), allocatable, intent(out) :: operator

    type(operator_swm_t), allocatable :: swm_op

    integer(kind=4) :: halo_width_xy, t

    !WORKAROUND
    halo_width_xy = 8

    allocate(swm_op)

    swm_op%div_op = create_div_operator(domain, swm_config%div_op_name)

    swm_op%grad_op =  create_grad_operator(domain, swm_config%grad_op_name)

    call create_coriolis(swm_op%coriolis_op, swm_config%coriolis_op_name, domain)

    call create_curl_operator_div_based(swm_op%curl_op, swm_config%div_op_name, domain)

    call create_KE_operator(swm_op%KE_op, swm_config%KE_op_name, domain)

    swm_op%massflux_op = create_massflux_operator(domain, swm_config%massflux_op_name)

    swm_op%co2contra_op = create_co2contra_operator(domain, swm_config%co2contra_op_name)

    call create_grid_field(swm_op%KE,  halo_width_xy, 0, domain%mesh_p)
    call create_grid_field(swm_op%div, halo_width_xy, 0, domain%mesh_p)


    call create_grid_field(swm_op%curl, halo_width_xy, 0, domain%mesh_p)

    call create_grid_field(swm_op%hu, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%hv, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%cor_u, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%cor_v, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%ut, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%vt, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%grad_x, halo_width_xy, 0, domain%mesh_u)
    call create_grid_field(swm_op%grad_y, halo_width_xy, 0, domain%mesh_v)

    call create_grid_field(swm_op%h_surf,    halo_width_xy, 0, domain%mesh_p)

    !WORKAROUND
    swm_op%grav = 1.0_8
    do t = domain%partition%ts, domain%partition%te
        swm_op%h_surf%tile(t)%p = 0.0_8
    end do

    call move_alloc(swm_op, operator)

end subroutine create_swm_operator

end module operator_swm_factory_mod
