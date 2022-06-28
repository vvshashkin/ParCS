module orography_factory_mod

use orography_mod,            only : orography_t
use config_mod,               only : config_t
use domain_mod,               only : domain_t
use grid_field_factory_mod,   only : create_grid_field
use parcomm_mod,              only : parcomm_global
use orography_test_field_mod, only : orography_test_field_t, orography_test_grad_t
use test_fields_3d_mod,       only : scalar_field3d_t, vector_field3d_t

implicit none

contains

subroutine create_orography(orography,orography_name,config,domain,halo_width)
    type(orography_t), intent(out) :: orography
    character(len=*),  intent(in)  :: orography_name
    class(config_t),   intent(in)  :: config
    type(domain_t),    intent(in)  :: domain
    integer(kind=4),   intent(in)  :: halo_width

    call create_grid_field(orography%o%h,          halo_width,0,domain%mesh_o)
    call create_grid_field(orography%o%dh_alpha,   halo_width,0,domain%mesh_o)
    call create_grid_field(orography%o%dh_beta,    halo_width,0,domain%mesh_o)

    call create_grid_field(orography%x%h,          halo_width,0,domain%mesh_x)
    call create_grid_field(orography%x%dh_alpha,   halo_width,0,domain%mesh_x)
    call create_grid_field(orography%x%dh_beta,    halo_width,0,domain%mesh_x)

    call create_grid_field(orography%y%h,          halo_width,0,domain%mesh_y)
    call create_grid_field(orography%y%dh_alpha,   halo_width,0,domain%mesh_y)
    call create_grid_field(orography%y%dh_beta,    halo_width,0,domain%mesh_y)

    call create_grid_field(orography%xy%h,         halo_width,0,domain%mesh_xy)
    call create_grid_field(orography%xy%dh_alpha,  halo_width,0,domain%mesh_xy)
    call create_grid_field(orography%xy%dh_beta,   halo_width,0,domain%mesh_xy)

    if(orography_name == "zero_orography") then
        call orography%o%h       %assign(0.0_8,domain%mesh_o,halo_width)
        call orography%o%dh_alpha%assign(0.0_8,domain%mesh_o,halo_width)
        call orography%o%dh_beta %assign(0.0_8,domain%mesh_o,halo_width)

        call orography%x%h       %assign(0.0_8,domain%mesh_x,halo_width)
        call orography%x%dh_alpha%assign(0.0_8,domain%mesh_x,halo_width)
        call orography%x%dh_beta %assign(0.0_8,domain%mesh_x,halo_width)

        call orography%y%h       %assign(0.0_8,domain%mesh_y,halo_width)
        call orography%y%dh_alpha%assign(0.0_8,domain%mesh_y,halo_width)
        call orography%y%dh_beta %assign(0.0_8,domain%mesh_y,halo_width)

        call orography%xy%h       %assign(0.0_8,domain%mesh_xy,halo_width)
        call orography%xy%dh_alpha%assign(0.0_8,domain%mesh_xy,halo_width)
        call orography%xy%dh_beta %assign(0.0_8,domain%mesh_xy,halo_width)
    else if(orography_name == "test_orography") then
        call create_analytic_orography(orography,orography_test_field_t(1.0_8), &
                                       orography_test_grad_t(1.0_8,domain%mesh_o%scale),&
                                       domain, halo_width)
    else
        call parcomm_global%abort("orography factory mod error - unknown orography name:"//&
                                   orography_name)
    end if

end subroutine create_orography

subroutine create_analytic_orography(orography,orography_gen,orography_grad_gen,&
                                                               domain,halo_width)
    type(orography_t),       intent(inout) :: orography
    class(scalar_field3d_t), intent(in)    :: orography_gen
    class(vector_field3d_t), intent(in)    :: orography_grad_gen
    type(domain_t),          intent(in)    :: domain
    integer(kind=4),         intent(in)    :: halo_width

    call orography_gen%get_scalar_field(orography%o%h, domain%mesh_o, halo_width)
    call orography_gen%get_scalar_field(orography%x%h, domain%mesh_x, halo_width)
    call orography_gen%get_scalar_field(orography%y%h, domain%mesh_y, halo_width)
    call orography_gen%get_scalar_field(orography%xy%h,domain%mesh_xy,halo_width)

    call orography_grad_gen%get_x_component(orography%o%dh_alpha, domain%mesh_o, halo_width,'covariant')
    call orography_grad_gen%get_y_component(orography%o%dh_beta , domain%mesh_o, halo_width,'covariant')
    call orography_grad_gen%get_x_component(orography%x%dh_alpha, domain%mesh_x, halo_width,'covariant')
    call orography_grad_gen%get_y_component(orography%x%dh_beta , domain%mesh_x, halo_width,'covariant')
    call orography_grad_gen%get_x_component(orography%y%dh_alpha, domain%mesh_y, halo_width,'covariant')
    call orography_grad_gen%get_y_component(orography%y%dh_beta , domain%mesh_y, halo_width,'covariant')
    call orography_grad_gen%get_x_component(orography%xy%dh_alpha,domain%mesh_xy,halo_width,'covariant')
    call orography_grad_gen%get_y_component(orography%xy%dh_beta ,domain%mesh_xy,halo_width,'covariant')

end subroutine create_analytic_orography

end module orography_factory_mod
