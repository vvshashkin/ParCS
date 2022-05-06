module test_diffops_3d_mod

use domain_mod,             only : domain_t
use domain_factory_mod,     only : create_domain
use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field
use parcomm_mod,            only : parcomm_global
use vec_math_mod,           only : l2norm
use config_mod,             only : config_t

use key_value_mod,          only : key_value_r8_t

implicit none

private
public :: test_w2uv_interp, test_uv2w_interp, test_scalar_advection_3d, &
          test_grad_3d, test_div_3d, test_co2contra_3d

contains

function test_grad_3d(Nh, Nz, hor_grad_name, diff_eta_name, &
                      horizontal_staggering, vertical_staggering) result(errs)

    use const_mod, only : Earth_radii

    use grad_3d_factory_mod,   only : create_grad_3d_operator
    use abstract_grad_3d_mod,  only : grad_3d_operator_t
    use config_domain_mod,     only : config_domain_t

    use test_fields_3d_mod,    only : scalar_field3d_t, vector_field3d_t
    use grad3d_test_field_mod, only : grad3d_test_input_t, grad3d_test_out_t

    integer(kind=4), intent(in)  :: Nh, Nz
    character(len=*), intent(in) :: hor_grad_name, diff_eta_name
    character(len=*), intent(in) :: horizontal_staggering, vertical_staggering

    type(key_value_r8_t) :: errs

    class(grad_3d_operator_t), allocatable :: grad_3d
    type(config_domain_t) :: config_domain
    type(domain_t)        :: domain

    type(grad3d_test_input_t) :: scalar_generator
    type(grad3d_test_out_t)   :: grad_generator

    type(grid_field_t) :: f, fx, fy, fz, fx_true, fy_true, fz_true

    integer(kind=4) :: halo_width = 8
    real(kind=8), parameter :: h_top = 30e3, T0(2) = [300._8,1e10_8]
    integer(kind=4) :: itemp

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

    call create_grad_3d_operator(grad_3d, domain, hor_grad_name, diff_eta_name)

    call create_grid_field(f,  halo_width, 0, domain%mesh_p)

    call create_grid_field(fx, 0, 0, domain%mesh_u)
    call create_grid_field(fy, 0, 0, domain%mesh_v)
    call create_grid_field(fz, 0, 0, domain%mesh_w)

    call create_grid_field(fx_true, 0, 0, domain%mesh_u)
    call create_grid_field(fy_true, 0, 0, domain%mesh_v)
    call create_grid_field(fz_true, 0, 0, domain%mesh_w)

    allocate(errs%keys(4), errs%values(4))
    errs%keys(1)%str = "verthor_ExnerP C_norm"
    errs%keys(2)%str = "verthor_ExnerP L2_norm"
    errs%keys(3)%str = "quasihor_ExnerP C_norm"
    errs%keys(4)%str = "quasihor_ExnerP L2_norm"

    do itemp = 1, size(T0)
        scalar_generator = grad3d_test_input_t(h_top=h_top,T0=T0(itemp))
        grad_generator   = grad3d_test_out_t(h_top=h_top,T0=T0(itemp))
        call scalar_generator%get_scalar_field(f, domain%mesh_p, 0)
        call grad_generator%get_vector_field(fx_true, fy_true, fz_true, &
                                             domain%mesh_u, domain%mesh_v, domain%mesh_w, &
                                             0,"covariant")

        call grad_3d%calc_grad(fx, fy, fz, f, domain)

        call fx%update(-1.0_8, fx_true, domain%mesh_u)
        call fy%update(-1.0_8, fy_true, domain%mesh_v)
        call fz%update(-1.0_8, fz_true, domain%mesh_w)
        errs%values(2*itemp-1) = fx%maxabs(domain%mesh_u, domain%parcomm)  + &
                                 fy%maxabs(domain%mesh_v, domain%parcomm)  + &
                                 fz%maxabs(domain%mesh_w, domain%parcomm)
        errs%values(2*itemp) = l2norm(fx, domain%mesh_u, domain%parcomm) + &
                               l2norm(fy, domain%mesh_v, domain%parcomm) + &
                               l2norm(fz, domain%mesh_w, domain%parcomm)
    end do
end function test_grad_3d

function test_div_3d(Nh, Nz, hor_div_name, diff_eta_name, &
                                  horizontal_staggering, vertical_staggering) result(errs)

    use const_mod, only : Earth_radii

    use div_3d_factory_mod,   only : create_div_3d_operator
    use abstract_div_3d_mod,  only : div_3d_operator_t
    use config_domain_mod,    only : config_domain_t
    use div3d_test_field_mod, only : div3d_test_wind_t, div3d_test_field_t

    integer(kind=4),  intent(in) :: Nh, Nz
    character(len=*), intent(in) :: hor_div_name, diff_eta_name
    character(len=*), intent(in) :: horizontal_staggering, vertical_staggering

    real(kind=8), parameter :: h_top = 30e3

    type(key_value_r8_t) :: errs

    class(div_3d_operator_t), allocatable :: div_3d

    type(config_domain_t) :: config_domain
    type(domain_t)        :: domain

    type(grid_field_t) :: div, u, v, w, div_true
    type(div3d_test_wind_t)  :: wind_generator
    type(div3d_test_field_t) :: div_generator

    integer(kind=4) :: halo_width = 8

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

    call create_div_3d_operator(div_3d, domain, hor_div_name, diff_eta_name)

    call create_grid_field(u, halo_width, 0, domain%mesh_u)
    call create_grid_field(v, halo_width, 0, domain%mesh_v)
    call create_grid_field(w, halo_width, 0, domain%mesh_w)

    call create_grid_field(div,      0, 0, domain%mesh_p)
    call create_grid_field(div_true, 0, 0, domain%mesh_p)

    wind_generator = div3d_test_wind_t(h_top=h_top)
    div_generator  = div3d_test_field_t(h_top=h_top)

    call wind_generator%get_vector_field(u, v, w, &
                                         domain%mesh_u, domain%mesh_v, domain%mesh_w, &
                                         0,"contravariant")
    call div_generator%get_scalar_field(div_true, domain%mesh_p, 0)

    call div_3d%calc_div(div, u, v, w, domain)

    allocate(errs%keys(2), errs%values(2))
    errs%keys(1)%str = "l_inf"
    errs%keys(2)%str = "l2"

    call div%update(-1.0_8, div_true, domain%mesh_p)
    errs%values(1) = div%maxabs(domain%mesh_p, domain%parcomm) / div_true%maxabs(domain%mesh_p, domain%parcomm)
    errs%values(2) = l2norm(div, domain%mesh_p, domain%parcomm) / l2norm(div_true, domain%mesh_p, domain%parcomm)
end function test_div_3d

type(key_value_r8_t) function test_co2contra_3d(Nh, nz, co2contra_3d_oper_name, &
                            horizontal_staggering, vertical_staggering) result(errs)

    use test_fields_mod,    only : set_vector_test_field, set_scalar_test_field, &
                                   vec_field_gen => cross_polar_flow_generator

    use abstract_co2contra_3d_mod, only : co2contra_3d_operator_t
    use co2contra_3d_factory_mod,  only : create_co2contra_3d_operator
    use mesh_mod,                  only : mesh_t
    use config_domain_mod,         only : config_domain_t
    use const_mod,                 only : Earth_radii
    use grad3d_test_field_mod,     only : grad3d_test_out_t

    integer(kind=4),  intent(in) :: Nh, nz
    character(len=*), intent(in) :: co2contra_3d_oper_name, horizontal_staggering, &
                                    vertical_staggering
    !locals:
    integer(kind=4), parameter  :: ex_halo_width = 8
    type(grid_field_t)          :: u_cov,  v_cov,  w_cov,  &
                                   u_test, v_test, w_test, &
                                   u_true, v_true, w_true

    type(domain_t)        :: domain
    type(config_domain_t) :: config_domain

    class(co2contra_3d_operator_t), allocatable :: co2contra_3d_op
    type(grad3d_test_out_t)   :: vec_generator
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

    call create_grid_field(u_cov, ex_halo_width, 0, domain%mesh_u)
    call create_grid_field(v_cov, ex_halo_width, 0, domain%mesh_v)
    call create_grid_field(w_cov, ex_halo_width, 0, domain%mesh_w)

    call create_grid_field(u_test, ex_halo_width, 0, domain%mesh_u)
    call create_grid_field(v_test, ex_halo_width, 0, domain%mesh_v)
    call create_grid_field(w_test, ex_halo_width, 0, domain%mesh_w)

    call create_grid_field(u_true, ex_halo_width, 0, domain%mesh_u)
    call create_grid_field(v_true, ex_halo_width, 0, domain%mesh_v)
    call create_grid_field(w_true, ex_halo_width, 0, domain%mesh_w)

    vec_generator = grad3d_test_out_t(h_top=h_top,T0=T0)
    call vec_generator%get_vector_field(u_cov, v_cov, w_cov,                          &
                                         domain%mesh_u, domain%mesh_v, domain%mesh_w, &
                                         0,"covariant")
    call vec_generator%get_vector_field(u_true, v_true, w_true,                      &
                                        domain%mesh_u, domain%mesh_v, domain%mesh_w, &
                                        0,"contravariant")

    call create_co2contra_3d_operator(co2contra_3d_op, domain, co2contra_3d_oper_name)

    call co2contra_3d_op%transform(u_test, v_test, w_test, &
                                   u_cov,  v_cov,  w_cov, domain)

    allocate(errs%keys(2), errs%values(2))
    errs%keys(1)%str = "l_inf"
    errs%keys(2)%str = "l2"

    call u_test%update(-1.0_8, u_true, domain%mesh_u)
    call v_test%update(-1.0_8, v_true, domain%mesh_v)
    call w_test%update(-1.0_8, w_true, domain%mesh_w)

    errs%values(1) = (u_test%maxabs(domain%mesh_u, domain%parcomm) + &
                      v_test%maxabs(domain%mesh_v, domain%parcomm) + &
                      w_test%maxabs(domain%mesh_w, domain%parcomm)) / &
                     (u_true%maxabs(domain%mesh_u, domain%parcomm) + &
                      v_true%maxabs(domain%mesh_v, domain%parcomm) + &
                      w_true%maxabs(domain%mesh_w, domain%parcomm))
    errs%values(2) = (l2norm(u_test, domain%mesh_u, domain%parcomm) + &
                      l2norm(v_test, domain%mesh_v, domain%parcomm) + &
                      l2norm(w_test, domain%mesh_w, domain%parcomm)) / &
                     (l2norm(u_true, domain%mesh_u, domain%parcomm) + &
                      l2norm(v_true, domain%mesh_v, domain%parcomm) + &
                      l2norm(w_true, domain%mesh_w, domain%parcomm))
end function test_co2contra_3d

function test_w2uv_interp(Nh, Nz, w2uv_interpolator_name, w2uv_hor_part_name, w2uv_vert_part_name,&
                                  horizontal_staggering, vertical_staggering)   result(errs)

    use interpolator_w2uv_factory_mod,  only : create_w2uv_interpolator
    use abstract_interpolators3d_mod,   only : interpolator_w2uv_t
    use config_domain_mod,              only : config_domain_t

    use grad3d_test_field_mod,       only : grad3d_test_input_t

    real(kind=8), parameter :: h_top = 30e3_8

    integer(kind=4),  intent(in) :: Nh, Nz
    character(len=*), intent(in) :: w2uv_interpolator_name, w2uv_hor_part_name, w2uv_vert_part_name
    character(len=*), intent(in) :: horizontal_staggering, vertical_staggering

    type(key_value_r8_t) :: errs

    class(interpolator_w2uv_t), allocatable :: w2uv_operator

    type(config_domain_t) :: config_domain
    type(domain_t)        :: domain

    type(grid_field_t) :: w, wu, wv, wu_true, wv_true

    type(grad3d_test_input_t) :: scalar_field = grad3d_test_input_t(h_top=h_top,T0=300.0_8)

    integer(kind=4) :: halo_width = 8

    config_domain%N  = nh
    config_domain%Nz = nz
    config_domain%staggering_type     = horizontal_staggering
    config_domain%vertical_staggering = vertical_staggering
    config_domain%metric_type         = "shallow_atmosphere_metric"
    config_domain%topology_type       = "cube"
    config_domain%h_top = h_top
    call config_domain%config_metric%set_defaults()
    config_domain%config_metric%vertical_scale = h_top

    call create_domain(domain, config_domain)

    call create_w2uv_interpolator(w2uv_operator,w2uv_interpolator_name, &
                                  w2uv_hor_part_name, w2uv_vert_part_name, domain)

    call create_grid_field(w, halo_width, 0,  domain%mesh_w)
    call create_grid_field(wu, halo_width, 0, domain%mesh_u)
    call create_grid_field(wv, halo_width, 0, domain%mesh_v)
    call create_grid_field(wu_true, halo_width, 0, domain%mesh_u)
    call create_grid_field(wv_true, halo_width, 0, domain%mesh_v)

    call scalar_field%get_scalar_field(w, domain%mesh_w, 0)
    call scalar_field%get_scalar_field(wu_true, domain%mesh_u, 0)
    call scalar_field%get_scalar_field(wv_true, domain%mesh_v, 0)

    call w2uv_operator%interp_w2uv(wu,wv,w,domain)

    allocate(errs%keys(2), errs%values(2))
    errs%keys(1)%str = "vertical_ExnerP C_norm"
    errs%keys(2)%str = "vertical_ExnerP L2_norm"

    call wu%update(-1.0_8, wu_true, domain%mesh_u)
    call wv%update(-1.0_8, wv_true, domain%mesh_v)
    errs%values(1) = (wu%maxabs(domain%mesh_u, domain%parcomm)+&
                     wv%maxabs(domain%mesh_v, domain%parcomm)) / &
                     (wu_true%maxabs(domain%mesh_u, domain%parcomm)+&
                     wv_true%maxabs(domain%mesh_v, domain%parcomm))
    errs%values(2) = (l2norm(wu, domain%mesh_u, domain%parcomm)+&
                     l2norm(wv, domain%mesh_v, domain%parcomm)) / &
                     (l2norm(wu_true, domain%mesh_u, domain%parcomm)+&
                     l2norm(wv_true, domain%mesh_v, domain%parcomm))

end function test_w2uv_interp

function test_uv2w_interp(Nh, Nz, uv2w_interpolator_name, uv2w_hor_part_name, uv2w_vert_part_name, &
                                  horizontal_staggering, vertical_staggering) result(errs)

    use interpolator_uv2w_factory_mod,  only : create_uv2w_interpolator
    use abstract_interpolators3d_mod,   only : interpolator_uv2w_t
    use config_domain_mod,              only : config_domain_t

    use grad3d_test_field_mod,   only : grad3d_test_out_t

    real(kind=8), parameter :: h_top = 30e3_8

    integer(kind=4),  intent(in) :: Nh, Nz
    character(len=*), intent(in) :: uv2w_interpolator_name, uv2w_hor_part_name, &
                                    uv2w_vert_part_name
    character(len=*), intent(in) :: horizontal_staggering, vertical_staggering

    type(key_value_r8_t) :: errs

    class(interpolator_uv2w_t), allocatable :: uv2w_operator

    type(config_domain_t) :: config_domain
    type(domain_t)        :: domain

    type(grid_field_t) :: u, v, uw, vw, uw_true, vw_true

    type(grad3d_test_out_t) :: vec_field_gen = grad3d_test_out_t(h_top=h_top,T0=300.0_8)

    integer(kind=4) :: halo_width = 8

    config_domain%N  = nh
    config_domain%Nz = nz
    config_domain%staggering_type     = horizontal_staggering
    config_domain%vertical_staggering = vertical_staggering
    config_domain%metric_type         = "shallow_atmosphere_metric"
    config_domain%topology_type       = "cube"
    config_domain%h_top = h_top
    call config_domain%config_metric%set_defaults()
    config_domain%config_metric%vertical_scale = h_top

    call create_domain(domain, config_domain)

    call create_uv2w_interpolator(uv2w_operator,uv2w_interpolator_name, &
                                  uv2w_hor_part_name, uv2w_vert_part_name, domain)

    call create_grid_field(u, halo_width, 0, domain%mesh_u)
    call create_grid_field(v, halo_width, 0, domain%mesh_v)
    call create_grid_field(uw, halo_width, 0,  domain%mesh_w)
    call create_grid_field(vw, halo_width, 0,  domain%mesh_w)
    call create_grid_field(uw_true, halo_width, 0, domain%mesh_w)
    call create_grid_field(vw_true, halo_width, 0, domain%mesh_w)

    call vec_field_gen%get_vector_field(u,v,uw,domain%mesh_u,domain%mesh_v,&
                                        domain%mesh_w,0,"contravariant")
    call vec_field_gen%get_vector_field(uw_true,vw_true,uw,domain%mesh_w,&
                                        domain%mesh_w,domain%mesh_w,0,"contravariant")

    call uv2w_operator%interp_uv2w(uw,vw,u,v,domain)

    allocate(errs%keys(2), errs%values(2))
    errs%keys(1)%str = "vertical_ExnerP C_norm"
    errs%keys(2)%str = "vertical_ExnerP L2_norm"

    call uw%update(-1.0_8, uw_true, domain%mesh_w)
    call vw%update(-1.0_8, vw_true, domain%mesh_w)
    errs%values(1) = (uw%maxabs(domain%mesh_w, domain%parcomm)+&
                      vw%maxabs(domain%mesh_w, domain%parcomm)) / &
                     (uw_true%maxabs(domain%mesh_w, domain%parcomm)+&
                      vw_true%maxabs(domain%mesh_w, domain%parcomm))
    errs%values(2) = (l2norm(uw, domain%mesh_w, domain%parcomm)+&
                      l2norm(vw, domain%mesh_w, domain%parcomm)) / &
                     (l2norm(uw_true, domain%mesh_w, domain%parcomm)+&
                      l2norm(vw_true, domain%mesh_w, domain%parcomm))

end function test_uv2w_interp

function test_scalar_advection_3d(Nh, Nz, advection_oper_name, config_str, points_type, &
                                          horizontal_staggering, vertical_staggering) result(errs)
    use const_mod, only : Earth_radii

    use abstract_scalar_advection3d_mod, only : scalar_advection3d_t
    use scalar_advection_factory_mod,    only : create_scalar_advection3d_operator
    use config_advection_3d_mod,         only : get_advection_3d_config
    use config_domain_mod,               only : config_domain_t

    use mesh_mod,                        only : mesh_t

    use test_fields_3d_mod,              only : scalar_field3d_t, vector_field3d_t
    use grad3d_test_field_mod,           only : grad3d_test_input_t, grad3d_test_out_t

    use solid_rotation_wind_field_mod,   only : solid_rotation_wind_field_t

    integer(kind=4),  intent(in) :: Nh, Nz
    character(len=*), intent(in) :: advection_oper_name, config_str, points_type, &
                                    horizontal_staggering, vertical_staggering

    type(key_value_r8_t) :: errs

    class(config_t), allocatable             :: config
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
    real(kind=8), parameter :: h_top = 30e3, T0(2) = [300._8,1e10_8]
    integer(kind=4) :: itemp

    call get_advection_3d_config(config,advection_oper_name)
    call config%parse(config_str)

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

    allocate(errs%keys(4), errs%values(4))
    errs%keys(1)%str = "scalar advection 3d L_inf"
    errs%keys(2)%str = "scalar advection 3d L2"
    errs%keys(3)%str = "scalar advection 3d quasi-hor L_inf"
    errs%keys(4)%str = "scalar advection 3d quasi-hor L2"

    wind_generator   = solid_rotation_wind_field_t(w_max = 0.1_8, u0 = 1.0_8, t= 0.0_8)
    call wind_generator%get_vector_field(u,v,eta_dot,domain%mesh_u,domain%mesh_v, &
                                         domain%mesh_w,0, "contravariant")
    call wind_generator%get_vector_field(u_p,v_p,eta_dot_p, mesh, mesh, mesh, &
                                         0, "contravariant")

    call create_scalar_advection3d_operator(advection_op,advection_oper_name, &
                                            config,domain)

    do itemp = 1, size(T0)

        scalar_generator = grad3d_test_input_t(h_top=h_top,T0=T0(itemp))
        grad_generator   = grad3d_test_out_t(h_top=h_top,T0=T0(itemp))

        call scalar_generator%get_scalar_field(f, mesh, 0)
        call grad_generator%get_vector_field(fx, fy, fz, mesh, mesh, mesh, &
                                            0,"covariant")
        call advection_op%calc_adv3d(v_nabla_f,f,u,v,eta_dot,domain)

        call fx%assign_prod(-1.0_8,fx,u_p,mesh)
        call fy%assign_prod(-1.0_8,fy,v_p,mesh)
        call fz%assign_prod(-1.0_8,fz,eta_dot_p,mesh)

        call f%assign(1.0_8,fx,mesh)
        call f%update(1.0_8,fy,1.0_8,fz,mesh)
        call v_nabla_f%update(-1.0_8,f,mesh)

        errs%values(2*itemp-1) = v_nabla_f%maxabs(mesh, domain%parcomm) / f%maxabs(mesh, domain%parcomm)
        errs%values(2*itemp)   = l2norm(v_nabla_f, mesh, domain%parcomm) / l2norm(f, mesh, domain%parcomm)
    end do

end function test_scalar_advection_3d

end module test_diffops_3d_mod
