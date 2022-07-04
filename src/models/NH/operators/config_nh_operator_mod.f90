module config_nh_operator_mod

use config_mod,                  only : config_t
use config_advection_3d_mod,     only : get_advection_3d_config
use config_mixvec_transform_mod, only : config_mixvec_transform_t
use parcomm_mod,                 only : parcomm_global

implicit none

type, extends(config_t) :: config_Ptheta_linear_t
    real(kind=8) :: Nb   = 0.01_8
    real(kind=8) :: T0   = 300.0_8
    character(:), allocatable :: background_type, &
                                 grad_hor_part_name, grad_vert_part_name, &
                                 div_hor_part_name,  div_vert_part_name, &
                                 co2contra_operator_name, w2p_operator_name
    contains
    procedure :: parse => parse_Ptheta_linear_config
end type config_Ptheta_linear_t

type, extends(config_t) :: config_advection3d_t
    character(:),    allocatable :: p_advection_oper_name
    character(:),    allocatable :: theta_advection_oper_name
    class(config_t), allocatable :: config_p_advec, config_theta_advec
    character(:),    allocatable :: wind_field
    contains
    procedure :: parse => parse_advection3d_config
end type config_advection3d_t

type, extends(config_t) :: config_nonlin_nh_operator_t
    character(:),    allocatable :: grad_hor_part_name, grad_vert_part_name,        &
                                    div_hor_part_name,  div_vert_part_name,         &
                                    mixvec_transform_name,                          &
                                    theta2uv_operator_name, theta2uv_hor_part_name, &
                                    theta2uv_vert_part_name
    character(:),    allocatable :: p_advection_oper_name
    character(:),    allocatable :: theta_advection_oper_name
    class(config_t), allocatable :: config_mixvec_transform, config_p_advec, &
                                    config_theta_advec
    character(:),    allocatable :: vec_adv_op_name
    character(:),    allocatable :: uv_hor_adv_op_name, uv_ver_adv_op_name
    character(:),    allocatable :: w_adv_op_name, w_adv_hor_part_name, w_adv_ver_part_name
    class(config_t), allocatable :: config_momentum_advection_operator
    character(:),    allocatable :: coriolis_op_name

    contains
    procedure :: parse => parse_nonlinear_nh_operator_config
end type config_nonlin_nh_operator_t

contains

subroutine get_nh_operator_config(operator_config, operator_type)
    character(len=*), intent(in) :: operator_type
    class(config_t), allocatable, intent(out) :: operator_config

    select case(operator_type)
    case("nonlinear_nh")
        allocate(config_nonlin_nh_operator_t :: operator_config)
    case("Ptheta_linear")
        allocate(config_Ptheta_linear_t :: operator_config)
    case("advection_3d")
        allocate(config_advection3d_t :: operator_config)
    case default
        call parcomm_global%abort("get_nh_operator_config, unknown nh operator type: "// operator_type)
    end select
end subroutine get_nh_operator_config

subroutine parse_Ptheta_linear_config(this,config_string)
    class(config_Ptheta_linear_t), intent(inout) :: this
    character(len=*), intent(in) :: config_string

    namelist /Ptheta_linear_nh_operator/ background_type, Nb, T0,&
                                         grad_hor_part_name, grad_vert_part_name, &
                                         div_hor_part_name,  div_vert_part_name, &
                                         co2contra_operator_name, w2p_operator_name

    character(len=256) :: background_type, grad_hor_part_name, grad_vert_part_name, &
                          div_hor_part_name, div_vert_part_name, &
                          co2contra_operator_name, w2p_operator_name
    real(kind=8) :: Nb=0.01, T0=300.0

    read(config_string,Ptheta_linear_nh_operator)

    this%Nb = Nb
    this%T0 = T0

    this%background_type          = trim(background_type)
    this%grad_hor_part_name       = trim(grad_hor_part_name)
    this%grad_vert_part_name      = trim(grad_vert_part_name)
    this%div_hor_part_name        = trim(div_hor_part_name)
    this%div_vert_part_name       = trim(div_vert_part_name)
    this%co2contra_operator_name  = trim(co2contra_operator_name)
    this%w2p_operator_name        = trim(w2p_operator_name)

end subroutine parse_Ptheta_linear_config

subroutine parse_advection3d_config(this,config_string)
    class(config_advection3d_t), intent(inout) :: this
    character(len=*), intent(in) :: config_string

    namelist /advection3d_operator/  p_advection_oper_name,      &
                                     p_advection_config_str,     &
                                     theta_advection_oper_name,  &
                                     theta_advection_config_str, &
                                     wind_field

    character(len=512) ::  p_advection_oper_name,      &
                           p_advection_config_str,     &
                           theta_advection_oper_name,  &
                           theta_advection_config_str, &
                           wind_field

    read(config_string,advection3d_operator)

    this%p_advection_oper_name         = trim(p_advection_oper_name)
    call get_advection_3d_config(this%config_p_advec, this%p_advection_oper_name)
    call this%config_p_advec%parse(trim(p_advection_config_str))
    this%theta_advection_oper_name     = trim(theta_advection_oper_name)
    call get_advection_3d_config(this%config_theta_advec, this%theta_advection_oper_name)
    call this%config_theta_advec%parse(trim(theta_advection_config_str))
    this%wind_field                    = trim(wind_field)

end subroutine parse_advection3d_config

subroutine parse_nonlinear_nh_operator_config(this,config_string)
    class(config_nonlin_nh_operator_t), intent(inout) :: this
    character(len=*), intent(in) :: config_string

    namelist /nonlin_nh_operator/ grad_hor_part_name, grad_vert_part_name, &
                                  div_hor_part_name,  div_vert_part_name,  &
                                  mixvec_transform_name,                   &
                                  mixvec_transform_config_str,             &
                                  theta2uv_operator_name,                  &
                                  theta2uv_hor_part_name,                  &
                                  theta2uv_vert_part_name,                 &
                                  p_advection_oper_name,                   &
                                  p_advection_config_str,                  &
                                  theta_advection_oper_name,               &
                                  theta_advection_config_str,              &
                                  vec_adv_op_name,                         &
                                  vec_adv_oper_config_str,                 &
                                  coriolis_op_name

    character(len=512) :: grad_hor_part_name, grad_vert_part_name, &
                          div_hor_part_name, div_vert_part_name,   &
                          mixvec_transform_name,                   &
                          mixvec_transform_config_str,             &
                          theta2uv_operator_name,                  &
                          theta2uv_hor_part_name,                  &
                          theta2uv_vert_part_name,                 &
                          p_advection_oper_name,                   &
                          p_advection_config_str,                  &
                          theta_advection_oper_name,               &
                          theta_advection_config_str,              &
                          vec_adv_op_name,                         &
                          coriolis_op_name
    character(len=1024) :: vec_adv_oper_config_str

    read(config_string,nonlin_nh_operator)

    this%grad_hor_part_name            = trim(grad_hor_part_name)
    this%grad_vert_part_name           = trim(grad_vert_part_name)
    this%div_hor_part_name             = trim(div_hor_part_name)
    this%div_vert_part_name            = trim(div_vert_part_name)
    this%mixvec_transform_name         = trim(mixvec_transform_name)
    this%config_mixvec_transform       = config_mixvec_transform_t()
    call this%config_mixvec_transform%parse(mixvec_transform_config_str)
    this%theta2uv_operator_name        = trim(theta2uv_operator_name)
    this%theta2uv_hor_part_name        = trim(theta2uv_hor_part_name)
    this%theta2uv_vert_part_name       = trim(theta2uv_vert_part_name)
    this%p_advection_oper_name         = trim(p_advection_oper_name)
    call get_advection_3d_config(this%config_p_advec, this%p_advection_oper_name)
    call this%config_p_advec%parse(trim(p_advection_config_str))
    this%theta_advection_oper_name     = trim(theta_advection_oper_name)
    call get_advection_3d_config(this%config_theta_advec, this%theta_advection_oper_name)
    call this%config_theta_advec%parse(trim(theta_advection_config_str))
    this%vec_adv_op_name               = trim(vec_adv_op_name)
    call get_advection_3d_config(this%config_momentum_advection_operator, this%vec_adv_op_name)
    call this%config_momentum_advection_operator%parse(trim(vec_adv_oper_config_str))
    this%coriolis_op_name              = trim(coriolis_op_name)

end subroutine parse_nonlinear_nh_operator_config

end module config_nh_operator_mod
