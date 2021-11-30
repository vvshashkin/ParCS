module config_nh_operator_mod

use config_mod,  only : config_t
use parcomm_mod, only : parcomm_global

implicit none

type, extends(config_t) :: config_Ptheta_linear_t
    real(kind=8) :: Nb   = 0.01_8
    real(kind=8) :: T0   = 300.0_8
    character(:), allocatable :: background_type, &
                                 grad_hor_part_name, grad_vert_part_name, &
                                 div_hor_part_name,  div_vert_part_name, &
                                 co2contra_operator_name, &
                                 w2uv_operator_name, w2p_operator_name
    contains
    procedure :: parse => parse_Ptheta_linear_config
end type config_Ptheta_linear_t

type, extends(config_t) :: config_advection3d_t
    character(:), allocatable :: p_advection_oper_name
    character(:), allocatable :: p_hor_advection_oper_name
    character(:), allocatable :: p_z_advection_oper_name
    character(:), allocatable :: wind_field
    contains
    procedure :: parse => parse_advection3d_config
end type config_advection3d_t

contains

function get_nh_operator_config(operator_type) result(operator_config)
    character(len=*), intent(in) :: operator_type

    class(config_t), allocatable :: operator_config

    select case(operator_type)
    case("Ptheta_linear")
        operator_config = config_Ptheta_linear_t()
    case("advection_3d")
        operator_config = config_advection3d_t()
    case default
        call parcomm_global%abort("get_nh_operator_config, unknown nh operator type: "// operator_type)
    end select
end function get_nh_operator_config

subroutine parse_Ptheta_linear_config(this,config_string)
    class(config_Ptheta_linear_t), intent(inout) :: this
    character(len=*), intent(in) :: config_string

    namelist /Ptheta_linear_nh_operator/ background_type, Nb, T0,&
                                         grad_hor_part_name, grad_vert_part_name, &
                                         div_hor_part_name,  div_vert_part_name, &
                                         co2contra_operator_name, &
                                         w2uv_operator_name, w2p_operator_name

    character(len=256) :: background_type, grad_hor_part_name, grad_vert_part_name, &
                          div_hor_part_name, div_vert_part_name, &
                          co2contra_operator_name, w2uv_operator_name, w2p_operator_name
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
    this%w2uv_operator_name       = trim(w2uv_operator_name)
    this%w2p_operator_name        = trim(w2p_operator_name)

end subroutine parse_Ptheta_linear_config

subroutine parse_advection3d_config(this,config_string)
    class(config_advection3d_t), intent(inout) :: this
    character(len=*), intent(in) :: config_string

    namelist /advection3d_operator/  p_advection_oper_name,     &
                                     p_hor_advection_oper_name, &
                                     p_z_advection_oper_name,   &
                                     wind_field

    character(len=256) ::  p_advection_oper_name,     &
                           p_hor_advection_oper_name, &
                           p_z_advection_oper_name,   &
                           wind_field

    read(config_string,advection3d_operator)

    this%p_advection_oper_name     = trim(p_advection_oper_name)
    this%p_hor_advection_oper_name = trim(p_hor_advection_oper_name)
    this%p_z_advection_oper_name   = trim(p_z_advection_oper_name)
    this%wind_field                = trim(wind_field)

end subroutine parse_advection3d_config

end module config_nh_operator_mod
