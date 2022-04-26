module config_advection_3d_mod

use config_mod,  only : config_t
use parcomm_mod, only : parcomm_global

implicit none

type, extends(config_t) :: config_p_advection_t
    character(len=:), allocatable :: hor_advection_oper_name
    character(len=:), allocatable :: z_advection_oper_name
    character(len=:), allocatable :: w2p_operator_name
    character(len=:), allocatable :: uv2p_operator_name
    character(len=:), allocatable :: p_halo
    contains
    procedure :: parse => parse_p_advection_config
end type

type, extends(config_t) :: config_w_advection_t
    character(len=:), allocatable :: hor_advection_oper_name
    character(len=:), allocatable :: z_advection_oper_name
    character(len=:), allocatable :: uv2w_operator_name
    character(len=:), allocatable :: uv2w_hor_part_name
    character(len=:), allocatable :: uv2w_vert_part_name
    character(len=:), allocatable :: w_halo
    contains
    procedure :: parse => parse_w_advection_config
end type

contains

function get_advection_3d_config(advection_3d_operator_name) result(config)
    character(len=*), intent(in) :: advection_3d_operator_name
    class(config_t), allocatable :: config

    select case(advection_3d_operator_name)
    case("advection_p_staggered")
        config = config_p_advection_t()
    case("advection_w_staggered")
        config = config_w_advection_t()
    case default
        call parcomm_global%abort("get_advection_3d_config error: "     // &
                                  "unknown advection_3d_operator name " // &
                                   advection_3d_operator_name)
    end select
end function

subroutine parse_p_advection_config(this, config_string)

    class(config_p_advection_t),  intent(inout) :: this
    character(len=*),             intent(in)    :: config_string

    character(len=256) :: hor_advection_oper_name, z_advection_oper_name, &
                          w2p_operator_name, uv2p_operator_name, p_halo

    namelist /p_advection_conf/ hor_advection_oper_name, z_advection_oper_name, &
                                w2p_operator_name, uv2p_operator_name, p_halo

    read(config_string, p_advection_conf)

    this%hor_advection_oper_name = trim(hor_advection_oper_name)
    this%z_advection_oper_name   = trim(z_advection_oper_name)
    this%w2p_operator_name       = trim(w2p_operator_name)
    this%uv2p_operator_name      = trim(uv2p_operator_name)
    this%p_halo                  = trim(p_halo)

end subroutine parse_p_advection_config

subroutine parse_w_advection_config(this, config_string)

    class(config_w_advection_t),  intent(inout) :: this
    character(len=*),             intent(in)    :: config_string

    character(len=256) :: hor_advection_oper_name, z_advection_oper_name, &
                          uv2w_operator_name, uv2w_hor_part_name, uv2w_vert_part_name, &
                          w_halo

    namelist /w_advection_conf/ hor_advection_oper_name, z_advection_oper_name, &
                                uv2w_operator_name, uv2w_hor_part_name, uv2w_vert_part_name, &
                                w_halo

    read(config_string, w_advection_conf)

    this%hor_advection_oper_name = trim(hor_advection_oper_name)
    this%z_advection_oper_name   = trim(z_advection_oper_name)
    this%uv2w_operator_name      = trim(uv2w_operator_name)
    this%uv2w_hor_part_name      = trim(uv2w_hor_part_name)
    this%uv2w_vert_part_name     = trim(uv2w_vert_part_name)
    this%w_halo                  = trim(w_halo)

end subroutine parse_w_advection_config

end module config_advection_3d_mod
