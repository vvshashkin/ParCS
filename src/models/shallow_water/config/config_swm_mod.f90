module config_swm_mod

use config_mod,        only : config_t
use config_domain_mod, only : config_domain_t

implicit none

type, public, extends(config_t) :: config_swm_t

    type(config_domain_t) :: config_domain

    character(:), allocatable :: div_op_name
    character(:), allocatable :: grad_op_name
    character(:), allocatable :: curl_op_name
    character(:), allocatable :: coriolis_op_name
    character(:), allocatable :: KE_op_name
    character(:), allocatable :: massflux_op_name
    character(:), allocatable :: co2contra_op_name

contains
    procedure, public :: parse
end type config_swm_t

contains

subroutine parse(this, config_string)

    use parcomm_mod,       only : parcomm_global

    class(config_swm_t),  intent(inout) :: this
    character(len=*),     intent(in)    :: config_string

    character(len=255) :: div_op_name
    character(len=255) :: grad_op_name
    character(len=255) :: curl_op_name
    character(len=255) :: coriolis_op_name
    character(len=255) :: KE_op_name
    character(len=255) :: massflux_op_name
    character(len=255) :: co2contra_op_name

    namelist /shallow_water_model/ div_op_name, grad_op_name, curl_op_name, &
                                   coriolis_op_name, KE_op_name, massflux_op_name, &
                                   co2contra_op_name

    read(config_string, shallow_water_model)

    call this%config_domain%parse(config_string)

    this%div_op_name       = trim(div_op_name)
    this%grad_op_name      = trim(grad_op_name)
    this%curl_op_name      = trim(curl_op_name)
    this%coriolis_op_name  = trim(coriolis_op_name)
    this%KE_op_name        = trim(KE_op_name)
    this%massflux_op_name  = trim(massflux_op_name)
    this%co2contra_op_name = trim(co2contra_op_name)

end subroutine parse

end module config_swm_mod
