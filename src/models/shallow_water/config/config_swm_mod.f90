module config_swm_mod

use config_mod, only : config_t

implicit none

type, public, extends(config_t) :: config_swm_t
    integer(kind=4) :: N, Nz

    character(:), allocatable :: staggering_type
    character(:), allocatable :: topology_type

    character(:), allocatable :: div_op_name
    character(:), allocatable :: grad_op_name
    character(:), allocatable :: curl_op_name
    character(:), allocatable :: coriolis_op_name
    character(:), allocatable :: KE_op_name
    character(:), allocatable :: massflux_op_name

contains
    procedure, public :: parse
end type config_swm_t

contains

subroutine parse(this, config_filename)

    use namelist_read_mod, only : read_namelist_as_str
    use parcomm_mod,       only : parcomm_global

    class(config_swm_t),  intent(inout) :: this
    character(len=*),     intent(in)    :: config_filename

    integer(kind=4) :: N, Nz

    character(len=255) :: staggering_type
    character(len=255) :: topology_type

    character(len=255) :: div_op_name
    character(len=255) :: grad_op_name
    character(len=255) :: curl_op_name
    character(len=255) :: coriolis_op_name
    character(len=255) :: massflux_op_name
    character(len=255) :: KE_op_name

    namelist /shallow_water_model/ N, Nz, staggering_type, topology_type,   &
                                   div_op_name, grad_op_name, curl_op_name, &
                                   coriolis_op_name, KE_op_name, massflux_op_name
    character(:), allocatable :: namelist_string

    call read_namelist_as_str(namelist_string, config_filename, parcomm_global%myid)

    read(namelist_string, shallow_water_model)

    this%staggering_type = trim(staggering_type)
    this%topology_type   = trim(topology_type)

    this%div_op_name      = trim(div_op_name)
    this%grad_op_name     = trim(grad_op_name)
    this%curl_op_name     = trim(curl_op_name)
    this%coriolis_op_name = trim(coriolis_op_name)
    this%KE_op_name       = trim(KE_op_name)
    this%massflux_op_name = trim(massflux_op_name)

    this%N = N
    this%Nz = Nz

end subroutine parse

end module config_swm_mod
