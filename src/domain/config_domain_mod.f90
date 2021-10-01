module config_domain_mod

use config_mod, only : config_t

implicit none

type, public, extends(config_t) :: config_domain_t

    integer(kind=4) :: N, Nz

    character(len=:), allocatable :: staggering_type
    character(len=:), allocatable :: topology_type
    character(len=:), allocatable :: metric_type

contains
    procedure, public :: parse
end type config_domain_t

contains

subroutine parse(this, config_string)

    use parcomm_mod, only : parcomm_global

    class(config_domain_t), intent(inout) :: this
    character(len=*),       intent(in)    :: config_string

    integer(kind=4) :: N, Nz

    character(len=255) :: staggering_type
    character(len=255) :: topology_type
    character(len=255) :: metric_type

    namelist /domain/ N, Nz, staggering_type, topology_type, metric_type

    read(config_string, domain)

    this%N = N
    this%Nz = Nz

    this%staggering_type = trim(staggering_type)
    this%topology_type   = trim(topology_type)
    this%metric_type     = trim(metric_type)

end subroutine parse

end module config_domain_mod
