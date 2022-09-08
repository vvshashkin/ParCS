module config_domain_mod

use config_mod,           only : config_t
use config_metric_mod,    only : config_metric_t
use config_orography_mod, only : config_test_orography_t

implicit none

type, public, extends(config_t) :: config_domain_t

    integer(kind=4) :: N, Nz

    character(len=:), allocatable :: staggering_type
    character(len=:), allocatable :: vertical_staggering
    character(len=:), allocatable :: topology_type
    character(len=:), allocatable :: metric_type
    real(kind=8)                  :: h_top = 1.0_8
    logical                       :: is_orographic_curvilinear = .false.
    character(len=:), allocatable :: orography_name
    class(config_t),  allocatable :: config_orography
    integer(kind=4) :: halo_width = 8

    type(config_metric_t) :: config_metric

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
    character(len=255) :: vertical_staggering = "None"
    real(kind=8)       :: h_top = 1.0_8
    logical            :: is_orographic_curvilinear
    character(len=255) :: orography_name, orography_config_str

    namelist /domain/ N, Nz, staggering_type, topology_type, metric_type,&
                             vertical_staggering, h_top, is_orographic_curvilinear, &
                             orography_name, orography_config_str

    read(config_string, domain)

    this%N = N
    this%Nz = Nz

    this%staggering_type = trim(staggering_type)
    this%topology_type   = trim(topology_type)
    this%metric_type     = trim(metric_type)
    this%vertical_staggering     = trim(vertical_staggering)

    this%h_top = h_top
    this%is_orographic_curvilinear = is_orographic_curvilinear
    this%config_orography = config_test_orography_t()
    this%orography_name   = trim(orography_name)
    if(is_orographic_curvilinear) &
        call this%config_orography%parse(orography_config_str)

    call this%config_metric%parse(config_string)

end subroutine parse

end module config_domain_mod
