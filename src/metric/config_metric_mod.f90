module config_metric_mod

use config_mod, only : config_t

implicit none

type, public, extends(config_t) :: config_metric_t

    integer(kind=4) :: N, Nz

    real(kind=8) :: scale = 1.0_8
    real(kind=8) :: omega = 1.0_8
    real(kind=8) :: rotation_matrix(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    real(kind=8) :: rotation_axis(3)     = [0, 0, 1]

contains
    procedure, public :: parse
end type config_metric_t

contains

subroutine parse(this, config_string)

    use parcomm_mod, only : parcomm_global

    class(config_metric_t), intent(inout) :: this
    character(len=*),       intent(in)    :: config_string

    real(kind=8) :: scale = 1.0_8
    real(kind=8) :: omega = 1.0_8
    real(kind=8) :: rotation_matrix(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    real(kind=8) :: rotation_axis(3)     = [0, 0, 1]

    namelist /metric/ scale, omega, rotation_matrix, rotation_axis

    read(config_string, metric)

    this%scale = scale
    this%omega = omega
    this%rotation_matrix = rotation_matrix
    this%rotation_axis = rotation_axis

end subroutine parse

end module config_metric_mod
