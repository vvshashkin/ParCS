module nh_model_mod

use domain_mod,     only : domain_t
use timescheme_mod, only : timescheme_t
use stvec_mod,      only : stvec_t

use abstract_postprocessing_mod, only : postprocessing_t

implicit none

type nh_model_t

    type(domain_t) :: domain
    class(timescheme_t), allocatable :: timescheme
    class(stvec_t),      allocatable :: stvec

    class(postprocessing_t), allocatable :: postproc

    real(kind=8) :: dt
    real(kind=8) :: tau_write
    real(kind=8) :: tau_diagnostics
    real(kind=8) :: simulation_time

    contains
    procedure, public :: run

end type nh_model_t

contains

subroutine run(this)
    class(nh_model_t), intent(inout) :: this

    call this%postproc%write_fields(1,this%stvec,this%domain)
end subroutine run

end module nh_model_mod
