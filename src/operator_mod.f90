module operator_mod

use stvec_mod,      only : stvec_t
use domain_mod,     only : domain_t
use parcomm_mod,    only : parcomm_global

implicit none

private

type, abstract, public :: operator_t

contains

    procedure(apply_to),  deferred :: apply_to !vout=A*vin
    procedure, public              :: solve    !vout=inverse(I-dt*A)*rhs

end type operator_t

abstract interface
    subroutine apply_to(this, vout, vin, domain)
        import stvec_t
        import operator_t
        import domain_t
        class(operator_t),         intent(inout) :: this
        class(stvec_t),            intent(inout) :: vout !inout to enable preallocated bectors
        class(stvec_t),            intent(in)    :: vin
        class(domain_t),           intent(in)    :: domain
    end subroutine apply_to
end interface

contains

subroutine solve(this, vout, rhs, dt, domain)
    class(operator_t),         intent(inout) :: this
    class(stvec_t),            intent(inout) :: vout !inout to enable preallocated bectors
    class(stvec_t),            intent(in)    :: rhs
    real(kind=8),              intent(in)    :: dt
    class(domain_t),           intent(in)    :: domain

    call parcomm_global%abort("solv function not implemented for specific operator class")
end subroutine solve

end module operator_mod
