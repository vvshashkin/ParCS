module operator_abstract_mod

use stvec_abstract_mod,     only : stvec_abstract_t
use container_abstract_mod, only : model_parameters_abstract_t

implicit none

type, abstract, public :: operator_abstract_t

contains

    procedure(act),  deferred :: act
    procedure                 :: solv ! solves (I+dt*A)f = r for f

end type operator_abstract_t

abstract interface
    subroutine act(this, vout, vin, model_params)
        import stvec_abstract_t
        import operator_abstract_t
        import model_parameters_abstract_t
        class(operator_abstract_t),         intent(inout) :: this
        class(stvec_abstract_t),            intent(inout) :: vout !inout to enable preallocated bectors
        class(stvec_abstract_t),            intent(in)    :: vin
        class(model_parameters_abstract_t), intent(in)    :: model_params
    end subroutine act
end interface

contains

subroutine solv(this, dt, vout, rhs, model_params)
    use stvec_abstract_mod,     only : stvec_abstract_t
    use container_abstract_mod, only : model_parameters_abstract_t

    class(operator_abstract_t),         intent(inout) :: this
    real(kind=8),                       intent(in)    :: dt
    class(stvec_abstract_t),            intent(inout) :: vout !inout to enable preallocated bectors
    class(stvec_abstract_t),            intent(in)    :: rhs
    class(model_parameters_abstract_t), intent(in)    :: model_params

    call avost("solv function not implemented for specific operator class")
end subroutine solv

end module operator_abstract_mod
