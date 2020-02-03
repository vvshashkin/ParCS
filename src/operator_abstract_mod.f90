module operator_abstract_mod

use stvec_abstract_mod, only : stvec_abstract_t

implicit none

type, abstract, public :: operator_abstract_t

contains

    procedure(act),  deferred :: act

end type operator_abstract_t

abstract interface
    subroutine act(this, vout, vin)
        !calculates linear combination alpha*this+beta*other
        import stvec_abstract_t
        import operator_abstract_t
        class(operator_abstract_t),    intent(inout) :: this
        class(stvec_abstract_t),       intent(inout) :: vout !inout to enable preallocated bectors
        class(stvec_abstract_t),       intent(in)    :: vin
    end subroutine act
end interface

contains

end module operator_abstract_mod
