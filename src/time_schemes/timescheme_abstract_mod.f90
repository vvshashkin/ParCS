module timescheme_abstract_mod

use stvec_abstract_mod,    only : stvec_abstract_t
use operator_abstract_mod, only : operator_abstract_t

implicit none

type, abstract, public :: timescheme_abstract_t
    contains
    procedure(step),  deferred :: step
end type timescheme_abstract_t

abstract interface
    subroutine step(this, operator, v0, dt)
        import stvec_abstract_t
        import operator_abstract_t
        import timescheme_abstract_t

        class(timescheme_abstract_t),    intent(inout) :: this
        class(operator_abstract_t),      intent(in)    :: operator
        class(stvec_abstract_t), target, intent(inout) :: v0
        real(kind=8),                    intent(in)    :: dt
    end subroutine step
end interface

contains

end module timescheme_abstract_mod
