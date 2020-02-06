module timescheme_abstract_mod

use state_abstract_mod,    only : state_abstract_t

implicit none

type, abstract, public :: timescheme_abstract_t
    contains
    procedure(step),  deferred :: step
end type timescheme_abstract_t

abstract interface
    subroutine step(this, v0, dt)
        import state_abstract_t
        import timescheme_abstract_t

        class(timescheme_abstract_t),    intent(inout) :: this
        class(state_abstract_t), target, intent(inout) :: v0
        real(kind=8),                    intent(in)    :: dt
    end subroutine step
end interface

contains

end module timescheme_abstract_mod
