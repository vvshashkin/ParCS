module explicit_Eul1_mod

use stvec_abstract_mod,      only : stvec_abstract_t
use timescheme_abstract_mod, only : timescheme_abstract_t

implicit none

type, extends(timescheme_abstract_t), public :: explicit_Eul1_t

contains

    procedure, public :: step => step_expl_Eul1

end type explicit_Eul1_t

contains

subroutine step_expl_Eul1(this, operator, v0, dt)

    use operator_abstract_mod, only : operator_abstract_t

    class(explicit_Eul1_t),          intent(inout) :: this
    class(operator_abstract_t),      intent(inout) :: operator
    class(stvec_abstract_t), target, intent(inout) :: v0
    real(kind=8),                    intent(in)    :: dt

    !local
    class(stvec_abstract_t), allocatable :: v

    allocate(v,source=v0) !infer specific type of v0 into v
    call operator%act(v,v0)
    call v0.add(v,1._8,dt)
end subroutine step_expl_Eul1

end module explicit_Eul1_mod
