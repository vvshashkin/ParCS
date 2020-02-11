module explicit_Eul1_mod

use container_abstract_mod,  only : state_abstract_t, model_parameters_abstract_t
use stvec_abstract_mod,      only : stvec_abstract_t
use timescheme_abstract_mod, only : timescheme_abstract_t
use operator_abstract_mod,   only : operator_abstract_t

implicit none

type, extends(timescheme_abstract_t), public :: explicit_Eul1_t

    class(operator_abstract_t), allocatable :: operator

contains

    procedure, public :: step => step_expl_Eul1

end type explicit_Eul1_t

contains

function init_expl_Eul1_ts(operator) result(Eul1_ts)
    type(explicit_Eul1_t) Eul1_ts
    class(operator_abstract_t), intent(in) :: operator

    Eul1_ts%operator = operator
end function init_expl_Eul1_ts

subroutine step_expl_Eul1(this, v0, model_params, dt)

    class(explicit_Eul1_t),             intent(inout) :: this
    class(state_abstract_t),            intent(inout) :: v0
    class(model_parameters_abstract_t), intent(in)    :: model_params
    real(kind=8),                       intent(in)    :: dt

    !local
    class(stvec_abstract_t), allocatable :: v

    select type(v0)
        class is (stvec_abstract_t)
            allocate(v,source=v0) !infer specific type of v0 into v
            call this%operator%act(v,v0)
            call v0.add(v,1._8,dt)
        class default
            call avost("Explicit Eulerian 1st order scheme works only with class(stvec_abstract_t)")
    end select
end subroutine step_expl_Eul1

end module explicit_Eul1_mod
