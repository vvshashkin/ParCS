module rk4_mod

use state_abstract_mod,      only : state_abstract_t
use stvec_abstract_mod,      only : stvec_abstract_t
use timescheme_abstract_mod, only : timescheme_abstract_t
use operator_abstract_mod,   only : operator_abstract_t

implicit none

type, extends(timescheme_abstract_t), public :: rk4_t

    class(operator_abstract_t), allocatable :: operator
    class(stvec_abstract_t),    allocatable :: k1, k2, k3, k4, y

contains

    procedure, public :: step => step_rk4

end type rk4_t

contains

function init_rk4(operator, v) result(ts_rk4)
    type(rk4_t)                               :: ts_rk4
    class(operator_abstract_t), intent(in)    :: operator
    class(stvec_abstract_t),    intent(in)    :: v !example of model state vector

    ts_rk4%operator = operator

    !preallocate additional state vectors
    allocate(ts_rk4%k1, source=v)
    allocate(ts_rk4%k2, source=v)
    allocate(ts_rk4%k3, source=v)
    allocate(ts_rk4%k4, source=v)
    allocate(ts_rk4%y,  source=v)

end function init_rk4

subroutine step_rk4(this, v0, dt)

    class(rk4_t),                    intent(inout) :: this
    class(state_abstract_t), target, intent(inout) :: v0
    real(kind=8),                    intent(in)    :: dt

    select type(v0)
        class is (stvec_abstract_t)

    call this%operator%act(this%k1,v0)

    call this%y%copy(v0)
    call this%y%add(this%k1,1._8,0.5_8*dt)
    call this%operator%act(this%k2,this%y)

    call this%y%copy(v0)
    call this%y%add(this%k2,1._8,0.5_8*dt)
    call this%operator%act(this%k3,this%y)

    call this%y%copy(v0)
    call this%y%add(this%k3,1.0_8,dt)
    call this%operator%act(this%k4,this%y)

    call v0%add(this%k1,1._8,dt/6._8)
    call v0%add(this%k2,1._8,dt/3._8)
    call v0%add(this%k3,1._8,dt/3._8)
    call v0%add(this%k4,1._8,dt/6._8)

        class default
            call avost("RK4 scheme works only with class(stvec_abstract_t)")
    end select


end subroutine step_rk4

end module rk4_mod
