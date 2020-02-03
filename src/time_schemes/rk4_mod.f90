module rk4_mod

use stvec_abstract_mod,      only : stvec_abstract_t
use timescheme_abstract_mod, only : timescheme_abstract_t

implicit none

type, extends(timescheme_abstract_t), public :: rk4_t

    class(stvec_abstract_t), allocatable :: k1, k2, k3, k4, y

contains

    procedure, public :: step => step_rk4

end type rk4_t

contains

subroutine init_rk4(ts_rk4, v)
    !preallocate additional arrays
    type(rk4_t),              intent(inout) :: ts_rk4
    class(stvec_abstract_t),  intent(in)    :: v !example of model state vector

    allocate(ts_rk4%k1, source=v)
    allocate(ts_rk4%k2, source=v)
    allocate(ts_rk4%k3, source=v)
    allocate(ts_rk4%k4, source=v)
    allocate(ts_rk4%y,  source=v)

end subroutine init_rk4

subroutine step_rk4(this, operator, v0, dt)

    use operator_abstract_mod, only : operator_abstract_t

    class(rk4_t),                    intent(inout) :: this
    class(operator_abstract_t),      intent(inout) :: operator
    class(stvec_abstract_t), target, intent(inout) :: v0
    real(kind=8),                    intent(in)    :: dt

    !local
    class(stvec_abstract_t), allocatable :: v

    call operator%act(this%k1,v0)

    call this%y%copy(v0)
    call this%y%add(this%k1,1._8,0.5_8*dt)
    call operator%act(this%k2,this%y)

    call this%y%copy(v0)
    call this%y%add(this%k2,1._8,0.5_8*dt)
    call operator%act(this%k3,this%y)

    call this%y%copy(v0)
    call this%y%add(this%k3,1.0_8,dt)
    call operator%act(this%k4,this%y)

    call v0%add(this%k1,1._8,dt/6._8)
    call v0%add(this%k2,1._8,dt/3._8)
    call v0%add(this%k3,1._8,dt/3._8)
    call v0%add(this%k4,1._8,dt/6._8)

end subroutine step_rk4

end module rk4_mod
