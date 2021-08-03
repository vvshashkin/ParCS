module rk4_mod

use stvec_mod,      only : stvec_t
use timescheme_mod, only : timescheme_t
use operator_mod,   only : operator_t
use domain_mod,     only : domain_t

implicit none

private
public :: init_rk4

type, extends(timescheme_t), public :: rk4_t

    class(stvec_t), allocatable :: k1, k2, k3, k4, y

contains

    procedure, public :: step => step_rk4

end type rk4_t

contains

subroutine init_rk4(tscheme, v)
    class(timescheme_t), allocatable, intent(out) :: tscheme
    class(stvec_t), intent(in) :: v !example of model state vector

    type(rk4_t), allocatable :: rk4

    allocate(rk4)
    !preallocate additional state vectors
    call v%create_similar(rk4%k1)
    call v%create_similar(rk4%k2)
    call v%create_similar(rk4%k3)
    call v%create_similar(rk4%k4)
    call v%create_similar(rk4%y)

    call move_alloc(rk4, tscheme)

end subroutine init_rk4

subroutine step_rk4(this, v0, operator, domain, dt)

    class(rk4_t),       intent(inout) :: this
    class(stvec_t),     intent(inout) :: v0
    class(operator_t),  intent(inout) :: operator
    class(domain_t),    intent(in)    :: domain
    real(kind=8),       intent(in)    :: dt

    call operator%apply_to(this%k1,v0,domain)

    call this%y%assign_s1v1s2v2(1.0_8,v0,0.5_8*dt,this%k1,domain)
    call operator%apply_to(this%k2,this%y,domain)

    call this%y%assign_s1v1s2v2(1.0_8,v0,0.5_8*dt,this%k2,domain)
    call operator%apply_to(this%k3,this%y,domain)

    call this%y%assign_s1v1s2v2(1.0_8,v0,dt,this%k3,domain)
    call operator%apply_to(this%k4,this%y,domain)

    call v0%update_s1v1s2v2(dt/6.0_8,this%k1,dt/3.0_8,this%k2,domain)
    call v0%update_s1v1s2v2(dt/3.0_8,this%k3,dt/6.0_8,this%k4,domain)

end subroutine step_rk4

end module rk4_mod
