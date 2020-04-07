module ars232_mod

use container_abstract_mod,  only : state_abstract_t, model_parameters_abstract_t
use stvec_abstract_mod,      only : stvec_abstract_t
use timescheme_abstract_mod, only : timescheme_abstract_t
use operator_abstract_mod,   only : operator_abstract_t

implicit none

real(kind=8), parameter :: ae(4,4) = [ [0.0         , 0.0         , 0.0         , 0.0], &
                                       [0.4358665215, 0.0         , 0.0         , 0.0], &
                                       [0.3212788860, 0.3966543747, 0.0         , 0.0], &
                                       [-0.105858296, 0.5529291479, 0.5529291479, 0.0]]
real(kind=8), parameter :: b(4)    =   [0.0,          1.208496649,  -0.644363171, 0.4358665215]
real(kind=8), parameter :: ai(4,4) = [ [0.0,         0.0,          0.0,          0.0],&
                                       [0.0,         0.4358665215, 0.0,          0.0],&
                                       [0.0,         0.2820667392, 0.4358665215, 0.0],&
                                       [0.0,         1.208496649,  -0.644363171, 0.4358665215]]

type, extends(timescheme_abstract_t), public :: ars232_t

    class(operator_abstract_t), allocatable :: oper_e, oper_i
    class(stvec_abstract_t),    allocatable :: y1,y2,y3,y4
    class(stvec_abstract_t),    allocatable :: q2,q3,q4
    class(stvec_abstract_t),    allocatable :: r,s

contains

    procedure, public :: step => step_ars232

end type ars232_t

contains

subroutine init_ars232(ts_ars232, oper_e, oper_i, v)
    class(timescheme_abstract_t),   allocatable, &
                                    intent(inout) :: ts_ars232
    class(operator_abstract_t),     intent(in)    :: oper_e, oper_i
    class(stvec_abstract_t),        intent(in)    :: v !example of model state vector

    allocate(ars232_t :: ts_ars232)

    select type(ts_ars232)
        class is (ars232_t)

    ts_ars232%oper_e = oper_e
    ts_ars232%oper_i = oper_i

    !preallocate additional state vectors
    allocate(ts_ars232%y1, source=v)
    allocate(ts_ars232%y2, source=v)
    allocate(ts_ars232%y3, source=v)
    allocate(ts_ars232%y4, source=v)
    allocate(ts_ars232%q2, source=v)
    allocate(ts_ars232%q3, source=v)
    allocate(ts_ars232%q4, source=v)
    allocate(ts_ars232%r,  source=v)
    allocate(ts_ars232%s,  source=v)

    class default
        call avost("type mismatch in init_ars232")
    end select

end subroutine init_ars232

subroutine step_ars232(this, v0, model_params, dt)

    class(ars232_t),                    intent(inout) :: this
    class(state_abstract_t),            intent(inout) :: v0
    class(model_parameters_abstract_t), intent(in)    :: model_params
    real(kind=8),                       intent(in)    :: dt

    select type(v0)
        class is (stvec_abstract_t)

!    !solv test
!    call this%oper_e%act(this%y1,v0,model_params)
!    call this%y1%add(v0,3.0_8*dt,1._8)
!    call this%r%copy(this%y1)
!
!    call this%oper_i%act(this%y2,this%y1,model_params)
!    call this%r%add(this%y2,1._8,-2.5_8*dt)
!
!    call this%oper_i%solv(2.5_8*dt,this%y3,this%r,model_params)

    !call test_solution(this%y3, this%y1, this%oper_i, 0.5_8*dt)
    !call setw0(this%y1)
    !call setw0(this%y3)

    call this%y3%add(this%y1,1._8,-1._8)
    print *, this%y3%norm() / this%y1%norm()

    stop

!    call this%oper_e%act(this%y1,v0,model_params)
!
!    call this%r%copy(v0); call this%r%add(this%y1,1._8,ae(1,2))
!    call this%s%copy(this%r)
!    call this%oper_e%act(this%y2,this%s,model_params)
!
!    call this%r%copy(v0); call this%r%add(this%y1,1._8,ae(1,3)); call this%r%add(this%y2,1._8,ae(2,3))
!    call this%s%copy(this%r)
!    call this%oper_e%act(this%y3,this%s,model_params)
!
!    call this%r%copy(v0); call this%r%add(this%y1,1._8,ae(1,4)); call this%r%add(this%y2,1._8,ae(2,4))
!                          call this%r%add(this%y3,1._8,ae(3,4))
!    call this%s%copy(this%r)
!    call this%oper_e%act(this%y4,this%s,model_params)
!
!    call v0%add(this%y2,1._8,b(2)*dt)
!    call v0%add(this%y3,1._8,b(3)*dt)
!    call v0%add(this%y4,1._8,b(4)*dt)
!
        class default
            call avost("ARS232 scheme works only with class(stvec_abstract_t)")
    end select

end subroutine step_ars232

!subroutine test_solution(y, r, oper, dt)
!    use stvec_abstract_mod, only: stvec_abstract_t
!    use stvec_NHlin_mod, only :stvec_NHlin_t
!    use operator_NHlin_split_mod, only : operator_NHlin_implicit_t
!
!    class(stvec_abstract_t) , intent(in) :: y, r
!    class(operator_abstract_t), intent(in) :: oper
!    real(kind=8)             ,  intent(in) :: dt
!
!    integer k,i,j, ind
!
!    select type(oper)
!    class is (operator_NHlin_implicit_t)
!    select type(y)
!    class is(stvec_NHlin_t)
!        select type(r)
!        class is(stvec_NHlin_t)
!            do k=1,9
!                print *, "rhs", oper%d2z(-1,k)*r%w(1)%p(1,1,k-1)+&
!                               (oper%d2z( 0,k)+oper%cwz(k)/dt**2)*r%w(1)%p(1,1,k)+&
!                               oper%d2z( 1,k)*r%w(1)%p(1,1,k+1)
!            end do
!
!            i=1; j=1; ind=1
!            do k=0,10
!                print *, "W", y%w(ind)%p(i,j,k) -r%w(ind)%p(i,j,k)
!            end do
!        end select
!    end select
!end select
!end subroutine test_solution
!
!subroutine setw0(y)
!    use stvec_abstract_mod
!    use stvec_NHlin_mod
!    class(stvec_abstract_t) y
!        select type(y)
!    class is (stvec_NHlin_t)
!        y%prex(1)%p=0._8
!        y%prex(2)%p=0._8
!        y%prex(3)%p=0._8
!        y%prex(4)%p=0._8
!        y%prex(5)%p=0._8
!        y%prex(6)%p=0._8
!        print *, "set!"
!    end select
!end subroutine setw0

end module ars232_mod
