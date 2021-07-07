module exp_taylor_mod

use container_abstract_mod,  only : state_abstract_t, model_parameters_abstract_t
use stvec_abstract_mod,      only : stvec_abstract_t
use timescheme_abstract_mod, only : timescheme_abstract_t
use operator_abstract_mod,   only : operator_abstract_t

implicit none

type, extends(timescheme_abstract_t), public :: exp_taylor_t

    class(operator_abstract_t), allocatable :: operator
    class(stvec_abstract_t),    allocatable :: v1, v2
    real(kind=8)                            :: tolerance = 1e-3_8 !norm(residual)/norm(f0)
    integer(kind=4)                         :: maxiter = 1000

contains

    procedure, public :: step => step_exp

end type exp_taylor_t

contains

function init_exp_taylor(operator, v, tolerance, maxiter) result(ts_exp)
    type(exp_taylor_t)                        :: ts_exp
    class(operator_abstract_t), intent(in)    :: operator
    class(stvec_abstract_t),    intent(in)    :: v !example of model state vector
    real(kind=8),     optional, intent(in)    :: tolerance
    integer(kind=4),  optional, intent(in)    :: maxiter


    if(present(tolerance)) ts_exp%tolerance = tolerance
    if(present(maxiter)) ts_exp%maxiter = maxiter

    ts_exp%operator = operator

    !preallocate additional state vectors
    allocate(ts_exp%v1, mold=v)
    allocate(ts_exp%v2, mold=v)

end function init_exp_taylor

subroutine step_exp(this, v0, model_params, dt)

    class(exp_taylor_t),                intent(inout) :: this
    class(state_abstract_t),            intent(inout) :: v0
    class(model_parameters_abstract_t), intent(in)    :: model_params
    real(kind=8),                       intent(in)    :: dt

    real(kind=8)     :: norm0, norm_res
    integer(kind=4)  :: iter

    select type(v0)
        class is (stvec_abstract_t)

    call this%v1%copy(v0)
    call this%v2%copy(v0)

    iter = 1
    norm0 = v0%norm()
    norm_res = norm0
    do while(norm_res/norm0 > this%tolerance .and. iter <= this%maxiter)
        call this%operator%act(this%v2,this%v1,model_params)
        call this%v2%smult(dt/dble(iter))
        call v0%add(this%v2,1._8,1._8)
        call this%v1%copy(this%v2)
        norm_res = this%v1%norm()
        !print *, iter, norm_res/norm0
        iter = iter+1
    end do

    class default
        call avost("EXP_Taylor scheme works only with class(stvec_abstract_t)")
    end select


end subroutine step_exp

end module exp_taylor_mod
