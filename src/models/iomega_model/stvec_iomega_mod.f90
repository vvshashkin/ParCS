module stvec_iomega_mod

use stvec_abstract_mod, only: stvec_abstract_t

implicit none

type, extends(stvec_abstract_t) :: stvec_iomega_t
    integer(kind=4) N
    complex(kind=8), allocatable :: f(:)

    contains
    procedure, public :: add  => add
    procedure, public :: copy => copy
    procedure, public :: dot  => dot
    procedure, public :: norm
    procedure, public :: smult
end type stvec_iomega_t

contains

subroutine init_stvec_iomega(new_stvec, N, f)
    type(stvec_iomega_t),      intent(out) :: new_stvec
    integer(kind=4),           intent(in)  :: N
    complex(kind=8), optional, intent(in)  :: f(1:N)

    if(allocated(new_stvec%f) .and. new_stvec%N /= N) then
        deallocate(new_stvec%f)
    end if
    if (.not. allocated(new_stvec%f)) then
        allocate(new_stvec%f(1:N))
    end if

    new_stvec%N = N
    if(present(f)) then
        new_stvec%f(1:N) = f(1:N)
    else
        new_stvec%f(1:N) = 0._8
    end if

end subroutine init_stvec_iomega

subroutine add(this,other,alpha,beta)
    class(stvec_iomega_t),   intent(inout) :: this
    class(stvec_abstract_t), intent(in)    :: other
    real(kind=8),         intent(in)    :: alpha, beta

    select type (other)
    class is (stvec_iomega_t)
        this%f(1:this%N) = alpha*this%f(1:this%N)+beta*other%f(1:this%N)
    class default
        stop
    end select
end subroutine add

subroutine copy(this,source_stvec)
    class(stvec_iomega_t),   intent(inout) :: this
    class(stvec_abstract_t), intent(in)    :: source_stvec

    select type (source_stvec)
    class is (stvec_iomega_t)
        call init_stvec_iomega(this,source_stvec%N,source_stvec%f)
    class default
        stop
    end select
end subroutine copy

real(kind=8) function dot(this, other) result(dot_prod)
    class(stvec_iomega_t),   intent(in)    :: this
    class(stvec_abstract_t), intent(in)    :: other

    select type (other)
    class is (stvec_iomega_t)
        !N-vector of complex is considered as 2N-vector of real values
        dot_prod = sum(real(this%f(1:this%N))*real(other%f(1:this%N))+ &
                       imag(this%f(1:this%N))*imag(other%f(1:this%N)))
    class default
        stop
    end select

end function dot

real(kind=8) function norm(this) result(l2)
    class(stvec_iomega_t), intent(in) :: this
    !calculates norm of state vector
    l2 = sqrt(sum(real(this%f(1:this%n))**2+imag(this%f(1:this%n))**2))
end function norm

subroutine smult(this, alpha)
    !multiplicates state vector by scalar
    class(stvec_iomega_t), intent(inout) :: this
    real(kind=8),          intent(in)    :: alpha

    this%f(1:this%N) = alpha*this%f(1:this%N)

end subroutine smult

end module stvec_iomega_mod
