module operator_iomega_mod

use operator_abstract_mod, only: operator_abstract_t
use stvec_abstract_mod,    only: stvec_abstract_t
use stvec_iomega_mod,      only: stvec_iomega_t

implicit none

type, extends(operator_abstract_t) :: operator_iomega_t
    integer(kind=4) N
    complex(kind=8), allocatable :: omega(:) !eigen values, imag == oscillation frequency,
                                             !              real == amplification/decay
    contains
    procedure, public :: act => act
end type operator_iomega_t

contains

subroutine init_operator_iomega(new_operator, N, omega)
    type(operator_iomega_t), intent(out) :: new_operator
    integer(kind=4),         intent(in)  :: N
    complex(kind=8),         intent(in)  :: omega(1:N)

    if(allocated(new_operator%omega) .and. new_operator%N /= N) then
        deallocate(new_operator%omega)
    end if
    if (.not. allocated(new_operator%omega)) then
        allocate(new_operator%omega(1:N))
    end if

    new_operator%N = N
    new_operator%omega(1:N) = omega(1:N)

end subroutine init_operator_iomega

subroutine act(this,vout,vin)
    class(operator_iomega_t), intent(in)    :: this
    class(stvec_abstract_t),  intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_abstract_t),  intent(in)    :: vin

    integer i

    select type (vout)
    class is (stvec_iomega_t)
        select type (vin)
            class is (stvec_iomega_t)
                !vout%N = this%N
                !allocate(vout%f(1:this%N)) !thanks to intent(in)!
                do i=1,this%N
                    vout%f(i) = this%omega(i)*vin%f(i)
                end do
        class default
            stop
        end select
    class default
        stop
    end select
end subroutine act

end module operator_iomega_mod
