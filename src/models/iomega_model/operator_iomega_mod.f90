module operator_iomega_mod

use operator_abstract_mod,  only: operator_abstract_t
use stvec_abstract_mod,     only: stvec_abstract_t
use stvec_iomega_mod,       only: stvec_iomega_t
use container_abstract_mod, only: model_parameters_abstract_t
use parameters_iomega_mod,  only: parameters_iomega_t

implicit none

type, extends(operator_abstract_t) :: operator_iomega_t
    integer(kind=4) N
    complex(kind=8), allocatable :: omega(:) !eigen values, imag == oscillation frequency,
                                             !              real == amplification/decay
    contains
    procedure, public :: act => act
end type operator_iomega_t

contains

subroutine act(this,vout,vin,model_params)
    class(operator_iomega_t),           intent(inout) :: this
    class(stvec_abstract_t),            intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_abstract_t),            intent(in)    :: vin
    class(model_parameters_abstract_t), intent(in)    :: model_params

    integer i

    select type (model_params)
    class is (parameters_iomega_t)

    select type (vout)
    class is (stvec_iomega_t)
        select type (vin)
            class is (stvec_iomega_t)
                do i=1,model_params%N
                    vout%f(i) = model_params%omega(i)*vin%f(i)
                end do
        class default
            print *, "iomega operator failure: vin of wrong type"
            stop
        end select
    class default
        print *, "iomega operator failure: vout of wrong type"
        stop
    end select

    class default
        print *, "non-iomega parameters are passed to iomega model operator"
        stop
    end select
end subroutine act

end module operator_iomega_mod
