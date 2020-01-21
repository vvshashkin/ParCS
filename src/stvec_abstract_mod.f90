module stvec_abstract_mod

implicit none

type, abstract, public :: stvec_abstract_t

contains

    procedure(add),  deferred :: add
    procedure(copy), deferred :: copy
    procedure(dot),  deferred :: dot

end type stvec_abstract_t

abstract interface
    subroutine add(this, other, alpha, beta)
        !calculates linear combination alpha*this+beta*other
        import stvec_abstract_t
        class(stvec_abstract_t), intent(inout) :: this
        class(stvec_abstract_t), intent(in)    :: other
        real(kind=8),  intent(in)    :: alpha, beta
    end subroutine add

    subroutine copy(this, source_stvec)
        !copies information from source_stvec to this
        import stvec_abstract_t
        class(stvec_abstract_t), intent(inout) :: this
        class(stvec_abstract_t), intent(in)    :: source_stvec
    end subroutine copy

    real(kind=8) function dot(this,other) result(dot_prod)
        !calculates dot(this,other)
        import stvec_abstract_t
        class(stvec_abstract_t), intent(in) :: this
        class(stvec_abstract_t), intent(in)    :: other
    end function dot

end interface

contains

end module stvec_abstract_mod
