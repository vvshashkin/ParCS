module exchange_abstract_mod

use grid_function_mod, only : grid_function_t

implicit none

type, abstract, public :: exchange_t

contains

    procedure(exchange_procedure), deferred :: do

end type exchange_t

abstract interface
    subroutine exchange_procedure(this, f, ts, te)
        import exchange_t, grid_function_t
        class(exchange_t),     intent(inout) :: this
        integer(kind=4),       intent(in)    :: ts, te
        type(grid_function_t), intent(inout) :: f(ts:te)
    end subroutine
end interface

contains


end module exchange_abstract_mod
