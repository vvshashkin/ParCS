module exchange_vec_abstract_mod

use grid_function_mod, only : grid_function_t

implicit none

type, abstract, public :: exchange_vec_t

contains

    procedure(exchange_vec_procedure), deferred :: do

end type exchange_vec_t

abstract interface
    subroutine exchange_vec_procedure(this, u, v, ts, te)
        import exchange_vec_t, grid_function_t
        class(exchange_vec_t),     intent(inout) :: this
        integer(kind=4),       intent(in)    :: ts, te
        type(grid_function_t), intent(inout) :: u(ts:te), v(ts:te)
    end subroutine
end interface

contains


end module exchange_vec_abstract_mod
