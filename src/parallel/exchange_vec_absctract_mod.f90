module exchange_vec_abstract_mod

use grid_field_mod, only : grid_field_t

implicit none

type, abstract, public :: exchange_vec_t
contains
    procedure(exchange_vec_procedure), deferred :: do
end type exchange_vec_t

abstract interface
    subroutine exchange_vec_procedure(this, u, v)
        import exchange_vec_t, grid_field_t
        class(exchange_vec_t), intent(inout) :: this
        type(grid_field_t),    intent(inout) :: u, v
    end subroutine
end interface

contains


end module exchange_vec_abstract_mod
