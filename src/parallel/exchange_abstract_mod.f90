module exchange_abstract_mod

use grid_field_mod, only : grid_field_t

implicit none

type, abstract, public :: exchange_t
contains
    procedure(exchange_procedure), deferred :: do
end type exchange_t

abstract interface
    subroutine exchange_procedure(this, f)
        import exchange_t, grid_field_t
        class(exchange_t),  intent(inout) :: this
        type(grid_field_t), intent(inout) :: f
    end subroutine
end interface

contains

end module exchange_abstract_mod
