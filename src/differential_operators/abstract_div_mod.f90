module abstract_div_mod

use grid_function_mod, only : grid_function_t
use domain_mod,        only : domain_t

implicit none

type, abstract, public :: div_operator_t

contains

procedure(div_calc_procedure), deferred :: calc_div

end type div_operator_t

abstract interface
    subroutine div_calc_procedure(this, div, u, v, domain, multiplier)
        import div_operator_t, grid_function_t, domain_t
        class(div_operator_t),  intent(inout) :: this
        type(domain_t),         intent(in)    :: domain
        type(grid_function_t),  intent(inout) :: u, v
        real(kind=8), optional, intent(in)    :: multiplier
        !out put
        type(grid_function_t),  intent(inout) :: div
    end subroutine div_calc_procedure
end interface

    contains

end module abstract_div_mod
