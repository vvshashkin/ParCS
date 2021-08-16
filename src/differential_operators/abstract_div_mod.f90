module abstract_div_mod

use grid_function_mod, only : grid_function_t
use partition_mod,     only : partition_t
use mesh_mod,          only : mesh_t

implicit none

type, abstract, public :: div_operator_t

contains

procedure(div_calc_procedure), deferred :: calc_div

end type div_operator_t

abstract interface
    subroutine div_calc_procedure(this, u, v, div, mesh, partition, multiplier)
        import div_operator_t, grid_function_t, mesh_t, partition_t
        class(div_operator_t),  intent(inout) :: this
        type(partition_t),      intent(in)    :: partition
        type(grid_function_t),  intent(inout) ::   u(partition%ts:partition%te), &
                                                   v(partition%ts:partition%te), &
                                                 div(partition%ts:partition%te)
        type(mesh_t),           intent(in)    :: mesh(partition%ts:partition%te)
        real(kind=8), optional, intent(in)    :: multiplier
    end subroutine div_calc_procedure
end interface

    contains

end module abstract_div_mod
