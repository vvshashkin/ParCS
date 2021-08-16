module abstract_grad_mod

use grid_function_mod, only : grid_function_t
use partition_mod,     only : partition_t
use mesh_mod,          only : mesh_t

implicit none

type, abstract, public :: grad_operator_t
contains
    procedure(grad_calc_procedure), deferred :: calc_grad
end type grad_operator_t

abstract interface
    subroutine grad_calc_procedure(this, f, gx, gy, mesh, partition, multiplier)
        import grad_operator_t, grid_function_t, mesh_t, partition_t
        class(grad_operator_t), intent(inout) :: this
        type(partition_t),      intent(in)    :: partition
        type(grid_function_t),  intent(inout) ::  f(partition%ts:partition%te), &
                                                 gx(partition%ts:partition%te), &
                                                 gy(partition%ts:partition%te)
        type(mesh_t),           intent(in)    :: mesh(partition%ts:partition%te)
        real(kind=8), optional, intent(in)    :: multiplier
    end subroutine grad_calc_procedure
end interface

    contains

end module abstract_grad_mod
