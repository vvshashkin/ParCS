!Abstract interfaces for horizontal differential (and also non-differential) operators
module hor_difops_abstract_mod

use grid_function_mod, only : grid_function_t
use mesh_mod,          only : mesh_t

implicit none

abstract interface
    subroutine gradient(gx, gy, f, mesh, multiplier)
        import grid_function_t, mesh_t
        type(grid_function_t),  intent(inout) :: gx, gy
        type(grid_function_t),  intent(in)    :: f
        type(mesh_t),           intent(in)    :: mesh
        real(kind=8), optional, intent(in)    :: multiplier
    end subroutine gradient
    subroutine divergence(div, u, v, mesh, multiplier)
        import grid_function_t, mesh_t
        type(grid_function_t),  intent(inout) :: div
        type(grid_function_t),  intent(in)    :: u, v
        type(mesh_t),           intent(in)    :: mesh
        real(kind=8), optional, intent(in)    :: multiplier
    end subroutine divergence
end interface

end module hor_difops_abstract_mod
