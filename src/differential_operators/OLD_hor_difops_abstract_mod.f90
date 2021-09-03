!Abstract interfaces for horizontal differential (and also non-differential) operators
module hor_difops_abstract_mod

use grid_field_mod, only : block_t
use mesh_mod,       only : mesh_t

implicit none

abstract interface
    subroutine gradient(gx, gy, f, mesh, multiplier)
        import block_t, mesh_t
        type(block_t),          intent(inout) :: gx, gy
        type(block_t),          intent(in)    :: f
        type(mesh_t),           intent(in)    :: mesh
        real(kind=8), optional, intent(in)    :: multiplier
    end subroutine gradient
    subroutine divergence(div, u, v, mesh, multiplier)
        import block_t, mesh_t
        type(block_t),          intent(inout) :: div
        type(block_t),          intent(in)    :: u, v
        type(mesh_t),           intent(in)    :: mesh
        real(kind=8), optional, intent(in)    :: multiplier
    end subroutine divergence
end interface

end module hor_difops_abstract_mod
