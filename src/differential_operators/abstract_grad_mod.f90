module abstract_grad_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t
    
implicit none
    
type, abstract, public :: grad_operator_t
    
contains
    
procedure(grad_calc_procedure), deferred :: calc_grad
    
end type grad_operator_t
    
abstract interface
    subroutine grad_calc_procedure(this, gx, gy, f, domain, multiplier)
        import grad_operator_t, grid_field_t, domain_t
        class(grad_operator_t), intent(inout) :: this
        type(domain_t),         intent(in)    :: domain
        type(grid_field_t),     intent(inout) :: f
        real(kind=8), optional, intent(in)    :: multiplier
        !output:
        type(grid_field_t),     intent(inout) :: gx, gy
    end subroutine grad_calc_procedure
end interface
    
    contains
   
end module abstract_grad_mod
