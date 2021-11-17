module sbp_vertical_gradient_mod

use abstract_vertical_operator_mod, only : vertical_operator_t
use domain_mod,                     only : domain_t
use grid_field_mod,                 only : grid_field_t
use sbp_operator_mod,               only : sbp_operator_t

implicit none

type, extends(vertical_operator_t) :: sbp_vertical_grad_t
    class(sbp_operator_t), allocatable :: sbp_op
    contains
    procedure :: apply
end type sbp_vertical_grad_t

contains

subroutine apply(this, f_out, f_in, domain)
    class(sbp_vertical_grad_t), intent(in)    :: this
    type(grid_field_t),         intent(inout) :: f_out
    type(grid_field_t),         intent(in)    :: f_in
    type(domain_t),             intent(in)    :: domain

    call f_out%assign(0.0_8,domain%mesh_n)
end subroutine apply

end module sbp_vertical_gradient_mod
