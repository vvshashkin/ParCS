module vertical_operator_factory_mod

use abstract_vertical_operator_mod, only : vertical_operator_t
use sbp_vertical_gradient_mod,      only : sbp_vertical_grad_t
use parcomm_mod,                    only : parcomm_global

implicit none

contains

subroutine create_vertical_operator(vertical_op, vertical_op_name)
    class(vertical_operator_t), allocatable, intent(out) :: vertical_op
    character(len=*), intent(in) :: vertical_op_name

    vertical_op = sbp_vertical_grad_t()
end subroutine create_vertical_operator

end module vertical_operator_factory_mod
