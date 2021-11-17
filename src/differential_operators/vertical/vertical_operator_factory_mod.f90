module vertical_operator_factory_mod

use abstract_vertical_operator_mod, only : vertical_operator_t
use parcomm_mod,                    only : parcomm_global

implicit none

contains

subroutine create_vertical_operator(vertical_op, vertical_op_name)
    class(vertical_operator_t), allocatable, intent(out) :: vertical_op
    character(len=*), intent(in) :: vertical_op_name

    select case(vertical_op_name)
    case("vertical_grad_staggered_sbp21")
        call create_sbp_vertical_op(vertical_op, "D21_staggered_c2i")
    case("vertical_grad_staggered_sbp42")
        call create_sbp_vertical_op(vertical_op, "D42_staggered_c2i")
    case("vertical_grad_sbp21")
        call create_sbp_vertical_op(vertical_op, "d21")
    case("vertical_grad_sbp42")
        call create_sbp_vertical_op(vertical_op, "d42")
    case("vertical_grad_sbp63")
        call create_sbp_vertical_op(vertical_op, "d63")
    case default
        call parcomm_global%abort("create_vertical_operator error, unknown operator "//&
                                  vertical_op_name)
    end select
end subroutine create_vertical_operator

subroutine create_sbp_vertical_op(vertical_op, sbp_op_name)

    use sbp_vertical_operator_mod, only : sbp_vertical_op_t
    use sbp_factory_mod,           only : create_sbp_operator

    class(vertical_operator_t), allocatable, intent(out) :: vertical_op
    character(len=*), intent(in) :: sbp_op_name

    type(sbp_vertical_op_t), allocatable :: sbp_vertical_op

    allocate(sbp_vertical_op)
    sbp_vertical_op%sbp_op = create_sbp_operator(sbp_op_name)

    call move_alloc(sbp_vertical_op, vertical_op)
end subroutine create_sbp_vertical_op

end module vertical_operator_factory_mod
