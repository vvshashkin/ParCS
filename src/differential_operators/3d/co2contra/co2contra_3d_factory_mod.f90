module co2contra_3d_factory_mod

use domain_mod,                 only : domain_t
use abstract_co2contra_3d_mod,  only : co2contra_3d_operator_t
use co2contra_3d_colocated_mod, only : co2contra_3d_colocated_t
use parcomm_mod,                only : parcomm_global

implicit none

contains

subroutine create_co2contra_3d_operator(co2contra_3d_op, domain, co2contra_3d_operator_name)
    class(co2contra_3d_operator_t), allocatable , intent(out) :: co2contra_3d_op
    type(domain_t),   intent(in) :: domain
    character(len=*), intent(in) :: co2contra_3d_operator_name

    select case(co2contra_3d_operator_name)
    case("co2contra_3d_colocated")
        co2contra_3d_op = co2contra_3d_colocated_t()
    case default
        call parcomm_global%abort("unknown co2contra_3d operator: "//co2contra_3d_operator_name)
    end select

end subroutine create_co2contra_3d_operator

end module co2contra_3d_factory_mod
