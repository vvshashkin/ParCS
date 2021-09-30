module KE_factory_mod

use abstract_KE_mod, only : KE_operator_t
use domain_mod,      only : domain_t
use parcomm_mod,     only : parcomm_global

implicit none

contains

subroutine create_KE_operator(ke_operator, ke_operator_name, domain)

    use ke_unstaggered_mod, only : ke_unstaggered_t

    class(KE_operator_t), allocatable, intent(out) :: ke_operator
    character(len=*),                  intent(in)  :: ke_operator_name
    type(domain_t),                    intent(in)  :: domain

    select case(ke_operator_name)

    case("KE_A_Ah")
        ke_operator = ke_unstaggered_t()
    case default
        call parcomm_global%abort("Unknown KE operator: "//ke_operator_name)
    end select

end subroutine create_KE_operator

end module KE_factory_mod
