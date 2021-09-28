module massflux_factory_mod

use domain_mod,             only : domain_t
use abstract_massflux_mod,  only : massflux_operator_t
use massflux_colocated_mod, only : massflux_colocated_t
use parcomm_mod,            only : parcomm_global

implicit none

contains

function create_massflux_operator(domain, massflux_operator_name) result(massflux)
    type(domain_t),    intent(in)    :: domain
    character(len=*),  intent(in)    :: massflux_operator_name

    class(massflux_operator_t), allocatable :: massflux

    if(massflux_operator_name == 'massflux_colocated') then
        massflux = massflux_colocated_t()
    else
        call parcomm_global%abort("unknown massflux operator: "//massflux_operator_name)
    end if

end function create_massflux_operator

end module massflux_factory_mod
