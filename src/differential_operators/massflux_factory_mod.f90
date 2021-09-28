module massflux_factory_mod

use domain_mod,             only : domain_t
use abstract_massflux_mod,  only : massflux_operator_t
use massflux_colocated_mod, only : massflux_colocated_t
use massflux_Cgrid_mod,     only : massflux_chalo_t
use parcomm_mod,            only : parcomm_global
use halo_factory_mod,       only : create_halo_procedure


implicit none

contains

function create_massflux_operator(domain, massflux_operator_name) result(massflux)
    type(domain_t),    intent(in)    :: domain
    character(len=*),  intent(in)    :: massflux_operator_name

    class(massflux_operator_t), allocatable :: massflux

    if(massflux_operator_name == 'massflux_colocated') then
        massflux = massflux_colocated_t()
    elseif(massflux_operator_name == 'massflux_c2') then
        massflux = massflux_chalo_t(order=2)
        select type(massflux)
        type is (massflux_chalo_t)
            call create_halo_procedure(massflux%halo,domain,2,"ECS_O")
        end select
    elseif(massflux_operator_name == 'massflux_c4') then
        massflux = massflux_chalo_t(order=4)
        select type(massflux)
        type is (massflux_chalo_t)
            call create_halo_procedure(massflux%halo,domain,2,"ECS_O")
        end select
    else
        call parcomm_global%abort("unknown massflux operator: "//massflux_operator_name)
    end if

end function create_massflux_operator

end module massflux_factory_mod
