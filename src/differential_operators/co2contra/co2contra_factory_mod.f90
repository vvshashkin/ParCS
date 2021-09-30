module co2contra_factory_mod

use domain_mod,              only : domain_t
use abstract_co2contra_mod,  only : co2contra_operator_t
use co2contra_colocated_mod, only : co2contra_colocated_t
use co2contra_Cgrid_mod,     only : co2contra_c_sbp_t
use parcomm_mod,             only : parcomm_global
use exchange_factory_mod,    only : create_symmetric_halo_vec_exchange_C

implicit none

contains

function create_co2contra_operator(domain, co2contra_operator_name) result(co2contra)
    type(domain_t),    intent(in)    :: domain
    character(len=*),  intent(in)    :: co2contra_operator_name

    class(co2contra_operator_t), allocatable :: co2contra

    if(co2contra_operator_name == 'co2contra_colocated') then
        co2contra = co2contra_colocated_t()
    else if(co2contra_operator_name == "co2contra_c_sbp21" .or. &
            co2contra_operator_name == "co2contra_c_sbp42") then
        co2contra = create_co2contra_c_sbp21_operator(domain, co2contra_operator_name)
    else
        call parcomm_global%abort("unknown co2contra operator: "//co2contra_operator_name)
    end if

end function create_co2contra_operator

function create_co2contra_c_sbp21_operator(domain, co2contra_operator_name) result(co2contra)

    type(domain_t),          intent(in) :: domain
    character(len=*),        intent(in) :: co2contra_operator_name
    type(co2contra_c_sbp_t)             :: co2contra

    integer(kind=4) :: halo_width

    co2contra%operator_name = co2contra_operator_name

    select case(co2contra_operator_name)
    case("co2contra_c_sbp21")
        halo_width = 1
    case("co2contra_c_sbp42")
        halo_width = 3
    case default
        call parcomm_global%abort("unknown co2contra_c_sbp operator "// co2contra_operator_name)
    end select

    co2contra%exchange_inner = &
          create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                                     domain%topology, halo_width, 'full')
end function

end module co2contra_factory_mod
