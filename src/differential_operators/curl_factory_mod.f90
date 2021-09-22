module curl_factory_mod

use domain_mod,         only : domain_t
use abstract_curl_mod,  only : curl_operator_t
use parcomm_mod,        only : parcomm_global

implicit none

contains

subroutine create_curl_operator_div_based(curl, div_operator_name, domain)

    use curl_div_based_mod,     only : curl_div_based_t
    use div_factory_mod,        only : create_div_operator
    use grid_field_factory_mod, only : create_grid_field

    class(curl_operator_t), allocatable, intent(out) :: curl
    character(len=*),                    intent(in)  :: div_operator_name
    type(domain_t),                      intent(in)  :: domain

    type(curl_div_based_t), allocatable :: curl_div_based

    integer(kind=4) :: halo_width

    !This is temporary solution
    !We need a way to get information about required halo width from operator
    halo_width = 5

    allocate(curl_div_based)

    curl_div_based%div_op = create_div_operator(domain, div_operator_name)

    call create_grid_field(curl_div_based%ut, halo_width, 0, domain%mesh_u)
    call create_grid_field(curl_div_based%vt, halo_width, 0, domain%mesh_v)

    call move_alloc(curl_div_based, curl)

end subroutine create_curl_operator_div_based

end module curl_factory_mod