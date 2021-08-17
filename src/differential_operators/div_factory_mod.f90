module div_factory_mod

use domain_mod,        only : domain_t
use mesh_mod,          only : mesh_t
use abstract_div_mod,  only : div_operator_t
use parcomm_mod,       only : parcomm_global

implicit none

contains

function create_div_operator(domain, div_operator_name) result(div)
    type(domain_t),    intent(in)    :: domain
    character(len=*),  intent(in)    :: div_operator_name

    class(div_operator_t), allocatable :: div

    if(div_operator_name == 'divergence_c2') then
        call parcomm_global%abort("not implemented: "//div_operator_name)
        !div = div_c2_t()
    elseif(div_operator_name == 'divergence_a2_ecs') then
        div = create_div_a2_ecs_operator(domain)
    else
        call parcomm_global%abort("unknown divergence operator: "//div_operator_name)
    end if
end

function create_div_a2_ecs_operator(domain) result(div)

    use div_a2_mod,       only : div_a2_t
    use halo_factory_mod, only : create_vector_halo_procedure

    type(domain_t), intent(in)  :: domain
    type(div_a2_t) :: div

    integer(kind=4), parameter :: halo_width=2

    div = div_a2_t()
    call create_vector_halo_procedure(div%halo_procedure,domain,halo_width,"ecs_A_vec")
    
end function create_div_a2_ecs_operator

end module div_factory_mod
