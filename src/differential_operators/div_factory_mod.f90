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
    elseif(div_operator_name == 'divergence_a2_ecs'  .or. &
           div_operator_name == 'divergence_a2_cons' .or. &
           div_operator_name == 'divergence_a2_fv') then
        div = create_div_a2_operator(domain,div_operator_name)
    elseif(div_operator_name == 'divergence_ah2') then
        div = create_div_ah2_operator(domain)
    elseif(div_operator_name == 'divergence_ah42_sbp') then
        div = create_div_ah42_sbp_operator(domain)
    else
        call parcomm_global%abort("unknown divergence operator: "//div_operator_name)
    end if
end

function create_div_a2_operator(domain, div_operator_name) result(div)

    use div_a2_mod,       only : div_a2_t
    use halo_factory_mod, only : create_vector_halo_procedure

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: div_operator_name
    type(div_a2_t) :: div

    integer(kind=4), parameter :: ecs_halo_width=2, default_halo_width=1

    if(div_operator_name == "divergence_a2_ecs") then
        call create_vector_halo_procedure(div%halo_procedure,domain,ecs_halo_width,"ecs_A_vec")
        div%subtype="default"
    else if(div_operator_name == "divergence_a2_cons") then
        call create_vector_halo_procedure(div%halo_procedure,domain,default_halo_width,"A_vec_default")
        div%subtype="cons"
    else if(div_operator_name == "divergence_a2_fv") then
        call create_vector_halo_procedure(div%halo_procedure,domain,default_halo_width,"A_vec_default")
        div%subtype="fv"
    else
        call parcomm_global%abort("div_factory_mod, create_div_a2_operator, "// &
                                  "unknown divergence_a2 operator subtype: "// div_operator_name)
    end if

end function create_div_a2_operator

function create_div_ah2_operator(domain) result(div)

    use div_ah2_mod,          only : div_ah2_t
    use exchange_factory_mod, only : create_symm_halo_exchange_Ah

    type(domain_t),   intent(in)  :: domain
    type(div_ah2_t)               :: div

    integer(kind=4), parameter :: halo_width=2

    div%exch_halo = create_symm_halo_exchange_Ah(domain%partition, domain%parcomm, &
                                                 domain%topology,  halo_width, 'full')

end function create_div_ah2_operator

function create_div_ah42_sbp_operator(domain) result(div)

    use div_ah42_sbp_mod,     only : div_ah42_sbp_t
    use exchange_factory_mod, only : create_symm_halo_exchange_Ah

    type(domain_t),   intent(in)  :: domain
    type(div_ah42_sbp_t)          :: div

    integer(kind=4), parameter :: halo_width_interior=3
    integer(kind=4), parameter :: halo_width_edges=1

    div%exch_uv_interior =  &
                    create_symm_halo_exchange_Ah(domain%partition, domain%parcomm, &
                                                 domain%topology,  halo_width_interior, 'full')
    div%exch_div_edges =  &
                    create_symm_halo_exchange_Ah(domain%partition, domain%parcomm, &
                                                 domain%topology,  halo_width_edges, 'full')

end function create_div_ah42_sbp_operator

end module div_factory_mod
