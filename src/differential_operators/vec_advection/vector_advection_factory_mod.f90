module vector_advection_factory_mod

use domain_mod,                     only : domain_t
use abstract_vector_advection_mod,  only : vector_advection_operator_t
use parcomm_mod,                    only : parcomm_global
use grid_field_mod,                 only : grid_field_t

implicit none

private

public :: create_vector_advection_operator

contains

subroutine create_vector_advection_operator(vec_advection_op, vec_advection_op_name, domain)

    class(vector_advection_operator_t), allocatable, intent(out) :: vec_advection_op
    character(len=*),                                intent(in)  :: vec_advection_op_name
    type(domain_t),                                  intent(in)  :: domain

    select case(vec_advection_op_name)

    case("vector_advection_Ah21_covariant")
        call create_vector_advection_Ah_covariant(vec_advection_op, "d21", 1, domain)
    case("vector_advection_Ah42_covariant")
        call create_vector_advection_Ah_covariant(vec_advection_op, "d42", 3, domain)
    case("vector_advection_Ah63_covariant")
        call create_vector_advection_Ah_covariant(vec_advection_op, "d63", 5, domain)
    case default
        call parcomm_global%abort("Unknown vector advection operator: "//vec_advection_op_name)
    end select

end subroutine create_vector_advection_operator

subroutine create_vector_advection_Ah_covariant(vec_advection_op,sbp_operator_name, &
                                                halo_width, domain)

    use vector_advection_Ah_mod,   only : vector_advection_Ah_t
    use sbp_factory_mod,           only : create_sbp_operator
    use exchange_factory_mod,      only : create_symm_halo_exchange_Ah
    use halo_factory_mod,          only : create_vector_halo_procedure

    class(vector_advection_operator_t), allocatable, intent(out) :: vec_advection_op
    character(len=*),                                intent(in)  :: sbp_operator_name
    integer(kind=4),                                 intent(in)  :: halo_width
    type(domain_t),                                  intent(in)  :: domain

    type(vector_advection_Ah_t), allocatable :: vec_advection_Ah_op

    allocate(vec_advection_Ah_op)

    vec_advection_Ah_op%sbp_op = create_sbp_operator(sbp_operator_name)
    call create_vector_halo_procedure(vec_advection_Ah_op%sync_edges,domain,0, &
                                                    "ecs_Ah_vec_sync_covariant")
    vec_advection_Ah_op%exch_uv_interior =  &
                    create_symm_halo_exchange_Ah(domain%partition, domain%parcomm, &
                                                 domain%topology,  halo_width, 'full')

    call move_alloc(vec_advection_Ah_op, vec_advection_op)

end subroutine create_vector_advection_Ah_covariant

end module vector_advection_factory_mod
