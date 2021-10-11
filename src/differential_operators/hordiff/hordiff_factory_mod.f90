module hordiff_factory_mod

use domain_mod,             only : domain_t
use grid_field_factory_mod, only : create_grid_field
use abstract_hordiff_mod,   only : hordiff_operator_t
use hordiff_colocated_mod,  only : hordiff_colocated_t

implicit none

contains

subroutine create_hordiff_operator(hordiff_op, hordiff_op_name, hordiff_coeff, domain)

    class(hordiff_operator_t), allocatable, intent(out) :: hordiff_op
    character(len=*),                       intent(in)  :: hordiff_op_name
    real(kind=8),                           intent(in)  :: hordiff_coeff
    type(domain_t),                         intent(in)  :: domain

    select case(hordiff_op_name)
    case("hordiff_c_biharm_div")
        call create_Cgrid_hordiff_div_operator(hordiff_op, hordiff_coeff, domain)
    case("hordiff_colocated")
        hordiff_op = hordiff_colocated_t()
    case default
        call domain%parcomm%abort("Unknown hordiff_op_name")
    end select

end subroutine create_hordiff_operator


subroutine create_Cgrid_hordiff_div_operator(hordiff_op, hordiff_coeff, domain)

    use hordiff_Cgrid_mod,     only : hordiff_c_div_t
    use div_factory_mod,       only : create_div_operator
    use grad_factory_mod,      only : create_grad_operator
    use co2contra_factory_mod, only : create_co2contra_operator

    class(hordiff_operator_t), allocatable, intent(out) :: hordiff_op
    real(kind=8),                           intent(in)  :: hordiff_coeff
    type(domain_t),                         intent(in)  :: domain

    type(hordiff_c_div_t), allocatable :: hordiff_div

    integer(kind=4) :: halo_width
    real(kind=8)    :: hx

    allocate(hordiff_div)

    !WORKAROUND
    halo_width = 4

    call create_grid_field(hordiff_div%div, halo_width, 0, domain%mesh_p)
    call create_grid_field(hordiff_div%ut,  halo_width, 0, domain%mesh_u)
    call create_grid_field(hordiff_div%vt,  halo_width, 0, domain%mesh_v)

    hordiff_div%div_op = create_div_operator(domain, "divergence_c_sbp21")
    hordiff_div%grad_op = create_grad_operator(domain, "gradient_c_sbp21")
    hordiff_div%co2contra_op = create_co2contra_operator(domain, "co2contra_c_sbp21")

    hx = domain%mesh_p%tile(domain%mesh_p%ts)%hx

    hordiff_div%diff_coeff = hordiff_coeff*domain%mesh_p%scale*hx
    hordiff_div%diff_order = 2

    call move_alloc(hordiff_div, hordiff_op)

end subroutine create_Cgrid_hordiff_div_operator

end module hordiff_factory_mod
