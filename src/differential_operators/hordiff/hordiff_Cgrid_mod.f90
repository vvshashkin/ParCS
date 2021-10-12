module hordiff_Cgrid_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

use abstract_hordiff_mod,   only : hordiff_operator_t
use abstract_div_mod,       only : div_operator_t
use abstract_grad_mod,      only : grad_operator_t
use abstract_co2contra_mod, only : co2contra_operator_t

implicit none

type, public, extends(hordiff_operator_t) :: hordiff_c_div_t

    integer(kind=4) :: diff_order
    real(kind=8)    :: diff_coeff

    class(div_operator_t),       allocatable :: div_op
    class(grad_operator_t),      allocatable :: grad_op
    class(co2contra_operator_t), allocatable :: co2contra_op

    type(grid_field_t) :: div, ut, vt

contains
    procedure, public :: calc_diff_vec => calc_diff_vec_div_part
end type hordiff_c_div_t

type, public, extends(hordiff_operator_t) :: hordiff_c_curl_dump_t

    integer(kind=4) :: diff_order
    real(kind=8)    :: diff_coeff

    class(div_operator_t),       allocatable :: curl_op
    class(grad_operator_t),      allocatable :: grad_perp_op
    class(co2contra_operator_t), allocatable :: co2contra_op

    type(grid_field_t) :: curl, ut, vt

contains
    procedure, public :: calc_diff_vec => calc_diff_vec_curl_part
end type hordiff_c_curl_dump_t

contains

subroutine calc_diff_vec_div_part(this, u_tend, v_tend, u, v, domain)

    class(hordiff_c_div_t), intent(inout) :: this
    type(grid_field_t),     intent(inout) :: u_tend, v_tend
    type(grid_field_t),     intent(inout) :: u, v
    type(domain_t),         intent(in)    :: domain

    integer(kind=4) :: p

    real(kind=8) :: coeff

    call this%co2contra_op%transform(this%ut, this%vt, u, v, domain)

    call this%div_op%calc_div(this%div, this%ut, this%vt, domain)

    do p = 2, this%diff_order

        call this%grad_op%calc_grad(u_tend, v_tend, this%div, domain)

        call this%co2contra_op%transform(this%ut, this%vt, u_tend, v_tend, domain)

        call this%div_op%calc_div(this%div, this%ut, this%vt, domain)

    end do

    call this%grad_op%calc_grad(this%ut, this%vt, this%div, domain)

    coeff = (-1.0_8)**(this%diff_order+1)*this%diff_coeff**(2*this%diff_order)

    call u_tend%assign(coeff, this%ut, domain%mesh_u)
    call v_tend%assign(coeff, this%vt, domain%mesh_v)

end subroutine calc_diff_vec_div_part

subroutine calc_diff_vec_curl_part(this, u_tend, v_tend, u, v, domain)

    class(hordiff_c_div_t), intent(inout) :: this
    type(grid_field_t),          intent(inout) :: u_tend, v_tend
    type(grid_field_t),          intent(inout) :: u, v
    type(domain_t),              intent(in)    :: domain

    integer(kind=4) :: p

    real(kind=8) :: coeff

    call this%co2contra_op%transform(this%ut, this%vt, u, v, domain)

    call this%curl_op%calc_curl(this%curl, u, v, domain)

    do p = 2, this%diff_order

        call this%grad_perp_op%calc_grad_perp(u_tend, v_tend, this%div, domain)

        call this%co2contra_op%transform(this%ut, this%vt, u_tend, v_tend, domain)

        call this%div_op%calc_div(this%div, this%ut, this%vt, domain)

    end do

    call this%grad_op%calc_grad(this%ut, this%vt, this%div, domain)

    coeff = (-1.0_8)**(this%diff_order+1)*this%diff_coeff**(2*this%diff_order)

    call u_tend%assign(coeff, this%ut, domain%mesh_u)
    call v_tend%assign(coeff, this%vt, domain%mesh_v)

end subroutine calc_diff_vec_curl_part

end module hordiff_Cgrid_mod
