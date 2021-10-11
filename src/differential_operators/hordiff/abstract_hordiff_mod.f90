module abstract_hordiff_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

implicit none

type, public, abstract :: horidff_operator_t

contains
    procedure, public :: calc_diff
    procedure, public :: calc_diff_vec
end type horidff_operator_t

contains

subroutine calc_diff(this, f_tend, f, domain)

    class(horidff_operator_t), intent(inout) :: this
    type(grid_field_t),        intent(inout) :: f_tend, f
    type(domain_t),            intent(in)    :: domain


    call domain%parcomm%abort("Calc_diff is not implemented!")

end subroutine

subroutine calc_diff_vec(this, u_tend, v_tend, u, v, domain)

    class(horidff_operator_t), intent(inout) :: this
    type(grid_field_t),        intent(inout) :: u_tend, v_tend, u, v
    type(domain_t),            intent(in)    :: domain


    call domain%parcomm%abort("Calc_diff is not implemented!")

end subroutine

end module abstract_hordiff_mod
