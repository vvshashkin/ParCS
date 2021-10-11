module hordiff_colocated_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

use abstract_hordiff_mod, only : hordiff_operator_t

implicit none

type, public, extends(hordiff_operator_t) :: hordiff_colocated_t
contains
    procedure, public :: calc_diff_vec
end type hordiff_colocated_t

contains

subroutine calc_diff_vec(this, u_tend, v_tend, u, v, domain)

    class(hordiff_colocated_t), intent(inout) :: this
    type(grid_field_t),     intent(inout) :: u_tend, v_tend
    type(grid_field_t),     intent(inout) :: u, v
    type(domain_t),         intent(in)    :: domain

    call u_tend%assign(0.0_8, domain%mesh_u)
    call v_tend%assign(0.0_8, domain%mesh_v)

end subroutine calc_diff_vec

end module hordiff_colocated_mod
