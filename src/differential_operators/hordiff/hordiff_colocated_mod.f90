module hordiff_colocated_mod

use grid_field_mod,     only : grid_field_t
use halo_mod,           only : halo_vec_t
use domain_mod,         only : domain_t

use abstract_hordiff_mod, only : hordiff_operator_t

implicit none

type, public, extends(hordiff_operator_t) :: hordiff_colocated_t
    class(hordiff_operator_t), allocatable :: hordiff_1comp
    class(halo_vec_t),         allocatable :: edge_sync
contains
    procedure, public :: calc_diff_vec
end type hordiff_colocated_t

contains

subroutine calc_diff_vec(this, u_tend, v_tend, u, v, domain)

    class(hordiff_colocated_t), intent(inout) :: this
    type(grid_field_t),     intent(inout) :: u_tend, v_tend
    type(grid_field_t),     intent(inout) :: u, v
    type(domain_t),         intent(in)    :: domain

    call this%hordiff_1comp%calc_diff(u_tend, u, domain%mesh_u, domain)
    call this%hordiff_1comp%calc_diff(v_tend, v, domain%mesh_v, domain)

    call this%edge_sync%get_halo_vector(u_tend,v_tend,domain,1)

end subroutine calc_diff_vec

end module hordiff_colocated_mod
