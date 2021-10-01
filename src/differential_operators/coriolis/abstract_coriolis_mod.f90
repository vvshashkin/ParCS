module abstract_coriolis_mod

use grid_field_mod,    only : grid_field_t
use domain_mod,        only : domain_t

implicit none

type, abstract, public :: coriolis_operator_t
contains
    procedure(calc_coriolis_i), deferred :: calc_coriolis
end type coriolis_operator_t

abstract interface
    subroutine calc_coriolis_i(this, cor_u, cor_v, ut, vt, domain)
        import coriolis_operator_t, grid_field_t, domain_t
        class(coriolis_operator_t), intent(inout) :: this
        type(domain_t),             intent(in)    :: domain
        type(grid_field_t),         intent(inout) :: ut, vt!contravariant components
        type(grid_field_t),         intent(inout) :: cor_u, cor_v!covariant components
    end subroutine calc_coriolis_i
end interface

contains

end module abstract_coriolis_mod
