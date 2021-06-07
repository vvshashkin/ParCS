module halo_A_default_mod

use halo_mod,          only : halo_t
use exchange_halo_mod, only : exchange_t

implicit none

type, extends(halo_t) :: halo_A_default_t

    class(exchange_t), allocatable  :: exch_halo

    contains

    procedure :: get_halo_scalar => get_A_default_scalar_halo

end type

contains

subroutine get_A_default_scalar_halo(this,f,parcomm,halo_width)
    use grid_field_mod, only : grid_field_t
    use parcomm_mod,    only : parcomm_t

    class(halo_A_default_t),  intent(inout) :: this
    class(grid_field_t),      intent(inout) :: f
    type(parcomm_t),          intent(in)    :: parcomm
    integer(kind=4),          intent(in)    :: halo_width

    call this%exch_halo%do(f, parcomm)

end subroutine get_A_default_scalar_halo

end module halo_A_default_mod
