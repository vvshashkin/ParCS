module halo_mod

use grid_field_mod, only : grid_field_t
use parcomm_mod,    only : parcomm_t

implicit none

type, abstract :: halo_t
    contains
    procedure(get_halo_scalar),  deferred, public :: get_halo_scalar  !scalar halo procedure
end type halo_t

type, abstract :: halo_vec_t
    contains
    procedure(get_halo_vector),  deferred, public :: get_halo_vector  !scalar halo procedure
end type halo_vec_t

interface
    subroutine get_halo_scalar(this,f,parcomm,halo_width)
        import halo_t, grid_field_t, parcomm_t
        class(halo_t),       intent(inout) :: this
        class(grid_field_t), intent(inout) :: f
        type(parcomm_t),     intent(in)    :: parcomm
        integer(kind=4),     intent(in)    :: halo_width
    end subroutine get_halo_scalar

    subroutine get_halo_vector(this,u,v,parcomm,halo_width)
        import halo_vec_t, grid_field_t, parcomm_t
        class(halo_vec_t),   intent(inout) :: this
        class(grid_field_t), intent(inout) :: u, v
        type(parcomm_t),     intent(in)    :: parcomm
        integer(kind=4),     intent(in)    :: halo_width
    end subroutine get_halo_vector
end interface

!Old code
!type, abstract :: halo_t
!    contains
!    procedure(halo_interp),  deferred, public :: interp  !scalar halo procedure
!    !procedure(halo_interpv), deferred, public :: interpv !vector halo procedure
!end type halo_t
!
!type, abstract :: halo_vec_t
!    contains
!    procedure(halo_interpv), deferred, public :: interpv !vector halo procedure
!end type halo_vec_t
!
!interface
!    subroutine halo_interp(this,f,halo_width)
!        import halo_t, block_t
!        class(halo_t),   intent(in)    :: this
!        type(block_t),   intent(inout) :: f
!        integer(kind=4), intent(in)    :: halo_width
!    end subroutine halo_interp
!
!    subroutine halo_interpv(this,u,v,halo_width)
!        import halo_vec_t, block_t
!        class(halo_vec_t), intent(in)    :: this
!        type(block_t),     intent(inout) :: u, v
!        integer(kind=4),   intent(in)    :: halo_width
!    end subroutine halo_interpv
!end interface

end module halo_mod
