module halo_mod

use grid_function_mod, only : grid_function_t

implicit none

type, abstract :: halo_t
    contains
    procedure(halo_interp),  deferred, public :: interp  !scalar halo procedure
    !procedure(halo_interpv), deferred, public :: interpv !vector halo procedure
end type halo_t

type, abstract :: halo_vec_t
    contains
    procedure(halo_interpv), deferred, public :: interpv !vector halo procedure
end type halo_vec_t

interface
    subroutine halo_interp(this,f,halo_width)
        import halo_t
        import grid_function_t
        class(halo_t),         intent(in)    :: this
        type(grid_function_t), intent(inout) :: f
        integer(kind=4),       intent(in)    :: halo_width
    end subroutine halo_interp

    subroutine halo_interpv(this,u,v,halo_width)
        import halo_vec_t
        import grid_function_t
        class(halo_vec_t),    intent(in)    :: this
        type(grid_function_t),intent(inout) :: u, v
        integer(kind=4),       intent(in)    :: halo_width
    end subroutine halo_interpv
end interface

end module halo_mod
