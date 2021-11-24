module abstract_interpolators3d_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

implicit none

type, abstract :: interpolator_w2uv_t
contains
    procedure(interp_w2uv), deferred :: interp_w2uv
end type interpolator_w2uv_t

abstract interface
    subroutine interp_w2uv(this, wu, wv, w, domain)
        import interpolator_w2uv_t, grid_field_t, domain_t
        class(interpolator_w2uv_t), intent(inout) :: this
        type(grid_field_t),         intent(inout) :: w
        type(domain_t),             intent(in)    :: domain
        !output:
        type(grid_field_t),         intent(inout) :: wu, wv
    end subroutine interp_w2uv
end interface

end module abstract_interpolators3d_mod
