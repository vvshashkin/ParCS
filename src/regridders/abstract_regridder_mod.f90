module abstract_regrid_mod

use grid_field_mod, only : grid_field_t


implicit none



type, abstract :: regrid_t
    contains
        procedure(do_regrid_scalar_interface), deferred :: do_regrid
        procedure(do_regrid_vector_interface), deferred :: do_regrid_vec
end type regrid_t



interface
    subroutine do_regrid_scalar_interface(this,fout,f)
        import grid_field_t, regrid_t
        class(regrid_t),     intent(inout) :: this
        real(kind=8),        intent(inout) :: fout(:,:,:)
        type(grid_field_t),  intent(in)    :: f
    end subroutine do_regrid_scalar_interface

    subroutine do_regrid_vector_interface(this,uout,vout,u,v)
        import grid_field_t, regrid_t
        class(regrid_t),     intent(inout) :: this
        real(kind=8),        intent(inout) :: uout(:,:,:), vout(:,:,:)
        type(grid_field_t),  intent(in)    :: u, v
    end subroutine do_regrid_vector_interface
end interface

end module abstract_regrid_mod
