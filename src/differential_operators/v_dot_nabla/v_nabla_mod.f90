module v_nabla_mod

use abstract_v_nabla_mod, only : v_nabla_operator_t
use grid_field_mod,       only : grid_field_t, tile_field_t
use mesh_mod,             only : mesh_t, tile_mesh_t

implicit none

type, public, extends(v_nabla_operator_t) :: v_nabla_up4_operator_t
contains
    procedure, public  :: calc_v_nabla => calc_v_nabla_up4
    procedure, private, nopass :: calc_v_nabla_tile => calc_v_nabla_up4_tile
end type v_nabla_up4_operator_t

contains

subroutine calc_v_nabla_up4(this, f_tend, f, ut, vt, mesh)

    class(v_nabla_up4_operator_t), intent(inout) :: this
    type(grid_field_t),            intent(in)    :: f      !advected field
    type(grid_field_t),            intent(in)    :: ut, vt !contravariant components
    type(mesh_t),                  intent(in)    :: mesh
    type(grid_field_t),            intent(inout) :: f_tend !advective tendency

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%calc_v_nabla_tile(f_tend%tile(t), f%tile(t), ut%tile(t), vt%tile(t), &
                                    mesh%tile(t), mesh%scale)
    end do
end subroutine calc_v_nabla_up4

subroutine calc_v_nabla_up4_tile(f_tend, f, ut, vt, mesh, scale)

    type(tile_field_t),            intent(in)    :: f      !advected field
    type(tile_field_t),            intent(in)    :: ut, vt !contravariant components
    type(tile_mesh_t),             intent(in)    :: mesh
    real(kind=8),                  intent(in)    :: scale
    type(tile_field_t),            intent(inout) :: f_tend !advective tendency

    integer(kind=4) :: i, j, k
    integer(kind=4) :: is, ie, js, je, ks, ke
    real(kind=8)    :: hx
    real(kind=8)    :: dx, dy, zl, zr

    hx = mesh%hx

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k = ks, ke
        do j = js, je
            do i = is, ie
                zl = .5_8+sign(.5_8,ut%p(i,j,k))
                zr = 1._8-zl
                dx = (zl*( 3._8*f%p(i+1,j,k)+10._8*f%p(i,j,k)- &
                          18._8*f%p(i-1,j,k)+6._8*f%p(i-2,j,k)-f%p(i-3,j,k))-&
                      zr*( 3._8*f%p(i-1,j,k)+10._8*f%p(i,j,k)- &
                          18._8*f%p(i+1,j,k)+6._8*f%p(i+2,j,k)-f%p(i+3,j,k))) / 12.0_8
                ! dx = (zl*( 2._8*uvt%p(i+1,j,k)+3._8*uvt%p(i,j,k)- &
                !           6._8*uvt%p(i-1,j,k)+1._8*uvt%p(i-2,j,k))+&
                !       zr*( 2._8*uvt%p(i-1,j,k)+3._8*uvt%p(i,j,k)- &
                !           6._8*uvt%p(i+1,j,k)+1._8*uvt%p(i+2,j,k))) / 6.0_8
                ! dx = (zl*(uvt%p(i,j,k)-uvt%p(i-1,j,k))+zr*(uvt%p(i+1,j,k)-uvt%p(i,j,k)))
                ! dx = 0.5_8*(uvt%p(i+1,j,k)-uvt%p(i-1,j,k))
                ! dx = (-uvt%p(i+2,j,k)+8._8*uvt%p(i+1,j,k)-8._8*uvt%p(i-1,j,k)+uvt%p(i-2,j,k)) / 12._8
                ! dx = (uvt%p(i+3,j,k)/60._8-3._8/20._8*uvt%p(i+2,j,k)+3._8/4._8*uvt%p(i+1,j,k)-&
                !       3._8/4._8*uvt%p(i-1,j,k)+3._8/20._8*uvt%p(i-2,j,k)-uvt%p(i-3,j,k)/60._8)
                zl = .5_8+sign(.5_8,vt%p(i,j,k))
                zr = 1._8-zl
                dy = (zl*( 3._8*f%p(i,j+1,k)+10._8*f%p(i,j,k)- &
                          18._8*f%p(i,j-1,k)+6._8*f%p(i,j-2,k)-f%p(i,j-3,k))-&
                      zr*( 3._8*f%p(i,j-1,k)+10._8*f%p(i,j,k)- &
                          18._8*f%p(i,j+1,k)+6._8*f%p(i,j+2,k)-f%p(i,j+3,k))) / 12.0_8
                ! dy = (zl*( 2._8*uvt%p(i,j+1,k)+3._8*uvt%p(i,j,k)- &
                !           6._8*uvt%p(i,j-1,k)+1._8*uvt%p(i,j-2,k))+&
                !       zr*( 2._8*uvt%p(i,j-1,k)+3._8*uvt%p(i,j,k)- &
                !           6._8*uvt%p(i,j+1,k)+1._8*uvt%p(i,j+2,k))) / 6.0_8
                ! dy = (zl*(uvt%p(i,j,k)-uvt%p(i,j-1,k))+zr*(uvt%p(i,j+1,k)-uvt%p(i,j,k)))
                ! dy = 0.5_8*(uvt%p(i,j+1,k)-uvt%p(i,j-1,k))
                !dy = (-uvt%p(i,j+2,k)+8._8*uvt%p(i,j+1,k)-8._8*uvt%p(i,j-1,k)+uvt%p(i,j-2,k)) / 12._8
                ! dy = (uvt%p(i,j+3,k)/60._8-3._8/20._8*uvt%p(i,j+2,k)+3._8/4._8*uvt%p(i,j+1,k)-&
                !       3._8/4._8*uvt%p(i,j-1,k)+3._8/20._8*uvt%p(i,j-2,k)-uvt%p(i,j-3,k)/60._8)
                f_tend%p(i,j,k) =-(ut%p(i,j,k)*dx+vt%p(i,j,k)*dy) / (hx*scale)
            end do
        end do
    end do

end subroutine calc_v_nabla_up4_tile

end module v_nabla_mod
