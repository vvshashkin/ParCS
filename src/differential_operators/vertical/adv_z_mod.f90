module adv_z_mod

use abstract_adv_z_mod, only : adv_z_t
use grid_field_mod,     only : grid_field_t, tile_field_t
use mesh_mod,           only : mesh_t, tile_mesh_t

implicit none

type, extends(adv_z_t) :: adv_z_c2_t
    contains
    procedure :: calc_z_adv_tile => calc_z_adv_c2_tile
end type adv_z_c2_t

contains

subroutine calc_z_adv_c2_tile(this, f_tend, f, eta_dot, mesh,scale)
    class(adv_z_c2_t),  intent(in)    :: this
    type(tile_field_t), intent(in)    :: f, eta_dot
    type(tile_mesh_t),  intent(in)    :: mesh
    real(kind=8),       intent(in)    :: scale
    !output
    type(tile_field_t), intent(inout) :: f_tend

    integer(kind=4) :: is, ie, js, je, ks, ke
    integer(kind=4) :: i, j, k
    integer(kind=4) :: km1, kp1

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks, ke
        !constant extension of field above/below the the upper/lower boundary
        km1 = max(k-1,ks)
        kp1 = min(k+1,ke)
        do j=js, je
            do i=is, ie
                f_tend%p(i,j,k) =-eta_dot%p(i,j,k)*(f%p(i,j,kp1)-f%p(i,j,km1)) / (2.0_8*mesh%hz*scale)
            end do
        end do
    end do
end subroutine calc_z_adv_c2_tile

end module adv_z_mod
