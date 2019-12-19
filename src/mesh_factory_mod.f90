module mesh_factory_mod

use mesh_mod,      only : mesh_t
use topology_mod,  only : ex, ey, n

implicit none


contains


subroutine create_equiangular_mesh(mesh, is, ie, js, je, ks, ke, nh, halo_width, panel_ind)

    type(mesh_t),    intent(inout) :: mesh
    integer(kind=4), intent(in)    :: is, ie, js, je, ks, ke
    integer(kind=4), intent(in)    :: nh, halo_width
    integer(kind=4), intent(in)    :: panel_ind

    integer(kind=4) :: i, j, k
    real(kind=8) :: alpha, beta, x, y, z, r
    real(kind=8) :: pi

    pi = acos(-1.0_8)

    call mesh%init(is, ie, js, je, ks, ke, halo_width)

    mesh%panel_ind = panel_ind
    mesh%hx = 0.5_8 * pi / (real(nh, 8))
    mesh%nx = nh

    do j = js - halo_width, je + halo_width
        beta = -0.25_8*pi + (j-0.5_8)*mesh.hx
        do i = is - halo_width, ie + halo_width
            alpha = -0.25_8*pi+(i-0.5_8)*mesh%hx

            !grid point at cube face: (assumes x = tan(alpha), y = tan(beta), z = 1)
            r = sqrt(1.0_8 + tan(alpha)**2 + tan(beta)**2)

            !coordinates at prototype spherical face: vec(r)/||r||, r = (x, y, z)
            x = tan(alpha)/r; y = tan(beta)/r; z = 1.0_8/r

            !transform to real spherical face
            mesh%rhx(i,j) = ex(1,panel_ind)*x + ey(1,panel_ind)*y - n(1,panel_ind)*z
            mesh%rhy(i,j) = ex(2,panel_ind)*x + ey(2,panel_ind)*y - n(2,panel_ind)*z
            mesh%rhz(i,j) = ex(3,panel_ind)*x + ey(3,panel_ind)*y - n(3,panel_ind)*z

        end do
    end do

end subroutine create_equiangular_mesh


end module mesh_factory_mod
