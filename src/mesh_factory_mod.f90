module mesh_factory_mod

use mesh_mod,      only : mesh_t
use topology_mod,  only : ex, ey, n

implicit none


contains


subroutine create_equiangular_mesh(mesh, is, ie, js, je, ks, ke, nh, halo_width, panel_ind)
    use const_mod,            only: pi
    use ecs_geometry_mod,     only: ecs_ab2xyz_proto, ecs_proto2face,       &
                                    ecs_acov_proto, ecs_bcov_proto,         &
                                    ecs_actv_proto, ecs_bctv_proto,         &
                                    ecs_metric_tensor, ecs_invmetric_tensor,&
                                    ecs_G
    use ecs_halo_factory_mod, only: init_ecs_halo

    type(mesh_t),    intent(inout) :: mesh
    integer(kind=4), intent(in)    :: is, ie, js, je, ks, ke
    integer(kind=4), intent(in)    :: nh, halo_width
    integer(kind=4), intent(in)    :: panel_ind

    integer(kind=4) :: i, j, k
    real(kind=8) :: alpha, beta, x, y, z, r
    real(kind=8) :: xyz(3)

    call mesh%init(is, ie, js, je, ks, ke, halo_width)

    mesh%panel_ind = panel_ind
    mesh%hx = 0.5_8 * pi / (real(nh, 8))
    mesh%nx = nh

    do j = js - halo_width, je + halo_width
        beta = -0.25_8*pi + (j-0.5_8)*mesh.hx
        do i = is - halo_width, ie + halo_width
            alpha = -0.25_8*pi+(i-0.5_8)*mesh%hx
            !cartesian coordinates of (alpha,beta) point
            xyz = ecs_ab2xyz_proto(alpha,beta)
            xyz = ecs_proto2face(xyz,panel_ind)
            mesh%rhx(i,j) = xyz(1)
            mesh%rhy(i,j) = xyz(2)
            mesh%rhz(i,j) = xyz(3)
            !curvilinear system basis vectors
            xyz = ecs_acov_proto(alpha,beta)
            mesh%acov(:,i,j) = ecs_proto2face(xyz,panel_ind)
            xyz = ecs_bcov_proto(alpha,beta)
            mesh%bcov(:,i,j) = ecs_proto2face(xyz,panel_ind)
            xyz = ecs_actv_proto(alpha,beta)
            mesh%actv(:,i,j) = ecs_proto2face(xyz,panel_ind)
            xyz = ecs_bctv_proto(alpha,beta)
            mesh%bctv(:,i,j) = ecs_proto2face(xyz,panel_ind)

            mesh%Q(:,i,j)  = ecs_metric_tensor(alpha,beta)
            mesh%Qi(:,i,j) = ecs_invmetric_tensor(alpha,beta)
            mesh%G(i,j)    = ecs_G(alpha, beta)
        end do
    end do
    do j=js-halo_width, je+halo_width
        beta = -0.25_8*pi + (j-0.5_8)*mesh.hx
        do i=is-halo_width-1, ie+halo_width
            alpha = -0.25_8*pi+i*mesh%hx
            mesh%Gu(i,j)    = ecs_G(alpha, beta)
        end do
    end do
    do j=js-halo_width-1, je+halo_width
        beta = -0.25_8*pi + j*mesh.hx
        do i=is-halo_width, ie+halo_width
            alpha = -0.25_8*pi+(i-0.5_8)*mesh%hx
            mesh%Gv(i,j)    = ecs_G(alpha, beta)
        end do
    end do
    !mesh%halo = init_ecs_halo(mesh%is, mesh%ie,           &
    !                          mesh%js, mesh%je,           &
    !                          mesh%nx, halo_width,   &
    !                          mesh%hx)

end subroutine create_equiangular_mesh


end module mesh_factory_mod
