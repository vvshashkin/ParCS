module mesh_factory_mod

use mesh_mod, only : mesh_t
use tile_mod, only : tile_t

implicit none

contains

subroutine create_equiangular_mesh(mesh, partition, halo_width, staggering_type, points_type)
    use const_mod,            only: pi
    use ecs_geometry_mod,     only: ecs_ab2xyz_proto, ecs_proto2face,       &
                                    ecs_acov_proto, ecs_bcov_proto,         &
                                    ecs_actv_proto, ecs_bctv_proto,         &
                                    ecs_metric_tensor, ecs_invmetric_tensor,&
                                    ecs_G

    use partition_mod, only : partition_t

    type(partition_t), target, intent(in)  :: partition
    type(mesh_t),              intent(out) :: mesh

    integer(kind=4),  intent(in)   :: halo_width
    character(len=*), intent(in)   :: staggering_type, points_type

    integer(kind=4) :: t, i, j, k, ts, te, is, ie, js, je, ks, ke, nh, nx, ny
    real(kind=8) :: alpha, beta, x, y, z, r
    real(kind=8) :: xyz(3)

    type(tile_t), pointer :: tile(:)

    real(kind=8) :: i_0, j_0, alpha_0, beta_0

    alpha_0 = -0.25_8*pi
    beta_0  = -0.25_8*pi

    i_0 = 0.5_8
    j_0 = 0.5_8

    select case(points_type)
    case('p')
        tile => partition%tile
    case('u')
        if (staggering_type=='C') i_0 = 0.0_8
        tile => partition%tile_u
    case('v')
        if (staggering_type=='C') j_0 = 0.0_8
        tile => partition%tile_v
    case default
        print*, 'Error! Wrong points_type! Abort!'
        stop
    end select

    ts = partition%ts
    te = partition%te

    allocate(mesh%tile(ts:te))
    mesh%ts = ts
    mesh%te = te

    do t = ts, te

        ks = tile(t)%ks; ke = tile(t)%ke
        js = tile(t)%js; je = tile(t)%je;
        is = tile(t)%is; ie = tile(t)%ie;

        nh = partition%nh

        call mesh%tile(t)%init(is, ie, js, je, ks, ke, halo_width)

        mesh%tile(t)%hx = 0.5_8 * pi / (real(nh, 8))

        mesh%tile(t)%i_0 = i_0
        mesh%tile(t)%j_0 = j_0

        mesh%tile(t)%alpha_0 = alpha_0
        mesh%tile(t)%beta_0  = beta_0

        do j = js - halo_width, je + halo_width
            beta = mesh%tile(t)%get_beta(j)
            do i = is - halo_width, ie + halo_width
                alpha = mesh%tile(t)%get_beta(i)
                !cartesian coordinates of (alpha,beta) point
                xyz = ecs_ab2xyz_proto(alpha, beta)
                xyz = ecs_proto2face(xyz, partition%panel_map(t))
                mesh%tile(t)%rhx(i,j) = xyz(1)
                mesh%tile(t)%rhy(i,j) = xyz(2)
                mesh%tile(t)%rhz(i,j) = xyz(3)
                !curvilinear system basis vectors
                xyz = ecs_acov_proto(alpha,beta)
                mesh%tile(t)%acov(1:3,i,j) = ecs_proto2face(xyz, partition%panel_map(t))

                xyz = ecs_bcov_proto(alpha,beta)
                mesh%tile(t)%bcov(1:3,i,j) = ecs_proto2face(xyz, partition%panel_map(t))

                xyz = ecs_actv_proto(alpha,beta)
                mesh%tile(t)%actv(1:3,i,j) = ecs_proto2face(xyz, partition%panel_map(t))

                xyz = ecs_bctv_proto(alpha,beta)
                mesh%tile(t)%bctv(1:3,i,j) = ecs_proto2face(xyz, partition%panel_map(t))

                mesh%tile(t)%Q (1:3,i,j) = ecs_metric_tensor(alpha, beta)
                mesh%tile(t)%Qi(1:3,i,j) = ecs_invmetric_tensor(alpha, beta)
                mesh%tile(t)%G(i,j)      = ecs_G(alpha, beta)
            end do
        end do
    end do


end subroutine create_equiangular_mesh


end module mesh_factory_mod
