module mesh_factory_mod

use mesh_mod,      only : mesh_t
use topology_mod,  only : ex, ey, n

implicit none

contains

! subroutine create_equiangular_mesh(mesh, is, ie, js, je, ks, ke, nh, halo_width, panel_ind)
!     use const_mod,            only: pi
!     use ecs_geometry_mod,     only: ecs_ab2xyz_proto, ecs_proto2face,       &
!                                     ecs_acov_proto, ecs_bcov_proto,         &
!                                     ecs_actv_proto, ecs_bctv_proto,         &
!                                     ecs_metric_tensor, ecs_invmetric_tensor,&
!                                     ecs_G
!     use ecs_halo_factory_mod, only: init_ecs_halo
!
!     type(mesh_t),    intent(inout) :: mesh
!     integer(kind=4), intent(in)    :: is, ie, js, je, ks, ke
!     integer(kind=4), intent(in)    :: nh, halo_width
!     integer(kind=4), intent(in)    :: panel_ind
!
!     integer(kind=4) :: i, j, k
!     real(kind=8) :: alpha, beta, x, y, z, r
!     real(kind=8) :: xyz(3)
!
!     call mesh%init(is, ie, js, je, ks, ke, halo_width)
!
!     mesh%panel_ind = panel_ind
!     mesh%hx = 0.5_8 * pi / (real(nh, 8))
!     mesh%nx = nh
!
!     do j = js - halo_width, je + halo_width
!         beta = -0.25_8*pi + (j-0.5_8)*mesh.hx
!         do i = is - halo_width, ie + halo_width
!             alpha = -0.25_8*pi+(i-0.5_8)*mesh%hx
!             !cartesian coordinates of (alpha,beta) point
!             xyz = ecs_ab2xyz_proto(alpha,beta)
!             xyz = ecs_proto2face(xyz,panel_ind)
!             mesh%rhx(i,j) = xyz(1)
!             mesh%rhy(i,j) = xyz(2)
!             mesh%rhz(i,j) = xyz(3)
!             !curvilinear system basis vectors
!             xyz = ecs_acov_proto(alpha,beta)
!             mesh%acov(:,i,j) = ecs_proto2face(xyz,panel_ind)
!             xyz = ecs_bcov_proto(alpha,beta)
!             mesh%bcov(:,i,j) = ecs_proto2face(xyz,panel_ind)
!             xyz = ecs_actv_proto(alpha,beta)
!             mesh%actv(:,i,j) = ecs_proto2face(xyz,panel_ind)
!             xyz = ecs_bctv_proto(alpha,beta)
!             mesh%bctv(:,i,j) = ecs_proto2face(xyz,panel_ind)
!
!             mesh%Q(:,i,j)  = ecs_metric_tensor(alpha,beta)
!             mesh%Qi(:,i,j) = ecs_invmetric_tensor(alpha,beta)
!             mesh%G(i,j)    = ecs_G(alpha, beta)
!         end do
!     end do
!     do j=js-halo_width, je+halo_width
!         beta = -0.25_8*pi + (j-0.5_8)*mesh.hx
!         do i=is-halo_width-1, ie+halo_width
!             alpha = -0.25_8*pi+i*mesh%hx
!             mesh%Gu(i,j)    = ecs_G(alpha, beta)
!             mesh%Qiu(:,i,j) = ecs_invmetric_tensor(alpha,beta)
!         end do
!     end do
!     do j=js-halo_width-1, je+halo_width
!         beta = -0.25_8*pi + j*mesh.hx
!         do i=is-halo_width, ie+halo_width
!             alpha = -0.25_8*pi+(i-0.5_8)*mesh%hx
!             mesh%Gv(i,j)    = ecs_G(alpha, beta)
!             mesh%Qiv(:,i,j) = ecs_invmetric_tensor(alpha,beta)
!         end do
!     end do
!     !mesh%halo = init_ecs_halo(mesh%is, mesh%ie,           &
!     !                          mesh%js, mesh%je,           &
!     !                          mesh%nx, halo_width,   &
!     !                          mesh%hx)
!
! end subroutine create_equiangular_mesh


subroutine create_equiangular_mesh(mesh, partition, halo_width, staggering_type, points_type)
    use const_mod,            only: pi
    use ecs_geometry_mod,     only: ecs_ab2xyz_proto, ecs_proto2face,       &
                                    ecs_acov_proto, ecs_bcov_proto,         &
                                    ecs_actv_proto, ecs_bctv_proto,         &
                                    ecs_metric_tensor, ecs_invmetric_tensor,&
                                    ecs_G

    use partition_mod, only : partition_t

    type(mesh_t), allocatable, intent(out) :: mesh(:)
    type(partition_t),         intent(in)  :: partition

    integer(kind=4), intent(in)            :: halo_width
    character(len=*), intent(in)           :: staggering_type, points_type

    integer(kind=4) :: t, i, j, k, ts, te, is, ie, js, je, ks, ke, nh, nx, ny
    real(kind=8) :: alpha, beta, x, y, z, r
    real(kind=8) :: xyz(3)

    real(kind=8) :: i_0, j_0, alpha_0, beta_0

    alpha_0 = -0.25_8*pi
    beta_0  = -0.25_8*pi

    select case(staggering_type)
    case('A')
        i_0 = 0.5_8
        j_0 = 0.5_8
        nx = partition%nh
        ny = partition%nh
    case('C')
        select case(points_type)
        case('p')
            i_0 = 0.5_8
            j_0 = 0.5_8
            nx = partition%nh
            ny = partition%nh
        case('u')
            i_0 = 0.0_8
            j_0 = 0.5_8
            nx = partition%nh+1
            ny = partition%nh
        case('v')
            i_0 = 0.5_8
            j_0 = 0.0_8
            nx = partition%nh
            ny = partition%nh+1
        case default
            print*, 'Error! Wrong points_type! Abort!'
            stop
        end select
    case default
        print*, 'Error! Wrong staggering_type! Abort!'
        stop
    end select


    ts = partition%ts
    te = partition%te

    allocate(mesh(ts:te))

    do t = ts, te

        ks = partition%tile(t)%ks; ke = partition%tile(t)%ke;

        js = partition%tile(t)%js; je = partition%tile(t)%je;
        if (je == partition%nh .and. ny == partition%nh+1) je = ny

        is = partition%tile(t)%is; ie = partition%tile(t)%ie;
        if (ie == partition%nh .and. nx == partition%nh+1) ie = nx

        nh = partition%nh

        call mesh(t)%init(is, ie, js, je, ks, ke, halo_width)

        ! mesh(t)%panel_ind = partition%tile(t)%panel_number
        mesh(t)%hx = 0.5_8 * pi / (real(nh, 8))

        mesh(t)%i_0 = i_0
        mesh(t)%j_0 = j_0

        mesh(t)%alpha_0 = alpha_0
        mesh(t)%beta_0  = beta_0

        do j = js - halo_width, je + halo_width
            beta = mesh(t)%get_beta(j)
            do i = is - halo_width, ie + halo_width
                alpha = mesh(t)%get_beta(i)
                !cartesian coordinates of (alpha,beta) point
                xyz = ecs_ab2xyz_proto(alpha, beta)
                xyz = ecs_proto2face(xyz, partition%tile(t)%panel_number)
                mesh(t)%rhx(i,j) = xyz(1)
                mesh(t)%rhy(i,j) = xyz(2)
                mesh(t)%rhz(i,j) = xyz(3)
                !curvilinear system basis vectors
                xyz = ecs_acov_proto(alpha,beta)
                mesh(t)%acov(1:3,i,j) = ecs_proto2face(xyz, partition%tile(t)%panel_number)
                xyz = ecs_bcov_proto(alpha,beta)
                mesh(t)%bcov(1:3,i,j) = ecs_proto2face(xyz, partition%tile(t)%panel_number)
                xyz = ecs_actv_proto(alpha,beta)
                mesh(t)%actv(1:3,i,j) = ecs_proto2face(xyz, partition%tile(t)%panel_number)
                xyz = ecs_bctv_proto(alpha,beta)
                mesh(t)%bctv(1:3,i,j) = ecs_proto2face(xyz, partition%tile(t)%panel_number)

                mesh(t)%Q (1:3,i,j) = ecs_metric_tensor(alpha, beta)
                mesh(t)%Qi(1:3,i,j) = ecs_invmetric_tensor(alpha, beta)
                mesh(t)%G(i,j)      = ecs_G(alpha, beta)
            end do
        end do
    end do


end subroutine create_equiangular_mesh


end module mesh_factory_mod
