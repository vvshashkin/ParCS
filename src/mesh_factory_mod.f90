module mesh_factory_mod

use mesh_mod, only : mesh_t
use tile_mod, only : tile_t

implicit none

contains

subroutine create_mesh(mesh, partition, metric, halo_width, points_type)

    use partition_mod, only : partition_t
    use metric_mod,    only : metric_t

    type(partition_t), target, intent(in)  :: partition
    class(metric_t),           intent(in)  :: metric
    type(mesh_t),              intent(out) :: mesh

    integer(kind=4),  intent(in)   :: halo_width
    character(len=*), intent(in)   :: points_type

    integer(kind=4) :: t, pind, i, j, k, ts, te, is, ie, js, je, ks, ke, nh, nx, ny
    real(kind=8) :: alpha, beta, vec(3)

    type(tile_t), pointer :: tile(:)

    real(kind=8) :: shift_i, shift_j, alpha_0, beta_0

    alpha_0 = metric%alpha0
    beta_0  = metric%beta0

    shift_i = 0.5_8
    shift_j = 0.5_8

    select case(points_type)
    case('c')
        tile => partition%tile
    case('x')
        shift_i = 0.0_8
        tile => partition%tile_x
    case('y')
        shift_j = 0.0_8
        tile => partition%tile_y
    case('xy')
        shift_i = 0.0_8
        shift_j = 0.0_8
        tile => partition%tile_xy
    case default
        print*, 'Error! Wrong points_type! Abort!'
        stop
    end select

    ts = partition%ts
    te = partition%te

    allocate(mesh%tile(ts:te))
    mesh%ts = ts
    mesh%te = te

    mesh%scale         = metric%scale
    mesh%omega         = metric%omega
    mesh%rotation_axis = metric%rotation_axis

    nh = partition%nh

    do t = ts, te

        ks = tile(t)%ks; ke = tile(t)%ke
        js = tile(t)%js; je = tile(t)%je;
        is = tile(t)%is; ie = tile(t)%ie;
        pind = partition%panel_map(t)

        call mesh%tile(t)%init(is, ie, js, je, ks, ke, halo_width)

        !WORKAROUND
        mesh%tile(t)%nx = nh
        mesh%tile(t)%ny = nh

        mesh%tile(t)%hx = (metric%alpha1 - metric%alpha0)/real(nh,8)
        mesh%tile(t)%hy = (metric%beta1  - metric%beta0 )/real(nh,8)

        mesh%tile(t)%shift_i = shift_i
        mesh%tile(t)%shift_j = shift_j

        mesh%tile(t)%alpha_0 = alpha_0
        mesh%tile(t)%beta_0  = beta_0

        do j = js - halo_width, je + halo_width
            beta = mesh%tile(t)%get_beta(j)

            do i = is - halo_width, ie + halo_width
                alpha = mesh%tile(t)%get_alpha(i)

                vec(1:3) = metric%point_r(pind, alpha, beta)

                mesh%tile(t)%rx(i,j) = vec(1)
                mesh%tile(t)%ry(i,j) = vec(2)
                mesh%tile(t)%rz(i,j) = vec(3)

                mesh%tile(t)%a1(1:3,i,j) = metric%a1(pind, alpha, beta)
                mesh%tile(t)%a2(1:3,i,j) = metric%a2(pind, alpha, beta)
                mesh%tile(t)%b1(1:3,i,j) = metric%b1(pind, alpha, beta)
                mesh%tile(t)%b2(1:3,i,j) = metric%b2(pind, alpha, beta)

                mesh%tile(t)%Q (1:3,i,j) = metric%Q (pind, alpha, beta)
                mesh%tile(t)%Qi(1:3,i,j) = metric%QI(pind, alpha, beta)

                mesh%tile(t)%G(i,j)= metric%G(pind, alpha, beta)
            end do
        end do
    end do


end subroutine create_mesh


end module mesh_factory_mod
