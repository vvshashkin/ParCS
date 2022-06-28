module mesh_factory_mod

use mesh_mod,       only : mesh_t, tile_mesh_t
use tiles_mod,      only : tiles_t
use parcomm_mod,    only : parcomm_global
use orography_mod,  only : orography_1mesh_t
use partition_mod,  only : partition_t
use metric_mod,     only : metric_t
use grid_field_mod, only : tile_field_t

implicit none

contains

subroutine create_mesh(mesh, partition, metric, halo_width, h_top, tiles, &
                       shift_xyz, orography)

    type(partition_t),       intent(in)   :: partition
    class(metric_t),         intent(in)   :: metric
    type(mesh_t),            intent(out)  :: mesh

    integer(kind=4),         intent(in)   :: halo_width
    real(kind=8),            intent(in)   :: h_top
    type(tiles_t),           intent(in)   :: tiles
    real(kind=8),            intent(in)   :: shift_xyz(3)
    type(orography_1mesh_t), intent(in), optional :: orography

    integer(kind=4) :: t, pind, i, j, k, ts, te, is, ie, js, je, ks, ke, nh, nx, ny, nz
    real(kind=8)    :: hx, hy, hz


    real(kind=8) :: shift_i, shift_j, shift_k

    ! define horizontal grid parameters
    shift_i = shift_xyz(1); shift_j = shift_xyz(2); shift_k = shift_xyz(3)

    nx = tiles%Nx; ny = tiles%Ny; nz = tiles%Nz
    hx = (metric%alpha1 - metric%alpha0)/(real(nx,8)-1+2.0_8*shift_i)
    hy = (metric%beta1  - metric%beta0 )/(real(ny,8)-1+2.0_8*shift_j)
    hz = 1.0_8 / max(1.0_8,real(nz,8)-1+2.0_8*shift_k)

    ts = partition%ts
    te = partition%te

    allocate(mesh%tile(ts:te))
    mesh%ts = ts
    mesh%te = te

    mesh%scale          = metric%scale
    mesh%vertical_scale = metric%vertical_scale
    mesh%omega          = metric%omega
    mesh%rotation_axis  = metric%rotation_axis

    do t = ts, te

        call tiles%tile(t)%getind(is, ie, js, je, ks, ke)
        pind = partition%panel_map(t)

        call mesh%tile(t)%init(is, ie, js, je, ks, ke, halo_width)

        mesh%tile(t)%nx = nx
        mesh%tile(t)%ny = ny
        mesh%tile(t)%nz = nz

        mesh%tile(t)%hx = hx
        mesh%tile(t)%hy = hy
        mesh%tile(t)%hz = hz

        mesh%tile(t)%shift_i = shift_i
        mesh%tile(t)%shift_j = shift_j
        mesh%tile(t)%shift_k = shift_k

        mesh%tile(t)%alpha_0 = metric%alpha0
        mesh%tile(t)%beta_0  = metric%beta0

        mesh%tile(t)%panel_ind = pind

        if(present(orography)) then
            call calc_mesh_metrics_orog(mesh%tile(t),metric,h_top,halo_width,&
                                        orography%h%tile(t),orography%dh_alpha%tile(t),&
                                        orography%dh_beta%tile(t))
        else
            call calc_mesh_metrics_2d(mesh%tile(t),metric,h_top,halo_width)
        end if

    end do

end subroutine create_mesh

subroutine calc_mesh_metrics_2d(mesh,metric,h_top,halo_width)

    type(tile_mesh_t), intent(inout) :: mesh
    class(metric_t),   intent(in)    :: metric
    real(kind=8),      intent(in)    :: h_top
    integer(kind=4),   intent(in)    :: halo_width

    integer(kind=4) :: i, j, k, pind
    real(kind=8) :: alpha, beta, eta, vec(3)

    pind = mesh%panel_ind

    do k = mesh%ks, mesh%ke
        eta = mesh%get_eta(k)
        do j = mesh%js - halo_width, mesh%je + halo_width
            beta = mesh%get_beta(j)
            do i = mesh%is - halo_width, mesh%ie + halo_width

                alpha = mesh%get_alpha(i)
                vec(1:3) = metric%calculate_r(pind, alpha, beta)

                mesh%rx(i,j,k) = vec(1)
                mesh%ry(i,j,k) = vec(2)
                mesh%rz(i,j,k) = vec(3)
                mesh%h(i,j,k)  = metric%calculate_h(pind,alpha,beta,eta,0.0_8,h_top)

                mesh%Q (1:6,i,j,k) = metric%calculate_Q(pind, alpha, beta)
                mesh%Qi(1:6,i,j,k) = metric%calculate_Qi(pind, alpha, beta)

                mesh%J(i,j,k)= metric%calculate_J(pind, alpha, beta)

                mesh%a1(1:4,i,j,k) = metric%calculate_a1(pind,alpha,beta)
                mesh%a2(1:4,i,j,k) = metric%calculate_a2(pind,alpha,beta)
                mesh%a3(1:4,i,j,k) = metric%calculate_a3(pind,alpha,beta)
                mesh%b1(1:3,i,j,k) = metric%calculate_b1(pind,alpha,beta)
                mesh%b2(1:3,i,j,k) = metric%calculate_b2(pind,alpha,beta)
                mesh%b3(1:4,i,j,k) = metric%calculate_b3(pind,alpha,beta)
            end do
        end do
    end do

end subroutine calc_mesh_metrics_2d

subroutine calc_mesh_metrics_orog(mesh,metric,h_top,halo_width,h_surf,dhs_alpha,dhs_beta)

    type(tile_mesh_t),       intent(inout) :: mesh
    class(metric_t),         intent(in)    :: metric
    real(kind=8),            intent(in)    :: h_top
    integer(kind=4),         intent(in)    :: halo_width
    type(tile_field_t),      intent(in)    :: h_surf, dhs_alpha, dhs_beta

    integer(kind=4) :: i, j, k, pind
    real(kind=8) :: alpha, beta, eta, vec(3)

    pind = mesh%panel_ind

    do k = mesh%ks, mesh%ke
        eta = mesh%get_eta(k)
        do j = mesh%js - halo_width, mesh%je + halo_width
            beta = mesh%get_beta(j)
            do i = mesh%is - halo_width, mesh%ie + halo_width

                alpha = mesh%get_alpha(i)
                vec(1:3) = metric%calculate_r(pind, alpha, beta,eta,h_surf%p(i,j,1),h_top)

                mesh%rx(i,j,k) = vec(1)
                mesh%ry(i,j,k) = vec(2)
                mesh%rz(i,j,k) = vec(3)
                mesh%h(i,j,k)  = metric%calculate_h(pind,alpha,beta,eta,h_surf%p(i,j,1),h_top)

                mesh%Q (1:6,i,j,k) = metric%calculate_Q(pind, alpha, beta, eta, &
                                     h_surf%p(i,j,1), dhs_alpha%p(i,j,1), dhs_beta%p(i,j,1),h_top)
                mesh%Qi(1:6,i,j,k) = metric%calculate_Qi(pind, alpha, beta, eta, &
                                     h_surf%p(i,j,1), dhs_alpha%p(i,j,1), dhs_beta%p(i,j,1),h_top)

                mesh%J(i,j,k)= metric%calculate_J(pind, alpha, beta, eta, h_surf%p(i,j,1), h_top)

                mesh%a1(1:4,i,j,k) = metric%calculate_a1(pind,alpha,beta,eta,&
                                                  h_surf%p(i,j,1),dhs_alpha%p(i,j,1),h_top)
                mesh%a2(1:4,i,j,k) = metric%calculate_a2(pind,alpha,beta,eta,&
                                                  h_surf%p(i,j,1),dhs_beta%p(i,j,1),h_top)
                mesh%a3(1:4,i,j,k) = metric%calculate_a3(pind,alpha,beta,eta,h_surf%p(i,j,1),h_top)
                mesh%b1(1:3,i,j,k) = metric%calculate_b1(pind,alpha,beta,eta,h_surf%p(i,j,1),&
                                                        dhs_alpha%p(i,j,1),dhs_beta%p(i,j,1),h_top)
                mesh%b2(1:3,i,j,k) = metric%calculate_b2(pind,alpha,beta,eta,h_surf%p(i,j,1),&
                                                        dhs_alpha%p(i,j,1),dhs_beta%p(i,j,1),h_top)
                mesh%b3(1:4,i,j,k) = metric%calculate_b3(pind,alpha,beta,eta,h_surf%p(i,j,1),&
                                                        dhs_alpha%p(i,j,1),dhs_beta%p(i,j,1),h_top)
            end do
        end do
    end do

end subroutine calc_mesh_metrics_orog

end module mesh_factory_mod
