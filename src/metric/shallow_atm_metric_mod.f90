module shallow_atm_metric_mod

use metric_mod,                      only : metric_t
use abstract_vertical_transform_mod, only : vertical_transform_t
use mesh_mod,                        only : mesh_t
use orography_mod,                   only : orography_1mesh_t

implicit none

type, extends(metric_t) :: shallow_atm_metric_t
    class(metric_t), allocatable             :: metric_2d
    class(vertical_transform_t), allocatable :: vertical_transform

    contains

    procedure :: calculate_r_orog
    procedure :: calculate_r_2d
    procedure :: calculate_h
    procedure :: calculate_a1_orog
    procedure :: calculate_a1_2d
    procedure :: calculate_a2_orog
    procedure :: calculate_a2_2d
    procedure :: calculate_a3_orog
    procedure :: calculate_a3_2d
    procedure :: calculate_b1_orog
    procedure :: calculate_b1_2d
    procedure :: calculate_b2_orog
    procedure :: calculate_b2_2d
    procedure :: calculate_b3_orog
    procedure :: calculate_b3_2d
    procedure :: calculate_Q_orog
    procedure :: calculate_Q_2d
    procedure :: calculate_Qi_orog
    procedure :: calculate_Qi_2d
    procedure :: calculate_J_orog
    procedure :: calculate_J_2d
    procedure :: calculate_G_orog
    procedure :: calculate_G_2d

    procedure :: transform_cartesian_to_native

    procedure :: set_curvilinear_mesh

end type shallow_atm_metric_t

contains

pure function calculate_r_orog(this, panel_ind, alpha, beta, eta, h_surf, h_top) result(r)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, h_top
    real(kind=8)                    :: r(3)

    r = this%metric_2d%calculate_r(panel_ind,alpha,beta)
end function calculate_r_orog

pure function calculate_r_2d(this, panel_ind, alpha, beta) result(r)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: r(3)

    r = this%metric_2d%calculate_r(panel_ind,alpha,beta)
end function calculate_r_2d
pure function calculate_h(this,panel_ind,alpha,beta,eta,h_surf,h_top) result(h)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),             intent(in) :: panel_ind
    real(kind=8),                intent(in) :: alpha, beta, eta
    real(kind=8),                intent(in) :: h_surf, h_top
    real(kind=8)                            :: h

    h = this%vertical_transform%calc_z(h_surf,h_top,eta)
end function calculate_h
pure function calculate_a1_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dcov_h_surf, h_top) result(a)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dcov_h_surf, h_top
    real(kind=8)                    :: a(4)

    a = this%metric_2d%calculate_a1(panel_ind,alpha,beta)
    a(4) = this%vertical_transform%calc_dz_dh_surf(eta)*dcov_h_surf

end function calculate_a1_orog

pure function calculate_a1_2d(this, panel_ind, alpha, beta) result(a)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: a(4)

    a = this%metric_2d%calculate_a1(panel_ind,alpha,beta)

end function calculate_a1_2d

pure function calculate_a2_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dcov_h_surf, h_top) result(a)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dcov_h_surf, h_top
    real(kind=8)                    :: a(4)

    a = this%metric_2d%calculate_a2(panel_ind,alpha,beta)
    a(4) = this%vertical_transform%calc_dz_dh_surf(eta)*dcov_h_surf

end function calculate_a2_orog

pure function calculate_a2_2d(this, panel_ind, alpha, beta) result(a)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: a(4)

    a = this%metric_2d%calculate_a2(panel_ind,alpha,beta)

end function calculate_a2_2d

pure function calculate_a3_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, h_top) result(a)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, h_top
    real(kind=8)                    :: a(4)

    a(1:3) = [0.0_8, 0.0_8, 0.0_8]
    a(4) = this%vertical_transform%calc_dz_deta(h_surf,h_top,eta) / this%vertical_scale

end function calculate_a3_orog

pure function calculate_a3_2d(this, panel_ind, alpha, beta) result(a)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: a(4)

    a = [0.0_8, 0.0_8, 0.0_8, 1.0_8]

end function calculate_a3_2d

pure function calculate_b1_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(b)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: b(3)

    b = this%metric_2d%calculate_b1(panel_ind,alpha,beta)

end function calculate_b1_orog

pure function calculate_b1_2d(this, panel_ind, alpha, beta) result(b)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: b(3)

    b = this%metric_2d%calculate_b1(panel_ind,alpha,beta)
end function calculate_b1_2d

pure function calculate_b2_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(b)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: b(3)

    b = this%metric_2d%calculate_b2(panel_ind,alpha,beta)

end function calculate_b2_orog

pure function calculate_b2_2d(this, panel_ind, alpha, beta) result(b)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: b(3)

    b = this%metric_2d%calculate_b2(panel_ind,alpha,beta)
end function calculate_b2_2d

pure function calculate_b3_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(b)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: b(4)

    real(kind=8) :: b1(3), b2(3), dh_dalpha, dh_dbeta, dh_deta, dh_dhs

    b1 = this%metric_2d%calculate_b1(panel_ind,alpha,beta)
    b2 = this%metric_2d%calculate_b2(panel_ind,alpha,beta)

    dh_dhs = this%vertical_transform%calc_dz_dh_surf(eta)
    dh_dalpha  = dh_dhs*dh_surf_dalpha
    dh_dbeta   = dh_dhs*dh_surf_dbeta
    dh_deta = this%vertical_transform%calc_dz_deta(h_surf,h_top,eta) / this%vertical_scale

    b(1) =-(dh_dalpha*b1(1)+dh_dbeta*b2(1)) / dh_deta
    b(2) =-(dh_dalpha*b1(2)+dh_dbeta*b2(2)) / dh_deta
    b(3) =-(dh_dalpha*b1(3)+dh_dbeta*b2(3)) / dh_deta
    b(4) = 1.0_8 / dh_deta

end function calculate_b3_orog

pure function calculate_b3_2d(this, panel_ind, alpha, beta) result(b)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: b(4)

    b = [0.0_8, 0.0_8, 0.0_8, 1.0_8]
end function calculate_b3_2d

pure function calculate_Q_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(Q)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: Q(6)

    real(kind=8) :: dh_dhs, dh_deta, dh_dalpha, dh_dbeta

    Q = this%metric_2d%calculate_Q(panel_ind,alpha,beta)

    dh_dhs = this%vertical_transform%calc_dz_dh_surf(eta)
    dh_dalpha  = dh_dhs*dh_surf_dalpha
    dh_dbeta   = dh_dhs*dh_surf_dbeta
    dh_deta = this%vertical_transform%calc_dz_deta(h_surf,h_top,eta) / this%vertical_scale

    Q(1) = Q(1) + dh_dalpha**2
    Q(2) = Q(2) + dh_dalpha*dh_dbeta
    Q(3) = Q(3) + dh_dbeta**2
    Q(4) = dh_dalpha*dh_deta
    Q(5) = dh_dbeta*dh_deta
    Q(6) = dh_deta**2

end function calculate_Q_orog

pure function calculate_Q_2d(this, panel_ind, alpha, beta) result(Q)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: Q(6)

    Q = this%metric_2d%calculate_Q(panel_ind,alpha,beta)

end function calculate_Q_2d

pure function calculate_Qi_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(Qi)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: Qi(6)

    real(kind=8) :: dh_dhs, dh_deta, dh_dalpha, dh_dbeta

    dh_dhs = this%vertical_transform%calc_dz_dh_surf(eta)
    dh_dalpha  = dh_dhs*dh_surf_dalpha
    dh_dbeta   = dh_dhs*dh_surf_dbeta
    dh_deta = this%vertical_transform%calc_dz_deta(h_surf,h_top,eta) / this%vertical_scale

    Qi = this%metric_2d%calculate_Qi(panel_ind,alpha,beta)
    Qi(4) = -(dh_dalpha*Qi(1)+dh_dbeta*Qi(2)) / dh_deta
    Qi(5) = -(dh_dalpha*Qi(2)+dh_dbeta*Qi(3)) / dh_deta
    Qi(6) = (1+dh_dalpha*(dh_dalpha*Qi(1)+dh_dbeta*Qi(2))+&
               dh_dbeta *(dh_dalpha*Qi(2)+dh_dbeta*Qi(3))) / dh_deta**2

end function calculate_Qi_orog

pure function calculate_Qi_2d(this, panel_ind, alpha, beta) result(Qi)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: Qi(6)

    Qi = this%metric_2d%calculate_Qi(panel_ind,alpha,beta)

end function calculate_Qi_2d

pure function calculate_J_orog(this, panel_ind, alpha, beta, eta, h_surf, h_top) result(J)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, h_top
    real(kind=8)                    :: J

    J = this%metric_2d%calculate_J(panel_ind,alpha,beta)*&
        this%vertical_transform%calc_dz_deta(h_surf,h_top,eta) / this%vertical_scale
end function calculate_J_orog

pure function calculate_J_2d(this, panel_ind, alpha, beta) result(J)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: J

    J = this%metric_2d%calculate_J(panel_ind,alpha,beta)

end function calculate_J_2d

pure function calculate_G_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(G)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: G(3,3,3)

    G = this%metric_2d%calculate_G(panel_ind,alpha,beta)

end function calculate_G_orog

pure function calculate_G_2d(this, panel_ind, alpha, beta) result(G)
    class(shallow_atm_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: G(3,3,3)

    G = this%metric_2d%calculate_G(panel_ind,alpha,beta)

end function calculate_G_2d

subroutine transform_cartesian_to_native(this,panel_ind, alpha, beta, r)
    class(shallow_atm_metric_t), intent(in)  :: this
    integer(kind=4),     intent(out) :: panel_ind
    real(kind=8),        intent(out) :: alpha, beta
    real(kind=8),        intent(in)  :: r(3)

    call this%metric_2d%transform_cartesian_to_native(panel_ind,alpha,beta,r)

end subroutine transform_cartesian_to_native

subroutine set_curvilinear_mesh(this,mesh,h_top,orog)
    class(shallow_atm_metric_t), intent(in)    :: this
    type(mesh_t),                intent(inout) :: mesh
    real(kind=8),                intent(in)    :: h_top
    type(orography_1mesh_t),     intent(in), optional :: orog

    if(present(orog)) then
        call set_curvilinear_mesh_orog(this,mesh,h_top,orog)
    else
        call set_curvilinear_mesh_norog(this,mesh,h_top)
    end if

end subroutine

subroutine set_curvilinear_mesh_orog(this,mesh,h_top,orog)
    class(shallow_atm_metric_t), intent(in)    :: this
    type(mesh_t),                intent(inout) :: mesh
    real(kind=8),                intent(in)    :: h_top
    type(orography_1mesh_t),     intent(in)    :: orog

    integer(kind=4) :: t, i, j, k, halo_width
    real(kind=8)    :: eta, dh_dhs, dh_dalpha, dh_dbeta, dh_deta, h_surf
    real(kind=8)    :: qi(3), b1(3), b2(3)

    call this%metric_2d%set_curvilinear_mesh(mesh,h_top)

    do t = mesh%ts, mesh%te
        halo_width = mesh%tile(t)%halo_width
        do k=mesh%tile(t)%ks, mesh%tile(t)%ke
            eta = mesh%tile(t)%get_eta(k)
            dh_dhs = this%vertical_transform%calc_dz_dh_surf(eta)
            do j = mesh%tile(t)%js-halo_width, mesh%tile(t)%je+halo_width
                do i = mesh%tile(t)%is-halo_width, mesh%tile(t)%ie+halo_width
                    dh_dalpha  = dh_dhs*orog%dh_alpha%tile(t)%p(i,j,1)
                    dh_dbeta   = dh_dhs*orog%dh_beta%tile(t)%p(i,j,1)
                    h_surf     = orog%h%tile(t)%p(i,j,1)
                    dh_deta = this%vertical_transform%calc_dz_deta(h_surf,h_top,eta) &
                                                                / this%vertical_scale

                    mesh%tile(t)%h(i,j,k) = this%vertical_transform%calc_z(h_surf,h_top,eta)

                    mesh%tile(t)%a1(4,i,j,k) = dh_dalpha
                    mesh%tile(t)%a2(4,i,j,k) = dh_dbeta
                    mesh%tile(t)%a3(1:4,i,j,k) = [0.0_8, 0.0_8, 0.0_8, dh_deta]
                    !b1 and b2 remains the same as initialized by metric_2d
                    b1(1:3) = mesh%tile(t)%b1(1:3,i,j,k)
                    b2(1:3) = mesh%tile(t)%b2(1:3,i,j,k)

                    mesh%tile(t)%b3(1,i,j,k) =-(dh_dalpha*b1(1)+dh_dbeta*b2(1)) / dh_deta
                    mesh%tile(t)%b3(2,i,j,k) =-(dh_dalpha*b1(2)+dh_dbeta*b2(2)) / dh_deta
                    mesh%tile(t)%b3(3,i,j,k) =-(dh_dalpha*b1(3)+dh_dbeta*b2(3)) / dh_deta
                    mesh%tile(t)%b3(4,i,j,k) = 1.0_8 / dh_deta

                    mesh%tile(t)%Q(1,i,j,k) = mesh%tile(t)%Q(1,i,j,k) + dh_dalpha**2
                    mesh%tile(t)%Q(2,i,j,k) = mesh%tile(t)%Q(2,i,j,k) + dh_dalpha*dh_dbeta
                    mesh%tile(t)%Q(3,i,j,k) = mesh%tile(t)%Q(3,i,j,k) + dh_dbeta**2
                    mesh%tile(t)%Q(4,i,j,k) = dh_dalpha*dh_deta
                    mesh%tile(t)%Q(5,i,j,k) = dh_dbeta*dh_deta
                    mesh%tile(t)%Q(6,i,j,k) = dh_deta**2

                    qi(1:3) = mesh%tile(t)%Qi(1:3,i,j,k)
                    mesh%tile(t)%Qi(4,i,j,k) = -(dh_dalpha*qi(1)+dh_dbeta*qi(2)) / dh_deta
                    mesh%tile(t)%Qi(5,i,j,k) = -(dh_dalpha*qi(2)+dh_dbeta*qi(3)) / dh_deta
                    mesh%tile(t)%Qi(6,i,j,k) = (1+dh_dalpha*(dh_dalpha*qi(1)+dh_dbeta*qi(2))+&
                                                dh_dbeta *(dh_dalpha*qi(2)+dh_dbeta*qi(3))) / dh_deta**2

                    mesh%tile(t)%J(i,j,k) = mesh%tile(t)%J(i,j,k)*dh_deta
                end do
            end do
        end do
    end do

end subroutine

subroutine set_curvilinear_mesh_norog(this,mesh,h_top)
    class(shallow_atm_metric_t), intent(in)    :: this
    type(mesh_t),                intent(inout) :: mesh
    real(kind=8),                intent(in)    :: h_top

    integer(kind=4) :: t, i, j, k, halo_width
    real(kind=8)    :: alpha, beta, eta

    call this%metric_2d%set_curvilinear_mesh(mesh,h_top)

    do t = mesh%ts, mesh%te
        halo_width = mesh%tile(t)%halo_width
        do k=mesh%tile(t)%ks, mesh%tile(t)%ke
            eta = mesh%tile(t)%get_eta(k)
            do j = mesh%tile(t)%js-halo_width, mesh%tile(t)%je+halo_width
                do i = mesh%tile(t)%is-halo_width, mesh%tile(t)%ie+halo_width
                    mesh%tile(t)%h(i,j,k) = this%vertical_transform%calc_z(0.0_8,h_top,eta)
                    mesh%tile(t)%Q(4,i,j,k) = 0.0_8
                    mesh%tile(t)%Q(5,i,j,k) = 0.0_8
                    mesh%tile(t)%Q(6,i,j,k) = 1._8
                    mesh%tile(t)%Qi(4,i,j,k) = 0.0_8
                    mesh%tile(t)%Qi(5,i,j,k) = 0.0_8
                    mesh%tile(t)%Qi(6,i,j,k) = 1._8
                    mesh%tile(t)%J(i,j,k) = mesh%tile(t)%J(i,j,k)*1.0_8
                end do
            end do
        end do
    end do
end subroutine

end module shallow_atm_metric_mod
