module nonlin_nh_oper_mod

use operator_mod,                    only : operator_t
use grid_field_mod,                  only : grid_field_t
use abstract_grad_3d_mod,            only : grad_3d_operator_t
use abstract_div_3d_mod,             only : div_3d_operator_t
! use abstract_co2contra_mod,          only : co2contra_operator_t
use abstract_mixvec_transform_mod,   only : mixvec_transform_t
use abstract_interpolators3d_mod,    only : interpolator_w2uv_t
use abstract_scalar_advection3d_mod, only : scalar_advection3d_t
use abstract_vector_advection3d_mod, only : vector_advection3d_t
use abstract_coriolis_mod,           only : coriolis_operator_t
use stvec_mod,                       only : stvec_t
use stvec_nh_mod,                    only : stvec_nh_t
use parcomm_mod,                     only : parcomm_global
use domain_mod,                      only : domain_t

implicit none

type, extends(operator_t) :: nonlin_nh_operator_t

    class(grad_3d_operator_t),   allocatable :: grad_op
    class(div_3d_operator_t),    allocatable :: div_op
    ! class(co2contra_operator_t), allocatable :: co2contra_op
    class(mixvec_transform_t),   allocatable :: mixvec_transform
    class(interpolator_w2uv_t),  allocatable :: theta2uv_op
    class(scalar_advection3d_t), allocatable :: p_adv_oper, theta_adv_oper
    class(vector_advection3d_t), allocatable :: momentum_adv_op
    class(coriolis_operator_t),  allocatable :: coriolis_op

    type(grid_field_t) :: theta_u, theta_v
    type(grid_field_t) :: grad_x, grad_y, grad_z
    type(grid_field_t) :: grad_x_contra, grad_y_contra, grad_z_real
    type(grid_field_t) :: div3, eta_dot
    type(grid_field_t) :: Phi !geopotential field
    type(grid_field_t) :: grad_Phi_x, grad_Phi_y, grad_Phi_z ! geopotential grad covariant components

    contains
    procedure, public :: apply

end type nonlin_nh_operator_t

contains

subroutine apply(this, vout, vin, domain)

    use const_mod, only : rgaz, Cv, Cp

    class(nonlin_nh_operator_t), intent(inout) :: this
    class(stvec_t),                     intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),                     intent(inout) :: vin
    type(domain_t),                     intent(in)    :: domain

    select type(vout)
    class is (stvec_nh_t)
    select type(vin)
    class is (stvec_nh_t)

        call this%mixvec_transform%transform_mix2contra(this%eta_dot,vin%u,vin%v,vin%w,domain)

        call this%momentum_adv_op%calc_vec_adv3d(vout%u,vout%v,vout%w, &
                                                vin%u,vin%v,vin%w,this%eta_dot,domain)

        call this%grad_op%calc_grad(this%grad_x, this%grad_y, this%grad_z,vin%P,domain)

        call this%theta2uv_op%interp_w2uv(this%theta_u,this%theta_v,vin%theta,domain)
        call this%grad_x%assign_prod(-Cp,this%grad_x,this%theta_u,domain%mesh_u)
        call this%grad_y%assign_prod(-Cp,this%grad_y,this%theta_v,domain%mesh_v)
        call this%grad_z%assign_prod(-Cp,this%grad_z,vin%theta,domain%mesh_w)
        call this%grad_x%update(-1.0_8,this%grad_phi_x,domain%mesh_u)
        call this%grad_y%update(-1.0_8,this%grad_phi_y,domain%mesh_v)
        call this%grad_z%update(-1.0_8,this%grad_phi_z,domain%mesh_w)

        call this%mixvec_transform%transform_co2mix(this%grad_x_contra, this%grad_y_contra,   &
                                                    this%grad_z_real,this%grad_x,this%grad_y, &
                                                    this%grad_z,domain)
        call vout%u%update(1.0_8,this%grad_x_contra,domain%mesh_u)
        call vout%v%update(1.0_8,this%grad_y_contra,domain%mesh_v)
        call vout%w%update(1.0_8,this%grad_z_real,  domain%mesh_w)

        call add_top_damping(vout%w,vin%w,domain%mesh_w)

        call this%coriolis_op%calc_coriolis_contra(this%grad_x, this%grad_y, &
                                                           vin%u, vin%v, domain)
        call vout%u%update(1.0_8,this%grad_x,domain%mesh_u)
        call vout%v%update(1.0_8,this%grad_y,domain%mesh_v)

        call this%theta_adv_oper%calc_adv3d(vout%theta,vin%theta,vin%u,vin%v,this%eta_dot,domain)
        call this%p_adv_oper%calc_adv3d(vout%P,vin%P,vin%u,vin%v,this%eta_dot,domain)

        call this%div_op%calc_div(this%div3,vin%u,vin%v,this%eta_dot,domain)
        call this%div3%assign_prod(-rgaz/Cv,vin%P,this%div3,domain%mesh_p)
        call vout%P%update(1.0_8,this%div3,domain%mesh_p)

        vout%model_time = 1.0_8 != dt / dt

        ! call apply_0boundary_conds(vout%w,domain%mesh_w)
    class default
        call parcomm_global%abort("Ptheta_linear_nh_oper_t, vin type error")
    end select
    class default
        call parcomm_global%abort("Ptheta_linear_nh_oper_t, vout type error")
    end select


end subroutine apply

subroutine apply_0boundary_conds(eta_dot, mesh)

    use mesh_mod, only : mesh_t

    type(grid_field_t), intent(inout) :: eta_dot
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t, i, j, ks, ke, is,ie, js, je

    do t = mesh%ts, mesh%te
        is = mesh%tile(t)%is; ie = mesh%tile(t)%ie
        js = mesh%tile(t)%js; je = mesh%tile(t)%je
        ks = mesh%tile(t)%ks; ke = mesh%tile(t)%ke

        do j = js, je
            do i=is,ie
                eta_dot%tile(t)%p(i,j,ks) = 0.0_8
                eta_dot%tile(t)%p(i,j,ke) = 0.0_8
            end do
        end do
    end do
end subroutine apply_0boundary_conds

subroutine add_top_damping(w1,w0,mesh)

    use mesh_mod, only : mesh_t

    type(grid_field_t), intent(inout) :: w1
    type(grid_field_t), intent(in)    :: w0
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: i, j, k, t
    real(kind=8)    :: sigma, sigma0 = 1.0_8/25.0_8
    real(kind=8)    :: h0 = 20e3_8, dh = 10e3_8

    do t=mesh%ts, mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie
                    sigma = sigma0*min(1.0_8,max(0.0_8,(mesh%tile(t)%h(i,j,k)-h0)/dh))
                    w1%tile(t)%p(i,j,k) = w1%tile(t)%p(i,j,k)-sigma*w0%tile(t)%p(i,j,k)
                end do
            end do
        end do
    end do
end subroutine add_top_damping

end module nonlin_nh_oper_mod
