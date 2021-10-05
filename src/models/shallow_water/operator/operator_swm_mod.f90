module operator_swm_mod

use stvec_mod,      only : stvec_t
use domain_mod,     only : domain_t
use operator_mod,   only : operator_t
use grid_field_mod, only : grid_field_t

use abstract_div_mod,       only : div_operator_t
use abstract_grad_mod,      only : grad_operator_t
use abstract_coriolis_mod,  only : coriolis_operator_t
use abstract_curl_mod,      only : curl_operator_t
use abstract_KE_mod,        only : KE_operator_t
use abstract_massflux_mod,  only : massflux_operator_t
use abstract_co2contra_mod, only : co2contra_operator_t

use stvec_swm_mod, only : stvec_swm_t
use parcomm_mod,   only : parcomm_global

use vec_math_mod, only : mass, l2norm

implicit none

type, public, extends(operator_t) :: operator_swm_t

    class(div_operator_t),       allocatable :: div_op
    class(grad_operator_t),      allocatable :: grad_op
    class(coriolis_operator_t),  allocatable :: coriolis_op
    class(curl_operator_t),      allocatable :: curl_op
    class(KE_operator_t),        allocatable :: KE_op
    class(massflux_operator_t),  allocatable :: massflux_op
    class(co2contra_operator_t), allocatable :: co2contra_op

    type(grid_field_t) :: h_surf !orography
    type(grid_field_t) :: div, grad_x, grad_y, curl
    type(grid_field_t) :: cor_u, cor_v
    type(grid_field_t) :: KE !kinetic energy
    type(grid_field_t) :: KE_diag_u, KE_diag_v !kinetic energy
    type(grid_field_t) :: ut, vt !contravariant components
    type(grid_field_t) :: hu, hv !mass fluxes in continuty eq.

!gravity acceleration. Should be moved to mesh?
    real(kind=8) :: grav

contains
    procedure, public :: apply
end type operator_swm_t

contains

subroutine apply(this, vout, vin, domain)
    class(operator_swm_t), intent(inout) :: this
    class(stvec_t),        intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),        intent(inout) :: vin
    type(domain_t),        intent(in)    :: domain

    real(kind=8) :: ke_u, ke_v, pe

    select type (vout)
    class is (stvec_swm_t)
        select type (vin)
        class is (stvec_swm_t)

            call this%co2contra_op%transform(this%ut, this%vt, vin%u, vin%v, domain)
            call this%massflux_op%calc_massflux(this%hu, this%hv, &
                         vin%h, this%ut, this%vt, domain)

            !!ENERGY DIAGNOSTICS
            ! call this%KE_diag_u%assign_prod(1.0_8, this%hu, vin%u, domain%mesh_u)
            ! call this%KE_diag_v%assign_prod(1.0_8, this%hv, vin%v, domain%mesh_v)
            !
            ! ke_u = mass(this%KE_diag_u, domain%mesh_u, domain%parcomm)
            ! ke_v = mass(this%KE_diag_v, domain%mesh_v, domain%parcomm)
            !
            ! pe = l2norm(vin%h, domain%mesh_p, domain%parcomm)
            !
            ! print*, ke_u+ke_v, pe

            !momentum eq part
            call this%KE_op%calc_KE(this%KE, vin%u, vin%v, this%ut, this%vt, domain)

            !store grav*grad(h+hs+KE) in KE array
            call this%KE%update(this%grav, vin%h, this%grav, this%h_surf, domain%mesh_p)

            call this%grad_op%calc_grad(this%grad_x, this%grad_y, this%KE, domain)

            call this%curl_op%calc_curl(this%curl, vin%u, vin%v, domain)

            call this%coriolis_op%calc_coriolis_vec_inv(this%cor_u, this%cor_v, &
                                            this%hu, this%hv, vin%h, this%curl, domain)

            call vout%u%assign(-1.0_8, this%grad_x, 1.0_8, this%cor_u, domain%mesh_u)
            call vout%v%assign(-1.0_8, this%grad_y, 1.0_8, this%cor_v, domain%mesh_v)

            !continuty eq part
            call this%div_op%calc_div(this%div, this%hu, this%hv, domain)

            call vout%h%assign(-1.0_8, this%div, domain%mesh_p)
        class default
            call parcomm_global%abort("swm operator failure: vin of wrong type")
        end select
    class default
        call parcomm_global%abort("swm operator failure: vout of wrong type")
    end select
end subroutine apply

end module operator_swm_mod
