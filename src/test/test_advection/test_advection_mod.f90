module test_advection_mod

use domain_mod,                     only : domain_t
use domain_factory_mod,             only : create_domain
use stvec_mod,                      only : stvec_t
use stvec_advection_mod,            only : stvec_advection_t
use stvec_advection_factory_mod,    only : create_stvec_advection
use operator_mod,                   only : operator_t
use operator_advection_factory_mod, only : create_advection_operator
use timescheme_mod,                 only : timescheme_t
use timescheme_factory_mod,         only : create_timescheme
use outputer_abstract_mod,          only : outputer_t
use outputer_factory_mod,           only : create_master_paneled_outputer

implicit none

contains

subroutine test_advection()

    type(domain_t) :: domain
    integer(kind=4), parameter  :: N = 32, nz = 1, halo_width = 5

    real(kind=8) :: dt = 0.05_8

    class(stvec_t),      allocatable :: state
    class(operator_t),   allocatable :: operator
    class(timescheme_t), allocatable :: timescheme
    class(outputer_t),   allocatable :: outputer

    integer(kind=4) :: it

    call create_domain(domain, "cube", "A", N, nz)

    call create_stvec_advection(state, domain, halo_width, 0)

    call create_advection_operator(operator, "divergence_a2_ecs", domain)

    call create_timescheme(timescheme, state, 'rk4')

    call create_master_paneled_outputer(outputer, "p", domain)

    call init_initial_conditions(state, domain)

    do it = 1, 100

        call timescheme%step(state, operator, domain, dt)

        select type(state)
        class is (stvec_advection_t)
            call outputer%write(state%h, domain, 'h_adv.dat', it)
        end select

    end do


end subroutine test_advection

subroutine init_initial_conditions(state, domain)

    use test_fields_mod, only : set_vector_test_field, &
                                vector_field => solid_rotation_field_generator, &
                                set_scalar_test_field, &
                                scalar_field=>gaussian_hill_scalar_field_generator

    class(stvec_t), intent(inout) :: state
    type(domain_t), intent(in)    :: domain

    select type(state)
    class is (stvec_advection_t)

        call set_scalar_test_field(state%h, scalar_field, domain%mesh_p, 0)
        call set_vector_test_field(state%u, state%v, vector_field, &
              domain%mesh_u, domain%mesh_v, 0, "contravariant")

    end select


end subroutine init_initial_conditions


end module test_advection_mod
