module test_swm_mod

use domain_mod,               only : domain_t
use domain_factory_mod,       only : create_domain
use stvec_mod,                only : stvec_t
use stvec_swm_mod,            only : stvec_swm_t
use stvec_swm_factory_mod,    only : create_stvec_swm
use operator_mod,             only : operator_t
use operator_swm_factory_mod, only : create_swm_operator
use timescheme_mod,           only : timescheme_t
use timescheme_factory_mod,   only : create_timescheme
use outputer_abstract_mod,    only : outputer_t
use outputer_factory_mod,     only : create_master_paneled_outputer,&
                                     create_latlon_outputer
use parcomm_mod,              only : parcomm_global
use config_swm_mod,           only : config_swm_t

use test_fields_mod, only : solid_rotation_t, solid_rotated_scalar_field_t

implicit none

type(solid_rotation_t),             allocatable :: velocity_field
type(solid_rotated_scalar_field_t), allocatable :: advected_field

contains

subroutine test_swm()

    use const_mod,    only : pi
    use vec_math_mod, only : l2norm

    use namelist_read_mod, only : read_namelist_as_str

    type(config_swm_t) :: config
    type(domain_t) :: domain

    class(stvec_t),      allocatable :: state, state_err
    class(operator_t),   allocatable :: operator
    class(timescheme_t), allocatable :: timescheme
    class(outputer_t),   allocatable :: outputer

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: namelist_string

    real(kind=8),    parameter :: rotation_period = 1.0_8, rotation_axis_angle = pi/8
    integer(kind=4), parameter :: N_periods=2

    real(kind=8), parameter    :: dt = 0.01_8
    real(kind=8), parameter    :: tau_write = 0.01_8
    integer(kind=4), parameter :: nstep_write = nint(tau_write/dt)

    real(kind=8)    :: time, l2err
    integer(kind=4) :: it

    call read_namelist_as_str(namelist_string, "namelist_swm", parcomm_global%myid)

    call config%parse(namelist_string)

    call create_domain(domain, config%config_domain)

    call create_stvec_swm(state,     domain, halo_width, 0)
    call create_stvec_swm(state_err, domain, 0         , 0)

    call create_swm_operator(operator, config, domain)

    call create_timescheme(timescheme, state, 'rk4')

    do it = 1, nint(N_periods*rotation_period/dt)

        call timescheme%step(state, operator, domain, dt)

        time = it*dt

        ! call get_exact_solution(state_err, domain, time)
        ! call state_err%update(-1.0_8, state, domain)

        ! if(mod(it,nstep_write) == 0) then
        !
        !     select type(state)
        !     class is (stvec_advection_t)
        !         call outputer%write(state%h, domain, 'h_adv.dat', int(it/nstep_write))
        !     end select
        !     select type(state_err)
        !     class is (stvec_advection_t)
        !         call outputer%write(state_err%h, domain, 'h_adv_err.dat', int(it/nstep_write))
        !         l2err = l2norm(state_err%h, domain%mesh_p, domain%parcomm)
        !         if (parcomm_global%myid==0) print*, "Periods = ", real(time ,4), &
        !                                             "L2err =", real(l2err,4)
        !     end select
        ! end if

    end do

end subroutine test_swm

end module test_swm_mod
