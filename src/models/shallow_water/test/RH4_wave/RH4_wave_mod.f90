module RH4_wave_mod

use domain_mod,               only : domain_t
use domain_factory_mod,       only : create_domain
use stvec_mod,                only : stvec_t
use stvec_swm_mod,            only : stvec_swm_t
use stvec_swm_factory_mod,    only : create_stvec_swm
use operator_mod,             only : operator_t
use operator_swm_factory_mod, only : create_swm_operator
use timescheme_mod,           only : timescheme_t
use timescheme_factory_mod,   only : create_timescheme
use outputer_abstract_mod,    only : outputer_t, outputer_vector_t
use outputer_factory_mod,     only : create_master_paneled_outputer,&
                                     create_latlon_outputer, create_latlon_vec_outputer
use parcomm_mod,              only : parcomm_global
use config_swm_mod,           only : config_swm_t

use const_mod,  only : Earth_grav, Earth_omega, Earth_radii, pi, Earth_sidereal_T

use test_fields_mod, only : rh4_wave_height_generator_t, rh4_wave_wind_generator_t

implicit none

real(kind=8) :: grav   = Earth_grav
real(kind=8) :: omega  = Earth_omega
real(kind=8) :: a      = Earth_radii
real(kind=8) :: h_mean = 8.0_8*10**3

! type(config_RH4_wave_t) :: config_RH4_wave

type(rh4_wave_height_generator_t) :: height_field
type(rh4_wave_wind_generator_t)   :: velocity_field


contains

subroutine run_RH4_wave()

    use const_mod,    only : pi
    use vec_math_mod, only : l2norm

    use namelist_read_mod, only : read_namelist_as_str

    type(config_swm_t) :: config
    type(domain_t)     :: domain

    class(stvec_t),           allocatable :: state, state_err, state_ex
    class(operator_t),        allocatable :: operator
    class(timescheme_t),      allocatable :: timescheme
    class(outputer_t),        allocatable :: outputer
    class(outputer_vector_t), allocatable :: outputer_vec

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: namelist_string

    ! real(kind=8),    parameter :: rotation_period = 1.0_8, rotation_axis_angle = 0.0_8*pi/8
    ! integer(kind=4), parameter :: N_periods=2

    real(kind=8), parameter    :: dt = 900._8
    real(kind=8), parameter    :: tau_write = 1800._8
    integer(kind=4), parameter :: nstep_write = nint(tau_write/dt)

    real(kind=8)    :: time, l2err, l2_ex
    integer(kind=4) :: it

    call read_namelist_as_str(namelist_string, "namelist_swm", parcomm_global%myid)

    call config%parse(namelist_string)

    config%config_domain%config_metric%omega = omega
    config%config_domain%config_metric%scale = a

    height_field = rh4_wave_height_generator_t(h_mean = h_mean, omega = omega, &
           a = a, grav = grav)

    velocity_field = rh4_wave_wind_generator_t(a = a, omega = omega)
    !
    ! if (parcomm_global%myid==0) print*, "Advective CFL = ", &
    !         real(dt*config_ts2%u0/(2*pi/4/config%config_domain%N)/config_ts2%a,4)
    !
    if (parcomm_global%myid==0) print*, "Gravity Wave CFL = ", &
            real(dt*sqrt(h_mean*grav)    &
            /(2*pi/4/config%config_domain%N)/a,4)

    call create_domain(domain, config%config_domain)

    call create_stvec_swm(state,     domain, halo_width, 0)
    call create_stvec_swm(state_ex,  domain, 0         , 0)
    call create_stvec_swm(state_err, domain, 0         , 0)

    call create_swm_operator(operator, grav, config, domain)

    call create_timescheme(timescheme, state, 'rk4')

    print*, 4*domain%partition%Nh, 2*domain%partition%Nh+1

    call create_latlon_outputer(outputer, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
    call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", &
                                   "covariant", domain)

    call get_exact_solution(state,    domain)
    call get_exact_solution(state_ex, domain)

    select type(state)
    class is (stvec_swm_t)
        call outputer%write(state%h, domain, 'h_swm.dat',1)
        call outputer_vec%write(state%u,state%v,domain,"u.dat","v.dat",1)
    end select
    ! select type(state_ex)
    ! class is (stvec_swm_t)
    !     l2_ex = l2norm(state_ex%h,  domain%mesh_p, domain%parcomm)
    ! end select
    !
    do it = 1, int(20*24*3600/dt)

        call timescheme%step(state, operator, domain, dt)

        time = it*dt

        ! call state_err%assign(1.0_8, state_ex, -1.0_8, state, domain)

        if(mod(it,nstep_write) == 0) then

            select type(state)
            class is (stvec_swm_t)
                call outputer%write(state%h, domain, 'h.dat', int(it/nstep_write))
                l2err = l2norm(state%h, domain%mesh_p, domain%parcomm)
                if (parcomm_global%myid==0) print*, "Hours = ", real(time/3600 ,4), &
                                                    "L2err =", real(l2err,4)
            end select
        end if
    end do

end subroutine run_RH4_wave

subroutine get_exact_solution(state, domain)

    use test_fields_mod, only : set_vector_test_field, &
                                set_scalar_test_field

    class(stvec_t), intent(inout) :: state
    type(domain_t), intent(in)    :: domain

    select type(state)
    class is (stvec_swm_t)

        call set_scalar_test_field(state%h, height_field, domain%mesh_p, 0)
        call set_vector_test_field(state%u, state%v, velocity_field, &
              domain%mesh_u, domain%mesh_v, 0, "covariant")

    end select


end subroutine get_exact_solution

end module RH4_wave_mod
