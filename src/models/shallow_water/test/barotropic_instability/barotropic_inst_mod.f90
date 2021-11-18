module barotropic_inst_mod

use domain_mod,                 only : domain_t
use domain_factory_mod,         only : create_domain
use stvec_mod,                  only : stvec_t
use stvec_swm_mod,              only : stvec_swm_t
use stvec_swm_factory_mod,      only : create_stvec_swm
use operator_mod,               only : operator_t
use operator_swm_factory_mod,   only : create_swm_operator
use timescheme_mod,             only : timescheme_t
use timescheme_factory_mod,     only : create_timescheme
use outputer_abstract_mod,      only : outputer_t, outputer_vector_t
use outputer_factory_mod,       only : create_master_paneled_outputer,&
                                     create_latlon_outputer, create_latlon_vec_outputer
use parcomm_mod,                only : parcomm_global
use config_swm_mod,             only : config_swm_t
use config_barotropic_inst_mod, only : config_barotropic_inst_t

use operator_swm_mod,           only : operator_swm_t

use operator_swm_diff_mod,         only : operator_swm_diff_t
use operator_swm_diff_factory_mod, only : create_swm_diff_operator

use const_mod,  only : Earth_grav, Earth_omega, Earth_radii, pi, Earth_sidereal_T

use test_fields_mod, only : barotropic_instability_height_generator_t, &
                            barotropic_instability_wind_generator_t
use key_value_mod, only : key_value_r8_t

use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field

implicit none

! type(config_RH4_wave_t) :: config_RH4_wave

type(barotropic_instability_height_generator_t) :: height_field
type(barotropic_instability_wind_generator_t)   :: velocity_field


contains

subroutine run_barotropic_inst()

    use const_mod,    only : pi
    use vec_math_mod, only : l2norm

    use namelist_read_mod, only : read_namelist_as_str

    type(config_swm_t) :: config
    type(domain_t)     :: domain

    class(stvec_t),           allocatable :: state, state_err, state_ex
    class(operator_t),        allocatable :: operator
    class(timescheme_t),      allocatable :: timescheme
    class(outputer_t),        allocatable :: outputer, outputer_curl
    class(outputer_vector_t), allocatable :: outputer_vec

    type(config_barotropic_inst_t) :: test_config

    type(grid_field_t) :: curl

    type(operator_swm_diff_t),   allocatable :: operator_diff
    class(timescheme_t),         allocatable :: timescheme_diff

    integer(kind=4) :: halo_width = 8

    character(:), allocatable :: namelist_string
    type(key_value_r8_t) :: diagnostics

    real(kind=8)      :: dt
    real(kind=8)      :: tau_write
    integer(kind=4)   :: nstep_write, nstep_diagnostics

    real(kind=8)    :: time, l2err, l2_ex
    integer(kind=4) :: it

    call read_namelist_as_str(namelist_string, "namelist_swm", parcomm_global%myid)

    call config%parse(namelist_string)

    test_config%Nq = 10*config%config_domain%N

    config%config_domain%config_metric%omega = test_config%omega
    config%config_domain%config_metric%scale = test_config%a

    dt = config%dt
    tau_write = config%tau_write
    nstep_write = nint(tau_write/dt)
    nstep_diagnostics = nint(config%tau_diagnostics/dt)

    height_field = create_barotropic_instability_height_field_generator(H0 = test_config%H0, &
                              omega = test_config%omega, grav = test_config%grav, &
                              h_pert = test_config%h_pert, u0 = test_config%u0,   &
                              a = test_config%a, Nq = test_config%Nq)

    velocity_field = barotropic_instability_wind_generator_t(u0=test_config%u0)

    if (parcomm_global%myid==0) print*, "Gravity Wave CFL = ", &
            real(dt*sqrt(test_config%H0*test_config%grav)    &
            /(2*pi/4/config%config_domain%N)/test_config%a,4)

    call create_domain(domain, config%config_domain)

    call create_stvec_swm(state,     domain, halo_width, 0)
    call create_stvec_swm(state_ex,  domain, 0         , 0)
    call create_stvec_swm(state_err, domain, 0         , 0)

    call create_swm_operator(operator, test_config%grav, config, domain)
    call create_timescheme(timescheme, state, 'rk4')

    call create_swm_diff_operator(operator_diff, config, domain)
    call create_timescheme(timescheme_diff, state, config%diff_time_scheme)

    if(config%config_domain%staggering_type == "Ah") then
        call create_latlon_outputer(outputer,          2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
        call create_latlon_outputer(outputer_curl,     2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", &
                                   config%v_components_type, domain)
    else if(config%config_domain%staggering_type == "C") then
        call create_latlon_outputer(outputer,      2*domain%partition%Nh+1, 4*domain%partition%Nh, "A", domain)
        call create_latlon_outputer(outputer_curl, 2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "C", &
                                       config%v_components_type, domain)
    else
        call parcomm_global%abort("This staggering is not implemented in"//&
                                  " barotropic instability swe test output:"//&
                                  config%config_domain%staggering_type)
    end if

    call get_exact_solution(state,    domain, config%v_components_type)
    call get_exact_solution(state_ex, domain, config%v_components_type)

    select type(state)
    class is (stvec_swm_t)
        call outputer%write(state%h, domain, 'h.dat',1)
        call outputer_vec%write(state%u,state%v,domain,"u.dat","v.dat",1)
    end select

    call create_grid_field(curl, halo_width, 0, domain%mesh_q)


    do it = 1, int(config%simulation_time/dt)

        !print *, "tstep", it
        call timescheme%step(state, operator, domain, dt)
        call timescheme_diff%step(state, operator_diff, domain, dt)

        time = it*dt

        if(mod(it, nstep_diagnostics) == 0) then
            diagnostics = operator%get_diagnostics(state, domain)
            call diagnostics%print()
        end if

        if(mod(it,nstep_write) == 0) then

            select type(state)
            class is (stvec_swm_t)

                select type(operator)
                class is (operator_swm_t)
                    call operator%curl_op%calc_curl(curl, state%u, state%v, domain)
                end select
                call outputer_curl%write(curl, domain, 'curl.dat', int(it/nstep_write)+1)
                call outputer%write(state%h, domain, 'h.dat', int(it/nstep_write)+1)
                call outputer_vec%write(state%u,state%v,domain,"u.dat","v.dat",int(it/nstep_write)+1)
                l2err = l2norm(state%h, domain%mesh_p, domain%parcomm)
                if (parcomm_global%myid==0) print*, "Hours = ", real(time/3600 ,4), &
                                                    "L2err =", real(l2err,4), int(it/nstep_write)+1
            end select
        end if
    end do

end subroutine run_barotropic_inst

subroutine get_exact_solution(state, domain, v_components_type)

    use test_fields_mod, only : set_vector_test_field, &
                                set_scalar_test_field

    class(stvec_t), intent(inout) :: state
    type(domain_t), intent(in)    :: domain
    character(*),   intent(in)    :: v_components_type

    select type(state)
    class is (stvec_swm_t)

        call set_scalar_test_field(state%h, height_field, domain%mesh_p, 0)
        call set_vector_test_field(state%u, state%v, velocity_field, &
              domain%mesh_u, domain%mesh_v, 0, v_components_type)

    end select


end subroutine get_exact_solution

function create_barotropic_instability_height_field_generator(H0, &
                              omega, grav, h_pert, u0, a, Nq) result (height_gen)

    use const_mod, only : pi
    use barotropic_instability_u_mod, only : barotropic_instability_u

    real(kind=8), intent(in) :: H0
    real(kind=8), intent(in) :: omega
    real(kind=8), intent(in) :: grav
    real(kind=8), intent(in) :: h_pert
    real(kind=8), intent(in) :: u0
    real(kind=8), intent(in) :: a
    integer(kind=4), intent(in) :: Nq

    type(barotropic_instability_height_generator_t) :: height_gen

    real(kind=8), parameter :: phi0 = pi/7d0, phi1 = .5d0*pi-phi0
    real(kind=8)            :: dphi
    real(kind=8)            :: phi, u, dh
    real(kind=8)            :: urhs
    integer(kind=4)         :: j, iq

    real(kind=8), parameter :: xq(3) = [-sqrt(0.6_8), 0.0_8, sqrt(0.6_8)]
    real(kind=8), parameter :: wq(3) = [5.0_8/9.0_8, 8.0_8/9.0_8, 5.0_8/9.0_8]

    height_gen%H0 = H0
    height_gen%Nq = Nq
    allocate(height_gen%H_zonal(0:Nq))
    dphi = (phi1 - phi0) / real(Nq,8)
    height_gen%dphi = dphi
    height_gen%h_pert = h_pert

    height_gen%H_zonal(0) = H0!0.0_8
    do j=1, Nq
        dh = 0.0_8
        do iq=1, size(wq)
            phi   = phi0+(j-0.5_8)*dphi+dphi*xq(iq)*0.5_8
            u = barotropic_instability_u(u0, phi)
            urhs  =-(u**2*tan(phi)+2._8*omega*a*sin(phi)*u) / grav
            dh = dh+wq(iq)*urhs
        end do
        height_gen%H_zonal(j) = height_gen%H_zonal(j-1)+0.5_8*dh*dphi
    end do

    height_gen%H_north = height_gen%H_zonal(Nq)

end function

end module barotropic_inst_mod
