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

implicit none

contains

subroutine test_swm()

    type(config_swm_t) :: config

    type(domain_t) :: domain

    class(stvec_t),      allocatable :: state, state_err
    class(operator_t),   allocatable :: operator
    class(timescheme_t), allocatable :: timescheme
    class(outputer_t),   allocatable :: outputer

    integer(kind=4) :: halo_width = 8

    call config%parse("namelist_swm")

    call create_domain(domain, config%topology_type, config%staggering_type, &
                               config%N, config%nz)

    call create_stvec_swm(state,     domain, halo_width, 0)
    call create_stvec_swm(state_err, domain, 0         , 0)

    call create_swm_operator(operator, config, domain)

    call create_timescheme(timescheme, state, 'rk4')

end subroutine test_swm

end module test_swm_mod
