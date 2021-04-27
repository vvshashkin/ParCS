module operator_swlin_mod

use operator_abstract_mod,       only : operator_abstract_t
use stvec_abstract_mod,          only : stvec_abstract_t
use stvec_swlin_mod,             only : stvec_swlin_t
use container_abstract_mod,      only : model_parameters_abstract_t
use parameters_swlin_mod,        only : parameters_swlin_t
use mesh_mod,                    only : mesh_t
use exchange_abstract_mod,       only : exchange_t
use ecs_halo_mod,                only : ecs_halo_t
use hor_difops_abstract_mod,     only : gradient, divergence

implicit none

type, extends(operator_abstract_t) :: operator_swlin_t

    integer(kind=4)                :: op_halo_width = 3

    class(exchange_t), allocatable :: exch_halo
    type(ecs_halo_t),  allocatable :: halo(:)

    procedure(gradient),   pointer, nopass :: grad_contra
    procedure(divergence), pointer, nopass :: div

    contains

    procedure, public :: act      => act
    procedure, public :: ext_halo => ext_halo

end type operator_swlin_t

contains

subroutine act(this, vout, vin, model_params)

    use const_mod, only : grav

    class(operator_swlin_t),             intent(inout) :: this
    class(stvec_abstract_t),             intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_abstract_t),             intent(in)    :: vin
    class(model_parameters_abstract_t),  intent(in)    :: model_params

    integer(kind=4) ind, is, ie, js, je, nvi, nvj
    integer(kind=4) isv, iev, jsv, jev
    integer(kind=4) mesh_isv, mesh_jsv, mesh_iev, mesh_jev, hw

    select type (model_params)
    class is (parameters_swlin_t)

    select type (vout)
    class is (stvec_swlin_t)
        select type (vin)
        class is (stvec_swlin_t)

            do ind = model_params%ts, model_params%te

                !d vec{u}/dt = -grav * nabla(h)
                call this%grad_contra(vout%u%block(ind), vout%v%block(ind), vin%h%block(ind), model_params%mesh(ind), -grav)
                !dh/dt = -H0 * nabla*u
                call this%div(vout%h%block(ind),vin%u%block(ind),vin%v%block(ind), model_params%mesh(ind),-model_params%H0)

            end do

            call this%ext_halo(vout, model_params%ts, model_params%te)

        class default
            call avost("swlin operator: non-swlin state vector on input. Stop")
        end select
    class default
        call avost("swlin operator: non-swlin state vector on output. Stop")
    end select

    class default
        call avost("sqlin operator: inconsistent model parameters type")
    end select

    !print *, type(vout)

end subroutine act

subroutine ext_halo(this, v, ts, te)
    class(operator_swlin_t),    intent(inout) :: this
    class(stvec_swlin_t),       intent(inout) :: v
    integer(kind=4),            intent(in)    :: ts, te

    integer(kind=4) ind

    call this%exch_halo%do(v%h)
    do ind = ts, te
        call this%halo(ind)%interp(v%h%block(ind),this%op_halo_width)
    end do

end subroutine ext_halo

end module operator_swlin_mod
