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

character(64) :: div_op_name = "divergence2", grad_op_name = "gradient2"
namelist /oper_ini/ div_op_name, grad_op_name

integer, parameter :: halo_width = 3

type, extends(operator_abstract_t) :: operator_swlin_t

    integer(kind=4)                :: ts, te
    type(mesh_t), pointer          :: mesh(:)
    class(exchange_t), allocatable :: exch_halo
    type(ecs_halo_t),  allocatable :: halo(:)
    real(kind=8)                   :: H0
    procedure(gradient),   pointer, nopass :: grad_contra
    procedure(divergence), pointer, nopass :: div

    contains

    procedure, public :: act      => act
    procedure, public :: ext_halo => ext_halo

end type operator_swlin_t

contains

function init_swlin_operator(ts, te, mesh, partition, ex_halo_width, &
                             master_id, myid, np, H0, namelist_str)    result(oper)

    use partition_mod,        only : partition_t
    use exchange_factory_mod, only : create_2d_halo_exchange
    use ecs_halo_factory_mod, only : init_ecs_halo
    use hor_difops_basic_mod, only : cl_gradient_contra_c2, cl_divergence_cgr2, &
                                     cl_gradient_0, cl_divergence_0

    type(operator_swlin_t)                 :: oper
    integer(kind=4),           intent(in)  :: ts, te
    type(mesh_t),     target,  intent(in)  :: mesh(ts:te)
    type(partition_t),         intent(in)  :: partition
    integer(kind=4),           intent(in)  :: ex_halo_width
    integer(kind=4),           intent(in)  :: master_id, myid, np
    real(kind=8),              intent(in)  :: H0
    character(:), allocatable, intent(in)  :: namelist_str

    integer(kind=4) ind

    oper%ts = ts; oper%te = te
    oper%mesh(ts:te) => mesh(ts:te)

    oper%exch_halo = create_2d_halo_exchange(partition, ex_halo_width, 'full', myid, np)

    allocate(oper%halo(ts:te))
    do ind=ts,te
        oper%halo(ind) = init_ecs_halo(mesh(ind)%is, mesh(ind)%ie, &
                                       mesh(ind)%js, mesh(ind)%je, &
                                       mesh(ind)%nx, halo_width,   &
                                       mesh(ind)%hx)
    end do

    oper%H0 = H0

    if(allocated(namelist_str)) then
        read(namelist_str, oper_ini)
    end if

    if(trim(grad_op_name) == "gradient2") then
        oper%grad_contra => cl_gradient_contra_c2
    else if(trim(grad_op_name) == "gradient0") then
        oper%grad_contra => cl_gradient_0
    else
        call avost("SWLIN operator init: unknown gradient operator -- " // &
                                                          trim(grad_op_name))
    end if

    if(trim(div_op_name) == "divergence2") then
        oper%div         => cl_divergence_cgr2
    else if(trim(div_op_name) == "divergence0") then
        oper%div         => cl_divergence_0
    else
        call avost("SWLIN operator init: unknown gradient operator -- " // &
                                                          trim(div_op_name))
    end if

    if(myid == master_id) then
        print *, "---Operator initialization"
        print *, "gradient operator: ", grad_op_name
        print *, "divergence operator: ", div_op_name
        print *, "--------------------------"
    end if

end function init_swlin_operator

subroutine act(this, vout, vin, model_params)

    use const_mod, only : grav

    class(operator_swlin_t),             intent(inout) :: this
    class(stvec_abstract_t),             intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_abstract_t),             intent(in)    :: vin
    class(model_parameters_abstract_t),  intent(in)    :: model_params

    integer(kind=4) ind, is, ie, js, je, nvi, nvj
    integer(kind=4) isv, iev, jsv, jev
    integer(kind=4) mesh_isv, mesh_jsv, mesh_iev, mesh_jev, hw

    select type (vout)
    class is (stvec_swlin_t)
        select type (vin)
        class is (stvec_swlin_t)

            do ind = this%ts, this%te

                !d vec{u}/dt = -grav * nabla(h)
                call this%grad_contra(vout%u(ind), vout%v(ind), vin%h(ind), this%mesh(ind), -grav)
                !dh/dt = -H0 * nabla*u
                call this%div(vout%h(ind),vin%u(ind),vin%v(ind), this%mesh(ind),-this%H0)

            end do

            call this%ext_halo(vout)

        class default
            call avost("swlin operator: non-swlin state vector on input. Stop")
        end select
    class default
        call avost("swlin operator: non-swlin state vector on output. Stop")
    end select

end subroutine act

subroutine ext_halo(this, v)
    class(operator_swlin_t),    intent(inout) :: this
    class(stvec_swlin_t),       intent(inout) :: v

    integer(kind=4) ind

    call this%exch_halo%do(v%h(this%ts:this%te), this%ts, this%te)
    do ind = this%ts, this%te
        call this%halo(ind)%interp(v%h(ind),halo_width)
    end do

end subroutine ext_halo

end module operator_swlin_mod
