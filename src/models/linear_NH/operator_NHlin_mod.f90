module operator_NHlin_mod

use operator_abstract_mod,       only : operator_abstract_t
use stvec_abstract_mod,          only : stvec_abstract_t
use stvec_NHlin_mod,             only : stvec_NHlin_t
use container_abstract_mod,      only : model_parameters_abstract_t
use parameters_NHlin_mod,        only : parameters_NHlin_t
use mesh_mod,                    only : mesh_t
use exchange_abstract_mod,       only : exchange_t
use ecs_halo_mod,                only : ecs_halo_t
use hor_difops_abstract_mod,     only : gradient, divergence

implicit none

character(64) :: div_op_name = "divergence2", grad_op_name = "gradient2"
namelist /oper_ini/ div_op_name, grad_op_name

type, extends(operator_abstract_t) :: operator_NHlin_t

    integer, private               :: op_halo_width = 1

    class(exchange_t), allocatable :: exch_halo
    type(ecs_halo_t),  allocatable :: halo(:)

    procedure(gradient),   pointer, nopass :: grad_contra
    procedure(divergence), pointer, nopass :: div

    contains

    procedure, public :: act      => act
    procedure, public :: ext_halo => ext_halo

end type operator_NHlin_t

contains

subroutine init_NHlin_operator(oper, model_params, master_id, myid, np, namelist_str)

    use partition_mod,        only : partition_t
    use exchange_factory_mod, only : create_Agrid_halo_exchange
    use ecs_halo_factory_mod, only : init_ecs_halo
    use hor_difops_basic_mod, only : cl_gradient_contra_c2, cl_divergence_cgr2, &
                                     cl_gradient_0, cl_divergence_0

    class(operator_NHlin_t)                :: oper
    type(parameters_NHlin_t),  intent(in)  :: model_params
    integer(kind=4),           intent(in)  :: master_id, myid, np
    character(:), allocatable, intent(in)  :: namelist_str

    integer(kind=4) ind

    oper%exch_halo = create_Agrid_halo_exchange(model_params%partition,    &
                                                model_params%halo_width,   &
                                                'full', myid, np)
    allocate(oper%halo(model_params%ts:model_params%te))
    do ind = model_params%ts, model_params%te
        oper%halo(ind) = init_ecs_halo(model_params%mesh(ind)%is, model_params%mesh(ind)%ie, &
                                       model_params%mesh(ind)%js, model_params%mesh(ind)%je, &
                                       model_params%mesh(ind)%nx, max(oper%op_halo_width,2),&
                                       model_params%mesh(ind)%hx)
    end do

    if(allocated(namelist_str)) then
        read(namelist_str, oper_ini)
    end if

    if(trim(grad_op_name) == "gradient2") then
        oper%grad_contra => cl_gradient_contra_c2
    else if(trim(grad_op_name) == "gradient0") then
        oper%grad_contra => cl_gradient_0
    else
        call avost("NHlin operator init: unknown gradient operator -- " // &
                                                          trim(grad_op_name))
    end if

    if(trim(div_op_name) == "divergence2") then
        oper%div         => cl_divergence_cgr2
    else if(trim(div_op_name) == "divergence0") then
        oper%div         => cl_divergence_0
    else
        call avost("NHlin operator init: unknown gradient operator -- " // &
                                                          trim(div_op_name))
    end if

    if(myid == master_id) then
        print *, "---Operator initialization"
        print *, "gradient operator: ", grad_op_name
        print *, "divergence operator: ", div_op_name
        print *, "--------------------------"
    end if

end subroutine init_NHlin_operator

subroutine act(this, vout, vin, model_params)

    use const_mod, only : grav, rgaz, Cp, Cv

    class(operator_NHlin_t),             intent(inout) :: this
    class(stvec_abstract_t),             intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_abstract_t),             intent(in)    :: vin
    class(model_parameters_abstract_t),  intent(in)    :: model_params

    integer(kind=4) :: ind, i, j, k
    integer(kind=4) :: is, ie, js, je, nvi, nvj
    integer(kind=4) :: isv, iev, jsv, jev
    integer(kind=4) :: mesh_isv, mesh_jsv, mesh_iev, mesh_jev, hw
    real(kind=8)    :: dwdz, w_at_p, th0p, dpdz, dp0dz

    select type (model_params)
    class is (parameters_NHlin_t)

    select type (vout)
    class is (stvec_NHlin_t)
        select type (vin)
        class is (stvec_NHlin_t)

            do ind = model_params%ts, model_params%te

                call model_params%partition%tile(ind)%getind(is=is,ie=ie,js=js,je=je)

                call this%grad_contra(vout%u(ind), vout%v(ind), vin%prex(ind), model_params%mesh(ind), model_params%radx)
                call this%div(vout%prex(ind),vin%u(ind),vin%v(ind), model_params%mesh(ind), model_params%radx)

                do k = 1, model_params%nz
                    do j=js,je; do i=is-1,ie
                        vout%u(ind)%p(i,j,k) =-Cp*0.5_8*(model_params%theta0(k-1)+model_params%theta0(k))*vout%u(ind)%p(i,j,k)
                    end do; end do
                    do j=js-1,je; do i=is,ie
                        vout%v(ind)%p(i,j,k) =-Cp*0.5_8*(model_params%theta0(k-1)+model_params%theta0(k))*vout%v(ind)%p(i,j,k)
                    end do; end do
                    do j=js,je; do i=is,ie
                        dwdz   = (vin%w(ind)%p(i,j,k)-vin%w(ind)%p(i,j,k-1)) / model_params%dz
                        w_at_p = 0.5_8*(vin%w(ind)%p(i,j,k)+vin%w(ind)%p(i,j,k-1))
                        vout%prex(ind)%p(i,j,k) = -model_params%prex0(k)*rgaz/Cv*(vout%prex(ind)%p(i,j,k)+dwdz) - w_at_p*model_params%prex0dz(k)
                    end do; end do
                end do

                vout%theta(ind)%p(is:ie,js:je,0) = 0._8
                vout%w(ind)%p(is:ie,js:je,0) = 0._8
                do k = 1, model_params%nz-1
                    dp0dz = (model_params%prex0(k+1)-model_params%prex0(k)) / model_params%dz
                    do j=js,je; do i=is,ie

                        vout%theta(ind)%p(i,j,k) = -vin%w(ind)%p(i,j,k)*model_params%theta0dz(k)

                        dpdz = (vin%prex(ind)%p(i,j,k+1)-vin%prex(ind)%p(i,j,k) ) / model_params%dz
                        vout%w(ind)%p(i,j,k) =-Cp*(model_params%theta0(k)*dpdz + vin%theta(ind)%p(i,j,k)*dp0dz)

                    end do; end do
                end do
                vout%theta(ind)%p(is:ie,js:je,model_params%nz) = 0._8
                vout%w(ind)%p(is:ie,js:je,model_params%nz) = 0._8

                !vout%prex(ind)%p = 0._8
                !vout%u(ind)%p = 0._8
                !vout%v(ind)%p = 0._8
                !vout%w(ind)%p = 0._8
                !vout%theta(ind)%p = 0._8
            end do

            call this%ext_halo(vout, model_params%ts, model_params%te)

        class default
            call avost("NHlin operator: non-NHlin state vector on input. Stop")
        end select
    class default
        call avost("NHlin operator: non-NHlin state vector on output. Stop")
    end select

    class default
        call avost("sqlin operator: inconsistent model parameters type")
    end select

    !print *, type(vout)

end subroutine act

subroutine ext_halo(this, v, ts, te)
    class(operator_NHlin_t),    intent(inout) :: this
    class(stvec_NHlin_t),       intent(inout) :: v
    integer(kind=4),            intent(in)    :: ts, te

    integer(kind=4) ind

    call this%exch_halo%do(v%prex(ts:te), ts, te)
    do ind = ts, te
        call this%halo(ind)%interp(v%prex(ind),this%op_halo_width)
    end do

end subroutine ext_halo

end module operator_NHlin_mod
