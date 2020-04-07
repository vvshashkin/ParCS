module operator_NHlin_split_mod

use operator_NHlin_mod,          only : operator_NHlin_t
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

type, extends(operator_NHlin_t) :: operator_NHlin_explicit_t

    contains

    procedure, public :: act => act_explicit

end type operator_NHlin_explicit_t

type, extends(operator_NHlin_t) :: operator_NHlin_implicit_t

    real(kind=8), allocatable :: d2z(:,:) !coefficients of d/dz*(P0*d/dz) for w-eq
    real(kind=8), allocatable :: cwz(:)    !coefficients of diag part of w-eq

    contains

    procedure, public :: act  => act_implicit
    procedure, public :: solv => solv_implicit

end type operator_NHlin_implicit_t

contains

subroutine init_NHlin_explicit_operator(oper, model_params, master_id, myid, np, namelist_str)
    use operator_NHlin_mod, only : init_NHlin_operator

    type(operator_NHlin_explicit_t)        :: oper
    type(parameters_NHlin_t),  intent(in)  :: model_params
    integer(kind=4),           intent(in)  :: master_id, myid, np
    character(:), allocatable, intent(in)  :: namelist_str

    !fall back to full operator inititalization
    call init_NHlin_operator(oper, model_params, master_id, myid, np, namelist_str)

end subroutine init_NHlin_explicit_operator

subroutine init_NHlin_implicit_operator(oper, model_params, master_id, myid, np, namelist_str)
    use operator_NHlin_mod, only : init_NHlin_operator
    use const_mod,          only : Cp, Cv, rgaz

    type(operator_NHlin_implicit_t)        :: oper
    type(parameters_NHlin_t),  intent(in)  :: model_params
    integer(kind=4),           intent(in)  :: master_id, myid, np
    character(:), allocatable, intent(in)  :: namelist_str

    integer(kind=4) k

    !fall back to full operator inititalization
    call init_NHlin_operator(oper, model_params, master_id, myid, np, namelist_str)
    !3-diag system initialize
    allocate(oper%d2z(-1:1,1:model_params%nz-1))
    allocate(oper%cwz(1:model_params%nz-1))

    oper%d2z(-1,1) = 0._8
    oper%d2z(0,1) =-(model_params%prex0(1)+model_params%prex0(2)) / model_params%dz**2
    oper%d2z(1,1) = (model_params%prex0(2)) / model_params%dz**2
    oper%cwz(1)   =-Cv/(Cp*rgaz*model_params%theta0(1))
    do k=2, model_params%nz-2
        oper%d2z(-1,k) = (model_params%prex0(k)) / model_params%dz**2
        oper%d2z( 0,k) =-(model_params%prex0(k)+model_params%prex0(k+1)) / model_params%dz**2
        oper%d2z( 1,k) = (model_params%prex0(k+1)) / model_params%dz**2
        oper%cwz(k)    =-Cv/(Cp*rgaz*model_params%theta0(k))
    end do
    k = model_params%nz-1
    oper%d2z(-1,k) = (model_params%prex0(k)) / model_params%dz**2
    oper%d2z( 0,k) =-(model_params%prex0(k)+model_params%prex0(k+1)) / model_params%dz**2
    oper%d2z( 1,k) = 0._8
    oper%cwz(k)    =-Cv/(Cp*rgaz*model_params%theta0(k))

end subroutine init_NHlin_implicit_operator

subroutine act_explicit(this, vout, vin, model_params)

    use const_mod, only : grav, rgaz, Cp, Cv

    class(operator_NHlin_explicit_t),    intent(inout) :: this
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
                w_at_p = 0.5_8*(vin%w(ind)%p(i,j,k)+vin%w(ind)%p(i,j,k-1))
                vout%prex(ind)%p(i,j,k) = -model_params%prex0(k)*rgaz/Cv*(vout%prex(ind)%p(i,j,k)) - w_at_p*model_params%prex0dz(k)
            end do; end do
        end do

        vout%theta(ind)%p(is:ie,js:je,0) = 0._8
        vout%w(ind)%p(is:ie,js:je,0) = 0._8
        do k = 1, model_params%nz-1
            dp0dz = (model_params%prex0(k+1)-model_params%prex0(k)) / model_params%dz
            do j=js,je; do i=is,ie
                vout%theta(ind)%p(i,j,k) = -vin%w(ind)%p(i,j,k)*model_params%theta0dz(k)
                vout%w(ind)%p(i,j,k) =-Cp*(vin%theta(ind)%p(i,j,k)*dp0dz)
            end do; end do
        end do
        vout%theta(ind)%p(is:ie,js:je,model_params%nz) = 0._8
        vout%w(ind)%p(is:ie,js:je,model_params%nz) = 0._8

    end do

    call this%ext_halo(vout, model_params%ts, model_params%te)

    class default
        call avost("NHlin operator: non-NHlin state vector on input. Stop")
    end select
    class default
        call avost("NHlin operator: non-NHlin state vector on output. Stop")
    end select

    class default
        call avost("NHlin operator: inconsistent model parameters type")
    end select

end subroutine act_explicit

subroutine act_implicit(this, vout, vin, model_params)

    use const_mod, only : grav, rgaz, Cp, Cv

    class(operator_NHlin_implicit_t),    intent(inout) :: this
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

        vout%u(ind)%p(is-1:ie,js:je,1:model_params%nz) = 0.
        vout%v(ind)%p(is:ie,js-1:je,1:model_params%nz) = 0.

        do k = 1, model_params%nz
            do j=js,je; do i=is,ie
                dwdz   = (vin%w(ind)%p(i,j,k)-vin%w(ind)%p(i,j,k-1)) / model_params%dz
                vout%prex(ind)%p(i,j,k) = -model_params%prex0(k)*rgaz/Cv*dwdz
            end do; end do
        end do

        vout%theta(ind)%p(is:ie,js:je,0:model_params%nz) = 0._8
        vout%w(ind)%p(is:ie,js:je,0) = 0._8
        do k = 1, model_params%nz-1
            do j=js,je; do i=is,ie
                dpdz = (vin%prex(ind)%p(i,j,k+1)-vin%prex(ind)%p(i,j,k) ) / model_params%dz
                vout%w(ind)%p(i,j,k) =-Cp*model_params%theta0(k)*dpdz
            end do; end do
        end do
        vout%w(ind)%p(is:ie,js:je,model_params%nz) = 0._8

    end do

    call this%ext_halo(vout, model_params%ts, model_params%te)

    class default
        call avost("NHlin operator: non-NHlin state vector on input. Stop")
    end select
    class default
        call avost("NHlin operator: non-NHlin state vector on output. Stop")
    end select

    class default
        call avost("NHlin operator: inconsistent model parameters type")
    end select

end subroutine act_implicit

subroutine solv_implicit(this, dt, vout, rhs, model_params)

    use const_mod, only : grav, rgaz, Cp, Cv

    class(operator_NHlin_implicit_t),    intent(inout) :: this
    real(kind=8),                        intent(in)    :: dt
    class(stvec_abstract_t),             intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_abstract_t),             intent(in)    :: rhs
    class(model_parameters_abstract_t),  intent(in)    :: model_params

    integer(kind=4) :: ind, i, j, k
    integer(kind=4) :: is, ie, js, je, nvi, nvj
    integer(kind=4) :: isv, iev, jsv, jev
    integer(kind=4) :: mesh_isv, mesh_jsv, mesh_iev, mesh_jev, hw
    real(kind=8)    :: dwdz, w_at_p, th0p, dpdz, dp0dz
    real(kind=8)    :: cf1, cf2

    select type (model_params)
    class is (parameters_NHlin_t)

    select type (vout)
    class is (stvec_NHlin_t)
        select type (rhs)
        class is (stvec_NHlin_t)

    do ind = model_params%ts, model_params%te

        call model_params%partition%tile(ind)%getind(is=is,ie=ie,js=js,je=je)

        vout%u(ind)%p(is-1:ie,js:je,1:model_params%nz) = rhs%u(ind)%p(is-1:ie,js:je,1:model_params%nz)
        vout%v(ind)%p(is:ie,js-1:je,1:model_params%nz) = rhs%v(ind)%p(is:ie,js-1:je,1:model_params%nz)
        vout%theta(ind)%p(is:ie,js:je,0:model_params%nz) = rhs%theta(ind)%p(is:ie,js:je,0:model_params%nz)

        vout%w(ind)%p(is:ie,js:je,0) = 0._8
        do k = 1, model_params%nz-1
            do j=js,je; do i=is,ie
                dpdz = (rhs%prex(ind)%p(i,j,k+1)-rhs%prex(ind)%p(i,j,k) ) / model_params%dz
                cf1 = Cv/(rgaz*dt)
                cf2 =-Cv/(rgaz*Cp*model_params%theta0(k)*dt**2)
                vout%w(ind)%p(i,j,k) = cf1*dpdz+cf2*rhs%w(ind)%p(i,j,k)
            end do; end do
        end do
        vout%w(ind)%p(is:ie,js:je,model_params%nz) = 0._8

        call solve_w(vout%w(ind),is,ie,js,je,model_params%nz,this%d2z,this%cwz,dt)

        do k = 1, model_params%nz
            do j=js,je; do i=is,ie
                dwdz   = (vout%w(ind)%p(i,j,k)-vout%w(ind)%p(i,j,k-1)) / model_params%dz
                vout%prex(ind)%p(i,j,k) = -model_params%prex0(k)*rgaz/Cv*dwdz*dt+rhs%prex(ind)%p(i,j,k)
            end do; end do
        end do

    end do

    call this%ext_halo(vout, model_params%ts, model_params%te)

    class default
        call avost("NHlin operator: non-NHlin state vector on input. Stop")
    end select
    class default
        call avost("NHlin operator: non-NHlin state vector on output. Stop")
    end select

    class default
        call avost("NHlin operator: inconsistent model parameters type")
    end select

end subroutine solv_implicit

subroutine solve_w(w,is,ie,js,je,nz,d2z,cwz,dt)
    use grid_function_mod, only : grid_function_t

    type(grid_function_t), intent(inout) :: w
    integer(kind=4),       intent(in)    :: is,ie,js,je,nz
    real(kind=8),          intent(in)    :: d2z(-1:1,1:nz-1)
    real(kind=8),          intent(in)    :: cwz(1:nz-1)
    real(kind=8),          intent(in)    :: dt

    integer(kind=4) :: i,j,k
    !real(kind=8)    :: wsol(is:ie,js:je,1:nz-1)
    real(kind=8)    :: wdiag(1:nz-1)

    wdiag(1) = d2z(0,1)+cwz(1)/dt**2
    do k=2, nz-1
        do j=js,je; do i=is,ie
            w%p(i,j,k) = w%p(i,j,k) - w%p(i,j,k-1)*d2z(-1,k)/wdiag(k-1)
            wdiag(k) = (d2z(0,k)+cwz(k)/dt**2)-d2z(-1,k)/wdiag(k-1)*d2z(1,k-1)
        end do; end do
    end do

    w%p(is:ie,js:je,nz-1) = w%p(is:ie,js:je,nz-1) / wdiag(nz-1)
    do k = nz-2,1,-1
        do j=js,je; do i=is,ie
            w%p(i,j,k) = (w%p(i,j,k) - w%p(i,j,k+1)*d2z(1,k))/wdiag(k)
        end do; end do
    end do

end subroutine solve_w

end module operator_NHlin_split_mod
