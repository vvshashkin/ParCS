module operator_swlin_mod

use operator_abstract_mod, only : operator_abstract_t
use stvec_abstract_mod,    only : stvec_abstract_t
use stvec_swlin_mod,       only : stvec_swlin_t
use mesh_mod,              only : mesh_t
use exchange_mod,          only : exchange_t
use ecs_halo_mod,          only : ecs_halo_t

implicit none

integer, parameter :: halo_width = 3

type, extends(operator_abstract_t) :: operator_swlin_t

    integer(kind=4)                :: ts, te
    type(mesh_t), pointer          :: mesh(:)
    type(exchange_t)               :: exch_halo
    type(ecs_halo_t), allocatable  :: halo(:)
    real(kind=8)                   :: H0

    contains

    procedure, public :: act      => act
    procedure, public :: ext_halo => ext_halo

end type operator_swlin_t

contains

subroutine init_swlin_operator(oper, ts, te, mesh, partition, ex_halo_width, myid, np, H0)

    use partition_mod,        only : partition_t
    use exchange_factory_mod, only : create_2d_full_halo_exchange
    use ecs_halo_factory_mod, only : init_ecs_halo

    type(operator_swlin_t), intent(out) :: oper
    integer(kind=4),        intent(in)  :: ts, te
    type(mesh_t),           intent(in), &
                            target      :: mesh(ts:te)
    type(partition_t),      intent(in)  :: partition
    integer(kind=4),        intent(in)  :: ex_halo_width
    integer(kind=4),        intent(in)  :: myid, np
    real(kind=8),           intent(in)  :: H0

    integer(kind=4) ind

    oper%ts = ts; oper%te = te
    oper%mesh(ts:te) => mesh(ts:te)

    call create_2d_full_halo_exchange(oper%exch_halo, partition, ex_halo_width, myid, np)

    allocate(oper%halo(ts:te))
    do ind=ts,te
        oper%halo(ind) = init_ecs_halo(mesh(ind)%is, mesh(ind)%ie, &
                                       mesh(ind)%js, mesh(ind)%je, &
                                       mesh(ind)%nx, halo_width,   &
                                       mesh(ind)%hx)
    end do

    oper%H0 = H0

end subroutine init_swlin_operator

subroutine act(this, vout, vin)

    class(operator_swlin_t),    intent(inout) :: this
    class(stvec_abstract_t),    intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_abstract_t),    intent(in)    :: vin

    integer(kind=4) ind, is, ie, js, je, nvi, nvj
    integer(kind=4) isv, iev, jsv, jev
    integer(kind=4) mesh_isv, mesh_jsv, mesh_iev, mesh_jev, hw

    select type (vout)
    class is (stvec_swlin_t)
        select type (vin)
        class is (stvec_swlin_t)
                do ind = this%ts, this%te
                    nvi = vin%h(ind)%nvi;    nvj = vin%h(ind)%nvj
                    js  = vin%h(ind)%js;     jsv = js - nvj
                    je  = vin%h(ind)%je;     jev = je + nvj
                    is  = vin%h(ind)%is;     isv = is - nvi
                    ie  = vin%h(ind)%ie;     iev = ie + nvi

                    hw = this%mesh(ind)%halo_width
                    mesh_isv = this%mesh(ind)%is - hw
                    mesh_iev = this%mesh(ind)%ie + hw
                    mesh_jsv = this%mesh(ind)%js - hw
                    mesh_jev = this%mesh(ind)%je + hw

                    call apply_swlin_oper(vout%h(ind)%p, vout%u(ind)%p, vout%v(ind)%p,  &
                                          vin%h(ind)%p,  vin%u(ind)%p,  vin%v(ind)%p,   &
                                          is, ie, js, je, isv, iev, jsv, jev,           &
                                          this%mesh(ind)%G,  this%mesh(ind)%Gu,         &
                                          this%mesh(ind)%Gv, this%mesh(ind)%QI,         &
                                          this%mesh(ind)%QIu, this%mesh(ind)%QIv,       &
                                          mesh_isv, mesh_iev, mesh_jsv, mesh_jev,       &
                                          this%mesh(ind)%hx, this%H0)

                end do
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

subroutine apply_swlin_oper(h1, u1, v1, h, u, v,                &
                            is, ie, js, je, isv, iev, jsv, jev, &
                            G, Gu, Gv, QI, QIu, QIv,            &
                            misv, miev, mjsv, mjev, hx, H0)

    use const_mod, only : grav, radz

    integer(kind=4), intent(in)  :: is, ie, js, je, isv, iev, jsv, jev !dimensions of grid arrays
    real(kind=8),    intent(out) :: h1(isv:iev,jsv:jev), u1(isv:iev,jsv:jev), v1(isv:iev,jsv:jev) !output tendencies
    real(kind=8),    intent(in)  :: h(isv:iev,jsv:jev), u(isv:iev,jsv:jev), v(isv:iev,jsv:jev)    !input variables
    integer(kind=4), intent(in)  :: misv, miev, mjsv, mjev
    !mesh parameters:
    real(kind=8),    intent(in)  :: G(misv:miev,mjsv:mjev)    !sqrt(det(Q)), where Q == metric tensor
    real(kind=8),    intent(in)  :: Gu(misv-1:miev,mjsv:mjev) !the same at u-flux points
    real(kind=8),    intent(in)  :: Gv(misv:miev,mjsv-1:mjev) !the same at v-flux points
    real(kind=8),    intent(in)  :: QI(3,misv:miev,mjsv:mjev)    !inversed metric tensor components
    real(kind=8),    intent(in)  :: QIu(3,misv-1:miev,mjsv:mjev) !the same at u-flux points
    real(kind=8),    intent(in)  :: QIv(3,misv:miev,mjsv-1:mjev) !the same at v-flux points
    real(kind=8),    intent(in)  :: hx !mesh grid spacing
    real(kind=8),    intent(in)  :: H0 !mean background water height
    !local
    integer(kind=4) i, j
    real(kind=8) gx(is-1:ie,js-1:je+1), gy(is-1:ie+1,js-1:je)

    !gradient covariant components
    do j=js-1, je+1
        do i=is-1, ie
            gx(i,j) =-grav*(h(i+1,j)-h(i,j))/(hx*radz)
        end do
    end do

    do j=js-1, je
        do i=is-1, ie+1
            gy(i,j) =-grav*(h(i,j+1)-h(i,j))/(hx*radz)
        end do
    end do

    !transform to contravariant components
    do j=js,je
        do i=is-1,ie
            u1(i,j) = Qiu(1,i,j)*gx(i,j)+Qiu(2,i,j)*0.25_8*(gy(i,j)+gy(i+1,j)+gy(i,j-1)+gy(i+1,j-1))
        end do
    end do
    do j=js-1,je
        do i=is,ie
            v1(i,j) = Qiv(3,i,j)*gy(i,j)+Qiv(2,i,j)*0.25_8*(gx(i,j)+gx(i-1,j)+gx(i,j+1)+gx(i-1,j+1))
        end do
    end do

    do j=js, je
        do i=is, ie
            h1(i,j) =-H0*(Gu(i,j)*u(i,j)-Gu(i-1,j)*u(i-1,j)+  &
                          Gv(i,j)*v(i,j)-Gv(i,j-1)*v(i,j-1))/(G(i,j)*radz*hx)
        end do
    end do
end subroutine apply_swlin_oper

end module operator_swlin_mod
