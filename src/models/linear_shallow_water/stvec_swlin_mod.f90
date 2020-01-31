module stvec_swlin_mod

use stvec_abstract_mod, only : stvec_abstract_t
use grid_function_mod,  only : grid_function_t

implicit none

type, extends(stvec_abstract_t) :: stvec_swlin_t
    integer(kind=4) ts, te
    type(grid_function_t), allocatable :: h(:)
    type(grid_function_t), allocatable :: u(:)
    type(grid_function_t), allocatable :: v(:)
    logical                            :: init_and_alloc = .false.

    contains
    procedure, public :: add  => add
    procedure, public :: copy => copy
    procedure, public :: dot  => dot
end type stvec_swlin_t

contains

subroutine init_stvec_swlin(new_stvec, ts, te, panel_ind, is, ie, js, je, ks, ke, halo_width)

    type(stvec_swlin_t), intent(out) :: new_stvec
    integer(kind=4),     intent(in)  :: ts, te
    integer(kind=4),     intent(in)  :: panel_ind(ts:te)
    integer(kind=4),     intent(in)  :: is(ts:te), ie(ts:te)
    integer(kind=4),     intent(in)  :: js(ts:te), je(ts:te)
    integer(kind=4),     intent(in)  :: ks(ts:te), ke(ts:te)
    integer(kind=4),     intent(in)  :: halo_width

    integer ind

    allocate(new_stvec%h(ts:te))
    allocate(new_stvec%u(ts:te))
    allocate(new_stvec%v(ts:te))

    new_stvec%ts = ts; new_stvec%te = te

    do ind = ts, te
        call new_stvec%h(ind)%init(panel_ind(ind),is(ind), ie(ind),    &
                                   js(ind), je(ind), ks(ind), ke(ind), &
                                   halo_width, halo_width, 0)
        call new_stvec%u(ind)%init(panel_ind(ind),is(ind), ie(ind),    &
                                   js(ind), je(ind), ks(ind), ke(ind), &
                                   halo_width, halo_width, 0)
        call new_stvec%v(ind)%init(panel_ind(ind),is(ind), ie(ind),    &
                                   js(ind), je(ind), ks(ind), ke(ind), &
                                   halo_width, halo_width, 0)

    end do

end subroutine init_stvec_swlin

subroutine add(this,other,alpha,beta)
    class(stvec_swlin_t),    intent(inout)  :: this
    class(stvec_abstract_t), intent(in)     :: other
    real(kind=8),            intent(in)     :: alpha, beta

    !select type (other)
    !class is (stvec_iomega_t)
    !    this%f(1:this%N) = alpha*this%f(1:this%N)+beta*other%f(1:this%N)
    !class default
    !    stop
    !end select
end subroutine add

subroutine copy(this,source_stvec)
    class(stvec_swlin_t),   intent(inout)  :: this
    class(stvec_abstract_t), intent(in)    :: source_stvec

    !select type (source_stvec)
    !class is (stvec_iomega_t)
    !    call init_stvec_iomega(this,source_stvec%N,source_stvec%f)
    !class default
    !    stop
    !end select
end subroutine copy

real(kind=8) function dot(this, other) result(dot_prod)
    class(stvec_swlin_t),   intent(in)     :: this
    class(stvec_abstract_t), intent(in)    :: other

    dot_prod = 0._8

    !select type (other)
    !class is (stvec_iomega_t)
    !    dot_prod = sum(this%f(1:this%N)*other%f(1:this%N))
    !class default
    !    stop
    !end select

end function dot

end module stvec_swlin_mod
