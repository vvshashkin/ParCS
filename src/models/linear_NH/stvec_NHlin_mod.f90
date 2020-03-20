module stvec_NHlin_mod

use stvec_abstract_mod, only : stvec_abstract_t
use grid_function_mod,  only : grid_function_t

implicit none

type, extends(stvec_abstract_t) :: stvec_NHlin_t
    integer(kind=4) ts, te
    type(grid_function_t), allocatable :: h(:)
    type(grid_function_t), allocatable :: u(:)
    type(grid_function_t), allocatable :: v(:)
    logical                            :: init_and_alloc = .false.

    contains
    procedure, public :: add  => add
    procedure, public :: copy => copy
    procedure, public :: dot  => dot
end type stvec_NHlin_t

interface init_stvec_NHlin
    module procedure init_stvec_NHlin_tiles
    module procedure init_stvec_NHlin_ints
end interface
contains

subroutine init_stvec_NHlin_tiles(new_stvec, ts, te, tiles, halo_width)
    use tile_mod, only : tile_t

    type(stvec_NHlin_t), intent(out) :: new_stvec
    integer(kind=4),     intent(in)  :: ts, te
    type(tile_t),        intent(in)  :: tiles(ts:te)
    integer(kind=4),     intent(in)  :: halo_width

    integer(kind=4), allocatable :: is(:), ie(:), js(:), je(:)
    integer(kind=4), allocatable :: ks(:), ke(:), panel_ind(:)
    integer(kind=4) ind

    panel_ind = tiles(ts:te)%panel_number
    is = tiles(ts:te)%is; ie = tiles(ts:te)%ie
    js = tiles(ts:te)%js; je = tiles(ts:te)%je
    ks = tiles(ts:te)%ks; ke = tiles(ts:te)%ke

    call init_stvec_NHlin_ints(new_stvec, ts, te, panel_ind, is, ie, js, je, ks, ke,  &
                               halo_width)

end subroutine init_stvec_NHlin_tiles

subroutine init_stvec_NHlin_ints(new_stvec, ts, te, panel_ind, is, ie, js, je, ks, ke, halo_width)

    type(stvec_NHlin_t), intent(out) :: new_stvec
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
        new_stvec%h(ind)%p(:,:,:) = 0._8
        call new_stvec%u(ind)%init(panel_ind(ind),is(ind), ie(ind),    &
                                   js(ind), je(ind), ks(ind), ke(ind), &
                                   halo_width, halo_width, 0)
        new_stvec%u(ind)%p(:,:,:) = 0._8
        call new_stvec%v(ind)%init(panel_ind(ind),is(ind), ie(ind),    &
                                   js(ind), je(ind), ks(ind), ke(ind), &
                                   halo_width, halo_width, 0)
        new_stvec%v(ind)%p(:,:,:) = 0._8

    end do

    new_stvec%init_and_alloc = .true.

end subroutine init_stvec_NHlin_ints

subroutine add(this,other,alpha,beta)
    class(stvec_NHlin_t),    intent(inout)  :: this
    class(stvec_abstract_t), intent(in)     :: other
    real(kind=8),            intent(in)     :: alpha, beta

    integer(kind=4) ts, te
    integer(kind=4) i1, i2, j1, j2, k1, k2
    integer(kind=4) ind


    select type (other)
    class is (stvec_NHlin_t)
        ts = this%ts; te = this%te
        do ind = ts,te
            j1 = this%h(ind)%js-this%h(ind)%nvj
            j2 = this%h(ind)%je+this%h(ind)%nvj
            i1 = this%h(ind)%is-this%h(ind)%nvi
            i2 = this%h(ind)%ie+this%h(ind)%nvi
            k1 = this%h(ind)%ks
            k2 = this%h(ind)%ke
            this%h(ind)%p(i1:i2,j1:j2,k1:k2) = alpha*this%h(ind)%p(i1:i2,j1:j2,k1:k2) + &
                                               beta*other%h(ind)%p(i1:i2,j1:j2,k1:k2)
            this%u(ind)%p(i1:i2,j1:j2,k1:k2) = alpha*this%u(ind)%p(i1:i2,j1:j2,k1:k2) + &
                                               beta*other%u(ind)%p(i1:i2,j1:j2,k1:k2)
            this%v(ind)%p(i1:i2,j1:j2,k1:k2) = alpha*this%v(ind)%p(i1:i2,j1:j2,k1:k2) + &
                                               beta*other%v(ind)%p(i1:i2,j1:j2,k1:k2)
        end do
    class default
        call avost("NHlin_stvec_t%add types mismatch. stop!")
    end select
end subroutine add

subroutine copy(this,source_stvec)
    class(stvec_NHlin_t),   intent(inout)  :: this
    class(stvec_abstract_t), intent(in)    :: source_stvec

    integer(kind=4), allocatable :: is(:), ie(:), js(:), je(:)
    integer(kind=4), allocatable :: ks(:), ke(:), panel_ind(:)
    integer(kind=4) ts, te, halo_width

    integer(kind=4) i1, i2, j1, j2, k1, k2
    integer(kind=4) ind

    select type (source_stvec)
    class is (stvec_NHlin_t)
        ts = source_stvec%ts; te = source_stvec%te
        if(.not. this%init_and_alloc) then
            panel_ind = source_stvec%h(ts:te)%panel_ind
            is = source_stvec%h(ts:te)%is; ie = source_stvec%h(ts:te)%ie
            js = source_stvec%h(ts:te)%js; je = source_stvec%h(ts:te)%je
            ks = source_stvec%h(ts:te)%ks; ke = source_stvec%h(ts:te)%ke
            halo_width = max(source_stvec%h(ts)%nvi, source_stvec%h(ts)%nvj)
            call init_stvec_NHlin(this, ts, te, panel_ind, is, ie, js,   &
                                  je, ks, ke, halo_width)
        end if
        do ind = ts,te
            j1 = source_stvec%h(ind)%js-source_stvec%h(ind)%nvj
            j2 = source_stvec%h(ind)%je+source_stvec%h(ind)%nvj
            i1 = source_stvec%h(ind)%is-source_stvec%h(ind)%nvi
            i2 = source_stvec%h(ind)%ie+source_stvec%h(ind)%nvi
            k1 = source_stvec%h(ind)%ks
            k2 = source_stvec%h(ind)%ke
            this%h(ind)%p(i1:i2,j1:j2,k1:k2) = source_stvec%h(ind)%p(i1:i2,j1:j2,k1:k2)
            this%u(ind)%p(i1:i2,j1:j2,k1:k2) = source_stvec%u(ind)%p(i1:i2,j1:j2,k1:k2)
            this%v(ind)%p(i1:i2,j1:j2,k1:k2) = source_stvec%v(ind)%p(i1:i2,j1:j2,k1:k2)
        end do
    class default
        call avost("NHlin_stvec_t%copy types mismatch. stop!")
    end select
end subroutine copy

real(kind=8) function dot(this, other) result(dot_prod)
    class(stvec_NHlin_t),   intent(in)     :: this
    class(stvec_abstract_t), intent(in)    :: other

    dot_prod = 0._8

    !select type (other)
    !class is (stvec_iomega_t)
    !    dot_prod = sum(this%f(1:this%N)*other%f(1:this%N))
    !class default
    !    stop
    !end select

end function dot

end module stvec_NHlin_mod
