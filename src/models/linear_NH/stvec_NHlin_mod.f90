module stvec_NHlin_mod

use stvec_abstract_mod, only : stvec_abstract_t
use grid_function_mod,  only : grid_function_t

implicit none

type, extends(stvec_abstract_t) :: stvec_NHlin_t
    integer(kind=4) ts, te
    type(grid_function_t), allocatable :: prex(:)
    type(grid_function_t), allocatable :: u(:)
    type(grid_function_t), allocatable :: v(:)
    type(grid_function_t), allocatable :: w(:)
    type(grid_function_t), allocatable :: theta(:)
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

    allocate(new_stvec%prex(ts:te))
    allocate(new_stvec%u(ts:te))
    allocate(new_stvec%v(ts:te))
    allocate(new_stvec%w(ts:te))
    allocate(new_stvec%theta(ts:te))

    new_stvec%ts = ts; new_stvec%te = te

    do ind = ts, te
        call new_stvec%prex(ind)%init(panel_ind(ind),is(ind), ie(ind),    &
                                      js(ind), je(ind), ks(ind), ke(ind), &
                                      halo_width, halo_width, 0)
        new_stvec%prex(ind)%p(:,:,:) = 0._8

        call new_stvec%u(ind)%init(panel_ind(ind),is(ind), ie(ind),    &
                                   js(ind), je(ind), ks(ind), ke(ind), &
                                   halo_width, halo_width, 0)
        new_stvec%u(ind)%p(:,:,:) = 0._8

        call new_stvec%v(ind)%init(panel_ind(ind),is(ind), ie(ind),    &
                                   js(ind), je(ind), ks(ind), ke(ind), &
                                   halo_width, halo_width, 0)

        new_stvec%v(ind)%p(:,:,:) = 0._8
        call new_stvec%w(ind)%init(panel_ind(ind),is(ind), ie(ind),     &
                                  js(ind), je(ind), ks(ind)-1, ke(ind), &
                                  halo_width, halo_width, 0)
        new_stvec%w(ind)%p(:,:,:) = 0._8

        call new_stvec%theta(ind)%init(panel_ind(ind),is(ind), ie(ind),     &
                                   js(ind), je(ind), ks(ind)-1, ke(ind), &
                                   halo_width, halo_width, 0)
        new_stvec%theta(ind)%p(:,:,:) = 0._8
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
            j1 = this%prex(ind)%js-1!-this%prex(ind)%nvj
            j2 = this%prex(ind)%je+1!+this%prex(ind)%nvj
            i1 = this%prex(ind)%is-1!-this%prex(ind)%nvi
            i2 = this%prex(ind)%ie+1!+this%prex(ind)%nvi
            k1 = this%prex(ind)%ks
            k2 = this%prex(ind)%ke
            this%prex(ind)%p(i1-1:i2+1,j1-1:j2+1,k1:k2)    = alpha*this%prex(ind)%p(i1-1:i2+1,j1-1:j2+1,k1:k2)    + &
                                                     beta*other%prex(ind)%p(i1-1:i2+1,j1-1:j2+1,k1:k2)
            this%u(ind)%p(i1-1:i2,j1:j2,k1:k2)       = alpha*this%u(ind)%p(i1-1:i2,j1:j2,k1:k2)       + &
                                                     beta*other%u(ind)%p(i1-1:i2,j1:j2,k1:k2)
            this%v(ind)%p(i1:i2,j1-1:j2,k1:k2)       = alpha*this%v(ind)%p(i1:i2,j1-1:j2,k1:k2)       + &
                                                     beta*other%v(ind)%p(i1:i2,j1-1:j2,k1:k2)
            this%w(ind)%p(i1:i2,j1:j2,k1-1:k2)     = alpha*this%w(ind)%p(i1:i2,j1:j2,k1-1:k2)     + &
                                                     beta*other%w(ind)%p(i1:i2,j1:j2,k1-1:k2)
            this%theta(ind)%p(i1:i2,j1:j2,k1-1:k2) = alpha*this%theta(ind)%p(i1:i2,j1:j2,k1-1:k2) + &
                                                     beta*other%theta(ind)%p(i1:i2,j1:j2,k1-1:k2)
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
            panel_ind = source_stvec%prex(ts:te)%panel_ind
            is = source_stvec%prex(ts:te)%is; ie = source_stvec%prex(ts:te)%ie
            js = source_stvec%prex(ts:te)%js; je = source_stvec%prex(ts:te)%je
            ks = source_stvec%prex(ts:te)%ks; ke = source_stvec%prex(ts:te)%ke
            halo_width = max(source_stvec%prex(ts)%nvi, source_stvec%prex(ts)%nvj)
            call init_stvec_NHlin(this, ts, te, panel_ind, is, ie, js,   &
                                  je, ks, ke, halo_width)
        end if
        do ind = ts,te
            j1 = source_stvec%prex(ind)%js!-1!-source_stvec%prex(ind)%nvj
            j2 = source_stvec%prex(ind)%je!+1!+source_stvec%prex(ind)%nvj
            i1 = source_stvec%prex(ind)%is!-1!-source_stvec%prex(ind)%nvi
            i2 = source_stvec%prex(ind)%ie!+1!+source_stvec%prex(ind)%nvi
            k1 = source_stvec%prex(ind)%ks
            k2 = source_stvec%prex(ind)%ke
            this%prex(ind)%p(i1-1:i2+1,j1-1:j2+1,k1:k2)    = source_stvec%prex(ind)%p(i1-1:i2+1,j1-1:j2+1,k1:k2)
            this%u(ind)%p(i1-1:i2,j1:j2,k1:k2)       = source_stvec%u(ind)%p(i1-1:i2,j1:j2,k1:k2)
            this%v(ind)%p(i1:i2,j1-1:j2,k1:k2)       = source_stvec%v(ind)%p(i1:i2,j1-1:j2,k1:k2)
            this%w(ind)%p(i1:i2,j1:j2,k1-1:k2)     = source_stvec%w(ind)%p(i1:i2,j1:j2,k1-1:k2)
            this%theta(ind)%p(i1:i2,j1:j2,k1-1:k2) = source_stvec%theta(ind)%p(i1:i2,j1:j2,k1-1:k2)
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
