module stvec_swlin_mod

use stvec_abstract_mod, only : stvec_abstract_t
use grid_field_mod,     only : grid_field_t

implicit none

type, extends(stvec_abstract_t) :: stvec_swlin_t
    integer(kind=4) ts, te
    type(grid_field_t) :: h
    type(grid_field_t) :: u
    type(grid_field_t) :: v
    logical            :: init_and_alloc = .false.

    contains
    procedure, public :: add  => add
    procedure, public :: copy => copy
    procedure, public :: dot  => dot
end type stvec_swlin_t

contains

subroutine add(this,other,alpha,beta)
    class(stvec_swlin_t),    intent(inout)  :: this
    class(stvec_abstract_t), intent(in)     :: other
    real(kind=8),            intent(in)     :: alpha, beta

    integer(kind=4) t, i1, i2, j1, j2

    select type (other)
    class is (stvec_swlin_t)
        do t = this%ts, this%te
            j1 = this%h%block(t)%js-this%h%block(t)%nvj
            j2 = this%h%block(t)%je+this%h%block(t)%nvj
            i1 = this%h%block(t)%is-this%h%block(t)%nvi
            i2 = this%h%block(t)%ie+this%h%block(t)%nvi
            this%h%block(t)%p(i1:i2,j1:j2,1) = alpha*this%h%block(t)%p(i1:i2,j1:j2,1) + &
                                               beta*other%h%block(t)%p(i1:i2,j1:j2,1)
            this%u%block(t)%p(i1:i2,j1:j2,1) = alpha*this%u%block(t)%p(i1:i2,j1:j2,1) + &
                                               beta*other%u%block(t)%p(i1:i2,j1:j2,1)
            this%v%block(t)%p(i1:i2,j1:j2,1) = alpha*this%v%block(t)%p(i1:i2,j1:j2,1) + &
                                               beta*other%v%block(t)%p(i1:i2,j1:j2,1)
        end do
    class default
        call avost("swlin_stvec_t%add types mismatch. stop!")
    end select
end subroutine add

subroutine copy(this,source_stvec)
    class(stvec_swlin_t),   intent(inout)  :: this
    class(stvec_abstract_t), intent(in)    :: source_stvec

    integer(kind=4), allocatable :: is(:), ie(:), js(:), je(:)
    integer(kind=4), allocatable :: ks(:), ke(:), panel_ind(:)
    integer(kind=4) ts, te, halo_width

    integer(kind=4) t, i1, i2, j1, j2

    select type (source_stvec)
    class is (stvec_swlin_t)
        ! if(.not. this%init_and_alloc) then
        !     panel_ind = source_stvec%h(ts:te)%panel_ind
        !     is = source_stvec%h(ts:te)%is; ie = source_stvec%h(ts:te)%ie
        !     js = source_stvec%h(ts:te)%js; je = source_stvec%h(ts:te)%je
        !     ks = source_stvec%h(ts:te)%ks; ke = source_stvec%h(ts:te)%ke
        !     halo_width = max(source_stvec%h(ts)%nvi, source_stvec%h(ts)%nvj)
        !     call init_stvec_swlin(this, ts, te, panel_ind, is, ie, js,   &
        !                           je, ks, ke, halo_width)
        ! end if
        do t = source_stvec%ts, source_stvec%te
            j1 = source_stvec%h%block(t)%js-source_stvec%h%block(t)%nvj
            j2 = source_stvec%h%block(t)%je+source_stvec%h%block(t)%nvj
            i1 = source_stvec%h%block(t)%is-source_stvec%h%block(t)%nvi
            i2 = source_stvec%h%block(t)%ie+source_stvec%h%block(t)%nvi
            this%h%block(t)%p(i1:i2,j1:j2,1) = source_stvec%h%block(t)%p(i1:i2,j1:j2,1)
            this%u%block(t)%p(i1:i2,j1:j2,1) = source_stvec%u%block(t)%p(i1:i2,j1:j2,1)
            this%v%block(t)%p(i1:i2,j1:j2,1) = source_stvec%v%block(t)%p(i1:i2,j1:j2,1)
        end do
    class default
        call avost("swlin_stvec_t%copy types mismatch. stop!")
    end select
end subroutine copy

real(kind=8) function dot(this, other) result(dot_prod)
    class(stvec_swlin_t),   intent(in)     :: this
    class(stvec_abstract_t), intent(in)    :: other

    print*, 'Dot is not implemented for stvec_swlin_t!!!'

    dot_prod = 0._8

end function dot

end module stvec_swlin_mod
