module stvec_NHlin_mod

use stvec_abstract_mod, only : stvec_abstract_t
use grid_function_mod,  only : grid_function_t

implicit none

type :: stencil_dep_t
    integer(kind=4) pl, pr, pd, pu, pb, pt !needed halo_zones for prex field
    integer(kind=4) ul, ur, ud, uu, ub, ut !needed halo_zones for u field
    integer(kind=4) vl, vr, vd, vu, vb, vt !needed halo_zones for v field
    integer(kind=4) tl, tr, td, tu, tb, tt !for theta field
    integer(kind=4) wl, wr, wd, wu, wb, wt !for theta field
end type stencil_dep_t

type, extends(stvec_abstract_t) :: stvec_NHlin_t
    integer(kind=4) ts, te
    type(grid_function_t), allocatable :: prex(:)
    type(grid_function_t), allocatable :: u(:)
    type(grid_function_t), allocatable :: v(:)
    type(grid_function_t), allocatable :: w(:)
    type(grid_function_t), allocatable :: theta(:)
    type(stencil_dep_t)                :: dep
    logical                            :: init_and_alloc = .false.

    contains
    procedure, public  :: add  => add
    procedure, public  :: copy => copy
    procedure, public  :: dot  => dot
    procedure, public  :: smult
    procedure, public  :: norm
    procedure, private :: get_dep_zone
end type stvec_NHlin_t

interface init_stvec_NHlin
    module procedure init_stvec_NHlin_tiles
    module procedure init_stvec_NHlin_ints
end interface
contains

subroutine init_stvec_NHlin_tiles(new_stvec, ts, te, tiles, halo_width)
    use tile_mod,           only : tile_t

    type(stvec_NHlin_t),    intent(out) :: new_stvec
    integer(kind=4),        intent(in)  :: ts, te
    type(tile_t),           intent(in)  :: tiles(ts:te)
    integer(kind=4),        intent(in)  :: halo_width

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

    new_stvec%dep = stencil_dep_t(pl=-1, pr=1, pd=-1, pu=1, pb= 0, pt=0,  &
                                  ul=-1, ur=0, ud= 0, uu=0, ub= 0, ut=0,  &
                                  vl= 0, vr=0, vd=-1, vu=0, vb= 0, vt=0,  &
                                  tl= 0, tr=0, td= 0, tu=0, tb=-1, tt=0,  &
                                  wl= 0, wr=0, wd= 0, wu=0, wb=-1, wt=0 )

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

            call this%get_dep_zone("prex",ind,i1,i2,j1,j2,k1,k2)
            this%prex(ind)%p(i1:i2,j1:j2,k1:k2)  = alpha*this%prex(ind)%p(i1:i2,j1:j2,k1:k2)  + &
                                                   beta*other%prex(ind)%p(i1:i2,j1:j2,k1:k2)
            call this%get_dep_zone("u",ind,i1,i2,j1,j2,k1,k2)
            this%u(ind)%p(i1:i2,j1:j2,k1:k2)     = alpha*this%u(ind)%p(i1:i2,j1:j2,k1:k2)     + &
                                                   beta*other%u(ind)%p(i1:i2,j1:j2,k1:k2)
            call this%get_dep_zone("v",ind,i1,i2,j1,j2,k1,k2)
            this%v(ind)%p(i1:i2,j1:j2,k1:k2)     = alpha*this%v(ind)%p(i1:i2,j1:j2,k1:k2)     + &
                                                   beta*other%v(ind)%p(i1:i2,j1:j2,k1:k2)
            call this%get_dep_zone("w",ind,i1,i2,j1,j2,k1,k2)
            this%w(ind)%p(i1:i2,j1:j2,k1:k2)     = alpha*this%w(ind)%p(i1:i2,j1:j2,k1:k2)     + &
                                                   beta*other%w(ind)%p(i1:i2,j1:j2,k1:k2)
            call this%get_dep_zone("theta",ind,i1,i2,j1,j2,k1,k2)
            this%theta(ind)%p(i1:i2,j1:j2,k1:k2) = alpha*this%theta(ind)%p(i1:i2,j1:j2,k1:k2) + &
                                                   beta*other%theta(ind)%p(i1:i2,j1:j2,k1:k2)
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
            !this%dep = source%dep
            call this%get_dep_zone("prex",ind,i1,i2,j1,j2,k1,k2)
            this%prex(ind)%p(i1:i2,j1:j2,k1:k2)    = source_stvec%prex(ind)%p(i1:i2,j1:j2,k1:k2)
            call this%get_dep_zone("u",ind,i1,i2,j1,j2,k1,k2)
            this%u(ind)%p(i1:i2,j1:j2,k1:k2)       = source_stvec%u(ind)%p(i1:i2,j1:j2,k1:k2)
            call this%get_dep_zone("v",ind,i1,i2,j1,j2,k1,k2)
            this%v(ind)%p(i1:i2,j1:j2,k1:k2)       = source_stvec%v(ind)%p(i1:i2,j1:j2,k1:k2)
            call this%get_dep_zone("w",ind,i1,i2,j1,j2,k1,k2)
            this%w(ind)%p(i1:i2,j1:j2,k1:k2)     = source_stvec%w(ind)%p(i1:i2,j1:j2,k1:k2)
            call this%get_dep_zone("theta",ind,i1,i2,j1,j2,k1,k2)
            this%theta(ind)%p(i1:i2,j1:j2,k1:k2) = source_stvec%theta(ind)%p(i1:i2,j1:j2,k1:k2)
        end do
    class default
        call avost("NHlin_stvec_t%copy types mismatch. stop!")
    end select
end subroutine copy

subroutine smult(this, alpha)
    class(stvec_NHlin_t), intent(inout) :: this
    real(kind=8),         intent(in)    :: alpha

    integer(kind=4) ts, te!, halo_width
    integer(kind=4) i1, i2, j1, j2, k1, k2
    integer(kind=4) ind

    ts = this%ts; te = this%te

    do ind = ts, te
        call this%get_dep_zone("prex",ind,i1,i2,j1,j2,k1,k2)
        this%prex(ind)%p(i1:i2,j1:j2,k1:k2)  = alpha*this%prex(ind)%p(i1:i2,j1:j2,k1:k2)
        call this%get_dep_zone("u",ind,i1,i2,j1,j2,k1,k2)
        this%u(ind)%p(i1:i2,j1:j2,k1:k2)     = alpha*this%u(ind)%p(i1:i2,j1:j2,k1:k2)
        call this%get_dep_zone("v",ind,i1,i2,j1,j2,k1,k2)
        this%v(ind)%p(i1:i2,j1:j2,k1:k2)     = alpha*this%v(ind)%p(i1:i2,j1:j2,k1:k2)
        call this%get_dep_zone("w",ind,i1,i2,j1,j2,k1,k2)
        this%w(ind)%p(i1:i2,j1:j2,k1:k2)     = alpha*this%w(ind)%p(i1:i2,j1:j2,k1:k2)
        call this%get_dep_zone("theta",ind,i1,i2,j1,j2,k1,k2)
        this%theta(ind)%p(i1:i2,j1:j2,k1:k2) = alpha*this%theta(ind)%p(i1:i2,j1:j2,k1:k2)
    end do

end subroutine smult

real(kind=8) function dot(this, other) result(dot_prod)

    use mpi

    class(stvec_NHlin_t),   intent(in)     :: this
    class(stvec_abstract_t), intent(in)    :: other

    integer(kind=4) :: ind, ts, te, ierr
    integer(kind=4) :: i1, i2, j1, j2, k1, k2
    real(kind=8)    :: dot_prod_loc

    dot_prod_loc = 0._8
    select type (other)
    class is (stvec_NHlin_t)

    ts = this%ts; te = this%te
    do ind = ts, te
        call this%get_dep_zone("prex",ind,i1,i2,j1,j2,k1,k2)
        dot_prod_loc = dot_prod_loc + &
            sum(this%prex(ind)%p(i1:i2,j1:j2,k1:k2)*other%prex(ind)%p(i1:i2,j1:j2,k1:k2))

        call this%get_dep_zone("u",ind,i1,i2,j1,j2,k1,k2)
        dot_prod_loc = dot_prod_loc + &
            sum(this%u(ind)%p(i1+1:i2-1,j1:j2,k1:k2)*other%u(ind)%p(i1+1:i2-1,j1:j2,k1:k2)) +&
            0.5_8*sum(this%u(ind)%p(i1,j1:j2,k1:k2)*other%u(ind)%p(i1,j1:j2,k1:k2))         +&
            0.5_8*sum(this%u(ind)%p(i2,j1:j2,k1:k2)*other%u(ind)%p(i2,j1:j2,k1:k2))

        call this%get_dep_zone("v",ind,i1,i2,j1,j2,k1,k2)
        dot_prod_loc = dot_prod_loc + &
            sum(this%v(ind)%p(i1:i2,j1+1:j2-1,k1:k2)*other%v(ind)%p(i1:i2,j1+1:j2-1,k1:k2)) +&
            0.5_8*sum(this%v(ind)%p(i1:i2,j1,k1:k2)*other%v(ind)%p(i1:i2,j1,k1:k2))  +&
            0.5_8*sum(this%v(ind)%p(i1:i2,j2,k1:k2)*other%v(ind)%p(i1:i2,j2,k1:k2))

        call this%get_dep_zone("w",ind,i1,i2,j1,j2,k1,k2)
        dot_prod_loc = dot_prod_loc + &
            sum(this%w(ind)%p(i1:i2,j1:j2,k1:k2)*other%w(ind)%p(i1:i2,j1:j2,k1:k2))

        call this%get_dep_zone("theta",ind,i1,i2,j1,j2,k1,k2)
        dot_prod_loc = dot_prod_loc + &
            sum(this%theta(ind)%p(i1:i2,j1:j2,k1:k2)*other%theta(ind)%p(i1:i2,j1:j2,k1:k2))
    end do

    call mpi_allreduce(dot_prod_loc, dot_prod, 1, mpi_double, mpi_sum, mpi_comm_world, ierr)

    class default
        call avost("type mismatch in stvec_NHlin%dot")
    end select

end function dot

real(kind=8) function norm(this) result(l2)

    use mpi

    class(stvec_NHlin_t),    intent(in)    :: this

    integer(kind=4) :: ind, ts, te, ierr
    integer(kind=4) :: i1, i2, j1, j2, k1, k2
    real(kind=8)    :: l2_loc

    l2_loc = 0._8

    ts = this%ts; te = this%te
    do ind = ts, te
        call this%get_dep_zone("prex",ind,i1,i2,j1,j2,k1,k2)
        l2_loc = l2_loc + &
            sum(this%prex(ind)%p(i1:i2,j1:j2,k1:k2)**2)

        call this%get_dep_zone("u",ind,i1,i2,j1,j2,k1,k2)
        l2_loc = l2_loc + &
            sum(this%u(ind)%p(i1+1:i2-1,j1:j2,k1:k2)**2) +&
            0.5_8*sum(this%u(ind)%p(i1,j1:j2,k1:k2)**2)         +&
            0.5_8*sum(this%u(ind)%p(i2,j1:j2,k1:k2)**2)

        call this%get_dep_zone("v",ind,i1,i2,j1,j2,k1,k2)
        l2_loc = l2_loc + &
            sum(this%v(ind)%p(i1:i2,j1+1:j2-1,k1:k2)**2) +&
            0.5_8*sum(this%v(ind)%p(i1:i2,j1,k1:k2)**2)  +&
            0.5_8*sum(this%v(ind)%p(i1:i2,j2,k1:k2)**2)

        call this%get_dep_zone("w",ind,i1,i2,j1,j2,k1,k2)
        l2_loc = l2_loc + &
            sum(this%w(ind)%p(i1:i2,j1:j2,k1:k2)**2)

        call this%get_dep_zone("theta",ind,i1,i2,j1,j2,k1,k2)
        l2_loc = l2_loc + &
            sum(this%theta(ind)%p(i1:i2,j1:j2,k1:k2)**2)
    end do

    call mpi_allreduce(l2_loc, l2, 1, mpi_double, mpi_sum, mpi_comm_world, ierr)

    l2 = sqrt(l2)
    
end function norm

subroutine get_dep_zone(this,fld,ind,i1,i2,j1,j2,k1,k2)
    class(stvec_NHlin_t), intent(in)  :: this
    character(*),         intent(in)  :: fld
    integer(kind=4),      intent(in)  :: ind
    integer(kind=4),      intent(out) :: i1,i2,j1,j2,k1,k2

    if(fld == "prex") then
        i1 = this%prex(ind)%is+this%dep%pl; i2 = this%prex(ind)%ie+this%dep%pr
        j1 = this%prex(ind)%js+this%dep%pd; j2 = this%prex(ind)%je+this%dep%pu
        k1 = this%prex(ind)%ks+this%dep%pb; k2 = this%prex(ind)%ke+this%dep%pt
    else if(fld == "u") then
        i1 = this%prex(ind)%is+this%dep%ul; i2 = this%prex(ind)%ie+this%dep%ur
        j1 = this%prex(ind)%js+this%dep%ud; j2 = this%prex(ind)%je+this%dep%uu
        k1 = this%prex(ind)%ks+this%dep%ub; k2 = this%prex(ind)%ke+this%dep%ut
    else if(fld == "v") then
        i1 = this%prex(ind)%is+this%dep%vl; i2 = this%prex(ind)%ie+this%dep%vr
        j1 = this%prex(ind)%js+this%dep%vd; j2 = this%prex(ind)%je+this%dep%vu
        k1 = this%prex(ind)%ks+this%dep%vb; k2 = this%prex(ind)%ke+this%dep%vt
    else if(fld == "w") then
        i1 = this%prex(ind)%is+this%dep%wl; i2 = this%prex(ind)%ie+this%dep%wr
        j1 = this%prex(ind)%js+this%dep%wd; j2 = this%prex(ind)%je+this%dep%wu
        k1 = this%prex(ind)%ks+this%dep%wb; k2 = this%prex(ind)%ke+this%dep%wt

    else if(fld == "theta") then
        i1 = this%prex(ind)%is+this%dep%tl; i2 = this%prex(ind)%ie+this%dep%tr
        j1 = this%prex(ind)%js+this%dep%td; j2 = this%prex(ind)%je+this%dep%tu
        k1 = this%prex(ind)%ks+this%dep%tb; k2 = this%prex(ind)%ke+this%dep%tt
    else
        call avost("stvec_NHlin_t%get_dep_zone unknown field" // fld // " . Exit!")
    end if
end subroutine get_dep_zone

end module stvec_NHlin_mod
