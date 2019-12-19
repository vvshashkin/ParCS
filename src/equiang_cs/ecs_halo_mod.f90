!Module to interpolate values to virtual points beyond face edge
!Here we use two prototype faces to precompute weights:
!    __________
!   / source  /
!  / /^b     /
! / +-->a   /
!/_________/     ^z  
!|         |     |  /^y
!|  ^b     |     | /
!|  | targ |     |/
!|  +-->a  |     +------->x
!|_________|
module ecs_halo_mod

use halo_mod, only : halo_t

implicit none

!type-container for precomputed interpolation weights to have
!multiple halo-interpolators for multiple grids (for geometric MG-solver)

type, extends(halo_t) :: ecs_halo_t
  integer(kind=4) n          !number of real grid points along cubed-sphere face edge
  integer(kind=4) halo_width !number of rows in halo-zone
  integer(kind=4) ws, we     !starting and ending index of interpolation weights
  real(kind=8)   , allocatable :: w(:,:,:)
  integer(kind=4), allocatable :: ind(:,:)
  contains
  procedure, public :: interp  => ecs_ext_halo
  procedure, public :: interpv => ecs_ext_halo_vect
end type ecs_halo_t

contains

!subroutine ecs_halo_mod_init(nn,halo_width)
!!input
!integer nn(:)      !dimension of each cub-sph grid (6*nn(i)**2)
!integer halo_width !number of rows in halo-zone
!!locals:
!integer nmg   !number of grids
!integer i
!
!nmg = size(nn)
!
!allocate(whalo(nmg))
!
!do i=1,nmg
!    call init_ecs_halo_t(nn(i),halo_width,whalo(i))
!end do
!
!end subroutine ecs_halo_mod_init

type(ecs_halo_t) function init_ecs_halo(n,halo_width) result(halo)
use const_mod, only: pi

integer(kind=4), intent(in) :: n, halo_width

!local:
integer i,j
real(kind=8) zda !angle grid-step
real(kind=8) za, zb !angle coordinates at target grid face
real(kind=8) zx, zz !cartesian coordinates, zy is not needed
real(kind=8) zxi ! normalized displacement inside interpolation stencil

halo%n = n
halo%halo_width = halo_width
halo%ws = -1
halo%we = n+2
allocate(halo%w(-1:2,halo%ws:halo%we,halo_width))
allocate(halo%ind(halo%ws:halo%we,halo_width))


zda = 0.5_8*pi/n
do j=1, halo_width
    zb = 0.25_8*pi+(j-0.5_8)*zda
    do i = halo%ws, halo%we
        za = -0.25_8*pi+(i-0.5_8)*zda
        !cartesian coordinates of target cube face with alpha = za, beta = zb:
        zx = tan(za)
        zz = tan(zb)
        !implicitly project to sphere r = r/||r|| and then to source cube face zx,zy = zx/zz,zy/zz,zz=1
        zx = zx / zz
        !get alpha on source face:
        za = atan(zx)
        !i-index of source face point immediately to the left of projected target point
        halo%ind(i,j) = floor((za+0.25_8*pi-0.5_8*zda)/zda)+1
        zxi = (za+0.25_8*pi)/zda - (halo%ind(i,j)-0.5_8)
        halo%w(-1,i,j) =- zxi      *(zxi-1._8)*(zxi-2._8) / 6._8
        halo%w( 0,i,j) = (zxi+1._8)*(zxi-1._8)*(zxi-2._8) / 2._8
        halo%w( 1,i,j) =-(zxi+1._8)* zxi      *(zxi-2._8) / 2._8
        halo%w( 2,i,j) = (zxi+1._8)* zxi      *(zxi-1._8) / 6._8
    end do
end do

end function init_ecs_halo

subroutine ecs_ext_halo(this, f)
!interpolate source face values at halo zones to target face virtual points
use grid_function_mod, only: grid_function_t

class(ecs_halo_t),     intent(in)    :: this
type(grid_function_t), intent(inout) :: f

!local
integer(kind=4) n, hw
integer(kind=4) is, ie, js, je
integer(kind=4) nvi, nvj, nvk
integer(kind=4) isv, iev, jsv,jev, ksv, kev
integer(kind=4) klev
integer(kind=4) i,j,k
logical lhalo(4) !halo-procedure at edge
logical lcorn(4) !corner halo-procedure  for numeration of edges and corners see below
logical lfail_hw, lfail_corn, lfail_halo_long
real(kind=8) zf_csp(6,f%ks-f%nvk:f%ke+f%nvk,4) !values at corner-points and near-corner points@first halo-row
!local real(kind=8) zbufc(1-f%nvi:f%nvi,1-f%nvj:f%nvj,f%ks-f%nvk:f%ke+f%nvk) !store values for corner procedure

!short names for needed params
n = this%n
hw = this%halo_width

is = f%is;      ie = f%ie
js = f%js;      je = f%je
nvi = f%nvi;    nvj = f%nvj;     nvk = f%nvk

isv = is-nvi;   iev = ie+nvi
jsv = js-nvj;   jev = je+nvj
ksv = f%ks-nvk; kev = f%ke+nvk

klev = kev-ksv+1


!check if we need ecs edge-halo procedures for each specific tile edge
                                     !4__________3
lhalo(1) = js == 1                   !| edge2    |
lhalo(2) = je == n                   !|e        e|
lhalo(3) = is == 1                   !|d        d|
lhalo(4) = ie == n                   !|g        g|
lcorn(1) = lhalo(1) .and. lhalo(3)   !|3        4|
lcorn(2) = lhalo(1) .and. lhalo(4)   !1----------2
lcorn(3) = lhalo(4) .and. lhalo(2)   !   edge1
lcorn(4) = lhalo(2) .and. lhalo(3)

if(lcorn(1).or.lcorn(2).or.lcorn(3).or.lcorn(4)) then
    !interpolate f%p to special points in corners, store obtained values in zf_csp
    !zf_csp -> f%p after edge halo procedure
    call ecs_halo_corners(zf_csp,                                       & !special oints interpolated values
                          f%p,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                          this%w,this%ws,this%we,hw,                    & !weight array and its bounds
                          lcorn,n,hw)                                   !operating parameters
end if

if(lhalo(1).or.lhalo(2).or.lhalo(3).or.lhalo(4)) then
    call ecs_halo_edges(f%p,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                        this%w,this%ws,this%we,hw,this%ind,           & !weight array and its bounds
                        lhalo,n,hw)                                   !operating parameters
end if

call ecs_halo_corners_put(f%p,isv,iev,jsv,jev,klev,zf_csp,n,lcorn)

end subroutine ecs_ext_halo

!corner procedures for halo-zones
!currently interpolates only left/right values in first halo row
!planned: halo-corner interpolation
subroutine ecs_halo_corners(pfcsp,                                       &
                            pf,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                            w,ws,we,whw,                                 & !weight array and its bounds
                            lcorn,n,hw)                                    !global operating parameters

!output
real(kind=8),   intent(out) :: pfcsp(6,klev,4) !values at corner-points and near-corner points in first halo-row
!pfcsp placement for corner#1:
!  ha-|real points
!  lo |
!  __1|____
!   43|2 halo
!   65|
!  halo
!  corner
!input
!function values:
integer(kind=4), intent(in) :: is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj
real(kind=8),    intent(in) :: pf(isv:iev,jsv:jev,klev)
!interpolation weights
integer(kind=4), intent(in) :: ws,we,whw
real(kind=8),    intent(in) :: w(-1:2,ws:we,whw) !interpolation weights
!operting parameters
logical,         intent(in) :: lcorn(4)
integer(kind=4), intent(in) :: n, hw
!local
real(kind=8)    zbufc(1-nvi:nvi,1-nvj:nvj,klev) !store values for corner procedure
integer(kind=4) k, j, i, icor
logical lfail_corn
integer(kind=4) ist(4),jst(4),jdir(4),idir(4)

ist = [0,n+1,n+1,0  ]
idir= [1, -1, -1, 1]
jst = [0, 0, n+1,n+1]
jdir= [1, 1,  -1, -1]

!check if we have enough data for corner-halo-procedure
lfail_corn = (lcorn(1) .and. (isv> -2 .or. iev<  3 .or. jsv> -2 .or. jev<  3)) .or.& !corner 1
             (lcorn(2) .and. (isv>n-2 .or. iev<n+3 .or. jsv> -2 .or. jev<  3)) .or.& !corner 2
             (lcorn(3) .and. (isv>n-2 .or. iev<n+3 .or. jsv>n-2 .or. jev<n+3)) .or.& !corner 3
             (lcorn(4) .and. (isv> -2 .or. iev<  3 .or. jsv>n-2 .or. jev<n+3))       !corner 4
if(lfail_corn) call halo_avost("Must be at least 3 halo points for interpolation at corners")

do icor = 1,4
    if(lcorn(icor)) then
        do k=1, klev
            do j=1-nvj,nvj
                do i=1-nvi,nvi
                    zbufc(i,j,k) = pf(ist(icor)+idir(icor)*i,jst(icor)+jdir(icor)*j,k)
                end do
            end do
        end do
        call ecs_halo_1corner(pfcsp(:,:,icor),zbufc,nvi,nvj,klev,w,ws,we,whw)
    end if
end do

end subroutine ecs_halo_corners

subroutine ecs_halo_1corner(pfcsp,fc,nvi,nvj,klev,pw,ws,we,whw)
!output
real(kind=8),   intent(out) :: pfcsp(6,klev)
!input:
real(kind=8),    intent(in) :: fc(1-nvi:nvi,1-nvj:nvj,klev)
integer(kind=4), intent(in) :: nvi, nvj, klev
integer(kind=4), intent(in) :: ws, we, whw
real(kind=8),    intent(in) :: pw(-1:2,ws:we,whw)
!locals:
integer(kind=4) k
real(kind=8) zpw(-1:2)
real(kind=8) zx, zy, zz
real(kind=8) zw, zw2, zw3
real(kind=8) zfy1, zfy2, zff1, zff2, zff11, zff22

!pfcsp = 0._8

zpw= pw(:,1,1)
zw = pw(-1,1,1);   zw2 = zw*zw;   zw3 = zw2*zw

do k=1, klev
    zx = zpw(0)*fc(0,1,k)+zpw(1)*fc(0, 2,k)+zpw(2)*fc(0, 3,k)
    zy = zpw(0)*fc(1,0,k)+zpw(1)*fc(1,-1,k)+zpw(2)*fc(1,-2,k)
    zz = zpw(0)*fc(1,1,k)+zpw(1)*fc(2, 1,k)+zpw(2)*fc(3, 1,k)
    pfcsp(1,k) = (zx+zw*zy+zw2*zz)/(1-zw3)
    zfy1       = (zy+zw*zz+zw2*zx)/(1-zw3)

    zx = zpw(0)*fc(1,0,k)+zpw(1)*fc( 2,0,k)+zpw(2)*fc( 3,0,k)
    zy = zpw(0)*fc(0,1,k)+zpw(1)*fc(-1,1,k)+zpw(2)*fc(-2,1,k)
    zz = zpw(0)*fc(1,1,k)+zpw(1)*fc( 1,2,k)+zpw(2)*fc( 1,3,k)
    pfcsp(2,k) = (zx+zw*zy+zw2*zz)/(1-zw3)
    zfy2       = (zy+zw*zz+zw2*zx)/(1-zw3)

    zff1 = pw(-1,1,2)*fc(2,0,k)+pw(0,1,2)*fc(2,-1,k)+pw(1,1,2)*fc(2,-2,k)+pw(2,1,2)*fc(2,-3,k)
    zff2 = pw(-1,1,2)*fc(0,2,k)+pw(0,1,2)*fc(-1,2,k)+pw(1,1,2)*fc(-2,2,k)+pw(2,1,2)*fc(-3,2,k)
    pfcsp(3,k) = 0.5_8*(pw(-1,0,1)*zff1+pw(0,0,1)*zfy1+pw(1,0,1)*fc(0,1,k)+pw(2,0,1)*fc(0,2,k)+&
                        pw(-1,0,1)*zff2+pw(0,0,1)*zfy2+pw(1,0,1)*fc(1,0,k)+pw(2,0,1)*fc(2,0,k))

    zff1  = pw(-1,2,1)*fc(1,0,k)+pw(0,2,1)*fc(1,-1,k)+pw(1,2,1)*fc(1,-2,k)+pw(2,2,1)*fc(1,-3,k)
    zff11 = pw(-1,2,2)*fc(2,-1,k)+pw(0,2,2)*fc(2,-2,k)+pw(1,2,2)*fc(2,-3,k)+pw(2,2,2)*fc(2,-4,k)
    zff2  = pw(-1,2,1)*fc(0,1,k)+pw(0,2,1)*fc(-1,1,k)+pw(1,2,1)*fc(-2,1,k)+pw(2,2,1)*fc(-3,1,k)
    zff22 = pw(-1,2,2)*fc(-1,2,k)+pw(0,2,2)*fc(-2,2,k)+pw(1,2,2)*fc(-3,2,k)+pw(2,2,2)*fc(-4,2,k)

    pfcsp(4,k) = pw(-1,0,2)*zff1+pw(0,0,2)*fc(-1,1,k)+pw(1,0,2)*fc(-1,2,k)+pw(2,0,2)*fc(-1,3,k)
    pfcsp(5,k) = pw(-1,0,2)*zff2+pw(0,0,2)*fc(1,-1,k)+pw(1,0,2)*fc(2,-1,k)+pw(2,0,2)*fc(3,-1,k)

    pfcsp(6,k) = 0.5_8*(pw(-1,-1,2)*zff11+pw(0,-1,2)*zff1+pw(1,-1,2)*fc(-1,1,k)+pw(2,-1,2)*fc(-1,2,k)+&
                        pw(-1,-1,2)*zff22+pw(0,-1,2)*zff2+pw(1,-1,2)*fc(1,-1,k)+pw(2,-1,2)*fc(2,-1,k))
    
end do

end subroutine ecs_halo_1corner

subroutine ecs_halo_corners_put(pf,isv,iev,jsv,jev,klev,pfcsp,n,lcorn)
!output
real(kind=8), intent(inout) :: pf(isv:iev,jsv:jev,klev) !output field
!input
integer(kind=4), intent(in) :: isv,iev,jsv,jev,klev   !dimensions of pf
real(kind=8),    intent(in) :: pfcsp(6,klev,4)          !previously interpolated values at 'special' corner points
integer(kind=4), intent(in) :: n
logical,         intent(in) :: lcorn(4)                 !if pf corner is cubedsphere corner

if(lcorn(1)) then
    pf( 0, 1,:) = pfcsp(1,:,1)
    pf( 1, 0,:) = pfcsp(2,:,1)
    pf( 0, 0,:) = pfcsp(3,:,1)
    pf(-1, 0,:) = pfcsp(4,:,1)
    pf( 0,-1,:) = pfcsp(5,:,1)
    pf(-1,-1,:) = pfcsp(6,:,1)
end if
if(lcorn(2)) then
    pf(n+1, 1,:) = pfcsp(1,:,2)
    pf(n  , 0,:) = pfcsp(2,:,2)
    pf(n+1, 0,:) = pfcsp(3,:,2)
    pf(n+2, 0,:) = pfcsp(4,:,2)
    pf(n+1,-1,:) = pfcsp(5,:,2)
    pf(n+2,-1,:) = pfcsp(6,:,2)
end if
if(lcorn(3)) then
    pf(n+1,n  ,:) = pfcsp(1,:,3)
    pf(n  ,n+1,:) = pfcsp(2,:,3)
    pf(n+1,n+1,:) = pfcsp(3,:,3)
    pf(n+2,n+1,:) = pfcsp(4,:,3)
    pf(n+1,n+2,:) = pfcsp(5,:,3)
    pf(n+2,n+2,:) = pfcsp(6,:,3)
end if
if(lcorn(4)) then
    pf( 0,n  ,:) = pfcsp(1,:,4)
    pf( 1,n+1,:) = pfcsp(2,:,4)
    pf( 0,n+1,:) = pfcsp(3,:,4)
    pf(-1,n+1,:) = pfcsp(4,:,4)
    pf( 0,n+2,:) = pfcsp(5,:,4)
    pf(-1,n+2,:) = pfcsp(6,:,4)
end if
end subroutine ecs_halo_corners_put

subroutine ecs_halo_edges(pf,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                          w,ws,we,whw,ind,                             & !weights and indices  arrays and their bounds
                          lhalo,n,hw)                                    !global operating parameters

!in-output
real(kind=8), intent(inout) :: pf(isv:iev,jsv:jev,klev)
!input
integer(kind=4), intent(in) :: is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj
!interpolation weights
integer(kind=4), intent(in) :: ws,we,whw
real(kind=8),    intent(in) :: w(-1:2,ws:we,whw) !interpolation weights
integer(kind=4), intent(in) :: ind(ws:we,whw)
!operting parameters
logical,         intent(in) :: lhalo(4)
integer(kind=4), intent(in) :: n, hw
!local
integer(kind=4) k, j, i
integer(kind=4) ish, ieh, jsh, jeh
logical lfail_hw
logical lfail_halo_long
real(kind=8) zbufy(jsv:jev, hw, klev)
real(kind=8) zbufx(isv:iev, hw, klev)

!Internal checks to avoid (at least some part of) out-of-array errors
!check if we have enough data for edge halo-procedure
lfail_hw = ((lhalo(1) .or. lhalo(2)) .and. nvi< hw) .or. ((lhalo(3) .or. lhalo(4)) .and. nvj< hw)
if(lfail_hw) call halo_avost("cubed sphere halo_width > grid_function_t halo width, can't continue")
!check if we have enough data along edges to perform interpolations
ish = max(1,is-hw); ieh = min(n,ie+hw)
jsh = max(1,js-hw); jeh = min(n,je+hw)
lfail_halo_long = (lhalo(1) .and. (minval(ind(ish,:))-1<isv .or. maxval(ind(ieh,:))+2>iev)) .or.& !edge 1
                  (lhalo(2) .and. (minval(ind(ish,:))-1<isv .or. maxval(ind(ieh,:))+2>iev)) .or.& !edge 2
                  (lhalo(3) .and. (minval(ind(jsh,:))-1<jsv .or. maxval(ind(jeh,:))+2>jev)) .or.& !edge 3
                  (lhalo(4) .and. (minval(ind(jsh,:))-1<jsv .or. maxval(ind(jeh,:))+2>jev))       !edge 4
if(lfail_halo_long) call halo_avost("not enough points along cub.sph edge to perform halo-interpolations (nvi=5 needed)")

!calculate values at corner special points 
!store them in sepparate arrays to not to spoil original f%p

!Halo-procedures along edges
if(lhalo(1)) then
    zbufx = pf(isv:iev, 0:1-hw:-1,:)
    call ecs_ext_halo_1e(zbufx,isv,iev,is,ie,hw,klev,w,ws,we,whw,ind,n)
    pf(isv:iev,1-hw:0,:) = zbufx(:,hw:1:-1,:)
end if
if(lhalo(2)) then
    zbufx = pf(isv:iev, n+1:n+hw, :)
    call ecs_ext_halo_1e(zbufx,isv,iev,is,ie,hw,klev,w,ws,we,whw,ind,n)
    pf(isv:iev, n+1:n+hw, :) = zbufx
end if
if(lhalo(3)) then
    do k=1,klev
        do j=jsv, jev
            zbufy(j,1:hw,k) = pf(0:1-hw:-1,j,k)
        end do
    end do
    call ecs_ext_halo_1e(zbufy,jsv,jev,js,je,hw,klev,w,ws,we,whw,ind,n)
    do k=1,klev
        do j=1,hw
            pf(1-j,jsv:jev,k) = zbufy(jsv:jev,j,k)
        end do
    end do
end if
if(lhalo(4)) then
    do k=1, klev
        do j=jsv,jev
            zbufy(j,1:hw,k) = pf(n+1:n+hw,j,k)
        end do
    end do
    call ecs_ext_halo_1e(zbufy,jsv,jev,js,je,hw,klev,w,ws,we,whw,ind,n)
    do k=1,klev
        do j=1,hw
            pf(n+j,jsv:jev,k) = zbufy(jsv:jev,j,k)
        end do
    end do
end if


end subroutine ecs_halo_edges

subroutine ecs_ext_halo_1e(zf,i1v,i2v,i1,i2,hw,klev,w,ws,we,whw,ind,n)
!input
integer(kind=4), intent(in) :: i1v,i2v,i1,i2,hw,klev
!in-output:
real(kind=8), intent(inout) :: zf(i1v:i2v,hw,klev)!input: source face values, output: interpolated target face values
!input
integer(kind=4), intent(in) :: ws,we,whw
real(kind=8),    intent(in) :: w(-1:2,ws:we,whw)
integer(kind=4), intent(in) :: ind(ws:we,whw)
integer(kind=4), intent(in) :: n
!locals
real(kind=8) zh(i1-hw:i2+hw)!buffer for interpolated values
integer i, j, k, ii
integer ihs(hw), ihe(hw)

!exclude "special" points at first row, where we initially have not enough surrounding values to use cubic-interpolation
ihs(2:) = 1; ihs(1) = 2
ihe(2:) = n; ihe(1) = n-1

do k=1, klev
    do j=1, hw
        do i = max(i1-hw,ihs(j)),min(i2+hw,ihe(j))
            ii = ind(i,j)
            zh(i) = sum(w(:,i,j)*zf(ii-1:ii+2,j,k))
        end do
        zf(max(i1-hw,ihs(j)):min(i2+hw,ihe(j)),j,k) = zh(max(i1-hw,ihs(j)):min(i2+hw,ihe(j)))
    end do
end do

end subroutine ecs_ext_halo_1e

subroutine halo_avost(str)
use mpi
integer ierr
character(*) str
print *, str
!print '(6(A,i8,1x))', "is=", is, "ie=", ie, "js=", js, "je=", je, "nvi=", iev-ie, "nvj=", jev-je
print *, "exit"
call mpi_finalize(ierr)
stop
end subroutine halo_avost


subroutine ecs_ext_halo_vect(this, u, v)
!interpolate source face values at halo zones to target face virtual points
use grid_function_mod, only: grid_function_t

class(ecs_halo_t),     intent(in)    :: this
type(grid_function_t), intent(inout) :: u, v

end subroutine ecs_ext_halo_vect

end module ecs_halo_mod
