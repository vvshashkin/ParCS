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

implicit none

!type-container for precomputed interpolation weights to have
!multiple halo-interpolators for multiple grids (for geometric MG-solver)

type ecs_halo_t
  integer n          !number of real grid points along cubed-sphere face edge
  integer halo_width !number of rows in halo-zone
  real(kind=8)   , allocatable :: w(:,:,:)
  integer(kind=4), allocatable :: ind(:,:)
  contains
  procedure, public :: ext_halo=>ecs_ext_halo
end type ecs_halo_t

type(ecs_halo_t), allocatable :: whalo(:)

contains

subroutine ecs_halo_mod_init(nn,halo_width)
!input
integer nn(:)      !dimension of each cub-sph grid (6*nn(i)**2)
integer halo_width !number of rows in halo-zone
!locals:
integer nmg   !number of grids
integer i

nmg = size(nn)

allocate(whalo(nmg))

do i=1,nmg
    call init_ecs_halo_t(nn(i),halo_width,whalo(i))
end do

end subroutine ecs_halo_mod_init

subroutine init_ecs_halo_t(n,halo_width,wh)
use const_mod, only: pi

integer n, halo_width
type(ecs_halo_t) wh

!local:
integer i,j
real(kind=8) zda !angle grid-step
real(kind=8) za, zb !angle coordinates at target grid face
real(kind=8) zx, zz !cartesian coordinates, zy is not needed
real(kind=8) zxi ! normalized displacement inside interpolation stencil

wh%n = n
wh%halo_width = halo_width
allocate(wh%w(-1:2,n,halo_width))
allocate(wh%ind(n,halo_width))


zda = 0.5_8*pi/n
do j=1, halo_width
    zb = 0.25_8*pi+(j-0.5_8)*zda
    do i = 1, n
        za = -0.25_8*pi+(i-0.5_8)*zda
        !cartesian coordinates of target cube face with alpha = za, beta = zb:
        zx = tan(za)
        zz = tan(zb)
        !implicitly project to sphere r = r/||r|| and then to source cube face zx,zy = zx/zz,zy/zz,zz=1
        zx = zx / zz
        !get alpha on source face:
        za = atan(zx)
        !i-index of source face point immediately to the left of projected target point
        wh%ind(i,j) = floor((za+0.25_8*pi-0.5_8*zda)/zda)+1
        zxi = (za+0.25_8*pi)/zda - (wh%ind(i,j)-0.5_8)
        wh%w(-1,i,j) =- zxi      *(zxi-1._8)*(zxi-2._8) / 6._8
        wh%w( 0,i,j) = (zxi+1._8)*(zxi-1._8)*(zxi-2._8) / 2._8
        wh%w( 1,i,j) =-(zxi+1._8)* zxi      *(zxi-2._8) / 2._8
        wh%w( 2,i,j) = (zxi+1._8)* zxi      *(zxi-1._8) / 6._8
    end do
end do

end subroutine init_ecs_halo_t

subroutine ecs_ext_halo(this, f)
!interpolate source face values at halo zones to target face virtual points
use grid_function_mod, only: grid_function_t

class(ecs_halo_t) this
type(grid_function_t) f

!local
real(kind=8) zbufy(f%js-f%nvj : f%je+f%nvj, this%halo_width, f%ks-f%nvk:f%ke+f%nvk)
real(kind=8) zbufx(f%is-f%nvi : f%ie+f%nvi, this%halo_width, f%ks-f%nvk:f%ke+f%nvk)
integer(kind=4) n, hw, klev
integer(kind=4) is,ie,js,je,isv,iev,jsv,jev, nvi, nvj
integer(kind=8) i,j,k
logical lhalo(4) !halo-procedurea at edge
logical lcorn(4) !corner halo-procedure  for numeration of edges and corners see below
logical lfail_hw, lfail_corn, lfail_halo_long
real(kind=8) zf_csp(2,f%ks-f%nvk:f%ke+f%nvk,4) !values at special points in corners
real(kind=8) zbufc(1-f%nvi:f%nvi,1-f%nvj:f%nvj,f%ks-f%nvk:f%ke+f%nvk) !store values for corner procedure
!short names for needed params
n = this%n; hw = this%halo_width; klev = f%ke-f%ks+2*f%nvk+1
is = f%is; ie = f%ie; js = f%js; je = f%je
nvi = f%nvi; nvj = f%nvj
isv = is-nvi; iev = ie+nvi
jsv = js-nvj; jev = je+nvj

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

!Internal checks to avoid (at least some part of) out-of-array errors
!check if we have enough data for edge halo-procedure
lfail_hw = ((lhalo(1) .or. lhalo(2)) .and. nvi< hw) .or. ((lhalo(3) .or. lhalo(4)) .and. nvj< hw)
if(lfail_hw) call halo_avost("cubed sphere halo_width > grid_function_t halo width, can't continue")
!check if we have enough data for corner-halo-procedure
lfail_corn = (lcorn(1) .and. (isv> -2 .or. iev<  3 .or. jsv> -2 .or. jev<  3)) .or.& !corner 1
             (lcorn(2) .and. (isv>n-2 .or. iev<n+3 .or. jsv> -2 .or. jev<  3)) .or.& !corner 2
             (lcorn(3) .and. (isv>n-2 .or. iev<n+3 .or. jsv>n-2 .or. jev<n+3)) .or.& !corner 3
             (lcorn(4) .and. (isv> -2 .or. iev<  3 .or. jsv>n-2 .or. jev<n+3))       !corner 4
if(lfail_corn) call halo_avost("Must be at least 3 halo points for interpolation at corners")
!check if we have enough data along edges to perform interpolations
lfail_halo_long = (lhalo(1) .and. (minval(this%ind(is,:))-1<isv .or. maxval(this%ind(ie,:))+2>iev)) .or.& !edge 1
                  (lhalo(2) .and. (minval(this%ind(is,:))-1<isv .or. maxval(this%ind(ie,:))+2>iev)) .or.& !edge 2
                  (lhalo(3) .and. (minval(this%ind(js,:))-1<jsv .or. maxval(this%ind(je,:))+2>jev)) .or.& !edge 3
                  (lhalo(4) .and. (minval(this%ind(js,:))-1<jsv .or. maxval(this%ind(je,:))+2>jev))       !edge 4
if(lfail_halo_long) call halo_avost("not enough points along cub.sph edge to perform halo-interpolations (nvi=5 needed)")

!calculate values at corner special points 
!store them in sepparate arrays to not to spoil original f%p
if(lcorn(1)) then
    do k=f%ks-f%nvk, f%ke+f%nvk
        do j=1-nvj,nvj
            do i=1-nvi,nvi
                zbufc(i,j,k) = f%p(i,j,k)
            end do
        end do
    end do
    call ecs_corner_halo(zf_csp(:,:,1),zbufc)
end if
if(lcorn(2)) then
    do k=f%ks-f%nvk, f%ke+f%nvk
        do j=1-nvj,nvj
            do i=1-nvi,nvi
                zbufc(i,j,k) = f%p(n-i+1,j,k)
            end do
        end do
    end do
    call ecs_corner_halo(zf_csp(:,:,2),zbufc)
end if
if(lcorn(3)) then
    do k=f%ks-f%nvk, f%ke+f%nvk
        do j=1-nvj,nvj
            do i=1-nvi,nvi
                zbufc(i,j,k) = f%p(n-i+1,n-j+1,k)
            end do
        end do
    end do
    call ecs_corner_halo(zf_csp(:,:,3),zbufc)
end if
if(lcorn(4)) then
    do k=f%ks-f%nvk, f%ke+f%nvk
        do j=1-nvj,nvj
            do i=1-nvi,nvi
                zbufc(i,j,k) = f%p(i,n-j+1,k)
            end do
        end do
    end do
    call ecs_corner_halo(zf_csp(:,:,4),zbufc)
end if

!Halo-procedures along edges
if(lhalo(1)) then
    zbufx = f%p(isv:iev, 0:1-hw:-1,:)
    call ecs_ext_halo_1e(isv,iev,is,ie,hw,zbufx)
    f%p(isv:iev,1-hw:0,:) = zbufx(:,hw:1:-1,:)
end if
if(lhalo(2)) then
    zbufx = f%p(isv:iev, n+1:n+hw, :)
    call ecs_ext_halo_1e(isv,iev,is,ie,hw,zbufx)
    f%p(isv:iev, n+1:n+hw, :) = zbufx
end if
if(lhalo(3)) then
    do k=f%ks-f%nvk, f%ke+f%nvk
        do j=jsv, jev
            zbufy(j,1:hw,k) = f%p(0:1-hw:-1,j,k)
        end do
    end do
    call ecs_ext_halo_1e(jsv,jev,js,je,hw,zbufy)
    do k=f%ks-f%nvk, f%ke+f%nvk
        do j=1,hw
            f%p(1-j,jsv:jev,k) = zbufy(jsv:jev,j,k)
        end do
    end do
end if
if(lhalo(4)) then
    do k=f%ks-f%nvk, f%ke+f%nvk
        do j=jsv,jev
            zbufy(j,1:hw,k) = f%p(n+1:n+hw,j,k)
        end do
    end do
    call ecs_ext_halo_1e(jsv,jev,js,je,hw,zbufy)
    do k=f%ks-f%nvk, f%ke+f%nvk
        do j=1,hw
            f%p(n+j,jsv:jev,k) = zbufy(jsv:jev,j,k)
        end do
    end do
end if

if(lcorn(1)) then
    f%p(0,1,:) = zf_csp(1,:,1)
    f%p(1,0,:) = zf_csp(2,:,1)
end if
if(lcorn(2)) then
    f%p(n+1,1,:) = zf_csp(1,:,2)
    f%p(n,0,:) = zf_csp(2,:,2)
end if
if(lcorn(3)) then
    f%p(n+1,n,:) = zf_csp(1,:,3)
    f%p(n,n+1,:) = zf_csp(2,:,3)
end if
if(lcorn(4)) then
    f%p(0,n,:) = zf_csp(1,:,4)
    f%p(1,n+1,:) = zf_csp(2,:,4)
end if

contains

subroutine ecs_corner_halo(fcsp,fc)
real(kind=8) fcsp(2,klev)
real(kind=8) fc(1-nvi:nvi,1-nvj:nvj,klev)
integer k
real(kind=8) zx, zy, zz

do k=1, klev
    zx = this%w(0,1,1)*fc(0,1,k)+this%w(1,1,1)*fc(0, 2,k)+this%w(2,1,1)*fc(0, 3,k)
    zy = this%w(0,1,1)*fc(1,0,k)+this%w(1,1,1)*fc(1,-1,k)+this%w(2,1,1)*fc(1,-2,k)
    zz = this%w(0,1,1)*fc(1,1,k)+this%w(1,1,1)*fc(2, 1,k)+this%w(2,1,1)*fc(3, 1,k)
    fcsp(1,k) = (zx+this%w(-1,1,1)*zy+this%w(-1,1,1)**2*zz)/(1-this%w(-1,1,1)**3)

    zx = this%w(0,1,1)*fc(1,0,k)+this%w(1,1,1)*fc( 2,0,k)+this%w(2,1,1)*fc( 3,0,k)
    zy = this%w(0,1,1)*fc(0,1,k)+this%w(1,1,1)*fc(-1,1,k)+this%w(2,1,1)*fc(-2,1,k)
    zz = this%w(0,1,1)*fc(1,1,k)+this%w(1,1,1)*fc( 1,2,k)+this%w(2,1,1)*fc( 1,3,k)
    fcsp(2,k) = (zx+this%w(-1,1,1)*zy+this%w(-1,1,1)**2*zz)/(1-this%w(-1,1,1)**3)
end do

end subroutine ecs_corner_halo

subroutine ecs_ext_halo_1e(i1v,i2v,i1,i2,hw,zf)
integer(kind=4) i1v,i2v,i1,i2,hw
real(kind=8) zf(i1v:i2v,hw,klev)!input: source face values, output: interpolated target face values
!locals
real(kind=8) zh(max(1,i1):min(n,i2))!buffer for interpolated values
integer i, j, k, ii
integer ihs(hw), ihe(hw)

!exclude "special" points at first row, where we initially have not enough surrounding values to use cubic-interpolation
ihs(2:) = 1; ihs(1) = 2
ihe(2:) = n; ihe(1) = n-1

do k=1, klev
    do j=1, hw
        do i = max(i1,ihs(j)),min(i2,ihe(j))
            ii = this%ind(i,j)
            zh(i) = sum(this%w(:,i,j)*zf(ii-1:ii+2,j,k))
        end do
        zf(max(i1,ihs(j)):min(i2,ihe(j)),j,k) = zh(max(i1,ihs(j)):min(i2,ihe(j)))
    end do
end do
end subroutine ecs_ext_halo_1e

subroutine halo_avost(str)
use mpi
integer ierr
character(*) str
print *, str
print '(6(A,i8,1x))', "is=", is, "ie=", ie, "js=", js, "je=", je, "nvi=", iev-ie, "nvj=", jev-je
print *, "exit"
call mpi_finalize(ierr)
stop
end subroutine halo_avost

end subroutine ecs_ext_halo

end module ecs_halo_mod