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
integer(kind=4) n, hw, j, klev, k
integer(kind=4) is,ie,js,je,isv,iev,jsv,jev
logical lhalo(4)
logical lfail_hw

n = this%n; hw = this%halo_width; klev = f%ke-f%ks+2*f%nvk+1
is = f%is; ie = f%ie; js = f%js; je = f%je
isv = is-f%nvi; iev = ie+f%nvi
jsv = js-f%nvj; jev = je+f%nvj
!check if we need ecs edge-halo procedures for each specific tile edge
                     ! __________
lhalo(1) = f%js == 1 !| edge2    |
lhalo(2) = f%je == n !|e        e|
lhalo(3) = f%is == 1 !|d        d|
lhalo(4) = f%ie == n !|g        g|
                     !|3        4|
                     ! ----------
                     !   edge1

!check if we have enough data
lfail_hw = ((lhalo(1) .or. lhalo(2)) .and. f%nvi< hw) .or. ((lhalo(3) .or. lhalo(4)) .and. f%nvj< hw)
if(lfail_hw) call halo_width_avost()

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

contains
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

subroutine halo_width_avost()
use mpi
integer ierr
print *, "cubed sphere halo_width > grid_function_t halo width, can't continue"
call mpi_finalize(ierr)
end subroutine halo_width_avost

end subroutine ecs_ext_halo

end module ecs_halo_mod
