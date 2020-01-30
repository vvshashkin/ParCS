!Precompute information to interpolate values to virtual points beyond face edge
!Here we use two prototype face (all other cases for scalars are symmetric):
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
module ecs_halo_factory_mod
use ecs_halo_mod, only : ecs_halo_t

implicit none

private
public   :: init_ecs_halo

contains

type(ecs_halo_t) function init_ecs_halo(is,ie,js,je,nx,halo_width,hx) result(halo)
use mesh_mod,  only : mesh_t

integer(kind=4), intent(in) :: is,ie,js,je,nx
integer(kind=4), intent(in) :: halo_width
real(kind=8),    intent(in) :: hx

halo%n = nx
halo%halo_width = halo_width

!check if we need ecs edge-halo procedures for each specific tile edge
                                                    !4__________3
halo%lhalo(1) = js == 1                             !| edge2    |
halo%lhalo(2) = je == nx                            !|e        e|
halo%lhalo(3) = is == 1                             !|d        d|
halo%lhalo(4) = ie == nx                            !|g        g|
halo%lcorn(1) = halo%lhalo(1) .and. halo%lhalo(3)   !|3        4|
halo%lcorn(2) = halo%lhalo(1) .and. halo%lhalo(4)   !1----------2
halo%lcorn(3) = halo%lhalo(4) .and. halo%lhalo(2)   !   edge1
halo%lcorn(4) = halo%lhalo(2) .and. halo%lhalo(3)

!initialize interpolation weights
if(halo%lhalo(1) .or. halo%lhalo(2)) then
    halo%wsx = max(is-halo_width,-1)
    halo%wex = min(ie+halo_width,nx+2)
    if(halo%lcorn(1).or.halo%lcorn(4)) halo%wsx = min(halo%wsx,-1)        !approve that we will have weights needed for corners
    if(halo%lcorn(2).or.halo%lcorn(3)) halo%wex = max(halo%wex,nx+2) !x-y and nx/2 symetri will be used
    call init_halo_interp(halo%wx, halo%indx, hx, halo%wsx, halo%wex, halo_width)
end if
if(halo%lhalo(3) .or. halo%lhalo(4)) then
    halo%wsy = max(js-halo_width,-1)
    halo%wey = min(je+halo_width,nx+2)
    call init_halo_interp(halo%wy, halo%indy, hx, halo%wsy, halo%wey, halo_width)
end if

end function init_ecs_halo

subroutine init_halo_interp(w, ind, dxa, i1, i2, hw)
use const_mod, only : pi
!output
real(kind=8),    allocatable, intent(out) :: w(:,:,:)
integer(kind=4), allocatable, intent(out) :: ind(:,:)
!input
real(kind=8),    intent(in) :: dxa
integer(kind=4), intent(in) :: i1, i2, hw
!local:
integer i,j
real(kind=8) za, zb !angle coordinates at target grid face
real(kind=8) zx, zz !cartesian coordinates, zy is not needed
real(kind=8) zxi ! normalized displacement inside interpolation stencil

allocate(w(-1:2,i1:i2,hw))
allocate(ind(i1:i2,hw))

do j=1, hw
    zb = 0.25_8*pi+(j-0.5_8)*dxa
    do i = i1, i2
        za = -0.25_8*pi+(i-0.5_8)*dxa
        !cartesian coordinates of target cube face with alpha = za, beta = zb:
        zx = tan(za)
        zz = tan(zb)
        !implicitly project to sphere r = r/||r|| and then to source cube face zx,zy = zx/zz,zy/zz,zz=1
        zx = zx / zz
        !get alpha on source face:
        za = atan(zx)
        !i-index of source face point immediately to the left of projected target point
        ind(i,j) = floor((za+0.25_8*pi-0.5_8*dxa)/dxa)+1
        zxi = (za+0.25_8*pi)/dxa - (ind(i,j)-0.5_8)
        w(-1,i,j) =- zxi      *(zxi-1._8)*(zxi-2._8) / 6._8
        w( 0,i,j) = (zxi+1._8)*(zxi-1._8)*(zxi-2._8) / 2._8
        w( 1,i,j) =-(zxi+1._8)* zxi      *(zxi-2._8) / 2._8
        w( 2,i,j) = (zxi+1._8)* zxi      *(zxi-1._8) / 6._8
    end do
end do

end subroutine init_halo_interp

end module ecs_halo_factory_mod