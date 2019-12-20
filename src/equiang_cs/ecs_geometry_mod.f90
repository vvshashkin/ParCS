!Module for equiangular cubed sphere geometrical parameters
module ecs_geometry_mod
implicit none
!horizontal dim
integer nx, halo_w
!angle step size:
real(kind=8) dxa
!normalized cartesian coordinates of cubed sphere grid points (2d):
!indices x, y, face
real(kind=8), allocatable :: rhx(:,:,:), rhy(:,:,:), rhz(:,:,:) !scalar points (cell centers)

contains
subroutine ecs_geometry_mod_init(nh, halo_width)
use const_mod,    only : pi
use topology_mod, only : ex, ey, n
!input:
integer nh         !number of points in one horizontal direction (panel size is nh*nh)
integer halo_width !halo zones width

!locals:
real(kind=8) za, zb     !cubed sphere panel angles
real(kind=8) zx, zy, zz !grid point cartesian coord at prototype spherical-cube face (with ex=[1,0,0], ey=[0,1,0], n=[0,0,1])
real(kind=8) zr         !cube-point radius vector length
integer i, j, ifc

nx = nh
halo_w = halo_width
dxa = 0.5_8*pi/(1._8*nh)

allocate(rhx(1-halo_width:nh+halo_width,1-halo_width:nh+halo_width,6))
allocate(rhy(1-halo_width:nh+halo_width,1-halo_width:nh+halo_width,6))
allocate(rhz(1-halo_width:nh+halo_width,1-halo_width:nh+halo_width,6))

do ifc = 1,6
    do j  =1-halo_width,nh+halo_width
        zb = -0.25_8*pi+(j-0.5_8)*dxa
        do i  =1-halo_width,nh+halo_width
            za = -0.25_8*pi+(i-0.5_8)*dxa
            !grid point at cube face: (assumes zx = tan(za), zy = tan(zb), zz = 1)
            zr = sqrt(1._8+tan(za)**2+tan(zb)**2)
            !coordinates at prototype spherical face: vec(r)/||r||, r = (zx, zy, zz)
            zx = tan(za)/zr; zy = tan(zb)/zr; zz = 1._8/zr
            !transform to real spherical face
            rhx(i,j,ifc) = ex(1,ifc)*zx+ey(1,ifc)*zy-n(1,ifc)*zz
            rhy(i,j,ifc) = ex(2,ifc)*zx+ey(2,ifc)*zy-n(2,ifc)*zz
            rhz(i,j,ifc) = ex(3,ifc)*zx+ey(3,ifc)*zy-n(3,ifc)*zz
        end do
    end do
end do

end subroutine ecs_geometry_mod_init

subroutine ecs_geometry_mod_check()
use topology_mod, only : ex, ey, n
logical lsym_check
integer(kind=4) ifc, is, ie, js, je
real(kind=8) ee(3)

lsym_check = .true.

is = 1-halo_w; ie = nx+halo_w
js = 1-halo_w; je = nx+halo_w
do ifc = 1,6
    ee = real(ex(:,ifc),8)
    lsym_check = lsym_check .and. symmetric(is,ie,js,js,ifc,ee)
    lsym_check = lsym_check .and. symmetric(is,ie,je,je,ifc,ee)
    ee = real(ey(:,ifc),8)
    lsym_check = lsym_check .and. symmetric(is,is,js,je,ifc,ee)
    lsym_check = lsym_check .and. symmetric(ie,ie,js,je,ifc,ee)
end do
if(lsym_check) then
    print *, "x,y,z symmetricity check passed"
else
    print *, "x,y,z symmetricity check failed"
end if
contains
    logical function symmetric(i1,i2,j1,j2,ifc,ee) result(lsym)
        integer(kind=4) i1,i2,j1,j2,ifc
        real(kind=8) ee(3)
        real(kind=8) a(3), b(3), zpr
        real(kind=8), parameter :: zeps = 1e-16_8
    
        a = [rhx(i1,j1,ifc), rhy(i1,j1,ifc), rhz(i1,j1,ifc)]
        b = [rhx(i2,j2,ifc), rhy(i2,j2,ifc), rhz(i2,j2,ifc)]
        zpr = sum(a*ee)
        a = a - 2._8*zpr*ee
        lsym = sum(abs(a-b))<zeps
    end function symmetric

end subroutine ecs_geometry_mod_check

end module ecs_geometry_mod
