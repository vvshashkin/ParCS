!Module for equiangular cubed sphere geometrical parameters
module ecs_geometry_mod
implicit none
!!horizontal dim
!integer nx, halo_w
!!angle step size:
!real(kind=8) dxa
!!normalized cartesian coordinates of cubed sphere grid points (2d):
!!indices x, y, face
!real(kind=8), allocatable :: rhx(:,:,:), rhy(:,:,:), rhz(:,:,:) !scalar points (cell centers)
!real(kind=8), allocatable :: acov(:,:,:,:), bcov(:,:,:,:) !covariant vectors
!real(kind=8), allocatable :: actv(:,:,:,:), bctv(:,:,:,:) !contravariant vectors
!
contains
!subroutine ecs_geometry_mod_init(nh, halo_width)
!use const_mod,    only : pi
!!input:
!integer nh         !number of points in one horizontal direction (panel size is nh*nh)
!integer halo_width !halo zones width
!
!!locals:
!real(kind=8)    za, zb     !cubed sphere panel angles
!real(kind=8)    zx, zy, zz !grid point cartesian coord at prototype spherical-cube face (with ex=[1,0,0], ey=[0,1,0], n=[0,0,1])
!real(kind=8)    zr         !cube-point radius vector length
!real(kind=8)    zav(3), zbv(3)
!integer(kind=4) i, j, ifc
!
!nx = nh
!halo_w = halo_width
!dxa = 0.5_8*pi/(1._8*nh)
!
!allocate(rhx(1-halo_width:nh+halo_width,1-halo_width:nh+halo_width,6))
!allocate(rhy(1-halo_width:nh+halo_width,1-halo_width:nh+halo_width,6))
!allocate(rhz(1-halo_width:nh+halo_width,1-halo_width:nh+halo_width,6))
!allocate(acov(3,1-halo_width:nh+halo_width,1-halo_width:nh+halo_width,6))
!allocate(bcov(3,1-halo_width:nh+halo_width,1-halo_width:nh+halo_width,6))
!allocate(actv(3,1-halo_width:nh+halo_width,1-halo_width:nh+halo_width,6))
!allocate(bctv(3,1-halo_width:nh+halo_width,1-halo_width:nh+halo_width,6))
!
!do ifc = 1,6
!    do j  =1-halo_width,nh+halo_width
!        zb = -0.25_8*pi+(j-0.5_8)*dxa
!        do i  =1-halo_width,nh+halo_width
!            za = -0.25_8*pi+(i-0.5_8)*dxa
!            call ecs_ab2xyz_proto(zx,zy,zz,za,zb)
!            call ecs_proto2face(rhx(i,j,ifc),rhy(i,j,ifc),rhz(i,j,ifc),zx,zy,zz,ifc)
!            call ecs_cov_proto(zav,zbv,za,zb)
!            call ecs_proto2face(acov(1,i,j,ifc),acov(2,i,j,ifc),acov(3,i,j,ifc),zav(1),zav(2),zav(3),ifc)
!            call ecs_proto2face(bcov(1,i,j,ifc),bcov(2,i,j,ifc),bcov(3,i,j,ifc),zbv(1),zbv(2),zbv(3),ifc)
!        end do
!    end do
!end do
!
!end subroutine ecs_geometry_mod_init

function ecs_ab2xyz_proto(alpha,beta) result(xyz)
    real(kind=8)                   xyz(3)
    real(kind=8), intent(in)    :: alpha, beta
    !locals:
    real(kind=8) r
    integer(kind=4) i

    !grid point at proto cube face: (assumes x = tan(alpha), y = tan(beta), z = 1)
    r = sqrt(1._8+tan(alpha)**2+tan(beta)**2)
    !coordinates at prototype spherical face: vec(r)/||r||, r = (x, y, z)
    xyz(1) = tan(alpha)/r; xyz(2) = tan(beta)/r; xyz(3) = 1._8/r

end function ecs_ab2xyz_proto

subroutine ecs_cov_proto(pacov,pbcov,pa,pb)
!covariant vectors at prototype face given alpha&beta
real(kind=8), intent(out) :: pacov(3), pbcov(3)
real(kind=8), intent(in)  :: pa, pb
!local
real(kind=8) zta, ztb

zta = tan(pa);   ztb = tan(pb)
!d (xyz)^T/ d alpha
pacov(1) = (1d0+zta**2d0)*(1d0+ztb**2d0)/(1d0+zta**2d0+ztb**2d0)**(3d0/2d0)
pacov(2) = -zta*ztb*(1d0+zta**2d0)/(1d0+zta**2d0+ztb**2d0)**(3d0/2d0)
pacov(3) = -zta*(1d0+zta**2d0) / (1d0+zta**2d0+ztb**2d0)**(3d0/2d0)
!d (xyz)^T/ d beta
pbcov(1) = -zta*ztb*(1d0+ztb**2d0)/(1d0+zta**2d0+ztb**2d0)**(3d0/2d0)
pbcov(2) = (1d0+zta**2d0)*(1d0+ztb**2d0)/(1d0+zta**2d0+ztb**2d0)**(3d0/2d0)
pbcov(3) = -ztb*(1d0+ztb**2d0) / (1d0+zta**2d0+ztb**2d0)**(3d0/2d0)
end subroutine ecs_cov_proto

function ecs_proto2face(xyz,panel_ind) result(xyz1)
!transform prototype-face xyz to real spherical face
use topology_mod, only : ex, ey, n
real(kind=8)                :: xyz1(3)
real(kind=8), intent(in)    :: xyz(3)
integer(kind=4), intent(in) :: panel_ind

xyz1(1) = ex(1,panel_ind)*xyz(1)+ey(1,panel_ind)*xyz(2)-n(1,panel_ind)*xyz(3)
xyz1(2) = ex(2,panel_ind)*xyz(1)+ey(2,panel_ind)*xyz(2)-n(2,panel_ind)*xyz(3)
xyz1(3) = ex(3,panel_ind)*xyz(1)+ey(3,panel_ind)*xyz(2)-n(3,panel_ind)*xyz(3)

end function ecs_proto2face

subroutine ecs_xyz2ab(pa,pb,px,py,pz,ifc)
use topology_mod, only : ex, ey, n
!transform x,y,z to a,b of panel number ifc
!potentially dangerous: ifc-a,b are not defined for some xyz
real(kind=8), intent(out)   :: pa, pb
real(kind=8), intent(in)    :: px, py, pz
integer(kind=4), intent(in) :: ifc
!locals
real(kind=8) zx, zy, zz

!transform to prototype-face coordinates
zx = ex(1,ifc)*px+ex(2,ifc)*py+ex(3,ifc)*pz
zy = ey(1,ifc)*px+ey(2,ifc)*py+ey(3,ifc)*pz
zz = -n(1,ifc)*px- n(2,ifc)*py- n(3,ifc)*pz

pa = atan(zx/zz);   pb = atan(zy/zz)

end subroutine ecs_xyz2ab

end module ecs_geometry_mod
