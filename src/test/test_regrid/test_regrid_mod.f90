module test_regrid_mod

use abstract_regrid_mod,    only : regrid_t
use regrid_factory_mod,     only : create_latlon_regrid
use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field
use domain_mod,             only : domain_t
use domain_factory_mod,     only : create_domain
use test_fields_mod,        only : set_scalar_test_field, xyz_f=>xyz_scalar_field_generator

use const_mod,              only : pi

implicit none

contains

subroutine test_regrid(staggering)
    character(len=*), intent(in)    :: staggering

    class(regrid_t),    allocatable :: regrid
    type(grid_field_t)              :: f
    type(domain_t)                  :: domain
    integer(kind=4),    parameter   :: N=32, nz=3
    integer(kind=4),    parameter   :: Nlon = 4*N, Nlat=2*N+1

    real(kind=8),       parameter   :: hlon = 2._8*pi / real(Nlon,8)
    real(kind=8),       parameter   :: hlat = pi / (Nlat-1)

    real(kind=8), dimension(Nlon,Nlat)    :: x_latlon, y_latlon, z_latlon
    real(kind=8), dimension(Nlon,Nlat,nz) :: ptest, ptest_cub, ptrue
    real(kind=8) :: lon, lat
    integer(kind=4) :: i,j

    call create_domain(domain, "cube", staggering, N, nz)
    call create_grid_field(f,0,0,domain%mesh_p)
    call set_scalar_test_field(f,xyz_f,domain%mesh_p,0)

    do j=1, Nlat
        lat = -0.5_8*pi+(j-1)*hlat
        do i=1, Nlon
            lon = (i-1)*hlon
            x_latlon(i,j) = cos(lat)*cos(lon)
            y_latlon(i,j) = cos(lat)*sin(lon)
            z_latlon(i,j) = sin(lat)
        end do
    end do
    call xyz_f%get_scalar_field(ptrue,Nlat*Nlon,nz,x_latlon,y_latlon,z_latlon)

    call create_latlon_regrid(regrid,domain,Nlat=Nlat,Nlon=Nlon,interp_type="linear",&
                              scalar_grid_type=staggering)
    call regrid%do_regrid(ptest,f)

    call create_latlon_regrid(regrid,domain,Nlat=Nlat,Nlon=Nlon,interp_type="cubic",&
                              scalar_grid_type=staggering)
    call regrid%do_regrid(ptest_cub,f)

    print *, "Interpolation error (linear)", maxval(abs(ptest-ptrue))
    print *, "Interpolation error (cubic) ", maxval(abs(ptest_cub-ptrue))

end subroutine test_regrid

end module test_regrid_mod
