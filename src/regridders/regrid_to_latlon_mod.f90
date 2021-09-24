module latlon_regrid_mod

use abstract_regrid_mod,    only : regrid_t
use grid_field_mod,         only : grid_field_t
use halo_mod,               only : halo_t, halo_vec_t

implicit none

type, public :: latlon_interp_t
    integer(kind=4)               :: Nlat, Nlon
    integer(kind=4)               :: ks, ke
    character(len=:), allocatable :: interp_type
    integer(kind=4)               :: stencil_width
    real(kind=8),     allocatable :: wx(:,:,:), wy(:,:,:)
    integer(kind=4),  allocatable :: tile_ind(:,:)
    integer(kind=4),  allocatable :: x_ind(:,:)
    integer(kind=4),  allocatable :: y_ind(:,:)

    contains
    procedure, public :: interpolate => latlon_interpolate
end type latlon_interp_t

type, public, extends(regrid_t) :: latlon_regrid_t
    integer(kind=4) :: Nlon, Nlat
    !interpolation weights
    character(len=:), allocatable :: scalar_grid_type
    type(latlon_interp_t)         :: interp_scalar
    class(halo_t), allocatable    :: scalar_halo

    contains
    procedure :: do_regrid     => do_latlon_scalar_regrid
    procedure :: do_regrid_vec => do_latlon_vector_regrid
end type latlon_regrid_t

contains

subroutine do_latlon_scalar_regrid(this,fout,f)
    class(latlon_regrid_t),    intent(inout) :: this
    real(kind=8),              intent(inout) :: fout(:,:,:)
    type(grid_field_t),        intent(in)    :: f

    integer i, j

    call this%interp_scalar%interpolate(fout,f)

end subroutine do_latlon_scalar_regrid

subroutine do_latlon_vector_regrid(this,uout,vout,u,v)
    class(latlon_regrid_t),    intent(inout) :: this
    real(kind=8),              intent(inout) :: uout(:,:,:), vout(:,:,:)
    class(grid_field_t),       intent(in)    :: u, v
end subroutine do_latlon_vector_regrid

subroutine latlon_interpolate(this,fout,f)
    class(latlon_interp_t), intent(in)    :: this
    real(kind=8),           intent(inout) :: fout(this%Nlon,this%Nlat,this%ks:this%ke)
    type(grid_field_t),     intent(in)    :: f

    integer(kind=4) :: i, j, k, i1, j1
    real(kind=8)    :: p

    do k=this%ks, this%ke
        do j = 1, this%Nlat
            do i = 1, this%Nlon

                fout(i,j,k) = 0.0_8
                do j1=1, this%stencil_width
                    do i1 = 1, this%stencil_width
                        p = f%tile(this%tile_ind(i,j))%p(this%x_ind(i,j)+i1-1,this%y_ind(i,j)+j1-1,k)
                        fout(i,j,k) = fout(i,j,k)+this%wx(i1,i,j)*this%wy(j1,i,j)*p
                    end do
                end do

            end do
        end do
    end do
end subroutine latlon_interpolate

end module latlon_regrid_mod
