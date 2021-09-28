module massflux_Cgrid_mod

use abstract_massflux_mod, only : massflux_operator_t
use grid_field_mod,        only : grid_field_t
use domain_mod,            only : domain_t
use halo_mod,              only : halo_t
use parcomm_mod,           only : parcomm_global

implicit none

type, extends(massflux_operator_t), public :: massflux_chalo_t

    integer(kind=4)            :: order
    class(halo_t), allocatable :: halo

    contains

    procedure :: calc_massflux => calc_c2_massflux

end type massflux_chalo_t

contains

subroutine calc_c2_massflux(this, fx, fy, f, u, v, domain)
    class(massflux_chalo_t), intent(inout) :: this
    type(domain_t),          intent(in)    :: domain
    type(grid_field_t),      intent(inout) :: f, u, v
    !output:
    type(grid_field_t),      intent(inout) :: fx, fy

    integer(kind=4) :: t


    select case(this%order)
    case(2)
        call this%halo%get_halo_scalar(f,domain,1)
        do t = domain%mesh_p%ts, domain%mesh_p%te
            call calc_c2_massflux_tile(fx%tile(t), fy%tile(t), &
                                       f%tile(t), u%tile(t), v%tile(t), &
                                       domain%mesh_x%tile(t), domain%mesh_y%tile(t))
        end do
    case(4)
        call this%halo%get_halo_scalar(f,domain,1)
        do t = domain%mesh_p%ts, domain%mesh_p%te
            call calc_c4_massflux_tile(fx%tile(t), fy%tile(t), &
                                       f%tile(t), u%tile(t), v%tile(t), &
                                       domain%mesh_x%tile(t), domain%mesh_y%tile(t))
        end do
    case default
        call parcomm_global%abort("massflux_chalo_t is currently implemented only for orders=2,4")
    end select

end subroutine calc_c2_massflux

subroutine calc_c2_massflux_tile(fx,fy,f,u,v,mesh_x,mesh_y)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_x, mesh_y

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke

    ks = mesh_x%ks; ke = mesh_x%ke

    do k=ks,ke
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je

        do j=js,je; do i=is,ie
            fx%p(i,j,k) = u%p(i,j,k)*0.5_8*(f%p(i,j,k)+f%p(i-1,j,k))
        end do; end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        do j=js,je; do i=is,ie
            fy%p(i,j,k) = v%p(i,j,k)*0.5_8*(f%p(i,j,k)+f%p(i,j-1,k))
        end do; end do

    end do
end subroutine calc_c2_massflux_tile

subroutine calc_c4_massflux_tile(fx,fy,f,u,v,mesh_x,mesh_y)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_x, mesh_y

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke

    ks = mesh_x%ks; ke = mesh_x%ke

    do k=ks,ke
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je

        do j=js,je; do i=is,ie
            fx%p(i,j,k) = u%p(i,j,k)*(-f%p(i-2,j,k)+7._8*f%p(i-1,j,k)+7._8*f%p(i,j,k)-f%p(i+1,j,k))/12._8
        end do; end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        do j=js,je; do i=is,ie
            fy%p(i,j,k) = v%p(i,j,k)*(-f%p(i,j-2,k)+7._8*f%p(i,j-1,k)+7._8*f%p(i,j,k)-f%p(i,j+1,k))/12._8
        end do; end do

    end do
end subroutine calc_c4_massflux_tile

end module massflux_Cgrid_mod
