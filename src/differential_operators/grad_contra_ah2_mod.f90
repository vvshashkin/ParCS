module grad_contra_ah2_mod

use domain_mod,             only : domain_t
use abstract_grad_mod,      only : grad_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use exchange_abstract_mod,  only : exchange_t
use parcomm_mod,            only : parcomm_global

implicit none

private

type, public, extends(grad_operator_t) :: grad_contra_ah2_t
    class(exchange_t), allocatable     :: exch_halo
    character(:),      allocatable     :: subtype
    !quantities to calculate the components in the edge points (gx,gy)
    !using components from neighbouring faces (gx, gy):
    !gx = gx', gy = g=q*gx'+gy'
    type(tile_edge_transform_t), allocatable     :: q(:)
contains
    procedure, public :: calc_grad => calc_grad_ah2
end type grad_contra_ah2_t

type, public :: tile_edge_transform_t
    real(kind=8), dimension(:), allocatable :: qb, qr, qt, ql
end type tile_edge_transform_t

contains

subroutine calc_grad_ah2(this, gx, gy, f, domain)
    class(grad_contra_ah2_t), intent(inout) :: this
    type(domain_t),           intent(in)    :: domain
    type(grid_field_t),       intent(inout) :: f
    !output:
    type(grid_field_t),     intent(inout) :: gx, gy

    integer(kind=4) :: t

    call this%exch_halo%do(f,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        gx%tile(t)%p = 0.0_8
        gy%tile(t)%p = 0.0_8
        call calc_grad_on_tile(gx%tile(t), gy%tile(t), f%tile(t),  &
                               domain%mesh_xy%tile(t), this%q(t),  &
                               domain%mesh_xy%scale)
    end do

end subroutine calc_grad_ah2

subroutine calc_grad_on_tile(gx, gy, f, mesh, q, scale)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),           intent(inout) :: gx, gy
    type(tile_field_t),           intent(in)    :: f
    type(tile_mesh_t),            intent(in)    :: mesh
    type(tile_edge_transform_t),  intent(in)    :: q
    real(kind=8),                 intent(in)    :: scale

    real(kind=8)    :: hx, mult_loc
    real(kind=8)    :: fdx, fdy, fdx1, fdy1, fdx0, fdy0
    integer(kind=4) :: ks, ke, js, je, is, ie, n, i, j, k

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx

    do k = ks, ke
        !Corner cases:
        if(js==1 .and. is == 1) then
            i = is; j = js
            fdx0 = (f%p(i+1,j,k)-f%p(i,j,k))/(hx*scale)
            fdy0 = (f%p(i,j+1,k)-f%p(i,j,k))/(hx*scale)
            fdx1 = (f%p(i-1,j,k)-f%p(i-2,j,k))/(hx*scale)
            fdy1 = (f%p(i,j-1,k)-f%p(i,j-2,k))/(hx*scale)
            fdx = (2.0_8*fdx0+fdx1+q%ql(j)*fdy0) / 3.0_8
            fdy = (2.0_8*fdy0+fdy1+q%qb(i)*fdx0) / 3.0_8
            gx%p(i,j,k) = (mesh%Qi(1,i,j)*fdx + mesh%Qi(2,i,j)*fdy)
            gy%p(i,j,k) = (mesh%Qi(3,i,j)*fdy + mesh%Qi(2,i,j)*fdx)
        end if
        if(js==1 .and. ie == mesh%nx+1) then
            i = ie; j = js
            fdx0 = (f%p(i,j,k)-f%p(i-1,j,k))/(hx*scale)
            fdy0 = (f%p(i,j+1,k)-f%p(i,j,k))/(hx*scale)
            fdx1 = (f%p(i+2,j,k)-f%p(i+1,j,k))/(hx*scale)
            fdy1 = (f%p(i,j-1,k)-f%p(i,j-2,k))/(hx*scale)
            fdx = (2.0_8*fdx0+fdx1+q%qr(j)*fdy0) / 3.0_8
            fdy = (2.0_8*fdy0+fdy1+q%qb(i)*fdx0) / 3.0_8
            gx%p(i,j,k) = (mesh%Qi(1,i,j)*fdx + mesh%Qi(2,i,j)*fdy)
            gy%p(i,j,k) = (mesh%Qi(3,i,j)*fdy + mesh%Qi(2,i,j)*fdx)
        end if
        if(je==mesh%ny+1 .and. ie == mesh%nx+1) then
            i = ie; j = je
            fdx0 = (f%p(i,j,k)-f%p(i-1,j,k))/(hx*scale)
            fdy0 = (f%p(i,j,k)-f%p(i,j-1,k))/(hx*scale)
            fdx1 = (f%p(i+2,j,k)-f%p(i+1,j,k))/(hx*scale)
            fdy1 = (f%p(i,j+2,k)-f%p(i,j+1,k))/(hx*scale)
            fdx = (2.0_8*fdx0+fdx1+q%qr(j)*fdy0) / 3.0_8
            fdy = (2.0_8*fdy0+fdy1+q%qt(i)*fdx0) / 3.0_8
            gx%p(i,j,k) = (mesh%Qi(1,i,j)*fdx + mesh%Qi(2,i,j)*fdy)
            gy%p(i,j,k) = (mesh%Qi(3,i,j)*fdy + mesh%Qi(2,i,j)*fdx)
        end if
        if(je==mesh%ny+1 .and. is == 1) then
            i = is; j = je
            fdx0 = (f%p(i+1,j,k)-f%p(i,j,k))/(hx*scale)
            fdy0 = (f%p(i,j,k)-f%p(i,j-1,k))/(hx*scale)
            fdx1 = (f%p(i-1,j,k)-f%p(i-2,j,k))/(hx*scale)
            fdy1 = (f%p(i,j+2,k)-f%p(i,j+1,k))/(hx*scale)
            fdx = (2.0_8*fdx0+fdx1+q%ql(j)*fdy0) / 3.0_8
            fdy = (2.0_8*fdy0+fdy1+q%qt(i)*fdx0) / 3.0_8
            gx%p(i,j,k) = (mesh%Qi(1,i,j)*fdx + mesh%Qi(2,i,j)*fdy)
            gy%p(i,j,k) = (mesh%Qi(3,i,j)*fdy + mesh%Qi(2,i,j)*fdx)
        end if
        !Edge cases:
        if(js == 1) then
            j=1
            do i = max(is,2), min(ie,mesh%nx)
                fdx  = (f%p(i+1,j,k) - f%p(i-1,j,k)) / (2.0_8*(hx*scale))
                fdy  = (f%p(i,j+1,k) - f%p(i,j,k)) / (hx*scale)
                fdy1 = (f%p(i,0,k) - f%p(i,-1,k)) / (hx*scale)
                fdy = 0.5_8*fdy+0.5_8*(q%qb(i)*fdx+fdy1)
                gx%p(i,j,k) = (mesh%Qi(1,i,j)*fdx + mesh%Qi(2,i,j)*fdy)
                gy%p(i,j,k) = (mesh%Qi(3,i,j)*fdy + mesh%Qi(2,i,j)*fdx)
            end do
        end if
        if(je == mesh%ny+1) then
            j=je
            do i = max(is,2), min(ie,mesh%nx)
                fdx  = (f%p(i+1,j,k) - f%p(i-1,j,k)) / (2.0_8*(hx*scale))
                fdy  = (f%p(i,j,k) - f%p(i,j-1,k)) / (hx*scale)
                fdy1 = (f%p(i,j+2,k) - f%p(i,j+1,k)) / (hx*scale)
                fdy = 0.5_8*fdy+0.5_8*(q%qt(i)*fdx+fdy1)
                gx%p(i,j,k) = (mesh%Qi(1,i,j)*fdx + mesh%Qi(2,i,j)*fdy)
                gy%p(i,j,k) = (mesh%Qi(3,i,j)*fdy + mesh%Qi(2,i,j)*fdx)
            end do
        end if
        if(is == 1) then
            i=1
            do j = max(js,2), min(je,mesh%ny)
                fdx  = (f%p(i+1,j,k) - f%p(i,j,k)) / (hx*scale)
                fdx1 = (f%p(i-1,j,k) - f%p(i-2,j,k)) / (hx*scale)
                fdy  = (f%p(i,j+1,k) - f%p(i,j-1,k)) / (2.0_8*(hx*scale))
                fdx = 0.5_8*fdx+0.5_8*(q%ql(j)*fdy+fdx1)
                gx%p(i,j,k) = (mesh%Qi(1,i,j)*fdx + mesh%Qi(2,i,j)*fdy)
                gy%p(i,j,k) = (mesh%Qi(3,i,j)*fdy + mesh%Qi(2,i,j)*fdx)
            end do
        end if
        if(ie == mesh%nx+1) then
            i=ie
            do j = max(js,2), min(je,mesh%ny)
                fdx  = (f%p(i,j,k) - f%p(i-1,j,k)) / (hx*scale)
                fdx1 = (f%p(i+2,j,k) - f%p(i+1,j,k)) / (hx*scale)
                fdy  = (f%p(i,j+1,k) - f%p(i,j-1,k)) / (2.0_8*(hx*scale))
                fdx = 0.5_8*fdx+0.5_8*(q%qr(j)*fdy+fdx1)
                gx%p(i,j,k) = (mesh%Qi(1,i,j)*fdx + mesh%Qi(2,i,j)*fdy)
                gy%p(i,j,k) = (mesh%Qi(3,i,j)*fdy + mesh%Qi(2,i,j)*fdx)
            end do
        end if
        !Regular points:
        do j = max(js,2), min(je,mesh%ny)
            do i = max(is,2), min(ie,mesh%nx)
                fdx = (f%p(i+1,j,k) - f%p(i-1,j,k)) / (2.0_8*(hx*scale))
                fdy = (f%p(i,j+1,k) - f%p(i,j-1,k)) / (2.0_8*(hx*scale))
                gx%p(i,j,k) = (mesh%Qi(1,i,j)*fdx + mesh%Qi(2,i,j)*fdy)
                gy%p(i,j,k) = (mesh%Qi(3,i,j)*fdy + mesh%Qi(2,i,j)*fdx)
            end do
        end do
    end do

end subroutine calc_grad_on_tile

end module grad_contra_ah2_mod
