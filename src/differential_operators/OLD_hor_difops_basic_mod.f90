module hor_difops_basic_mod

use grid_field_mod, only : block_t
use mesh_mod,       only : mesh_t
use const_mod,      only : pi, radz

implicit none

contains

subroutine cl_gradient_contra_c2(gx, gy, f, mesh, multiplier)
    type(block_t),  intent(inout) :: gx, gy
    type(block_t),  intent(in)    :: f
    type(mesh_t),           intent(in)    :: mesh
    real(kind=8), optional, intent(in)    :: multiplier

    integer i, j, k
    integer is, ie, js, je, ks, ke
    real(kind=8), allocatable :: fdx(:,:), fdy(:,:)
    real(kind=8) mult, hx, fdy_at_x, fdx_at_y

    if(present(multiplier)) then
        mult = multiplier
    else
        mult = 1._8
    end if

    hx = mesh%hx

    js  = f%js;    je  = f%je
    is  = f%is;    ie  = f%ie
    ks  = f%ks;    ke  = f%ke

    allocate(fdx(is-1:ie,js-1:je+1), fdy(is-1:ie+1,js-1:je))

    do k = ks, ke

        do j= js-1, je+1
            do i= is-1, ie
                fdx(i,j) = mult*(f%p(i+1,j,k)-f%p(i,j,k))/(hx*radz)
            end do
        end do

        do j= js-1, je
            do i= is-1, ie+1
                fdy(i,j) = mult*(f%p(i,j+1,k)-f%p(i,j,k))/(hx*radz)
            end do
        end do

        !transform to contravariant components
        do j=js,je
            do i=is-1,ie
                fdy_at_x = 0.25_8*(fdy(i,j)+fdy(i+1,j)+fdy(i,j-1)+fdy(i+1,j-1))
                gx%p(i,j,k) = mesh%Qiu(1,i,j)*fdx(i,j) + mesh%Qiu(2,i,j)*fdy_at_x
            end do
        end do
        do j=js-1,je
            do i=is,ie
                fdx_at_y = 0.25_8*(fdx(i,j)+fdx(i-1,j)+fdx(i,j+1)+fdx(i-1,j+1))
                gy%p(i,j,k) = mesh%Qiv(3,i,j)*fdy(i,j) + mesh%Qiv(2,i,j)*fdx_at_y
            end do
        end do

    end do

end subroutine cl_gradient_contra_c2

subroutine cl_divergence_cgr2(div, u, v, mesh, multiplier)
    type(block_t),          intent(inout) :: div
    type(block_t),          intent(in)    :: u, v
    type(mesh_t),           intent(in)    :: mesh
    real(kind=8), optional, intent(in)    :: multiplier

    integer i, j, k
    integer is, ie, js, je, ks, ke
    real(kind=8), allocatable :: fdx(:,:), fdy(:,:)
    real(kind=8) mult, hx

    if(present(multiplier)) then
        mult = multiplier
    else
        mult = 1._8
    end if

    js  = div%js;    je  = div%je
    is  = div%is;    ie  = div%ie
    ks  = div%ks;    ke  = div%ke

    hx = mesh%hx

    do k=ks, ke
        do j=js, je
            do i=is, ie
                div%p(i,j,k) = mult*(mesh%Gu(i,j)*u%p(i,j,k)-mesh%Gu(i-1,j)*u%p(i-1,j,k) +  &
                                     mesh%Gv(i,j)*v%p(i,j,k)-mesh%Gv(i,j-1)*v%p(i,j-1,k))/  &
                                     (mesh%G(i,j)*radz*hx)
            end do
        end do
    end do

end subroutine cl_divergence_cgr2

!fake test operator: 0-grad on any field f
subroutine cl_gradient_0(gx, gy, f, mesh, multiplier)
    type(block_t),          intent(inout) :: gx, gy
    type(block_t),          intent(in)    :: f
    type(mesh_t),           intent(in)    :: mesh
    real(kind=8), optional, intent(in)    :: multiplier

    gx%p = 0._8; gy%p = 0._8

end subroutine cl_gradient_0

!Fake test operator: 0-divergence on any u,v, field
subroutine cl_divergence_0(div, u, v, mesh, multiplier)
    type(block_t),          intent(inout) :: div
    type(block_t),          intent(in)    :: u, v
    type(mesh_t),           intent(in)    :: mesh
    real(kind=8), optional, intent(in)    :: multiplier

    div%p = 0._8

end subroutine cl_divergence_0

end module hor_difops_basic_mod
