module grad_contra_c_sbp42_mod

use abstract_grad_mod,      only : grad_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use mesh_mod,               only : mesh_t, tile_mesh_t
use domain_mod,             only : domain_t
use exchange_abstract_mod,  only : exchange_t

implicit none

type, public, extends(grad_operator_t) :: grad_contra_c_sbp42_t
    class(exchange_t), allocatable  :: exch_f
    class(exchange_t), allocatable  :: exch_covariant
    type(grid_field_t)              :: dx, dy
contains
    procedure, public :: calc_grad => calc_grad_contra_c_sbp42
end type grad_contra_c_sbp42_t

!Differentiation constants
integer(kind=4), parameter :: n_edge = 4
real(kind=8), parameter :: dedge(5,4) = reshape( &
                           [-2.0_8,       3.0_8,     -1.0_8,       0.0_8,       0.0_8,       &
                            -1.0_8,       1.0_8,      0.0_8,       0.0_8,       0.0_8,       &
                             1._8/24._8, -9._8/8._8,  9._8/8._8,  -1._8/24._8,  0.0_8,       &
                            -1._8/71._8,  6._8/71._8,-83._8/71._8, 81._8/71._8,-3._8/71._8], &
                            [5,4])
integer(kind=4), parameter :: d_last_nonzero(4) = [3,2,4,5]
real(kind=8), parameter :: din(4) = [1.0_8/24.0_8, -9.0_8/8.0_8, 9.0_8/8.0_8, -1.0/24.0_8]
real(kind=8), parameter :: proj_op(3) = [15.0_8/8.0_8, -5.0_8/4.0_8, 3.0_8/8.0_8]
real(kind=8), parameter :: sbp_quad(4)= [7._8/18._8, 9._8/8._8, 1._8, 71._8/72._8]

!Interpolation weights
!interpolation from vector points to p-points
real(kind=8), parameter :: wvec_edge(6,4) = reshape( &
[0.3585442775006853_8,  0.6119832845967544_8,  0.20040059830444001_8, -0.1709281604018784_8,  0.0_8,    0.0_8,    &
 0.16724022138832148_8, 0.36635673412335484_8, 0.2655658675883291_8,   0.20083717689999595_8, 0.0_8,    0.0_8,    &
-0.09295505106784016_8, 0.1751328387308906_8,  0.36859947574175345_8,  0.609222736595201_8,  -0.06_8,   0.0_8,    &
-0.04904109392262206_8,-0.04097407434909711_8, 0.16657143046607598_8,  0.42344373780564915_8, 0.5625_8,-0.0625_8],&
   [6,4])
integer(kind=4), parameter :: wvec_last_nonzero(4) = [4,4,5,6]
!interpolation from p-points to vector points
real(kind=8), parameter :: wp_edge(5,4) = reshape( &
[0.9988019158947662_8,  0.37629049812372334_8, -0.24898674393171474_8, -0.12610567008674245_8,  0.0_8, &
 0.5893172370190968_8,  0.2849441265403871_8,   0.16216003586193575_8, -0.03642139942141965_8,  0.0_8, &
 0.21710064816314334_8, 0.23237013413978796_8,  0.3839577872309932_8,   0.16657143046607598_8,  0.0_8, &
-0.18778023255417622_8, 0.17820763584084146_8,  0.6435451442907053_8,   0.42940773411277094_8, -0.06338028169014084_8],&
   [5,4])
integer(kind=4), parameter :: wp_last_nonzero(4) = [4,4,4,5]
!inner interpolation stencil
real(kind=8), parameter :: w_in(4) = [-1._8/16._8, 9._8/16._8, 9._8/16._8, -1._8/16._8]

contains

subroutine calc_grad_contra_c_sbp42(this, gx, gy, f, domain)
    class(grad_contra_c_sbp42_t), intent(inout) :: this
    type(grid_field_t),           intent(inout) :: gx
    type(grid_field_t),           intent(inout) :: gy
    type(grid_field_t),           intent(inout) :: f
    type(domain_t),               intent(in)    :: domain

    integer(kind=4) :: t

    call this%exch_f%do(f,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_covariant_grad_on_tile(this%dx%tile(t), this%dy%tile(t), f%tile(t),            &
                                         domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                                         domain%mesh_o%tile(t),domain%mesh_o%scale)
    end do

    call this%exch_covariant%do_vec(this%dx,this%dy,domain%parcomm)

    do t=domain%partition%ts, domain%partition%te
        call transform_to_contravariant(gx%tile(t),gy%tile(t),           &
                                        this%dx%tile(t),this%dy%tile(t), &
                                        domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                                        domain%mesh_o%tile(t))
    end do

end subroutine calc_grad_contra_c_sbp42

subroutine calc_covariant_grad_on_tile(gx, gy, f, mesh_x, mesh_y, mesh_o, scale)

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(in)    :: f
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_o
    real(kind=8),           intent(in)    :: scale

    integer(kind=4) :: i, i1, i2, j, j1, j2, k, ks, ke, is, ie, js, je, n
    real(kind=8)    :: dx1, dh, dhy(mesh_y%is:mesh_y%ie)

    dx1 = 1.0_8/(scale*mesh_o%hx)

    ks = mesh_o%ks; ke = mesh_o%ke

    do k = ks,ke
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je
        do j=js,je

            do i = is,min(ie,n_edge)
                i1 = d_last_nonzero(i)
                gx%p(i,j,k) =sum(dedge(1:i1,i)*f%p(1:i1,j,k)) *dx1
            end do
            if(is == 1) then
                dh = (proj_op(1)*(f%p(1,j,k)-f%p( 0,j,k))+&
                       proj_op(2)*(f%p(2,j,k)-f%p(-1,j,k))+&
                       proj_op(3)*(f%p(3,j,k)-f%p(-2,j,k)))*dx1
                gx%p(1,j,k) = gx%p(1,j,k)+0.5_8*dh/sbp_quad(1)
            end if

            do i=max(is,n_edge+1),min(ie,mesh_o%nx+1-n_edge)
                gx%p(i,j,k) = sum(f%p(i-2:i+1,j,k)*din(1:4)) *dx1
            end do

            n = mesh_o%nx+1
            do i = max(is,n-n_edge+1),ie
                i2 = n-i+1
                i1 = d_last_nonzero(i2)
                gx%p(i,j,k) =-sum(dedge(1:i1,i2)*f%p(n-1:n-i1:-1,j,k)) *dx1
            end do
            if(ie == n) then
                dh = (proj_op(1)*(f%p(n  ,j,k)-f%p(n-1,j,k))+&
                       proj_op(2)*(f%p(n+1,j,k)-f%p(n-2,j,k))+&
                       proj_op(3)*(f%p(n+2,j,k)-f%p(n-3,j,k)))*dx1
                gx%p(n,j,k) = gx%p(n,j,k)+0.5_8*dh/sbp_quad(1)
            end if
        end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        do j = js,min(je,n_edge)
            j1 = d_last_nonzero(j)
            do i=is,ie
                gy%p(i,j,k) =sum(dedge(1:j1,j)*f%p(i,1:j1,k)) *dx1
            end do
        end do
        if(js==1) then
            do i=is,ie
                dh = (proj_op(1)*(f%p(i,1,k)-f%p(i, 0,k))+&
                      proj_op(2)*(f%p(i,2,k)-f%p(i,-1,k))+&
                      proj_op(3)*(f%p(i,3,k)-f%p(i,-2,k)))*dx1
                gy%p(i,1,k) = gy%p(i,1,k)+0.5_8*dh/sbp_quad(1)
            end do
        end if

        do j=max(js,n_edge+1),min(je,mesh_o%ny+1-n_edge)
            do i=is,ie
                gy%p(i,j,k) = sum(f%p(i,j-2:j+1,k)*din(1:4)) *dx1
            end do
        end do

        n = mesh_o%ny+1
        do j = max(js,n-n_edge+1),je
            j2 = n-j+1
            j1 = d_last_nonzero(j2)
            do i=is,ie
                gy%p(i,j,k) =-sum(dedge(1:j1,j2)*f%p(i,n-1:n-j1:-1,k)) *dx1
            end do
        end do
        if(je==n) then
            do i=is,ie
                dh = (proj_op(1)*(f%p(i,n  ,k)-f%p(i,n-1,k))+&
                      proj_op(2)*(f%p(i,n+1,k)-f%p(i,n-2,k))+&
                      proj_op(3)*(f%p(i,n+2,k)-f%p(i,n-3,k)))*dx1
                gy%p(i,n,k) = gy%p(i,n,k)+0.5_8*dh/sbp_quad(1)
            end do
        end if
    end do
end subroutine calc_covariant_grad_on_tile

subroutine transform_to_contravariant(gx, gy, dx, dy, mesh_x, mesh_y, mesh_o)

    use sbp_mod, only : sbp_apply

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(inout)    :: dx, dy
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_o

    integer(kind=4) :: i, j, k, l, nsrc, ntarg
    integer(kind=4) :: is, ie, js, je, ks, ke
    !indices for grad components:
    integer(kind=4) :: isx, iex, jsx, jex
    integer(kind=4) :: isy, iey, jsy, jey
    !indices for grad components interpolated to p-points:
    integer(kind=4) :: ispx, iepx, jspx, jepx
    integer(kind=4) :: ispy, iepy, jspy, jepy

    integer(kind=4), parameter :: hw = 3 !maximum of interpolation stencil displacement
    real(kind=8)    :: dx_at_p(mesh_y%is:mesh_y%ie,mesh_y%js-hw:mesh_y%je+hw)
    real(kind=8)    :: dy_at_p(mesh_x%is-hw:mesh_x%ie+hw,mesh_x%js:mesh_x%je)

    ks = mesh_o%ks; ke = mesh_o%ke

    !bounding box for required values of gx
    isx = mesh_x%is; iex = mesh_x%ie
    jsx = mesh_x%js; jex = mesh_x%je
    !bounding box for needed values of dy_at_p:
    ispy = interp_p_stencil_start(isx,mesh_o%nx)
    iepy = interp_p_stencil_end  (iex,mesh_o%nx)
    jspy = jsx
    jepy = jex

    !bounding box for required values of gy
    isy = mesh_y%is; iey = mesh_y%ie
    jsy = mesh_y%js; jey = mesh_y%je
    !bounding box for needed values of dx_at_p:
    ispx = isy
    iepx = iey
    jspx = interp_p_stencil_start(jsx,mesh_o%ny)
    jepx = interp_p_stencil_end  (jex,mesh_o%ny)

    do k = ks,ke

       !interpolate dy to p-points
       ntarg = mesh_o%ny
       nsrc  = mesh_o%ny+1
        call sbp_apply(dy_at_p,isx-hw,iex+hw,jsx,jex,ispy,iepy,jspy,jepy,           &
                      dy%p(dy%is:dy%ie,dy%js:dy%je,k), dy%is, dy%ie, dy%js, dy%je, &
                      mesh_o%ny+1,mesh_o%ny,   &
                      wvec_edge,wvec_last_nonzero,w_in,-1,1.0_8,'y')
        !mutiplicate by metric terms in p-points
        do j = jspy,jepy
            do i = ispy, iepy
                dy_at_p(i,j) = mesh_o%G(i,j)*mesh_o%Qi(2,i,j)*dy_at_p(i,j)
            end do
        end do
        !interpolate to u-points
        call sbp_apply(gx%p(gx%is:gx%ie,gx%js:gx%je,k),gx%is,gx%ie,gx%js,gx%je, &
                      isx,iex,jsx,jex,            &
                      dy_at_p, isx-hw,iex+hw,jsx,jex, &
                      mesh_o%nx,mesh_o%nx+1,   &
                      wp_edge,wp_last_nonzero,w_in,-2,1.0_8,'x')
        !interpolate dx to p-points
        call sbp_apply(dx_at_p,isy,iey,jsx-hw,jex+hw,ispx,iepx,jspx,jepx,           &
                      dx%p(dx%is:dx%ie,dx%js:dx%je,k), dx%is, dx%ie, dx%js, dx%je, &
                      mesh_o%nx+1,mesh_o%nx,   &
                      wvec_edge,wvec_last_nonzero,w_in,-1,1.0_8,'x')
        !multiplicate by metric terms in p-points
        do j = jspx,jepx
            do i = ispx, iepx
                dx_at_p(i,j) = mesh_o%G(i,j)*mesh_o%Qi(2,i,j)*dx_at_p(i,j)
            end do
        end do
        !interpolate dx to v-points
        call sbp_apply(gy%p(gy%is:gy%ie,gy%js:gy%je,k),gy%is,gy%ie,gy%js,gy%je, &
                      isy,iey,jsy,jey,            &
                      dx_at_p, isy,iey,jsy-hw,jey+hw, &
                      mesh_o%ny,mesh_o%ny+1,   &
                      wp_edge,wp_last_nonzero,w_in,-2,1.0_8,'y')

        do j=jsx,jex
            do i=isx,iex
                gx%p(i,j,k) = mesh_x%Qi(1,i,j)*dx%p(i,j,k)+gx%p(i,j,k)/mesh_x%G(i,j)
            end do
        end do

        do j=jsy,jey
            do i=isy,iey
                gy%p(i,j,k) = mesh_y%Qi(3,i,j)*dy%p(i,j,k)+gy%p(i,j,k)/mesh_y%G(i,j)
            end do
        end do
    end do
end subroutine transform_to_contravariant

!from u/v points to p points
integer(kind=4) function interp_v_stencil_start(i, n) result(is)
    integer(kind=4) :: i, n
    if(i <= 4) then
        is = 1
    else
        is = min(i-1,n-3)
    end if
end

integer(kind=4) function interp_v_stencil_end(i, n) result(ie)
    integer(kind=4) :: i, n
    if(i <= n-5) then
        ie = max(4,i+2)
    else
        ie = n
    end if
end

!from p-points to u/v points
integer(kind=4) function interp_p_stencil_start(i, n) result(is)
    integer(kind=4) :: i, n
    if(i <= 4) then
        is = 1
    else
        is = min(i-2,n-3)
    end if
end

integer(kind=4) function interp_p_stencil_end(i, n) result(ie)
    integer(kind=4) :: i, n
    if(i <= n-3) then
        ie = max(4,i+1)
    else
        ie = n
    end if
end

end module grad_contra_c_sbp42_mod
