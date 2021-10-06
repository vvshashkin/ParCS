module div_c_sbp42_mod

use abstract_div_mod,   only : div_operator_t
use domain_mod,         only : domain_t
use mesh_mod,           only : mesh_t, tile_mesh_t
use grid_field_mod,     only : grid_field_t, tile_field_t
use exchange_halo_mod,  only : exchange_t
use halo_mod,           only : halo_vec_t

implicit none

type, public, extends(div_operator_t) :: div_c_sbp42_t
    class(exchange_t), allocatable :: exch_halo
    class(halo_vec_t), allocatable :: halo_procedure
contains
    procedure, public :: calc_div => calc_div_c_sbp42
end type div_c_sbp42_t

real(kind=8), parameter :: d1(4) = [-79.0_8/78.0_8, 27.0_8/26.0_8, -1.0/26.0_8, 1.0/78.0_8]
real(kind=8), parameter :: d2(4) = [2.0_8/21.0_8, -9.0_8/7.0_8, 9.0_8/7.0_8, -2.0_8/21.0_8]
real(kind=8), parameter :: d3(5) = [1.0_8/75.0_8, 0.0_8, -27.0_8/25.0_8, 83.0_8/75.0_8, -1.0_8/25.0_8]
real(kind=8), parameter :: din(4) = [1.0_8/24.0_8, -9.0_8/8.0_8, 9.0_8/8.0_8, -1.0/24.0_8]
real(kind=8), parameter :: proj_op(3) = [15.0_8/8.0_8, -5.0_8/4.0_8, 3.0_8/8.0_8]
real(kind=8), parameter :: sbp_quad(3)= [13.0_8/12.0_8, 7.0_8/8.0_8, 25.0_8/24.0_8]
real(kind=8), parameter :: penalty(3) = proj_op(1:3)/sbp_quad(1:3)

contains

subroutine calc_div_c_sbp42(this, div, u, v, domain)

    class(div_c_sbp42_t),   intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: u, v
    !out put
    type(grid_field_t),     intent(inout) :: div

    integer(kind=4) :: i, j, k, t
    real(kind=8) hx

    call this%exch_halo%do_vec(u,v,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_div_on_tile_sbp42(div%tile(t), u%tile(t), v%tile(t),            &
                                    domain%mesh_x%tile(t), domain%mesh_y%tile(t), &
                                    domain%mesh_o%tile(t),domain%mesh_o%scale)
    end do

end subroutine calc_div_c_sbp42

subroutine calc_div_on_tile_sbp42(div, u, v, mesh_u, mesh_v, mesh_p, scale)

    type(tile_field_t),  intent(inout) :: div
    type(tile_field_t),  intent(in)    :: u, v
    type(tile_mesh_t),   intent(in)    :: mesh_u, mesh_v, mesh_p
    real(kind=8),        intent(in)    :: scale

    real(kind=8)    :: hx
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k
    real(kind=8)    :: dx(mesh_p%is:mesh_p%ie,mesh_p%js:mesh_p%je)
    real(kind=8)    :: dy(mesh_p%is:mesh_p%ie,mesh_p%js:mesh_p%je)
    real(kind=8)    :: du, dv(mesh_p%is:mesh_p%ie)

    is = mesh_p%is; ie = mesh_p%ie
    js = mesh_p%js; je = mesh_p%je
    ks = mesh_p%ks; ke = mesh_p%ke

    hx = mesh_p%hx

    do k = ks, ke

        do j=js,je
            if(is<=3) du = mesh_u%G(1,j)*(u%p(1,j,k)-u%p(0,j,k))

            if(is == 1) then
                dx(1,j) = d1(1)*mesh_u%G(1,j)*u%p(1,j,k)+&
                          d1(2)*mesh_u%G(2,j)*u%p(2,j,k)+&
                          d1(3)*mesh_u%G(3,j)*u%p(3,j,k)+&
                          d1(4)*mesh_u%G(4,j)*u%p(4,j,k)+&
                          0.5_8*penalty(1)*du
            end if
            if(is<=2 .and. ie >= 2) then
                dx(2,j) = d2(1)*mesh_u%G(1,j)*u%p(1,j,k)+&
                          d2(2)*mesh_u%G(2,j)*u%p(2,j,k)+&
                          d2(3)*mesh_u%G(3,j)*u%p(3,j,k)+&
                          d2(4)*mesh_u%G(4,j)*u%p(4,j,k)+&
                          0.5_8*penalty(2)*du
            end if
            if(is<=3 .and. ie >= 3) then
                dx(3,j) = d3(1)*mesh_u%G(1,j)*u%p(1,j,k)+&
                          !zero coefficient:!d3(2)*mesh_u%G(2,j)*u%p(2,j,k)+&
                          d3(3)*mesh_u%G(3,j)*u%p(3,j,k)+&
                          d3(4)*mesh_u%G(4,j)*u%p(4,j,k)+&
                          d3(5)*mesh_u%G(5,j)*u%p(5,j,k)+&
                          0.5_8*penalty(3)*du
            end if

            do i=max(4,is),min(ie,mesh_p%globnx-3)
                dx(i,j) = (din(1)*mesh_u%G(i-1,j)*u%p(i-1,j,k)+ &
                           din(2)*mesh_u%G(i  ,j)*u%p(i  ,j,k)+ &
                           din(3)*mesh_u%G(i+1,j)*u%p(i+1,j,k)+ &
                           din(4)*mesh_u%G(i+2,j)*u%p(i+2,j,k))
            end do

            if(ie>=mesh_p%globnx-2) then
                du = mesh_u%G(mesh_p%globnx+1,j)*(u%p(mesh_p%globnx+2,j,k)-u%p(mesh_p%globnx+1,j,k))
            end if
            if(is<=mesh_p%globnx-2 .and. ie >= mesh_p%globnx-2) then
                i = mesh_p%globnx-2
                dx(i,j) =-d3(1)*mesh_u%G(i+3,j)*u%p(i+3,j,k)-&
                          !zero coefficient:!d3(2)*mesh_u%G(i+2,j)*u%p(i+2,j,k)+&
                          d3(3)*mesh_u%G(i+1,j)*u%p(i+1,j,k)-&
                          d3(4)*mesh_u%G(i  ,j)*u%p(i  ,j,k)-&
                          d3(5)*mesh_u%G(i-1,j)*u%p(i-1,j,k)+&
                          0.5_8*penalty(3)*du
            end if
            if(is<=mesh_p%globnx-1 .and. ie >= mesh_p%globnx-1) then
                i = mesh_p%globnx-1
                dx(i,j) =-d2(1)*mesh_u%G(i+2,j)*u%p(i+2,j,k)-&
                          d2(2)*mesh_u%G(i+1,j)*u%p(i+1,j,k)-&
                          d2(3)*mesh_u%G(i  ,j)*u%p(i  ,j,k)-&
                          d2(4)*mesh_u%G(i-1,j)*u%p(i-1,j,k)+&
                          0.5_8*penalty(2)*du
            end if
            if(ie == mesh_p%globnx) then
                i = mesh_p%globnx
                dx(i,j) =-d1(1)*mesh_u%G(i+1,j)*u%p(i+1,j,k)-&
                          d1(2)*mesh_u%G(i  ,j)*u%p(i  ,j,k)-&
                          d1(3)*mesh_u%G(i-1,j)*u%p(i-1,j,k)-&
                          d1(4)*mesh_u%G(i-2,j)*u%p(i-2,j,k)+&
                          0.5_8*penalty(1)*du
            end if
        end do

        if(js<=3) then
            do i=is,ie
                dv(i) = mesh_v%G(i,1)*(v%p(i,1,k)-v%p(i,0,k))
            end do
        end if

        if(js == 1) then
            do i=is,ie
                dy(i,1) = d1(1)*mesh_v%G(i,1)*v%p(i,1,k)+&
                          d1(2)*mesh_v%G(i,2)*v%p(i,2,k)+&
                          d1(3)*mesh_v%G(i,3)*v%p(i,3,k)+&
                          d1(4)*mesh_v%G(i,4)*v%p(i,4,k)+&
                          0.5_8*penalty(1)*dv(i)
            end do
        end if
        if(js<=2 .and. je >= 2) then
            do i=is,ie
                dy(i,2) = d2(1)*mesh_v%G(i,1)*v%p(i,1,k)+&
                          d2(2)*mesh_v%G(i,2)*v%p(i,2,k)+&
                          d2(3)*mesh_v%G(i,3)*v%p(i,3,k)+&
                          d2(4)*mesh_v%G(i,4)*v%p(i,4,k)+&
                          0.5_8*penalty(2)*dv(i)
            end do
        end if
        if(js<=3 .and. je >= 3) then
            do i=is,ie
                dy(i,3) = d3(1)*mesh_v%G(i,1)*v%p(i,1,k)+&
                          !zero coefficient:!d3(2)*mesh_u%G(i,2)*u%p(i,2,k)+&
                          d3(3)*mesh_v%G(i,3)*v%p(i,3,k)+&
                          d3(4)*mesh_v%G(i,4)*v%p(i,4,k)+&
                          d3(5)*mesh_v%G(i,5)*v%p(i,5,k)+&
                          0.5_8*penalty(3)*dv(i)
            end do
        end if
        do j=max(4,js),min(je,mesh_p%globny-3)
            do i=is,ie
                dy(i,j) = (din(1)*mesh_v%G(i,j-1)*v%p(i,j-1,k)+ &
                           din(2)*mesh_v%G(i,j  )*v%p(i,j  ,k)+ &
                           din(3)*mesh_v%G(i,j+1)*v%p(i,j+1,k)+ &
                           din(4)*mesh_v%G(i,j+2)*v%p(i,j+2,k))
            end do
        end do

        if(je>=mesh_p%globny-2) then
            do i=is,ie
                dv(i) = mesh_v%G(i,mesh_p%globny+1)*(v%p(i,mesh_p%globny+2,k)-v%p(i,mesh_p%globny+1,k))
            end do
        end if

        if(js<=mesh_p%globny-2 .and. je >= mesh_p%globny-2) then
            do i=is,ie
                j = mesh_p%globny-2
                dy(i,j) =-d3(1)*mesh_v%G(i,j+3)*v%p(i,j+3,k)-&
                          !zero coefficient:!d3(2)*mesh_v%G(i,j+2)*v%p(i,j+2,k)+&
                          d3(3)*mesh_v%G(i,j+1)*v%p(i,j+1,k)-&
                          d3(4)*mesh_v%G(i,j  )*v%p(i,j  ,k)-&
                          d3(5)*mesh_v%G(i,j-1)*v%p(i,j-1,k)+&
                          0.5_8*penalty(3)*dv(i)
            end do
        end if
        if(js<=mesh_p%globny-1 .and. je >= mesh_p%globny-1) then
            do i=is,ie
                j = mesh_p%globny-1
                dy(i,j) =-d2(1)*mesh_v%G(i,j+2)*v%p(i,j+2,k)-&
                          d2(2)*mesh_v%G(i,j+1)*v%p(i,j+1,k)-&
                          d2(3)*mesh_v%G(i,j  )*v%p(i,j  ,k)-&
                          d2(4)*mesh_v%G(i,j-1)*v%p(i,j-1,k)+&
                          0.5_8*penalty(2)*dv(i)
            end do
        end if
        if(je == mesh_p%globny) then
            do i=is,ie
                j = mesh_p%globny
                dy(i,j) =-d1(1)*mesh_v%G(i,j+1)*v%p(i,j+1,k)-&
                          d1(2)*mesh_v%G(i,j  )*v%p(i,j  ,k)-&
                          d1(3)*mesh_v%G(i,j-1)*v%p(i,j-1,k)-&
                          d1(4)*mesh_v%G(i,j-2)*v%p(i,j-2,k)+&
                          0.5_8*penalty(1)*dv(i)
            end do
        end if

        do j = js,je
            do i = is, ie
                div%p(i,j,k) =  (dx(i,j)+dy(i,j))/(mesh_p%G(i,j)*hx*scale)
            end do
        end do
    end do

end subroutine calc_div_on_tile_sbp42

end module div_c_sbp42_mod
