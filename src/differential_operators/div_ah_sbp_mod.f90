module div_ah_sbp_mod

use domain_mod,             only : domain_t
use abstract_div_mod,       only : div_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use halo_mod,               only : halo_vec_t
use exchange_abstract_mod,  only : exchange_t
use parcomm_mod,            only : parcomm_global
use sbp_mod,                only : sbp_diff

implicit none

type, public, extends(div_operator_t) :: div_ah_sbp_t
    class(exchange_t), allocatable     :: exch_uv_interior
    class(exchange_t), allocatable     :: exch_div_edges
    character(:),      allocatable     :: sbp_operator_name
contains
    procedure, public :: calc_div => calc_div_ah_sbp
end type div_ah_sbp_t

contains

subroutine calc_div_ah_sbp(this, div, u, v, domain)
    class(div_ah_sbp_t),    intent(inout) :: this
    type(domain_t),         intent(in)    :: domain
    type(grid_field_t),     intent(inout) :: u, v
    !output
    type(grid_field_t),     intent(inout) :: div

    integer(kind=4), parameter :: halo_width = 1
    integer(kind=4) :: t

    call this%exch_uv_interior%do_vec(u,v,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_div_on_tile(div%tile(t), u%tile(t), v%tile(t),  &
                              domain%mesh_xy%tile(t), this%sbp_operator_name,&
                              domain%mesh_xy%scale)
    end do

    call this%exch_div_edges%do(div,domain%parcomm)
    do t = domain%partition%ts, domain%partition%te
        call sync_div_on_edges(div%tile(t), domain%mesh_xy%tile(t))
    end do

end subroutine calc_div_ah_sbp

subroutine calc_div_on_tile(div, u, v, mesh, sbp_oper_name,scale)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),     intent(inout) :: div
    type(tile_field_t),     intent(in)    :: u, v
    type(tile_mesh_t),      intent(in)    :: mesh
    character(len=*),       intent(in)    :: sbp_oper_name
    real(kind=8),           intent(in)    :: scale

    real(kind=8)    :: hx
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k
    real(kind=8)    :: Gu(u%is:u%ie,u%js:u%je)
    real(kind=8)    :: Gv(v%is:v%ie,v%js:v%je)
    real(kind=8)    :: Dx(mesh%is:mesh%ie,mesh%js:mesh%je)
    real(kind=8)    :: Dy(mesh%is:mesh%ie,mesh%js:mesh%je)

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx

    do k = ks, ke

        do j = max(js-3,1), min(je+3,mesh%ny+1)
            do i = max(is-3,1),min(ie+3,mesh%nx+1)
                Gu(i,j) = mesh%G(i,j)*u%p(i,j,k)
            end do
        end do
        call sbp_diff(Dx,is,ie,js,je,is,ie,js,je,"x",sbp_oper_name, &
                      Gu,u%is,u%ie,u%js,u%je,           &
                      mesh%nx+1,mesh%nx+1)

        do j = max(js-3,1), min(je+3,mesh%ny+1)
            do i = max(is-3,1),min(ie+3,mesh%nx+1)
                Gv(i,j) = mesh%G(i,j)*v%p(i,j,k)
            end do
        end do
        call sbp_diff(Dy,is,ie,js,je,is,ie,js,je,"y",sbp_oper_name, &
                      Gv,v%is,v%ie,v%js,v%je,           &
                      mesh%ny+1,mesh%ny+1)

        do j = js, je
            do i = is,ie
                div%p(i,j,k) = (Dx(i,j)+Dy(i,j)) / (mesh%G(i,j)*hx*scale)
            end do
        end do
    end do

end subroutine calc_div_on_tile

subroutine sync_div_on_edges(div,mesh)
    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),     intent(inout) :: div
    type(tile_mesh_t),      intent(in)    :: mesh

    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k, ip1,im1,jp1,jm1,mx,my

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    if(is == 1) then
        do k=ks,ke
            if(js == 1) then
                div%p(1,1,k) = (div%p(1,1,k)+div%p(0,1,k)+div%p(1,0,k))/3.0_8
            end if
            do j = max(2,js),min(mesh%ny,je)
                div%p(1,j,k) = 0.5_8*(div%p(0,j,k)+div%p(1,j,k))
            end do
            if(je == mesh%ny+1) then
                div%p(1,je,k) = (div%p(1,je,k)+div%p(0,je,k)+div%p(1,je+1,k))/3.0_8
            end if
        end do
    end if
    if(ie == mesh%nx+1) then
        do k=ks,ke
            if(js == 1) then
                div%p(ie,1,k) = (div%p(ie,1,k)+div%p(ie+1,1,k)+div%p(ie,0,k))/3.0_8
            end if
            do j = max(2,js),min(mesh%ny,je)
                div%p(ie,j,k) = 0.5_8*(div%p(ie,j,k)+div%p(ie+1,j,k))
            end do
            if(je == mesh%ny+1) then
                div%p(ie,je,k) = (div%p(ie,je,k)+div%p(ie+1,je,k)+div%p(ie,je+1,k))/3.0_8
            end if
        end do
    end if
    if(js == 1) then
        do k=ks,ke
            do i = max(2,is),min(mesh%nx,ie)
                div%p(i,1,k) = 0.5_8*(div%p(i,1,k)+div%p(i,0,k))
            end do
        end do
    end if
    if(je == mesh%ny+1) then
        do k=ks,ke
            do i = max(2,is),min(mesh%nx,ie)
                div%p(i,je,k) = 0.5_8*(div%p(i,je,k)+div%p(i,je+1,k))
            end do
        end do
    end if
end subroutine sync_div_on_edges

end module div_ah_sbp_mod
