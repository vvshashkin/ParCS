module grad_contra_ah_sbp_mod

use domain_mod,             only : domain_t
use abstract_grad_mod,      only : grad_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use exchange_abstract_mod,  only : exchange_t
use parcomm_mod,            only : parcomm_global
use sbp_mod,                only : sbp_diff
use grad_contra_ah2_mod,    only : tile_edge_transform_t

implicit none

private

type, public, extends(grad_operator_t) :: grad_contra_ah_sbp_t
    class(exchange_t), allocatable     :: exch_scalar_interior
    class(exchange_t), allocatable     :: exch_vector_edges
    character(:),      allocatable     :: subtype
    character(:),      allocatable     :: sbp_operator_name
    !quantities to calculate the components in the edge points (gx,gy)
    !using components from neighbouring faces:
    type(tile_edge_transform_t), allocatable     :: q(:)
contains
    procedure, public :: calc_grad => calc_grad_ah_sbp
end type grad_contra_ah_sbp_t

contains

subroutine calc_grad_ah_sbp(this, gx, gy, f, domain, multiplier)
    class(grad_contra_ah_sbp_t),   intent(inout) :: this
    type(domain_t),                intent(in)    :: domain
    type(grid_field_t),            intent(inout) :: f
    real(kind=8), optional,        intent(in)    :: multiplier
    !output:
    type(grid_field_t),     intent(inout) :: gx, gy

    integer(kind=4) :: t
    real(kind=8)    :: mult_loc

    mult_loc = 1.0_8
    if (present(multiplier)) mult_loc = multiplier

    call this%exch_scalar_interior%do(f,domain%parcomm)

    do t = domain%partition%ts, domain%partition%te
        call calc_grad_on_tile(gx%tile(t), gy%tile(t), f%tile(t),  &
                               domain%mesh_xy%tile(t), this%sbp_operator_name, mult_loc)
    end do

    call this%exch_vector_edges%do_vec(gx,gy,domain%parcomm)
    do t = domain%partition%ts, domain%partition%te
        call syncronize_grad_on_edges(gx%tile(t), gy%tile(t), this%q(t), &
                                      domain%mesh_xy%tile(t))
    end do

end subroutine calc_grad_ah_sbp

subroutine calc_grad_on_tile(gx, gy, f, mesh, sbp_oper, multiplier)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(in)    :: f
    type(tile_mesh_t),      intent(in)    :: mesh
    character(len=*),       intent(in)    :: sbp_oper
    real(kind=8),           intent(in)    :: multiplier

    real(kind=8)    :: hx, mult_loc
    real(kind=8)    :: fdx, fdy, fdx1, fdy1, fdx0, fdy0
    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k

    real(kind=8)    :: Dx(mesh%is:mesh%ie,mesh%js:mesh%je)
    real(kind=8)    :: Dy(mesh%is:mesh%ie,mesh%js:mesh%je)

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx

    do k = ks, ke
        call sbp_diff("x",sbp_oper,f%p(f%is:f%ie,f%js:f%je,k),f%is,f%ie,f%js,f%je,mesh%nx+1,is,ie,js,je,Dx)
        call sbp_diff("y",sbp_oper,f%p(f%is:f%ie,f%js:f%je,k),f%is,f%ie,f%js,f%je,mesh%ny+1,is,ie,js,je,Dy)

        do j=js,je
            do i=is,ie
                gx%p(i,j,k) = (mesh%Qi(1,i,j)*Dx(i,j) + mesh%Qi(2,i,j)*Dy(i,j))*multiplier/hx
                gy%p(i,j,k) = (mesh%Qi(3,i,j)*Dy(i,j) + mesh%Qi(2,i,j)*Dx(i,j))*multiplier/hx
            end do
        end do
    end do

end subroutine calc_grad_on_tile

subroutine syncronize_grad_on_edges(gx,gy,q,mesh)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),          intent(inout) :: gx, gy
    type(tile_edge_transform_t), intent(in)    :: q
    type(tile_mesh_t),           intent(in)    :: mesh

    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    if(is==1 .and. js==1) then
        do k=ks,ke
            gx%p(1,1,k) = (gx%p(1,1,k)+gx%p(0,1,k)+gx%p(1,0,k)+q%qb(1)*gy%p(1,0,k))/3.0_8
            gy%p(1,1,k) = (gy%p(1,1,k)+gy%p(0,1,k)+gy%p(1,0,k)+q%ql(1)*gx%p(0,1,k))/3.0_8
        end do
    end if
    if(ie==mesh%nx+1 .and. js==1) then
        do k=ks,ke
            gx%p(ie,1,k) = (gx%p(ie,1,k)+gx%p(ie+1,1,k)+gx%p(ie,0,k)+q%qb(ie)*gy%p(ie,0,k))/3.0_8
            gy%p(ie,1,k) = (gy%p(ie,1,k)+gy%p(ie+1,1,k)+gy%p(ie,0,k)+q%qr(1)*gx%p(ie+1,1,k))/3.0_8
        end do
    end if
    if(is==1 .and. je==mesh%ny+1) then
        do k=ks,ke
            gx%p(1,je,k) = (gx%p(1,je,k)+gx%p(0,je,k)+gx%p(1,je+1,k)+q%qt(1)*gy%p(1,je+1,k))/3.0_8
            gy%p(1,je,k) = (gy%p(1,je,k)+gy%p(0,je,k)+gy%p(1,je+1,k)+q%ql(je)*gx%p(0,je,k))/3.0_8
        end do
    end if
    if(ie==mesh%nx+1 .and. je==mesh%ny+1) then
        do k=ks,ke
            gx%p(ie,je,k) = (gx%p(ie,je,k)+gx%p(ie+1,je,k)+gx%p(ie,je+1,k)+q%qt(ie)*gy%p(ie,je+1,k))/3.0_8
            gy%p(ie,je,k) = (gy%p(ie,je,k)+gy%p(ie+1,je,k)+gy%p(ie,je+1,k)+q%qr(je)*gx%p(ie+1,je,k))/3.0_8
        end do
    end if
    if(is == 1) then
        do k=ks,ke; do j=max(js,2),min(je,mesh%ny)
            gx%p(1,j,k) = 0.5_8*(gx%p(0,j,k)+gx%p(1,j,k))
            gy%p(1,j,k) = 0.5_8*(gy%p(1,j,k)+gy%p(0,j,k)+q%ql(j)*gx%p(0,j,k))
        end do; end do
    end if
    if(ie == mesh%nx+1) then
        do k=ks,ke; do j=max(js,2),min(je,mesh%ny)
            gx%p(ie,j,k) = 0.5_8*(gx%p(ie,j,k)+gx%p(ie+1,j,k))
            gy%p(ie,j,k) = 0.5_8*(gy%p(ie,j,k)+gy%p(ie+1,j,k)+q%qr(j)*gx%p(ie+1,j,k))
        end do; end do
    end if
    if(js==1) then
        do k=ks,ke; do i=max(is,2),min(ie,mesh%nx)
            gy%p(i,1,k) = 0.5_8*(gy%p(i,0,k)+gy%p(i,1,k))
            gx%p(i,1,k) = 0.5_8*(gx%p(i,0,k)+gx%p(i,1,k)+q%qb(i)*gy%p(i,0,k))
        end do; end do
    end if
    if(je == mesh%ny+1) then
        do k=ks,ke; do i=max(is,2),min(ie,mesh%nx)
            gy%p(i,je,k) = 0.5_8*(gy%p(i,je,k)+gy%p(i,je+1,k))
            gx%p(i,je,k) = 0.5_8*(gx%p(i,je,k)+gx%p(i,je+1,k)+q%qt(i)*gy%p(i,je+1,k))
        end do; end do
    end if
end subroutine syncronize_grad_on_edges

end module grad_contra_ah_sbp_mod
