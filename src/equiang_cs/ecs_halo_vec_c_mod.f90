!Object to interpolate vector values to virtual points beyond face edge
!C-staggered grid.

module ecs_halo_vec_c_mod

use halo_mod,          only : halo_vec_t
use exchange_halo_mod, only : exchange_t

implicit none

type, extends(halo_vec_t) :: ecs_halo_vec_c_t

    integer(kind=4)                         :: ts, te
    class(exchange_t),          allocatable :: exch_halo
    type(ecs_tile_halo_cvec_t), allocatable :: tile(:)

    contains

    procedure :: get_halo_vector => get_ecs_c_vector_halo

end type

type ecs_tile_halo_cvec_t
    logical is_left_edge, is_right_edge, is_bottom_edge, is_top_edge
    !interpolation weights and stencils for normal components
    real(kind=8),    allocatable :: wn_left(:,:,:), wn_right(:,:,:), wn_bottom(:,:,:), wn_top(:,:,:)
    integer(kind=4), allocatable :: in_left(:,:), in_right(:,:), in_bottom(:,:), in_top(:,:)
    contains

    procedure :: interpv => interp_ecs_tile_halo_cvec
end type ecs_tile_halo_cvec_t

contains

subroutine get_ecs_c_vector_halo(this,u,v,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(ecs_halo_vec_c_t),  intent(inout) :: this
    class(grid_field_t),      intent(inout) :: u,v
    type(domain_t),           intent(in)    :: domain
    integer(kind=4),          intent(in)    :: halo_width

    integer(kind=4) t

    call this%exch_halo%do_vec(u,v, domain%parcomm)

    do t=this%ts,this%te
        call this%tile(t)%interpv(u%tile(t),v%tile(t),domain%partition%tile_u(t), &
                                  domain%partition%tile_v(t),halo_width)
    end do
end subroutine get_ecs_c_vector_halo

subroutine interp_ecs_tile_halo_cvec(this,u,v,tile_u,tile_v,halo_width)
    use grid_field_mod, only: tile_field_t
    use tile_mod,       only: tile_t

    class(ecs_tile_halo_cvec_t), intent(in)    :: this
    type(tile_field_t),          intent(inout) :: u, v
    type(tile_t),                intent(in)    :: tile_u, tile_v
    integer(kind=4),             intent(in)    :: halo_width

    integer(kind=4) :: i,is,ie,js,je,ks,ke

    is = tile_u%is
    ie = tile_u%ie
    js = tile_u%js
    je = tile_u%je
    ks = tile_u%ks
    ke = tile_u%ke

    if(this%is_left_edge) then
        u%p(1,js:je,ks:ke) = 0.5_8*(u%p(1,js:je,ks:ke)+u%p(0,js:je,ks:ke))
        do i=1,halo_width
            u%p(1-i,js:je,ks:ke) = u%p(-i,js:je,ks:ke)
        end do
        call interp_edge(u,this%wn_left,this%in_left,js,je,ks,ke,halo_width,ie,'left')
    end if
    if(this%is_right_edge) then
        u%p(ie,js:je,ks:ke) = 0.5_8*(u%p(ie,js:je,ks:ke)+u%p(ie+1,js:je,ks:ke))
        do i=1,halo_width
            u%p(ie+i,js:je,ks:ke) = u%p(ie+i+1,js:je,ks:ke)
        end do
        call interp_edge(u,this%wn_right,this%in_right,js,je,ks,ke,halo_width,ie,'right')
    end if

    is = tile_v%is
    ie = tile_v%ie
    js = tile_v%js
    je = tile_v%je
    ks = tile_v%ks
    ke = tile_v%ke

    if(this%is_bottom_edge) then
        v%p(is:ie,1,ks:ke) = 0.5_8*(v%p(is:ie,1,ks:ke)+v%p(is:ie,0,ks:ke))
        do i=1,halo_width
            v%p(is:ie,1-i,ks:ke) = v%p(is:ie,-i,ks:ke)
        end do
        call interp_edge(v,this%wn_bottom,this%in_bottom,is,ie,ks,ke,halo_width,je,'bottom')
    end if
    if(this%is_top_edge) then
        v%p(is:ie,je,ks:ke) = 0.5_8*(v%p(is:ie,je,ks:ke)+v%p(is:ie,je+1,ks:ke))
        do i=1,halo_width
            v%p(is:ie,je+i,ks:ke) = v%p(is:ie,je+i+1,ks:ke)
        end do
        call interp_edge(v,this%wn_top,this%in_top,is,ie,ks,ke,halo_width,je,'top')
    end if

end subroutine interp_ecs_tile_halo_cvec

subroutine interp_edge(f, w, ind, is, ie, ks, ke, halo_width, nh, edge)
    use grid_field_mod, only : tile_field_t
    use parcomm_mod,    only : parcomm_global

    type(tile_field_t), intent(inout) :: f
    integer(kind=4),    intent(in)    :: is, ie, ks, ke, halo_width, nh
    real(kind=8),       intent(in)    :: w(-1:2,is:ie,1:halo_width)
    integer(kind=4),    intent(in)    :: ind(is:ie,1:halo_width)
    character(len=*),   intent(in)    :: edge

    integer(kind=4) :: i, j, k
    real(kind=8)    :: f1(is:ie)

    do k=ks,ke
        do j=1,halo_width
            if(edge=="left") then
                do i=is,ie
                    f1(i) = w(0,i,j)*f%p(1-j,ind(i,j),k)+w(1,i,j)*f%p(1-j,ind(i,j)+1,k)
                end do
                f%p(1-j,is:ie,k) = f1(is:ie)
            elseif(edge=="right") then
                do i=is,ie
                    f1(i) = w(0,i,j)*f%p(nh+j,ind(i,j),k)+w(1,i,j)*f%p(nh+j,ind(i,j)+1,k)
                end do
                f%p(nh+j,is:ie,k) = f1(is:ie)
            elseif(edge=="bottom") then
                do i=is,ie
                    f1(i) = w(0,i,j)*f%p(ind(i,j),1-j,k)+w(1,i,j)*f%p(ind(i,j)+1,1-j,k)
                end do
                f%p(is:ie,1-j,k) = f1(is:ie)
            elseif(edge=="top") then
                do i=is,ie
                    f1(i) = w(0,i,j)*f%p(ind(i,j),nh+j,k)+w(1,i,j)*f%p(ind(i,j)+1,nh+j,k)
                end do
                f%p(is:ie,nh+j,k) = f1(is:ie)
            else
                call parcomm_global%abort("unknown edge "//edge)
            end if
        end do
    end do

end subroutine interp_edge

end module ecs_halo_vec_c_mod
