!Initialization routine for halo_vec object on eq.cubsph grid
!C-type staggering:
module ecs_halo_vec_c_factory_mod
use ecs_halo_mod,       only : ecs_halo_t
use ecs_halo_vec_c_mod, only : ecs_halo_vec_c_t
use parcomm_mod,        only : parcomm_global

implicit none

!integer, parameter :: corner_halo_width = 5!minimum halo-width to compute 2x2 corner-halo-areas

private
public   :: create_ecs_C_vec_halo_procedure

contains

subroutine create_ecs_C_vec_halo_procedure(halo_out,domain,halo_width)
    use halo_mod,               only : halo_vec_t
    use domain_mod,             only : domain_t
    use exchange_factory_mod,   only : create_symmetric_halo_vec_exchange_C

    class(halo_vec_t), allocatable, intent(out) :: halo_out
    class(domain_t),                intent(in)  :: domain
    integer(kind=4),                intent(in)  :: halo_width

    !locals
    type(ecs_halo_vec_c_t), allocatable :: halo
    integer(kind=4)      :: ex_halo_width = 8
    integer(kind=4)      :: ts, te, is,ie, js, je, nh, t
    real(kind=8)         :: hx
    character(:), allocatable :: normal_comp_interp_type

    allocate(halo)
    ts = domain%partition%ts
    te = domain%partition%te
    halo%ts = ts
    halo%te = te
    nh = domain%partition%nh

    halo%exch_halo = create_symmetric_halo_vec_exchange_C(domain%partition, &
                              domain%parcomm, domain%topology, halo_width, 'full')

    allocate(halo%tile(ts:te))

    normal_comp_interp_type = "cubic_lag"
    do t=ts,te
        is = domain%partition%tile_u(t)%is
        ie = domain%partition%tile_u(t)%ie
        js = domain%partition%tile_u(t)%js
        je = domain%partition%tile_u(t)%je

        halo%tile(t)%is_left_edge   = (is == 1)
        if(halo%tile(t)%is_left_edge) then
            call init_tile_normal_interp(halo%tile(t)%wn_left,halo%tile(t)%in_left,    &
                                         js,je,halo_width,nh,domain%mesh_u%tile(t)%hx, &
                                         normal_comp_interp_type)
        end if
        
        halo%tile(t)%is_right_edge  = (ie == nh+1)
        if(halo%tile(t)%is_right_edge) then
            call init_tile_normal_interp(halo%tile(t)%wn_right,halo%tile(t)%in_right,  &
                                         js,je,halo_width,nh,domain%mesh_u%tile(t)%hx, &
                                         normal_comp_interp_type)
        end if

        is = domain%partition%tile_v(t)%is
        ie = domain%partition%tile_v(t)%ie
        js = domain%partition%tile_v(t)%js
        je = domain%partition%tile_v(t)%je

        halo%tile(t)%is_bottom_edge = (js == 1)
        if(halo%tile(t)%is_bottom_edge) then
            call init_tile_normal_interp(halo%tile(t)%wn_bottom,halo%tile(t)%in_bottom, &
                                         is,ie,halo_width,nh,domain%mesh_v%tile(t)%hx,  &
                                         normal_comp_interp_type)
        end if

        halo%tile(t)%is_top_edge    = (je == nh+1)
        if(halo%tile(t)%is_top_edge) then
            call init_tile_normal_interp(halo%tile(t)%wn_top,halo%tile(t)%in_top,      &
                                         is,ie,halo_width,nh,domain%mesh_v%tile(t)%hx, &
                                         normal_comp_interp_type)
        end if

    end do

    call move_alloc(halo, halo_out)

end

subroutine init_tile_normal_interp(w,ind,is,ie,halo_width,nx,hx,interp_type)
    use const_mod,   only : pi
    use parcomm_mod, only : parcomm_global

    real(kind=8),    allocatable, intent(out) :: w(:,:,:)
    integer(kind=4), allocatable, intent(out) :: ind(:,:)
    integer(kind=4),              intent(in)  :: is, ie, halo_width,nx
    real(kind=8),                 intent(in)  :: hx
    character(len=*),             intent(in)  :: interp_type

    integer(kind=4) :: i, j
    real(kind=8)    :: alpha, beta, xi, xi1
    real(kind=8)    :: x,y,z

    allocate(ind(is:ie,1:halo_width))
    select case(interp_type)
    case("linear")
        allocate(w(0:1,is:ie,1:halo_width))
    case("cubic_lag")
        allocate(w(-1:2,is:ie,1:halo_width))
    case default
        call parcomm_global%abort("ecs_halo_vec_c_factory_mod, init_tile_normal_interp - unknown interp type: "//interp_type)
    end select

    do j=1, halo_width
        beta = 0.25_8*pi+j*hx
        select case(interp_type)
        case("linear")
            do i=is,ie
                alpha = -0.25_8*pi+(i-0.5_8)*hx
                x = tan(alpha)
                y = tan(beta)
                z = 1.0_8 / sqrt(1.0_8+x**2+y**2)
                x = x*z
                y = y*z
                alpha = atan(x/y)
                xi = (alpha+0.25_8*pi-0.5_8*hx) / hx
                ind(i,j) = int(xi)+1
                xi = xi - ind(i,j)+1
                w(0:1,i,j) = [1._8-xi, xi]
            end do
        case("cubic_lag")
            do i=is,ie
                alpha = -0.25_8*pi+(i-0.5_8)*hx
                x = tan(alpha)
                y = tan(beta)
                z = 1.0_8 / sqrt(1.0_8+x**2+y**2)
                x = x*z
                y = y*z
                alpha = atan(x/y)
                xi = (alpha+0.25_8*pi-0.5_8*hx) / hx
                !ind(i,j) = int(xi)+1
                !xi = xi - ind(i,j)+1
                !w(-1:2,i,j) = [0.0_8, 1._8-xi, xi, 0.0_8]
                ind(i,j) = max(2,min(int(xi)+1,nx-2))
                xi = xi - ind(i,j)+1
                w(-1:2,i,j) = [-xi*(xi-1.0_8)*(xi-2.0_8) / 6.0_8,         &
                                (xi+1.0_8)*(xi-1.0_8)*(xi-2.0_8) / 2.0_8, &
                               -(xi+1.0_8)*xi*(xi-2.0_8) / 2.0_8,         &
                                (xi+1.0_8)*xi*(xi-1.0_8) / 6.0_8]
            end do
        end select
    end do
    
end subroutine init_tile_normal_interp

end module ecs_halo_vec_c_factory_mod
