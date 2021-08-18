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
    integer(kind=4)      :: ex_halo_width = 7
    integer(kind=4)      :: ts, te, is, ie, js, je, nh, t
    integer(kind=4)      :: is1, ie1, js1, je1
    real(kind=8)         :: hx
    character(:), allocatable :: normal_comp_interp_type, tangential_comp_interp_type
    character(:), allocatable :: normal_tang_interp_type

    allocate(halo)
    ts = domain%partition%ts
    te = domain%partition%te
    halo%ts = ts
    halo%te = te
    nh = domain%partition%nh

    halo%exch_halo = create_symmetric_halo_vec_exchange_C(domain%partition, &
                              domain%parcomm, domain%topology, ex_halo_width, 'full')

    halo%exchange_width = ex_halo_width
    allocate(halo%tile(ts:te))

    normal_comp_interp_type = "cubic_lag"
    tangential_comp_interp_type = "cubic_lag"
    normal_tang_interp_type = "cubic_lag"
    do t=ts,te
        is = domain%partition%tile_x(t)%is
        ie = domain%partition%tile_x(t)%ie
        js = domain%partition%tile_x(t)%js
        je = domain%partition%tile_x(t)%je
        js1 = domain%partition%tile_y(t)%js
        je1 = domain%partition%tile_y(t)%je

        halo%tile(t)%is_left_edge   = (is == 1)
        if(halo%tile(t)%is_left_edge) then
            call init_1d_tile_interp(halo%tile(t)%wtn_left,halo%tile(t)%itn_left,    &
                                         js1,je1,halo_width,nh,domain%mesh_v%tile(t)%hx, &
                                         normal_tang_interp_type,"y","x")
            call init_basis_transform(halo%tile(t)%qtt_left,halo%tile(t)%qtn_left,js1,je1,&
                                      halo_width,domain%mesh_v%tile(t)%hx)
            call init_1d_tile_interp(halo%tile(t)%wt_left,halo%tile(t)%it_left,    &
                                         js1,je1,halo_width,nh+1,domain%mesh_v%tile(t)%hx, &
                                         tangential_comp_interp_type,"x","x")
            call init_1d_tile_interp(halo%tile(t)%wn_left,halo%tile(t)%in_left,    &
                                         js,je,halo_width,nh,domain%mesh_u%tile(t)%hx, &
                                         normal_comp_interp_type,"y","y")
        end if
        
        halo%tile(t)%is_right_edge  = (ie == nh+1)
        if(halo%tile(t)%is_right_edge) then
            call init_1d_tile_interp(halo%tile(t)%wtn_right,halo%tile(t)%itn_right,    &
                                         js1,je1,halo_width,nh,domain%mesh_v%tile(t)%hx, &
                                         normal_tang_interp_type,"y","x")
            call init_basis_transform(halo%tile(t)%qtt_right,halo%tile(t)%qtn_right,js1,je1,&
                                      halo_width,domain%mesh_v%tile(t)%hx)
            call init_1d_tile_interp(halo%tile(t)%wt_right,halo%tile(t)%it_right,    &
                                     js1,je1,halo_width,nh+1,domain%mesh_v%tile(t)%hx, &
                                     tangential_comp_interp_type,"x","x")
            call init_1d_tile_interp(halo%tile(t)%wn_right,halo%tile(t)%in_right,  &
                                         js,je,halo_width,nh,domain%mesh_u%tile(t)%hx, &
                                         normal_comp_interp_type,"y","y")
        end if

        is = domain%partition%tile_y(t)%is
        ie = domain%partition%tile_y(t)%ie
        js = domain%partition%tile_y(t)%js
        je = domain%partition%tile_y(t)%je
        is1 = domain%partition%tile_x(t)%is
        ie1 = domain%partition%tile_x(t)%ie

        halo%tile(t)%is_bottom_edge = (js == 1)
        if(halo%tile(t)%is_bottom_edge) then
            call init_1d_tile_interp(halo%tile(t)%wtn_bottom,halo%tile(t)%itn_bottom,    &
                                         is1,ie1,halo_width,nh,domain%mesh_u%tile(t)%hx, &
                                         normal_tang_interp_type,"y","x")
            call init_basis_transform(halo%tile(t)%qtt_bottom,halo%tile(t)%qtn_bottom,is1,ie1,&
                                      halo_width,domain%mesh_u%tile(t)%hx)
            call init_1d_tile_interp(halo%tile(t)%wt_bottom,halo%tile(t)%it_bottom,    &
                                     is1,ie1,halo_width,nh+1,domain%mesh_u%tile(t)%hx, &
                                     tangential_comp_interp_type,"x","x")
            call init_1d_tile_interp(halo%tile(t)%wn_bottom,halo%tile(t)%in_bottom, &
                                         is,ie,halo_width,nh,domain%mesh_v%tile(t)%hx,  &
                                         normal_comp_interp_type,"y","y")
        end if

        halo%tile(t)%is_top_edge    = (je == nh+1)
        if(halo%tile(t)%is_top_edge) then
            call init_1d_tile_interp(halo%tile(t)%wtn_top,halo%tile(t)%itn_top,    &
                                        is1,ie1,halo_width,nh,domain%mesh_u%tile(t)%hx, &
                                        normal_tang_interp_type,"y","x")
            call init_basis_transform(halo%tile(t)%qtt_top,halo%tile(t)%qtn_top,is1,ie1,&
                                      halo_width,domain%mesh_u%tile(t)%hx)
            call init_1d_tile_interp(halo%tile(t)%wt_top,halo%tile(t)%it_top,    &
                                     is1,ie1,halo_width,nh+1,domain%mesh_u%tile(t)%hx, &
                                     tangential_comp_interp_type,"x","x")
            call init_1d_tile_interp(halo%tile(t)%wn_top,halo%tile(t)%in_top,      &
                                         is,ie,halo_width,nh,domain%mesh_v%tile(t)%hx, &
                                         normal_comp_interp_type,"y","y")
        end if

    end do

    call move_alloc(halo, halo_out)

end

subroutine init_1d_tile_interp(w,ind,is,ie,halo_width,nx,hx,&
                               interp_type,source_point_type,target_point_type)
    use const_mod,   only : pi
    use parcomm_mod, only : parcomm_global

    real(kind=8),    allocatable, intent(out) :: w(:,:,:)
    integer(kind=4), allocatable, intent(out) :: ind(:,:)
    integer(kind=4),              intent(in)  :: is, ie, halo_width,nx
    real(kind=8),                 intent(in)  :: hx !grid spacing
    character(len=*),             intent(in)  :: interp_type
    character(len=*),             intent(in)  :: source_point_type
    character(len=*),             intent(in)  :: target_point_type

    integer(kind=4) :: i, j
    real(kind=8)    :: alpha, beta, xi, xi1
    real(kind=8)    :: sx, sy, tx, ty !source and target grid-shift
    real(kind=8)    :: x,y,z

    select case(source_point_type)
    case("y")
        sx = 0.5_8
        sy = 0.0_8
    case("x")
        sx = 0.0_8
        sy =-0.5_8
    case default
        call parcomm_global%abort("ecs_halo_vec_c_factory_mod, init_tile_interp - unknown point type: "//source_point_type)
    end select

    select case(target_point_type)
    case("y")
        tx = 0.5_8
        ty = 0.0_8
    case("x")
        tx = 0.0_8
        ty =-0.5_8
    case default
        call parcomm_global%abort("ecs_halo_vec_c_factory_mod, init_tile_interp - unknown point type: "//source_point_type)
    end select

    allocate(ind(is:ie,1:halo_width))
    select case(interp_type)
    case("linear")
        allocate(w(0:1,is:ie,1:halo_width))
    case("cubic_lag")
        allocate(w(-1:2,is:ie,1:halo_width))
    case default
        call parcomm_global%abort("ecs_halo_vec_c_factory_mod, init_tile_interp - unknown interp type: "//interp_type)
    end select

    do j=1, halo_width
        beta = 0.25_8*pi+(j+ty)*hx
        select case(interp_type)
        case("linear")
            do i=is,ie
                alpha = -0.25_8*pi+(i-1.0_8+tx)*hx
                x = tan(alpha)
                y = tan(beta)
                z = 1.0_8 / sqrt(1.0_8+x**2+y**2)
                x = x*z
                y = y*z
                alpha = atan(x/y)
                xi = (alpha+0.25_8*pi-sx*hx) / hx
                ind(i,j) = int(xi)+1
                xi = xi - ind(i,j)+1
                w(0:1,i,j) = [1._8-xi, xi]
            end do
        case("cubic_lag")
            do i=is,ie
                alpha = -0.25_8*pi+(i-1.0_8+tx)*hx
                x = tan(alpha)
                y = tan(beta)
                z = 1.0_8 / sqrt(1.0_8+x**2+y**2)
                x = x*z
                y = y*z
                alpha = atan(x/y)
                xi = (alpha+0.25_8*pi-sx*hx) / hx
                ind(i,j) = max(2,min(int(xi)+1,nx-2))
                xi = xi - ind(i,j)+1
                w(-1:2,i,j) = [-xi*(xi-1.0_8)*(xi-2.0_8) / 6.0_8,         &
                                (xi+1.0_8)*(xi-1.0_8)*(xi-2.0_8) / 2.0_8, &
                               -(xi+1.0_8)*xi*(xi-2.0_8) / 2.0_8,         &
                                (xi+1.0_8)*xi*(xi-1.0_8) / 6.0_8]
            end do
        end select
    end do
    
end subroutine init_1d_tile_interp

subroutine init_basis_transform(qtt,qtn,si,ei,halo_width,hx)
!coefficients for the following transform:
!target face tangential component = qtt * source face tangential + qtn * source face normal

use const_mod,      only : pi
use ecs_metric_mod, only : ecs_b1_proto, ecs_a1_proto, ecs_a2_proto

real(kind=8), allocatable, intent(out) :: qtt(:,:), qtn(:,:)
integer(kind=4),           intent(in)  :: si, ei, halo_width
real(kind=8),              intent(in)  :: hx

integer(kind=4)    :: i, j
real(kind=8)       :: alpha, beta, alpha1, beta1
real(kind=8)       :: x, y, z
real(kind=8)       :: bt(3), an(3), at(3), v(3)

allocate(qtt(si:ei,halo_width), qtn(si:ei,halo_width))

do j=1, halo_width
    beta = 0.25_8*pi+(j-0.5_8)*hx
    do i=si,ei
        alpha = -0.25_8*pi+(i-1.0_8)*hx
        x = tan(alpha)
        y = tan(beta)
        z = 1.0_8 / sqrt(1.0_8+x**2+y**2)
        x = x*z
        y = y*z
        alpha1 = atan(x/y)
        beta1  = -0.25_8*pi+(j-0.5_8)*hx
        bt = ecs_b1_proto(alpha, beta)
        v  = ecs_a1_proto(alpha1,beta1)
        at = [v(1),v(3),-v(2)] 
        v  = ecs_a2_proto(alpha1,beta1)
        an = [v(1),v(3),-v(2)]
        qtt(i,j) = sum(at(1:3)*bt(1:3)) 
        qtn(i,j) = sum(an(1:3)*bt(1:3))
    end do
end do

end subroutine init_basis_transform

end module ecs_halo_vec_c_factory_mod
