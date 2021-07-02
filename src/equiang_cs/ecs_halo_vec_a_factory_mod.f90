!Initialization routine for halo_vec object on eq.cubsph grid
!no staggering:
!1)save information about scalar halo-procedure (interpolation of
!  individual components)
!2)precompute vector-component transformation matrices from source panel
!  to target panel
module ecs_halo_vec_a_factory_mod
use ecs_halo_mod,       only : ecs_halo_t
use ecs_halo_vec_a_mod, only : ecs_halo_vec_t, ecs_tile_halo_vec_t

implicit none

integer, parameter :: corner_halo_width = 5!minimum halo-width to compute 2x2 corner-halo-areas

private
public   :: create_ecs_A_vec_halo_procedure

contains

subroutine create_ecs_A_vec_halo_procedure(halo_out,domain,halo_width)
    use halo_mod,               only : halo_vec_t
    use domain_mod,             only : domain_t
    use exchange_factory_mod,   only : create_symm_halo_exchange_A

    class(halo_vec_t), allocatable, intent(out) :: halo_out
    class(domain_t),                intent(in)  :: domain
    integer(kind=4),                intent(in)  :: halo_width

    !locals
    type(ecs_halo_vec_t), allocatable :: halo
    integer(kind=4)      :: ex_halo_width = 8
    integer(kind=4)      :: ts, te, is,ie, js, je, nh, t
    real(kind=8)         :: hx

    allocate(halo)
    ts = domain%partition%ts
    te = domain%partition%te
    nh = domain%partition%nh
    halo%ts = ts
    halo%te = te
    allocate(halo%tile(ts:te))

    halo%exch_halo = create_symm_halo_exchange_A(domain%partition, domain%parcomm, ex_halo_width, 'full')
    do t=ts,te
        hx = domain%mesh_p%tile(t)%hx
        call domain%partition%tile(t)%getind(is,ie,js,je)
        call init_ecs_tile_halo_vect(halo%tile(t),domain%partition%panel_map(t), &
                                     is,ie,js,je,nh,halo_width,hx,domain%topology, &
                                     domain%metric)
    end do

    call move_alloc(halo, halo_out)

end

subroutine init_ecs_tile_halo_vect(tile_halo, panel_ind,is,ie,js,je, &
                              nx,halo_width,hx,topology,metric)

    use ecs_halo_factory_mod, only : init_ecs_tile_halo
    use topology_mod,         only : topology_t
    use metric_mod,           only : metric_t

    integer(kind=4),   intent(in) :: panel_ind,is,ie,js,je,nx
    integer(kind=4),   intent(in) :: halo_width
    real(kind=8),      intent(in) :: hx
    class(topology_t), intent(in) :: topology
    class(metric_t),   intent(in) :: metric

    type(ecs_tile_halo_vec_t), intent(inout) :: tile_halo

    !locals
    !direction to edge in local alpha-beta coordinates:
    integer(kind=4), parameter :: edge_dir(2,4) = [[0,-1], [0,1], [-1,0], [1,0]]
    integer(kind=4) ish, ieh, jsh, jeh
    integer(kind=4) jst, jinc(2), ist, iinc(2) !start points and increments for adjacent panel
                                               ! inc(1) - along edge, inc(2) accross
    integer(kind=4) adjacent_panel_ind, cor_hw

    allocate(tile_halo%scalar_halo)
    call init_ecs_tile_halo(tile_halo%scalar_halo,is,ie,js,je,nx,halo_width,hx)

    tile_halo%n = nx
    tile_halo%lhalo = tile_halo%scalar_halo%lhalo
    tile_halo%panel_ind = panel_ind
    tile_halo%halo_width = max(halo_width,corner_halo_width)
    tile_halo%corner_halo_width = corner_halo_width
    cor_hw = corner_halo_width

    if(tile_halo%lhalo(1) .or. tile_halo%lhalo(2)) then
        ish = max(1,is-halo_width)
        ish = max(1,minval(tile_halo%scalar_halo%indx(ish,1:halo_width)-1))
        ieh = min(nx,ie+halo_width)
        ieh = min(nx,maxval(tile_halo%scalar_halo%indx(ieh,1:halo_width))+2)
        tile_halo%ish = ish; tile_halo%ieh = ieh
    end if
    if(tile_halo%lhalo(3) .or. tile_halo%lhalo(4)) then
        jsh = max(1,js-halo_width)
        jsh = max(1,minval(tile_halo%scalar_halo%indy(jsh,1:halo_width)-1))
        jeh = min(nx,je+halo_width)
        jeh = min(nx,maxval(tile_halo%scalar_halo%indy(jeh,1:halo_width))+2)
        tile_halo%jsh = jsh; tile_halo%jeh = jeh
    end if

    if(tile_halo%lhalo(1)) then
        allocate(tile_halo%TM1(4,ish:ieh,1:max(halo_width,cor_hw)))
        call find_adjacent_panel(adjacent_panel_ind, ist, iinc, jst, jinc, &
                                 panel_ind, edge_dir(1:3,1), nx, topology)
        !print '(8I5)',panel_ind, adjacent_panel_ind, ist, iinc, jst, jinc
        call init_transform_matrix(tile_halo%TM1,ish,ieh,max(halo_width,cor_hw), panel_ind, &
                              adjacent_panel_ind, ist, iinc, jst, jinc, hx,metric)
    end if
    if(tile_halo%lhalo(2)) then
        allocate(tile_halo%TM2(4,ish:ieh,1:max(halo_width,cor_hw)))
        call find_adjacent_panel(adjacent_panel_ind, ist, iinc, jst, jinc, &
                                 panel_ind, edge_dir(1:3,2), nx,topology)
        !print '(8I5)',panel_ind, adjacent_panel_ind, ist, iinc, jst, jinc
        call init_transform_matrix(tile_halo%TM2,ish,ieh,max(halo_width,cor_hw), panel_ind,&
                              adjacent_panel_ind, ist, iinc, jst, jinc, hx,metric)
    end if
    if(tile_halo%lhalo(3)) then
        allocate(tile_halo%TM3(4,jsh:jeh,1:max(halo_width,cor_hw)))
        call find_adjacent_panel(adjacent_panel_ind, ist, iinc, jst, jinc, &
                                 panel_ind, edge_dir(1:3,3), nx,topology)
        !print '(8I5)',panel_ind, adjacent_panel_ind, ist, iinc, jst, jinc
        call init_transform_matrix(tile_halo%TM3,jsh,jeh,max(halo_width,cor_hw), panel_ind, &
                              adjacent_panel_ind, ist, iinc, jst, jinc, hx,metric)
    end if
    if(tile_halo%lhalo(4)) then
        allocate(tile_halo%TM4(4,jsh:jeh,1:max(halo_width,cor_hw)))
        call find_adjacent_panel(adjacent_panel_ind, ist, iinc, jst, jinc, &
                                 panel_ind, edge_dir(1:3,4), nx,topology)
        !print '(8I5)',panel_ind, adjacent_panel_ind, ist, iinc, jst, jinc
        call init_transform_matrix(tile_halo%TM4,jsh,jeh,max(halo_width,cor_hw), panel_ind,&
                              adjacent_panel_ind, ist, iinc, jst, jinc, hx,metric)
    end if
end subroutine init_ecs_tile_halo_vect

subroutine find_adjacent_panel(adjacent_panel_ind, ist, iinc, jst, jinc, &
                               panel_ind, edge_dir, nx, topology)

    use topology_mod, only : topology_t

    integer(kind=4), intent(out) :: adjacent_panel_ind
    integer(kind=4), intent(out) :: ist, iinc(2), jst, jinc(2)
    integer(kind=4), intent(in)  :: panel_ind, edge_dir(2), nx
    class(topology_t), intent(in), target :: topology
    !local
    integer(kind=8) edge_along(2) !direction along edge in local alpha-beta
    integer(kind=4) edge_xyz(3), edge_along_xyz(3)
    integer(kind=4) ee !for scalar prods
    integer(kind=4) pind
    integer(kind=4), dimension(:,:), pointer :: ex, ey, n

    ex => topology%ex
    ey => topology%ey
    n  => topology%n

    edge_xyz(1:3) = edge_dir(1)*ex(1:3,panel_ind)+edge_dir(2)*ey(1:3,panel_ind)
    edge_along(1:2) = [1,1]-abs(edge_dir(1:2))
    edge_along_xyz(1:3) = edge_along(1)*ex(1:3,panel_ind)+edge_along(2)*ey(1:3,panel_ind)

    do pind=1,6
        if(sum(abs(edge_xyz(1:3)+n(1:3,pind)))==0) then
            adjacent_panel_ind = pind
            ee = min(1,max(-1,sum(edge_along_xyz(1:3)*ex(1:3,pind))))
            if(abs(ee)==1) then
                iinc = [ee,0]
                ist  = (1+ee)/2 + (1-ee)/2*nx
                ee = min(1,max(-1,sum(n(1:3,panel_ind)*ey(1:3,pind))))
                jinc = [0, ee]
                jst = (1+ee)/2 + (1-ee)/2*nx
            else
                ee = min(1,max(-1,sum(edge_along_xyz(1:3)*ey(1:3,pind))))
                jinc = [ee,0]
                jst  = (1+ee)/2 + (1-ee)/2*nx
                ee = min(1,max(-1,sum(n(1:3,panel_ind)*ex(1:3,pind))))
                iinc = [0, ee]
                ist = (1+ee)/2 + (1-ee)/2*nx
            end if
            return
        end if
    end do
    call halo_avost("cannot find adjacent face! (init_ecs_htile_alo_vect)")
end subroutine find_adjacent_panel

subroutine init_transform_matrix(TM, i1, i2, hw, panel_ind,    &
                                 adjacent_panel_ind, ist, iinc, jst, jinc, hx,metric)
    use const_mod,        only : pi
    use metric_mod,       only : metric_t
    use ecs_geometry_mod, only : ecs_acov_proto, ecs_bcov_proto, ecs_proto2face,&
                                 ecs_actv_proto, ecs_bctv_proto, ecs_xyz2ab,    &
                                 ecs_ab2xyz_proto

    integer(kind=4), intent(in)  :: i1, i2, hw
    real(kind=8),    intent(out) :: TM(4,i1:i2,1:hw)
    integer(kind=4), intent(in)  :: panel_ind, adjacent_panel_ind
    integer(kind=4), intent(in)  :: ist, iinc(2), jst, jinc(2)
    real(kind=8),    intent(in)  :: hx
    class(metric_t), intent(in)  :: metric
    !locals
    integer(kind=4) i, j, ia, ja
    real(kind=8) as(3), bs(3), at(3), bt(3)
    real(kind=8) alpha, beta, r(3) !alpha, beta and xyz-cartesian

    do j=1, hw
        do i=i1,i2
            !indices at source panel to get source-alpha|beta
            ia = ist+iinc(1)*(i-1)+iinc(2)*(j-1)
            ja = jst+jinc(1)*(i-1)+jinc(2)*(j-1)
            alpha = -0.25_8*pi+(ia-0.5_8)*hx
            beta  = -0.25_8*pi+(ja-0.5_8)*hx
            !covariant basis vectors at source panel
            as = ecs_acov_proto(alpha, beta)
            as = ecs_proto2face(as,adjacent_panel_ind)
            bs = ecs_bcov_proto(alpha, beta)
            bs = ecs_proto2face(bs,adjacent_panel_ind)
            !xyz of source panel point to get alpha/beta at target panel
            r = ecs_ab2xyz_proto(alpha, beta)
            r = ecs_proto2face(r,adjacent_panel_ind)
            call ecs_xyz2ab(alpha,beta,r(1),r(2),r(3), panel_ind)
            at = ecs_actv_proto(alpha, beta)
            at = ecs_proto2face(at,panel_ind)
            bt = ecs_bctv_proto(alpha, beta)
            bt = ecs_proto2face(bt,panel_ind)
            !transform matrix <- dot-products of target_contra_vec*source_cov_vec
            !TM = |1  2|
            !     |3  4|
            TM(1,i,j) = sum(as(1:3)*at(1:3))
            TM(2,i,j) = sum(bs(1:3)*at(1:3))
            TM(3,i,j) = sum(as(1:3)*bt(1:3))
            TM(4,i,j) = sum(bs(1:3)*bt(1:3))
        end do
    end do

end subroutine init_transform_matrix

subroutine halo_avost(str)
use mpi
integer ierr
character(*) str
print *, str
print *, "exit"
call mpi_finalize(ierr)
stop
end subroutine halo_avost

end module ecs_halo_vec_a_factory_mod
