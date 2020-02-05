!Initialization routine for halo_vec object on eq.cubsph grid
!no staggering:
!1)save information about scalar halo-procedure (interpolation of
!  individual components)
!2)precompute vector-component transformation matrices from source panel
!  to target panel
module ecs_halo_vec_a_factory_mod
use ecs_halo_mod,       only : ecs_halo_t
use ecs_halo_vec_a_mod, only : ecs_halo_vec_t

implicit none

integer, parameter :: corner_halo_width = 5!minimum halo-width to compute 2x2 corner-halo-areas

private
public   :: init_ecs_halo_vect

contains

type(ecs_halo_vec_t) function init_ecs_halo_vect(panel_ind,is,ie,js,je, &
                                            nx,halo_width,hx,halo) result(halo_vec)

    integer(kind=4),   intent(in) :: panel_ind,is,ie,js,je,nx
    integer(kind=4),   intent(in) :: halo_width
    real(kind=8),      intent(in) :: hx
    class(ecs_halo_t), intent(in) :: halo

    !locals
    !direction to edge in local alpha-beta coordinates:
    integer(kind=4), parameter :: edge_dir(2,4) = [[0,-1], [0,1], [-1,0], [1,0]]
    integer(kind=4) ish, ieh, jsh, jeh
    integer(kind=4) jst, jinc(2), ist, iinc(2) !start points and increments for adjacent panel
                                               ! inc(1) - along edge, inc(2) accross
    integer(kind=4) adjacent_panel_ind, cor_hw

    halo_vec%n = nx
    halo_vec%halo = halo !implicit allocate here
    halo_vec%lhalo = halo%lhalo
    halo_vec%panel_ind = panel_ind
    halo_vec%halo_width = max(halo_width,corner_halo_width)
    halo_vec%corner_halo_width = corner_halo_width
    cor_hw = corner_halo_width

    if(halo_vec%lhalo(1) .or. halo_vec%lhalo(2)) then
        ish = max(1,is-halo_width)
        ish = max(1,minval(halo%indx(ish,1:halo_width)-1))
        ieh = min(nx,ie+halo_width)
        ieh = min(nx,maxval(halo%indx(ieh,1:halo_width))+2)
        halo_vec%ish = ish; halo_vec%ieh = ieh
    end if
    if(halo_vec%lhalo(3) .or. halo_vec%lhalo(4)) then
        jsh = max(1,js-halo_width)
        jsh = max(1,minval(halo%indy(jsh,1:halo_width)-1))
        jeh = min(nx,je+halo_width)
        jeh = min(nx,maxval(halo%indy(jeh,1:halo_width))+2)
        halo_vec%jsh = jsh; halo_vec%jeh = jeh
    end if

    if(halo_vec%lhalo(1)) then
        allocate(halo_vec%TM1(4,ish:ieh,1:max(halo_width,cor_hw)))
        call find_adjacent_panel(adjacent_panel_ind, ist, iinc, jst, jinc, &
                                 panel_ind, edge_dir(1:3,1), nx)
        !print '(8I5)',panel_ind, adjacent_panel_ind, ist, iinc, jst, jinc
        call init_transform_matrix(halo_vec%TM1,ish,ieh,max(halo_width,cor_hw), panel_ind, &
                              adjacent_panel_ind, ist, iinc, jst, jinc, hx)
    end if
    if(halo_vec%lhalo(2)) then
        allocate(halo_vec%TM2(4,ish:ieh,1:max(halo_width,cor_hw)))
        call find_adjacent_panel(adjacent_panel_ind, ist, iinc, jst, jinc, &
                                 panel_ind, edge_dir(1:3,2), nx)
        !print '(8I5)',panel_ind, adjacent_panel_ind, ist, iinc, jst, jinc
        call init_transform_matrix(halo_vec%TM2,ish,ieh,max(halo_width,cor_hw), panel_ind,&
                              adjacent_panel_ind, ist, iinc, jst, jinc, hx)
    end if
    if(halo_vec%lhalo(3)) then
        allocate(halo_vec%TM3(4,jsh:jeh,1:max(halo_width,cor_hw)))
        call find_adjacent_panel(adjacent_panel_ind, ist, iinc, jst, jinc, &
                                 panel_ind, edge_dir(1:3,3), nx)
        !print '(8I5)',panel_ind, adjacent_panel_ind, ist, iinc, jst, jinc
        call init_transform_matrix(halo_vec%TM3,jsh,jeh,max(halo_width,cor_hw), panel_ind, &
                              adjacent_panel_ind, ist, iinc, jst, jinc, hx)
    end if
    if(halo_vec%lhalo(4)) then
        allocate(halo_vec%TM4(4,jsh:jeh,1:max(halo_width,cor_hw)))
        call find_adjacent_panel(adjacent_panel_ind, ist, iinc, jst, jinc, &
                                 panel_ind, edge_dir(1:3,4), nx)
        !print '(8I5)',panel_ind, adjacent_panel_ind, ist, iinc, jst, jinc
        call init_transform_matrix(halo_vec%TM4,jsh,jeh,max(halo_width,cor_hw), panel_ind,&
                              adjacent_panel_ind, ist, iinc, jst, jinc, hx)
    end if
end function init_ecs_halo_vect

subroutine find_adjacent_panel(adjacent_panel_ind, ist, iinc, jst, jinc, &
                               panel_ind, edge_dir, nx)
    use topology_mod, only : ex, ey, n
    integer(kind=4), intent(out) :: adjacent_panel_ind
    integer(kind=4), intent(out) :: ist, iinc(2), jst, jinc(2)
    integer(kind=4), intent(in)  :: panel_ind, edge_dir(2), nx
    !local
    integer(kind=8) edge_along(2) !direction along edge in local alpha-beta
    integer(kind=4) edge_xyz(3), edge_along_xyz(3)
    integer(kind=4) ee !for scalar prods
    integer(kind=4) pind

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
    call halo_avost("cannot find adjacent face! (init_ecs_halo_vect)")
end subroutine find_adjacent_panel

subroutine init_transform_matrix(TM, i1, i2, hw, panel_ind,    &
                                 adjacent_panel_ind, ist, iinc, jst, jinc, hx)
    use const_mod,        only : pi
    use ecs_geometry_mod, only : ecs_acov_proto, ecs_bcov_proto, ecs_proto2face,&
                                 ecs_actv_proto, ecs_bctv_proto, ecs_xyz2ab,    &
                                 ecs_ab2xyz_proto

    integer(kind=4), intent(in)  :: i1, i2, hw
    real(kind=8),    intent(out) :: TM(4,i1:i2,1:hw)
    integer(kind=4), intent(in)  :: panel_ind, adjacent_panel_ind
    integer(kind=4), intent(in)  :: ist, iinc(2), jst, jinc(2)
    real(kind=8),    intent(in)  :: hx
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
!print '(6(A,i8,1x))', "is=", is, "ie=", ie, "js=", js, "je=", je, "nvi=", iev-ie, "nvj=", jev-je
print *, "exit"
call mpi_finalize(ierr)
stop
end subroutine halo_avost

end module ecs_halo_vec_a_factory_mod
