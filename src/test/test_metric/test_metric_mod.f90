module test_metric_mod
implicit none
private
public :: test_metric
contains

subroutine test_metric()

use mpi

use grid_function_mod,     only : grid_function_t
use exchange_abstract_mod, only : exchange_t
use partition_mod,         only : partition_t
use exchange_factory_mod,  only : create_Agrid_halo_exchange
use mesh_factory_mod,      only : create_equiangular_mesh
use mesh_mod,              only : mesh_t

class(exchange_t),     allocatable :: exch_halo
type(partition_t)                  :: partition
type(grid_function_t), allocatable :: f1(:)
type(grid_function_t), allocatable :: f2(:)
type(mesh_t),          allocatable :: mesh(:)

integer(kind=4), parameter         :: nh=16, nz=3, halo_width=1
real(kind=8),    parameter         :: zeps = 1e-16_8
integer(kind=4)                    :: myid, np, ierr, code

integer(kind=4) :: ts, te
integer(kind=4) :: ind, i, j, k, ifc, err_sum, gl_err_sum
integer(kind=4) :: ie, is, js, je

integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number

character(:), allocatable :: frmt
logical lcross_edge_xyz, lsymm_check
real(kind=8) zq(1:nh,3), zp(1:nh,3)
real(kind=8) err_ort_rab, err_ort_abt, err_q, err_qi, err_g

call MPI_comm_rank(mpi_comm_world , myid, ierr)
call MPI_comm_size(mpi_comm_world , Np  , ierr)

if (myid==0) print *, 'equiangular cubed-sphere test'

if(Np >1) then
    print *, "this test currently runs only with 1 mpi-process. exiting"
    return
end if

call mpi_barrier(mpi_comm_world, ierr)

call partition%init(nh, nz, max(1,Np/6), myid, Np, strategy = 'default')

!find start and end index of tiles belonging to the current proccesor
ts = findloc(partition%proc_map, myid, dim=1)
te = findloc(partition%proc_map, myid, back = .true., dim=1)

allocate(mesh(ts:te))
lsymm_check = .true.
do ind = ts, te
    call create_equiangular_mesh(mesh(ind), partition%tile(ind)%is, partition%tile(ind)%ie, &
                                            partition%tile(ind)%js, partition%tile(ind)%je, &
                                            partition%tile(ind)%ks, partition%tile(ind)%ke, &
                                            nh, halo_width, partition%tile(ind)%panel_number)
    lsymm_check = lsymm_check .and. symmetricity_check_h(mesh(ind)%rhx,mesh(ind)%rhy,mesh(ind)%rhz, &
                                                         nh+2*halo_width,mesh(ind)%panel_ind)
end do

if(lsymm_check) then
    print *, "symmetricity test passed"
else
    print *, "symmetricity test failed"
end if

!Init arrays

allocate(f1(ts:te))
allocate(f2(ts:te))

do i = ts, te

    call f1(i)%init(partition%tile(i)%panel_number,             &
                    partition%tile(i)%is, partition%tile(i)%ie, &
                    partition%tile(i)%js, partition%tile(i)%je, &
                    partition%tile(i)%ks, partition%tile(i)%ke, &
                    halo_width, halo_width, 0)
    call f2(i)%init(partition%tile(i)%panel_number,             &
                    partition%tile(i)%is, partition%tile(i)%ie, &
                    partition%tile(i)%js, partition%tile(i)%je, &
                    partition%tile(i)%ks, partition%tile(i)%ke, &
                    halo_width, halo_width, 0)
end do

do ind = ts, te
     ifc = partition%tile(ind)%panel_number
     f1(ind).p(:,:,1) = mesh(ind)%rhx(:,:)
     f1(ind).p(:,:,2) = mesh(ind)%rhy(:,:)
     f1(ind).p(:,:,3) = mesh(ind)%rhz(:,:)
     f2(ind).p(:,:,:) = f1(ind).p(:,:,:)
end do

!Init exchange
exch_halo = create_Agrid_halo_exchange(partition, halo_width, 'full', myid, np)

!Perform exchange
call exch_halo%do(f1, lbound(f1, 1), ubound(f1, 1))

lcross_edge_xyz = .true.
do ind = ts, te
    zq = f1(ind)%p(1:nh,0,:); zp = f2(ind)%p(1:nh,0,:)
    lcross_edge_xyz = lcross_edge_xyz .and. cross_edge_xyz_check(zq,zp,nh)
    zq = f1(ind)%p(1:nh,nh+1,:); zp = f2(ind)%p(1:nh,nh+1,:)
    lcross_edge_xyz = lcross_edge_xyz .and. cross_edge_xyz_check(zq,zp,nh)
    zq = f1(ind)%p(0,1:nh,:); zp = f2(ind)%p(0,1:nh,:)
    lcross_edge_xyz = lcross_edge_xyz .and. cross_edge_xyz_check(zq,zp,nh)
    zq = f1(ind)%p(nh+1,1:nh,:); zp = f2(ind)%p(nh+1,1:nh,:)
    lcross_edge_xyz = lcross_edge_xyz .and. cross_edge_xyz_check(zq,zp,nh)
end do

if(lcross_edge_xyz) then
    print *, "cross edge xyz test passed"
else
    print *, "cross edge xyz test failed"
end if

!check the following facts
!1)cubsph basis vectors tangent to sphere
!2)acov orthogonal to bctv and vice versa
!3)metric tensor is matrix of covariant vectors dot-prods
!4)actv = qi11*acov+qi12*bcov and same for bctv
!5)G = det(Q) is determinant of Q
err_ort_abt = 0._8; err_ort_rab = 0._8
err_q = 0._8; err_qi = 0._8; err_g = 0._8
do ind=ts,te
    js = mesh(ind)%js-mesh(ind)%halo_width
    je = mesh(ind)%je+mesh(ind)%halo_width
    is = mesh(ind)%is-mesh(ind)%halo_width
    ie = mesh(ind)%ie+mesh(ind)%halo_width
    do j=mesh(ind)%js-mesh(ind)%halo_width,mesh(ind)%je+mesh(ind)%halo_width
        do i=mesh(ind)%is-mesh(ind)%halo_width,mesh(ind)%ie+mesh(ind)%halo_width
            err_ort_rab = err_ort_rab+&
                              abs(mesh(ind)%acov(1,i,j)*mesh(ind)%rhx(i,j)+&
                                  mesh(ind)%acov(2,i,j)*mesh(ind)%rhy(i,j)+&
                                  mesh(ind)%acov(3,i,j)*mesh(ind)%rhz(i,j))
            err_ort_rab = err_ort_rab+&
                              abs(mesh(ind)%bcov(1,i,j)*mesh(ind)%rhx(i,j)+&
                                  mesh(ind)%bcov(2,i,j)*mesh(ind)%rhy(i,j)+&
                                  mesh(ind)%bcov(3,i,j)*mesh(ind)%rhz(i,j))
            err_ort_rab = err_ort_rab+&
                              abs(mesh(ind)%actv(1,i,j)*mesh(ind)%rhx(i,j)+&
                                  mesh(ind)%actv(2,i,j)*mesh(ind)%rhy(i,j)+&
                                  mesh(ind)%actv(3,i,j)*mesh(ind)%rhz(i,j))
            err_ort_rab = err_ort_rab+&
                              abs(mesh(ind)%bctv(1,i,j)*mesh(ind)%rhx(i,j)+&
                                  mesh(ind)%bctv(2,i,j)*mesh(ind)%rhy(i,j)+&
                                  mesh(ind)%bctv(3,i,j)*mesh(ind)%rhz(i,j))
            err_ort_abt = err_ort_abt+&
                              abs(sum(mesh(ind)%acov(:,i,j)*mesh(ind)%bctv(:,i,j)))
            err_ort_abt = err_ort_abt+&
                              abs(sum(mesh(ind)%actv(:,i,j)*mesh(ind)%bcov(:,i,j)))
            err_q = err_q+abs(sum(mesh(ind)%acov(:,i,j)**2)-mesh(ind)%Q(1,i,j))+&
                          abs(sum(mesh(ind)%acov(:,i,j)*mesh(ind)%bcov(:,i,j))-mesh(ind)%Q(2,i,j))+&
                          abs(sum(mesh(ind)%bcov(:,i,j)**2)-mesh(ind)%Q(3,i,j))
            err_qi = err_qi+sum(abs(mesh(ind)%actv(:,i,j)                    -&
                                    mesh(ind)%Qi(1,i,j)*mesh(ind)%acov(:,i,j)-&
                                    mesh(ind)%Qi(2,i,j)*mesh(ind)%bcov(:,i,j)))
            err_qi = err_qi+sum(abs(mesh(ind)%bctv(:,i,j)                    -&
                                    mesh(ind)%Qi(2,i,j)*mesh(ind)%acov(:,i,j)-&
                                    mesh(ind)%Qi(3,i,j)*mesh(ind)%bcov(:,i,j)))
            err_g = err_g + abs(sqrt(mesh(ind)%Q(1,i,j)*mesh(ind)%Q(3,i,j)-mesh(ind)%Q(2,i,j)**2)-&
                                     mesh(ind)%G(i,j))
        end do
    end do
    err_ort_rab = err_ort_rab / (4*(ie-is+1)*(je-js+1))
    err_ort_abt = err_ort_abt / (2*(ie-is+1)*(je-js+1))
    err_q = err_q / ((ie-is+1)*(je-js+1))
    err_qi = err_qi / (2*(ie-is+1)*(je-js+1))
    err_g = err_g / ((ie-is+1)*(je-js+1))
end do
    err_ort_rab = err_ort_rab / (te-ts+1)
    err_ort_abt = err_ort_abt / (te-ts+1)
    err_q = err_q / (te-ts+1)
    err_qi = err_qi / (te-ts+1)
    err_g = err_g / (te-ts+1)
!print *, "metric errors", err_ort_rab, err_ort_abt, err_q, err_qi, err_g
if(max(err_ort_rab, err_ort_abt, err_q, err_qi, err_g) < zeps) then
    print *, "ecs metric test passed"
end if

contains
logical function cross_edge_xyz_check(q,p,nx) result(lpass)
real(kind=8) p(1:nx,3), q(1:nx,3)
integer nx
real(kind=8) znorm(3), zpr
real(kind=8) za(3), zb(3)

za = [q(nx,1)-q(1,1),q(nx,2)-q(1,2),q(nx,3)-q(1,3)]
zb = [p(nx,1)-p(1,1),p(nx,2)-p(1,2),p(nx,3)-p(1,3)]
lpass = sum(za**2) > sum(zb**2)

za = za / sqrt(sum(za**2))
zb = zb / sqrt(sum(zb**2))

lpass = lpass .and. sum(abs(za-zb))<zeps

znorm = [q(1,2)*q(nx,3)-q(1,3)*q(nx,2), q(1,3)*q(nx,1)-q(1,1)*q(nx,3), q(1,1)*q(nx,2)-q(1,2)*q(nx,1)]

zpr = abs(sum(p(1,:)*znorm))+abs(sum(p(nx,:)*znorm))

lpass = lpass .and. zpr<zeps

!middle connection check
if(mod(nx,2) == 0) then
    za = [q(nx/2,1)+q(nx/2+1,1),q(nx/2,2)+q(nx/2+1,2),q(nx/2,3)+q(nx/2+1,3)]
    zb = [p(nx/2,1)+p(nx/2+1,1),p(nx/2,2)+p(nx/2+1,2),p(nx/2,3)+p(nx/2+1,3)]
else
    za = [q(nx/2+1,1),q(nx/2+1,2),q(nx/2+1,3)]
    zb = [p(nx/2+1,1),p(nx/2+1,2),p(nx/2+1,3)]
end if

za = za / sqrt(sum(za**2))
zb = zb / sqrt(sum(zb**2))

zpr = acos(sum(za*zb))

lpass = lpass .and. abs(zpr) < zeps

end function cross_edge_xyz_check

end subroutine test_metric

!check symmetricity of pressure(h)-points on ecs grid face
logical function symmetricity_check_h(phx,phy,phz,nn,panel_ind) result(lsym_check)
use topology_mod, only : n, ex, ey
integer(kind=4), intent(in) :: nn, panel_ind !number of h-points incl halo, ind of cs-panel
real(kind=8),    intent(in) :: phx(nn, nn), phy(nn, nn), phz(nn, nn)
!local
real(kind=8) ee(3)

lsym_check = .true.

ee = real(ex(:,panel_ind),8)
lsym_check = lsym_check .and. symmetric( 1,nn, 1, 1,panel_ind,ee)
lsym_check = lsym_check .and. symmetric( 1,nn,nn,nn,panel_ind,ee)
ee = real(ey(:,panel_ind),8)
lsym_check = lsym_check .and. symmetric( 1, 1, 1,nn,panel_ind,ee)
lsym_check = lsym_check .and. symmetric(nn,nn, 1,nn,panel_ind,ee)

contains
logical function symmetric(i1,i2,j1,j2,panel_ind,ee) result(lsym)
    integer(kind=4), intent(in) :: i1,i2,j1,j2,panel_ind
    real(kind=8),    intent(in) :: ee(3)
    !local
    real(kind=8) a(3), b(3), zpr
    real(kind=8), parameter :: zeps = 1e-16_8

    a = [phx(i1,j1), phy(i1,j1), phz(i1,j1)]
    b = [phx(i2,j2), phy(i2,j2), phz(i2,j2)]
    zpr = sum(a*ee)
    a = a - 2._8*zpr*ee !flip a-vector
    lsym = sum(abs(a-b))<zeps
end function symmetric

end function symmetricity_check_h

end module test_metric_mod
