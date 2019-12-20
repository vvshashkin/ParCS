module test_metric_mod
implicit none

contains

subroutine test_metric()

use mpi

use grid_function_mod,    only : grid_function_t
use exchange_mod,         only : exchange_t
use partition_mod,        only : partition_t
use exchange_factory_mod, only : create_2d_full_halo_exchange, create_2d_cross_halo_exchange
use mesh_factory_mod,     only : create_equiangular_mesh
use mesh_mod,             only : mesh_t

type(exchange_t)                   :: exch_halo
type(partition_t)                  :: partition
type(grid_function_t), allocatable :: f1(:)
type(grid_function_t), allocatable :: f2(:)
type(mesh_t),          allocatable :: mesh(:)

integer(kind=4), parameter         :: nh=16, nz=3, halo_width=1
integer(kind=4)                    :: myid, np, ierr, code

integer(kind=4) :: ts, te
integer(kind=4) :: ind, i, j, k, ifc, err_sum, gl_err_sum
integer(kind=4) :: ie, is

integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number

character(:), allocatable :: frmt
logical lcross_edge_xyz
real(kind=8) zq(1:nh,3), zp(1:nh,3)


call MPI_comm_rank(mpi_comm_world , myid, ierr)
call MPI_comm_size(mpi_comm_world , Np  , ierr)

if (myid==0) print *, 'equiangular cubed-sphere test'

if(Np >1) then
    print *, "this test currently runs only with 1 mpi-process. exiting"
    return
end if

!if(myid == 0) then
!    call ecs_geometry_mod_check()
!end if

call mpi_barrier(mpi_comm_world, ierr)

call partition%init(nh, nz, max(1,Np/6), Np, strategy = 'default')

!find start and end index of tiles belonging to the current proccesor
ts = findloc(partition%proc_map, myid, dim=1)
te = findloc(partition%proc_map, myid, back = .true., dim=1)

allocate(mesh(ts:te))
do ind = ts, te
    call create_equiangular_mesh(mesh(ind), partition%tile(ind)%is, partition%tile(ind)%ie, &
                                            partition%tile(ind)%js, partition%tile(ind)%je, &
                                            partition%tile(ind)%ks, partition%tile(ind)%ke, &
                                            nh, halo_width, partition%tile(ind)%panel_number)
!    mesh(ind)%halo = init_ecs_halo(mesh(ind)%is, mesh(ind)%ie, &
!                                   mesh(ind)%js, mesh(ind)%je, &
!                                   mesh(ind)%nx, halo_width,   &
!                                   mesh(ind)%hx)
end do

!Init arrays

allocate(f1(ts:te))
allocate(f2(ts:te))

do i = ts, te

    call f1(i)%init(partition%tile(i)%is, partition%tile(i)%ie, &
                    partition%tile(i)%js, partition%tile(i)%je, &
                    partition%tile(i)%ks, partition%tile(i)%ke, &
                    halo_width, halo_width, 0)
    call f2(i)%init(partition%tile(i)%is, partition%tile(i)%ie, &
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
call create_2d_full_halo_exchange(exch_halo, partition, halo_width, myid, np)

!Perform exchange
call exch_halo%do(f1, ts, te)
call mpi_barrier(mpi_comm_world, ierr)

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

contains
logical function cross_edge_xyz_check(q,p,nx) result(lpass)
real(kind=8) p(1:nx,3), q(1:nx,3)
integer nx
real(kind=8) znorm(3), zpr
real(kind=8) za(3), zb(3)
real(kind=8), parameter :: zeps = 1e-15_8

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

end module test_metric_mod
