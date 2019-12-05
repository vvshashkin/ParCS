module test_halo_mod
implicit none

contains

subroutine test_ecs_halo()

use mpi

use grid_function_mod,    only : grid_function_t
use exchange_mod,         only : exchange_t
use partition_mod,        only : partition_t
use exchange_factory_mod, only : create_2d_full_halo_exchange, create_2d_cross_halo_exchange

use ecs_geometry_mod,     only : ecs_geometry_mod_init, ecs_geometry_mod_check, rhx, rhy, rhz
use ecs_halo_mod,         only : ecs_halo_mod_init, whalo

type(exchange_t)                   :: exch_halo
type(partition_t)                  :: partition
type(grid_function_t), allocatable :: f1(:)
type(grid_function_t), allocatable :: f2(:)

integer(kind=4), parameter         :: nh=128, nz=3, halo_width=3, ex_halo_width=5
integer(kind=4), parameter         :: nn(3) = [nh, nh/2, nh/4]
integer(kind=4)                    :: myid, np, ierr, code

integer(kind=4) :: ts, te
integer(kind=4) :: ind, i, j, k, ifc, err_sum, gl_err_sum
integer(kind=4) :: iev, isv, jev, jsv
integer(kind=4) :: ie, is, js, je, klev

real(kind=8) err, inface_err, cross_edge_err
real(kind=8) err_max, inface_err_max, cross_edge_err_max
real(kind=8) gl_inface_err, gl_cross_edge_err
real(kind=8) gl_inface_err_max, gl_cross_edge_err_max

integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number

character(:), allocatable :: frmt
logical lpass, glpass

call MPI_comm_rank(mpi_comm_world , myid, ierr)
call MPI_comm_size(mpi_comm_world , Np  , ierr)

call ecs_geometry_mod_init(nh, ex_halo_width)

call mpi_barrier(mpi_comm_world, ierr)
if (myid==0) print *, 'equiangular cubed-sphere halo-zone interpolation test'

call ecs_halo_mod_init(nn,halo_width)

lpass = ( size(whalo) == size(nn) )
do k = 1, size(nn)
    lpass = lpass .and. (whalo(k).n == nn(k))
end do

glpass = .false.

call mpi_allreduce(lpass, glpass, 1, mpi_logical, mpi_land, mpi_comm_world, ierr)

if(glpass) then
   if(myid == 0) print *, "ecs halo module initialization consistent"
else
   if(myid == 0) print *, "ecs halo module initialization inconsistent"
   return
end if

call partition%init(nh, nz, max(1,Np/6), Np, strategy = 'default')
!call partition%init(nh, nz, 4, Np, strategy = 'default')

!find start and end index of tiles belonging to the current proccesor
ts = findloc(partition%proc_map, myid, dim=1)
te = findloc(partition%proc_map, myid, back = .true., dim=1)

!Init arrays

allocate(f1(ts:te))
allocate(f2(ts:te))

do i = ts, te

    call f1(i)%init(partition%tile(i)%is, partition%tile(i)%ie, &
                    partition%tile(i)%js, partition%tile(i)%je, &
                    partition%tile(i)%ks, partition%tile(i)%ke, &
                    ex_halo_width, ex_halo_width, 0)
    call f2(i)%init(partition%tile(i)%is, partition%tile(i)%ie, &
                    partition%tile(i)%js, partition%tile(i)%je, &
                    partition%tile(i)%ks, partition%tile(i)%ke, &
                    ex_halo_width, ex_halo_width, 0)
end do

do ind = ts, te
     ifc = partition%tile(ind)%panel_number
     isv = f1(ind)%is-f1(ind)%nvi
     iev = f1(ind)%ie+f1(ind)%nvi
     jsv = f1(ind)%js-f1(ind)%nvj
     jev = f1(ind)%je+f1(ind)%nvj
     f1(ind).p(isv:iev,jsv:jev,1) = rhx(isv:iev,jsv:jev,ifc)
     f1(ind).p(isv:iev,jsv:jev,2) = rhy(isv:iev,jsv:jev,ifc)
     f1(ind).p(isv:iev,jsv:jev,3) = rhz(isv:iev,jsv:jev,ifc)
     f2(ind).p(:,:,:) = f1(ind).p(:,:,:)
end do



!Init exchange
call create_2d_full_halo_exchange(exch_halo, partition, ex_halo_width, myid, np)

!Perform exchange
call exch_halo%do(f1, ts, te)
call mpi_barrier(mpi_comm_world, ierr)
do ind = ts, te
    call whalo(1)%ext_halo(f1(ind))
end do

!if(myid == 0) then
!    do i = 1, nh
!        !print *, i, f1(ts)%p(0,i,1), f2(ts)%p(0,i,1)
!        print '(i4,3e15.7)', i, f1(ts)%p(i,0,1)-f2(ts)%p(i,0,1), f1(ts)%p(i,0,2)-f2(ts)%p(i,0,2), f1(ts)%p(i,0,3)-f2(ts)%p(i,0,3)
!        !print '(i4,3e15.7)', i, f1(ts)%p(i,0,1)-f2(ts)%p(i,0,1), f1(ts)%p(i,-1,1)-f2(ts)%p(i,-1,1), f1(ts)%p(i,-2,1)-f2(ts)%p(i,-2,1)
!    end do
!end if

inface_err = 0._8; inface_err_max = 0._8
cross_edge_err = 0._8; cross_edge_err_max = 0._8
do ind = ts, te

    is = partition%tile(ind)%is; ie = partition%tile(ind)%ie
    js = partition%tile(ind)%js; je = partition%tile(ind)%je
    klev = partition%tile(ind)%ke-partition%tile(ind)%ks+1

    err     = sum(abs(f1(ind)%p(is-halo_width:is-1,js:je,:)-f2(ind)%p(is-halo_width:is-1,js:je,:)))/nh
    err_max = maxval(abs(f1(ind)%p(is-halo_width:is-1,js:je,:)-f2(ind)%p(is-halo_width:is-1,js:je,:)))
    if(is == 1) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else 
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     = sum(abs(f1(ind)%p(ie+1:ie+halo_width,js:je,:)-f2(ind)%p(ie+1:ie+halo_width,js:je,:)))/nh
    err_max = maxval(abs(f1(ind)%p(ie+1:ie+halo_width,js:je,:)-f2(ind)%p(ie+1:ie+halo_width,js:je,:)))
    if(ie == nh) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else 
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     = sum(abs(f1(ind)%p(is:ie,js-halo_width:js-1,:)-f2(ind)%p(is:ie,js-halo_width:js-1,:)))/nh
    err_max = maxval(abs(f1(ind)%p(is:ie,js-halo_width:js-1,:)-f2(ind)%p(is:ie,js-halo_width:js-1,:)))
    if(js == 1) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else 
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     = sum(abs(f1(ind)%p(is:ie,je+1:je+halo_width,:)-f2(ind)%p(is:ie,je+1:je+halo_width,:)))/nh
    err_max = maxval(abs(f1(ind)%p(is:ie,je+1:je+halo_width,:)-f2(ind)%p(is:ie,je+1:je+halo_width,:)))
    if(js == 1) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else 
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if
end do

call mpi_allreduce(cross_edge_err, gl_cross_edge_err, 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(inface_err,     gl_inface_err    , 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(cross_edge_err_max, gl_cross_edge_err_max, 1, mpi_double, mpi_max, mpi_comm_world, ierr)
call mpi_allreduce(inface_err_max,     gl_inface_err_max    , 1, mpi_double, mpi_max, mpi_comm_world, ierr)

gl_cross_edge_err = gl_cross_edge_err/(6*4*klev*halo_width) !scale to 1 element error

if(myid == 0) then
  print *, "halo errors:"
  print *, "inface mean abs, inface max, cross edge mean abs, cross edge max"
  print '(4e15.7)', gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max
end if


end subroutine test_ecs_halo

end module test_halo_mod
