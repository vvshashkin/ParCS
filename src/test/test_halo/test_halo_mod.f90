module test_halo_mod
implicit none

contains

subroutine test_ecs_halo()

use mpi

use grid_function_mod,    only : grid_function_t
use exchange_mod,         only : exchange_t
use partition_mod,        only : partition_t
use exchange_factory_mod, only : create_2d_full_halo_exchange, create_2d_cross_halo_exchange
use mesh_factory_mod,     only : create_equiangular_mesh
use mesh_mod,             only : mesh_t

use ecs_halo_mod,         only : ecs_halo_mod_init, whalo

type(exchange_t)                   :: exch_halo
type(partition_t)                  :: partition
type(mesh_t),          allocatable :: mesh(:)
type(grid_function_t), allocatable :: f1(:)
type(grid_function_t), allocatable :: f2(:)

integer(kind=4), parameter         :: nh=128, nz=3, halo_width=3, ex_halo_width=8
integer(kind=4), parameter         :: nn(3) = [nh, nh/2, nh/4]
integer(kind=4)                    :: myid, np, ierr, code

integer(kind=4) :: ts, te
integer(kind=4) :: ind, i, j, k, ifc, err_sum, gl_err_sum
integer(kind=4) :: iev, isv, jev, jsv
integer(kind=4) :: ie, is, js, je, klev

real(kind=8) err, inface_err, cross_edge_err
real(kind=8) err_max, inface_err_max, cross_edge_err_max
real(kind=8) inface_corner_err, inface_corner_err_max
real(kind=8) inedge_corner_err, inedge_corner_err_max
real(kind=8) halo_corner_err, halo_corner_err_max
real(kind=8) gl_inface_err, gl_cross_edge_err
real(kind=8) gl_inface_err_max, gl_cross_edge_err_max
real(kind=8) gl_inface_corner_err, gl_inface_corner_err_max
real(kind=8) gl_inedge_corner_err, gl_inedge_corner_err_max
real(kind=8) gl_halo_corner_err, gl_halo_corner_err_max
integer(kind=4) num_inedge_corners, gl_num_inedge_corners

logical lpass, glpass

call MPI_comm_rank(mpi_comm_world , myid, ierr)
call MPI_comm_size(mpi_comm_world , Np  , ierr)

call partition%init(nh, nz, max(1,Np/6), Np, strategy = 'default')
!call partition%init(nh, nz, 64, Np, strategy = 'default')

!find start and end index of tiles belonging to the current proccesor
ts = findloc(partition%proc_map, myid, dim=1)
te = findloc(partition%proc_map, myid, back = .true., dim=1)

allocate(mesh(ts:te))
do ind = ts, te
    call create_equiangular_mesh(mesh(ind), partition%tile(ind)%is, partition%tile(ind)%ie, &
                                            partition%tile(ind)%js, partition%tile(ind)%je, &
                                            partition%tile(ind)%ks, partition%tile(ind)%ke, &
                                            nh, ex_halo_width, partition%tile(ind)%panel_number)
end do

call ecs_halo_mod_init(nn,halo_width)
call mpi_barrier(mpi_comm_world, ierr)
if (myid==0) print *, 'equiangular cubed-sphere halo-zone interpolation test'

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
     isv = f1(ind)%is-f1(ind)%nvi
     iev = f1(ind)%ie+f1(ind)%nvi
     jsv = f1(ind)%js-f1(ind)%nvj
     jev = f1(ind)%je+f1(ind)%nvj
     f1(ind).p(isv:iev,jsv:jev,1) = mesh(ind)%rhx(isv:iev,jsv:jev)
     f1(ind).p(isv:iev,jsv:jev,2) = mesh(ind)%rhy(isv:iev,jsv:jev)
     f1(ind).p(isv:iev,jsv:jev,3) = mesh(ind)%rhz(isv:iev,jsv:jev)
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

inface_err = 0._8; inface_err_max = 0._8
cross_edge_err = 0._8; cross_edge_err_max = 0._8
inface_corner_err = 0._8; inface_corner_err_max = 0._8
inedge_corner_err = 0._8; inedge_corner_err_max = 0._8
halo_corner_err = 0._8; halo_corner_err_max = 0._8
num_inedge_corners = 0
do ind = ts, te

    is = partition%tile(ind)%is; ie = partition%tile(ind)%ie
    js = partition%tile(ind)%js; je = partition%tile(ind)%je
    klev = partition%tile(ind)%ke-partition%tile(ind)%ks+1

    !tile edge errors
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
    if(je == nh) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else 
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     = sum(abs(f1(ind)%p(is-halo_width:is-1,js-halo_width:js-1,:)-f2(ind)%p(is-halo_width:is-1,js-halo_width:js-1,:))) / halo_width**2
    err_max = maxval(abs(f1(ind)%p(is-halo_width:is-1,js-halo_width:js-1,:)-f2(ind)%p(is-halo_width:is-1,js-halo_width:js-1,:)))
    if(is == 1 .and. js == 1) then
        halo_corner_err = halo_corner_err+err
        halo_corner_err_max = max(halo_corner_err_max,err_max)
    else if(is == 1 .or. js == 1) then
        inedge_corner_err = inedge_corner_err+err
        inedge_corner_err_max = max(inedge_corner_err_max,err_max)
        num_inedge_corners = num_inedge_corners+1
    else
        inface_corner_err = inface_corner_err+err
        inface_corner_err_max = max(inface_corner_err_max,err_max)
    end if

    err     = sum(abs(f1(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,:)-f2(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,:))) / halo_width**2
    err_max = maxval(abs(f1(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,:)-f2(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,:)))
    if(ie == nh .and. js == 1) then
        halo_corner_err = halo_corner_err+err
        halo_corner_err_max = max(halo_corner_err_max,err_max)
    else if(ie == nh .or. js == 1) then
        inedge_corner_err = inedge_corner_err+err
        inedge_corner_err_max = max(inedge_corner_err_max,err_max)
        num_inedge_corners = num_inedge_corners+1
    else
        inface_corner_err = inface_corner_err+err
        inface_corner_err_max = max(inface_corner_err_max,err_max)
    end if

    err     = sum(abs(f1(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,:)-f2(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,:))) / halo_width**2
    err_max = maxval(abs(f1(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,:)-f2(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,:)))
    if(ie == nh .and. je == nh) then
        halo_corner_err = halo_corner_err+err
        halo_corner_err_max = max(halo_corner_err_max,err_max)
    else if(ie == nh .or. je == nh) then
        inedge_corner_err = inedge_corner_err+err
        inedge_corner_err_max = max(inedge_corner_err_max,err_max)
        num_inedge_corners = num_inedge_corners+1
    else
        inface_corner_err = inface_corner_err+err
        inface_corner_err_max = max(inface_corner_err_max,err_max)
    end if

    err     = sum(abs(f1(ind)%p(is-halo_width:is-1,je+1:je+halo_width,:)-f2(ind)%p(is-halo_width:is-1,je+1:je+halo_width,:))) / halo_width**2
    err_max = maxval(abs(f1(ind)%p(is-halo_width:is-1,je+1:je+halo_width,:)-f2(ind)%p(is-halo_width:is-1,je+1:je+halo_width,:)))
    if(is == 1 .and. je == nh) then
        halo_corner_err = halo_corner_err+err
        halo_corner_err_max = max(halo_corner_err_max,err_max)
    else if(is == 1 .or. je == nh) then
        inedge_corner_err = inedge_corner_err+err
        inedge_corner_err_max = max(inedge_corner_err_max,err_max)
        num_inedge_corners = num_inedge_corners+1
    else
        inface_corner_err = inface_corner_err+err
        inface_corner_err_max = max(inface_corner_err_max,err_max)
    end if
end do

call mpi_allreduce(cross_edge_err, gl_cross_edge_err, 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(inface_err,     gl_inface_err    , 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(cross_edge_err_max, gl_cross_edge_err_max, 1, mpi_double, mpi_max, mpi_comm_world, ierr)
call mpi_allreduce(inface_err_max,     gl_inface_err_max    , 1, mpi_double, mpi_max, mpi_comm_world, ierr)

gl_cross_edge_err = gl_cross_edge_err/(6*4*klev*halo_width) !scale to 1 element error

call mpi_allreduce(inface_corner_err,     gl_inface_corner_err    , 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(inface_corner_err_max,     gl_inface_corner_err_max    , 1, mpi_double, mpi_max, mpi_comm_world, ierr)
call mpi_allreduce(inedge_corner_err,     gl_inedge_corner_err    , 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(inedge_corner_err_max,     gl_inedge_corner_err_max    , 1, mpi_double, mpi_max, mpi_comm_world, ierr)
call mpi_allreduce(num_inedge_corners,     gl_num_inedge_corners    , 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(halo_corner_err,     gl_halo_corner_err    , 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(halo_corner_err_max,     gl_halo_corner_err_max    , 1, mpi_double, mpi_max, mpi_comm_world, ierr)

gl_inedge_corner_err = gl_inedge_corner_err / max(gl_num_inedge_corners,1) / klev
gl_halo_corner_err = gl_halo_corner_err / 48 / klev

if(myid == 0) then
  print *, "halo errors:"
  print *, "inface mean abs, inface max, cross edge mean abs, cross edge max"
  print '(4e15.7)', gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max
  print *, "inface-corner mean abs, inface-corner max, inedge-corner mean abs, inedge-corner max, halo-corner mean abs, halo-corner max"
  print '(6e15.7)', gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err, gl_inedge_corner_err_max,&
                    gl_halo_corner_err, gl_halo_corner_err_max

end if


end subroutine test_ecs_halo

end module test_halo_mod
