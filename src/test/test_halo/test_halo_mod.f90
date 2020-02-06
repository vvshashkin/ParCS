module test_halo_mod
implicit none

private
public   :: test_ecs_halo

contains

subroutine test_ecs_halo()

use mpi

use grid_function_mod,     only : grid_function_t
use exchange_abstract_mod, only : exchange_t
use partition_mod,         only : partition_t
use exchange_factory_mod,  only : create_2d_halo_exchange
use mesh_factory_mod,      only : create_equiangular_mesh
use mesh_mod,              only : mesh_t
use ecs_halo_mod,          only : ecs_halo_t
use ecs_halo_factory_mod,  only : init_ecs_halo
use ecs_halo_vec_a_factory_mod, only : init_ecs_halo_vect


class(exchange_t),     allocatable :: exch_halo
type(partition_t)                  :: partition
type(mesh_t),          allocatable :: mesh(:)
type(grid_function_t), allocatable :: f_test(:), u_test(:), v_test(:)
type(grid_function_t), allocatable :: f_true(:), u_true(:), v_true(:)
type(ecs_halo_t),      allocatable :: halo

integer(kind=4), parameter         :: nh=128, nz=3, halo_width=3, ex_halo_width=8
integer(kind=4)                    :: myid, np, ierr, code

integer(kind=4) :: ts, te
integer(kind=4) :: ind, i, j, k, ifc, err_sum, gl_err_sum
integer(kind=4) :: iev, isv, jev, jsv
integer(kind=4) :: ie, is, js, je, klev, ke, ks

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
    halo = init_ecs_halo(mesh(ind)%is, mesh(ind)%ie, &
                         mesh(ind)%js, mesh(ind)%je, &
                         mesh(ind)%nx, halo_width,   &
                         mesh(ind)%hx)
    mesh(ind)%halo = halo
    mesh(ind)%halo_vec = init_ecs_halo_vect(mesh(ind)%panel_ind,&
                                   mesh(ind)%is, mesh(ind)%ie,  &
                                   mesh(ind)%js, mesh(ind)%je,  &
                                   mesh(ind)%nx, halo_width,    &
                                   mesh(ind)%hx, halo)
end do

call mpi_barrier(mpi_comm_world, ierr)
if (myid==0) print *, 'equiangular cubed-sphere halo-zone interpolation test'

call init_scalar_halo_test_fun(f_test,ts,te,partition,mesh,ex_halo_width)
call init_scalar_halo_test_fun(f_true,ts,te,partition,mesh,ex_halo_width)
call init_vector_halo_test_fun(u_test,v_test,ts,te,partition,mesh,ex_halo_width)
call init_vector_halo_test_fun(u_true,v_true,ts,te,partition,mesh,ex_halo_width)

!Init exchange
exch_halo = create_2d_halo_exchange(partition, ex_halo_width, 'full', myid, np)

!Perform exchange
call exch_halo%do(f_test, lbound(f_test, 1), ubound(f_test, 1))
call exch_halo%do(u_test, lbound(u_test, 1), ubound(u_test, 1))
call exch_halo%do(v_test, lbound(v_test, 1), ubound(v_test, 1))

do ind = ts, te
    call mesh(ind)%halo%interp(f_test(ind),halo_width)
    call mesh(ind)%halo_vec%interpv(u_test(ind),v_test(ind),halo_width)
end do

call halo_err(gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max, &
              gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err,       &
              gl_inedge_corner_err_max, gl_halo_corner_err, gl_halo_corner_err_max,       &
              nh, halo_width, f_test, f_true, te-ts+1)

if(myid == 0) then
  print *, "halo errors: F"
  print *, "inface mean abs, inface max, cross edge mean abs, cross edge max"
  print '(4e15.7)', gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max
  print *, "inface-corner mean abs, inface-corner max, inedge-corner mean abs, inedge-corner max, halo-corner mean abs, halo-corner max"
  print '(6e15.7)', gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err, gl_inedge_corner_err_max,&
                    gl_halo_corner_err, gl_halo_corner_err_max
end if
call halo_err(gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max, &
              gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err,       &
              gl_inedge_corner_err_max, gl_halo_corner_err, gl_halo_corner_err_max,       &
              nh, halo_width, u_test, u_true, te-ts+1)

if(myid == 0) then
  print *, "halo errors: U"
  print *, "inface mean abs, inface max, cross edge mean abs, cross edge max"
  print '(4e15.7)', gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max
  print *, "inface-corner mean abs, inface-corner max, inedge-corner mean abs, inedge-corner max, halo-corner mean abs, halo-corner max"
  print '(6e15.7)', gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err, gl_inedge_corner_err_max,&
                    gl_halo_corner_err, gl_halo_corner_err_max
end if

call halo_err(gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max, &
              gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err,       &
              gl_inedge_corner_err_max, gl_halo_corner_err, gl_halo_corner_err_max,       &
              nh, halo_width, v_test, v_true, te-ts+1)

if(myid == 0) then
  print *, "halo errors: V"
  print *, "inface mean abs, inface max, cross edge mean abs, cross edge max"
  print '(4e15.7)', gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max
  print *, "inface-corner mean abs, inface-corner max, inedge-corner mean abs, inedge-corner max, halo-corner mean abs, halo-corner max"
  print '(6e15.7)', gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err, gl_inedge_corner_err_max,&
                    gl_halo_corner_err, gl_halo_corner_err_max
end if

end subroutine test_ecs_halo



subroutine init_scalar_halo_test_fun(f,ts,te,partition,mesh,halo_width)

use grid_function_mod, only: grid_function_t
use partition_mod,     only: partition_t
use mesh_mod,          only: mesh_t

type(grid_function_t), allocatable, intent(out) :: f(:)

integer(kind=4),   intent(in) :: ts, te
type(partition_t), intent(in) :: partition
type(mesh_t),      intent(in) :: mesh(ts:te)
integer(kind=4),   intent(in) :: halo_width

!locals
integer(kind=4) ind, isv, iev, jsv, jev

!Init arrays
allocate(f(ts:te))

do ind = ts, te
     call f(ind)%init(partition%tile(ind)%panel_number,               &
                      partition%tile(ind)%is, partition%tile(ind)%ie, &
                      partition%tile(ind)%js, partition%tile(ind)%je, &
                      partition%tile(ind)%ks, partition%tile(ind)%ke, &
                      halo_width, halo_width, 0)
     isv = mesh(ind)%is-halo_width
     iev = mesh(ind)%ie+halo_width
     jsv = mesh(ind)%js-halo_width
     jev = mesh(ind)%je+halo_width
     f(ind).p(isv:iev,jsv:jev,1) = mesh(ind)%rhx(isv:iev,jsv:jev)
     f(ind).p(isv:iev,jsv:jev,2) = mesh(ind)%rhy(isv:iev,jsv:jev)
     f(ind).p(isv:iev,jsv:jev,3) = mesh(ind)%rhz(isv:iev,jsv:jev)
end do

end subroutine init_scalar_halo_test_fun



subroutine init_vector_halo_test_fun(u,v,ts,te,partition,mesh,halo_width)

use grid_function_mod, only: grid_function_t
use partition_mod,     only: partition_t
use mesh_mod,          only: mesh_t

type(grid_function_t), allocatable, intent(out) :: u(:), v(:)

integer(kind=4),   intent(in) :: ts, te
type(partition_t), intent(in) :: partition
type(mesh_t),      intent(in) :: mesh(ts:te)
integer(kind=4),   intent(in) :: halo_width

!locals
integer(kind=4) ind, isv, iev, jsv, jev, i, j
real(kind=8) vx, vy, vz
!Init arrays
allocate(u(ts:te), v(ts:te))

do ind = ts, te
     call u(ind)%init(partition%tile(ind)%panel_number,               &
                      partition%tile(ind)%is, partition%tile(ind)%ie, &
                      partition%tile(ind)%js, partition%tile(ind)%je, &
                      partition%tile(ind)%ks, partition%tile(ind)%ke, &
                      halo_width, halo_width, 0)
     call v(ind)%init(partition%tile(ind)%panel_number,               &
                      partition%tile(ind)%is, partition%tile(ind)%ie, &
                      partition%tile(ind)%js, partition%tile(ind)%je, &
                      partition%tile(ind)%ks, partition%tile(ind)%ke, &
                      halo_width, halo_width, 0)

     isv = mesh(ind)%is-halo_width
     iev = mesh(ind)%ie+halo_width
     jsv = mesh(ind)%js-halo_width
     jev = mesh(ind)%je+halo_width
     do j=jsv,jev
         do i=isv,iev
             !level 1: solid rotation around (1,0,0) axis
             vx = 0._8; vy = mesh(ind)%rhz(i,j); vz =-mesh(ind)%rhy(i,j)
             u(ind)%p(i,j,1) = sum([vx,vy,vz]*mesh(ind)%actv(:,i,j))
             v(ind)%p(i,j,1) = sum([vx,vy,vz]*mesh(ind)%bctv(:,i,j))
              !level 2: solid rotation around (0,1,0) axis
              vx =-mesh(ind)%rhz(i,j); vy = 0._8; vz = mesh(ind)%rhx(i,j)
              u(ind)%p(i,j,2) = sum([vx,vy,vz]*mesh(ind)%actv(:,i,j))
              v(ind)%p(i,j,2) = sum([vx,vy,vz]*mesh(ind)%bctv(:,i,j))
             !level 3: solid rotation around (0,0,1) axis
             vx = mesh(ind)%rhy(i,j); vy =-mesh(ind)%rhx(i,j); vz = 0._8
             u(ind)%p(i,j,3) = sum([vx,vy,vz]*mesh(ind)%actv(:,i,j))
             v(ind)%p(i,j,3) = sum([vx,vy,vz]*mesh(ind)%bctv(:,i,j))
         end do
     end do
end do

end subroutine init_vector_halo_test_fun


subroutine halo_err(&
              gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max, &
              gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err,       &
              gl_inedge_corner_err_max, gl_halo_corner_err, gl_halo_corner_err_max,       &
              nh, halo_width, f1, f2, ntiles)

use mpi
use grid_function_mod, only: grid_function_t

real(kind=8),          intent(out) :: gl_inface_err, gl_cross_edge_err
real(kind=8),          intent(out) :: gl_inface_err_max, gl_cross_edge_err_max
real(kind=8),          intent(out) :: gl_inface_corner_err, gl_inface_corner_err_max
real(kind=8),          intent(out) :: gl_inedge_corner_err, gl_inedge_corner_err_max
real(kind=8),          intent(out) :: gl_halo_corner_err, gl_halo_corner_err_max
integer(kind=4),       intent(in)  :: nh, halo_width, ntiles
type(grid_function_t), intent(in)  :: f1(ntiles), f2(ntiles)
!locals
real(kind=8)      err, inface_err, cross_edge_err
real(kind=8)      err_max, inface_err_max, cross_edge_err_max
real(kind=8)      inface_corner_err, inface_corner_err_max
real(kind=8)      inedge_corner_err, inedge_corner_err_max
real(kind=8)      halo_corner_err, halo_corner_err_max
integer(kind=4)   num_inedge_corners, gl_num_inedge_corners
integer(kind=4)   is, ie, js, je, ks, ke, klev, ind, ierr

inface_err = 0._8;          inface_err_max = 0._8
cross_edge_err = 0._8;      cross_edge_err_max = 0._8
inface_corner_err = 0._8;   inface_corner_err_max = 0._8
inedge_corner_err = 0._8;   inedge_corner_err_max = 0._8
halo_corner_err = 0._8;     halo_corner_err_max = 0._8
num_inedge_corners = 0

do ind = 1, ntiles

    is = f1(ind)%is; ie = f1(ind)%ie
    js = f1(ind)%js; je = f1(ind)%je
    ks = f1(ind)%ks; ke = f1(ind)%ke
    klev = ke-ks+1

    !tile edge errors
    err     = sum(abs(f1(ind)%p(is-halo_width:is-1,js:je,ks:ke)-f2(ind)%p(is-halo_width:is-1,js:je,ks:ke)))/nh
    err_max = maxval(abs(f1(ind)%p(is-halo_width:is-1,js:je,ks:ke)-f2(ind)%p(is-halo_width:is-1,js:je,ks:ke)))
    if(is == 1) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     = sum(abs(f1(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)-f2(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)))/nh
    err_max = maxval(abs(f1(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)-f2(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)))
    if(ie == nh) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     = sum(abs(f1(ind)%p(is:ie,js-halo_width:js-1,ks:ke)-f2(ind)%p(is:ie,js-halo_width:js-1,ks:ke)))/nh
    err_max = maxval(abs(f1(ind)%p(is:ie,js-halo_width:js-1,ks:ke)-f2(ind)%p(is:ie,js-halo_width:js-1,ks:ke)))
    if(js == 1) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     = sum(abs(f1(ind)%p(is:ie,je+1:je+halo_width,ks:ke)-f2(ind)%p(is:ie,je+1:je+halo_width,ks:ke)))/nh
    err_max = maxval(abs(f1(ind)%p(is:ie,je+1:je+halo_width,ks:ke)-f2(ind)%p(is:ie,je+1:je+halo_width,ks:ke)))
    if(je == nh) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     = sum(abs(f1(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke)-f2(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke))) / halo_width**2
    err_max = maxval(abs(f1(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke)-f2(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke)))
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

    err     = sum(abs(f1(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke)-f2(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke))) / halo_width**2
    err_max = maxval(abs(f1(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke)-f2(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke)))
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

    err     = sum(abs(f1(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke)-f2(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke))) / halo_width**2
    err_max = maxval(abs(f1(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke)-f2(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke)))
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

    err     = sum(abs(f1(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke)-f2(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke))) / halo_width**2
    err_max = maxval(abs(f1(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke)-f2(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke)))
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

end subroutine halo_err

end module test_halo_mod
