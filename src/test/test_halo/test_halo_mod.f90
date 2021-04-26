module test_halo_mod

use grid_field_mod, only : grid_field_t
use mesh_mod,       only : mesh_t
use partition_mod,  only : partition_t

implicit none

private
public   :: test_ecs_halo

contains

subroutine test_ecs_halo()

    use mpi

    use grid_field_factory_mod,     only : create_grid_field
    use exchange_abstract_mod,      only : exchange_t
    use exchange_factory_mod,       only : create_Agrid_halo_exchange
    use mesh_factory_mod,           only : create_equiangular_mesh
    use ecs_halo_mod,               only : ecs_halo_t
    use ecs_halo_vec_a_mod,         only : ecs_halo_vec_t
    use ecs_halo_factory_mod,       only : init_ecs_halo
    use ecs_halo_vec_a_factory_mod, only : init_ecs_halo_vect


    class(exchange_t),     allocatable :: exch_halo
    type(partition_t)                  :: partition
    type(mesh_t),          allocatable :: mesh(:)
    type(grid_field_t)                 :: f_test, u_test, v_test
    type(grid_field_t)                 :: f_true, u_true, v_true

    type(ecs_halo_t),     allocatable :: halo(:)
    type(ecs_halo_vec_t), allocatable :: halo_vec(:)

    integer(kind=4), parameter         :: nh=128, nz=3, halo_width=3, ex_halo_width=8
    integer(kind=4)                    :: myid, np, ierr, code

    integer(kind=4) :: ts, te, t
    integer(kind=4) :: i, j, k, ifc, err_sum, gl_err_sum
    integer(kind=4) :: iev, isv, jev, jsv
    integer(kind=4) :: ie, is, js, je, klev, ke, ks

    real(kind=8) gl_inface_err, gl_cross_edge_err
    real(kind=8) gl_inface_err_max, gl_cross_edge_err_max
    real(kind=8) gl_inface_corner_err, gl_inface_corner_err_max
    real(kind=8) gl_inedge_corner_err, gl_inedge_corner_err_max
    real(kind=8) gl_halo_corner_err, gl_halo_corner_err_max
    integer(kind=4) num_inedge_corners, gl_num_inedge_corners

    logical is_test_passed
    real(kind=8), parameter :: tolerance = 0.3e-7_8

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    call partition%init(nh, nz, max(1,Np/6), myid, Np, strategy = 'default')
    !call partition%init(nh, nz, 64, Np, strategy = 'default')

    !find start and end index of tiles belonging to the current proccesor
    ts = partition%ts
    te = partition%te

    allocate(mesh(ts:te))
    do t = ts, te
        call create_equiangular_mesh(mesh(t), partition%tile(t)%is, partition%tile(t)%ie, &
                                              partition%tile(t)%js, partition%tile(t)%je, &
                                              partition%tile(t)%ks, partition%tile(t)%ke, &
                                              nh, ex_halo_width, partition%tile(t)%panel_number)
    end do

    allocate(halo(ts:te), halo_vec(ts:te))
    do t = ts, te
        halo(t) = init_ecs_halo(mesh(t)%is, mesh(t)%ie, &
                                mesh(t)%js, mesh(t)%je, &
                                mesh(t)%nx, halo_width, &
                                mesh(t)%hx)
        halo_vec(t) = init_ecs_halo_vect(mesh(t)%panel_ind,&
                                         mesh(t)%is, mesh(t)%ie,  &
                                         mesh(t)%js, mesh(t)%je,  &
                                         mesh(t)%nx, halo_width,    &
                                         mesh(t)%hx, halo(t))
    end do

    call mpi_barrier(mpi_comm_world, ierr)
    if (myid==0) print *, 'equiangular cubed-sphere halo-zone interpolation test'

    call create_grid_field(f_test, ex_halo_width, 0, partition)
    call create_grid_field(f_true, ex_halo_width, 0, partition)

    call create_grid_field(u_test, ex_halo_width, 0, partition)
    call create_grid_field(u_true, ex_halo_width, 0, partition)

    call create_grid_field(v_test, ex_halo_width, 0, partition)
    call create_grid_field(v_true, ex_halo_width, 0, partition)

    call init_scalar_halo_test_fun(f_test,         partition, mesh, ex_halo_width)
    call init_scalar_halo_test_fun(f_true,         partition, mesh, ex_halo_width)
    call init_vector_halo_test_fun(u_test, v_test, partition, mesh, ex_halo_width)
    call init_vector_halo_test_fun(u_true, v_true, partition, mesh, ex_halo_width)

    !Init exchange
    exch_halo = create_Agrid_halo_exchange(partition, ex_halo_width, 'full', myid, np)

    !Perform exchange
    call exch_halo%do(f_test)
    call exch_halo%do(u_test)
    call exch_halo%do(v_test)

    do t = ts, te
        call halo(t)%interp(f_test%block(t),halo_width)
        call halo_vec(t)%interpv(u_test%block(t),v_test%block(t),halo_width)
    end do

    call halo_err(gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max, &
                  gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err,       &
                  gl_inedge_corner_err_max, gl_halo_corner_err, gl_halo_corner_err_max,       &
                  nh, halo_width, f_test, f_true, ts, te)

    if(myid == 0) then
      print *, "halo errors: F"
      print *, "inface mean abs, inface max, cross edge mean abs, cross edge max"
      print '(4e15.7)', gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max
      print *, "inface-corner mean abs, inface-corner max, inedge-corner mean abs, inedge-corner max, halo-corner mean abs, halo-corner max"
      print '(6e15.7)', gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err, gl_inedge_corner_err_max,&
                        gl_halo_corner_err, gl_halo_corner_err_max

      is_test_passed = max(gl_inedge_corner_err_max, gl_halo_corner_err_max, &
                           gl_cross_edge_err_max) < tolerance
    end if
    call halo_err(gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max, &
                  gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err,       &
                  gl_inedge_corner_err_max, gl_halo_corner_err, gl_halo_corner_err_max,       &
                  nh, halo_width, u_test, u_true, ts, te)

    if(myid == 0) then
      print *, "halo errors: U"
      print *, "inface mean abs, inface max, cross edge mean abs, cross edge max"
      print '(4e15.7)', gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max
      print *, "inface-corner mean abs, inface-corner max, inedge-corner mean abs, inedge-corner max, halo-corner mean abs, halo-corner max"
      print '(6e15.7)', gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err, gl_inedge_corner_err_max,&
                        gl_halo_corner_err, gl_halo_corner_err_max
      is_test_passed = is_test_passed .and. max(gl_inedge_corner_err_max, gl_halo_corner_err_max, &
                                           gl_cross_edge_err_max) < tolerance
    end if

    call halo_err(gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max, &
                  gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err,       &
                  gl_inedge_corner_err_max, gl_halo_corner_err, gl_halo_corner_err_max,       &
                  nh, halo_width, v_test, v_true, ts, te)

    if(myid == 0) then
      print *, "halo errors: V"
      print *, "inface mean abs, inface max, cross edge mean abs, cross edge max"
      print '(4e15.7)', gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max
      print *, "inface-corner mean abs, inface-corner max, inedge-corner mean abs, inedge-corner max, halo-corner mean abs, halo-corner max"
      print '(6e15.7)', gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err, gl_inedge_corner_err_max,&
                        gl_halo_corner_err, gl_halo_corner_err_max
      is_test_passed = is_test_passed .and. max(gl_inedge_corner_err_max, gl_halo_corner_err_max, &
                                                gl_cross_edge_err_max) < tolerance

      if (is_test_passed) then
          print *, "halo test passed"
      else
          print *, "halo test failed"
      end if
    end if

end subroutine test_ecs_halo

subroutine init_scalar_halo_test_fun(f, partition, mesh, halo_width)

type(grid_field_t), intent(inout) :: f

type(partition_t), intent(in) :: partition
type(mesh_t),      intent(in) :: mesh(partition%ts:partition%te)
integer(kind=4),   intent(in) :: halo_width

!locals
integer(kind=4) t, isv, iev, jsv, jev

do t = partition%ts, partition%te
     isv = mesh(t)%is-halo_width
     iev = mesh(t)%ie+halo_width
     jsv = mesh(t)%js-halo_width
     jev = mesh(t)%je+halo_width
     f%block(t)%p(isv:iev,jsv:jev,1) = mesh(t)%rhx(isv:iev,jsv:jev)
     f%block(t)%p(isv:iev,jsv:jev,2) = mesh(t)%rhy(isv:iev,jsv:jev)
     f%block(t)%p(isv:iev,jsv:jev,3) = mesh(t)%rhz(isv:iev,jsv:jev)
end do

end subroutine init_scalar_halo_test_fun

subroutine init_vector_halo_test_fun(u, v, partition, mesh, halo_width)

type(grid_field_t), intent(inout) :: u, v

type(partition_t), intent(in) :: partition
type(mesh_t),      intent(in) :: mesh(partition%ts:partition%te)
integer(kind=4),   intent(in) :: halo_width

!locals
integer(kind=4) :: isv, iev, jsv, jev, i, j, t
real(kind=8)    ::  vx, vy, vz

do t = partition%ts, partition%te

     isv = mesh(t)%is-halo_width
     iev = mesh(t)%ie+halo_width
     jsv = mesh(t)%js-halo_width
     jev = mesh(t)%je+halo_width
     do j = jsv, jev
         do i = isv, iev
             !level 1: solid rotation around (1,0,0) axis
             vx = 0._8; vy = mesh(t)%rhz(i,j); vz =-mesh(t)%rhy(i,j)
             u%block(t)%p(i,j,1) = sum([vx,vy,vz]*mesh(t)%actv(1:3,i,j))
             v%block(t)%p(i,j,1) = sum([vx,vy,vz]*mesh(t)%bctv(1:3,i,j))
              !level 2: solid rotation around (0,1,0) axis
              vx =-mesh(t)%rhz(i,j); vy = 0._8; vz = mesh(t)%rhx(i,j)
              u%block(t)%p(i,j,2) = sum([vx,vy,vz]*mesh(t)%actv(1:3,i,j))
              v%block(t)%p(i,j,2) = sum([vx,vy,vz]*mesh(t)%bctv(1:3,i,j))
             !level 3: solid rotation around (0,0,1) axis
             vx = mesh(t)%rhy(i,j); vy =-mesh(t)%rhx(i,j); vz = 0._8
             u%block(t)%p(i,j,3) = sum([vx,vy,vz]*mesh(t)%actv(1:3,i,j))
             v%block(t)%p(i,j,3) = sum([vx,vy,vz]*mesh(t)%bctv(1:3,i,j))
         end do
     end do
end do

end subroutine init_vector_halo_test_fun


subroutine halo_err(&
              gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max, &
              gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err,       &
              gl_inedge_corner_err_max, gl_halo_corner_err, gl_halo_corner_err_max,       &
              nh, halo_width, f1, f2, ts ,te)

use mpi
use grid_function_mod, only: grid_function_t

real(kind=8),          intent(out) :: gl_inface_err, gl_cross_edge_err
real(kind=8),          intent(out) :: gl_inface_err_max, gl_cross_edge_err_max
real(kind=8),          intent(out) :: gl_inface_corner_err, gl_inface_corner_err_max
real(kind=8),          intent(out) :: gl_inedge_corner_err, gl_inedge_corner_err_max
real(kind=8),          intent(out) :: gl_halo_corner_err, gl_halo_corner_err_max
integer(kind=4),       intent(in)  :: nh, halo_width, ts, te
type(grid_field_t),    intent(in)  :: f1, f2
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

do ind = ts, te

    is = f1%block(ind)%is; ie = f1%block(ind)%ie
    js = f1%block(ind)%js; je = f1%block(ind)%je
    ks = f1%block(ind)%ks; ke = f1%block(ind)%ke
    klev = ke-ks+1

    !tile edge errors
    err     = sum(abs(f1%block(ind)%p(is-halo_width:is-1,js:je,ks:ke)-f2%block(ind)%p(is-halo_width:is-1,js:je,ks:ke)))/nh
    err_max = maxval(abs(f1%block(ind)%p(is-halo_width:is-1,js:je,ks:ke)-f2%block(ind)%p(is-halo_width:is-1,js:je,ks:ke)))
    if(is == 1) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     =    sum(abs(f1%block(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)-f2%block(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)))/nh
    err_max = maxval(abs(f1%block(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)-f2%block(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)))
    if(ie == nh) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     =    sum(abs(f1%block(ind)%p(is:ie,js-halo_width:js-1,ks:ke)-f2%block(ind)%p(is:ie,js-halo_width:js-1,ks:ke)))/nh
    err_max = maxval(abs(f1%block(ind)%p(is:ie,js-halo_width:js-1,ks:ke)-f2%block(ind)%p(is:ie,js-halo_width:js-1,ks:ke)))
    if(js == 1) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     =    sum(abs(f1%block(ind)%p(is:ie,je+1:je+halo_width,ks:ke)-f2%block(ind)%p(is:ie,je+1:je+halo_width,ks:ke)))/nh
    err_max = maxval(abs(f1%block(ind)%p(is:ie,je+1:je+halo_width,ks:ke)-f2%block(ind)%p(is:ie,je+1:je+halo_width,ks:ke)))
    if(je == nh) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err/nh
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     = sum(abs(f1%block(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke)-f2%block(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke))) / halo_width**2
    err_max = maxval(abs(f1%block(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke)-f2%block(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke)))
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

    err     = sum(abs(f1%block(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke)-f2%block(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke))) / halo_width**2
    err_max = maxval(abs(f1%block(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke)-f2%block(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke)))
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

    err     = sum(abs(f1%block(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke)-f2%block(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke))) / halo_width**2
    err_max = maxval(abs(f1%block(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke)-f2%block(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke)))
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

    err     = sum(abs(f1%block(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke)-f2%block(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke))) / halo_width**2
    err_max = maxval(abs(f1%block(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke)-f2%block(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke)))
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
