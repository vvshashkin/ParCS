module test_mod
implicit none

contains

subroutine test_cross_halo_exchange()

use mpi

use grid_function_mod,    only : grid_function_t
use exchange_mod,         only : exchange_t
use partition_mod,        only : partition_t
use exchange_factory_mod, only : create_2d_cross_halo_exchange

type(exchange_t)                   :: exch_halo
type(partition_t)                  :: partition
type(grid_function_t), allocatable :: f1(:)

integer(kind=4)                    :: nh=100, nz=10, halo_width=50
integer(kind=4)                    :: myid, np, ierr, code

integer(kind=4) :: ts, te
integer(kind=4) :: ind, i, j, k, err_sum, gl_err_sum

integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number

call MPI_comm_rank(mpi_comm_world , myid, ierr)
call MPI_comm_size(mpi_comm_world , Np  , ierr)

if (myid==0) print*, 'Running cross_halo_exchange test!'

call partition%init(nh, nz, max(1,Np/6), Np, strategy = 'default')

!find start and end index of tiles belonging to the current proccesor
ts = findloc(partition%proc_map, myid, dim=1)
te = findloc(partition%proc_map, myid, back = .true., dim=1)

!Init arrays

allocate(f1(ts:te))

do i = ts, te

    call f1(i)%init(partition%tile(i)%panel_number,             &
                    partition%tile(i)%is, partition%tile(i)%ie, &
                    partition%tile(i)%js, partition%tile(i)%je, &
                    partition%tile(i)%ks, partition%tile(i)%ke, &
                    halo_width, halo_width, 0)
    f1(i).p(:,:,:) = huge(1.0_8)
end do

do ind = ts, te
    do k = partition%tile(ind)%ks, partition%tile(ind)%ke
        do j = partition%tile(ind)%js, partition%tile(ind)%je
            do i = partition%tile(ind)%is, partition%tile(ind)%ie
                f1(ind).p(i,j,k) =  (partition%tile(ind)%panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k
            end do
        end do
    end do
end do

!Init exchange
call create_2d_cross_halo_exchange(exch_halo, partition, halo_width, myid, np)

!Perform exchange
call exch_halo%do(f1, lbound(f1,1), ubound(f1,1))


call mpi_barrier(mpi_comm_world, ierr)

err_sum = 0

do ind = 1, exch_halo%profile%recv_exch_num

    local_tile_ind           = exch_halo%profile%recv_to_tile_ind(ind)
    remote_tile_ind          = exch_halo%profile%recv_from_tile_ind(ind)

    remote_tile_panel_number = partition%tile(remote_tile_ind)%panel_number
    local_tile_panel_number  = partition%tile(local_tile_ind )%panel_number

    do k = exch_halo%profile%recv_ks(ind), exch_halo%profile%recv_ke(ind)
        do j = exch_halo%profile%recv_js(ind), exch_halo%profile%recv_je(ind)
            do i = exch_halo%profile%recv_is(ind), exch_halo%profile%recv_ie(ind)

                if (local_tile_panel_number == remote_tile_panel_number) then
                    err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k)))
                else
                    if (exch_halo%profile%send_j_step(ind)==1 .and. exch_halo%profile%send_i_step(ind)==1 .and. exch_halo%profile%first_dim_index(ind)=='i') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(modulo(j-1,nh)+1-1) + nz*(modulo(i-1,nh)+1-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==1 .and. exch_halo%profile%send_i_step(ind)==1 .and. exch_halo%profile%first_dim_index(ind)=='j') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(modulo(i-1,nh)+1-1) + nz*(modulo(j-1,nh)+1-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==1 .and. exch_halo%profile%send_i_step(ind)==-1 .and. exch_halo%profile%first_dim_index(ind)=='i') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(modulo(j-1,nh)+1-1) + nz*(nh-modulo(i-1,nh)-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==1 .and. exch_halo%profile%send_i_step(ind)==-1 .and. exch_halo%profile%first_dim_index(ind)=='j') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*(modulo(j-1,nh)+1-1) + nz*nh*(nh-modulo(i-1,nh)-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==-1 .and. exch_halo%profile%send_i_step(ind)==1 .and. exch_halo%profile%first_dim_index(ind)=='i') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(nh-modulo(j-1,nh)-1) + nz*(modulo(i-1,nh)+1-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==-1 .and. exch_halo%profile%send_i_step(ind)==1 .and. exch_halo%profile%first_dim_index(ind)=='j') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*(nh-modulo(j-1,nh)-1) + nz*nh*(modulo(i-1,nh)+1-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==-1 .and. exch_halo%profile%send_i_step(ind)==-1 .and. exch_halo%profile%first_dim_index(ind)=='i') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(nh-modulo(j-1,nh)-1) + nz*(nh-modulo(i-1,nh)-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==-1 .and. exch_halo%profile%send_i_step(ind)==-1 .and. exch_halo%profile%first_dim_index(ind)=='j') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*(nh-modulo(j-1,nh)-1) + nz*nh*(nh-modulo(i-1,nh)-1) + k)))
                    else
                        print*, 'Error!!!', 'myid=', myid
                        call mpi_abort(mpi_comm_world, code, ierr)
                        stop
                    end if
                end if
            end do
        end do
    end do
end do

call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)

if (gl_err_sum==0) then
    if (myid==0) print*, 'Test passed!'
else
    if (myid==0) print*, 'Test not passed! Error! Abort!'
    call mpi_abort(mpi_comm_world, code, ierr)
    stop
end if

end subroutine test_cross_halo_exchange

subroutine test_full_halo_exchange()

use mpi

use grid_function_mod,    only : grid_function_t
use exchange_mod,         only : exchange_t
use partition_mod,        only : partition_t
use exchange_factory_mod, only : create_2d_full_halo_exchange

type(exchange_t)                   :: exch_halo
type(partition_t)                  :: partition
type(grid_function_t), allocatable :: f1(:)

integer(kind=4)                    :: nh=100, nz=10, halo_width=50
integer(kind=4)                    :: myid, np, ierr, code

integer(kind=4) :: ts, te
integer(kind=4) :: ind, i, j, k, err_sum, gl_err_sum

integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number

call MPI_comm_rank(mpi_comm_world , myid, ierr)
call MPI_comm_size(mpi_comm_world , Np  , ierr)


if (myid==0) print*, 'Running full_halo_exchange test!'

call partition%init(nh, nz, max(1,Np/6), Np, strategy = 'default')

!find start and end index of tiles belonging to the current proccesor
ts = findloc(partition%proc_map, myid, dim=1)
te = findloc(partition%proc_map, myid, back = .true., dim=1)

!Init arrays

allocate(f1(ts:te))

do i = ts, te
    call f1(i)%init(partition%tile(ind)%panel_number,           &
                    partition%tile(i)%is, partition%tile(i)%ie, &
                    partition%tile(i)%js, partition%tile(i)%je, &
                    partition%tile(i)%ks, partition%tile(i)%ke, &
                    halo_width, halo_width, 0)
    f1(i).p(:,:,:) = huge(1.0_8)
end do

do ind = ts, te
    do k = partition%tile(ind)%ks, partition%tile(ind)%ke
        do j = partition%tile(ind)%js, partition%tile(ind)%je
            do i = partition%tile(ind)%is, partition%tile(ind)%ie
                f1(ind).p(i,j,k) =  (partition%tile(ind)%panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k
            end do
        end do
    end do
end do

!Init exchange
call create_2d_full_halo_exchange(exch_halo, partition, halo_width, myid, np)

!Perform exchange
call exch_halo%do(f1, lbound(f1,1), ubound(f1,1))

call mpi_barrier(mpi_comm_world, ierr)

err_sum = 0

do ind = 1, exch_halo%profile%recv_exch_num

    local_tile_ind           = exch_halo%profile%recv_to_tile_ind(ind)
    remote_tile_ind          = exch_halo%profile%recv_from_tile_ind(ind)

    remote_tile_panel_number = partition%tile(remote_tile_ind)%panel_number
    local_tile_panel_number  = partition%tile(local_tile_ind )%panel_number

    do k = exch_halo%profile%recv_ks(ind), exch_halo%profile%recv_ke(ind)
        do j = exch_halo%profile%recv_js(ind), exch_halo%profile%recv_je(ind)
            do i = exch_halo%profile%recv_is(ind), exch_halo%profile%recv_ie(ind)
                if (local_tile_panel_number == remote_tile_panel_number) then
                    err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k)))
                else
                    if (exch_halo%profile%send_j_step(ind)==1 .and. exch_halo%profile%send_i_step(ind)==1 .and. exch_halo%profile%first_dim_index(ind)=='i') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(modulo(j-1,nh)+1-1) + nz*(modulo(i-1,nh)+1-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==1 .and. exch_halo%profile%send_i_step(ind)==1 .and. exch_halo%profile%first_dim_index(ind)=='j') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(modulo(i-1,nh)+1-1) + nz*(modulo(j-1,nh)+1-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==1 .and. exch_halo%profile%send_i_step(ind)==-1 .and. exch_halo%profile%first_dim_index(ind)=='i') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(modulo(j-1,nh)+1-1) + nz*(nh-modulo(i-1,nh)-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==1 .and. exch_halo%profile%send_i_step(ind)==-1 .and. exch_halo%profile%first_dim_index(ind)=='j') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*(modulo(j-1,nh)+1-1) + nz*nh*(nh-modulo(i-1,nh)-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==-1 .and. exch_halo%profile%send_i_step(ind)==1 .and. exch_halo%profile%first_dim_index(ind)=='i') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(nh-modulo(j-1,nh)-1) + nz*(modulo(i-1,nh)+1-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==-1 .and. exch_halo%profile%send_i_step(ind)==1 .and. exch_halo%profile%first_dim_index(ind)=='j') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*(nh-modulo(j-1,nh)-1) + nz*nh*(modulo(i-1,nh)+1-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==-1 .and. exch_halo%profile%send_i_step(ind)==-1 .and. exch_halo%profile%first_dim_index(ind)=='i') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(nh-modulo(j-1,nh)-1) + nz*(nh-modulo(i-1,nh)-1) + k)))
                    else if (exch_halo%profile%send_j_step(ind)==-1 .and. exch_halo%profile%send_i_step(ind)==-1 .and. exch_halo%profile%first_dim_index(ind)=='j') then
                        err_sum = err_sum + abs(int(f1(local_tile_ind).p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*(nh-modulo(j-1,nh)-1) + nz*nh*(nh-modulo(i-1,nh)-1) + k)))
                    else
                        print*, 'Error!!!', 'myid=', myid
                        call mpi_abort(mpi_comm_world, code, ierr)
                        stop
                    end if
                end if
            end do
        end do
    end do
end do

call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)

if (gl_err_sum==0) then
    if (myid==0) print*, 'Test passed!'
else
    if (myid==0) print*, 'Test not passed! Error! Abort!'
    call mpi_abort(mpi_comm_world, code, ierr)
    stop
end if

end subroutine test_full_halo_exchange

subroutine test_gather_exchange()

use mpi

use grid_function_mod,    only : grid_function_t
use exchange_mod,         only : exchange_t
use partition_mod,        only : partition_t
use exchange_factory_mod, only : create_gather_exchange

type(exchange_t)                   :: exch_gather
type(partition_t)                  :: partition
type(grid_function_t), allocatable :: f1(:)

integer(kind=4)                    :: nh=100, nz=10, halo_width=10
integer(kind=4)                    :: myid, np, ierr, code
integer(kind=4)                    :: master_id = 0

integer(kind=4) :: ts, te
integer(kind=4) :: ind, i, j, k

real(kind=8) :: err_sum

integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number

call MPI_comm_rank(mpi_comm_world , myid, ierr)
call MPI_comm_size(mpi_comm_world , Np  , ierr)


if (myid==0) print*, 'Running gather_exchange test!'

call partition%init(nh, nz, max(1,Np/6), Np, strategy = 'default')

!find start and end index of tiles belonging to the current proccesor
ts = findloc(partition%proc_map, myid, dim=1)
te = findloc(partition%proc_map, myid, back = .true., dim=1)

!Init arrays


if (myid == master_id) then
    allocate(f1(1:6*partition%num_tiles))
    do i = 1, 6*partition%num_tiles
        call f1(i)%init(partition%tile(ind)%panel_number,           &
                        partition%tile(i)%is, partition%tile(i)%ie, &
                        partition%tile(i)%js, partition%tile(i)%je, &
                        partition%tile(i)%ks, partition%tile(i)%ke, &
                        halo_width, halo_width, 0)
        f1(i).p = huge(1.0_8)
    end do
else
    allocate(f1(ts:te))
    do i = ts, te
        call f1(i)%init(partition%tile(ind)%panel_number,           &
                        partition%tile(i)%is, partition%tile(i)%ie, &
                        partition%tile(i)%js, partition%tile(i)%je, &
                        partition%tile(i)%ks, partition%tile(i)%ke, &
                        halo_width, halo_width, 0)
        f1(i).p = huge(1.0_8)
    end do
end if

do ind = ts, te
    do k = partition%tile(ind)%ks, partition%tile(ind)%ke
        do j = partition%tile(ind)%js, partition%tile(ind)%je
            do i = partition%tile(ind)%is, partition%tile(ind)%ie
                f1(ind).p(i,j,k) =  (partition%tile(ind)%panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k
            end do
        end do
    end do
end do

!Init exchange
call create_gather_exchange(exch_gather, partition, master_id, myid, np)

!Perform exchange
call exch_gather%do(f1, lbound(f1,1), ubound(f1,1))

call mpi_barrier(mpi_comm_world, ierr)

if (myid == master_id) then

    err_sum = 0

    do ind = 1, 6*partition%num_tiles
        do k = partition%tile(ind)%ks, partition%tile(ind)%ke
            do j = partition%tile(ind)%js, partition%tile(ind)%je
                do i = partition%tile(ind)%is, partition%tile(ind)%ie
                    err_sum = err_sum + abs(f1(ind).p(i,j,k) - ((partition%tile(ind)%panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k))
                end do
            end do
        end do
    end do

    if (int(err_sum)==0) then
        print*, 'Test passed!'
    else
        print*, 'Test not passed! Error! Abort!'
        call mpi_abort(mpi_comm_world, code, ierr)
        stop
    end if

end if

end subroutine test_gather_exchange


end module test_mod
