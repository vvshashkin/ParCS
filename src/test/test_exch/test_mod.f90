module test_mod
implicit none

contains

subroutine test_tile_mod()

    use tile_mod, only : tile_t

    type(tile_t) :: tile
    integer(kind=4) :: is, ie, js, je, ks, ke, panel_number

    is = 1; ie = 10
    js = 1; je = 10
    ks = 1; ke = 10;
    panel_number = 1;

    call tile%init(is, ie, js, je, ks, ke, panel_number)

end subroutine test_tile_mod

subroutine test_partition_mod()

    use partition_mod, only : partition_t

    type(partition_t) :: partition

    ! call partition%init(1000,10,50,'default')

end subroutine test_partition_mod

subroutine test_exchange()

use mpi

use grid_function_mod, only : grid_function_t
use exchange_mod, only : exchange_t
use partition_mod, only : partition_t
use exchange_factory_mod, only : create_2d_cross_halo_exchange

type(partition_t) :: partition
type(exchange_t) :: exch

type(grid_function_t), allocatable :: array(:)

integer(kind=4) :: myid, Np, tile_s, tile_e
integer(kind=4) :: i, ierr, j

real(kind=8), allocatable :: f(:,:,:)


call MPI_init(ierr)
call MPI_comm_rank(mpi_comm_world , myid, ierr)
call MPI_comm_size(mpi_comm_world , Np  , ierr)

if (mod(Np,6) /= 0) then
    if (myid ==0) print*, 'Error! Wrong Np! Abort!'
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize(ierr)
    stop
end if

call partition%init(4, 1, Np/6, Np, 'default')

print*, 'Mpi mode. Np = ', Np, 'myid = ', myid

    call mpi_barrier(mpi_comm_world, ierr)

do i = 1, 6*partition%num_tiles
    call mpi_barrier(mpi_comm_world, ierr)

    if (partition%proc_map(i) == myid) then
        print*, 'myid = ', myid
        print*, 'tile number = ', i
        call partition%tile(i)%print()
    end if
    call mpi_barrier(mpi_comm_world, ierr)
end do

call create_2d_cross_halo_exchange(exch, partition, 1, myid, Np)

call mpi_barrier(mpi_comm_world, ierr)
tile_s = findloc(partition%proc_map, myid, dim=1)
tile_e = findloc(partition%proc_map, myid, back = .true., dim=1)

print*, myid, tile_s, tile_e

allocate(array(tile_s:tile_e))

do i = tile_s, tile_e

    call array(i)%init(partition%tile(i)%is, partition%tile(i)%ie, &
                       partition%tile(i)%js, partition%tile(i)%je, &
                       partition%tile(i)%ks, partition%tile(i)%ke, &
                       1, 1, 0)
    array(i)%p = myid

end do

call exch%do(array, tile_s, tile_e)

if ( myid == 11) then
    do j = partition%tile(tile_s)%js-1, partition%tile(tile_s)%je+1
     do    i = partition%tile(tile_s)%is-1, partition%tile(tile_s)%ie+1
    print*, j, i, array(tile_s).p(i,j,1)
end do
end do
end if

    if (myid ==0) call partition%write_to_txt('partition.txt')

call mpi_barrier(mpi_comm_world, ierr)
call mpi_finalize(ierr)

end subroutine test_exchange

end module test_mod
