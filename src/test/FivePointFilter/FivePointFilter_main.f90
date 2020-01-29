program FivePointFilter

use mpi
use exchange_mod,         only : exchange_t
use partition_mod,        only : partition_t
use grid_function_mod,    only : grid_function_t
use exchange_factory_mod, only : create_2d_cross_halo_exchange
use FivePointFilter_mod,  only : calc_fivepointfilter

implicit none

type(exchange_t)                   :: exch_halo1
type(partition_t)                  :: partition
type(grid_function_t), allocatable :: f1(:), f2(:)

integer(kind=4)                    :: nh=4, nz=1
integer(kind=4)                    :: myid, np, ierr

integer(kind=4) :: ts, te
integer(kind=4) :: ind, i, j, k


call MPI_init(ierr)
call MPI_comm_rank(mpi_comm_world , myid, ierr)
call MPI_comm_size(mpi_comm_world , Np  , ierr)

! if (mod(Np,6) /= 0) then
!     if (myid ==0) print*, 'Error! Wrong Np! Abort!'
!     call mpi_barrier(mpi_comm_world, ierr)
!     call mpi_finalize(ierr)
!     stop
! end if

call partition%init(nh, nz, max(1,Np/6), Np, strategy = 'default')

print*, 'NUM TILES', partition%num_tiles

print*, 'Mpi mode. Np = ', Np, 'myid = ', myid

    call mpi_barrier(mpi_comm_world, ierr)

if (myid ==0) print*, partition%proc_map

do i = 1, 6*partition%num_tiles
    call mpi_barrier(mpi_comm_world, ierr)

    if (partition%proc_map(i) == myid) then
        print*, 'myid = ', myid
        print*, 'tile number = ', i
        call partition%tile(i)%print()
        print*, ''
    end if
    call mpi_barrier(mpi_comm_world, ierr)
end do

if (myid == 0) call partition%write_to_txt('partition.txt')

!find start and end index of tiles belonging to the current proccesor
ts = findloc(partition%proc_map, myid, dim=1)
te = findloc(partition%proc_map, myid, back = .true., dim=1)


!Init arrays

allocate(f1(ts:te), f2(ts:te))

do i = ts, te

    call f1(i)%init(partition%tile(i)%panel_number,             &
                    partition%tile(i)%is, partition%tile(i)%ie, &
                    partition%tile(i)%js, partition%tile(i)%je, &
                    partition%tile(i)%ks, partition%tile(i)%ke, &
                    1, 1, 0)

    call f2(i)%init(partition%tile(i)%panel_number,             &
                    partition%tile(i)%is, partition%tile(i)%ie, &
                    partition%tile(i)%js, partition%tile(i)%je, &
                    partition%tile(i)%ks, partition%tile(i)%ke, &
                    1, 1, 0)
end do

do ind = ts, te
    do k = partition%tile(ind)%ks, partition%tile(ind)%ke
        do j = partition%tile(ind)%js, partition%tile(ind)%je
            do i = partition%tile(ind)%is, partition%tile(ind)%ie
                f1(ind).p(i,j,k) = 1000*partition%tile(ind)%panel_number + i+j+k
            end do
        end do
    end do
end do

!Init exchange

call create_2d_cross_halo_exchange(exch_halo1, partition, 1, myid, np)

!Perform exchange

call exch_halo1%do(f1, ts, te)

do ind = ts, te
    call calc_fivepointfilter(f1(ind), f2(ind), partition%tile(ind))
end do

call mpi_barrier(mpi_comm_world, ierr)
call mpi_finalize(ierr)

end program FivePointFilter
