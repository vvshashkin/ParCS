module test_output_mod
contains

subroutine test_output()

    use mpi

    use grid_function_mod,           only : grid_function_t
    use partition_mod,               only : partition_t
    use exchange_factory_mod,        only : create_gather_exchange
    use outputer_abstract_mod,       only : outputer_t
    use outputer_factory_mod,        only : create_master_process_outputer
    use mesh_mod,                    only : mesh_t

    type(partition_t)                  :: partition
    class(outputer_t), allocatable     :: outputer_bin, outputer_txt
    type(grid_function_t), allocatable :: f1(:)
    type(mesh_t), allocatable          :: mesh(:)

    integer(kind=4)                    :: nh=20, nz=2, halo_width=10
    integer(kind=4)                    :: myid, np, ierr, code
    integer(kind=4)                    :: master_id = 0

    integer(kind=4) :: ts, te
    integer(kind=4) :: ind, i, j, k

    real(kind=8) :: err_sum

    integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)


    if (myid==0) print*, 'Running master process output test!'

    call partition%init(nh, nz, max(1,Np/6), Np, strategy = 'default')

    !find start and end index of tiles belonging to the current proccesor
    ts = findloc(partition%proc_map, myid, dim=1)
    te = findloc(partition%proc_map, myid, back = .true., dim=1)

    !Init arrays


    if (myid == master_id) then
        allocate( f1(size(partition%tile,1)) )
        allocate( mesh(size(partition%tile,1)) )
        do i = 1, size(partition%tile,1)
            call f1(i)%init(partition%tile(i)%panel_number,             &
                            partition%tile(i)%is, partition%tile(i)%ie, &
                            partition%tile(i)%js, partition%tile(i)%je, &
                            partition%tile(i)%ks, partition%tile(i)%ke, &
                            halo_width, halo_width, 0)
            f1(i).p = huge(1.0_8)
        end do
    else
        allocate(f1(ts:te))
        allocate(mesh(ts:te))
        do i = ts, te
            call f1(i)%init(partition%tile(i)%panel_number,             &
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


    outputer_bin = create_master_process_outputer(master_id = 0, &
                   gather_exch =create_gather_exchange(partition, master_id, myid, np), &
                   write_type = 'bin')

    outputer_txt = create_master_process_outputer(master_id = 0, &
                   gather_exch =create_gather_exchange(partition, master_id, myid, np), &
                   write_type = 'txt')

    do ind = 1, 5

        call outputer_bin%write(f1, mesh, lbound(f1, dim=1), ubound(f1, dim=1), 'test_output')
        call outputer_txt%write(f1, mesh, lbound(f1, dim=1), ubound(f1, dim=1), 'test_output')

    end do

    call mpi_barrier(mpi_comm_world, ierr)

end subroutine test_output

end module test_output_mod
