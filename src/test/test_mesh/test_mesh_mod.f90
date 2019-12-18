module test_mesh_mod
contains

subroutine test_mesh()

    use mpi

    use partition_mod,        only : partition_t
    use mesh_factory_mod,     only : create_equiangular_mesh
    use mesh_mod,             only : mesh_t

    type(partition_t)                  :: partition
    type(mesh_t),          allocatable :: mesh(:)

    integer(kind=4)                    :: nh=100, nz=10, halo_width=50
    integer(kind=4)                    :: myid, np, ierr, code

    integer(kind=4) :: ts, te
    integer(kind=4) :: ind

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    if (myid==0) print*, 'Running mesh test!'

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
    end do

end subroutine test_mesh

end module test_mesh_mod
