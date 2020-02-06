module test_paneled_output_mod
contains

subroutine test_paneled_output()

    use mpi

    use grid_function_mod,           only : grid_function_t
    use partition_mod,               only : partition_t
    use exchange_factory_mod,        only : create_gather_exchange
    use outputer_abstract_mod,       only : outputer_t
    use outputer_factory_mod,        only : create_master_paneled_outputer
    use mesh_mod,                    only : mesh_t
    use mesh_factory_mod,            only : create_equiangular_mesh

    type(partition_t)                  :: partition
    class(outputer_t), allocatable     :: outputer
    type(grid_function_t), allocatable :: f1(:)
    type(mesh_t), allocatable          :: mesh(:)
    type(mesh_t)                       :: mesh_global

    integer(kind=4)                    :: nh=128, nz=3, halo_width=3
    integer(kind=4)                    :: myid, np, ierr, code
    integer(kind=4)                    :: master_id = 0

    character(*), parameter            :: file_name = "h.dat"
    integer(kind=4)                    :: fdunit

    integer(kind=4) :: ts, te, ntiles
    integer(kind=4) :: ind, i, j, k

    real(kind=4), allocatable :: bufcheck(:,:,:,:), bufin(:,:,:,:)
    real(kind=4)  err
    real(kind=4), parameter :: tolerance = 1e-16

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    if (myid==0) print*, 'Running master process paneled output test!'

    call partition%init(nh, nz, max(1,Np/6), Np, strategy = 'default')

    !find start and end index of tiles belonging to the current proccesor
    ts = findloc(partition%proc_map, myid, dim=1)
    te = findloc(partition%proc_map, myid, back = .true., dim=1)
    ntiles = size(partition%tile,1)

    !Init arrays

    if (myid == master_id) then

        allocate( f1(ntiles) )
        allocate( mesh(ntiles) )

        do i = 1, ntiles
            call create_equiangular_mesh(mesh(i), partition%tile(i)%is, partition%tile(i)%ie, &
                                         partition%tile(i)%js, partition%tile(i)%je, &
                                         partition%tile(i)%ks, partition%tile(i)%ke, &
                                         nh, halo_width, partition%tile(i)%panel_number)
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
            call create_equiangular_mesh(mesh(i), partition%tile(i)%is, partition%tile(i)%ie, &
                                         partition%tile(i)%js, partition%tile(i)%je, &
                                         partition%tile(i)%ks, partition%tile(i)%ke, &
                                         nh, halo_width, partition%tile(i)%panel_number)
            call f1(i)%init(partition%tile(i)%panel_number,             &
                            partition%tile(i)%is, partition%tile(i)%ie, &
                            partition%tile(i)%js, partition%tile(i)%je, &
                            partition%tile(i)%ks, partition%tile(i)%ke, &
                            halo_width, halo_width, 0)
            f1(i).p = huge(1.0_8)
        end do
    end if

!    do ind = ts, te
!        do k = partition%tile(ind)%ks, partition%tile(ind)%ke
!            do j = partition%tile(ind)%js, partition%tile(ind)%je
!                do i = partition%tile(ind)%is, partition%tile(ind)%ie
!                    f1(ind).p(i,j,k) =  f1(ind)%panel_ind*10+k
!                end do
!            end do
!        end do
!    end do
    do ind = ts, te
        do j = partition%tile(ind)%js, partition%tile(ind)%je
            do i = partition%tile(ind)%is, partition%tile(ind)%ie
                f1(ind).p(i,j,1) =  mesh(ind)%rhx(i,j)
                f1(ind).p(i,j,2) =  mesh(ind)%rhy(i,j)
                f1(ind).p(i,j,3) =  mesh(ind)%rhz(i,j)
            end do
        end do
    end do

    outputer = create_master_paneled_outputer( master_id = 0, &
    gather_exch = create_gather_exchange(partition, master_id, myid, np) )

    call outputer%write(f1, mesh, lbound(f1, dim=1), ubound(f1, dim=1), file_name)

    call mpi_barrier(mpi_comm_world, ierr)

    if(myid == 0) then
        allocate(bufcheck(nh,nh,6,nz))
        allocate(   bufin(nh,nh,6,nz))

        do ind = 1, 6
            call create_equiangular_mesh(mesh_global, 1, nh, 1, nh, & !horizontal dimensions
                                                      1, 1, nh,0, ind)
            bufcheck(1:nh,1:nh,ind,1) = mesh_global%rhx(1:nh,1:nh)
            bufcheck(1:nh,1:nh,ind,2) = mesh_global%rhy(1:nh,1:nh)
            bufcheck(1:nh,1:nh,ind,3) = mesh_global%rhz(1:nh,1:nh)
        end do

        open(file="h.dat",newunit = fdunit, access="direct", recl = 6*nh*nh*nz)
        read(fdunit,rec=1) bufin

        err = maxval(abs(bufcheck-bufin))
        print *, "discrepancy", err
        if( err < tolerance) then
            print *, "paneled master-process output test passed"
        else
            print *, "paneled master-process output test failed"
        end if
    end if

end subroutine test_paneled_output

end module test_paneled_output_mod
