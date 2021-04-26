module test_paneled_output_mod

implicit none

contains

subroutine test_master_paneled_output()

    use exchange_factory_mod,        only : create_gather_exchange
    use outputer_abstract_mod,       only : outputer_t
    use outputer_factory_mod,        only : create_master_paneled_outputer
    use partition_mod,               only : partition_t
    use mpi

    class(outputer_t),allocatable :: outputer
    type(partition_t) :: partition

    integer(kind=4)                    :: nh=128, nz=3, halo_width=3
    integer(kind=4)                    :: myid, np, ierr, code
    integer(kind=4)                    :: master_id = 0

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    if (myid==0) print*, 'Running master process paneled output test!'

    call partition%init(nh, nz, max(1,Np/6), myid, Np, strategy = 'default')

    outputer = create_master_paneled_outputer( &
    master_id   = 0                                                      , &
    gather_exch = create_gather_exchange(partition, master_id, myid, np) , &
    partition   = partition)

    call test_case_1(outputer, partition, 'master_outputer')

end subroutine test_master_paneled_output

subroutine test_mpi_paneled_output()

    use outputer_abstract_mod,       only : outputer_t
    use outputer_factory_mod,        only : create_mpi_paneled_outputer
    use partition_mod,               only : partition_t
    use mpi

    class(outputer_t),allocatable :: outputer
    type(partition_t) :: partition

    integer(kind=4)                    :: nh=128, nz=3, halo_width=3
    integer(kind=4)                    :: myid, np, ierr, code
    integer(kind=4)                    :: master_id = 0

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    if (myid==0) print*, 'Running mpi paneled output test!'

    call partition%init(nh, nz, max(1,Np/6), myid, Np, strategy = 'default')

    outputer = create_mpi_paneled_outputer(partition)

    call test_case_1(outputer, partition, 'mpi_outputer')

end subroutine test_mpi_paneled_output

subroutine test_case_1(outputer, partition, class_name)

    use mpi

    use grid_field_mod,              only : grid_field_t
    use grid_field_factory_mod,      only : create_grid_field
    use partition_mod,               only : partition_t
    use outputer_abstract_mod,       only : outputer_t
    use mesh_mod,                    only : mesh_t
    use mesh_factory_mod,            only : create_equiangular_mesh

    class(outputer_t), intent(inout)   :: outputer
    type(partition_t), intent(in)      :: partition
    character(*),      intent(in)      :: class_name

    type(grid_field_t) :: f

    type(mesh_t), allocatable          :: mesh(:)
    type(mesh_t)                       :: mesh_global

    character(*), parameter            :: file_name = "h.dat"
    integer(kind=4)                    :: fdunit

    integer(kind=4)                    :: halo_width=3
    integer(kind=4)                    :: nh, nz

    integer(kind=4) :: myid, np, ierr, code
    integer(kind=4) :: ts, te, ntiles
    integer(kind=4) :: ind, t, i, j, k

    real(kind=4), allocatable :: bufcheck(:,:,:,:), bufin(:,:,:,:)
    real(kind=4)  err
    real(kind=4), parameter :: tolerance = 1e-16

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    !start and end index of tiles belonging to the current proccesor
    ts = partition%ts
    te = partition%te
    ntiles = partition%num_panels*partition%num_tiles

    !Init arrays

    nh=partition%Nh; nz=partition%Nz

    allocate(mesh(ts:te))
    do i = ts, te
        call create_equiangular_mesh(mesh(i), partition%tile(i)%is, partition%tile(i)%ie, &
                                     partition%tile(i)%js, partition%tile(i)%je, &
                                     partition%tile(i)%ks, partition%tile(i)%ke, &
                                     nh, halo_width, partition%tile(i)%panel_number)
    end do

    call create_grid_field(f, halo_width, 0, mesh, partition)

    do t = ts, te
        f%block(t)%p = huge(1.0_8)
        do j = partition%tile(t)%js, partition%tile(t)%je
            do i = partition%tile(t)%is, partition%tile(t)%ie
                f%block(t)%p(i,j,1) =  mesh(t)%rhx(i,j)
                f%block(t)%p(i,j,2) =  mesh(t)%rhy(i,j)
                f%block(t)%p(i,j,3) =  mesh(t)%rhz(i,j)
            end do
        end do
    end do

    call outputer%write(f, partition, file_name)

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
            print *, class_name, " passed output testcase 1"
        else
            print *, class_name, " failed output testcase 1"
        end if
    end if

end subroutine test_case_1

end module test_paneled_output_mod
