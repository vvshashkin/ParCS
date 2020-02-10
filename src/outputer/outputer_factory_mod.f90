module outputer_factory_mod

implicit none

contains

function create_master_paneled_outputer(master_id, gather_exch, partition) result(outputer)

    use master_paneled_outputer_mod, only : master_paneled_outputer_t
    use exchange_abstract_mod,       only : exchange_t
    use partition_mod,               only : partition_t
    use mpi

    integer(kind=4),   intent(in) :: master_id
    class(exchange_t), intent(in) :: gather_exch
    type(partition_t), intent(in) :: partition

    type(master_paneled_outputer_t) :: outputer

    integer(kind=4) :: i, myid, ierr

    call mpi_comm_rank(mpi_comm_world, myid, ierr)

    if (myid == master_id) then
        allocate(outputer%buf(1:partition%num_tiles*partition%num_panels))
        do i=1, partition%num_tiles*partition%num_panels
            call outputer%buf(i)%init( &
                 partition%tile(i)%panel_number,      &
                 partition%tile(i)%is,                &
                 partition%tile(i)%ie,                &
                 partition%tile(i)%js,                &
                 partition%tile(i)%je,                &
                 partition%tile(i)%ks,                &
                 partition%tile(i)%ke,                &
                 0,                                   &
                 0,                                   &
                 0)
        end do
    end if

    outputer%master_id = master_id
    outputer%gather_exch = gather_exch


end function create_master_paneled_outputer


function create_mpi_paneled_outputer(partition) result(outputer)

    use mpi_paneled_outputer_mod, only : mpi_paneled_outputer_t
    use partition_mod,            only : partition_t
    use mpi

    type(partition_t), intent(in)  :: partition
    type(mpi_paneled_outputer_t)   :: outputer

    integer(kind=mpi_offset_kind) :: disp
    integer(kind=4) :: type_local(partition%ts:partition%te)
    integer(kind=4) :: type_global(partition%ts:partition%te)
    integer(kind=4) :: type_domain
    integer(kind=4) :: dim, gsize(4), lsize(4), start(4), size
    integer(kind=4) :: myid, ierr
    integer(kind=4) :: t, is, ie, js, je, ks, ke

    call mpi_comm_rank(mpi_comm_world, myid, ierr)

    dim = 4

    gsize = [partition%Nh, partition%Nh, partition%num_panels, partition%Nz]
    lsize = [partition%Nh, partition%Nh, partition%num_panels, partition%Nz]
    start = [0,0,0,0]

    call MPI_Type_create_subarray(dim, gsize, lsize, start, MPI_ORDER_FORTRAN, mpi_real4, type_domain, ierr)
    call MPI_Type_commit(type_domain, ierr)
    call mpi_type_size(type_domain, size, ierr)
    call mpi_type_free(type_domain, ierr)

    disp = int(size,kind=mpi_offset_kind)

    do t = partition%ts, partition%te

        is = partition%tile(t)%is; ie = partition%tile(t)%ie
        js = partition%tile(t)%js; je = partition%tile(t)%je
        ks = partition%tile(t)%ks; ke = partition%tile(t)%ke

        gsize = [partition%Nh, partition%Nh, partition%num_panels, partition%Nz]
        lsize = [ie-is+1, je-js+1, 1, ke-ks+1]
        start = [is-1, js-1,partition%tile(t)%panel_number-1, ks-1]

        call MPI_Type_create_subarray(dim, gsize, lsize, start, MPI_ORDER_FORTRAN, mpi_real4, type_global(t), ierr)
        call MPI_Type_commit(type_global(t), ierr)

        gsize = [ie-is+1, je-js+1, 1, ke-ks+1]
        lsize = [ie-is+1, je-js+1, 1, ke-ks+1]
        start = [0, 0, 0, 0]

        call MPI_Type_create_subarray(dim, gsize, lsize, start, MPI_ORDER_FORTRAN, mpi_real4, type_local(t), ierr)
        call MPI_Type_commit(type_local(t), ierr)

    end do

    outputer = mpi_paneled_outputer_t(record_disp = disp,      &
                                      type_local = type_local, &
                                      type_global = type_global )

end function create_mpi_paneled_outputer

end module outputer_factory_mod
