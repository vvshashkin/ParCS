module test_paneled_output_mod

implicit none

contains

subroutine test_master_paneled_output()

    use exchange_factory_mod,  only : create_gather_exchange
    use outputer_abstract_mod, only : outputer_t
    use outputer_factory_mod,  only : create_master_paneled_outputer
    use partition_mod,         only : partition_t
    use domain_mod,            only : domain_t
    use domain_factory_mod,     only : create_domain

    type(domain_t) :: domain

    class(outputer_t), allocatable :: outputer

    integer(kind=4)                    :: nh=128, nz=3, halo_width=3
    integer(kind=4)                    :: myid, np, ierr, code
    integer(kind=4)                    :: master_id = 0

    call create_domain(domain, "cube", 'A', nh, nz)

    call domain%parcomm%print('Running master process paneled output test!')

    call create_master_paneled_outputer(outputer, 'p', domain)

    call test_case_1(outputer, domain, 'master_outputer')

end subroutine test_master_paneled_output

! subroutine test_mpi_paneled_output()
!
!     use outputer_abstract_mod,       only : outputer_t
!     use outputer_factory_mod,        only : create_mpi_paneled_outputer
!     use partition_mod,               only : partition_t
!     use mpi
!
!     class(outputer_t),allocatable :: outputer
!     type(partition_t) :: partition
!
!     integer(kind=4)                    :: nh=128, nz=3, halo_width=3
!     integer(kind=4)                    :: myid, np, ierr, code
!     integer(kind=4)                    :: master_id = 0
!
!     call MPI_comm_rank(mpi_comm_world , myid, ierr)
!     call MPI_comm_size(mpi_comm_world , Np  , ierr)
!
!     if (myid==0) print*, 'Running mpi paneled output test!'
!
!     call partition%init(nh, nz, max(1,Np/6), myid, Np, strategy = 'default')
!
!     outputer = create_mpi_paneled_outputer(partition)
!
!     call test_case_1(outputer, partition, 'mpi_outputer')
!
! end subroutine test_mpi_paneled_output

subroutine test_case_1(outputer, domain, class_name)

    use domain_mod, only : domain_t

    use grid_field_mod,              only : grid_field_t
    use grid_field_factory_mod,      only : create_grid_field
    use outputer_abstract_mod,       only : outputer_t

    class(outputer_t), intent(inout) :: outputer
    type(domain_t),    intent(in)    :: domain
    character(*),      intent(in)    :: class_name

    type(grid_field_t) :: f

    character(*), parameter            :: file_name = "h.dat"
    integer(kind=4)                    :: fdunit

    integer(kind=4)                    :: halo_width=3
    integer(kind=4)                    :: nh, nz

    integer(kind=4) :: myid, np, ierr, code
    integer(kind=4) :: ts, te
    integer(kind=4) :: pn, t, i, j, k

    real(kind=4), allocatable :: bufcheck(:,:,:,:), bufin(:,:,:,:)
    real(kind=4)  err
    real(kind=4), parameter :: tolerance = 1e-16

#define __fun(i,j,k,pn) ((pn-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k)

    !start and end index of tiles belonging to the current proccesor
    ts = domain%partition%ts
    te = domain%partition%te

    !Init arrays

    nh=domain%partition%Nh; nz=domain%partition%Nz

    call create_grid_field(f, halo_width, 0, domain%mesh_p)

    do t = ts, te
        f%tile(t)%p = huge(1.0_8)
        pn = domain%partition%panel_map(t)
        do k = domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
            do j = domain%mesh_p%tile(t)%js, domain%mesh_p%tile(t)%je
                do i = domain%mesh_p%tile(t)%is, domain%mesh_p%tile(t)%ie
                    f%tile(t)%p(i,j,k) =  __fun(i,j,k,pn)
                end do
            end do
        end do
    end do

    call outputer%write(f, domain, file_name)

    if(domain%parcomm%myid == 0) then
        allocate(bufcheck(nh,nh,6,nz))
        allocate(   bufin(nh,nh,6,nz))

        do k = 1, nz
            do pn = 1, 6
                do j = 1, nh
                    do i = 1, nh
                        bufcheck(i,j,pn,k) = __fun(i,j,k,pn)
                    end do
                end do
            end do
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

#undef __fun

end subroutine test_case_1

end module test_paneled_output_mod
