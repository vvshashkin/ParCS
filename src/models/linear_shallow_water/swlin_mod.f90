module swlin_mod

implicit none

integer(kind=4)        :: nx    = 64
namelist /dims/ nx

real(kind=8)           :: H0    = 1e3_8
real(kind=8)           :: dt    = 100._8
namelist /dyn/ H0, dt

integer(kind=4)        :: nstep = 1000
integer(kind=4)        :: nzap  = 100
namelist /ctr/ nstep, nzap

contains

subroutine swlin_model_main()

    use mpi
    use cmd_args_mod,             only : cmd_arg_t, get_cmd_args
    use namelist_read_mod,        only : read_namelist_as_str

    integer(kind=4) myid, Np, ierr

    type(cmd_arg_t), allocatable :: cmd_args(:)
    integer(kind=4) nargs

    character(:), allocatable :: namelist_str

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)


    call get_cmd_args(cmd_args, nargs)

    if(nargs > 1) then
        call read_namelist_as_str(namelist_str,cmd_args(2)%str, myid, master_id = 0)
    end if
    if(allocated(namelist_str)) then
        read(namelist_str, dims)
        read(namelist_str, dyn)
        read(namelist_str, ctr)
    end if
    print *, "nx=", nx
    print *, "H0=", H0
    print *, "dt=", dt
    print *, "nstep=", nstep
    print *, "nzap=", nzap

end subroutine swlin_model_main

end module swlin_mod
