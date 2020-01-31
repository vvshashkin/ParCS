module swlin_mod

use partition_mod,   only : partition_t
use stvec_swlin_mod, only : stvec_swlin_t, init_stvec_swlin
use mesh_mod,        only : mesh_t

implicit none

!dimensions
integer(kind=4)            :: nx    = 64
integer(kind=4), parameter :: nz    = 1 !shallow water!
namelist /dims/ nx

!model specific constants
real(kind=8)               :: H0    = 1e3_8
real(kind=8)               :: dt    = 100._8
namelist /dyn/ H0, dt

!run & output controls
integer(kind=4)            :: nstep = 1000
integer(kind=4)            :: nzap  = 100
integer(kind=4)            :: test_case_num = 1
namelist /ctr/ nstep, nzap, test_case_num

!internal controls
integer(kind=4)            :: master_id = 0

integer(kind=4)            :: halo_width = 8
type(partition_t)          :: partition
type(stvec_swlin_t)        :: stvec
type(mesh_t), allocatable  :: mesh(:)

contains

subroutine init_swlin_model()

    use mpi
    use cmd_args_mod,             only : cmd_arg_t, get_cmd_args
    use namelist_read_mod,        only : read_namelist_as_str
    use mesh_factory_mod,         only : create_equiangular_mesh
    use swlin_output_mod,         only : init_swlin_output
    use swlin_initial_cond_mod,   only : set_swlin_initial_conditions

    integer(kind=4) myid, Np, ierr

    type(cmd_arg_t), allocatable :: cmd_args(:)
    integer(kind=4) nargs

    character(:), allocatable :: namelist_str

    integer(kind=4) ts, te
    integer(kind=4), allocatable :: is(:), ie(:), js(:), je(:)
    integer(kind=4), allocatable :: ks(:), ke(:), panel_ind(:)
    integer(kind=4) ind

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    if(myid == master_id) print *, "linear shallow water model"

    call get_cmd_args(cmd_args, nargs)

    if(nargs > 1) then
        call read_namelist_as_str(namelist_str,cmd_args(2)%str, myid, master_id = master_id)
    end if
    if(allocated(namelist_str)) then
        read(namelist_str, dims)
        read(namelist_str, dyn)
        read(namelist_str, ctr)
    end if
    if(myid == master_id) then
        print *, "---Model parameters:"
        print *, "nx=", nx
        print *, "H0=", H0
        print *, "dt=", dt
        print *, "nstep=", nstep
        print *, "nzap=", nzap
        print *, "--------------------"
    end if

    call partition%init(nx, nz, max(1,Np/6), Np, strategy = 'default')

    ts = findloc(partition%proc_map, myid, dim=1)
    te = findloc(partition%proc_map, myid, back = .true., dim=1)

    panel_ind = partition%tile(ts:te)%panel_number
    is = partition%tile(ts:te)%is; ie = partition%tile(ts:te)%ie
    js = partition%tile(ts:te)%js; je = partition%tile(ts:te)%je
    ks = partition%tile(ts:te)%ks; ke = partition%tile(ts:te)%ke
    call init_stvec_swlin(stvec, ts, te, panel_ind, is, ie, js,   &
                          je, ks, ke, halo_width)

    if(myid == master_id) then
        allocate(mesh(partition%num_tiles*6))
    else
        allocate(mesh(ts:te))
    end if
    do ind = lbound(mesh, dim=1),ubound(mesh, dim=1)
        call create_equiangular_mesh(mesh(ind), partition%tile(ind)%is, partition%tile(ind)%ie, &
                                                partition%tile(ind)%js, partition%tile(ind)%je, &
                                                partition%tile(ind)%ks, partition%tile(ind)%ke, &
                                                nx, halo_width, partition%tile(ind)%panel_number)
    end do
    call init_swlin_output(myid, master_id, Np, partition,                &
                           lbound(mesh, dim=1),ubound(mesh, dim=1),mesh)

    call set_swlin_initial_conditions(stvec, test_case_num, ts,te, mesh(ts:te))

end subroutine init_swlin_model

subroutine run_swlin_model()
    use mpi
    use swlin_output_mod,         only : write_swlin

    integer(kind=4) myid, Np, ierr
    integer(kind=4) ind
    type(stvec_swlin_t) v2

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    call write_swlin(myid, master_id, stvec%ts, stvec%te, stvec,  &
                     lbound(mesh, dim=1),ubound(mesh, dim=1), mesh, 1)
    print *, "wr"
    call v2.copy(stvec)
    call v2.add(stvec,1._8, 0.5_8)

    do ind=v2%ts, v2%te
        print *, "hmax", ind, maxval(v2%h(ind)%p(:,:,:))
    end do

    call write_swlin(myid, master_id, stvec%ts, stvec%te, v2,  &
                     lbound(mesh, dim=1),ubound(mesh, dim=1), mesh, 2)
end subroutine run_swlin_model

end module swlin_mod
