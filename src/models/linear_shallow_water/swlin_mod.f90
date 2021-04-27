module swlin_mod

use parameters_swlin_mod,    only : parameters_swlin_t, init_swlin_parameters
use stvec_swlin_mod,         only : stvec_swlin_t
use operator_swlin_mod,      only : operator_swlin_t
use timescheme_abstract_mod, only : timescheme_abstract_t
use diag_swlin_mod,          only : init_swlin_diag_mod

implicit none

!run & output controls
integer(kind=4)            :: nstep = 10
integer(kind=4)            :: nzap  = 1
namelist /ctr/ nstep, nzap

!internal controls
integer(kind=4)            :: master_id = 0

!model setup
class(parameters_swlin_t), allocatable    :: params
type(stvec_swlin_t)                       :: stvec
type(operator_swlin_t)                    :: operator
class(timescheme_abstract_t), allocatable :: time_scheme

contains

subroutine init_swlin_model()

    use mpi
    use cmd_args_mod,               only : cmd_arg_t, get_cmd_args
    use namelist_read_mod,          only : read_namelist_as_str
    use swlin_output_mod,           only : init_swlin_output
    use swlin_initial_cond_mod,     only : set_swlin_initial_conditions
    use swlin_operator_factory_mod, only : create_swlin_operator
    use stvec_swlin_factory_mod,    only : create_stvec_swlin
    use rk4_mod,                    only : rk4_t, init_rk4

    integer(kind=4) myid, Np, ierr

    type(cmd_arg_t), allocatable :: cmd_args(:)
    integer(kind=4) nargs

    character(:), allocatable :: namelist_str

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    if(myid == master_id) print *, "linear shallow water model"

    call get_cmd_args(cmd_args, nargs)

    if(nargs > 1) then
        call read_namelist_as_str(namelist_str,cmd_args(2)%str, myid, master_id = master_id)
    end if
    if(allocated(namelist_str)) then
        read(namelist_str, ctr)
    end if

    call init_swlin_parameters(params, namelist_str, myid, Np, master_id)

    call create_stvec_swlin(stvec,  params%halo_width, params%partition)

    call init_swlin_output(myid, master_id, Np, params%partition)

    call init_swlin_diag_mod()

    call set_swlin_initial_conditions(stvec, namelist_str, params%ts, params%te, &
                                      params%mesh, myid, master_id)

    call create_swlin_operator(operator, params, master_id, myid, Np, namelist_str)
    call operator%ext_halo(stvec, params%ts, params%te)

    time_scheme = init_rk4(operator, stvec)

    print *, "-----------------------------------------"
    print *, "|", nstep, "time steps will be performed"
    print *, "|", "output each ", nzap, "steps"
    print *, "========================================"
    print *, "----------execution:--------------------"

end subroutine init_swlin_model

subroutine run_swlin_model()
    use mpi
    use swlin_output_mod, only : write_swlin!, print_swlin_diag
    use diag_swlin_mod,   only : print_swlin_diag

    integer(kind=4) myid, Np, ierr
    integer(kind=4) istep, ind, irec

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    call write_swlin(stvec, params%partition, 1)

    irec = 2
    do istep = 1, nstep
        call time_scheme%step(stvec, params, params%dt)
        if(mod(istep,nzap) == 0) then
            call write_swlin(stvec, params%partition,irec)
            call print_swlin_diag(istep, stvec, params, myid, master_id)
            irec = irec+1
        end if
    end do
end subroutine run_swlin_model

end module swlin_mod
