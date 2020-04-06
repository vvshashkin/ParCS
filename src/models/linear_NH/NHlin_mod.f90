module NHlin_mod

use parameters_NHlin_mod,    only : parameters_NHlin_t, init_NHlin_parameters
use stvec_NHlin_mod,         only : stvec_NHlin_t, init_stvec_NHlin
use operator_NHlin_mod,      only : operator_NHlin_t
use timescheme_abstract_mod, only : timescheme_abstract_t
use diag_NHlin_mod,          only : init_NHlin_diag_mod

implicit none

!run & output controls
integer(kind=4)            :: nstep = 10
integer(kind=4)            :: nzap  = 1
namelist /ctr/ nstep, nzap

!internal controls
integer(kind=4)            :: master_id = 0

!model setup
class(parameters_NHlin_t), allocatable    :: params
type(stvec_NHlin_t)                       :: stvec
type(operator_NHlin_t)                    :: operator
class(timescheme_abstract_t), allocatable :: time_scheme

contains

subroutine init_NHlin_model()

    use mpi
    use cmd_args_mod,             only : cmd_arg_t, get_cmd_args
    use namelist_read_mod,        only : read_namelist_as_str
    use NHlin_output_mod,         only : init_NHlin_output
    use NHlin_initial_cond_mod,   only : set_NHlin_initial_conditions
    use operator_NHlin_mod,       only : init_NHlin_operator
    use tscheme_NHlin_mod,        only : init_tscheme_NHlin

    integer(kind=4) myid, Np, ierr

    type(cmd_arg_t), allocatable :: cmd_args(:)
    integer(kind=4) nargs

    character(:), allocatable :: namelist_str

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    if(myid == master_id) print *, "linear non-hydrostatic model"

    call get_cmd_args(cmd_args, nargs)

    if(nargs > 1) then
        call read_namelist_as_str(namelist_str,cmd_args(2)%str, myid, master_id = master_id)
    end if
    if(allocated(namelist_str)) then
        read(namelist_str, ctr)
    end if

    call init_NHlin_parameters(params, namelist_str, myid, Np, master_id)

    call init_NHlin_operator(operator, params, master_id, myid, Np, namelist_str)

    call init_stvec_NHlin(stvec, params%ts, params%te, params%tiles,  &
                          params%halo_width)

    call init_NHlin_output(myid, master_id, Np, params%partition)

    call init_NHlin_diag_mod()

    call set_NHlin_initial_conditions(stvec, namelist_str, params, myid, master_id)

    call operator%ext_halo(stvec, params%ts, params%te)

    call init_tscheme_NHlin(time_scheme,operator,params,stvec,myid,master_id,np,namelist_str)

    if(myid == master_id) then
        print *, "-----------------------------------------"
        print *, "|", nstep, "time steps will be performed"
        print *, "|", "output each ", nzap, "steps"
        print *, "========================================"
        print *, "----------execution:--------------------"
    end if

end subroutine init_NHlin_model

subroutine run_NHlin_model()
    use mpi
    use NHlin_output_mod, only : write_NHlin!, print_NHlin_diag
    use diag_NHlin_mod,   only : print_NHlin_diag

    integer(kind=4) myid, Np, ierr
    integer(kind=4) istep, ind, irec

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    call write_NHlin(stvec, params%partition, 1)
    call print_NHlin_diag(0, stvec, params, myid, master_id)

    irec = 2
    do istep = 1, nstep
        call time_scheme%step(stvec, params, params%dt)
        if(mod(istep,nzap) == 0) then
            call write_NHlin(stvec, params%partition,irec)
            call print_NHlin_diag(istep, stvec, params, myid, master_id)
            irec = irec+1
        end if
    end do
end subroutine run_NHlin_model

end module NHlin_mod
