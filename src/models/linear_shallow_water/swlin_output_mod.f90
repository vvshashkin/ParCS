module swlin_output_mod

use grid_function_mod, only : grid_function_t
use outputer_abstract_mod,       only : outputer_t

implicit none

class(outputer_t),     allocatable :: outputer

contains

subroutine init_swlin_output(myid, master_id, np, partition)
    use partition_mod,         only : partition_t
    use exchange_factory_mod,  only : create_gather_exchange
    use outputer_factory_mod,  only : create_master_paneled_outputer


    integer(kind=4),     intent(in) :: myid, master_id, np
    type(partition_t),   intent(in) :: partition

    integer(kind=4) i

    outputer = create_master_paneled_outputer(master_id = master_id,  &
    gather_exch = create_gather_exchange(partition, master_id, myid, np), &
    partition = partition)
end subroutine init_swlin_output

subroutine write_swlin(stvec, partition, rec_num)

    use stvec_swlin_mod,      only : stvec_swlin_t
    use partition_mod,        only : partition_t

    type(stvec_swlin_t), intent(inout) :: stvec
    type(partition_t),   intent(in)    :: partition
    integer(kind=4),     intent(in)    :: rec_num

    integer(kind=4) :: ts, te

    ts = partition%ts; te = partition%te

    call outputer%write(stvec%h(ts:te), ts, te, partition, "h.dat", rec_num)
    call outputer%write(stvec%u(ts:te), ts, te, partition, "u.dat", rec_num)
    call outputer%write(stvec%v(ts:te), ts, te, partition, "v.dat", rec_num)

end subroutine write_swlin

subroutine print_swlin_diag(stvec,ts, te, mesh, myid, master_id, istep)
    use stvec_swlin_mod,      only : stvec_swlin_t
    use mesh_mod,             only : mesh_t
    use mpi

    type(stvec_swlin_t), intent(in) :: stvec
    integer(kind=4),     intent(in) :: ts, te
    type(mesh_t),        intent(in) :: mesh(ts:te)
    integer(kind=4),     intent(in) :: myid, master_id
    integer(kind=4),     intent(in) :: istep

    real(kind=8) hmin, hmax, gl_hmin, gl_hmax
    integer(kind=4) ierr
    integer(kind=4) is, ie, js, je, ind

    ind = stvec%ts
    is  = stvec%h(ind)%is
    js  = stvec%h(ind)%js
    hmin = stvec%h(ind)%p(is,js,1)
    hmax = hmin

    do ind = stvec%ts, stvec%te
        is  = stvec%h(ind)%is;    ie  = stvec%h(ind)%ie
        js  = stvec%h(ind)%js;    je  = stvec%h(ind)%je
        hmin = min(hmin,minval(stvec%h(ind)%p(is:ie,js:je,1)))
        hmax = max(hmax,maxval(stvec%h(ind)%p(is:ie,js:je,1)))
    end do

    call mpi_reduce(hmin,gl_hmin,1,MPI_DOUBLE,MPI_MIN,master_id,MPI_COMM_WORLD,ierr)
    call mpi_reduce(hmax,gl_hmax,1,MPI_DOUBLE,MPI_MAX,master_id,MPI_COMM_WORLD,ierr)

    if(myid == master_id) then
        print "(A,I7,A,2E15.7)", "T = ", istep, " H MIN/MAX = ", gl_hmin, gl_hmax
    end if

end subroutine print_swlin_diag

end module swlin_output_mod
