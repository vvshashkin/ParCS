module NHlin_output_mod

use partition_mod,         only : partition_t
use grid_function_mod,     only : grid_function_t
use outputer_abstract_mod, only : outputer_t

implicit none

class(outputer_t),     allocatable :: outputer, outputer_w
type(partition_t) partition_w


contains

subroutine init_NHlin_output(myid, master_id, np, partition)
    use partition_mod,         only : partition_t
    use exchange_factory_mod,  only : create_gather_exchange
    use outputer_factory_mod,  only : create_master_paneled_outputer


    integer(kind=4),     intent(in) :: myid, master_id, np
    type(partition_t),   intent(in) :: partition

    integer(kind=4) i

    partition_w = partition
    do i=1,partition%num_tiles*partition%num_panels
        partition_w%tile(i)%ks = partition_w%tile(i)%ks-1
    end do

    outputer = create_master_paneled_outputer(master_id = master_id,  &
    gather_exch = create_gather_exchange(partition, master_id, myid, np), &
    partition = partition)

    outputer_w = create_master_paneled_outputer(master_id = master_id,  &
    gather_exch = create_gather_exchange(partition_w, master_id, myid, np), &
    partition = partition_w)
end subroutine init_NHlin_output

subroutine write_NHlin(stvec, partition, rec_num)

    use stvec_NHlin_mod,      only : stvec_NHlin_t
    use partition_mod,        only : partition_t

    type(stvec_NHlin_t), intent(inout) :: stvec
    type(partition_t),   intent(in)    :: partition
    integer(kind=4),     intent(in)    :: rec_num

    integer(kind=4) :: ts, te

    ts = partition%ts; te = partition%te

    !call outputer%write(stvec%prex(ts:te), ts, te, partition, "prex.dat", rec_num)
    !call outputer%write(stvec%u(ts:te), ts, te, partition, "u.dat", rec_num)
    !call outputer%write(stvec%v(ts:te), ts, te, partition, "v.dat", rec_num)
    call outputer_w%write(stvec%theta(ts:te), ts, te, partition_w, "theta.dat", rec_num)

end subroutine write_NHlin

end module NHlin_output_mod
