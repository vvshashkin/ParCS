module NHlin_output_mod

use grid_function_mod, only : grid_function_t
use outputer_abstract_mod,       only : outputer_t

implicit none

class(outputer_t),     allocatable :: outputer

contains

subroutine init_NHlin_output(myid, master_id, np, partition)
    use partition_mod,         only : partition_t
    use exchange_factory_mod,  only : create_gather_exchange
    use outputer_factory_mod,  only : create_master_paneled_outputer


    integer(kind=4),     intent(in) :: myid, master_id, np
    type(partition_t),   intent(in) :: partition

    integer(kind=4) i

    outputer = create_master_paneled_outputer(master_id = master_id,  &
    gather_exch = create_gather_exchange(partition, master_id, myid, np), &
    partition = partition)
end subroutine init_NHlin_output

subroutine write_NHlin(stvec, partition, rec_num)

    use stvec_NHlin_mod,      only : stvec_NHlin_t
    use partition_mod,        only : partition_t

    type(stvec_NHlin_t), intent(inout) :: stvec
    type(partition_t),   intent(in)    :: partition
    integer(kind=4),     intent(in)    :: rec_num

    integer(kind=4) :: ts, te

    ts = partition%ts; te = partition%te

    call outputer%write(stvec%h(ts:te), ts, te, partition, "h.dat", rec_num)
    call outputer%write(stvec%u(ts:te), ts, te, partition, "u.dat", rec_num)
    call outputer%write(stvec%v(ts:te), ts, te, partition, "v.dat", rec_num)

end subroutine write_NHlin

end module NHlin_output_mod
